"""Interactive napari viewer for a Tessera paired dataset.

Browse nuclei on the spatial tissue map or the UMAP, colour by cell type or by
any gene / ATAC peak, and click a nucleus to see its image crop (when imaging is
paired) alongside its top expressed genes and modality depth.

    python -m tessera.viewer.napari_view data/processed/melanoma.tessera

or from Python::

    from tessera.core import PairedData
    from tessera.viewer.napari_view import view
    view(PairedData.load("data/processed/melanoma.tessera"))
"""
from __future__ import annotations

import sys
from typing import Optional

import numpy as np

from ..core import PairedData

_TAB20 = None


def _cell_type_colors(categories):
    import matplotlib.pyplot as plt

    cmap = plt.get_cmap("tab20")
    return {c: np.array(cmap(i % 20)) for i, c in enumerate(categories)}


def _feature_vector(paired: PairedData, name: str) -> Optional[np.ndarray]:
    """Dense per-nucleus vector for a gene (RNA) or peak (ATAC), else None."""
    a = paired.adata
    if name in a.var_names:
        col = a[:, name].X
        return np.asarray(col.todense()).ravel() if hasattr(col, "todense") else np.asarray(col).ravel()
    for m in paired.mods.values():
        if name in m.var_names:
            col = m[:, name].X
            return np.asarray(col.todense()).ravel() if hasattr(col, "todense") else np.asarray(col).ravel()
    return None


class _Panel:
    """Qt dock widget: controls + per-nucleus profile (matplotlib)."""

    def __init__(self, paired: PairedData, viewer, points):
        from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg
        from matplotlib.figure import Figure
        from qtpy.QtWidgets import (
            QComboBox, QLabel, QLineEdit, QPushButton, QVBoxLayout, QWidget,
        )

        self.paired = paired
        self.viewer = viewer
        self.points = points
        self.adata = paired.adata
        self.cats = list(self.adata.obs["cell_type"].cat.categories)
        self.colors = _cell_type_colors(self.cats)

        self.widget = QWidget()
        lay = QVBoxLayout(self.widget)

        lay.addWidget(QLabel("<b>Embedding</b>"))
        self.embed = QComboBox(); self.embed.addItems(["spatial", "umap"])
        self.embed.currentTextChanged.connect(self._set_embedding)
        lay.addWidget(self.embed)

        lay.addWidget(QLabel("<b>Colour by</b>"))
        self.feat = QLineEdit(); self.feat.setPlaceholderText("cell_type, or a gene / peak")
        self.feat.setText("cell_type")
        btn = QPushButton("Apply colour")
        btn.clicked.connect(lambda: self._color_by(self.feat.text().strip()))
        lay.addWidget(self.feat); lay.addWidget(btn)

        self.info = QLabel("Click a nucleus…"); self.info.setWordWrap(True)
        lay.addWidget(self.info)

        self.fig = Figure(figsize=(3.2, 4.4), facecolor="white")
        self.canvas = FigureCanvasQTAgg(self.fig)
        lay.addWidget(self.canvas)

    # -- controls ----------------------------------------------------------- #
    def _set_embedding(self, key: str):
        coords = self.paired.spatial if key == "spatial" else self.adata.obsm["X_umap"]
        self.points.data = coords[:, ::-1]  # napari is (row, col) = (y, x)
        self.viewer.reset_view()

    def _color_by(self, name: str):
        if name in ("", "cell_type"):
            self.points.face_color = "cell_type"
            self.points.face_color_cycle = [self.colors[c] for c in self.cats]
            self.points.refresh_colors()
            return
        vec = _feature_vector(self.paired, name)
        if vec is None:
            self.info.setText(f"'{name}' not found in RNA genes or ATAC peaks.")
            return
        self.points.features = {**dict(self.points.features), name: vec}
        self.points.face_color = name
        self.points.face_colormap = "viridis"
        self.points.refresh_colors()
        self.info.setText(f"Colouring by {name} (range {vec.min():.0f}–{vec.max():.0f}).")

    # -- selection ---------------------------------------------------------- #
    def show_nucleus(self, idx: int):
        obs = self.adata.obs.iloc[idx]
        cid = self.adata.obs_names[idx]
        lines = [f"<b>{cid}</b>", f"cell_type: {obs['cell_type']}"]
        for c in ("tcr_alpha", "tcr_beta"):
            if c in obs and str(obs[c]) not in ("nan", "NA"):
                lines.append(f"{c}: {obs[c]}")
        self.info.setText("<br>".join(lines))

        self.fig.clear()
        crop = self.paired.get_crop(cid)
        ax_img = self.fig.add_subplot(211)
        if crop is not None:
            ax_img.imshow(crop[0], cmap="gray"); ax_img.set_title("nucleus crop", fontsize=8)
        else:
            ax_img.text(0.5, 0.5, "no image\n(morphology slot empty)", ha="center",
                        va="center", fontsize=8, color="0.5")
            ax_img.set_title("morphology", fontsize=8)
        ax_img.set_xticks([]); ax_img.set_yticks([])

        ax_g = self.fig.add_subplot(212)
        row = self.adata.X[idx]
        row = np.asarray(row.todense()).ravel() if hasattr(row, "todense") else np.asarray(row).ravel()
        top = np.argsort(row)[::-1][:10][::-1]
        ax_g.barh(range(len(top)), row[top], color="#4C72B0")
        ax_g.set_yticks(range(len(top)))
        ax_g.set_yticklabels(self.adata.var_names[top], fontsize=7)
        atac = self.paired.mods.get("atac")
        depth = ""
        if atac is not None:
            arow = atac.X[idx]
            n_peaks = int((arow > 0).sum()) if not hasattr(arow, "nnz") else int(arow.nnz)
            depth = f"  |  ATAC peaks: {n_peaks}"
        ax_g.set_title(f"top RNA genes{depth}", fontsize=8)
        ax_g.set_facecolor("white")
        self.fig.tight_layout()
        self.canvas.draw_idle()


def view(paired: PairedData, point_size: float = 30.0):
    """Open the napari viewer on a PairedData object. Returns the viewer."""
    import napari

    coords = paired.spatial[:, ::-1]  # (y, x) for napari
    cats = list(paired.adata.obs["cell_type"].cat.categories)
    colors = _cell_type_colors(cats)
    ct = paired.adata.obs["cell_type"].astype(str).to_numpy()

    viewer = napari.Viewer(title="Tessera — paired morphology + multi-omics")
    points = viewer.add_points(
        coords,
        name="nuclei",
        size=point_size,
        features={"cell_type": ct, "index": np.arange(paired.n_cells)},
        face_color="cell_type",
        face_color_cycle=[colors[c] for c in cats],
        border_width=0,
    )
    viewer.reset_view()

    panel = _Panel(paired, viewer, points)
    viewer.window.add_dock_widget(panel.widget, area="right", name="nucleus")

    @points.mouse_drag_callbacks.append
    def _on_click(layer, event):
        start = np.asarray(event.position)
        yield
        # treat as a click only if the pointer barely moved
        if np.linalg.norm(np.asarray(event.position) - start) > point_size:
            return
        data = np.asarray(layer.data)
        d = np.linalg.norm(data - np.asarray(event.position)[-data.shape[1]:], axis=1)
        idx = int(np.argmin(d))
        if d[idx] <= point_size:
            layer.selected_data = {idx}
            panel.show_nucleus(idx)

    return viewer


def main():
    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit(1)
    import napari

    paired = PairedData.load(sys.argv[1])
    print(repr(paired))
    view(paired)
    napari.run()


if __name__ == "__main__":
    main()
