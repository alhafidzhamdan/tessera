"""Verify the viewer's non-GL logic without opening napari's canvas:
feature lookup, colour-by, and the per-nucleus profile panel. Saves a panel
screenshot for a T cell and a tumour cell.
"""
import os
import sys
import types

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

import matplotlib

matplotlib.use("Agg")
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
from tessera.core import PairedData  # noqa: E402
from tessera.viewer.napari_view import _Panel, _feature_vector  # noqa: E402


def main():
    from qtpy.QtWidgets import QApplication

    app = QApplication.instance() or QApplication([])

    p = PairedData.load("data/processed/melanoma.tessera")
    print(repr(p))

    # dummy viewer + points (no GL): capture attribute writes
    dummy_view = types.SimpleNamespace(reset_view=lambda: None)
    dummy_pts = types.SimpleNamespace(
        data=None, features={}, face_color=None, face_color_cycle=None,
        face_colormap=None, refresh_colors=lambda: None,
    )
    panel = _Panel(p, dummy_view, dummy_pts)

    # -- feature lookup ----------------------------------------------------- #
    genes = [g for g in ("CD8A", "CD3D", "MLANA", "PMEL", "PTPRC") if g in p.adata.var_names]
    print("marker genes present:", genes)
    vec = _feature_vector(p, genes[0])
    assert vec is not None and vec.shape == (p.n_cells,)
    print(f"_feature_vector({genes[0]}) range: {vec.min():.0f}-{vec.max():.0f}, "
          f"mean {vec.mean():.2f}")
    peak = p.mods["atac"].var_names[100]
    assert _feature_vector(p, peak) is not None
    print(f"_feature_vector(peak {peak}) OK")
    assert _feature_vector(p, "NOT_A_FEATURE") is None
    print("unknown feature -> None OK")

    # -- colour-by ---------------------------------------------------------- #
    panel._color_by(genes[0])
    assert dummy_pts.face_color == genes[0] and genes[0] in dummy_pts.features
    panel._color_by("cell_type")
    assert dummy_pts.face_color == "cell_type"
    print("colour-by gene and cell_type OK")

    # -- embedding swap ----------------------------------------------------- #
    panel._set_embedding("umap")
    assert dummy_pts.data.shape == (p.n_cells, 2)
    panel._set_embedding("spatial")
    print("embedding swap OK")

    # -- profile panel screenshots ----------------------------------------- #
    obs = p.adata.obs
    examples = {}
    for ct in ("T_CD8", "tumour_1"):
        idx = int(np.where(obs["cell_type"].to_numpy() == ct)[0][0])
        panel.show_nucleus(idx)
        out = f"data/processed/panel_{ct}.png"
        panel.fig.savefig(out, dpi=130, facecolor="white", bbox_inches="tight")
        examples[ct] = out
        print(f"saved {out}")

    print("ALL PANEL CHECKS PASSED")


if __name__ == "__main__":
    main()
