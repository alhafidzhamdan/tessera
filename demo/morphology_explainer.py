"""Annotated morphology cheat-sheet: what each feature measures, drawn on
geometric diagrams and on real 7-AAD melanoma nuclei with their measured values.

    python demo/morphology_explainer.py
"""
import os
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.patches import Ellipse, FancyArrowPatch, Polygon, Rectangle

sys.path.insert(0, "src")
from tessera.core import PairedData  # noqa: E402

OUT = "data/processed/morphology_features.pdf"
BUNDLE = "data/processed/melanoma_7aad.tessera"


def representative(obs, ct, feat="area"):
    """cell id nearest the median 'feat' for a cell type."""
    sub = obs[(obs["cell_type"] == ct) & obs["area"].notna()]
    med = sub[feat].median()
    return (sub[feat] - med).abs().idxmin()


def main():
    p = PairedData.load(BUNDLE)
    obs = p.adata.obs
    ink, muted, acc = "#14242A", "#697B81", "#0EA5B7"

    with PdfPages(OUT) as pdf:
        fig = plt.figure(figsize=(15, 9.2), facecolor="white")
        fig.suptitle("Tessera — nuclear morphology features: what they measure",
                     fontsize=15, y=0.98, color=ink)

        # ---- A. SIZE ---------------------------------------------------- #
        ax = fig.add_axes([0.03, 0.55, 0.29, 0.36]); _clean(ax)
        ax.set_title("Size", loc="left", fontsize=12, color=ink, weight="bold")
        el = Ellipse((0.5, 0.5), 0.6, 0.42, angle=22, facecolor="#cbe7e4",
                     edgecolor="#2a8f86", lw=1.6, alpha=.9)
        ax.add_patch(el)
        # major / minor axes
        _axis_line(ax, (0.5, 0.5), 22, 0.30, acc, "major axis")
        _axis_line(ax, (0.5, 0.5), 112, 0.21, "#C2543A", "minor axis")
        # equivalent-diameter circle (same area)
        r = np.sqrt(0.30 * 0.21) / 1.0  # geo-mean radius ~ equiv radius
        ax.add_patch(plt.Circle((0.5, 0.5), r, fill=False, ls=(0, (4, 3)),
                                ec=muted, lw=1.4))
        ax.text(0.5, 0.05,
                "area = filled pixels · equivalent diameter = circle of equal area",
                ha="center", fontsize=8.4, color=muted)

        # ---- B. SHAPE / CONTOUR ---------------------------------------- #
        ax = fig.add_axes([0.36, 0.55, 0.30, 0.36]); _clean(ax)
        ax.set_title("Shape / contour", loc="left", fontsize=12, color=ink, weight="bold")
        # a lobulated (concave) nucleus
        t = np.linspace(0, 2 * np.pi, 220)
        rr = 0.30 + 0.06 * np.sin(3 * t) - 0.05 * (np.cos(t) > 0.6)
        px, py = 0.5 + rr * np.cos(t), 0.52 + rr * 0.8 * np.sin(t)
        ax.fill(px, py, color="#e7d7ef", ec="#8455a6", lw=1.5, alpha=.9, zorder=2)
        # convex hull (solidity)
        from scipy.spatial import ConvexHull
        pts = np.column_stack([px, py]); hull = ConvexHull(pts)
        hp = pts[np.append(hull.vertices, hull.vertices[0])]
        ax.plot(hp[:, 0], hp[:, 1], ls=(0, (5, 3)), color="#8455a6", lw=1.3, zorder=3)
        # bounding box (extent)
        x0, x1, y0, y1 = px.min(), px.max(), py.min(), py.max()
        ax.add_patch(Rectangle((x0, y0), x1 - x0, y1 - y0, fill=False, ls=":",
                               ec=muted, lw=1.2, zorder=1))
        # fitted ellipse (eccentricity) + orientation
        ax.add_patch(Ellipse((0.5, 0.52), 0.52, 0.40, angle=8, fill=False,
                             ec=acc, lw=1.3, alpha=.8, zorder=3))
        ax.add_patch(FancyArrowPatch((0.5, 0.52), (0.5 + 0.26 * np.cos(np.radians(8)),
                     0.52 + 0.26 * np.sin(np.radians(8))), color=acc, lw=1.6,
                     arrowstyle="-|>", mutation_scale=11, zorder=4))
        ax.text(0.5, 0.05,
                "solidity = area ÷ convex hull (dashed) · extent = area ÷ box (dotted)\n"
                "eccentricity & orientation = fitted ellipse (teal)",
                ha="center", fontsize=8.4, color=muted)

        # ---- C. INTENSITY / CHROMATIN ---------------------------------- #
        ax = fig.add_axes([0.70, 0.55, 0.27, 0.36]); ax.axis("off")
        ax.set_title("Intensity / chromatin", loc="left", fontsize=12, color=ink, weight="bold")
        tcid = representative(obs, "tumour_1")
        crop = p.get_crop(tcid)[0]
        axim = fig.add_axes([0.70, 0.60, 0.155, 0.27])
        axim.imshow(crop, cmap="magma", vmin=0, vmax=255); axim.set_xticks([]); axim.set_yticks([])
        for s in axim.spines.values(): s.set_color("#a01646")
        axh = fig.add_axes([0.865, 0.60, 0.10, 0.27]); axh.set_facecolor("white")
        vals = crop[crop > 20]
        axh.hist(vals, bins=24, orientation="horizontal", color="#a01646", alpha=.8)
        axh.set_xticks([]); axh.set_ylim(0, 255)
        axh.set_ylabel("7-AAD signal", fontsize=8, color=muted)
        axh.tick_params(labelsize=7, colors=muted)
        for sp in ("top", "right", "bottom"): axh.spines[sp].set_visible(False)
        fig.text(0.70, 0.565,
                 "mean = chromatin density · std = texture / coarseness · max = nucleoli",
                 fontsize=8.4, color=muted)

        # ---- D. REAL NUCLEI: concept -> numbers ------------------------ #
        fig.text(0.03, 0.47, "The same features on real 7-AAD melanoma nuclei",
                 fontsize=12, color=ink, weight="bold")
        examples = ["tumour_1", "mono-mac", "T_CD8", "plasma"]
        feats = [("area", "Area", "px"), ("equivalent_diameter_area", "Equiv. diam.", "px"),
                 ("eccentricity", "Eccentricity", ""), ("solidity", "Solidity", ""),
                 ("intensity_mean", "Signal mean", ""), ("intensity_std", "Signal std", "")]
        n = len(examples); w = 0.225; gap = 0.018
        for j, ct in enumerate(examples):
            cid = representative(obs, ct)
            crop = p.get_crop(cid)[0]
            x = 0.03 + j * (w + gap)
            axc = fig.add_axes([x, 0.20, 0.11, 0.22])
            axc.imshow(crop, cmap="gray", vmin=0, vmax=255); axc.set_xticks([]); axc.set_yticks([])
            axc.set_title(ct, fontsize=10, color=ink)
            axt = fig.add_axes([x + 0.115, 0.20, w - 0.115, 0.22]); axt.axis("off")
            row = obs.loc[cid]
            for k, (f, lab, u) in enumerate(feats):
                v = row[f]
                s = "—" if v != v else (f"{v:.0f}" if abs(v) >= 20 else f"{v:.2f}")
                axt.text(0, 1 - k * 0.17, lab, fontsize=8, color=muted, va="top")
                axt.text(1, 1 - k * 0.17, f"{s} {u}".strip(), fontsize=8.5, color=ink,
                         va="top", ha="right", family="monospace")

        fig.text(0.03, 0.14,
                 "Reading across: tumour nuclei are large, less solid (irregular) and "
                 "have coarser chromatin (higher signal std); lymphocytes are small, round "
                 "and uniformly condensed — the classic cytopathology axes, computed per "
                 "nucleus and paired to its RNA + ATAC profile.",
                 fontsize=9.5, color=ink, wrap=True)

        pdf.savefig(fig, facecolor="white"); plt.close(fig)
    print("wrote", OUT)


def _clean(ax):
    ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.set_aspect("equal")
    ax.set_xticks([]); ax.set_yticks([])
    for s in ax.spines.values(): s.set_visible(False)


def _axis_line(ax, c, ang, half, color, label):
    a = np.radians(ang); dx, dy = half * np.cos(a), half * np.sin(a)
    ax.plot([c[0] - dx, c[0] + dx], [c[1] - dy, c[1] + dy], color=color, lw=2)
    ax.text(c[0] + dx, c[1] + dy, " " + label, fontsize=8, color=color, va="center")


if __name__ == "__main__":
    main()
