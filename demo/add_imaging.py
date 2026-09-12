"""Wire synthetic 7-AAD nuclear imaging into the melanoma bundle: register a
nuclear-stain image to each nucleus's slide-tags (x,y), extract real morphology,
save a new bundle, and render a QC figure.

    python demo/add_imaging.py
"""
import os
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages

sys.path.insert(0, "src")
from tessera.core import PairedData  # noqa: E402
from tessera.ingest import synthesize_nuclear_image  # noqa: E402

IN = "data/processed/melanoma.tessera"
OUT = "data/processed/melanoma_7aad.tessera"
FIG = "data/processed/tessera_7aad_qc.pdf"
STAIN = "7-AAD"


def main():
    p = PairedData.load(IN)
    p, img, labels = synthesize_nuclear_image(
        p, stain=STAIN, px_per_um=1.0, crop_size=64,
        crops_path="data/processed/_crops_7aad.zarr", seed=1, return_image=True,
    )
    print(repr(p))
    print("morphology features:", p.morphology_features)
    print("\nMean morphology by cell type (REAL, from 7-AAD image):")
    print(p.catalogue_table().round(2).to_string(index=False))
    p.save(OUT)
    print(f"\nsaved -> {OUT}")
    _figure(p, img, FIG)
    print(f"figure -> {FIG}")


def _figure(p, img, path):
    a = p.adata
    obs = a.obs
    cats = list(obs["cell_type"].cat.categories)
    cmap = plt.get_cmap("tab20")
    color = {c: cmap(i % 20) for i, c in enumerate(cats)}
    # a 7-AAD-like far-red lookup for the overview
    from matplotlib.colors import LinearSegmentedColormap
    aad = LinearSegmentedColormap.from_list("aad", ["#05010a", "#3a0a1e", "#a01646", "#ff5a7a", "#ffd7de"])

    with PdfPages(path) as pdf:
        fig = plt.figure(figsize=(15, 9.5), facecolor="white")

        # overview (downsampled) with a zoom inset on a dense tumour nest
        ax = fig.add_axes([0.04, 0.30, 0.44, 0.66]); ax.set_facecolor("white")
        ax.imshow(img[::4, ::4], cmap=aad, vmin=0, vmax=255, interpolation="nearest")
        ax.set_title(f"Synthetic {STAIN} nuclear image · registered to slide-tags (x,y)", fontsize=11)
        ax.set_xticks([]); ax.set_yticks([])
        # inset
        xy = a.obsm["spatial"]; pad = 64
        cx = int(np.median(xy[:, 0]) - xy[:, 0].min() + pad)
        cy = int(np.median(xy[:, 1]) - xy[:, 1].min() + pad)
        h = 260
        sub = img[max(0, cy - h):cy + h, max(0, cx - h):cx + h]
        axi = fig.add_axes([0.50, 0.62, 0.16, 0.34]); axi.imshow(sub, cmap=aad, vmin=0, vmax=255)
        axi.set_title("zoom", fontsize=9); axi.set_xticks([]); axi.set_yticks([])
        for s in axi.spines.values(): s.set_color("#a01646")

        # crop gallery: 3 nuclei per cell type
        axg = fig.add_axes([0.50, 0.30, 0.16, 0.30]); axg.axis("off")
        axg.set_title("Per-nucleus 7-AAD crops", fontsize=9, loc="left")
        rows = len(cats)
        for r, c in enumerate(cats):
            ids = obs.index[(obs["cell_type"] == c).to_numpy()][:3]
            for k, cid in enumerate(ids):
                crop = p.get_crop(cid)[0]
                s = axg.inset_axes([0.30 + k * 0.22, 1 - (r + 1) / rows, 0.2, 0.92 / rows])
                s.imshow(crop, cmap="gray", vmin=0, vmax=255); s.set_xticks([]); s.set_yticks([])
            axg.text(0, 1 - (r + 0.5) / rows, c, va="center", fontsize=6.5,
                     color=color[c], transform=axg.transAxes)

        # real morphology by cell type: area + eccentricity
        for j, feat in enumerate(["area", "eccentricity"]):
            ax = fig.add_axes([0.72, 0.56 - j * 0.48, 0.25, 0.38]); ax.set_facecolor("white")
            data = [obs.loc[obs["cell_type"] == c, feat].dropna().to_numpy() for c in cats]
            bp = ax.boxplot(data, vert=False, patch_artist=True, widths=0.6)
            for patch, c in zip(bp["boxes"], cats):
                patch.set_facecolor(color[c]); patch.set_alpha(.75)
            for med in bp["medians"]:
                med.set_color("black")
            ax.set_yticklabels(cats, fontsize=7)
            ax.set_xlabel(f"nucleus {feat}" + (" (px)" if feat == "area" else ""))
            ax.set_title(f"Real morphology · {feat}", fontsize=10)

        fig.suptitle("Tessera — 7-AAD imaging wired into melanoma slide-tags: morphology is now real",
                     fontsize=14, y=0.99)
        pdf.savefig(fig, facecolor="white"); plt.close(fig)


if __name__ == "__main__":
    main()
