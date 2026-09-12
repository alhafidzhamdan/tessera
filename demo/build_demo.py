"""Build the synthetic paired dataset, print the catalogue, and render a
first-look figure proving the morphology <-> multi-omic pairing end to end.

    python demo/build_demo.py [outdir]
"""
import os
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
from tessera.ingest import make_synthetic  # noqa: E402


def main(outdir="demo_out"):
    os.makedirs(outdir, exist_ok=True)
    pd_obj = make_synthetic(crops_path=os.path.join(outdir, "crops.zarr"), seed=0)
    print(repr(pd_obj))

    cat = pd_obj.catalogue()
    print("\nCatalogue summary")
    print("-" * 60)
    for k, v in cat.items():
        print(f"  {k}: {v}")
    print("\nPer-cell-type table")
    print(pd_obj.catalogue_table().to_string(index=False))

    pd_obj.save(os.path.join(outdir, "melanoma_demo.tessera"))
    print(f"\nSaved paired dataset -> {outdir}/melanoma_demo.tessera")

    _figure(pd_obj, os.path.join(outdir, "tessera_first_look.pdf"))
    print(f"Figure -> {outdir}/tessera_first_look.pdf")


def _figure(pd_obj, path):
    adata = pd_obj.adata
    obs = adata.obs
    cats = list(obs["cell_type"].cat.categories)
    cmap = plt.get_cmap("tab10")
    color = {c: cmap(i % 10) for i, c in enumerate(cats)}
    xy = adata.obsm["spatial"]

    with PdfPages(path) as pdf:
        fig = plt.figure(figsize=(14, 9), facecolor="white")

        # (1) spatial tissue map
        ax = fig.add_axes([0.04, 0.55, 0.42, 0.4])
        ax.set_facecolor("white")
        for c in cats:
            m = (obs["cell_type"] == c).to_numpy()
            ax.scatter(xy[m, 0], xy[m, 1], s=6, color=color[c], label=c, linewidths=0)
        ax.invert_yaxis()
        ax.set_title("Spatial tissue map (slide-tags coords)")
        ax.set_xticks([]); ax.set_yticks([])
        ax.legend(markerscale=2, fontsize=6, ncol=2, loc="upper right")

        # (2) UMAP
        ax = fig.add_axes([0.54, 0.55, 0.42, 0.4])
        ax.set_facecolor("white")
        um = adata.obsm["X_umap"]
        for c in cats:
            m = (obs["cell_type"] == c).to_numpy()
            ax.scatter(um[m, 0], um[m, 1], s=6, color=color[c], linewidths=0)
        ax.set_title("Transcriptome UMAP")
        ax.set_xticks([]); ax.set_yticks([])

        # (3) nucleus gallery: 2 example nuclei per cell type
        ax_g = fig.add_axes([0.04, 0.06, 0.5, 0.42])
        ax_g.axis("off")
        ax_g.set_title("Paired nucleus crops (DAPI) by cell type", loc="left")
        n_per = 2
        rows = len(cats)
        for r, c in enumerate(cats):
            ids = obs.index[(obs["cell_type"] == c).to_numpy()][:n_per]
            for k, cid in enumerate(ids):
                crop = pd_obj.get_crop(cid)[0]
                sub = ax_g.inset_axes(
                    [0.18 + k * 0.11, 1 - (r + 1) / rows, 0.1, 0.9 / rows]
                )
                sub.imshow(crop, cmap="gray", vmin=0, vmax=255)
                sub.set_xticks([]); sub.set_yticks([])
            ax_g.text(0.0, 1 - (r + 0.5) / rows, c, va="center", fontsize=8,
                      color=color[c], transform=ax_g.transAxes)

        # (4) morphology <-> omics link: nucleus area by cell type
        ax = fig.add_axes([0.62, 0.06, 0.34, 0.42])
        ax.set_facecolor("white")
        data = [obs.loc[obs["cell_type"] == c, "area"].to_numpy() for c in cats]
        bp = ax.boxplot(data, vert=False, patch_artist=True, widths=0.6)
        for patch, c in zip(bp["boxes"], cats):
            patch.set_facecolor(color[c]); patch.set_alpha(0.7)
        for med in bp["medians"]:
            med.set_color("black")
        ax.set_yticklabels(cats, fontsize=8)
        ax.set_xlabel("nucleus area (px)")
        ax.set_title("Morphology feature by transcriptomic cell type")

        fig.suptitle(
            "Tessera - synthetic slide-tags melanoma: morphology paired to multi-omics",
            fontsize=13, y=0.99,
        )
        pdf.savefig(fig, facecolor="white")
        plt.close(fig)


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "demo_out")
