"""Ingest the real SCP2176 slide-tags melanoma multiome, save a Tessera bundle,
and render a first-look figure.

    python demo/build_real.py
"""
import os
import sys
import time

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
from tessera.ingest import load_slidetags  # noqa: E402

RAW = "data/raw/SCP2176"
OUT = "data/processed/melanoma.tessera"
FIG = "data/processed/tessera_real_first_look.pdf"


def main():
    t0 = time.time()
    p = load_slidetags(RAW, load_atac=True)
    print(repr(p), f"[loaded in {time.time()-t0:.0f}s]")
    print("modality sizes:", p.catalogue()["modality_sizes"])

    p.save(OUT)
    print(f"saved bundle -> {OUT}")

    _figure(p, FIG)
    print(f"figure -> {FIG}")


def _figure(p, path):
    a = p.adata
    obs = a.obs
    cats = list(obs["cell_type"].cat.categories)
    cmap = plt.get_cmap("tab20")
    color = {c: cmap(i % 20) for i, c in enumerate(cats)}
    xy = a.obsm["spatial"]
    um = a.obsm["X_umap"]

    rna_counts = np.asarray(a.X.sum(1)).ravel()
    rna_genes = np.asarray((a.X > 0).sum(1)).ravel()
    atac = p.mods.get("atac")
    atac_acc = np.asarray((atac.X > 0).sum(1)).ravel() if atac is not None else None

    with PdfPages(path) as pdf:
        fig = plt.figure(figsize=(15, 9.5), facecolor="white")

        ax = fig.add_axes([0.04, 0.55, 0.42, 0.4]); ax.set_facecolor("white")
        for c in cats:
            m = (obs["cell_type"] == c).to_numpy()
            ax.scatter(xy[m, 0], xy[m, 1], s=8, color=color[c], label=c, linewidths=0)
        ax.invert_yaxis(); ax.set_aspect("equal")
        ax.set_title("Spatial tissue map — per-nucleus slide-tags coords (µm)")
        ax.set_xlabel("x (µm)"); ax.set_ylabel("y (µm)")
        ax.legend(markerscale=2, fontsize=6, ncol=2, loc="upper right")

        ax = fig.add_axes([0.54, 0.55, 0.42, 0.4]); ax.set_facecolor("white")
        for c in cats:
            m = (obs["cell_type"] == c).to_numpy()
            ax.scatter(um[m, 0], um[m, 1], s=8, color=color[c], linewidths=0)
        ax.set_title("Transcriptome UMAP"); ax.set_xticks([]); ax.set_yticks([])

        # cell-type counts
        ax = fig.add_axes([0.04, 0.08, 0.26, 0.36]); ax.set_facecolor("white")
        counts = obs["cell_type"].value_counts().reindex(cats)
        ax.barh(cats, counts.to_numpy(), color=[color[c] for c in cats])
        ax.invert_yaxis(); ax.set_xlabel("nuclei"); ax.set_title("Cell-type composition")
        for i, v in enumerate(counts.to_numpy()):
            ax.text(v, i, f" {v}", va="center", fontsize=7)

        # multi-omic depth: RNA counts vs ATAC accessible peaks, colored by type
        ax = fig.add_axes([0.38, 0.08, 0.28, 0.36]); ax.set_facecolor("white")
        if atac_acc is not None:
            for c in cats:
                m = (obs["cell_type"] == c).to_numpy()
                ax.scatter(rna_counts[m], atac_acc[m], s=6, color=color[c], linewidths=0)
            ax.set_xscale("log"); ax.set_yscale("log")
            ax.set_xlabel("RNA UMIs / nucleus"); ax.set_ylabel("ATAC accessible peaks / nucleus")
            ax.set_title("Paired RNA + ATAC depth")

        # morphology slot placeholder
        ax = fig.add_axes([0.72, 0.08, 0.24, 0.36]); ax.set_facecolor("white")
        ax.axis("off")
        ax.add_patch(plt.Rectangle((0.05, 0.05), 0.9, 0.9, fill=False, ls="--",
                                   ec="0.6", transform=ax.transAxes))
        ax.text(0.5, 0.62, "Morphology slot", ha="center", fontsize=11, weight="bold",
                transform=ax.transAxes)
        ax.text(0.5, 0.45,
                "awaiting a co-registered\nnuclear image.\nEach nucleus already has an\n(x, y) — the join key.",
                ha="center", fontsize=9, color="0.35", transform=ax.transAxes)

        fig.suptitle(
            "Tessera — SCP2176 slide-tags melanoma: 2,535 nuclei, RNA + ATAC + spatial",
            fontsize=14, y=0.99,
        )
        pdf.savefig(fig, facecolor="white")
        plt.close(fig)


if __name__ == "__main__":
    main()
