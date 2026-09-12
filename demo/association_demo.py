"""Morphology <-> multi-omics association on the paired melanoma bundle.

Correlates nuclear/cell morphology with gene expression + ATAC, and shows how
strongly each morphology feature separates transcriptomic cell types.

    python demo/association_demo.py
"""
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages

sys.path.insert(0, "src")
from tessera.analysis import (  # noqa: E402
    morphology_celltype_effect, morphology_feature_correlation, morphology_feature_gene,
)
from tessera.core import PairedData  # noqa: E402

BUNDLE = "data/processed/melanoma_7aad.tessera"
OUT = "data/processed/tessera_association.pdf"


def top_bar(ax, df, k, title, ink, muted):
    pos = df.head(k)[::-1]
    neg = df.tail(k)
    d = list(neg.itertuples()) + list(pos.itertuples())
    names = [r.name for r in d]
    rs = [r.r for r in d]
    colors = ["#C2543A" if v < 0 else "#2a8f86" for v in rs]
    ax.barh(range(len(rs)), rs, color=colors)
    ax.set_yticks(range(len(names))); ax.set_yticklabels(names, fontsize=7.5, family="monospace")
    ax.axvline(0, color=muted, lw=.8)
    ax.set_xlabel("Pearson r (log1p-CPM)", fontsize=8, color=muted)
    ax.set_title(title, fontsize=10.5, color=ink, weight="bold")
    ax.tick_params(colors=muted, labelsize=7.5)
    for s in ("top", "right"): ax.spines[s].set_visible(False)


def main():
    p = PairedData.load(BUNDLE)
    ink, muted = "#14242A", "#697B81"

    area_rna = morphology_feature_gene(p, "area", layer="rna", min_cells=120)
    nc_rna = morphology_feature_gene(p, "nc_ratio", layer="rna", min_cells=120)
    tex_rna = morphology_feature_gene(p, "intensity_std", layer="rna", min_cells=120)
    area_atac = morphology_feature_gene(p, "area", layer="atac", min_cells=200)

    print("Top genes + correlated with nuclear AREA:")
    print(area_rna.head(8).to_string(index=False))
    print("\nTop genes + correlated with N:C ratio:")
    print(nc_rna.head(8).to_string(index=False))
    print("\nTop ATAC peaks + correlated with nuclear AREA:")
    print(area_atac.head(5).to_string(index=False))

    feats = ["area", "circularity", "integrated_intensity", "intensity_std",
             "eccentricity", "solidity", "cell_area", "nc_ratio"]
    feats = [f for f in feats if f in p.adata.obs]
    effects = [morphology_celltype_effect(p, f) for f in feats]
    effects = sorted(effects, key=lambda e: e["eta2"], reverse=True)
    print("\nCell-type separation (eta^2):")
    for e in effects:
        print(f"  {e['feature']:22s} eta2={e['eta2']:.3f}  p={e['p']:.1e}")

    with PdfPages(OUT) as pdf:
        fig = plt.figure(figsize=(15, 9.5), facecolor="white")
        fig.suptitle("Tessera — morphology ↔ multi-omics association (melanoma slide-tags)",
                     fontsize=14, y=0.98, color=ink)

        ax = fig.add_axes([0.05, 0.56, 0.26, 0.34]); ax.set_facecolor("white")
        top_bar(ax, area_rna, 10, "Genes vs nuclear AREA (RNA)", ink, muted)
        ax = fig.add_axes([0.38, 0.56, 0.26, 0.34]); ax.set_facecolor("white")
        top_bar(ax, nc_rna, 10, "Genes vs N:C ratio (RNA)", ink, muted)
        ax = fig.add_axes([0.71, 0.56, 0.26, 0.34]); ax.set_facecolor("white")
        top_bar(ax, tex_rna, 10, "Genes vs chromatin texture (RNA)", ink, muted)

        # morphology feature correlation heatmap
        corr = morphology_feature_correlation(p, feats)
        ax = fig.add_axes([0.06, 0.09, 0.34, 0.36]); ax.set_facecolor("white")
        im = ax.imshow(corr.to_numpy(), cmap="RdBu_r", vmin=-1, vmax=1)
        ax.set_xticks(range(len(feats))); ax.set_yticks(range(len(feats)))
        ax.set_xticklabels(feats, rotation=45, ha="right", fontsize=7.5)
        ax.set_yticklabels(feats, fontsize=7.5)
        ax.set_title("Morphology feature correlation", fontsize=10.5, color=ink, weight="bold")
        for (i, j), v in np.ndenumerate(corr.to_numpy()):
            ax.text(j, i, f"{v:.1f}", ha="center", va="center", fontsize=6,
                    color="white" if abs(v) > 0.6 else ink)
        cb = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.03); cb.ax.tick_params(labelsize=7)

        # eta2 barh
        ax = fig.add_axes([0.50, 0.09, 0.30, 0.36]); ax.set_facecolor("white")
        names = [e["feature"] for e in effects][::-1]
        vals = [e["eta2"] for e in effects][::-1]
        ax.barh(range(len(vals)), vals, color="#0EA5B7")
        ax.set_yticks(range(len(names))); ax.set_yticklabels(names, fontsize=8)
        ax.set_xlabel("eta² (cell-type separation)", fontsize=8, color=muted)
        ax.set_title("Which morphology features track cell identity", fontsize=10.5,
                     color=ink, weight="bold")
        ax.tick_params(colors=muted, labelsize=8)
        for s in ("top", "right"): ax.spines[s].set_visible(False)

        fig.text(0.83, 0.44,
                 "Note: in this test the\nmorphology is synthesised\nfrom transcriptomic cell\n"
                 "type, so these associations\nrecover cell-type marker\nprograms (e.g. melanoma\n"
                 "genes ↑ with nuclear size,\nT-cell genes ↓). On real\nimaging the same call\n"
                 "surfaces genuine\nmorphology–expression\nlinks.",
                 fontsize=8.5, color=muted, va="top")

        pdf.savefig(fig, facecolor="white"); plt.close(fig)
    print("\nwrote", OUT)


if __name__ == "__main__":
    main()
