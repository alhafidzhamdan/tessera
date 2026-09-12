"""Synthetic slide-tags-like paired dataset.

Generates a :class:`~tessera.core.PairedData` that mirrors the SCP2176 human
melanoma multiome study: the same cell types and rough proportions, spatial
tissue coordinates, an RNA matrix, an ATAC gene-activity layer, a UMAP, and --
crucially -- a synthetic DAPI-like tissue image from which real morphology
features and per-nucleus crops are extracted. This lets the container, catalogue
and viewer be built and tested before any real imaging exists.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from skimage.draw import ellipse
from skimage.filters import gaussian

from ..core import CropStore, PairedData
from .morphology import MORPHOLOGY_FEATURES, cut_crops, extract_morphology

# (name, fraction, major_axis_mean, major_axis_sd, eccentricity, dapi_intensity)
CELL_TYPES = [
    ("tumour_1", 0.221, 20.0, 4.0, 0.55, 150),
    ("tumour_2", 0.107, 18.0, 3.5, 0.60, 140),
    ("T_CD8", 0.402, 9.0, 1.0, 0.30, 200),
    ("T_CD4", 0.026, 9.0, 1.0, 0.30, 200),
    ("T_reg", 0.028, 9.0, 1.0, 0.32, 195),
    ("mono-mac", 0.121, 13.0, 2.0, 0.45, 170),
    ("myeloid", 0.014, 13.0, 2.5, 0.50, 165),
    ("plasma", 0.065, 12.0, 1.5, 0.55, 185),
    ("mDC", 0.011, 11.0, 1.5, 0.48, 175),
    ("pDC", 0.005, 11.0, 1.5, 0.48, 175),
]


def make_synthetic(
    n_cells: int = 2535,
    n_genes: int = 200,
    canvas: int = 1400,
    crop_size: int = 64,
    crops_path: str = "synthetic_crops.zarr",
    seed: int = 0,
) -> PairedData:
    rng = np.random.default_rng(seed)
    names = [c[0] for c in CELL_TYPES]
    fracs = np.array([c[1] for c in CELL_TYPES])
    fracs = fracs / fracs.sum()

    # -- assign cell types -------------------------------------------------- #
    counts = np.round(fracs * n_cells).astype(int)
    counts[-1] += n_cells - counts.sum()  # fix rounding
    ct = np.concatenate([[names[i]] * counts[i] for i in range(len(names))])
    rng.shuffle(ct)
    n = len(ct)
    ct = pd.Categorical(ct, categories=names)
    type_idx = {name: i for i, name in enumerate(names)}
    ti = np.array([type_idx[c] for c in ct])

    # -- spatial layout: tumour nests + immune infiltrate ------------------- #
    margin = crop_size
    nests = rng.uniform(margin + 100, canvas - margin - 100, size=(3, 2))
    xy = np.zeros((n, 2))
    is_tumour = np.array([c.startswith("tumour") or c == "plasma" for c in ct])
    for i in range(n):
        if is_tumour[i]:
            nest = nests[rng.integers(0, 3)]
            xy[i] = rng.normal(nest, 130)
        else:  # immune: scattered across tissue, mild nest avoidance
            xy[i] = rng.uniform(margin, canvas - margin, size=2)
    xy = np.clip(xy, margin, canvas - margin)

    # -- render DAPI-like tissue image + label image ------------------------ #
    img = np.zeros((canvas, canvas), dtype=np.float32)
    labels = np.zeros((canvas, canvas), dtype=np.int32)
    majors, eccs = np.zeros(n), np.zeros(n)
    for i in range(n):
        _, _, maj_m, maj_sd, ecc, inten = CELL_TYPES[ti[i]]
        major = max(4.0, rng.normal(maj_m, maj_sd)) / 2.0
        minor = major * np.sqrt(max(0.05, 1 - ecc**2))
        majors[i], eccs[i] = major, ecc
        rot = rng.uniform(0, np.pi)
        rr, cc = ellipse(
            xy[i, 1], xy[i, 0], major, minor, shape=img.shape, rotation=rot
        )
        val = inten + rng.normal(0, 15, size=rr.shape)
        img[rr, cc] = np.clip(val, 0, 255)
        labels[rr, cc] = i + 1  # later cells overwrite on overlap (like seg)

    img = gaussian(img, sigma=1.2, preserve_range=True)
    img = np.clip(img, 0, 255).astype(np.uint8)

    # -- morphology features (real extraction path) ------------------------- #
    morph = extract_morphology(labels, intensity_image=img)
    feat_cols = [f for f in MORPHOLOGY_FEATURES if f in morph.columns]
    feat_cols += [c for c in ("intensity_mean", "intensity_std") if c in morph.columns]
    morph_aligned = morph.reindex(np.arange(1, n + 1))  # label i+1 -> cell i
    morph_aligned.index = np.arange(n)

    # crops centred on the intended positions (robust to overlap merges)
    crops = cut_crops(img, centroids=list(zip(xy[:, 1], xy[:, 0])), size=crop_size)

    # -- RNA + ATAC --------------------------------------------------------- #
    gene_names = [f"gene_{i:03d}" for i in range(n_genes)]
    base = rng.gamma(1.5, 1.0, size=n_genes)  # housekeeping baseline
    programs = np.tile(base, (len(names), 1))
    for t in range(len(names)):  # each type boosts ~12 marker genes
        markers = rng.choice(n_genes, size=12, replace=False)
        programs[t, markers] += rng.uniform(4, 9, size=markers.size)
    lib = rng.lognormal(0.0, 0.3, size=n)  # library-size variation
    rna_mean = programs[ti] * lib[:, None]
    rna = rng.poisson(rna_mean).astype(np.float32)
    # ATAC gene-activity: correlated with RNA program + independent noise
    atac_mean = 0.6 * rna_mean + rng.gamma(1.0, 0.8, size=(n, n_genes))
    atac = rng.poisson(atac_mean).astype(np.float32)

    # -- UMAP-ish embedding via PCA(2) on log-normalised RNA ---------------- #
    x = np.log1p(rna / (rna.sum(1, keepdims=True) + 1e-9) * 1e4)
    x = x - x.mean(0)
    u, s, vt = np.linalg.svd(x, full_matrices=False)
    umap = (u[:, :2] * s[:2]) + rng.normal(0, 0.15, size=(n, 2))

    # -- assemble AnnData + crop store -------------------------------------- #
    import anndata as ad

    cell_ids = [f"cell_{i:05d}" for i in range(n)]
    obs = pd.DataFrame(index=cell_ids)
    obs["cell_type"] = ct.astype(str)
    obs["cell_type"] = obs["cell_type"].astype("category")
    for col in feat_cols:
        obs[col] = morph_aligned[col].to_numpy()
    obs["centroid_x"] = xy[:, 0]
    obs["centroid_y"] = xy[:, 1]

    adata = ad.AnnData(
        X=rna,
        obs=obs,
        var=pd.DataFrame(index=gene_names),
        layers={"atac": atac},
    )
    adata.obsm["spatial"] = xy
    adata.obsm["X_umap"] = umap
    adata.uns[PairedData.MORPH_KEY] = feat_cols
    adata.uns["tissue_image_shape"] = list(img.shape)

    store = CropStore.create(
        crops_path, cell_ids=cell_ids, shape=(crop_size, crop_size),
        channels=1, channel_names=["DAPI"], dtype="uint8",
    )
    store._arr[:] = crops
    return PairedData(adata, store)
