"""Load *your own* data into a Tessera :class:`~tessera.core.PairedData`.

Two entry points:

- ``from_anndata`` — you already have an AnnData (the scanpy/scverse case). It just
  needs an expression matrix and per-nucleus spatial coordinates; cell types, a
  UMAP and an ATAC modality are used if present.
- ``build_paired`` — you have plain arrays/matrices (counts, coordinates, labels)
  and want them assembled without touching AnnData directly.

The only hard requirements are an omics matrix (one row per nucleus) and a
per-nucleus spatial (x, y). Everything else is optional and slots in when present.
"""
from __future__ import annotations

from typing import Optional, Sequence

import numpy as np
import pandas as pd
import scipy.sparse as sp

from ..core import PairedData


def from_anndata(
    adata,
    spatial_key: str = "spatial",
    umap_key: str = "X_umap",
    celltype_key: str = "cell_type",
    atac=None,
) -> PairedData:
    """Wrap an existing AnnData. Requires ``obsm[spatial_key]`` (n_obs x 2).

    ``atac`` may be a second AnnData (peaks in var) sharing obs order, and is
    attached as the 'atac' modality.
    """
    if spatial_key not in adata.obsm:
        raise ValueError(
            f"adata.obsm['{spatial_key}'] missing — spatial coordinates are required"
        )
    if celltype_key in adata.obs and not str(adata.obs[celltype_key].dtype).startswith("category"):
        adata.obs[celltype_key] = adata.obs[celltype_key].astype("category")
    if celltype_key != "cell_type" and celltype_key in adata.obs:
        adata.obs["cell_type"] = adata.obs[celltype_key]
    adata.obsm["spatial"] = np.asarray(adata.obsm[spatial_key], dtype=float)[:, :2]
    if umap_key in adata.obsm and umap_key != "X_umap":
        adata.obsm["X_umap"] = np.asarray(adata.obsm[umap_key], dtype=float)[:, :2]
    adata.uns.setdefault(PairedData.MORPH_KEY, [])
    mods = {}
    if atac is not None:
        if list(atac.obs_names) != list(adata.obs_names):
            atac = atac[adata.obs_names].copy()
        mods["atac"] = atac
    return PairedData(adata, mods=mods)


def build_paired(
    rna,
    var_names: Sequence[str],
    cell_ids: Sequence[str],
    spatial,
    cell_type: Optional[Sequence] = None,
    umap=None,
    atac=None,
    atac_peaks: Optional[Sequence[str]] = None,
    obs: Optional[pd.DataFrame] = None,
) -> PairedData:
    """Assemble a PairedData from arrays.

    Parameters
    ----------
    rna : (n_cells, n_genes) array or sparse matrix of counts.
    var_names : gene names (length n_genes).
    cell_ids : per-nucleus ids (length n_cells).
    spatial : (n_cells, 2) tissue coordinates — the join key to imaging.
    cell_type : optional per-nucleus labels.
    umap : optional (n_cells, 2) embedding.
    atac : optional (n_cells, n_peaks) matrix; names via ``atac_peaks``.
    obs : optional extra per-nucleus columns (indexed like cell_ids or aligned).
    """
    import anndata as ad

    ids = [str(c) for c in cell_ids]
    o = pd.DataFrame(index=pd.Index(ids, name="cell_id"))
    if cell_type is not None:
        o["cell_type"] = pd.Categorical([str(c) for c in cell_type])
    if obs is not None:
        obs = obs.reset_index(drop=True) if len(obs) == len(ids) else obs.reindex(ids)
        for col in obs.columns:
            o[col] = np.asarray(obs[col])

    adata = ad.AnnData(X=sp.csr_matrix(rna), obs=o, var=pd.DataFrame(index=list(var_names)))
    adata.var_names_make_unique()
    sp_arr = np.asarray(spatial, dtype=float)[:, :2]
    if sp_arr.shape[0] != adata.n_obs:
        raise ValueError(f"spatial has {sp_arr.shape[0]} rows, expected {adata.n_obs}")
    adata.obsm["spatial"] = sp_arr
    if umap is not None:
        adata.obsm["X_umap"] = np.asarray(umap, dtype=float)[:, :2]
    adata.uns[PairedData.MORPH_KEY] = []

    mods = {}
    if atac is not None:
        at = ad.AnnData(X=sp.csr_matrix(atac))
        at.obs_names = ids
        at.var_names = list(atac_peaks) if atac_peaks is not None else [
            f"peak_{i}" for i in range(at.n_vars)
        ]
        mods["atac"] = at
    return PairedData(adata, mods=mods)
