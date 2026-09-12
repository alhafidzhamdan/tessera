"""Ingest a slide-tags multiome study (SCP2176 human melanoma layout) into a
:class:`~tessera.core.PairedData`.

Expected directory layout (as downloaded from the Broad Single Cell Portal)::

    root/
      cluster/HumanMelanomaMultiome_cluster.csv    # UMAP X/Y + cell_type
      cluster/HumanMelanomaMultiome_spatial.csv    # tissue X/Y + cell_type
      metadata/HumanMelanomaMultiome_metadata.csv  # donor/disease/organ...
      expression/<hash>/{matrix.mtx.gz,features.tsv.gz,barcodes.tsv.gz}  # RNA
      expression/HumanMelanomaMultiome_atac.csv.gz # peaks x cells
      other/slidetags_multiome_tcr.csv             # TCR alpha/beta

The nucleus set is defined by the spatial file (the analysed nuclei). RNA/ATAC
are subset and reordered to match, so every modality shares nuclei order/ids.
The tissue coordinate becomes ``obsm['spatial']`` -- the join key to a future
registered nuclear image.
"""
from __future__ import annotations

import glob
import gzip
import os
from typing import Optional

import numpy as np
import pandas as pd
import scipy.io
import scipy.sparse as sp

from ..core import PairedData


# --------------------------------------------------------------------------- #
# helpers
# --------------------------------------------------------------------------- #
def _read_scp_csv(path: str) -> pd.DataFrame:
    """Read a Single Cell Portal CSV, dropping the second 'TYPE' header row.
    Indexed by the NAME column (the cell barcode)."""
    df = pd.read_csv(path)
    if len(df) and str(df.iloc[0, 0]).upper() == "TYPE":
        df = df.iloc[1:].reset_index(drop=True)
    df = df.set_index(df.columns[0])
    df.index.name = "cell_id"
    return df


def _find_rna_dir(root: str) -> str:
    for cand in glob.glob(os.path.join(root, "expression", "*")):
        if os.path.isdir(cand) and os.path.exists(os.path.join(cand, "matrix.mtx.gz")):
            return cand
    raise FileNotFoundError("RNA 10x triplet (matrix.mtx.gz) not found under expression/")


def _load_rna(rna_dir: str) -> "ad.AnnData":  # noqa: F821
    import anndata as ad

    with gzip.open(os.path.join(rna_dir, "matrix.mtx.gz"), "rb") as fh:
        m = scipy.io.mmread(fh)  # genes x cells
    x = sp.csr_matrix(m.T)  # cells x genes
    feats = pd.read_csv(
        os.path.join(rna_dir, "features.tsv.gz"), sep="\t", header=None,
        names=["gene_id", "symbol", "feature_type", "chrom", "start", "end"],
    )
    barcodes = pd.read_csv(
        os.path.join(rna_dir, "barcodes.tsv.gz"), header=None
    )[0].tolist()
    var = feats.set_index("symbol")
    adata = ad.AnnData(X=x, var=var)
    adata.obs_names = barcodes
    adata.var["gene_ids"] = feats["gene_id"].to_numpy()
    adata.var_names_make_unique()
    return adata


def _load_atac(path: str) -> "ad.AnnData":  # noqa: F821
    import anndata as ad

    with gzip.open(path, "rb") as fh:
        df = pd.read_csv(fh, index_col=0)  # peaks x cells (integer counts)
    peaks = df.index.astype(str)
    cells = df.columns.astype(str)
    # peaks x cells -> sparse -> transpose to cells x peaks (float32)
    x = sp.csr_matrix(df.to_numpy(dtype=np.float32)).T.tocsr()
    del df
    adata = ad.AnnData(X=x)
    adata.obs_names = cells
    adata.var_names = peaks
    return adata


# --------------------------------------------------------------------------- #
# main
# --------------------------------------------------------------------------- #
def load_slidetags(
    root: str,
    prefix: str = "HumanMelanomaMultiome",
    load_atac: bool = True,
) -> PairedData:
    import anndata as ad

    cl_dir = os.path.join(root, "cluster")
    spatial = _read_scp_csv(os.path.join(cl_dir, f"{prefix}_spatial.csv"))
    umap = _read_scp_csv(os.path.join(cl_dir, f"{prefix}_cluster.csv"))
    meta = _read_scp_csv(os.path.join(root, "metadata", f"{prefix}_metadata.csv"))

    cells = list(spatial.index)  # the analysed nuclei define order

    # --- RNA (primary) ----------------------------------------------------- #
    rna = _load_rna(_find_rna_dir(root))
    missing = [c for c in cells if c not in set(rna.obs_names)]
    if missing:
        raise ValueError(f"{len(missing)} spatial nuclei missing from RNA matrix")
    rna = rna[cells].copy()  # subset + reorder to the analysed nuclei

    # --- obs: cell_type, spatial, umap, metadata, TCR ---------------------- #
    obs = pd.DataFrame(index=pd.Index(cells, name="cell_id"))
    obs["cell_type"] = spatial["cell_type"].astype("category")
    for col in meta.columns:
        if col != "cluster":
            obs[col] = meta.reindex(cells)[col].to_numpy()
    tcr_path = os.path.join(root, "other", "slidetags_multiome_tcr.csv")
    if os.path.exists(tcr_path):
        tcr = pd.read_csv(tcr_path)
        if "CB" in tcr.columns:
            tcr = tcr.set_index("CB")
            for col in ("alpha", "beta"):
                if col in tcr.columns:
                    obs[f"tcr_{col}"] = tcr.reindex(cells)[col].to_numpy()
    rna.obs = obs

    rna.obsm["spatial"] = spatial.loc[cells, ["X", "Y"]].to_numpy(dtype=float)
    rna.obsm["X_umap"] = umap.reindex(cells)[["X", "Y"]].to_numpy(dtype=float)
    rna.uns["study"] = {"accession": "SCP2176", "assay": "slide-tags multiome",
                        "disease": "melanoma", "organ": "lymph node"}
    rna.uns[PairedData.MORPH_KEY] = []  # no imaging yet

    mods = {}
    if load_atac:
        atac = _load_atac(os.path.join(root, "expression", f"{prefix}_atac.csv.gz"))
        atac_missing = [c for c in cells if c not in set(atac.obs_names)]
        if atac_missing:
            raise ValueError(f"{len(atac_missing)} nuclei missing from ATAC matrix")
        atac = atac[cells].copy()
        atac.obs["cell_type"] = obs["cell_type"].to_numpy()
        mods["atac"] = atac

    return PairedData(rna, crops=None, mods=mods)
