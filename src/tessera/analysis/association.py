"""Associate per-nucleus morphology with the paired multi-omics.

The payoff of pairing: once every nucleus has both morphology and an RNA/ATAC
profile, ask which molecular programs track a morphological axis (nuclear size,
N:C ratio, chromatin texture…).

- ``morphology_feature_gene`` — correlate one morphology feature with every
  (sufficiently detected) gene or ATAC peak; returns a signed, sorted table.
- ``morphology_feature_correlation`` — correlation matrix among morphology
  features themselves.
- ``morphology_celltype_effect`` — how strongly a feature separates cell types
  (Kruskal–Wallis H, p, and eta^2 effect size).

Correlation is Pearson on log1p-CPM-normalised expression by default (the single-
cell standard for feature–gene correlation); pass ``method='spearman'`` for rank
correlation.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from ..core import PairedData


def _matrix(paired: PairedData, layer: str):
    if layer == "rna":
        return paired.adata.X, list(paired.adata.var_names)
    m = paired.mods[layer]
    return m.X, list(m.var_names)


def _pearson_vec(x: np.ndarray, Y: np.ndarray) -> np.ndarray:
    """Correlation of vector x (n,) with each column of Y (n, g)."""
    xc = x - x.mean()
    yc = Y - Y.mean(0)
    den = np.sqrt((yc ** 2).sum(0)) * np.sqrt((xc ** 2).sum())
    with np.errstate(divide="ignore", invalid="ignore"):
        r = (yc.T @ xc) / den
    return r


def morphology_feature_gene(
    paired: PairedData,
    feature: str,
    layer: str = "rna",
    method: str = "pearson",
    min_cells: int = 100,
    normalize: bool = True,
) -> pd.DataFrame:
    """Signed correlation of ``feature`` with each detected feature in ``layer``.

    Returns a DataFrame [name, r, n_detected] sorted by r (positive first).
    """
    a = paired.adata
    if feature not in a.obs:
        raise KeyError(f"{feature!r} not in obs; run imaging first")
    y = a.obs[feature].to_numpy(dtype=float)
    mask = ~np.isnan(y)
    yv = y[mask]

    X, names = _matrix(paired, layer)
    X = X.tocsr()[mask]
    det = np.asarray((X > 0).sum(0)).ravel()
    keep = np.where(det >= min_cells)[0]
    if keep.size == 0:
        raise ValueError("no features pass min_cells; lower it")
    Xk = X[:, keep].toarray().astype(float)
    names = [names[i] for i in keep]

    if normalize:
        rs = Xk.sum(1, keepdims=True)
        rs[rs == 0] = 1.0
        Xk = np.log1p(Xk / rs * 1e4)

    if method == "spearman":
        xv = _rank(yv)
        Xk = np.apply_along_axis(_rank, 0, Xk)
    else:
        xv = yv
    r = _pearson_vec(xv, Xk)

    out = pd.DataFrame({"name": names, "r": r, "n_detected": det[keep]}).dropna(subset=["r"])
    return out.sort_values("r", ascending=False).reset_index(drop=True)


def _rank(v: np.ndarray) -> np.ndarray:
    from scipy.stats import rankdata
    return rankdata(v)


def morphology_feature_correlation(paired: PairedData, features=None) -> pd.DataFrame:
    a = paired.adata
    if features is None:
        features = list(a.uns.get(PairedData.MORPH_KEY, [])) + list(a.uns.get("cell_features", []))
    features = [f for f in features if f in a.obs]
    return a.obs[features].astype(float).corr()


def morphology_celltype_effect(paired: PairedData, feature: str, groupby: str = "cell_type") -> dict:
    from scipy.stats import kruskal

    a = paired.adata
    df = a.obs[[feature, groupby]].dropna()
    groups = [g[feature].to_numpy(float) for _, g in df.groupby(groupby, observed=True) if len(g) > 1]
    if len(groups) < 2:
        return {"feature": feature, "H": np.nan, "p": np.nan, "eta2": np.nan}
    H, p = kruskal(*groups)
    k, N = len(groups), len(df)
    eta2 = (H - k + 1) / (N - k) if N > k else np.nan
    return {"feature": feature, "H": float(H), "p": float(p), "eta2": float(eta2), "n_groups": k}
