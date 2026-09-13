"""Write a test .h5ad from the melanoma bundle to develop/verify the browser
HDF5 reader. Subsets to a manageable gene panel (incl. SOX2 + markers)."""
import sys

import numpy as np

sys.path.insert(0, "src")
from tessera.core import PairedData  # noqa: E402

p = PairedData.load("data/processed/melanoma_7aad.tessera")
a = p.adata
X = a.X.tocsc()
det = np.asarray((X > 0).sum(0)).ravel()
keep = np.where(det >= 30)[0]
# make sure SOX2 + a few markers are present
for g in ["SOX2", "MLANA", "CD8A", "SOX10", "LYZ", "MS4A1"]:
    if g in a.var_names:
        gi = list(a.var_names).index(g)
        if gi not in keep:
            keep = np.append(keep, gi)
keep = np.unique(keep)
sub = a[:, keep].copy()
sub.X = sub.X.tocsr()
# keep only the light bits (drop big uns) for a clean test file
sub.uns = {}
for k in list(sub.obsm.keys()):
    if k not in ("spatial", "X_umap"):
        del sub.obsm[k]
out = "/private/tmp/claude-501/-Users-alhafidzhamdan-Claude-Projects-github-page-personal/3147b092-afb4-45d7-b201-fc762c67b4f1/scratchpad/test.h5ad"
sub.write_h5ad(out)
print("wrote", out, "shape", sub.shape, "SOX2 in var:", "SOX2" in sub.var_names)
