"""Export a compact web-viewer data payload from a Tessera bundle and inject it
into the HTML template to produce a self-contained index.html.

    python demo/export_web.py <bundle> <template.html> <out.html>
"""
import json
import sys

import numpy as np

sys.path.insert(0, "src")
from tessera.core import PairedData  # noqa: E402

MARKERS = [
    # melanoma
    "MLANA", "PMEL", "TYR", "MITF", "DCT", "TYRP1", "SOX10",
    # T / NK
    "CD3D", "CD3E", "CD8A", "CD4", "IL7R", "FOXP3", "CTLA4", "GZMB", "NKG7", "GNLY",
    # B / plasma
    "MS4A1", "MZB1", "IGHG1", "JCHAIN", "CD79A",
    # myeloid / DC
    "LYZ", "CD68", "CD14", "ITGAX", "C1QA", "FCGR3A", "CLEC9A", "LILRA4", "CLEC4C",
    "PTPRC",
]


def main(bundle, template, out):
    p = PairedData.load(bundle)
    a = p.adata
    obs = a.obs
    cats = list(obs["cell_type"].cat.categories)
    ct = obs["cell_type"].cat.codes.to_numpy().astype(int)

    xy = a.obsm["spatial"]
    um = a.obsm["X_umap"]

    X = a.X.tocsr()
    umis = np.asarray(X.sum(1)).ravel().astype(int)
    genes = np.asarray((X > 0).sum(1)).ravel().astype(int)
    atac = p.mods["atac"].X.tocsr()
    npk = np.asarray((atac > 0).sum(1)).ravel().astype(int)

    var_names = list(a.var_names)
    var_idx = {g: i for i, g in enumerate(var_names)}

    # top 6 genes per nucleus (from sparse rows)
    top = []
    for i in range(a.n_obs):
        s, e = X.indptr[i], X.indptr[i + 1]
        idx, dat = X.indices[s:e], X.data[s:e]
        if len(dat) == 0:
            top.append([])
            continue
        k = min(6, len(dat))
        sub = np.argpartition(dat, -k)[-k:]
        sub = sub[np.argsort(dat[sub])[::-1]]
        top.append([[var_names[int(idx[j])], int(dat[j])] for j in sub])

    # marker gene per-nucleus counts (present markers only)
    markers = {}
    for g in MARKERS:
        if g in var_idx:
            col = X[:, var_idx[g]]
            markers[g] = np.asarray(col.todense()).ravel().astype(int).tolist()

    def tcr(col):
        if col not in obs:
            return None
        out = []
        for v in obs[col]:
            sv = str(v)
            out.append(None if sv in ("nan", "NA", "None", "") else sv)
        return out

    data = {
        "categories": cats,
        "n": int(a.n_obs),
        "ids": list(a.obs_names),
        "sx": [round(float(v), 1) for v in xy[:, 0]],
        "sy": [round(float(v), 1) for v in xy[:, 1]],
        "ux": [round(float(v), 3) for v in um[:, 0]],
        "uy": [round(float(v), 3) for v in um[:, 1]],
        "ct": ct.tolist(),
        "umis": umis.tolist(),
        "genes": genes.tolist(),
        "npk": npk.tolist(),
        "top": top,
        "markers": markers,
        "tcr_a": tcr("tcr_alpha"),
        "tcr_b": tcr("tcr_beta"),
    }

    payload = json.dumps(data, separators=(",", ":"))
    html = open(template).read().replace("__TESSERA_DATA__", payload)
    open(out, "w").write(html)
    mb = len(html.encode()) / 1e6
    print(f"markers embedded: {len(markers)} ({', '.join(markers)})")
    print(f"wrote {out}  ({mb:.2f} MB)")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3])
