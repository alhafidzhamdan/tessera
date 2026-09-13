"""Export a compact web-viewer data payload from a Tessera bundle and inject it
into the HTML template to produce a self-contained index.html.

    python demo/export_web.py <bundle> <template.html> <out.html>
"""
import json
import os
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


def write_expression_files(out_dir, X, var_names, min_cells=3, chunk_chars=14_000_000):
    """Write the full sparse RNA matrix as companion JS files for gene search.

    Encoding: CSC over genes detected in >= min_cells; concatenated uint16 cell
    indices + uint8 counts (lossless for snRNA counts), plus a uint32 indptr, all
    base64. Big base64 arrays are split into <=chunk_chars pieces so no single file
    exceeds the artifact per-file cap. Returns {published_name: source_path}.
    """
    import base64

    Xc = X.tocsc()
    det = np.asarray((Xc > 0).sum(0)).ravel()
    keep = np.where(det >= min_cells)[0]
    genes = [var_names[i] for i in keep]
    idx_parts, cnt_parts, indptr = [], [], [0]
    for gi in keep:
        s, e = Xc.indptr[gi], Xc.indptr[gi + 1]
        idx_parts.append(Xc.indices[s:e].astype("<u2"))
        cnt_parts.append(np.clip(Xc.data[s:e], 0, 255).astype(np.uint8))
        indptr.append(indptr[-1] + int(e - s))
    idx = np.concatenate(idx_parts) if idx_parts else np.array([], "<u2")
    cnt = np.concatenate(cnt_parts) if cnt_parts else np.array([], np.uint8)
    indptr = np.asarray(indptr, dtype="<u4")

    def shard(b64, prefix, arrname):
        names = []
        for j in range(0, len(b64), chunk_chars):
            nm = f"{prefix}_{j // chunk_chars}.js"
            path = os.path.join(out_dir, nm)
            with open(path, "w") as fh:
                fh.write(f'window.EXPR.{arrname}.push("{b64[j:j+chunk_chars]}");')
            names.append(nm)
        return names

    idx_files = shard(base64.b64encode(idx.tobytes()).decode(), "expr_idx", "idxChunks")
    cnt_files = shard(base64.b64encode(cnt.tobytes()).decode(), "expr_cnt", "cntChunks")

    meta = {
        "n_cells": int(X.shape[0]), "genes": genes,
        "indptr_b64": base64.b64encode(indptr.tobytes()).decode(),
        "idx_files": idx_files, "cnt_files": cnt_files,
    }
    files = {"expr_meta.js": os.path.join(out_dir, "expr_meta.js")}
    with open(files["expr_meta.js"], "w") as fh:
        fh.write("window.EXPR=" + json.dumps(meta, separators=(",", ":")) +
                 ";window.EXPR.idxChunks=[];window.EXPR.cntChunks=[];")
    for nm in idx_files + cnt_files:
        files[nm] = os.path.join(out_dir, nm)
    total = sum(os.path.getsize(p) for p in files.values()) / 1e6
    print(f"expression companion: {len(genes)} genes, {idx.size} nonzeros, "
          f"{len(files)} files, {total:.1f} MB")
    return files


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

    # DE gene panel: curated markers + top-dispersion HVGs, raw counts per cell,
    # so the browser can run region-vs-rest differential expression on real genes.
    Xc = X.tocsc()
    nS = X.shape[0]
    csum = np.asarray(Xc.sum(0)).ravel()
    csq = np.asarray(Xc.multiply(Xc).sum(0)).ravel()
    gmean = csum / nS
    gvar = csq / nS - gmean ** 2
    gdet = np.asarray((Xc > 0).sum(0)).ravel()
    disp = np.divide(gvar, gmean, out=np.zeros_like(gvar), where=gmean > 0)
    elig = np.where(gdet >= 100)[0]
    hvg = list(elig[np.argsort(disp[elig])[::-1]][:160])
    mk = [var_idx[g] for g in MARKERS if g in var_idx]
    panel = list(dict.fromkeys(hvg + mk))
    deg = {}
    for gi in panel:
        deg[var_names[gi]] = np.asarray(Xc[:, gi].todense()).ravel().astype(int).tolist()
    print(f"DE panel: {len(deg)} genes")

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

    # morphology columns: schema + values (empty until a nuclear image is paired)
    morph_features = p.morphology_features
    morph = {}
    for f in morph_features:
        if f in obs:
            morph[f] = [
                None if (v != v) else round(float(v), 2) for v in obs[f].to_numpy()
            ]
    data["morph_features"] = morph_features
    data["morph"] = morph
    data["stain"] = a.uns.get("stain")
    data["morph_source"] = a.uns.get("morphology_source")

    # whole-cell morphology + N:C ratio
    cell_features = list(a.uns.get("cell_features", []))
    cellmorph = {}
    for f in cell_features:
        if f in obs:
            cellmorph[f] = [
                None if (v != v) else round(float(v), 3) for v in obs[f].to_numpy()
            ]
    data["cell_features"] = cell_features
    data["cellmorph"] = cellmorph
    data["deg_genes"] = [var_names[gi] for gi in panel]
    data["deg"] = deg

    # embed per-nucleus crops as one RGB PNG atlas: nucleus -> red, membrane -> green
    if p.crops is not None:
        import base64
        import io

        from PIL import Image

        cs = p.crops.crop_shape[0]
        cols = 51
        rows = (a.n_obs + cols - 1) // cols
        atlas = np.zeros((rows * cs, cols * cs, 3), np.uint8)
        for i, cid in enumerate(a.obs_names):
            cr = p.crops.get(cid)  # (channels, H, W)
            r, c = divmod(i, cols)
            tile = atlas[r * cs:(r + 1) * cs, c * cs:(c + 1) * cs]
            tile[..., 0] = cr[0]
            if cr.shape[0] > 1:
                tile[..., 1] = cr[1]  # membrane -> green
            else:
                tile[..., 1] = cr[0]
                tile[..., 2] = cr[0]  # grayscale
        buf = io.BytesIO()
        Image.fromarray(atlas, "RGB").save(buf, format="PNG", optimize=True)
        b64 = base64.b64encode(buf.getvalue()).decode()
        data["atlas"] = "data:image/png;base64," + b64
        data["atlas_cols"] = cols
        data["crop_px"] = cs
        data["crop_channels"] = list(p.crops.channel_names)
        print(f"crop atlas: {rows}x{cols} @ {cs}px RGB, {len(b64)/1e6:.2f} MB base64")

    payload = json.dumps(data, separators=(",", ":"))
    html = open(template).read().replace("__TESSERA_DATA__", payload)
    open(out, "w").write(html)
    mb = len(html.encode()) / 1e6
    print(f"markers embedded: {len(markers)} ({', '.join(markers)})")
    print(f"wrote {out}  ({mb:.2f} MB)")

    files = write_expression_files(os.path.dirname(os.path.abspath(out)), X, var_names)
    print("COMPANION_FILES " + " ".join(sorted(files)))


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3])
