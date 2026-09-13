/* Parse an AnnData .h5ad (HDF5) in the browser with jsfive (pure JS, no WASM),
   and build the Tessera DATA payload + an in-memory CSC matrix for gene search.
   Requires window.hdf5 (jsfive) to be loaded first. Exposes window.parseH5ad(buf,name)
   -> { DATA, expr }.  expr = { n_cells, genes, ptr:Uint32, idx:Uint32, cnt:Float32 }. */
(function () {
  const MARKERS = ["MLANA","PMEL","TYR","MITF","DCT","TYRP1","SOX10","CD3D","CD3E",
    "CD8A","CD4","IL7R","FOXP3","CTLA4","GZMB","NKG7","GNLY","MS4A1","MZB1","IGHG1",
    "JCHAIN","CD79A","LYZ","CD68","CD14","ITGAX","C1QA","FCGR3A","CLEC9A","LILRA4",
    "CLEC4C","PTPRC"];
  const MORPH_NUC = ["area","perimeter","eccentricity","solidity","extent","orientation",
    "axis_major_length","axis_minor_length","equivalent_diameter_area","circularity",
    "intensity_mean","intensity_std","integrated_intensity"];
  const MORPH_CELL = ["cell_area","cell_eccentricity","cell_solidity",
    "cell_axis_major_length","cell_circularity","nc_ratio"];

  const isGroup = (x) => x && typeof x.get === "function" && Array.isArray(x.keys);
  const toNum = (arr) => {
    if (arr instanceof BigInt64Array || arr instanceof BigUint64Array) {
      const o = new Float64Array(arr.length);
      for (let i = 0; i < arr.length; i++) o[i] = Number(arr[i]);
      return o;
    }
    return arr;
  };
  function readStr(x) {
    if (!x) return [];
    if (isGroup(x)) {
      if (x.keys.includes("values")) return Array.from(x.get("values").value);
      if (x.keys.includes("categories")) {  // categorical -> decode
        const cats = Array.from(x.get("categories").value);
        const codes = toNum(x.get("codes").value);
        return Array.from(codes, (c) => (c < 0 ? null : cats[c]));
      }
    }
    return Array.from(x.value);
  }
  function col2(ds) {  // 2-D dataset -> {x:[],y:[]}
    const v = ds.value, sh = ds.shape, k = sh[1] || 1, n = sh[0];
    const x = new Float64Array(n), y = new Float64Array(n);
    for (let i = 0; i < n; i++) { x[i] = Number(v[i * k]); y[i] = Number(v[i * k + 1] || 0); }
    return { x, y };
  }

  window.parseH5ad = function (buf, name) {
    const H = window.hdf5;
    if (!H) throw new Error("HDF5 reader (jsfive) not loaded");
    const f = new H.File(buf, name || "file.h5ad");
    const get = (p) => { try { return f.get(p); } catch (e) { return null; } };

    // --- var / obs indices ---
    const varG = get("var"); if (!varG) throw new Error("no /var");
    const vIdx = (varG.attrs && varG.attrs["_index"]) || "_index";
    const genes = readStr(get("var/" + vIdx)).map(String);
    const obsG = get("obs"); if (!obsG) throw new Error("no /obs");
    const oIdx = (obsG.attrs && obsG.attrs["_index"]) || "_index";
    const ids = readStr(get("obs/" + oIdx)).map(String);
    const n = ids.length;

    // --- cell_type ---
    let cats = ["cell"], ct = new Array(n).fill(0);
    const cc = get("obs/cell_type");
    if (cc) {
      if (isGroup(cc) && cc.keys.includes("categories")) {
        cats = Array.from(cc.get("categories").value).map(String);
        ct = Array.from(toNum(cc.get("codes").value), Number);
      } else {
        const v = readStr(cc);
        cats = [...new Set(v)]; const m = new Map(cats.map((c, i) => [c, i]));
        ct = v.map((x) => m.get(x));
      }
    }

    // --- coordinates (spatial required) ---
    const sp = get("obsm/spatial");
    if (!sp) throw new Error("obsm['spatial'] is required — no spatial coordinates found");
    const S = col2(sp);
    const um = get("obsm/X_umap");
    const U = um ? col2(um) : S;

    // --- X (CSR expected) ---
    const X = get("X");
    if (!isGroup(X)) throw new Error("dense X not supported here — provide sparse (CSR) X");
    const enc = String((X.attrs["encoding-type"] || X.attrs["h5sparse_format"] || "csr"));
    const shape = X.attrs["shape"] || X.attrs["h5sparse_shape"];
    const nVar = Number(shape[1]);
    let data = toNum(X.get("data").value);
    let indices = X.get("indices").value;   // gene index per nnz (CSR)
    let indptr = X.get("indptr").value;      // per cell, length n+1
    if (enc.indexOf("csc") >= 0) throw new Error("CSC X not yet supported — re-save as CSR");
    const nnz = data.length;

    // --- per-cell summaries + top genes (from CSR rows) ---
    const umis = new Float64Array(n), nGenes = new Int32Array(n);
    const top = new Array(n);
    for (let r = 0; r < n; r++) {
      const s = Number(indptr[r]), e = Number(indptr[r + 1]);
      let sum = 0; const items = [];
      for (let k = s; k < e; k++) { const d = data[k]; sum += d; items.push([indices[k], d]); }
      umis[r] = sum; nGenes[r] = e - s;
      items.sort((a, b) => b[1] - a[1]);
      top[r] = items.slice(0, 6).map(([gi, d]) => [genes[gi], Math.round(d)]);
    }

    // --- transpose CSR -> CSC (gene-major) for gene vectors ---
    const ptr = new Uint32Array(nVar + 1);
    for (let k = 0; k < nnz; k++) ptr[indices[k] + 1]++;
    for (let g = 0; g < nVar; g++) ptr[g + 1] += ptr[g];
    const idxC = new Uint32Array(nnz), cntC = new Float32Array(nnz);
    const fill = ptr.slice(0, nVar);
    for (let r = 0; r < n; r++) {
      const s = Number(indptr[r]), e = Number(indptr[r + 1]);
      for (let k = s; k < e; k++) { const g = indices[k]; const pos = fill[g]++; idxC[pos] = r; cntC[pos] = data[k]; }
    }
    const geneVec = (g) => { const v = new Float32Array(n); for (let k = ptr[g]; k < ptr[g + 1]; k++) v[idxC[k]] = cntC[k]; return v; };
    const geneRow = new Map(genes.map((g, i) => [g.toUpperCase(), i]));

    // --- markers present ---
    const markers = {};
    for (const g of MARKERS) { const gi = geneRow.get(g); if (gi != null) markers[g] = Array.from(geneVec(gi), (x) => Math.round(x)); }

    // --- DE panel: top-dispersion HVGs + markers ---
    const disp = new Float64Array(nVar);
    for (let g = 0; g < nVar; g++) {
      const s = ptr[g], e = ptr[g + 1]; let sum = 0, sq = 0;
      for (let k = s; k < e; k++) { const d = cntC[k]; sum += d; sq += d * d; }
      const mean = sum / n, varr = sq / n - mean * mean;
      disp[g] = mean > 0 && (e - s) >= Math.min(20, n * 0.02) ? varr / mean : 0;
    }
    const order = Array.from({ length: nVar }, (_, i) => i).sort((a, b) => disp[b] - disp[a]);
    const panel = new Set(order.slice(0, 160));
    for (const g of MARKERS) { const gi = geneRow.get(g); if (gi != null) panel.add(gi); }
    const deg = {}, deg_genes = [];
    for (const gi of panel) { deg_genes.push(genes[gi]); deg[genes[gi]] = Array.from(geneVec(gi), (x) => Math.round(x)); }

    // --- morphology columns from obs (if present) ---
    const readObsNum = (col) => { const d = get("obs/" + col); return d && !isGroup(d) ? Array.from(toNum(d.value), Number) : null; };
    const morph = {}, morph_features = [];
    for (const fld of MORPH_NUC) { const v = readObsNum(fld); if (v) { morph[fld] = v.map((x) => Math.round(x * 100) / 100); morph_features.push(fld); } }
    const cellmorph = {}, cell_features = [];
    for (const fld of MORPH_CELL) { const v = readObsNum(fld); if (v) { cellmorph[fld] = v.map((x) => Math.round(x * 1000) / 1000); cell_features.push(fld); } }

    // --- TCR (optional) ---
    const readObsCat = (col) => { const d = get("obs/" + col); if (!d) return null; const v = readStr(d).map((x) => (x == null || x === "NaN" || x === "" ? null : String(x))); return v.length ? v : null; };

    const round = (a, p) => Array.from(a, (v) => Math.round(v * p) / p);
    const DATA = {
      categories: cats, n, ids,
      sx: round(S.x, 10), sy: round(S.y, 10),
      ux: round(U.x, 1000), uy: round(U.y, 1000),
      ct, umis: Array.from(umis, Math.round), genes: Array.from(nGenes),
      npk: new Array(n).fill(0), top, markers,
      tcr_a: readObsCat("tcr_alpha"), tcr_b: readObsCat("tcr_beta"),
      morph_features, morph, cell_features, cellmorph,
      deg_genes, deg,
      source: name || "uploaded .h5ad",
    };
    const expr = { n_cells: n, genes, ptr, idx: idxC, cnt: cntC };
    return { DATA, expr };
  };
})();
