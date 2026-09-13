"""Generate the drop-a-file web app (web/viewer_app.html) from the baked viewer
template. The app shares the viewer's markup + logic but, instead of a baked-in
dataset, shows a drop zone and boots the viewer from an uploaded Tessera bundle
JSON (see demo/export_bundle.py).

    python demo/make_app.py
"""
import os

TPL = "web/viewer_template.html"
OUT = "web/viewer_app.html"

DROP_UI = """
<script src="https://cdn.jsdelivr.net/npm/jsfive@0.3.10/dist/browser/hdf5.js"></script>
<script src="tessera_h5ad.js"></script>
<div id="drop-overlay">
  <div id="dropzone">
    <div class="dz-title">Tessera viewer</div>
    <div class="dz-sub">Drop an <code>.h5ad</code> (AnnData) or a Tessera bundle <code>.json</code> to explore it — spatial &amp; UMAP, colour by any gene, cell table, lasso, differential expression, A/B compare. It just needs an expression matrix and <code>obsm['spatial']</code>.</div>
    <label class="dz-btn" for="drop-file">Choose file…</label>
    <input id="drop-file" type="file" accept=".h5ad,.h5,.hdf5,.json,application/json" hidden>
    <button class="dz-demo" id="drop-demo">Load demo dataset</button>
    <div class="dz-msg" id="drop-msg"></div>
    <div class="dz-foot">Reads <code>.h5ad</code> directly in your browser (no server, nothing uploaded). Very large files are read in memory, so a processed AnnData works best. Or make a lightweight bundle: <code>python demo/export_bundle.py in.tessera out.tessera.json</code>.</div>
  </div>
</div>
<style>
#drop-overlay{position:fixed;inset:0;z-index:100;display:flex;align-items:center;justify-content:center;
  background:var(--bg);padding:24px}
#dropzone{max-width:520px;width:100%;border:1.5px dashed var(--border);border-radius:16px;
  background:var(--panel);box-shadow:var(--shadow);padding:34px 30px;text-align:center}
#dropzone.over{border-color:var(--accent);background:var(--panel-2)}
.dz-title{font-size:20px;font-weight:600;letter-spacing:-.01em;color:var(--ink)}
.dz-sub{font-size:13px;color:var(--muted);margin:8px auto 20px;max-width:400px;line-height:1.5}
.dz-btn{display:inline-block;cursor:pointer;background:var(--accent);color:var(--accent-ink);
  font-weight:500;font-size:13px;padding:9px 18px;border-radius:9px}
.dz-btn:hover{filter:brightness(1.05)}
.dz-demo{margin-left:10px;cursor:pointer;background:var(--panel-2);color:var(--ink-2);
  border:1px solid var(--border);font-family:inherit;font-size:13px;padding:9px 16px;border-radius:9px}
.dz-demo:hover{color:var(--ink)}
.dz-msg{font-size:12px;color:#C2543A;margin-top:14px;min-height:16px}
.dz-foot{font-size:11px;color:var(--muted);margin-top:22px;line-height:1.5}
.dz-foot code{font-family:"IBM Plex Mono",monospace;background:var(--panel-2);padding:1px 5px;border-radius:4px}
</style>
<script>
(function(){
  const ov=document.getElementById('drop-overlay');
  const msg=document.getElementById('drop-msg');
  function boot(d){try{ov.style.display='none';bootTessera(d);}catch(e){ov.style.display='flex';msg.textContent='Could not open: '+e.message;}}
  function readFile(file){
    if(/\\.(h5ad|h5|hdf5)$/i.test(file.name)){
      if(!window.parseH5ad){msg.textContent='HDF5 reader still loading — try again in a moment';return;}
      msg.textContent='reading '+file.name+' …';
      const fr=new FileReader();
      fr.onerror=()=>msg.textContent='Could not read the file.';
      fr.onload=()=>setTimeout(()=>{try{const r=window.parseH5ad(fr.result,file.name);
        window.__PRE_EXPR=r.expr;boot(r.DATA);}
        catch(e){msg.textContent='.h5ad: '+e.message;}},10);
      fr.readAsArrayBuffer(file);return;
    }
    const fr=new FileReader();
    fr.onerror=()=>msg.textContent='Could not read the file.';
    fr.onload=()=>{try{window.__PRE_EXPR=null;boot(JSON.parse(fr.result));}catch(e){msg.textContent='Not a valid Tessera bundle JSON: '+e.message;}};
    fr.readAsText(file);}
  document.getElementById('drop-file').addEventListener('change',e=>{if(e.target.files[0])readFile(e.target.files[0]);});
  const dz=document.getElementById('dropzone');
  ['dragover','dragenter'].forEach(ev=>dz.addEventListener(ev,e=>{e.preventDefault();dz.classList.add('over');}));
  ['dragleave','dragend'].forEach(ev=>dz.addEventListener(ev,e=>{e.preventDefault();dz.classList.remove('over');}));
  dz.addEventListener('drop',e=>{e.preventDefault();dz.classList.remove('over');
    const f=e.dataTransfer.files[0];if(f)readFile(f);});
  document.getElementById('drop-demo').addEventListener('click',()=>{
    msg.textContent='loading demo…';window.__PRE_EXPR=null;
    if(window.DEMO_BUNDLE){boot(window.DEMO_BUNDLE);return;}
    const s=document.createElement('script');s.src='demo_bundle.js';
    s.onload=()=>window.DEMO_BUNDLE?boot(window.DEMO_BUNDLE):(msg.textContent='demo unavailable');
    s.onerror=()=>msg.textContent='demo unavailable';document.body.appendChild(s);});
})();
</script>
"""


def main():
    tpl = open(TPL).read()
    app = tpl.replace('<script src="expr_meta.js"></script>\n', "")
    app = app.replace("<title>Melanoma Slide-tags Atlas</title>", "<title>Tessera Viewer</title>")
    # genericise the baked header for an arbitrary uploaded dataset
    app = app.replace("<h1>Melanoma Slide-tags Atlas</h1>", "<h1>Tessera</h1>")
    app = app.replace(
        '<div class="sub"><b id="ncells">—</b> nuclei · single-nucleus <b>RNA</b> + <b>ATAC</b> + spatial · <span class="mono">SCP2176</span></div>',
        '<div class="sub"><b id="ncells">—</b> nuclei · paired morphology + multi-omics</div>')
    assert "const DATA = __TESSERA_DATA__;" in app, "DATA anchor not found"
    app = app.replace("const DATA = __TESSERA_DATA__;",
                      "let DATA=null;\nfunction bootTessera(d){DATA=d;")
    idx = app.rfind("</script>")
    app = app[:idx] + "}\n</script>" + app[idx + len("</script>"):]
    app += DROP_UI
    with open(OUT, "w") as fh:
        fh.write(app)
    print(f"wrote {OUT}  ({len(app.encode())/1e6:.2f} MB)")


if __name__ == "__main__":
    main()
