# Tessera

Pair single-nucleus **morphology** to **multi-omic** data — spatially.

Tessera is a paired-data container + viewer for experiments where each nucleus
has (or will have) both an image and a multi-omic profile. It is built first
around **slide-tags multiome** (snRNA-seq + snATAC-seq with spatial
coordinates), where the spatial coordinate is the bridge that lets a
co-registered nuclear image (DAPI / H&E / IF) be paired to each nucleus's
transcriptome and chromatin profile.

> Working name — easily renamed. "Tessera" = a tile in a mosaic; each nucleus
> is one tile of the tissue.

## Why this design

- **One row = one nucleus**, joined by a stable `cell_id`.
- **AnnData backbone** so the whole scverse/Bioconductor ecosystem is available
  for the later association / embedding / prediction phases.
- **Zarr crop store**, one chunk per nucleus → instant random-access to any
  nucleus's image crop, at 10⁵–10⁶ nuclei.
- **Everything serialises to Zarr** → readable from a web viewer and from R
  (`zellkonverter` / `arrow`) with no re-export. This is the hybrid split:
  Python for speed + imaging + viewer, R for stats + publication figures.
- The **morphology slot is designed-in but optional**: works today on
  slide-tags RNA+ATAC+spatial with the slot empty; fills in the moment a
  registered image arrives, via spatial coordinates.

## Data model

| Where | Holds |
|---|---|
| `adata.X` | RNA counts |
| `adata.layers['atac']` | ATAC gene activity / peaks |
| `adata.obs` | `cell_type`, morphology features, centroids |
| `adata.obsm['spatial']` | tissue (x, y) — the bridge to imaging |
| `adata.obsm['X_umap']` | embedding |
| crop store (Zarr) | per-nucleus image crop, keyed by `cell_id` |

## Quick start

```bash
pip install -e .
python demo/build_demo.py demo_out   # synthetic slide-tags melanoma + figure
```

```python
from tessera.core import PairedData
p = PairedData.load("demo_out/melanoma_demo.tessera")
p.catalogue()                 # what's paired
p.get_crop(p.cell_ids[0])     # instant crop access
p.catalogue_table()           # per-cell-type counts + mean morphology
```

## Loading your own data

### What data is needed

| | Field | Notes |
|---|---|---|
| **Required** | omics matrix | one row per nucleus (e.g. RNA counts, cells × genes) |
| **Required** | spatial (x, y) | per-nucleus tissue coordinates — the join key to imaging |
| Recommended | `cell_type` | per-nucleus labels (drives colouring/composition) |
| Recommended | embedding | a 2-D UMAP/embedding per nucleus |
| Optional | ATAC | second modality, peaks × cells (any extra modality works) |
| Optional | TCR / metadata | any extra per-nucleus columns |
| Optional | imaging | a registered nuclear image **+ segmentation mask** → real morphology |

Morphology is optional and designed-in: everything works without an image; when a
co-registered nuclear image arrives, morphology fills in per nucleus **via the
spatial coordinate**.

### Three ways to load

**1. From an AnnData** (scanpy / scverse) — just needs `obsm['spatial']`:

```python
from tessera.ingest import from_anndata
p = from_anndata(adata)                 # uses obs['cell_type'], obsm['X_umap'] if present
p = from_anndata(adata, atac=atac_adata)  # optional second modality
p.save("mydata.tessera")
```

**2. From plain arrays / matrices:**

```python
from tessera.ingest import build_paired
p = build_paired(rna, var_names, cell_ids, spatial,
                 cell_type=labels, umap=umap, atac=atac, atac_peaks=peak_ids)
p.save("mydata.tessera")
```

**3. From a slide-tags / SCP directory** (SCP2176 layout — 10x `matrix.mtx`
triplet, `*_spatial.csv`, `*_cluster.csv`, metadata, ATAC csv, TCR):

```python
from tessera.ingest import load_slidetags
p = load_slidetags("path/to/SCP_dir", prefix="HumanMelanomaMultiome")
p.save("mydata.tessera")
```

### Add real morphology from imaging

```python
from tessera.ingest import attach_morphology_from_mask
attach_morphology_from_mask(p, nuclear_image="dapi.tif", labels="mask.tif",
                            um_per_px=0.5, origin=(x0, y0))   # joins by tissue coord
p.save("mydata.tessera")
```

(no image yet? `synthesize_nuclear_image(p)` fabricates a 7-AAD test image to
exercise the pipeline.)

### View it

**Desktop (napari)** — reads any bundle directly; any gene/peak, full data in memory:

```bash
PYTHONPATH=src python -m tessera.viewer.napari_view mydata.tessera
```

**Drop-a-file web app** (`web/viewer_app.html`) — a generic hosted viewer with an
upload zone. Export your bundle to a single JSON and drop it in; nothing is sent
to a server (the file is read in the browser):

```bash
python demo/export_bundle.py mydata.tessera mydata.tessera.json
# open web/viewer_app.html (or the hosted app) and drop mydata.tessera.json
python demo/make_app.py            # regenerate web/viewer_app.html from the template
```

**Baked page** (`demo/export_web.py`) — one dataset compiled into its own
self-contained page + companion files (adds full any-gene search):

```bash
python demo/export_web.py mydata.tessera web/viewer_template.html out/index.html
```

The drop-a-file app and baked page cover the current in-browser bundle format;
native `.h5ad` upload would need an embedded HDF5/WASM reader (the browser sandbox
blocks fetching one), so for now convert to a bundle with `export_bundle.py`.

## Roadmap

- [x] **Phase 1** — paired container + catalogue + morphology extraction
- [x] Synthetic slide-tags-like generator (RNA + ATAC + spatial + image)
- [x] **Real ingest** — SCP2176 melanoma multiome (RNA 90k genes + ATAC 53k
      peaks + spatial + UMAP + TCR), 2,535 nuclei → `data/processed/melanoma.tessera`
- [x] Multi-modal container (RNA primary + ATAC peaks sharing obs)
- [x] **Phase 2a** — napari desktop viewer (`tessera.viewer.napari_view`):
      spatial/UMAP toggle, colour by cell_type/gene/peak, click → profile panel
- [x] **Phase 2b** — web viewer: `demo/export_web.py` builds a compact payload
      into `web/viewer_template.html` → a self-contained, shareable HTML page
      (pan/zoom scope, colour by lineage/subtype/any-gene/depth, cell table, lasso
      region stats + differential expression, saved regions, A/B compare)
- [x] **Drop-a-file web app** (`web/viewer_app.html`) — generic hosted viewer;
      drop a Tessera bundle JSON (`demo/export_bundle.py`) to explore, no re-export
- [ ] **Phase 3** — morphology ↔ omics association statistics
- [ ] **Phase 4** — joint embedding
- [ ] **Phase 5** — cross-modal prediction
- [ ] R-side readers for stats + figures

## Status

Pre-alpha. Real SCP2176 slide-tags melanoma ingested; morphology demonstrated via
synthetic 7-AAD imaging pending real microscopy.

## License

**Non-commercial.** Tessera is licensed under the
[PolyForm Noncommercial License 1.0.0](LICENSE.md): free to use, modify, and
share for any non-commercial purpose — including academic research, teaching, and
non-profit / government use — but not for commercial use. For a commercial
licence, contact the author.

Note: the SCP2176 melanoma data is **not** included in this repository and is
governed by its own Broad Single Cell Portal terms; obtain it directly from the
portal. Tessera is not affiliated with the slide-tags authors.
