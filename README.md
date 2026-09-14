# Tessera

Pair single-nucleus **morphology** to **multi-omic** data — spatially.

Each nucleus is one row, joined by a stable `cell_id`. Built first around
**slide-tags multiome** (snRNA + snATAC + spatial): the spatial coordinate is the
bridge that lets a co-registered nuclear image (DAPI / H&E / 7-AAD) attach
morphology to each nucleus's transcriptome and chromatin profile. Morphology is
optional and designed-in — everything works without an image and fills in, per
nucleus, the moment one arrives.

> Working name — a *tessera* is one tile of a mosaic; each nucleus is one tile of
> the tissue.

## Data model

AnnData backbone + a Zarr crop store, all serialising to Zarr (readable from the
web viewer and from R):

| Where | Holds |
|---|---|
| `adata.X` | RNA counts |
| `adata.mods['atac']` | ATAC peaks (any extra modality) |
| `adata.obs` | `cell_type`, morphology + cell (N:C) features |
| `adata.obsm['spatial']` | tissue (x, y) — the join key to imaging |
| `adata.obsm['X_umap']` | embedding |
| crop store (Zarr) | per-nucleus image crop, keyed by `cell_id` |

## Quick start

```bash
pip install -e .
python demo/build_demo.py demo_out          # synthetic demo bundle + figure
```

```python
from tessera.core import PairedData
p = PairedData.load("demo_out/melanoma_demo.tessera")
p.catalogue()          # what's paired
p.catalogue_table()    # per-cell-type counts + mean morphology
```

## Load your own data

Only an **omics matrix** and **per-nucleus `spatial` (x, y)** are required;
`cell_type`, a UMAP, ATAC, and imaging are used when present.

```python
from tessera.ingest import from_anndata, build_paired, load_slidetags

p = from_anndata(adata)                              # AnnData with obsm['spatial']
p = build_paired(rna, var_names, cell_ids, spatial,  # or plain arrays
                 cell_type=labels, umap=umap, atac=atac, atac_peaks=peaks)
p = load_slidetags("path/to/SCP_dir")                # or a slide-tags/SCP directory
p.save("mydata.tessera")
```

Add real morphology once you have imaging (joins by tissue coordinate):

```python
from tessera.ingest import attach_morphology_from_mask
attach_morphology_from_mask(p, nuclear_image="dapi.tif", labels="mask.tif", um_per_px=0.5)
# no image yet? synthesize_nuclear_image(p) fabricates a 7-AAD test image
```

## View it

- **Desktop (napari):** `PYTHONPATH=src python -m tessera.viewer.napari_view mydata.tessera`
- **Drop-a-file web app** (`web/viewer_app.html`): drop an **`.h5ad`** (parsed
  in-browser via [jsfive](https://github.com/usnistgov/jsfive) — nothing uploaded,
  with any-gene search) or a bundle from `demo/export_bundle.py`.
- **Baked page:** `python demo/export_web.py mydata.tessera web/viewer_template.html out/index.html`

The web viewer does spatial/UMAP, colour by lineage / subtype / any gene / depth,
a cell table, lasso region stats + differential expression, saved regions,
direct A/B compare, and a **tissue-image (H&E/DAPI) underlay** aligned to the
spatial map (load a downsampled overview, then opacity + offset/scale/rotate/flip).

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
      drop an `.h5ad` (parsed in-browser via jsfive, with any-gene search) or a
      Tessera bundle JSON — no re-export, nothing uploaded to a server
- [x] **Tissue-image underlay** — register an H&E/DAPI overview under the spatial
      nuclei (opacity + offset/scale/rotate/flip) to check slide-tags alignment
- [ ] **Phase 3** — morphology ↔ omics association statistics
- [ ] **Phase 4** — joint embedding
- [ ] **Phase 5** — cross-modal prediction
- [ ] R-side readers for stats + figures

## Status

Pre-alpha. Real SCP2176 slide-tags melanoma ingested; morphology demonstrated via
synthetic 7-AAD imaging pending real microscopy.

## License

**Non-commercial** — [PolyForm Noncommercial License 1.0.0](LICENSE.md): free for
research, teaching, and non-profit / government use; contact the author for a
commercial licence. The SCP2176 data is **not** included and is governed by its
own Broad Single Cell Portal terms; Tessera is not affiliated with the slide-tags
authors.
