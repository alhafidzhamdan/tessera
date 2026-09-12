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
      (pan/zoom scope, lineage/subtype/gene/depth colouring, per-nucleus inspector)
- [ ] **Phase 3** — morphology ↔ omics association statistics
- [ ] **Phase 4** — joint embedding
- [ ] **Phase 5** — cross-modal prediction
- [ ] R-side readers for stats + figures

## Status

Pre-alpha scaffold. Synthetic data only; real slide-tags ingest pending the
paper + SCP2176 download.
