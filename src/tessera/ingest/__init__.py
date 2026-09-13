"""Ingest: build PairedData from real or synthetic sources."""
from .generic import build_paired, from_anndata
from .imaging import attach_morphology_from_mask, synthesize_nuclear_image
from .morphology import MORPHOLOGY_FEATURES, cut_crops, extract_morphology
from .slidetags import load_slidetags
from .synthetic import make_synthetic

__all__ = [
    "from_anndata",
    "build_paired",
    "load_slidetags",
    "make_synthetic",
    "synthesize_nuclear_image",
    "attach_morphology_from_mask",
    "extract_morphology",
    "cut_crops",
    "MORPHOLOGY_FEATURES",
]
