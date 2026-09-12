"""Ingest: build PairedData from real or synthetic sources."""
from .morphology import MORPHOLOGY_FEATURES, cut_crops, extract_morphology
from .slidetags import load_slidetags
from .synthetic import make_synthetic

__all__ = [
    "make_synthetic",
    "load_slidetags",
    "extract_morphology",
    "cut_crops",
    "MORPHOLOGY_FEATURES",
]
