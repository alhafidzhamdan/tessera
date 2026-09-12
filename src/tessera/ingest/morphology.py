"""Morphology feature extraction and crop cutting.

This is the *real* path used once a co-registered nuclear image is available:
segment the image (Cellpose / StarDist / classical) to a label image, then

    features = extract_morphology(label_image, intensity_image)
    crops    = cut_crops(intensity_image, features[["centroid_y","centroid_x"]])

The same functions back the synthetic generator, so demo and real data flow
through identical code.
"""
from __future__ import annotations

from typing import Optional, Sequence

import numpy as np
import pandas as pd
from skimage.measure import regionprops_table

# Shape features (unitless or in pixels); intensity features appended when an
# intensity image is supplied.
_SHAPE_PROPS = (
    "label",
    "area",
    "perimeter",
    "eccentricity",
    "solidity",
    "extent",
    "orientation",
    "axis_major_length",
    "axis_minor_length",
    "equivalent_diameter_area",
    "centroid",
)
_INTENSITY_PROPS = ("intensity_mean", "intensity_max", "intensity_min")

MORPHOLOGY_FEATURES = [
    "area",
    "perimeter",
    "eccentricity",
    "solidity",
    "extent",
    "orientation",
    "axis_major_length",
    "axis_minor_length",
    "equivalent_diameter_area",
    "circularity",
]


def extract_morphology(
    label_image: np.ndarray,
    intensity_image: Optional[np.ndarray] = None,
) -> pd.DataFrame:
    """Per-label morphology table, indexed by integer label.

    Adds ``centroid_y`` / ``centroid_x`` and (if ``intensity_image`` given)
    intensity mean/max/min plus a derived ``intensity_std``.
    """
    props = list(_SHAPE_PROPS)
    kwargs = {}
    if intensity_image is not None:
        props += list(_INTENSITY_PROPS)
        kwargs["intensity_image"] = intensity_image
    tbl = regionprops_table(label_image, properties=props, **kwargs)
    df = pd.DataFrame(tbl)
    df = df.rename(
        columns={"centroid-0": "centroid_y", "centroid-1": "centroid_x"}
    )
    # circularity: 4*pi*area / perimeter^2 (1 = perfect circle, ->0 = irregular)
    with np.errstate(divide="ignore", invalid="ignore"):
        df["circularity"] = (4 * np.pi * df["area"] / (df["perimeter"] ** 2))
    df["circularity"] = df["circularity"].replace([np.inf, -np.inf], np.nan).clip(upper=1.0)
    if intensity_image is not None:
        # regionprops has no std; compute it cheaply per label.
        stds = {}
        for lbl in df["label"]:
            stds[lbl] = float(intensity_image[label_image == lbl].std())
        df["intensity_std"] = df["label"].map(stds)
        # integrated intensity = total signal in the region (area x mean) ~ DNA content
        df["integrated_intensity"] = (df["area"] * df["intensity_mean"]).round(1)
    return df.set_index("label")


def cut_crops(
    image: np.ndarray,
    centroids: Sequence[Sequence[float]],
    size: int = 64,
) -> np.ndarray:
    """Cut fixed-size (channels, size, size) crops centred on each centroid.

    ``image`` is (H, W) grayscale or (H, W, C). Out-of-bounds regions are
    zero-padded. Returns (n, channels, size, size), dtype matching ``image``.
    """
    if image.ndim == 2:
        image = image[..., None]
    h, w, c = image.shape
    half = size // 2
    out = np.zeros((len(centroids), c, size, size), dtype=image.dtype)
    for i, (cy, cx) in enumerate(centroids):
        cy, cx = int(round(cy)), int(round(cx))
        y0, y1 = cy - half, cy - half + size
        x0, x1 = cx - half, cx - half + size
        sy0, sx0 = max(0, y0), max(0, x0)
        sy1, sx1 = min(h, y1), min(w, x1)
        if sy1 <= sy0 or sx1 <= sx0:
            continue
        patch = image[sy0:sy1, sx0:sx1, :]
        oy0, ox0 = sy0 - y0, sx0 - x0
        out[i, :, oy0 : oy0 + patch.shape[0], ox0 : ox0 + patch.shape[1]] = (
            np.moveaxis(patch, -1, 0)
        )
    return out
