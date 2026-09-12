"""Synthesize a nuclear-stain image registered to existing nuclei, then run the
real morphology pipeline on it.

This is the bridge that turns the designed-in morphology slot into real columns.
Given a :class:`~tessera.core.PairedData` that already has per-nucleus tissue
coordinates (e.g. slide-tags), it paints a synthetic nuclear-stain image at each
nucleus's (x, y) with cell-type-dependent nuclear morphology and chromatin
texture, builds the segmentation label image, and then runs the *same*
``extract_morphology`` + ``cut_crops`` used for genuine microscopy. The result is
attached to the object: morphology feature columns in ``obs``, ``uns`` schema, and
a per-nucleus crop store.

Default stain is 7-AAD (7-aminoactinomycin D), a DNA intercalator. Swap in a real
registered image later and the downstream code is identical.
"""
from __future__ import annotations

import numpy as np
from skimage.draw import disk, ellipse
from skimage.filters import gaussian

from ..core import CropStore, PairedData
from .morphology import MORPHOLOGY_FEATURES, cut_crops, extract_morphology
from .synthetic import CELL_TYPES

# per-cell-type nuclear params: (major-axis mean µm, sd, eccentricity, stain intensity)
PARAMS = {c[0]: {"maj_m": c[2], "maj_sd": c[3], "ecc": c[4], "inten": c[5]} for c in CELL_TYPES}
_DEFAULT = {"maj_m": 11.0, "maj_sd": 2.0, "ecc": 0.42, "inten": 175}


def synthesize_nuclear_image(
    paired: PairedData,
    stain: str = "7-AAD",
    px_per_um: float = 1.0,
    crop_size: int = 64,
    crops_path: str = "crops_stain.zarr",
    seed: int = 0,
    return_image: bool = False,
):
    if not paired.has_spatial:
        raise ValueError("nuclei need spatial coordinates to register an image")
    rng = np.random.default_rng(seed)
    a = paired.adata
    ct = a.obs["cell_type"].astype(str).to_numpy()
    n = a.n_obs

    xy = paired.spatial[:, :2]
    pad = crop_size
    xmin, ymin = xy[:, 0].min(), xy[:, 1].min()
    xpx = (xy[:, 0] - xmin) * px_per_um + pad
    ypx = (xy[:, 1] - ymin) * px_per_um + pad
    Wp = int(xpx.max() + pad)
    Hp = int(ypx.max() + pad)

    img = np.zeros((Hp, Wp), np.float32)
    labels = np.zeros((Hp, Wp), np.int32)
    for i in range(n):
        p = PARAMS.get(ct[i], _DEFAULT)
        major = max(3.0, rng.normal(p["maj_m"], p["maj_sd"])) * 0.5 * px_per_um
        minor = major * np.sqrt(max(0.05, 1 - p["ecc"] ** 2))
        rot = rng.uniform(0, np.pi)
        rr, cc = ellipse(ypx[i], xpx[i], major, minor, shape=img.shape, rotation=rot)
        img[rr, cc] = np.clip(p["inten"] + rng.normal(0, 14, rr.shape), 0, 255)
        labels[rr, cc] = i + 1
        # chromatin: a few brighter heterochromatin / nucleolar puncta
        for _ in range(int(rng.integers(2, 5))):
            ny = int(np.clip(ypx[i] + rng.normal(0, major * 0.4), 0, Hp - 1))
            nx = int(np.clip(xpx[i] + rng.normal(0, minor * 0.4), 0, Wp - 1))
            dr, dc = disk((ny, nx), max(1.0, major * 0.18), shape=img.shape)
            img[dr, dc] = np.clip(img[dr, dc] + rng.uniform(25, 55), 0, 255)

    img = gaussian(img, sigma=0.8 * px_per_um, preserve_range=True)
    img = np.clip(img, 0, 255).astype(np.uint8)

    # real morphology pipeline on the (synthetic) segmented image
    morph = extract_morphology(labels, intensity_image=img)
    feat = [f for f in MORPHOLOGY_FEATURES if f in morph.columns]
    feat += [c for c in ("intensity_mean", "intensity_std") if c in morph.columns]
    ma = morph.reindex(np.arange(1, n + 1))
    ma.index = np.arange(n)
    for c in feat:
        a.obs[c] = ma[c].to_numpy()
    a.uns[PairedData.MORPH_KEY] = feat
    a.uns["stain"] = stain
    a.uns["morphology_source"] = "synthetic"
    a.uns["px_per_um"] = float(px_per_um)

    crops = cut_crops(img, centroids=list(zip(ypx, xpx)), size=crop_size)
    store = CropStore.create(
        crops_path, cell_ids=list(a.obs_names), shape=(crop_size, crop_size),
        channels=1, channel_names=[stain], dtype="uint8",
    )
    store._arr[:] = crops
    paired.crops = store

    if return_image:
        return paired, img, labels
    return paired
