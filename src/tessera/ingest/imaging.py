"""Turn images into real per-nucleus / per-cell morphology.

Two entry points, one shared attach step:

- ``synthesize_nuclear_image`` — paint a synthetic nuclear-stain image (and,
  optionally, a membrane/cytoplasm channel giving whole-cell masks) registered to
  each nucleus's spatial (x, y), then run the real morphology pipeline. Used to
  exercise and demo the whole path before real microscopy exists.
- ``attach_morphology_from_mask`` — the *real* loader: given a genuine nuclear
  image + a segmentation label image (Cellpose/StarDist/…) and the pixel size,
  join each segmented object to a nucleus by its tissue coordinate and attach the
  extracted morphology.

Both produce: nuclear morphology (incl. circularity + integrated intensity ~ DNA
content), optional whole-cell morphology + nucleus:cytoplasm (N:C) ratio, and a
per-nucleus crop store (multi-channel when a membrane image is present).
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from skimage.draw import ellipse
from skimage.filters import gaussian

from ..core import CropStore, PairedData
from .morphology import MORPHOLOGY_FEATURES, cut_crops, extract_morphology
from .synthetic import CELL_TYPES

# nuclear params: (major-axis mean µm, sd, eccentricity, stain intensity)
PARAMS = {c[0]: {"maj_m": c[2], "maj_sd": c[3], "ecc": c[4], "inten": c[5]} for c in CELL_TYPES}
_DEFAULT = {"maj_m": 11.0, "maj_sd": 2.0, "ecc": 0.42, "inten": 175}
# whole-cell size = nucleus x cyto (drives N:C = 1/cyto^2): lymphocytes ~thin rim,
# macrophages ~abundant cytoplasm.
CYTO = {"tumour_1": 1.5, "tumour_2": 1.5, "T_CD8": 1.12, "T_CD4": 1.12, "T_reg": 1.14,
        "mono-mac": 2.0, "myeloid": 1.9, "mDC": 1.7, "pDC": 1.5, "plasma": 1.35}

_CELL_BASE = ["area", "eccentricity", "solidity", "axis_major_length", "circularity"]
_NUC_INTENSITY = ["intensity_mean", "intensity_std", "integrated_intensity"]


# --------------------------------------------------------------------------- #
# synthetic
# --------------------------------------------------------------------------- #
def synthesize_nuclear_image(
    paired: PairedData,
    stain: str = "7-AAD",
    membrane_stain: str = "WGA",
    px_per_um: float = 1.0,
    crop_size: int = 64,
    with_membrane: bool = True,
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
    Wp, Hp = int(xpx.max() + pad), int(ypx.max() + pad)

    nuc_img = np.zeros((Hp, Wp), np.float32)
    nuc_lab = np.zeros((Hp, Wp), np.int32)
    mem_img = np.zeros((Hp, Wp), np.float32)
    cell_lab = np.zeros((Hp, Wp), np.int32)

    order = np.argsort(-np.array([PARAMS.get(c, _DEFAULT)["maj_m"] for c in ct]))  # big first
    majors = np.zeros(n)
    for i in order:
        p = PARAMS.get(ct[i], _DEFAULT)
        major = max(3.0, rng.normal(p["maj_m"], p["maj_sd"])) * 0.5 * px_per_um
        minor = major * np.sqrt(max(0.05, 1 - p["ecc"] ** 2))
        majors[i] = major
        rot = rng.uniform(0, np.pi)
        if with_membrane:
            cf = CYTO.get(ct[i], 1.4)
            crr, ccc = ellipse(ypx[i], xpx[i], major * cf, minor * cf, shape=nuc_img.shape, rotation=rot)
            cell_lab[crr, ccc] = i + 1
            mem_img[crr, ccc] = np.clip(70 + rng.normal(0, 12, crr.shape), 0, 255)
        rr, cc = ellipse(ypx[i], xpx[i], major, minor, shape=nuc_img.shape, rotation=rot)
        nuc_img[rr, cc] = np.clip(p["inten"] + rng.normal(0, 14, rr.shape), 0, 255)
        nuc_lab[rr, cc] = i + 1
        if with_membrane:
            mem_img[rr, cc] = np.clip(12 + rng.normal(0, 6, rr.shape), 0, 255)  # nucleus dark in membrane ch
        # chromatin puncta
        from skimage.draw import disk
        for _ in range(int(rng.integers(2, 5))):
            ny = int(np.clip(ypx[i] + rng.normal(0, major * 0.4), 0, Hp - 1))
            nx = int(np.clip(xpx[i] + rng.normal(0, minor * 0.4), 0, Wp - 1))
            dr, dc = disk((ny, nx), max(1.0, major * 0.18), shape=nuc_img.shape)
            nuc_img[dr, dc] = np.clip(nuc_img[dr, dc] + rng.uniform(25, 55), 0, 255)

    nuc_img = np.clip(gaussian(nuc_img, 0.8 * px_per_um, preserve_range=True), 0, 255).astype(np.uint8)
    ma = _morph_aligned(nuc_lab, nuc_img, n)
    if with_membrane:
        mem_img = np.clip(gaussian(mem_img, 0.8 * px_per_um, preserve_range=True), 0, 255).astype(np.uint8)
        cma = _morph_aligned(cell_lab, None, n)
        stack = np.stack([nuc_img, mem_img], axis=-1)
        chan = [stain, membrane_stain]
    else:
        cma = None
        stack = nuc_img
        chan = [stain]

    crops = cut_crops(stack, centroids=list(zip(ypx, xpx)), size=crop_size)
    _attach(paired, ma, cma, crops, chan, stain, crops_path, "synthetic", px_per_um)

    # a registered tissue overview (pseudo-H&E) for the viewer's image underlay
    nrm = nuc_img.astype(np.float32) / 255.0
    mrm = (mem_img.astype(np.float32) / 255.0) if with_membrane else np.zeros_like(nrm)
    hema = np.array([80, 45, 135]) / 255.0   # nuclei -> haematoxylin (purple)
    eos = np.array([240, 150, 190]) / 255.0  # cytoplasm -> eosin (pink)
    he = np.ones((*nuc_img.shape, 3), np.float32)
    he *= (1 - mrm[..., None] * (1 - eos))
    he *= (1 - nrm[..., None] * (1 - hema))
    he = np.clip(he * 255, 0, 255).astype(np.uint8)
    step = max(1, int(np.ceil(max(he.shape[:2]) / 1100)))
    a = paired.adata
    a.uns["tissue_image"] = np.ascontiguousarray(he[::step, ::step])
    a.uns["tissue_fy"] = True  # synthetic image row increases with spatial y

    if return_image:
        return paired, nuc_img, (mem_img if with_membrane else None)
    return paired


# --------------------------------------------------------------------------- #
# real loader
# --------------------------------------------------------------------------- #
def attach_morphology_from_mask(
    paired: PairedData,
    nuclear_image,
    labels,
    cell_labels=None,
    membrane_image=None,
    um_per_px: float = 1.0,
    origin: tuple = (0.0, 0.0),
    crop_size: int = 64,
    crops_path: str = "crops_real.zarr",
    stain: str = "DAPI",
    membrane_stain: str = "membrane",
    match_tol_px: float | None = None,
    return_unmatched: bool = False,
):
    """Attach morphology from a real nuclear image + segmentation mask.

    ``labels`` is an integer segmentation image (0 = background). Each nucleus's
    tissue (x, y) is mapped to pixels via ``um_per_px`` and ``origin`` (the tissue
    coordinate of pixel (0, 0)); the segment covering that pixel — else the nearest
    segment centroid within ``match_tol_px`` — supplies the morphology. Images and
    label arrays may be passed as arrays or file paths.
    """
    nuclear_image = _asarray(nuclear_image)
    labels = _asarray(labels).astype(np.int32)
    membrane_image = _asarray(membrane_image) if membrane_image is not None else None
    cell_labels = _asarray(cell_labels).astype(np.int32) if cell_labels is not None else None

    a = paired.adata
    n = a.n_obs
    xy = paired.spatial[:, :2]
    px = (xy[:, 0] - origin[0]) / um_per_px
    py = (xy[:, 1] - origin[1]) / um_per_px
    H, W = labels.shape

    nuc_morph = extract_morphology(labels, intensity_image=nuclear_image)
    cell_morph = extract_morphology(cell_labels, intensity_image=membrane_image) if cell_labels is not None else None

    if match_tol_px is None:
        match_tol_px = crop_size / 2
    nuc_lbl = _match_labels(labels, nuc_morph, px, py, H, W, match_tol_px)
    cell_lbl = (_match_labels(cell_labels, cell_morph, px, py, H, W, match_tol_px)
                if cell_labels is not None else None)

    ma = nuc_morph.reindex(nuc_lbl).reset_index(drop=True)
    cma = cell_morph.reindex(cell_lbl).reset_index(drop=True) if cell_morph is not None else None

    if membrane_image is not None:
        stack = np.stack([nuclear_image, membrane_image], axis=-1)
        chan = [stain, membrane_stain]
    else:
        stack = nuclear_image
        chan = [stain]
    crops = cut_crops(stack, centroids=list(zip(py, px)), size=crop_size)
    _attach(paired, ma, cma, crops, chan, stain, crops_path, "imaged", um_per_px and 1.0 / um_per_px)

    n_matched = int(np.sum(~pd.isna(nuc_lbl)))
    print(f"matched {n_matched}/{n} nuclei to segments")
    if return_unmatched:
        return paired, np.where(pd.isna(nuc_lbl))[0]
    return paired


# --------------------------------------------------------------------------- #
# shared helpers
# --------------------------------------------------------------------------- #
def _morph_aligned(labels, intensity, n):
    m = extract_morphology(labels, intensity_image=intensity)
    m = m.reindex(np.arange(1, n + 1))
    m.index = np.arange(n)
    return m


def _match_labels(labels, morph, px, py, H, W, tol):
    """Return, per nucleus, the label id covering its pixel (or nearest centroid)."""
    out = np.full(len(px), np.nan)
    cent = None
    for i in range(len(px)):
        x, y = px[i], py[i]
        xi, yi = int(round(x)), int(round(y))
        if 0 <= yi < H and 0 <= xi < W and labels[yi, xi] > 0:
            out[i] = labels[yi, xi]
            continue
        if cent is None:
            cent = morph[["centroid_x", "centroid_y"]].to_numpy()
            cent_lbl = morph.index.to_numpy()
        d = np.hypot(cent[:, 0] - x, cent[:, 1] - y)
        j = int(np.argmin(d))
        if d[j] <= tol:
            out[i] = cent_lbl[j]
    return out


def _attach(paired, ma, cma, crops, channel_names, stain, crops_path, source, px_per_um):
    a = paired.adata
    nuc_feats = [f for f in MORPHOLOGY_FEATURES if f in ma.columns]
    nuc_feats += [c for c in _NUC_INTENSITY if c in ma.columns]
    for c in nuc_feats:
        a.obs[c] = np.asarray(ma[c], dtype=float)
    a.uns[PairedData.MORPH_KEY] = nuc_feats

    cell_feats = []
    if cma is not None:
        for c in _CELL_BASE:
            if c in cma.columns:
                a.obs["cell_" + c] = np.asarray(cma[c], dtype=float)
                cell_feats.append("cell_" + c)
        with np.errstate(divide="ignore", invalid="ignore"):
            nc = np.asarray(ma["area"], float) / np.asarray(cma["area"], float)
        a.obs["nc_ratio"] = np.round(nc, 3)
        cell_feats.append("nc_ratio")
    a.uns["cell_features"] = cell_feats
    a.uns["stain"] = stain
    a.uns["morphology_source"] = source
    if px_per_um:
        a.uns["px_per_um"] = float(px_per_um)

    store = CropStore.create(
        crops_path, cell_ids=list(a.obs_names),
        shape=(crops.shape[-2], crops.shape[-1]),
        channels=crops.shape[1], channel_names=channel_names, dtype="uint8",
    )
    store._arr[:] = crops
    paired.crops = store


def _asarray(x):
    if isinstance(x, str):
        from skimage.io import imread
        return imread(x)
    return np.asarray(x)
