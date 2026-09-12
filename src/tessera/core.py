"""Core paired-data container for Tessera.

One row = one nucleus. Everything is joined by a stable ``cell_id``:

- ``adata.X``                 omic matrix (e.g. RNA)
- ``adata.layers[...]``       additional omic modalities (e.g. ATAC gene activity)
- ``adata.obs``               per-nucleus metadata + morphology features
- ``adata.obsm['spatial']``   tissue (x, y) per nucleus -- the canonical join
                              key (slide-tags bead coordinates); the bridge to
                              imaging
- ``adata.obs['centroid_x/y']`` pixel location of the nucleus *within a
                              registered image* (exists only once imaging is
                              paired; equals ``spatial`` in the synthetic demo,
                              where the image is the tissue)
- ``adata.obsm['X_umap']``    embeddings
- crop store (Zarr)           per-nucleus image crop, keyed by ``cell_id``

The image-crop store is a chunked Zarr array with one chunk per nucleus, so any
single nucleus's crop loads with a single-chunk read (fast random access). The
whole object serialises to Zarr, which keeps it readable from a browser viewer
and from R (via ``zellkonverter`` / ``arrow``) without re-export.
"""
from __future__ import annotations

import json
import os
from typing import Iterable, Optional, Sequence

import numpy as np
import pandas as pd
import zarr


# --------------------------------------------------------------------------- #
# Image-crop store
# --------------------------------------------------------------------------- #
class CropStore:
    """Chunked Zarr store of fixed-size per-nucleus image crops.

    Layout: array ``crops`` of shape (n_cells, channels, H, W), chunked
    (1, channels, H, W) so one crop == one chunk. ``cell_ids`` (row order) and
    ``channel_names`` live in the group attrs.
    """

    def __init__(self, group: "zarr.Group"):
        self._g = group
        self._arr = group["crops"]
        self.cell_ids = list(group.attrs["cell_ids"])
        self.channel_names = list(group.attrs.get("channel_names", []))
        self._index = {cid: i for i, cid in enumerate(self.cell_ids)}

    # -- construction ------------------------------------------------------- #
    @classmethod
    def create(
        cls,
        path: str,
        cell_ids: Sequence[str],
        shape: tuple[int, int] = (64, 64),
        channels: int = 1,
        channel_names: Optional[Sequence[str]] = None,
        dtype: str = "uint8",
    ) -> "CropStore":
        n = len(cell_ids)
        h, w = shape
        g = zarr.open_group(path, mode="w")
        _create_array(
            g,
            "crops",
            shape=(n, channels, h, w),
            chunks=(1, channels, h, w),
            dtype=dtype,
        )
        g.attrs["cell_ids"] = list(map(str, cell_ids))
        g.attrs["channel_names"] = list(channel_names) if channel_names else [
            f"ch{i}" for i in range(channels)
        ]
        return cls(g)

    @classmethod
    def open(cls, path: str, mode: str = "r") -> "CropStore":
        return cls(zarr.open_group(path, mode=mode))

    # -- access ------------------------------------------------------------- #
    @property
    def n_cells(self) -> int:
        return self._arr.shape[0]

    @property
    def crop_shape(self) -> tuple[int, int]:
        return tuple(self._arr.shape[-2:])

    def __contains__(self, cell_id: str) -> bool:
        return str(cell_id) in self._index

    def get(self, cell_id: str) -> np.ndarray:
        """Return one crop, shape (channels, H, W)."""
        return np.asarray(self._arr[self._index[str(cell_id)]])

    def set(self, cell_id: str, crop: np.ndarray) -> None:
        """Store one crop. ``crop`` is (H, W) or (channels, H, W)."""
        if crop.ndim == 2:
            crop = crop[None, ...]
        self._arr[self._index[str(cell_id)]] = crop

    def get_many(self, cell_ids: Iterable[str]) -> np.ndarray:
        idx = [self._index[str(c)] for c in cell_ids]
        return np.asarray(self._arr[idx])


def _create_array(group, name, shape, chunks, dtype):
    """zarr v2/v3-compatible array creation."""
    try:  # zarr v3
        return group.create_array(name, shape=shape, chunks=chunks, dtype=dtype)
    except AttributeError:  # zarr v2
        return group.create_dataset(name, shape=shape, chunks=chunks, dtype=dtype)


# --------------------------------------------------------------------------- #
# Paired container
# --------------------------------------------------------------------------- #
class PairedData:
    """Pairs per-nucleus morphology + spatial position to multi-omic profiles."""

    MORPH_KEY = "morphology_features"  # uns key listing morphology obs columns

    def __init__(self, adata, crops: Optional[CropStore] = None, mods: Optional[dict] = None):
        self.adata = adata  # primary modality (RNA); holds shared obs/obsm
        self.crops = crops
        # additional modalities (e.g. ATAC peaks) sharing obs_names with adata
        self.mods: dict = mods or {}
        self._check_mods()

    def _check_mods(self) -> None:
        for name, m in self.mods.items():
            if list(m.obs_names) != list(self.adata.obs_names):
                raise ValueError(
                    f"modality '{name}' obs_names are not aligned to the primary "
                    "modality; every modality must share nuclei order/ids"
                )

    # -- convenience -------------------------------------------------------- #
    @property
    def n_cells(self) -> int:
        return self.adata.n_obs

    @property
    def cell_ids(self) -> list[str]:
        return list(self.adata.obs_names)

    @property
    def spatial(self) -> Optional[np.ndarray]:
        """Per-nucleus tissue (x, y) as an (n, 2) array. The canonical join key
        between a nucleus's multi-omic profile and its (future) image crop."""
        return self.adata.obsm.get("spatial")

    @property
    def has_spatial(self) -> bool:
        return self.spatial is not None

    @property
    def spatial_df(self) -> Optional[pd.DataFrame]:
        """Per-nucleus tissue coordinates as a DataFrame indexed by ``cell_id``."""
        s = self.spatial
        if s is None:
            return None
        return pd.DataFrame(s[:, :2], columns=["x", "y"], index=self.adata.obs_names)

    def set_spatial(self, coords, x="x", y="y") -> None:
        """Attach per-nucleus tissue coordinates from an (n, 2) array or a
        DataFrame/dict aligned to ``cell_id`` (used by the real slide-tags
        ingest, where coordinates come from the bead barcodes)."""
        if isinstance(coords, pd.DataFrame):
            coords = coords.reindex(self.adata.obs_names)
            arr = coords[[x, y]].to_numpy(dtype=float)
        else:
            arr = np.asarray(coords, dtype=float)
        if arr.shape != (self.n_cells, 2):
            raise ValueError(
                f"spatial coords must be ({self.n_cells}, 2), got {arr.shape}"
            )
        if np.isnan(arr).any():
            raise ValueError("spatial coords contain NaN (unaligned cell_ids?)")
        self.adata.obsm["spatial"] = arr

    @property
    def modalities(self) -> list[str]:
        mods = ["rna"] if self.adata.X is not None else []
        mods += [k for k in self.adata.layers.keys() if k is not None]
        mods += list(self.mods.keys())
        return mods

    def modality(self, name: str):
        """Return the AnnData for a modality ('rna' -> primary, else self.mods)."""
        if name == "rna":
            return self.adata
        return self.mods[name]

    @property
    def morphology_features(self) -> list[str]:
        return list(self.adata.uns.get(self.MORPH_KEY, []))

    @property
    def has_morphology(self) -> bool:
        return bool(self.morphology_features) and self.crops is not None

    def get_crop(self, cell_id: str) -> Optional[np.ndarray]:
        if self.crops is None or cell_id not in self.crops:
            return None
        return self.crops.get(cell_id)

    # -- catalogue ---------------------------------------------------------- #
    def catalogue(self, groupby: str = "cell_type") -> dict:
        """Fast summary of what the paired catalogue contains."""
        obs = self.adata.obs
        summary = {
            "n_nuclei": int(self.n_cells),
            "n_genes": int(self.adata.n_vars),
            "modality_sizes": {
                "rna": int(self.adata.n_vars),
                **{k: int(m.n_vars) for k, m in self.mods.items()},
            },
            "modalities": self.modalities,
            "has_spatial": self.has_spatial,
            "has_morphology": self.has_morphology,
            "morphology_features": self.morphology_features,
            "crop_shape": self.crops.crop_shape if self.crops else None,
        }
        if self.has_spatial:
            s = self.spatial[:, :2]
            summary["spatial_extent"] = {
                "x": (float(s[:, 0].min()), float(s[:, 0].max())),
                "y": (float(s[:, 1].min()), float(s[:, 1].max())),
            }
        if groupby in obs.columns:
            summary["per_group"] = (
                obs[groupby].value_counts().sort_index().to_dict()
            )
            summary["groupby"] = groupby
        return summary

    def catalogue_table(self, groupby: str = "cell_type") -> pd.DataFrame:
        """Per-group counts + mean morphology, as a tidy DataFrame."""
        obs = self.adata.obs
        cols = [groupby] if groupby in obs.columns else []
        feats = [f for f in self.morphology_features if f in obs.columns]
        if not cols:
            return pd.DataFrame({"n_nuclei": [self.n_cells]})
        agg = obs.groupby(groupby, observed=True).size().to_frame("n_nuclei")
        if feats:
            agg = agg.join(obs.groupby(groupby, observed=True)[feats].mean())
        return agg.reset_index()

    def __repr__(self) -> str:
        c = self.catalogue()
        return (
            f"PairedData: {c['n_nuclei']} nuclei x {c['n_genes']} genes | "
            f"modalities={c['modalities']} | spatial={c['has_spatial']} | "
            f"morphology={c['has_morphology']}"
        )

    # -- persistence -------------------------------------------------------- #
    def save(self, path: str) -> None:
        """Write to ``path/`` as ``adata.zarr`` + ``crops.zarr`` + manifest."""
        os.makedirs(path, exist_ok=True)
        self.adata.write_zarr(os.path.join(path, "adata.zarr"))
        for name, m in self.mods.items():
            m.write_zarr(os.path.join(path, f"mod_{name}.zarr"))
        manifest = {
            "format": "tessera",
            "version": 1,
            "has_crops": self.crops is not None,
            "mods": list(self.mods.keys()),
        }
        if self.crops is not None:
            src = self.crops._g.store
            # crops already persisted on disk if created with a path; re-copy to be safe
            _copy_crops(self.crops, os.path.join(path, "crops.zarr"))
        with open(os.path.join(path, "manifest.json"), "w") as fh:
            json.dump(manifest, fh, indent=2)

    @classmethod
    def load(cls, path: str) -> "PairedData":
        import anndata as ad

        adata = ad.read_zarr(os.path.join(path, "adata.zarr"))
        manifest = {}
        mpath = os.path.join(path, "manifest.json")
        if os.path.exists(mpath):
            with open(mpath) as fh:
                manifest = json.load(fh)
        mods = {}
        for name in manifest.get("mods", []):
            mods[name] = ad.read_zarr(os.path.join(path, f"mod_{name}.zarr"))
        crops = None
        crops_path = os.path.join(path, "crops.zarr")
        if os.path.exists(crops_path):
            crops = CropStore.open(crops_path, mode="r")
        return cls(adata, crops, mods=mods)


def _copy_crops(src: CropStore, dst_path: str) -> None:
    dst = CropStore.create(
        dst_path,
        cell_ids=src.cell_ids,
        shape=src.crop_shape,
        channels=src._arr.shape[1],
        channel_names=src.channel_names,
        dtype=str(src._arr.dtype),
    )
    dst._arr[:] = src._arr[:]
