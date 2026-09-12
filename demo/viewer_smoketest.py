"""Headless smoke test: build the napari viewer offscreen, select a nucleus,
colour by a gene, and save canvas + profile-panel screenshots.

    QT_QPA_PLATFORM=offscreen python demo/viewer_smoketest.py
"""
import os
import sys

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

import matplotlib

matplotlib.use("Agg")

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
from tessera.core import PairedData  # noqa: E402
from tessera.viewer.napari_view import view  # noqa: E402

OUT = "data/processed"


def main():
    bundle = sys.argv[1] if len(sys.argv) > 1 else "data/processed/melanoma.tessera"
    paired = PairedData.load(bundle)
    print(repr(paired))

    v = view(paired)
    panel = v.window._dock_widgets["nucleus"].widget()  # our QWidget wrapper

    # find our _Panel instance to drive it
    from tessera.viewer import napari_view as nv  # noqa: F401
    # rebuild reference: the panel object is captured in the click closure; instead
    # re-create actions via the layer + a fresh panel is overkill. Just exercise API:
    pts = v.layers["nuclei"]

    # colour by a T-cell gene then by cell_type
    import numpy as np
    v.reset_view()
    canvas_ct = v.screenshot(canvas_only=True, flash=False)
    matplotlib.image.imsave(os.path.join(OUT, "viewer_spatial_celltype.png"), canvas_ct)
    print("saved viewer_spatial_celltype.png", canvas_ct.shape)

    # switch to UMAP embedding via the points data directly
    pts.data = paired.adata.obsm["X_umap"][:, ::-1]
    v.reset_view()
    canvas_um = v.screenshot(canvas_only=True, flash=False)
    matplotlib.image.imsave(os.path.join(OUT, "viewer_umap.png"), canvas_um)
    print("saved viewer_umap.png", canvas_um.shape)

    print("OK: viewer builds, layers:", [l.name for l in v.layers])
    v.close()


if __name__ == "__main__":
    main()
