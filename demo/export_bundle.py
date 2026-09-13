"""Export a .tessera bundle to a single self-contained JSON file that the
drop-a-file web app (web/viewer_app.html) can open — no HTML re-export needed.

    python demo/export_bundle.py mydata.tessera mydata.tessera.json

The JSON holds everything the interactive viewer needs (coordinates, cell types,
depth, morphology, crop atlas, DE panel) except the full any-gene expression
matrix (that stays a napari/baked-page feature). Share the .json with a
collaborator; they drop it into the app.
"""
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, "src")
from export_web import build_payload  # noqa: E402
from tessera.core import PairedData  # noqa: E402


def main(bundle, out):
    p = PairedData.load(bundle)
    data = build_payload(p)
    payload = json.dumps(data, separators=(",", ":"))
    with open(out, "w") as fh:
        fh.write(payload)
    print(f"wrote {out}  ({len(payload.encode())/1e6:.2f} MB)")
    # optional .js wrapper (used as the app's built-in demo, loadable under CSP)
    if len(sys.argv) > 3 and sys.argv[3] == "--js":
        js = os.path.splitext(out)[0] + ".js"
        with open(js, "w") as fh:
            fh.write("window.DEMO_BUNDLE=" + payload + ";")
        print(f"wrote {js}")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
