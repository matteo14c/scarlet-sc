"""The report bundle: TSVs, figures (PNG / SVG) and a small JSON manifest; the HTML's single source.

Layout: <out>/report.html and <out>/bundle/{manifest.json, tables/*.tsv, figures/*.png|svg}.
"""
import colorsys
import datetime
import hashlib
import json
import os

import numpy as np
import pandas as pd

SCHEMA = "scarlet.report.bundle/1"
REFUSED_COLOUR = "#BDBDBD"
# colour-blind-safe: Okabe & Ito (2008) without black, then Paul Tol's muted and light schemes (light without its
# pale grey, too close to the refusal grey): 24 identities; beyond them evenly spaced hues, not colour-blind-safe
OKABE_ITO = ["#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7"]
TOL_MUTED = ["#CC6677", "#332288", "#DDCC77", "#117733", "#88CCEE", "#882255", "#44AA99", "#999933", "#AA4499"]
TOL_LIGHT = ["#77AADD", "#EE8866", "#EEDD88", "#FFAABB", "#99DDFF", "#44BB99", "#BBCC33", "#AAAA00"]
BASE = OKABE_ITO + TOL_MUTED + TOL_LIGHT
CVD_SAFE_MAX = len(BASE)


def colours(identities):
    """One fixed colour per identity (sorted by name), the same in every panel: the colour-blind-safe list first."""
    ids = sorted(map(str, identities))
    out = {}
    for i, k in enumerate(ids):
        if i < len(BASE):
            out[k] = BASE[i]
        else:
            h = ((i - len(BASE)) * 0.618033988749895) % 1.0
            r, g, b = colorsys.hsv_to_rgb(h, 0.65, 0.80)
            out[k] = "#%02X%02X%02X" % (int(r * 255), int(g * 255), int(b * 255))
    return out


def palette_note(n):
    """The page's line on the palette for n identities."""
    if n <= CVD_SAFE_MAX:
        return f"colour-blind-safe palette (Okabe-Ito, then Paul Tol's muted and light; {n} identities)"
    return (f"{n} identities: the first {CVD_SAFE_MAX} in the colour-blind-safe palette (Okabe-Ito, Paul Tol's muted "
            f"and light), the other {n - CVD_SAFE_MAX} in evenly spaced hues that are not colour-blind-safe")


def file_hash(path, limit_bytes=2 * 1024 ** 3):
    """sha256 of an input file (above 2 GB: size and mtime)"""
    st = os.stat(path)
    if st.st_size > limit_bytes:
        return dict(path=str(path), size=st.st_size, mtime=int(st.st_mtime), sha256=None,
                    note="above 2 GB: size and mtime, not hashed")
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 20), b""):
            h.update(b)
    return dict(path=str(path), size=st.st_size, sha256=h.hexdigest())


def checked_hash(path, recorded):
    """An input file's sha256 for a report: the hash the run recorded (learn, transfer) checked against the file the
    report reads (file_hash), with a line saying which; None without a path (an input in memory)."""
    if not path:
        return None
    f = file_hash(path)
    if not recorded:
        f["check"] = "computed by the report; the run recorded no hash"
    elif f.get("sha256") is None:
        f.update(sha256=recorded, check="recorded by the run; the file, above 2 GB, was not hashed again")
    elif f["sha256"] == recorded:
        f["check"] = "recorded by the run; the file the report read is the same"
    else:
        f["check"] = (f"the run recorded {recorded}; the file the report read differs: it is not the run's input")
    return f


class Bundle:
    def __init__(self, out, kind):
        self.out = str(out)
        self.dir = os.path.join(self.out, "bundle")
        os.makedirs(os.path.join(self.dir, "tables"), exist_ok=True)
        os.makedirs(os.path.join(self.dir, "figures"), exist_ok=True)
        self.manifest = dict(schema=SCHEMA, kind=kind, tables={}, figures={}, meta={},
                             written_at=datetime.datetime.now().astimezone().isoformat(timespec="seconds"))

    def table(self, name, df, about):
        p = os.path.join("tables", f"{name}.tsv")
        df.to_csv(os.path.join(self.dir, p), sep="\t", index=False)
        self.manifest["tables"][name] = dict(path=p, rows=int(len(df)), about=about)

    def meta(self, name, obj):
        self.manifest["meta"][name] = obj

    def figure(self, name, fig, about):
        """a matplotlib figure, saved as PNG and SVG"""
        paths = []
        for ext in ("png", "svg"):
            p = os.path.join("figures", f"{name}.{ext}")
            fig.savefig(os.path.join(self.dir, p), dpi=150 if ext == "png" else None, bbox_inches="tight")
            paths.append(p)
        self.manifest["figures"][name] = dict(paths=paths, about=about)

    def save(self):
        with open(os.path.join(self.dir, "manifest.json"), "w") as fh:
            json.dump(self.manifest, fh, indent=1, default=_jsonable)


def _jsonable(o):
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, (np.floating,)):
        return float(o)
    if isinstance(o, np.ndarray):
        return o.tolist()
    return str(o)


def load(out):
    """(manifest, {table name: DataFrame}) of a written bundle"""
    d = os.path.join(str(out), "bundle")
    man = json.load(open(os.path.join(d, "manifest.json")))
    tabs = {k: pd.read_csv(os.path.join(d, v["path"]), sep="\t", dtype={"cell": str}, keep_default_na=False,
                           na_values=[""]) for k, v in man["tables"].items()}
    return man, tabs
