"""The SCARLET model file: written by learn (and transfer), read by assign and load_model.

One zip archive (extension ``.scarlet`` by convention) holding, under fixed names:
  manifest.json             format name and version, the member list with sha256 and size
  markers.tsv               the ranked markers per cell type
  programmes.json           per cell type, its marker genes in rank order (derived from markers.tsv)
  universe.tsv              per marker gene, the detection rate of every cell type, plus the background column
  background_universe.tsv   the detection rate of every gene across the reference
  marker_finder_params.json every setting of the marker search as used
  settings.json             the settings learn resolved and its options
  provenance.json           versions, input, labels, time
  positive_evidence.json    per cell type, the expected and typical number of its markers a cell detects
  summary.txt               the run's plain-language summary
"""

from __future__ import annotations

import hashlib
import io
import json
import zipfile
from dataclasses import dataclass, field
from pathlib import Path

import pandas as pd

FORMAT = "scarlet.model"
FORMAT_VERSION = 1
MEMBERS = ("markers.tsv", "programmes.json", "universe.tsv", "background_universe.tsv",
           "marker_finder_params.json", "settings.json", "provenance.json", "in_class.tsv")
_DATE = (1980, 1, 1, 0, 0, 0)


@dataclass
class Model:
    """A learned SCARLET model, as learn returns it and load_model reads it from a model file.

    A model file is a zip archive of members (tables and records). ``files`` holds each member's bytes; the
    properties below read them. See [Reading the output](../output.md) for what each member holds.

    Attributes:
        files (dict): Member name -> bytes.
        path (str | None): The file the model was read from or saved to, if any.

    Examples:
        >>> import scarlet
        >>> m = scarlet.load_model("reference.scarlet")
        >>> m.identities
        >>> m.markers.head()
    """
    files: dict = field(default_factory=dict)
    path: str | None = None

    @property
    def markers(self) -> pd.DataFrame:
        """The markers of every cell type, ranked (markers.tsv)."""
        return pd.read_csv(io.BytesIO(self.files["markers.tsv"]), sep="\t")

    @property
    def universe(self) -> pd.DataFrame:
        """Per marker gene, its detection rate in every cell type and in the background (universe.tsv)."""
        return pd.read_csv(io.BytesIO(self.files["universe.tsv"]), sep="\t", index_col=0)

    @property
    def in_class(self) -> pd.DataFrame | None:
        """Every gene the marker search scored (after the blocklist) x each identity: the share of the identity's cells
        detecting it, the bootstrap estimate the universe holds for the markers, with the background column bg (the
        universe is this table's marker rows). None for a model without it (written before 0.6.3)."""
        if "in_class.tsv" not in self.files:
            return None
        return pd.read_csv(io.BytesIO(self.files["in_class.tsv"]), sep="\t", index_col=0)

    @property
    def background(self) -> pd.DataFrame | None:
        """The detection rate of every gene across the reference (background_universe.tsv), or None."""
        b = self.files.get("background_universe.tsv")
        return None if b is None else pd.read_csv(io.BytesIO(b), sep="\t")

    @property
    def programmes(self) -> dict:
        """Cell type -> its marker genes, in rank order (programmes.json)."""
        return json.loads(self.files["programmes.json"])

    @property
    def identities(self) -> list:
        """The model's cell types (its identities), in the order of the model."""
        return list(self.programmes)

    @property
    def settings(self) -> dict:
        """The command that learned the model, the settings it used and its options (settings.json)."""
        return json.loads(self.files["settings.json"])

    @property
    def provenance(self) -> dict:
        """When and with which versions the model was learned, and from which input (provenance.json)."""
        return json.loads(self.files["provenance.json"])

    @property
    def marker_finder_params(self) -> dict:
        """Every setting of the marker search as used (marker_finder_params.json)."""
        return json.loads(self.files["marker_finder_params.json"])

    def extract(self, member: str, dest) -> Path:
        """Write one member of the model file, byte for byte, to ``dest``; returns its path."""
        dest = Path(dest)
        dest.parent.mkdir(parents=True, exist_ok=True)
        dest.write_bytes(self.files[member])
        return dest

    def save(self, path) -> Path:
        """Write the model file to ``path`` (a zip archive with a manifest of checksums); returns its path."""
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        manifest = {"format": FORMAT, "format_version": FORMAT_VERSION,
                    "members": {n: {"sha256": hashlib.sha256(b).hexdigest(), "bytes": len(b)}
                                for n, b in sorted(self.files.items())}}
        with zipfile.ZipFile(path, "w", compression=zipfile.ZIP_DEFLATED) as z:
            for name, data in [("manifest.json", _json_bytes(manifest))] + sorted(self.files.items()):
                info = zipfile.ZipInfo(name, date_time=_DATE)
                info.compress_type = zipfile.ZIP_DEFLATED
                z.writestr(info, data)
        self.path = str(path)
        return path

    def __repr__(self):
        ids = self.identities
        s = self.settings.get("resolved", {})
        z = s["sizes"]
        sizes = f"{z['n_markers']} markers, shortlist {z['marker_shortlist']}"
        return (f"<SCARLET model: {len(ids)} identities, technology {s['technology']['name']}, {sizes}"
                f"{', ' + self.path if self.path else ''}>")


def load_model(path) -> Model:
    """Read a model file written by learn or transfer.

    Args:
        path (str): The model file (``.scarlet``).

    Returns:
        Model: The model.

    Raises:
        FileNotFoundError: The file does not exist.
        ValueError: The file is not a SCARLET model, was written by a newer version, or a member does not match
            its checksum (the file is damaged).

    Examples:
        >>> import scarlet
        >>> m = scarlet.load_model("reference.scarlet")
    """
    with zipfile.ZipFile(path) as z:
        manifest = json.loads(z.read("manifest.json"))
        if manifest.get("format") != FORMAT:
            raise ValueError(f"{path} is not a SCARLET model file (format {manifest.get('format')!r})")
        if int(manifest.get("format_version", 0)) > FORMAT_VERSION:
            raise ValueError(f"{path}: model format version {manifest['format_version']} is newer than this SCARLET "
                             f"reads ({FORMAT_VERSION})")
        files = {}
        for name, meta in manifest["members"].items():
            data = z.read(name)
            if hashlib.sha256(data).hexdigest() != meta["sha256"]:
                raise ValueError(f"{path}: member {name} does not match its checksum")
            files[name] = data
    return Model(files=files, path=str(path))


def from_outputs(markers, universe, params, settings, provenance, background=None, in_class=None) -> Model:
    """Build a model from the marker finder's output files (read as bytes, unchanged)."""
    files = {"markers.tsv": Path(markers).read_bytes(), "universe.tsv": Path(universe).read_bytes(),
             "marker_finder_params.json": Path(params).read_bytes(),
             "settings.json": _json_bytes(settings), "provenance.json": _json_bytes(provenance)}
    if background is not None and Path(background).exists():
        files["background_universe.tsv"] = Path(background).read_bytes()
    if in_class is not None and Path(in_class).exists():
        files["in_class.tsv"] = Path(in_class).read_bytes()
    files["programmes.json"] = _json_bytes(programmes_of(pd.read_csv(io.BytesIO(files["markers.tsv"]), sep="\t")))
    return Model(files=files)


def programmes_of(markers: pd.DataFrame) -> dict:
    """{identity: [marker genes in rank order]} from a marker table (columns cell_type, gene)."""
    out = {}
    for ident, gene in zip(markers["cell_type"].astype(str), markers["gene"].astype(str)):
        out.setdefault(ident, []).append(gene)
    return out


def _json_bytes(obj) -> bytes:
    return (json.dumps(obj, indent=2, default=_jsonable) + "\n").encode()


def _jsonable(o):
    import numpy as np
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, (np.floating,)):
        return float(o)
    if isinstance(o, (np.ndarray, pd.Index)):
        return list(o)
    if isinstance(o, Path):
        return str(o)
    return str(o)
