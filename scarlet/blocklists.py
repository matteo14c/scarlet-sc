"""The blocklist option (auto, human, mouse, none or a file) and the identity of a blocklist file.

choose_blocklist() turns the option into the file a run uses and records how it was chosen. blocklist_identity()
records which file a run used (path, md5, species and, for a list SCARLET ships or shipped, its version label,
looked up by the file's md5 in KNOWN_BLOCKLISTS). The version label appears only in the run's provenance. Nothing
here changes which genes are filtered.
"""
from __future__ import annotations

import hashlib
import os
import re

KNOWN_BLOCKLISTS = {
    "c99d3c7e8b4b7418771505d2d8e86343": ("U1", "human"),
    "55a970417de47dbef310f88e124d4443": ("U2", "human"),
    "975e1c23e51d18901c5cd8dfe63effe9": ("U3", "human"),
    "b00199a00154ea4ab96be2896a9805f8": ("M1", "mouse"),
    "9a52b9b2b72337cae21c9e1a482cb8ba": ("U3.1", "human"),
    "bb1073b41dc9f19f460a40804356d1dc": ("M2", "mouse"),
    "610cb617c50dffa3fa3defac0c13e97f": ("M3", "mouse"),
    "9bc5fdbd9c8473ac50b75b6e9f1a230d": ("M4-draft", "mouse"),
    "68b329da9893e34099c7d8ad5cb9c940": ("E", "none"),
    "af4070d133cb6db821b7d0c6ebde14f8": ("E0", "none"),
    # the lists with gene groups (Category = group, scarlet/blocklist_modules.py)
    "3b7fddd9cb2dd8667dce9f9c63a6d938": ("U4", "human"),
    "d15c550f4fd8031325c5c97b7a01b745": ("M4", "mouse"),
}

_ENSG = re.compile(r"^ENSG\d{11}(\.\d+)?$")
_ENSMUSG = re.compile(r"^ENSMUSG\d{11}(\.\d+)?$")


def list_id_counts(blocklist_df) -> dict:
    """How many rows of a blocklist carry a human (ENSG), a mouse (ENSMUSG), another or no id."""
    if "ENSG_ID" not in blocklist_df.columns:
        return {"ensembl_human": 0, "ensembl_mouse": 0, "other": 0, "empty": int(len(blocklist_df))}
    ids = blocklist_df["ENSG_ID"].fillna("").astype(str).str.strip()
    hs = int(ids.str.match(_ENSG).sum())
    mm = int(ids.str.match(_ENSMUSG).sum())
    empty = int((ids == "").sum())
    return {"ensembl_human": hs, "ensembl_mouse": mm, "other": int(len(ids)) - hs - mm - empty, "empty": empty}


# Per-process memo of results that depend only on their inputs (a transfer learns many times on the same genes and
# the same blocklist file): the same values, computed once. Files are keyed by path, size and modification time.
_MEMO = {}


def _memo(key, compute):
    if key not in _MEMO:
        if len(_MEMO) > 64:
            _MEMO.clear()
        _MEMO[key] = compute()
    return _MEMO[key]


def _file_key(path):
    st = os.stat(path)
    return os.path.realpath(path), st.st_size, st.st_mtime_ns


def read_table(path):
    """pd.read_csv(path, sep='\\t') of a blocklist file, read once per process (a copy each time)."""
    import pandas as pd
    return _memo(("table",) + _file_key(path), lambda: pd.read_csv(path, sep="\t")).copy()


def file_md5(path) -> str:
    return _memo(("md5",) + _file_key(path), lambda: _file_md5(path))


def _file_md5(path) -> str:
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


# the blocklist option. 'auto' (the default) picks the bundled human or mouse list from the gene names;
# 'human', 'mouse' and 'none' are the bundled lists ('none' is a header-only file: nothing is removed); anything else
# is a file path, which must exist. The list's version is recorded in provenance only.
BLOCKLIST_CHOICES = ("auto", "human", "mouse", "none")
_PACKAGE_DIR = os.path.dirname(os.path.abspath(__file__))
BUNDLED_FILES = {"human": "scarlet_blocklist_universal.tsv", "mouse": "scarlet_blocklist_mouse.tsv",
                 "none": "scarlet_blocklist_none.tsv"}


class ChosenBlocklist(str):
    """The path of the blocklist a run uses, carrying how it was chosen (`record`: option, choice, basis), so the
    choice made at the start of a command reaches every record the command writes."""

    def __new__(cls, path, record):
        obj = super().__new__(cls, path)
        obj.record = dict(record)
        return obj

    def __reduce__(self):
        return ChosenBlocklist, (str(self), self.record)


def bundled_path(name):
    return os.path.join(_PACKAGE_DIR, BUNDLED_FILES[name])


def check_blocklist_option(value):
    """Refuse, before any work, a blocklist that is neither an option nor an existing file."""
    if value is None or isinstance(value, ChosenBlocklist) or str(value) in BLOCKLIST_CHOICES:
        return
    if not os.path.isfile(str(value)):
        raise FileNotFoundError(f"blocklist file {str(value)!r} not found. Give an existing file, or one of: "
                                f"{', '.join(BLOCKLIST_CHOICES)} (auto: chosen from the gene names)")


def choose_blocklist(value, gene_names=None):
    """The blocklist a run uses, as a ChosenBlocklist (a path with its record: option, choice, basis).

    value: None or 'auto' (the default: the bundled human or mouse list, from the gene names), 'human', 'mouse',
    'none' (nothing removed) or a path to an existing file. 'auto' reads the gene names: upper-case symbols (CD3E)
    and ENSG ids count as human, capitalised symbols (Cd3e) and ENSMUSG ids as mouse; at least 60 % mouse-style
    names decide mouse, at most 40 % human; otherwise, or with fewer than 20 usable names, it raises ValueError."""
    if isinstance(value, ChosenBlocklist):
        return value
    option = "auto" if value is None else str(value)
    if option == "auto":
        if gene_names is None:
            raise ValueError("blocklist='auto' needs the data's gene names")
        sp, frac, n = symbol_case_species(gene_names)
        if sp is None:
            how = (f"{n} usable gene names (at least 20 are needed)" if n < 20 else
                   f"{frac:.0%} of {n:,} gene names are mouse-style (capitalised, e.g. Cd3e) and the rest human-style "
                   f"(upper-case, e.g. CD3E); auto needs at least 60 % of one style")
            raise ValueError(f"blocklist='auto' cannot tell whether the data are human or mouse: {how}. Pass "
                             f"blocklist='human', 'mouse', 'none' or a file.")
        share = frac if sp == "mouse" else 1 - frac
        style = "capitalised symbols or ENSMUSG ids" if sp == "mouse" else "upper-case symbols or ENSG ids"
        rec = dict(option="auto", choice=sp,
                   basis=f"auto: {share:.1%} of {n:,} gene names are {sp}-style ({style}); the {sp} list")
        return ChosenBlocklist(bundled_path(sp), rec)
    if option in BUNDLED_FILES:
        basis = "given: nothing removed" if option == "none" else f"given: the {option} list"
        return ChosenBlocklist(bundled_path(option), dict(option=option, choice=option, basis=basis))
    check_blocklist_option(option)
    return ChosenBlocklist(option, dict(option="file", choice="file", basis=f"given: the file {option}"))


def blocklist_identity(path) -> dict:
    if path is None or str(path) == "" or not os.path.exists(str(path)):
        return _blocklist_identity(path)
    return dict(_memo(("identity", str(path)) + _file_key(str(path)), lambda: _blocklist_identity(path)))


def _blocklist_identity(path) -> dict:
    """{path, realpath, md5, version, species, n_rows} of a blocklist file (version/species None if unknown)."""
    if path is None or str(path) == "":
        return {"path": None, "md5": None, "version": None, "species": None, "n_rows": 0}
    p = str(path)
    out = {"path": p, "realpath": os.path.realpath(p)}
    if not os.path.exists(p):
        out.update(md5=None, version=None, species=None, n_rows=None, error="file not found")
        return out
    md5 = file_md5(p)
    version, species = KNOWN_BLOCKLISTS.get(md5, (None, None))
    with open(p) as fh:
        lines = [l for l in fh.read().splitlines() if l.strip()]
    n_rows = max(0, len(lines) - 1) if lines and lines[0].startswith("ENSG_ID") else len(lines)
    out.update(md5=md5, version=version, species=species, n_rows=n_rows)
    return out


def gene_id_type(names) -> str:
    names = [str(n) for n in names]
    return _memo(("id_type", tuple(names)), lambda: _gene_id_type(names))


def _gene_id_type(names) -> str:
    """'ensembl_human' | 'ensembl_mouse' | 'symbol' | 'mixed' (>= 90 % of one kind decides)."""
    names = [str(n).strip() for n in names]
    if not names:
        return "symbol"
    n = len(names)
    hs = sum(1 for g in names if _ENSG.match(g))
    mm = sum(1 for g in names if _ENSMUSG.match(g))
    if hs / n >= 0.9:
        return "ensembl_human"
    if mm / n >= 0.9:
        return "ensembl_mouse"
    if (hs + mm) / n <= 0.1:
        return "symbol"
    return "mixed"


def symbol_case_species(names, min_votes: int = 20):
    names = [str(n) for n in names]
    return _memo(("species", min_votes, tuple(names)), lambda: _symbol_case_species(names, min_votes))


def _symbol_case_species(names, min_votes: int = 20):
    """Species from the symbol case pattern: mouse 'Xkr4'/'Gm12345'/'0610005C13Rik', human 'XKR4'.

    Counts names with at least one letter; a name with a lower-case letter votes mouse, an all-upper
    name votes human (a '__chrN' suffix and the 'orf' of human open reading frames are ignored). Returns ('mouse'|'human'|None, fraction_mouse, n_votes). Ensembl ids decide by
    prefix (ENSG human, ENSMUSG mouse). None if fewer than ``min_votes`` usable names.
    """
    mouse = human = 0
    for g in names:
        g = str(g).strip()
        if _ENSMUSG.match(g):
            mouse += 1
            continue
        if _ENSG.match(g):
            human += 1
            continue
        # 'A1BG__chr19' (Muraro) -> 'A1BG'; human open reading frames 'C1orf112' are not mouse case
        g = g.split("__")[0].replace("orf", "")
        letters = [c for c in g if c.isalpha()]
        if not letters:
            continue
        # R-style mangling and clone names (RP11-..., AC012345.1, LINC...) are upper case in human data
        if any(c.islower() for c in letters):
            mouse += 1
        else:
            human += 1
    total = mouse + human
    if total < min_votes:
        return None, (mouse / total if total else None), total
    frac = mouse / total
    if frac >= 0.6:
        return "mouse", frac, total
    if frac <= 0.4:
        return "human", frac, total
    return None, frac, total
