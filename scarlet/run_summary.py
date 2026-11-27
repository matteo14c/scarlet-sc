"""The plain-language summary every run writes next to its outputs, and the report shows.

`summary()` turns what a command recorded (its inputs, the settings it resolved, its parameters and the package
versions) into a few lines a person can read without knowing SCARLET's records: what the command did, on which data,
which technology preset and blocklist it used and why, the main settings, the outputs and the versions.
"""
from __future__ import annotations

import time

FILE_NAME = "summary.txt"          # inside an output folder; next to a single output file: <stem>.summary.txt

WHAT = {
    "learn": "learned a model of the labelled cell identities of a reference",
    "assign": "annotated cells with a model",
    "transfer": "transferred the labels of a reference to a query",
}


def _what(command, params):
    """What the command did; a transfer says how (draws, refinement or not)."""
    w = WHAT.get(command, command)
    if command != "transfer":
        return w
    n = params.get("draws")
    how = ("a first pass, then a refinement on the query's own stably anchored cells" if params.get("refine") not in
           ("no", False) else "a first pass only, without refinement")
    if n:
        return f"{w} in {n} draws, each with its own seed ({how}); each cell's call is the consensus of the draws"
    return f"{w} ({how})"


def _num(x, fmt=",.0f"):
    return format(x, fmt) if isinstance(x, (int, float)) else str(x)


def _data_line(role, d):
    s = f"{role}: {d.get('path') or 'an AnnData object in memory'}"
    bits = []
    if d.get("n_cells") is not None:
        bits.append(f"{_num(d['n_cells'])} cells x {_num(d.get('n_genes'))} genes")
    if d.get("mean_detected_genes") is not None:
        bits.append(f"on average {_num(d['mean_detected_genes'])} genes detected per cell")
    if d.get("labels"):
        n = d.get("n_labels")
        bits.append(f"labels from column '{d['labels']}'" + (f" ({n} identities)" if n else ""))
    return s + (" (" + "; ".join(bits) + ")" if bits else "")


def technology_text(resolved, how=None):
    t = (resolved or {}).get("technology") or {}
    name = t.get("name")
    if name is None:
        return None
    kind = ("dense data (e.g. plate-based Smart-seq2): no background penalty in the marker search" if name == "dense"
            else "sparse data (e.g. droplet-based 10x): a background penalty in the marker search")
    kind += "; a gene counts as detected when its count is above 0"
    why = f"chosen automatically: {how}" if how else "given"
    return f"{name}, the preset for {kind}; {why}"


def blocklist_text(resolved):
    from scarlet.settings import blocklist_text as _bt
    b = (resolved or {}).get("blocklist")
    if not b:
        return None
    s = _bt(b)
    m = b.get("modules") or {}
    if m.get("modular_list") and m.get("off"):
        s += f"; gene groups left out for these labels: {', '.join(m['off'])}"
    return s


def cell_cycle_text(b):
    """Whether a learn kept the cell-cycle genes out of the markers, in words; None when the blocklist
    record has no module decision (a record made without the labels)."""
    if not b:
        return None
    if b.get("choice") == "none":
        return "not blocked (no blocklist applied)"
    m = b.get("modules")
    if not m:
        return None
    if not m.get("modular_list"):
        return "as the given list says (it has no gene groups)"
    if "cell_cycle" in (m.get("on") or []):
        return "blocked"
    labs = m.get("cell_cycle_off_by") or []
    if m.get("mode") == "auto" and labs:
        many = len(labs) > 1
        return (f"not blocked, because the label{'s' if many else ''} " + ", ".join(f"'{x}'" for x in labs)
                + f" name{'' if many else 's'} dividing cells")
    return f"not blocked (blocklist_modules {m.get('mode')})"


def _default_text(name, default):
    if default is None:
        return "the technology preset's" if name in ("n_markers", "marker_shortlist") else "none"
    return default


ADVANCED = ("n_markers", "marker_shortlist", "n_iterations", "n_cells", "blocklist_modules", "alpha")
BENCHMARK_NOTE = ("the published benchmarks used the defaults of the advanced options; with other values they do not "
                  "apply to this run")


def options_report_text(not_default):
    """The report's line on the run's options not at their defaults, with the benchmark note for an advanced one."""
    if not_default is None:
        return "not recorded"
    if not not_default:
        return "yes"
    s = "no: " + "; ".join(f"{k} {v['value']} (default: {_default_text(k, v['default'])})"
                           for k, v in not_default.items())
    return s + (f" ({BENCHMARK_NOTE})" if any(k in ADVANCED for k in not_default) else "")


def options_text(not_default):
    """The run's options that differ from their defaults (scarlet.api._options_not_default), in one line."""
    if not_default is None:
        return None
    if not not_default:
        return "All options at their defaults."
    return "Options not at their defaults: " + "; ".join(
        f"{k} {v['value']} (default: {_default_text(k, v['default'])})" for k, v in not_default.items()) + "."


def summary(command, *, inputs=(), resolved=None, params=None, versions=None, outputs=(), notes=(), when=None,
            not_default=None):
    """The summary text. inputs: [(role, {path, n_cells, n_genes, mean_detected_genes, labels, n_labels})];
    resolved: the command's settings record (scarlet.settings.resolve); params: its recorded parameters;
    versions: its recorded versions; outputs: [(path, what it is)]; notes: further plain lines; not_default: the
    run's options that differ from their defaults."""
    params, versions = params or {}, versions or {}
    when = when or time.strftime("%Y-%m-%d %H:%M")
    out = [f"SCARLET {versions.get('scarlet', '')} {command}, {when}".replace("  ", " "),
           f"What it did: {_what(command, params)}."]
    for role, d in inputs:
        out.append(_data_line(role, d))
    tr = params.get("technology_resolved") or {}
    tech = technology_text(resolved, tr.get("how") or params.get("technology_how"))
    if tech:
        out.append(f"Technology preset: {tech}.")
    bl = blocklist_text(resolved)
    if bl:
        out.append(f"Blocklist (genes left out of the markers): {bl}.")
    cc = cell_cycle_text((resolved or {}).get("blocklist"))
    if cc:
        out.append(f"Cell-cycle genes: {cc}.")
    sz = (resolved or {}).get("sizes") or {}
    if sz.get("n_markers") is not None:
        out.append(f"Markers: {sz['n_markers']} per cell type, ranked on a shortlist of {sz['marker_shortlist']}"
                   + (" (given)" if sz.get("explicit") else " (the technology preset's sizes)") + ".")
    im = (resolved or {}).get("imbalanced") or {}
    if im.get("handling") is not None:
        on = im["handling"] != "standard"
        out.append("Unequal identity sizes: " + ("corrected (the largest identity holds "
                   f"{im['f_max']:.0%} of the cells)" if on else "no correction needed (largest identity "
                   f"{im['f_max']:.0%} of the cells)") + ".")
    shown = dict(params)          # transfer records refine / depth matching by method: shown as the options' yes / no
    if shown.get("refine") is not None:
        shown["refine"] = "no" if shown["refine"] == "no" else "yes"
    if shown.get("depth_match") is not None:
        shown["depth_match"] = "no" if shown["depth_match"] == "none" else "yes"
    sets = [f"{k} {shown[k]}" for k in ("seed", "draws", "alpha", "refine", "depth_match", "flip")
            if shown.get(k) is not None]
    if sets:
        out.append("Settings: " + ", ".join(sets) + ".")
    opts = options_text(not_default)
    if opts:
        out.append(opts)
    out += list(notes)
    for path, what in outputs:
        out.append(f"Output: {path} ({what})")
    vs = [f"{k} {v}" for k, v in versions.items() if k not in ("scarlet", "git_commit") and v]
    if vs:
        out.append("Versions: SCARLET " + str(versions.get("scarlet")) + "; " + ", ".join(vs) + ".")
    return "\n".join(out) + "\n"


def write(path, text):
    with open(path, "w") as fh:
        fh.write(text)
    return path
