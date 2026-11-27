"""SCARLET's Python interface: learn, assign, transfer and report, the same commands as the command line.

    learn(reference, labels, model)              a labelled reference -> a model file
    assign(model, data, output)                  a model and new cells -> the cells annotated
    transfer(reference, labels, query, output)   learn on the reference, annotate the query, refine on the query
    report(run)                                  the HTML report of a finished run
"""

from __future__ import annotations

import argparse
import functools
import inspect
import io
import json
import os
import platform
import shutil
import sys
import tempfile
import time
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd

from scarlet import _verbosity
from scarlet import run_summary as _run_summary
from scarlet import blocklists as _blocklists
from scarlet import settings as _settings
from scarlet.model import Model, from_outputs, load_model

__all__ = ["learn", "assign", "transfer", "report", "load_model", "Model", "TransferResult"]

# --------------------------------------------------------------------------------------------------------------
# The options' types and allowed values (scarlet/cli.py reads them). Each option is described in its function's
# docstring (the command line shows the same text); the defaults are the function signatures'.
# --------------------------------------------------------------------------------------------------------------
PARAMS = {
    "reference": (str,), "labels": (str,), "model": (str,), "data": (str,), "query": (str,), "output": (str,),
    "predictions": (str,), "technology": (str, ("auto", "sparse", "dense")), "blocklist": (str,),
    "blocklist_modules": (str,), "imbalanced": (str, ("auto", "yes", "no")), "gene_col": (str,),
    "n_markers": (int,), "marker_shortlist": (int,), "n_iterations": (int,), "n_cells": (int,), "seed": (int,),
    "ncpu": (int,), "report": (bool,), "alpha": (float,), "flip": (str, ("refuse", "flag")),
    "save_probs": (bool,), "chunk_size": (int,), "origin_col": (str,), "refine": (bool,), "depth_match": (bool,),
    "qc": (bool,), "draws": (int,), "seeds": (str,), "run": (str,), "self_assessment": (bool,), "n_folds": (int,), "umap_cells": (int,),
    "query_labels": (str,), "umap_all_genes": (bool,),
    "verbosity": (str, _verbosity.LEVELS),
}
# An option whose allowed values differ by command (scarlet/cli.py reads these first)
PARAMS_BY_COMMAND = {"transfer": {"flip": (str, ("flag", "refuse"))}}

# Settings of learn and transfer that are fixed: part of the method, not options.
_FIXED_SETTINGS = dict(exclude_blocklist_categories=None, include_blocklist_categories=None, universe="compute",
                       bootstrap_order="keyed", normalize_active_genes=False, active_genes_mean=1000,
                       active_genes_sd=250, min_genes_expressed=0, marker_alpha=1.5, marker_alpha_rule=None,
                       marker_alpha_rho=0.05, classifier_bg_mode=None, universe_mode="natural",
                       universe_temper=0.5, marker_beta=1.0, marker_gamma=None, marker_eps=0.0001)

_SETTING_NAMES = ("technology", "blocklist", "blocklist_modules", "exclude_blocklist_categories",
                  "include_blocklist_categories", "imbalanced", "gene_col", "universe", "n_iterations", "n_cells",
                  "n_markers", "marker_shortlist", "seed", "bootstrap_order", "normalize_active_genes", "active_genes_mean",
                  "active_genes_sd", "min_genes_expressed", "marker_alpha", "marker_alpha_rule", "marker_alpha_rho",
                  "classifier_bg_mode", "universe_mode", "universe_temper", "marker_beta", "marker_gamma",
                  "marker_eps", "ncpu", "marker_extra")


# --------------------------------------------------------------------------------------------------------------
# helpers
# --------------------------------------------------------------------------------------------------------------
def _versions():
    from importlib import metadata
    from scarlet import __version__
    out = {"python": platform.python_version(), "scarlet": __version__}
    for dist in ("numpy", "pandas", "scipy", "anndata", "scikit-learn", "statsmodels", "networkx", "scanpy"):
        try:
            out[dist] = metadata.version(dist)
        except metadata.PackageNotFoundError:
            out[dist] = None
    from scarlet.marker_finder import _git_commit
    out["git_commit"] = _git_commit()
    return out


def _provenance(command, params, inputs):
    return {"schema": "scarlet.provenance/1", "command": command, "created_at": time.strftime("%Y-%m-%dT%H:%M:%S%z"),
            "host": platform.node(), "argv": sys.argv, "versions": _versions(),
            "scarlet_file": str(Path(__file__).resolve().parent), "parameters": params, "inputs": inputs}


class _InputPath:
    """An input given as a path, or an AnnData written to a temporary .h5ad (so both take the same code path)."""

    def __init__(self, obj, name):
        self.tmp = None
        if isinstance(obj, (str, os.PathLike)):
            self.path = str(obj)
            if not Path(self.path).exists():
                raise FileNotFoundError(f"{name}: {self.path} does not exist")
        else:
            self.tmp = tempfile.mkdtemp(prefix="scarlet_input_")
            self.path = os.path.join(self.tmp, f"{name}.h5ad")
            obj.write_h5ad(self.path)

    def close(self):
        if self.tmp:
            shutil.rmtree(self.tmp, ignore_errors=True)


def _file_sha256(obj):
    """The sha256 of an input file (hex), or None for an in-memory AnnData."""
    if not isinstance(obj, (str, os.PathLike)):
        return None
    import hashlib
    h = hashlib.sha256()
    with open(obj, "rb") as fh:
        for blk in iter(lambda: fh.read(1 << 24), b""):
            h.update(blk)
    return h.hexdigest()


def _input_source(obj, sha256=None):
    """What a learn's in-memory input came from, for its records (the marker finder's input and --in): an input file's
    path and sha256 (computed unless given), or an in-memory AnnData and its shape."""
    if isinstance(obj, (str, os.PathLike)):
        return f"{obj} (sha256 {sha256 or _file_sha256(obj)})"
    return f"in-memory AnnData ({obj.n_obs} cells x {obj.n_vars} genes)"


def _blocklist_given(b):
    """The blocklist as the user gave it (auto, human, mouse, none or the file), for the recorded parameters."""
    rec = getattr(b, "record", None)
    return str(b) if rec is None or rec["option"] == "file" else rec["option"]


def _summary_path(output):
    """Where a run's plain-language summary goes: <output>/summary.txt for a folder, <stem>.summary.txt beside a
    file."""
    p = Path(output)
    return p / _run_summary.FILE_NAME if p.is_dir() else p.with_name(p.stem + ".summary.txt")


def _summarise_learn(m, model, params, data, labels, resolved):
    text = _run_summary.summary(
        "learn", inputs=[("Reference", dict(path=params["reference"], n_cells=data["n_cells"], n_genes=data["n_genes"],
                                            mean_detected_genes=data.get("mean_detected_genes"), labels=labels,
                                            n_labels=data.get("n_labels")))],
        resolved=resolved, params=params, versions=m.provenance.get("versions"),
        not_default=params.get("options_not_default"),
        outputs=[(str(model) if model is not None else "the returned Model",
                  "the model; use it with scarlet assign, or open its report with scarlet report")])
    m.files["summary.txt"] = text.encode()


def _check_setting_values(kw):
    _blocklists.check_blocklist_option(kw.get("blocklist"))   # a missing file stops before any work
    if kw.get("technology") not in _settings.TECHNOLOGY_CHOICES:
        raise ValueError(f"technology must be one of {_settings.TECHNOLOGY_CHOICES} (got {kw.get('technology')!r})")
    for name in ("imbalanced",):
        if name in kw and kw[name] is not None and kw[name] not in PARAMS[name][1]:
            raise ValueError(f"{name} must be one of {PARAMS[name][1]} (got {kw[name]!r})")


def _resolve(kw, labels=None, gene_names=None):
    """The settings a learn uses (scarlet.settings.resolve). kw['blocklist'] becomes the chosen list
    (scarlet.blocklists.choose_blocklist on gene_names), so every later step of the command uses, and records, the
    same list and how it was chosen."""
    if gene_names is not None or str(kw["blocklist"]) not in ("None", "auto"):
        kw["blocklist"] = _blocklists.choose_blocklist(kw["blocklist"], gene_names)
    return _settings.resolve(
        technology=kw["technology"], blocklist=kw["blocklist"],
        blocklist_modules=kw["blocklist_modules"], exclude_blocklist_categories=kw["exclude_blocklist_categories"],
        include_blocklist_categories=kw["include_blocklist_categories"], imbalanced=kw["imbalanced"],
        n_markers=kw["n_markers"], marker_shortlist=kw["marker_shortlist"], marker_gamma=kw["marker_gamma"],
        universe=kw["universe"], marker_alpha=kw["marker_alpha"],
        marker_alpha_rule=kw["marker_alpha_rule"], classifier_bg_mode=kw["classifier_bg_mode"], labels=labels,
        gene_names=gene_names)


def _marker_namespace(command, data_path, labels, kw, *, output_markers, output_universe, output_background):
    """The namespace the engine's marker search reads (scarlet._engine.run_find)."""
    compute = kw["universe"] == "compute"
    ns = dict(active_genes_mean=kw["active_genes_mean"], active_genes_sd=kw["active_genes_sd"],
              active_genes_seed=kw["seed"], background_genes_N=kw["marker_shortlist"], blocklist=kw["blocklist"],
              blocklist_modules=kw["blocklist_modules"], bootstrap_order=kw["bootstrap_order"],
              classifier_bg_mode=kw["classifier_bg_mode"], command="marker", compute_universe=compute,
              dry_run=False, exclude_blocklist_categories=kw["exclude_blocklist_categories"],
              external_universe=None if compute else kw["universe"], gene_col=kw["gene_col"], group_key=labels,
              imbalanced=kw["imbalanced"], include_blocklist_categories=kw["include_blocklist_categories"],
              input=str(data_path), marker_alpha=kw["marker_alpha"], marker_alpha_rho=kw["marker_alpha_rho"],
              marker_alpha_rule=kw["marker_alpha_rule"], marker_beta=kw["marker_beta"], marker_command=command,
              marker_eps=kw["marker_eps"], marker_extra=list(kw["marker_extra"] or []),
              marker_gamma=kw["marker_gamma"],
              min_genes_expressed=kw["min_genes_expressed"], n_cells=kw["n_cells"],
              n_iterations=kw["n_iterations"], ncpu=kw["ncpu"], normalize_active_genes=bool(kw["normalize_active_genes"]),
              output_background_universe=str(output_background), output_markers=str(output_markers),
              output_universe=str(output_universe), technology=kw["technology"],
              top_final=kw["n_markers"], universe_mode=kw["universe_mode"], universe_temper=kw["universe_temper"])
    return argparse.Namespace(**ns)


def _model_from_files(stem_dir, stem, settings, provenance, universe_name=None, background_name=None):
    stem_dir = Path(stem_dir)
    return from_outputs(stem_dir / f"{stem}.tsv", stem_dir / (universe_name or "universe.tsv"),
                        stem_dir / f"{stem}.params.json", settings, provenance,
                        background=stem_dir / (background_name or "background_universe.tsv"),
                        in_class=stem_dir / "in_class.tsv")


def _settings_record(command, resolved, params):
    return {"schema": "scarlet.model.settings/1", "command": command, "resolved": resolved, "parameters": params}


def _labels_frame(annotated):
    obs = annotated.obs
    uniq = obs["is_scarlet_unique"].astype(str).str.lower().eq("true").to_numpy()
    best = obs["scarlet_best"].astype(str).to_numpy()
    out = pd.DataFrame({"scarlet_best": best, "scarlet_unique": np.where(uniq, best, "unassigned"),
                        "is_scarlet_unique": uniq, "scarlet_fdr": obs["scarlet_fdr"].to_numpy()},
                       index=pd.Index(obs.index.astype(str), name="cell"))
    for c in ("scarlet_flip", "scarlet_reason", "scarlet_pe_low", "scarlet_pe_note"):   # the flags, when computed
        if c in obs.columns:
            out[c] = obs[c].to_numpy()
    return out


def _learn_into(model_path, data_path, labels, kw, resolved, command, params, inputs, quiet_preflight=False,
                data=None, scoring_stats=False, in_class=False):
    """learn's work: the marker finder (find) on `data` (the input in memory; read from `data_path` when None; the
    records name `data_path`), packed into a Model (saved to model_path if given).
    scoring_stats (learn and transfer): the marker finder also records each identity's n_eff (scarlet.scoring)."""
    from scarlet import _engine
    with tempfile.TemporaryDirectory(prefix="scarlet_learn_") as work:
        ns = _marker_namespace("find", data_path, labels, kw, output_markers=Path(work) / "markers.tsv",
                               output_universe=Path(work) / "universe.tsv",
                               output_background=Path(work) / "background_universe.tsv")
        if scoring_stats:
            ns.marker_extra = list(ns.marker_extra) + ["--scoring_stats"]
        if in_class:   # the full in-class table, for the models that are kept (learn's; transfer's first draw)
            ns._output_in_class = str(Path(work) / "in_class.tsv")
        with _settings_errors():
            _engine.run_find(ns, adata=data)
        m = _model_from_files(work, "markers", _settings_record(command, resolved, params),
                              _provenance(command, params, inputs))
    if model_path is not None:
        m.save(model_path)
    return m


class _settings_errors:
    """The marker finder validates its settings with argparse (parser.error -> SystemExit); API callers get a
    ValueError instead (the message is printed above it)."""

    def __enter__(self):
        return self

    def __exit__(self, typ, exc, tb):
        if typ is SystemExit and exc.code not in (0, None):
            raise ValueError(f"invalid settings: {exc.code if isinstance(exc.code, str) else 'see the message above'}") \
                from exc
        return False


def _add_positive_evidence(m, data, labels, kw, resolved, command, params, inputs):
    """The model's positive-evidence record (scarlet.evidence; the ratio rule's E_k and each identity's in-sample
    median of markers detected; no extra learning)."""
    from scarlet import evidence as _ev
    rec = _ev.reference(data, labels, m, gene_col=kw["gene_col"])
    m.files[_ev.MEMBER] = (json.dumps(rec, indent=1) + "\n").encode()
    return m


# Options that change no call: left out of a run's options_not_default
_OPERATIONAL = ("ncpu", "chunk_size", "report", "predictions", "save_probs", "verbosity")


def _json_value(v):
    return v if isinstance(v, (bool, int, float, str, type(None))) else str(v)


def _options_not_default(fn, given):
    """{option: {value, default}} for each option of `fn` (scarlet.api) whose value in `given` differs from its
    default; inputs, outputs and _OPERATIONAL left out. Recorded by every run, stated by its summary and report."""
    out = {}
    for k, p in inspect.signature(fn).parameters.items():
        if p.default is inspect.Parameter.empty or k in _OPERATIONAL:
            continue
        if given[k] != p.default:
            out[k] = {"value": _json_value(given[k]), "default": _json_value(p.default)}
    return out


def _kw(local):
    return {k: local[k] for k in _SETTING_NAMES}


def _command(fn):
    """learn, assign, transfer and report run inside scarlet._verbosity.command: the outermost sets what is printed
    (its verbosity option); the signature and docstring are fn's (the command line and the documentation read them)."""
    sig = inspect.signature(fn)

    @functools.wraps(fn)
    def run(*args, **kwargs):
        bound = sig.bind(*args, **kwargs)
        bound.apply_defaults()
        with _verbosity.command(fn.__name__, bound.arguments["verbosity"]):
            return fn(*args, **kwargs)
    return run


# --------------------------------------------------------------------------------------------------------------
# learn
# --------------------------------------------------------------------------------------------------------------
@_command
def learn(reference, labels, model, *, technology="auto", blocklist="auto", blocklist_modules="auto",
          imbalanced="auto", gene_col="index", n_markers=None, marker_shortlist=None, n_iterations=100, n_cells=100,
          seed=42, ncpu=1, report=False, verbosity="normal") -> Model:
    """Learn a model of the cell types of a labelled reference.

    For each cell type of the reference, SCARLET finds marker genes: genes consistently detected in that cell type
    and rarely elsewhere. It writes them, with the detection rates it needs to score new cells, to a model file, and
    writes a plain-language summary of the run next to it. The guide: [Learn](../guide/learn.md).

    Args:
        reference (str | AnnData): The labelled reference: an ``.h5ad`` file of raw counts (cells x genes); in
            Python also an AnnData object.
        labels (str): The column of the reference's ``obs`` that holds the cell-type labels to learn.
        model (str | None): The model file to write (``.scarlet``). In Python, ``None`` returns the model without
            writing it.
        technology (str): The technology preset: ``auto`` chooses ``dense`` when the reference detects at least
            5,000 genes per cell on average, else ``sparse``; ``sparse`` suits droplet data (e.g. 10x), ``dense``
            plate data (e.g. Smart-seq2).
        blocklist (str): Genes kept out of the markers (technical and stress genes): ``auto`` picks the bundled
            human or mouse list from the gene names; ``human``, ``mouse``, ``none``, or the path of your own
            tab-separated file with a ``Symbol`` and/or ``ENSG_ID`` column.
        blocklist_modules (str): Advanced. Which gene groups of the bundled lists to apply: ``auto`` (a fixed rule
            on the labels), ``all``, or a comma-separated list of groups. See [Advanced
            options](../advanced.md).
        imbalanced (str): Unequal class sizes: ``auto`` corrects when the largest cell type holds more than a
            third of the cells, ``yes`` always, ``no`` never. See [Advanced options](../advanced.md).
        gene_col (str): The ``var`` column with the gene names; ``index`` uses the var names.
        n_markers (int | None): Advanced. Markers per cell type; by default the technology preset's (50).
        marker_shortlist (int | None): Advanced. The shortlist the markers are ranked on: in each round of the
            marker search a gene counts towards its stability when it is among the cell type's top
            marker_shortlist genes; by default the technology preset's (150). Keep it larger than n_markers.
        n_iterations (int): Advanced. Rounds of the marker search.
        n_cells (int): Advanced. Cells drawn per cell type in each round of the marker search.
        seed (int): Random seed; the same inputs, options and seed give the same model.
        ncpu (int): Number of CPU cores to use; it does not change the result.
        report (bool): Also write the run's HTML report when the model is written (as ``scarlet report`` on it).
        verbosity (str): What the command prints: ``normal``, its pre-flight once, one progress bar and one line
            per fold or draw; ``quiet``, nothing; ``full``, everything each step prints. It does not change the
            result.

    Returns:
        Model: The learned model (also written to ``model`` when a path is given).

    Writes:
        ``model``, the model file; ``<model stem>.summary.txt``, the run's summary, next to it; with ``report``,
        the folder ``<model stem>_report/`` (``report.html`` and its data).

    Raises:
        FileNotFoundError: The reference, or a blocklist file given by path, does not exist.
        KeyError: ``labels`` is not a column of the reference's ``obs``.
        ValueError: An option has a value outside its documented ones; with ``blocklist='auto'``, the gene names
            do not tell whether the data are human or mouse; the reference holds negative values (e.g.
            batch-corrected or scaled data; SCARLET needs counts).

    Examples:
        >>> import scarlet
        >>> model = scarlet.learn("reference.h5ad", "cell_type", "reference.scarlet")

        The same from the command line:

            scarlet learn --reference reference.h5ad --labels cell_type --model reference.scarlet
    """
    _stop_if_unguarded_import()
    not_default = _options_not_default(learn, locals())
    if report:   # before any work: without plotly (the 'report' extra) the run goes on and the report is skipped
        report = _report_available()
    if technology not in ("auto",) + _settings.TECHNOLOGY_CHOICES:
        raise ValueError(f"technology must be auto, sparse or dense (got {technology!r})")
    exclude_blocklist_categories, include_blocklist_categories, universe, bootstrap_order, normalize_active_genes, \
        active_genes_mean, active_genes_sd, min_genes_expressed, marker_alpha, marker_alpha_rule, marker_alpha_rho, \
        classifier_bg_mode, universe_mode, universe_temper, marker_beta, marker_gamma, marker_eps = \
        _FIXED_SETTINGS.values()
    marker_extra = None
    if technology != "auto":
        _check_setting_values(_kw(locals()))
    sha_f = None
    if isinstance(reference, (str, os.PathLike)):
        if not Path(reference).exists():
            raise FileNotFoundError(f"reference: {reference} does not exist")
        from concurrent.futures import ThreadPoolExecutor
        from scarlet._read import read_h5ad
        hasher = ThreadPoolExecutor(max_workers=1)   # the file's sha256 while it is read (hashlib releases the GIL)
        sha_f = hasher.submit(_file_sha256, reference)
        hasher.shutdown(wait=False)
        R_in, ref_name = read_h5ad(reference), str(reference)
    else:   # used as given (not changed); the learn's records name it with its shape
        R_in, ref_name = reference, _input_source(reference)
    try:
        data = _settings.data_summary(R_in, labels, gene_col)
        tech_resolved = None
        if technology == "auto":   # the documented rule, as transfer (no plate label)
            d = data["mean_detected_genes"]
            technology = "dense" if d is not None and d >= _settings.HIGH_DENS_MIN_GENES else "sparse"
            tech_resolved = {"technology": technology,
                             "how": (f"resolved from the reference ({d:,.0f} mean detected genes per cell; dense "
                                     f"from {_settings.HIGH_DENS_MIN_GENES:,})" if d is not None else
                                     "resolved from the reference (no cells: sparse)")}
        kw = _kw(locals())
        _check_setting_values(kw)
        resolved = _resolve(kw, labels=data["labels"], gene_names=data["gene_names"])
        extra = ({"resolved": f"technology {tech_resolved['technology']} ({tech_resolved['how']})"}
                 if tech_resolved else None)
        _verbosity.say(_settings.preflight_text("learn", resolved, data, extra), _verbosity.at_top("learn"))
        params = dict(reference=str(reference) if isinstance(reference, (str, os.PathLike)) else "<AnnData>",
                      labels=labels, model=str(model), positive_evidence=True, **kw)
        params["blocklist"] = _blocklist_given(kw["blocklist"])
        params["options_not_default"] = not_default
        if tech_resolved:
            params["technology_resolved"] = tech_resolved
        inputs = {"reference": {"path": params["reference"], "sha256": sha_f.result() if sha_f else None,
                                "n_cells": data["n_cells"], "n_genes": data["n_genes"], "labels": labels,
                                "class_counts": _counts(data["labels"])}}
        from scarlet.marker_finder import binary_core_enabled, positive_pattern
        if isinstance(reference, (str, os.PathLike)) and binary_core_enabled() and not kw["normalize_active_genes"]:
            # memory: the input read here (not a caller's object) keeps its detection only once its values are
            # summarised (data, above): the marker search and the positive evidence read X > 0 (the binary core;
            # the marker search's own check, then its pattern of this pattern: the same matrix)
            from scarlet.input_checks import check_non_negative
            check_non_negative(R_in.X, "the marker-finder input")
            R_in.X = positive_pattern(R_in.X)
        m = _learn_into(None, ref_name, labels, kw, resolved, "learn", params, inputs, scoring_stats=True, data=R_in,
                        in_class=True)
        _add_positive_evidence(m, R_in, labels, kw, resolved, "learn", params, inputs)
        _summarise_learn(m, model, params, data, labels, resolved)
        if model is not None:
            m.save(model)
            _run_summary.write(_summary_path(model), m.files["summary.txt"].decode())
            _verbosity.say(f"learn finished: {model} ({len(m.identities)} cell types)", _verbosity.at_top("learn"))
    finally:
        del R_in
    if report:   # after the model is written; read-only
        if model is None:
            raise ValueError("learn(report=True) needs a model path")
        _report_or_warn(model, reference=reference, ncpu=ncpu)
    return m


def _report_available():
    """learn / transfer with report=True: whether the report can be written (plotly installed); if not, a warning, and the
    run goes on without its report (a missing report never stops a run)."""
    from scarlet.html_report.core import require_plotly
    try:
        require_plotly()
        return True
    except Exception as e:
        print(f"SCARLET warning: the report will not be written ({e}); the run goes on.", file=sys.stderr)
        return False


def _report_or_warn(run, **kw):
    """The report of a finished learn / transfer run (report=True). A report that fails is a warning: the run's outputs
    are complete and kept, the command ends normally, and the report can be written later with scarlet report."""
    try:
        globals()["report"](run, **kw)
    except Exception as e:
        print(f"SCARLET warning: the report of {run} failed ({type(e).__name__}: {e}); the run's outputs are complete. "
              f"Write the report later with: scarlet report {run}", file=sys.stderr)


def _counts(labels):
    return {str(k): int(v) for k, v in pd.Series(np.asarray(labels)).astype(str).value_counts().items()}


# --------------------------------------------------------------------------------------------------------------
# assign
# --------------------------------------------------------------------------------------------------------------
@_command
def assign(model, data, output, *, predictions=None, gene_col="index", alpha=0.01, flip="refuse", save_probs=True,
           chunk_size=10000, ncpu=1, verbosity="normal") -> pd.DataFrame:
    """Annotate cells with a model.

    Each cell is scored against every cell type of the model and against the background, and receives three calls:
    the most likely cell type (forced), the best call (``unassigned`` when the background explains the cell as
    well), and the unique call (``unassigned`` unless the best cell type is clearly better than its competitor). The
    guide: [Assign](../guide/assign.md); the columns: [Reading the output](../output.md).

    Args:
        model (str | Model): The model: a ``.scarlet`` file written by learn; in Python also a Model.
        data (str | AnnData): The cells to annotate: an ``.h5ad`` file of raw counts (cells x genes); in Python
            also an AnnData object.
        output (str | None): The annotated ``.h5ad`` to write (the input with the calls in ``obs``). In Python,
            ``None`` returns the calls without writing them.
        predictions (str | None): Also write the calls per cell as a tab-separated table to this file.
        gene_col (str): The ``var`` column with the gene names; ``index`` uses the var names.
        alpha (float): Advanced. The threshold of a unique call: a cell's best call is unique when its adjusted
            p-value (``scarlet_fdr``) is below alpha. See [Advanced options](../advanced.md).
        flip (str): The flip test, for calls that depend on genes the cell type never showed: ``refuse`` (the
            default) leaves such a cell unassigned, ``flag`` keeps the call and marks it in ``scarlet_flip``.
        save_probs (bool): Store each cell's probability for every cell type in the annotated file
            (``obsm['scarlet_probs']``).
        chunk_size (int): Cells scored at a time; lower it to use less memory. It does not change the result.
        ncpu (int): Number of CPU cores to use; it does not change the result.
        verbosity (str): What the command prints: ``normal``, its pre-flight once, one progress bar and one line
            per fold or draw; ``quiet``, nothing; ``full``, everything each step prints. It does not change the
            result.

    Returns:
        pandas.DataFrame: One row per cell: ``scarlet_best``, ``scarlet_unique``, ``is_scarlet_unique``,
            ``scarlet_fdr`` and the flags ``scarlet_flip``, ``scarlet_reason``, ``scarlet_pe_low``,
            ``scarlet_pe_note``.

    Writes:
        ``output``, the annotated ``.h5ad``; ``<output stem>.summary.txt``, the run's summary, next to it;
        ``predictions``, when given.

    Raises:
        FileNotFoundError: The model or the data file does not exist.
        KeyError: ``gene_col`` is not a column of the data's ``var``.
        ValueError: An option has a value outside its documented ones; the model file is damaged or was not
            written by this version of SCARLET; the data hold negative values (SCARLET needs counts).

    Examples:
        >>> import scarlet
        >>> calls = scarlet.assign("reference.scarlet", "query.h5ad", "query_annotated.h5ad")

        The same from the command line:

            scarlet assign --model reference.scarlet --data query.h5ad --output query_annotated.h5ad
    """
    _stop_if_unguarded_import()
    not_default = _options_not_default(assign, locals())
    if flip not in PARAMS["flip"][1]:
        raise ValueError(f"flip must be one of {PARAMS['flip'][1]} (got {flip!r})")
    return _assign(model, data, output, predictions=predictions, gene_col=gene_col, alpha=alpha, eps=0.00001,
                   ncpu=ncpu, chunk_size=chunk_size, save_probs=save_probs, save_probs_to_obs=False, flip=flip,
                   options_not_default=not_default)


def _assign(model, data, output, *, predictions=None, gene_col="index", alpha=0.01, eps=0.00001, ncpu=1,
            chunk_size=10000, save_probs=True, save_probs_to_obs=False, flip="refuse",
            options_not_default=None, ranks=False) -> pd.DataFrame:
    """assign's steps. options_not_default: assign's record of its options that differ from their defaults (written
    with the output).
    ranks=True (transfer): the returned table also holds each cell's scarlet_top1 / scarlet_top2."""
    import anndata as ad
    from scarlet import scoring as _sc
    m = model if isinstance(model, Model) else load_model(model)
    flip_mode = _sc.resolve_flip(flip)
    jf = _sc.jeffreys(m)
    if jf is None:
        raise ValueError("this model has no scoring statistics (each cell type's effective number of cells): it was "
                         "not learned by this version of SCARLET; learn it again")
    used = dict(flip=flip_mode, positive_evidence=True)
    from scarlet import evidence as _ev
    if not 0 < float(alpha) < 1:
        raise ValueError(f"alpha must be in (0, 1) (got {alpha})")
    if int(chunk_size) < 1 or int(ncpu) < 1:
        raise ValueError("chunk_size and ncpu must be >= 1")
    A = ad.read_h5ad(data) if isinstance(data, (str, os.PathLike)) else data
    universe_df = pd.read_csv(io.BytesIO(m.files["universe.tsv"]), sep="\t", index_col=0)
    genes = pd.Index(A.var_names if gene_col == "index" else A.var[gene_col].astype(str))
    resolved = m.settings["resolved"]
    detected = _settings.detected_per_cell(A.X)
    summary = {"n_cells": int(A.n_obs), "n_genes": int(A.n_vars), "mean_detected_genes": float(np.mean(detected))
               if A.n_obs else None}
    extra = {"model": f"{len(m.identities)} identities; {len(universe_df)} programme genes, "
                      f"{int(genes.isin(universe_df.index.astype(str)).sum())} found in the data"}
    seed = int((m.settings.get("parameters") or {}).get("seed", 42))   # the model's learn seed (learn's default 42)
    _verbosity.say(_settings.preflight_text("assign", resolved, summary, extra), _verbosity.at_top("assign"))
    from scarlet.classifier import classify
    from scarlet.marker_finder import binary_core_enabled, positive_pattern
    from scarlet.input_checks import check_non_negative
    # binary core: the classifier reads detection only; for an input read here from a file, the
    # values are checked and released, and the annotated output receives the input file's X unchanged
    binary_core = binary_core_enabled() and isinstance(data, (str, os.PathLike))
    if binary_core:   # the output's X is copied from the input file: the values read here are released
        check_non_negative(A.X, "the classifier input")
        A.X = positive_pattern(A.X)
    annotated = classify(A, universe_df, alpha=alpha, eps=eps, ncpu=ncpu, chunk_size=chunk_size, gene_col=gene_col,
                         no_save_probs=not save_probs, save_probs_to_obs=save_probs_to_obs,
                         floor=jf, flip=jf if used["flip"] != "off" else None,
                         flip_mode=used["flip"] if used["flip"] != "off" else "refuse")
    # the scoring's floor per cell type (scarlet.scoring.jeffreys)
    annotated.uns["scarlet_floor"] = {"identities": list(jf[0]), "f_k": np.array(list(jf[0].values())),
                                      "f_bg": float(jf[1])}
    top = annotated.obs["scarlet_top1"].astype(str).to_numpy()
    pe = _ev.flags(A, m, top, gene_col=gene_col)
    for c in pe.columns:
        annotated.obs[c] = pe[c].to_numpy()
    _verbosity.say(f"Positive evidence: {int(pe.scarlet_pe_low.sum())} of {len(pe)} calls flagged 'low positive "
                   "evidence' (label only)", _verbosity.at_top("assign"))
    annotated.uns["scarlet_assign"] = {k: (str(x) if not isinstance(x, bool) else x) for k, x in used.items()}
    if output is not None:
        # the run's alpha and its options not at their defaults (JSON: uns does not store None)
        annotated.uns["scarlet_options"] = {"alpha": float(alpha),
                                            "options_not_default": json.dumps(options_not_default or {})}
        mp = m.provenance.get("parameters") or {}
        notes = [f"Model: {getattr(m, 'path', None) or 'a Model in memory'} ({extra['model']})."]
        if extra.get("note"):
            notes.append(f"Note: {extra['note']}.")
        text = _run_summary.summary(
            "assign", inputs=[("Data", dict(path=str(data) if isinstance(data, (str, os.PathLike)) else None,
                                            **summary))],
            resolved=resolved, params=dict(alpha=alpha, seed=seed,
                                           technology_how=(mp.get("technology_resolved") or {}).get("how")),
            not_default=options_not_default,
            versions=_versions(), notes=notes,
            outputs=[(str(output), "the annotated data: the calls are in obs, see Reading the output")]
            + ([(str(predictions), "the calls per cell as a table")] if predictions is not None else []))
        annotated.uns["scarlet_summary"] = text
        Path(output).parent.mkdir(parents=True, exist_ok=True)
        print(f"Saving to {output}...")
        annotated.write_h5ad(output)
        if binary_core:
            _copy_X(data, output)
        _run_summary.write(_summary_path(output), text)
    lab = _labels_frame(annotated)
    if ranks:
        for c in ("scarlet_top1", "scarlet_top2"):
            lab[c] = annotated.obs[c].astype(str).to_numpy() if c in annotated.obs.columns else ""
    n_unique = int(lab.is_scarlet_unique.sum())
    _verbosity.say(f"Done. Unique calls: {n_unique}/{len(lab)} ({(n_unique / max(1, len(lab))) * 100:.1f}%)",
                   _verbosity.at_top("assign"))
    if predictions is not None:
        Path(predictions).parent.mkdir(parents=True, exist_ok=True)
        lab.to_csv(predictions, sep="\t")
    return lab


# --------------------------------------------------------------------------------------------------------------
# transfer
# --------------------------------------------------------------------------------------------------------------
def _copy_X(src_path, dst_path):
    """Replace the X of the written .h5ad dst with the X of src (copied by HDF5, not loaded into memory)."""
    import h5py
    with h5py.File(src_path, "r") as src, h5py.File(dst_path, "r+") as dst:
        del dst["X"]
        src.copy(src["X"], dst, name="X")


@dataclass
class TransferResult:
    """What transfer returns.

    Attributes:
        reference_models (dict): Origin -> the Model learned on the reference (one model, 'all', without depth
            matching).
        query_models (dict): Origin -> the Model learned on the query in the refinement.
        first_pass (pandas.DataFrame): Per query cell, the calls before refinement (first_pass.tsv).
        labels (pandas.DataFrame): Per query cell, the final calls (labels.tsv).
        qc (pandas.DataFrame | None): The transfer check per cell type (qc.tsv), or None without it.
        record (dict): What was done (transfer_record.json).
    """
    reference_models: dict
    query_models: dict
    first_pass: pd.DataFrame
    labels: pd.DataFrame
    qc: pd.DataFrame | None
    record: dict = field(default_factory=dict)

    @property
    def reference_model(self):
        """The first (or only) reference model."""
        return next(iter(self.reference_models.values()))


DRAWS = 10   # draws of a transfer by default (accurate mode); draws=1 is quick mode


def _parse_seeds(seeds):
    """The seed list: None, a list of integers, or (command line) a comma-separated string."""
    if seeds is None:
        return None
    if isinstance(seeds, str):
        parts = [x for x in seeds.replace(" ", "").split(",") if x]
        try:
            seeds = [int(x) for x in parts]
        except ValueError:
            raise ValueError(f"seeds must be integers separated by commas (got {seeds!r})") from None
    out = []
    for x in seeds:
        if isinstance(x, bool) or not float(x).is_integer():
            raise ValueError(f"seeds must be integers (got {x!r})")
        out.append(int(x))
    if len(set(out)) != len(out):
        raise ValueError("seeds must not repeat: each draw needs its own seed")
    if any(x < 0 for x in out):
        raise ValueError("seeds must be non-negative integers")
    return out


def _draw_plan(draws, seeds):
    """(N, seed list) from the options, checked."""
    seeds = _parse_seeds(seeds)
    if seeds is not None and draws is not None and int(draws) != len(seeds):
        raise ValueError(f"draws ({draws}) and the number of seeds ({len(seeds)}) differ; give one of them")
    n = len(seeds) if seeds is not None else (DRAWS if draws is None else int(draws))
    if n < 1:
        raise ValueError(f"draws must be at least 1 (got {n})")
    seeds = list(range(1, n + 1)) if seeds is None else seeds
    return n, seeds


@_command
def transfer(reference, labels, query, output, *, origin_col=None, refine=True, depth_match=True, qc=True,
             technology="auto", blocklist="auto", blocklist_modules="auto", imbalanced="auto", gene_col="index",
             n_markers=None, marker_shortlist=None, n_iterations=100, n_cells=100, draws=None, seeds=None,
             alpha=0.01, flip="flag", chunk_size=10000, ncpu=1, report=False, verbosity="normal") -> TransferResult:
    """Transfer the labels of a reference to a query.

    A first pass learns the reference and annotates the query; the refinement then learns again on the query's own
    stably anchored cells and annotates the query once more. Depth matching first thins the deeper of the two
    datasets towards the other. The whole transfer is repeated in draws (10 by default), each with its own seed for
    depth matching and the marker search; a cell's final call is the outcome most draws give, when at least 8 of the
    10 draws give it (otherwise the cell is left unassigned: undetermined across draws). With ``draws=1`` (quick
    mode) the transfer runs once and its calls are that draw's. The transfer check flags cell types whose reference
    markers the query does not confirm. The guide: [Transfer](../guide/transfer.md).

    Args:
        reference (str | AnnData): The labelled reference: an ``.h5ad`` file of raw counts (cells x genes); in
            Python also an AnnData object.
        labels (str): The column of the reference's ``obs`` that holds the cell-type labels.
        query (str | AnnData): The cells to annotate: an ``.h5ad`` file of raw counts; in Python also an AnnData
            object.
        output (str): The output folder; it is created if it does not exist.
        origin_col (str | None): A column of the query's ``obs`` naming where each cell comes from (a batch, a
            technology or a dataset; never the labels). Depth matching, refinement and the transfer check are done
            per origin.
        refine (bool): Refine the first-pass calls on the query's own stably anchored cells.
        depth_match (bool): Match the sequencing depth of reference and query before the first pass.
        qc (bool): Run the transfer check: per cell type, the share of its reference markers that the query
            confirms; below 0.2 the cell type is flagged ``doubtful transfer``.
        technology (str): The technology preset: ``auto`` chooses ``dense`` when the shallower of reference and
            query (after depth matching) detects at least 5,000 genes per cell on average, else ``sparse``.
        blocklist (str): Genes kept out of the markers: ``auto`` picks the bundled human or mouse list from the
            gene names; ``human``, ``mouse``, ``none``, or the path of your own tab-separated file.
        blocklist_modules (str): Advanced. Which gene groups of the bundled lists to apply: ``auto``, ``all``, or a
            comma-separated list of groups. See [Advanced options](../advanced.md).
        imbalanced (str): Unequal class sizes: ``auto`` corrects when the largest cell type holds more than a
            third of the cells, ``yes`` always, ``no`` never.
        gene_col (str): The ``var`` column with the gene names (in both datasets); ``index`` uses the var names.
        n_markers (int | None): Advanced. Markers per cell type; by default the technology preset's (50).
        marker_shortlist (int | None): Advanced. The shortlist the markers are ranked on; by default the
            technology preset's (150).
        n_iterations (int): Advanced. Rounds of the marker search, in each learn of the transfer.
        n_cells (int): Advanced. Cells drawn per cell type in each round of the marker search.
        draws (int | None): How many times the transfer is repeated, each time with its own seed: by default 10
            (accurate mode), or the number of ``seeds``; ``1`` is quick mode (one draw, no consensus). A final call
            needs its outcome in at least 80 % of the draws, rounded up (8 of 10). See [Advanced
            options](../advanced.md).
        seeds (list | str | None): Advanced. The seeds of the draws, one per draw (on the command line separated
            by commas); by default 1, 2, ... up to ``draws``. The same inputs, options and seeds give the same
            result.
        alpha (float): Advanced. The threshold of a unique call; in a transfer it also selects the stable anchors
            the refinement learns from. See [Advanced options](../advanced.md).
        flip (str): The flip test on the first pass: ``flag`` records it per cell and keeps the calls, ``refuse``
            refuses a flipped cell.
        chunk_size (int): Cells scored at a time; lower it to use less memory. It does not change the result.
        ncpu (int): Number of CPU cores to use: the draws run side by side, up to ``ncpu`` at a time (faster, with
            more memory: each draw holds its own copy of the data). It does not change the result.
        report (bool): Also write the run's HTML report when the transfer ends.
        verbosity (str): What the command prints: ``normal``, its pre-flight once, one progress bar and one line
            per fold or draw; ``quiet``, nothing; ``full``, everything each step prints. It does not change the
            result.

    Returns:
        TransferResult: ``labels`` (the final calls per query cell), ``first_pass`` (before refinement), ``qc``
            (the transfer check, or ``None``), ``reference_models`` and ``query_models`` (the first draw's models,
            by origin) and ``record`` (what was done).

    Writes:
        In ``output``: ``labels.tsv``, ``first_pass.tsv``, ``qc.tsv`` (with ``qc``), ``transfer_record.json``,
        ``summary.txt``, the first draw's models ``reference.scarlet`` (or ``reference_<origin>.scarlet``) and
        ``refine_<origin>.scarlet``; with ``report``, the folder ``report/``.

    Raises:
        FileNotFoundError: The reference or the query file, or a blocklist file given by path, does not exist.
        KeyError: ``labels`` is not a column of the reference's ``obs``, or ``origin_col`` of the query's.
        ValueError: An option has a value outside its documented ones; the species cannot be told from the gene
            names with ``blocklist='auto'``; an input holds negative values (SCARLET needs counts).

    Examples:
        >>> import scarlet
        >>> result = scarlet.transfer("reference.h5ad", "cell_type", "query.h5ad", "transfer_out")
        >>> result.labels.head()

        The same from the command line:

            scarlet transfer --reference reference.h5ad --labels cell_type --query query.h5ad --output transfer_out
    """
    _stop_if_unguarded_import()
    not_default = _options_not_default(transfer, locals())
    return _transfer(reference, labels, query, output, origin_col=origin_col, refine=refine,
                     depth_match=depth_match, qc=qc, technology=technology, blocklist=blocklist,
                     blocklist_modules=blocklist_modules, imbalanced=imbalanced, gene_col=gene_col,
                     n_markers=n_markers, marker_shortlist=marker_shortlist, n_iterations=n_iterations,
                     n_cells=n_cells, draws=draws, seeds=seeds, alpha=alpha, flip=flip,
                     chunk_size=chunk_size, ncpu=ncpu, report=report, not_default=not_default)


# ---------------------------------------------------------------------------------------------------------------
# transfer's draws. Each draw is one call of _first_pass_draw (and, with refinement, of _refine_draw) on the run's
# context; the draws run in this process or side by side in worker processes (ncpu), each worker holding its own
# copy of the context. A draw depends only on its own seeds, so the results do not depend on how many run at once.
# ---------------------------------------------------------------------------------------------------------------
_DRAW_CONTEXT = None
_FLAG_COLS = ("scarlet_flip", "scarlet_pe_low", "scarlet_pe_note")
_RANK_COLS = ("scarlet_top1", "scarlet_top2")


_POOL_DIR_ENV = "SCARLET_DRAW_POOL_DIR"
_GUARD_MESSAGE = ("transfer with ncpu above 1 runs its draws in new Python processes, which import your script: call "
                  "transfer under `if __name__ == \"__main__\":` (or use ncpu=1, which gives the same result)")


_WORKER_GUARD_MESSAGE = ("with ncpu above 1, SCARLET runs work in new Python processes, which import your script: "
                         "call SCARLET under `if __name__ == \"__main__\":` (or use ncpu=1, which gives the same result)")


def _stop_if_unguarded_import():
    """The one guard against a script without the __main__ guard; every verb calls it first (learn, assign, transfer,
    report; merge and discover when they join the package). A pool started with spawn imports the caller's script in
    each new worker, and a script without the guard then calls a verb again in that worker, while the worker is still
    inheriting its state from the caller (multiprocessing's _inheriting). Such a call stops here, before any work and
    before any pool of its own: it marks the transfer pool's folder for the caller (_run_draws then raises
    _GUARD_MESSAGE, so the caller exits non-zero) and ends the worker with a non-zero exit and the message."""
    import multiprocessing as mp
    if not getattr(mp.current_process(), "_inheriting", False):
        return
    folder = os.environ.get(_POOL_DIR_ENV)
    if folder:
        try:
            open(os.path.join(folder, "unguarded"), "w").close()
        except OSError:
            pass
    raise SystemExit(f"SCARLET: {_WORKER_GUARD_MESSAGE}")


def _set_draw_context(ctx):
    global _DRAW_CONTEXT
    _DRAW_CONTEXT = ctx


def _run_draws(fn, ctx, n_draws, workers, pool=None, extra=None, done=None):
    """[fn(d) for d in range(n_draws)], in draw order: here (workers 1, or one draw) or in the pool's `workers`
    processes. extra: entries added to the workers' context for these draws (here they are already in ctx).
    done(d, result): called as each draw's result arrives, in draw order (the progress line)."""
    def collect(results):
        out = []
        for d, r in enumerate(results):
            out.append(r)
            if done is not None:
                done(d, r)
        return out
    if pool is None or not pool.active:
        _set_draw_context(ctx)
        try:
            return collect(fn(d) for d in range(n_draws))
        finally:
            _set_draw_context(None)
    from concurrent.futures.process import BrokenProcessPool
    try:
        return collect(_draw_results(pool, fn, extra, n_draws))
    except BrokenProcessPool:
        if pool.unguarded():   # the workers re-ran the caller's unguarded script (_stop_if_unguarded_import)
            raise RuntimeError(_GUARD_MESSAGE) from None
        raise


def _draw_results(pool, fn, extra, n_draws):
    """The draws' results in draw order, as the pool's map gives them (all draws submitted at once). While it waits, it
    looks every second for the mark of a worker that stopped importing the caller's unguarded script
    (_stop_if_unguarded_import): the pool itself may never report that worker (it starts workers lazily, and a worker
    that dies right after starting can be missed, so the caller would wait forever). Then it shuts the pool down and raises _GUARD_MESSAGE."""
    from concurrent.futures import TimeoutError as FuturesTimeout
    ex = pool.executor()
    futures = [ex.submit(_draw_task, fn, extra, d) for d in range(n_draws)]
    for f in futures:
        while True:
            try:
                r = f.result(timeout=1.0)
                break
            except FuturesTimeout:
                if pool.unguarded():
                    pool.abort()
                    raise RuntimeError(_GUARD_MESSAGE) from None
        yield r


def _draw_task(fn, extra, d):
    if extra:
        _DRAW_CONTEXT.update(extra)
    return fn(d)


def _learn_matrix(X):
    """A thinned side as the learn and assign of a draw take it: the matrix itself when it is a detection pattern (they
    read X > 0); otherwise a float32 copy, as before."""
    import scipy.sparse as sp
    X = sp.csr_matrix(X)
    return X if X.dtype == np.bool_ else X.astype(np.float32)


class _SharedCSR:
    """A CSR matrix saved to .npy files in the run's temporary folder and opened read-only (memory-mapped) by each draw
    worker: the workers share its pages instead of each unpickling its own copy. The same arrays."""

    def __init__(self, X, folder, name):
        import scipy.sparse as sp
        X = sp.csr_matrix(X)
        self.shape, self.files = X.shape, {}
        for a in ("data", "indices", "indptr"):
            f = os.path.join(folder, f"{name}.{a}.npy")
            np.save(f, getattr(X, a))
            self.files[a] = f

    def open(self):
        import scipy.sparse as sp
        d, i, p = (np.load(self.files[a], mmap_mode="r") for a in ("data", "indices", "indptr"))
        return sp.csr_matrix((d, i, p), shape=self.shape, copy=False)


class _SharedAnnData:
    """An AnnData for the draw workers: everything but X pickled (as AnnData.copy copies it), X shared (_SharedCSR)."""

    def __init__(self, A, X_handle):
        X = A.X
        try:
            A.X = None
            self.meta = A.copy()
        finally:
            A.X = X
        self.X = X_handle

    def open(self, opened):
        A = self.meta.copy()
        A.X = opened(self.X)
        return A


def _shared_context(ctx, folder):
    """ctx with the reference's, the query's and the prepared depth matching's matrices shared (_SharedCSR), each
    matrix once."""
    import scipy.sparse as sp
    handles = {}

    def handle(X):
        if id(X) not in handles:
            handles[id(X)] = _SharedCSR(X, folder, f"m{len(handles)}")
        return handles[id(X)]
    out = dict(ctx)
    for k in ("Rk", "Q"):
        if sp.issparse(ctx[k].X):
            out[k] = _SharedAnnData(ctx[k], handle(ctx[k].X))
    if ctx.get("dm_prep"):
        out["dm_prep"] = {g: (dict(p, X=handle(p["X"])) if sp.issparse(p.get("X")) else p)
                          for g, p in ctx["dm_prep"].items()}
    return out


def _open_context(ctx):
    """In a draw worker: the shared matrices opened (memory-mapped, each once)."""
    opened = {}

    def op(h):
        if id(h) not in opened:
            opened[id(h)] = h.open()
        return opened[id(h)]
    c = dict(ctx)
    for k in ("Rk", "Q"):
        if isinstance(c.get(k), _SharedAnnData):
            c[k] = c[k].open(op)
    if c.get("dm_prep"):
        c["dm_prep"] = {g: (dict(p, X=op(p["X"])) if isinstance(p.get("X"), _SharedCSR) else p)
                        for g, p in c["dm_prep"].items()}
    return c


def _init_draw_worker(ctx):
    _verbosity.enter_worker()
    _set_draw_context(_open_context(ctx))


def _init_draw_worker_from_file(path, lock=None):
    """A draw worker's start: the shared context from the file _DrawPool.executor wrote (one pickle, so the matrices
    shared by several entries stay one), then as _init_draw_worker. The workers load it one at a time (lock): started
    together, they would otherwise all load at once, with their loading peaks coinciding and competing reads; the draws
    themselves then run in parallel."""
    import contextlib
    import pickle
    with (lock if lock is not None else contextlib.nullcontext()):
        with open(path, "rb") as f:
            ctx = pickle.load(f)
        _init_draw_worker(ctx)


class _DrawPool:
    """The draw workers of one transfer, started once for the first pass and the refinement (spawned processes, as
    before), with the run's read-only matrices shared through memory-mapped files of a temporary folder (removed by
    close). Inactive with one worker or one draw: the draws then run in this process."""

    def __init__(self, ctx, n_draws, workers):
        self.ctx, self.n = ctx, min(workers, n_draws)
        self.active = workers > 1 and n_draws > 1
        self._ex = self._dir = None
        self._env = None

    def executor(self):
        if self._ex is None:
            import multiprocessing as mp
            from concurrent.futures import ProcessPoolExecutor
            self._dir = tempfile.mkdtemp(prefix="scarlet_draws_")
            self._env = os.environ.get(_POOL_DIR_ENV)   # the workers inherit the folder (_stop_if_unguarded_import)
            os.environ[_POOL_DIR_ENV] = self._dir
            shared = _shared_context(self.ctx, self._dir)
            # the workers read the context from a file: passed as initargs it would be written into each new
            # worker's start-up pipe, and a worker that stops while importing the caller's unguarded script never
            # reads it, so a context larger than the pipe's buffer blocks that write, and the caller, forever
            # (multiprocessing keeps the pipe's reading end open until the write ends)
            import pickle
            context_file = os.path.join(self._dir, "context.pkl")
            with open(context_file, "wb") as f:
                pickle.dump(shared, f, protocol=pickle.HIGHEST_PROTOCOL)
            spawn = mp.get_context("spawn")
            self._ex = ProcessPoolExecutor(max_workers=self.n, mp_context=spawn, initializer=_init_draw_worker_from_file,
                                           initargs=(context_file, spawn.Lock()))
        return self._ex

    def unguarded(self):
        return self._dir is not None and os.path.exists(os.path.join(self._dir, "unguarded"))

    def abort(self):
        """Stop the pool at once: its pending draws cancelled, any live worker ended (close removes the folder)."""
        if self._ex is not None:
            procs = list((getattr(self._ex, "_processes", None) or {}).values())
            self._ex.shutdown(wait=False, cancel_futures=True)
            for p in procs:
                if p.is_alive():
                    p.terminate()

    def close(self):
        if self._ex is not None:
            self._ex.shutdown()
            self._ex = None
        if self._dir is not None:
            if self._env is None:
                os.environ.pop(_POOL_DIR_ENV, None)
            else:
                os.environ[_POOL_DIR_ENV] = self._env
        if self._dir is not None:
            shutil.rmtree(self._dir, ignore_errors=True)
            self._dir = None


def _learn_step(c, A, lab_col, tag, inputs, kw_d, source, in_class=False):
    """One learn of a transfer (reference, refinement or transfer check) on the cells A; source: what A is, for the
    records (_input_source and how the cells were used). Returns (Model, (tag, cell-cycle text))."""
    res = _resolve(kw_d, labels=A.obs[lab_col].to_numpy(), gene_names=A.var_names)
    prm = dict(c["params"], step=tag, seed=kw_d["seed"])
    m = _learn_into(None, source, lab_col, kw_d, res, "transfer", prm, inputs, scoring_stats=True, data=A,
                    in_class=in_class)
    _add_positive_evidence(m, A, lab_col, kw_d, res, "transfer", prm, inputs)
    return m, (tag, _run_summary.cell_cycle_text(res["blocklist"]))


def _first_pass_draw(d):
    """Draw d's depth matching and first pass: its calls, flags and ranks per query cell (in the run's cell order),
    its depth-matching plan, and (first draw) its reference models and learns."""
    import anndata as ad
    import scipy.sparse as sp
    from scarlet.depth_match import METHODS as DM_METHODS, thin_draw
    c = _DRAW_CONTEXT
    Rk, Q, origin, cells, dm = c["Rk"], c["Q"], c["origin"], c["cells"], c["dm"]
    dm_entry, ls = c["draw_seeds"][d]
    kw_d = dict(c["kw"], seed=int(ls), ncpu=c["inner_ncpu"])
    assign_kw = dict(c["assign_kw"], ncpu=c["inner_ncpu"])
    ent = [int(x) for x in np.atleast_1d(dm_entry)]
    qparts, ref_parts, plan_d, learns, thin_seed = {}, {}, [], [], {}
    if dm in DM_METHODS:   # the run's prepared depth matching (_transfer: prepare_depth), this draw's keep draw
        for gi, g in enumerate(c["origin_order"]):
            sel = np.flatnonzero(origin == g)
            seed_g = int(np.random.SeedSequence(ent + [gi]).generate_state(1)[0])
            Y, r = thin_draw(c["dm_prep"][g], seed_g)
            Xr2 = Y if r["thinned"] == "reference" else None
            Xq2 = Y if r["thinned"] == "query" else None
            plan = {"origin": str(g), "cells": int(len(sel)), "seed": seed_g, **r}
            plan_d.append(plan)
            thin_seed[g] = seed_g
            print(f"[SCARLET] draw {d + 1}: depth matching ({dm}), origin {g}: {plan}", flush=True)
            if r["thinned"] == "reference":
                ref_parts[g] = ad.AnnData(_learn_matrix(Xr2), obs=Rk.obs.copy(), var=Rk.var.copy())
            if r["thinned"] == "query":
                qparts[g] = ad.AnnData(_learn_matrix(Xq2), obs=Q.obs.iloc[sel].copy(),
                                       var=(Q.var if c["same_genes"] else Rk.var).copy())
    models, lab_all = {}, []
    if dm != "none":   # per origin (each origin's own assign: its own BH adjustment)
        for g in c["origin_order"]:
            sel = origin == g
            key = g if g in ref_parts else "all"
            if key not in models:
                Rg = ref_parts[g] if g in ref_parts else Rk
                how = (f"in memory, thinned, draw {d + 1}, seed {thin_seed[g]} (origin {g})" if g in ref_parts
                       else f"in memory, not thinned, draw {d + 1}")
                models[key], cc = _learn_step(c, Rg, c["labels"], f"reference ({key})",
                                              {"reference": c["params"]["reference"], "origin": key, "draw": d + 1},
                                              kw_d, f"{c['sources']['reference']}: {how}", in_class=d == 0)
                learns.append(cc)
            Qg = qparts[g] if g in qparts else (Q if sel.all() else Q[sel].copy())   # assign reads it
            lab_all.append(_assign(models[key], Qg, None, flip=c["flip"], **assign_kw))
    else:
        models["all"], cc = _learn_step(c, Rk, c["labels"], "reference",
                                        {"reference": c["params"]["reference"], "draw": d + 1}, kw_d,
                                        f"{c['sources']['reference']}: in memory, not thinned, draw {d + 1}",
                                        in_class=d == 0)
        learns.append(cc)
        lab_all.append(_assign(models["all"], Q, None, flip=c["flip"], **assign_kw))
    lab = pd.concat(lab_all)
    lab.index = lab.index.astype(str)
    lab = lab.reindex(cells)
    out = {"best": lab.scarlet_best.astype(str).to_numpy(), "unique": lab.scarlet_unique.astype(str).to_numpy(),
           "ranks": {k: lab[k].fillna("").astype(str).to_numpy() for k in _RANK_COLS},
           "flags": {k: (lab[k].to_numpy() if k in lab.columns else np.full(len(cells), None, dtype=object))
                     for k in _FLAG_COLS},
           "depth": plan_d, "models": models if d == 0 else None, "learns": learns if d == 0 else []}
    return out


def _refine_draw(d):
    """Draw d's refinement, per origin of the plan: learn on the origin's stable anchors of the refined cell types,
    score the origin's cells against the pooled background; cells whose first-pass best call is a cell type not
    refined keep that draw's first-pass calls. Returns the draw's refined columns (run's cell order), the cells that
    kept their first-pass calls, the backgrounds, and (first draw) the refinement models and learns."""
    from scarlet import transfer_steps as T
    c = _DRAW_CONTEXT
    Q, origin, cells, anchor = c["Q"], c["origin"], c["cells"], c["anchor"]
    ls = int(c["draw_seeds"][d][1])
    kw_d = dict(c["kw"], seed=ls, ncpu=c["inner_ncpu"])
    assign_kw = dict(c["assign_kw"], ncpu=c["inner_ncpu"])
    Fb, Fu = c["F_best"][:, d], c["F_uniq"][:, d]
    b_all, q_all = Fb.copy(), Fu.copy()
    pe = {k: c["F_flags"][k][:, d].copy() for k in ("scarlet_pe_low", "scarlet_pe_note")}
    ranks = {k: c["F_rank"][k][:, d].copy() for k in _RANK_COLS}
    kept = np.zeros(len(cells), dtype=bool)
    bgs, models, learns = [], {}, []
    for rec in c["plan"]:
        g = rec["origin"]
        pos = np.flatnonzero(origin == g)
        if "skipped" in rec:
            continue
        lab_arr = np.where(np.isin(anchor[pos].astype(str), rec["refined"]), anchor[pos], T.UNASSIGNED).astype(object)
        train = np.flatnonzero(lab_arr != T.UNASSIGNED)
        S = _rows(Q, pos)   # the refinement learns on the query's own cells, unthinned
        S.obs["scarlet_anchor_label"] = lab_arr
        inp_ = {"query": c["params"]["query"], "origin": g, "n_train": int(len(train)), "draw": d + 1}
        m, cc = _learn_step(c, S[train].copy(), "scarlet_anchor_label", f"refine ({g})", inp_, kw_d,
                            f"{c['sources']['query']}: in memory, not thinned, draw {d + 1}, refinement (origin {g})",
                            in_class=d == 0)
        learns.append(cc)
        m, bg = T.pooled_background(m, S)
        bgs.append(dict(draw=d + 1, origin=g, **bg))
        if d == 0:
            models[g] = m
        lab_g = _assign(m, S, None, flip="off", **assign_kw)
        lab_g.index = lab_g.index.astype(str)
        lab_g = lab_g.reindex(cells[pos])
        b = lab_g.scarlet_best.astype(str).to_numpy().astype(object)
        q = lab_g.scarlet_unique.astype(str).to_numpy().astype(object)
        keep = T.keep_first(Fb[pos], rec["refined"])
        b[keep] = Fb[pos][keep]
        q[keep] = Fu[pos][keep]
        kept[pos] = keep
        b_all[pos], q_all[pos] = b, q
        for k in pe:
            if k in lab_g.columns:
                v = lab_g[k].to_numpy().astype(object)
                v[keep] = pe[k][pos][keep]
                pe[k][pos] = v
        for k in ranks:   # the refinement model's ranks where the cell was refined in this draw
            v = lab_g[k].fillna("").astype(str).to_numpy().astype(object)
            v[keep] = ranks[k][pos][keep]
            ranks[k][pos] = v
    return {"best": b_all, "unique": q_all, "kept": kept, "pe": pe, "ranks": ranks, "backgrounds": bgs,
            "models": models if d == 0 else None, "learns": learns if d == 0 else []}


def _transfer(reference, labels, query, output, *, origin_col, refine, depth_match, qc, technology, blocklist,
              blocklist_modules, imbalanced, gene_col, n_markers, marker_shortlist, n_iterations, n_cells, draws,
              seeds, alpha, flip, chunk_size, ncpu, report, not_default, draw_seeds=None, anchors=None,
              keep_draws=False):
    """transfer's work. Internal arguments, for comparisons with recorded runs only (not options):
    draw_seeds: per draw (depth-matching seed entry, an int or a list of ints; learn seed) instead of (s, s);
    anchors: {cell: cell type}, the stable anchors given instead of computed; keep_draws: also write every draw's
    calls to draws.tsv.gz."""
    import anndata as ad
    import scipy.sparse as sp
    from scarlet import transfer_steps as T
    from scarlet.depth_match import METHODS as DM_METHODS

    if report:   # before any work: without plotly (the 'report' extra) the run goes on and the report is skipped
        report = _report_available()
    if technology not in ("auto",) + _settings.TECHNOLOGY_CHOICES:
        raise ValueError(f"technology must be auto, sparse or dense (got {technology!r})")
    if flip not in ("flag", "refuse"):
        raise ValueError(f"flip must be flag or refuse in transfer (got {flip!r})")
    if int(ncpu) < 1:
        raise ValueError(f"ncpu must be at least 1 (got {ncpu})")
    n_draws, seed_list = _draw_plan(draws, seeds)
    quick = n_draws == 1
    if draw_seeds is None:
        draw_seeds = [(s, s) for s in seed_list]
    if len(draw_seeds) != n_draws:
        raise ValueError("draw_seeds must have one entry per draw")
    workers = max(1, min(int(ncpu), n_draws))
    inner_ncpu = int(ncpu) if workers == 1 else 1   # draws side by side: one CPU each
    refine, dm = ("anchored" if refine else "no"), ("capture" if depth_match else "none")
    exclude_blocklist_categories, include_blocklist_categories, universe, bootstrap_order, normalize_active_genes, \
        active_genes_mean, active_genes_sd, min_genes_expressed, marker_alpha, marker_alpha_rule, marker_alpha_rho, \
        classifier_bg_mode, universe_mode, universe_temper, marker_beta, marker_gamma, marker_eps = \
        _FIXED_SETTINGS.values()
    marker_extra, eps = None, 0.00001
    seed = int(draw_seeds[0][1])   # the first draw's learn seed (kw's seed; each draw sets its own)
    out = Path(output)
    out.mkdir(parents=True, exist_ok=True)
    # the learns' records (_input_source): the files' checksums in a background thread (hashlib releases the GIL)
    # while the data are read and prepared
    from concurrent.futures import ThreadPoolExecutor
    hasher = ThreadPoolExecutor(max_workers=1)
    sources_f = hasher.submit(lambda: {"reference": _file_sha256(reference), "query": _file_sha256(query)})
    hasher.shutdown(wait=False)
    from scarlet._read import read_h5ad   # int64 counts held as int32 (the same values; scarlet._read)
    R = read_h5ad(reference) if isinstance(reference, (str, os.PathLike)) else reference.copy()
    Q = read_h5ad(query) if isinstance(query, (str, os.PathLike)) else query.copy()
    R.obs_names = R.obs_names.astype(str)
    Q.obs_names = Q.obs_names.astype(str)
    shapes = {"reference": R.shape, "query": Q.shape}   # as read, for the record's inputs
    if labels not in R.obs.columns:
        raise KeyError(f"label column {labels!r} is not an obs column of the reference")
    if origin_col is not None and origin_col not in Q.obs.columns:
        raise KeyError(f"origin column {origin_col!r} is not an obs column of the query")
    rd = _settings.data_summary(R, labels, gene_col)
    qd = _settings.data_summary(Q, None, gene_col)
    origin = (Q.obs[origin_col].astype(str).to_numpy() if origin_col is not None
              else np.array(["all"] * Q.n_obs, dtype=object))
    origins = sorted(set(origin))
    origin_order = list(dict.fromkeys(origin.tolist()))   # each origin's seed stream: its order of first appearance
    keep_ref = R.obs[labels].notna().to_numpy()
    Rk = R if keep_ref.all() else R[keep_ref].copy()  # no copy of a fully labelled reference
    del R
    cells = Q.obs_names.to_numpy().astype(str)
    n_q = len(cells)

    # the gene alignment and the counts check of depth matching (once); the technology from the inputs' depths
    same_genes, col = True, None
    if dm in DM_METHODS:
        T.check_counts(Rk.X, "reference")
        T.check_counts(Q.X, "query")
        rg = pd.Index(Rk.var_names if gene_col == "index" else Rk.var[gene_col].astype(str))
        qg = pd.Index(Q.var_names if gene_col == "index" else Q.var[gene_col].astype(str))
        same_genes = rg.equals(qg)
        col = None if same_genes else qg.get_indexer(rg)   # query column of each reference gene (-1: absent)
    tech = technology
    tech_note = "as given"
    if tech == "auto":   # from the shallower side (after matching, both sides have about its density)
        depths = {"reference": rd["mean_detected_genes"]}
        det_q = _settings.detected_per_cell(Q.X)   # per cell, as for each origin's rows (no copy of them)
        for g in origins:
            depths["query" if g == "all" else f"query origin {g}"] = float(np.mean(det_q[origin == g]))
        del det_q
        side, dmin = min(((k, x) for k, x in depths.items() if x is not None), key=lambda kv: kv[1])
        tech = "dense" if dmin >= _settings.HIGH_DENS_MIN_GENES else "sparse"
        tech_note = (f"resolved from the shallower side ({side}, {dmin:,.0f} mean detected genes per cell; dense "
                     f"from {_settings.HIGH_DENS_MIN_GENES:,})")
    kw = _kw(dict(locals(), technology=tech))
    _check_setting_values(kw)
    resolved = _resolve(kw, labels=rd["labels"], gene_names=rd["gene_names"])
    pre = _settings.transfer_preflight(rd, qd, resolved["technology"]["name"], depth_match=dm)
    pre["technology_resolved"] = {"technology": tech, "how": tech_note}
    k_need = T.needed(n_draws)
    mode_text = (f"quick mode: 1 draw (seed {seed_list[0]}; no consensus)" if quick else
                 f"{n_draws} draws (a final call needs {k_need} of {n_draws})")
    _verbosity.say(_settings.preflight_text("transfer", resolved, rd, {
        "query": f"{qd['n_cells']:,} cells x {qd['n_genes']:,} genes; mean {qd['mean_detected_genes']:,.0f} detected "
                 f"genes per cell; {len(origins)} origin(s)",
        "overlap": f"{pre['shared_genes']:,} shared genes ({pre['share_of_reference_genes']:.0%} of the reference, "
                   f"{pre['share_of_query_genes']:.0%} of the query)",
        "technology": f"{tech} ({tech_note})",
        "advice": "; ".join(pre["recommendation"]),
        "steps": f"depth matching {'yes' if dm != 'none' else 'no'}; flip test {flip}; refine "
                 f"{'yes' if refine != 'no' else 'no'}; transfer check {'yes' if qc else 'no'}; {mode_text}; "
                 f"{workers} draw(s) at a time"}), _verbosity.at_top("transfer"))
    params = dict(reference=str(reference) if isinstance(reference, (str, os.PathLike)) else "<AnnData>",
                  labels=labels, query=str(query) if isinstance(query, (str, os.PathLike)) else "<AnnData>",
                  output=str(output), origin_col=origin_col, refine=refine, depth_match=dm, draws=n_draws,
                  seeds=seed_list, qc=qc, alpha=alpha, eps=eps, chunk_size=chunk_size, flip=flip,
                  positive_evidence=True, **kw)
    params.pop("seed", None)   # a transfer's seeds are the draws' (params["seeds"])
    params["blocklist"] = _blocklist_given(kw["blocklist"])
    internal = {k: v for k, v in dict(draw_seeds=None if draw_seeds == [(s, s) for s in seed_list] else draw_seeds,
                                      anchors_given=anchors is not None).items()
                if v not in (None, False)}
    consensus_rec = None if quick else {
        "rule": "the most frequent outcome across the draws (unassigned counted as an outcome) is the final call "
                "when it holds at least 80 % of the draws, rounded up; otherwise unassigned (undetermined across "
                "draws)", "share": T.AGREE_SHARE, "needed": k_need}
    record = {"schema": "scarlet.transfer.record/4", "parameters": params, "options_not_default": not_default,
              "settings": resolved, "preflight": pre, "origins": origins,
              "draws": {"n": n_draws, "seeds": seed_list, "mode": "quick" if quick else "accurate",
                        "consensus": consensus_rec},
              "depth_match": [], "refine": None, "qc_performance": T.QC_PERFORMANCE, "versions": _versions()}
    if internal:
        record["internal"] = {k: (str(v) if k == "draw_seeds" else v) for k, v in internal.items()}
    assign_kw = dict(gene_col=gene_col, alpha=alpha, eps=eps, ncpu=ncpu, chunk_size=chunk_size, ranks=True)
    from scarlet.marker_finder import binary_core_enabled, positive_pattern
    if binary_core_enabled() and dm == "none":
        # binary core: without depth matching, every step reads the query's
        # detection only; its values are checked as the classifier checks them, after the reference as the
        # marker finder checks it, and released.
        from scarlet.input_checks import check_non_negative
        check_non_negative(Rk.X, "the marker-finder input")
        check_non_negative(Q.X, "the classifier input")
        Q.X = positive_pattern(Q.X)

    # depth matching's seed-independent part, once per origin (each draw then makes only its keep draw)
    dm_prep = None
    if dm in DM_METHODS:
        from scarlet.depth_match import prepare_depth
        dm_prep, dm_cache = {}, {}
        Xq_all = Q.X if sp.issparse(Q.X) else sp.csr_matrix(Q.X)
        for g in origin_order:
            sel = np.flatnonzero(origin == g)
            Xs_ = Xq_all if len(sel) == Xq_all.shape[0] else Xq_all[sel]   # the whole query: no copy
            dm_prep[g] = prepare_depth(Rk.X, Xs_ if same_genes else _aligned(Xs_, col), dm, dm_cache)
        del dm_cache, Xq_all, Xs_
        if binary_core_enabled():
            # detection only from here on: the values were read for the counts check, the summaries, the technology
            # choice and the deep side's log-means (above); every later step reads X > 0 (the binary core)
            pattern = {}
            for A_ in (Rk, Q):
                old_X = A_.X
                A_.X = positive_pattern(old_X)
                pattern[id(old_X)] = A_.X
            for prep in dm_prep.values():
                if "X" in prep:
                    prep["X"] = pattern[id(prep["X"])] if id(prep["X"]) in pattern else positive_pattern(prep["X"])
            del pattern, old_X

    # the first pass, draw by draw
    sha = sources_f.result()
    sources = {"reference": _input_source(reference, sha["reference"]), "query": _input_source(query, sha["query"])}
    record["inputs"] = {k: {"path": str(v) if isinstance(v, (str, os.PathLike)) else "<AnnData>", "sha256": sha[k],
                            "n_cells": int(shapes[k][0]), "n_genes": int(shapes[k][1])}
                        for k, v in (("reference", reference), ("query", query))}
    ctx = dict(Rk=Rk, Q=Q, origin=origin, origin_order=origin_order, cells=cells, dm=dm, same_genes=same_genes,
               col=col, draw_seeds=draw_seeds, kw=kw, params=params, labels=labels, flip=flip, assign_kw=assign_kw,
               inner_ncpu=inner_ncpu, dm_prep=dm_prep, sources=sources)
    pool = _DrawPool(ctx, n_draws, workers)   # one pool of draw workers for the first pass and the refinement
    passes = 2 if refine != "no" else 1
    progress = _verbosity.bar(n_draws * passes, "transfer")
    t_start = time.time()

    def draw_done(step):
        def line(d, r):
            n_called = int((np.asarray(r["unique"]).astype(str) != T.UNASSIGNED).sum())
            progress.update(1)
            _verbosity.say(f"{step}, draw {d + 1} of {n_draws}: {n_called:,} of {n_q:,} query cells called uniquely "
                           f"({time.time() - t_start:.0f} s)", _verbosity.at_top("transfer"))
        return line
    try:
        first_draws = _run_draws(_first_pass_draw, ctx, n_draws, workers, pool, done=draw_done("first pass"))
        F_best = np.empty((n_q, n_draws), dtype=object)
        F_uniq = np.empty((n_q, n_draws), dtype=object)
        F_flags = {k: np.empty((n_q, n_draws), dtype=object) for k in _FLAG_COLS}
        F_rank = {k: np.full((n_q, n_draws), "", dtype=object) for k in _RANK_COLS}
        for d, r in enumerate(first_draws):
            F_best[:, d], F_uniq[:, d] = r["best"], r["unique"]
            for k in _FLAG_COLS:
                F_flags[k][:, d] = r["flags"][k]
            for k in _RANK_COLS:
                F_rank[k][:, d] = r["ranks"][k]
            if dm in DM_METHODS:
                record["depth_match"].append({"draw": d + 1, "origins": r["depth"]})
        ref_models = first_draws[0]["models"]
        cell_cycle = list(first_draws[0]["learns"])   # (learn, whether it blocked the cell-cycle genes), first draw
        del first_draws

        def flipped_refused(d):
            if flip != "refuse":
                return np.zeros(n_q, bool)
            return pd.Series(F_flags["scarlet_flip"][:, d]).astype(str).str.lower().eq("true").to_numpy()

        # the stable anchors (every draw of the run) and the plan of the refinement
        anchor = np.full(n_q, T.UNASSIGNED, dtype=object)
        plan = []
        if refine != "no":
            if anchors is not None:
                anchor = np.array([str(anchors.get(cl, T.UNASSIGNED)) for cl in cells], dtype=object)
            else:
                anchor = T.stable_anchors(F_uniq)
            plan = T.refine_plan(anchor, origin, F_best, F_uniq)

        # the refinement, draw by draw
        R_best, R_uniq = F_best.copy(), F_uniq.copy()
        R_pe = {k: F_flags[k].copy() for k in ("scarlet_pe_low", "scarlet_pe_note")}
        R_rank = {k: M.copy() for k, M in F_rank.items()}
        kept_first = np.zeros((n_q, n_draws), dtype=bool)
        bg_records, query_models = [], {}
        if refine != "no":
            extra = dict(anchor=anchor, plan=plan, F_best=F_best, F_uniq=F_uniq, F_flags=F_flags, F_rank=F_rank)
            ctx.update(extra)
            for d, r in enumerate(_run_draws(_refine_draw, ctx, n_draws, workers, pool, extra,
                                             done=draw_done("refinement"))):
                R_best[:, d], R_uniq[:, d], kept_first[:, d] = r["best"], r["unique"], r["kept"]
                for k in R_pe:
                    R_pe[k][:, d] = r["pe"][k]
                for k in R_rank:
                    R_rank[k][:, d] = r["ranks"][k]
                bg_records += r["backgrounds"]
                if d == 0:
                    query_models, cell_cycle = r["models"], cell_cycle + r["learns"]
            # the stable anchors keep their anchor label (every stable anchor, in every origin)
            a = anchor != T.UNASSIGNED
            R_best[a, :] = anchor[a][:, None]
            R_uniq[a, :] = anchor[a][:, None]
        del ctx
    finally:
        pool.close()
        progress.close()
    for d in range(n_draws):   # the first pass's flip refusals carried through
        fr = flipped_refused(d)
        if fr.any():
            R_best[fr, d] = T.UNASSIGNED
            R_uniq[fr, d] = T.UNASSIGNED

    # the consensus of the draws (quick mode: the one draw's calls)
    def frame(B, Uq, flags, ranks):
        rb = np.empty(B.shape, dtype=object)
        ru = np.empty(B.shape, dtype=object)
        for d in range(n_draws):
            rb[:, d], ru[:, d] = T.draw_reasons(B[:, d], Uq[:, d], flipped_refused(d))
        lb, kb, why_b = T.consensus(B, rb)
        lu, ku, why_u = T.consensus(Uq, ru)
        df = pd.DataFrame(index=pd.Index(cells, name="cell"))
        df["origin"] = origin
        df["scarlet_best"] = lb
        df["scarlet_unique"] = lu
        if not quick:
            df["scarlet_best_agreement"] = kb
            df["scarlet_unique_agreement"] = ku
        df["scarlet_best_reason"] = why_b
        df["scarlet_unique_reason"] = why_u
        t1, run_up, n_run_up = T.runner_up(lb, ranks["scarlet_top1"], ranks["scarlet_top2"])
        df["scarlet_top1"] = t1
        df["scarlet_runner_up"] = run_up
        if not quick:
            df["scarlet_runner_up_draws"] = n_run_up
        for k, M in flags.items():
            if k == "scarlet_pe_note" or M is None or all(x is None for x in M[:, 0]):
                continue
            hit = pd.DataFrame(M).astype(str).apply(lambda s: s.str.lower().eq("true")).to_numpy()
            df[k] = hit.sum(axis=1) * 2 >= n_draws
            if not quick:
                df[k + "_draws"] = hit.sum(axis=1)
            if k == "scarlet_pe_low" and "scarlet_pe_note" in flags:
                notes = np.asarray(flags["scarlet_pe_note"], dtype=object)
                first_hit = np.where(hit.any(axis=1), hit.argmax(axis=1), -1)
                df["scarlet_pe_note"] = [str(notes[i, j]) if (ok and j >= 0) else ""
                                         for i, (j, ok) in enumerate(zip(first_hit, df[k].to_numpy()))]
        if refine != "no":
            df["scarlet_stable_anchor"] = anchor != T.UNASSIGNED
        return df

    first = frame(F_best, F_uniq, F_flags, F_rank)
    for key, m in ref_models.items():
        m.save(out / ("reference.scarlet" if key == "all" else f"reference_{_safe(key)}.scarlet"))
    first.to_csv(out / "first_pass.tsv", sep="\t")
    final = frame(R_best, R_uniq, dict(F_flags, **R_pe) if refine != "no" else F_flags, R_rank)
    if refine != "no":
        final["scarlet_refine_note"] = ""
        for rec in plan:
            pos = np.flatnonzero(origin == rec["origin"])
            if "skipped" in rec:
                final.iloc[pos, final.columns.get_loc("scarlet_refine_note")] = T.ORIGIN_NOT_REFINED_NOTE
                continue
            nr = kept_first[pos].sum(axis=1) * 2 > n_draws
            idx = pos[nr]
            final.iloc[idx, final.columns.get_loc("scarlet_refine_note")] = T.NOT_STABLE_NOTE
            rec["not_refined_cells"] = int(nr.sum())
        record["refine"] = {"decision": "refine where stably anchored (stable anchors: the cells called uniquely to "
                                        "the same cell type in every draw of the run)",
                            "anchors": "stable", "learns_on": "values",
                            "min_cells": T.MIN_CELLS, "rounds": 1, "background": "pooled: all the origin's query "
                            "cells", "anchors_keep_their_label": True, "origin_col": origin_col, "groups": plan,
                            "stable_anchor_cells": int((anchor != T.UNASSIGNED).sum()),
                            "backgrounds": bg_records}
        for g, m in query_models.items():
            m.save(out / f"refine_{_safe(g)}.scarlet")
        for rec in plan:
            print(f"[SCARLET] refine origin {rec['origin']}: {rec}", flush=True)
    final.to_csv(out / "labels.tsv", sep="\t")
    per_draw = []
    for d in range(n_draws):
        per_draw.append(dict(draw=d + 1, seed=seed_list[d],
                             first_best_labelled=int((F_best[:, d] != T.UNASSIGNED).sum()),
                             first_unique_labelled=int((F_uniq[:, d] != T.UNASSIGNED).sum()),
                             best_labelled=int((R_best[:, d] != T.UNASSIGNED).sum()),
                             unique_labelled=int((R_uniq[:, d] != T.UNASSIGNED).sum())))
    record["draws"]["per_draw"] = per_draw
    record["draws"]["final"] = {o: {str(k): int(v) for k, v in final[f"scarlet_{o}_reason"].replace(
        "", "labelled").value_counts().items()} for o in ("best", "unique")}
    if keep_draws:
        rows = []
        for d in range(n_draws):
            rows.append(pd.DataFrame({"cell": cells, "draw": d + 1, "seed": seed_list[d], "first_best": F_best[:, d],
                                      "first_unique": F_uniq[:, d], "best": R_best[:, d],
                                      "unique": R_uniq[:, d]}))
        pd.concat(rows).to_csv(out / "draws.tsv.gz", sep="\t", index=False)

    # the transfer check on the final labels (the first draw's reference models; the first seed)
    qc_df = None
    if qc:
        rows = []
        kw_q = dict(kw, seed=int(draw_seeds[0][1]))
        qctx = dict(params=params)
        for g in origins:
            pos = np.where(origin == g)[0]
            lab_g = final.scarlet_best.values[pos].astype(str)
            vc, ok = T.qc_classes(lab_g)
            mq = None
            if len(ok) >= 2:
                S = _rows(Q, pos)
                S.obs["scarlet_qc_label"] = lab_g
                St = S[np.isin(lab_g, ok)].copy()
                res = _resolve(kw_q, labels=St.obs["scarlet_qc_label"].to_numpy(), gene_names=St.var_names)
                cell_cycle.append((f"qc ({g})", _run_summary.cell_cycle_text(res["blocklist"])))
                mq = _learn_into(None, f"{sources['query']}: in memory, not thinned, transfer check (origin {g})",
                                 "scarlet_qc_label", kw_q, res, "transfer",
                                 dict(qctx["params"], step=f"qc ({g})", seed=kw_q["seed"]),
                                 {"query": params["query"], "origin": g}, data=St).markers
            ref_m = ref_models[g] if g in ref_models else ref_models["all"]
            rows += T.qc_rows(g, ref_m.markers, mq, lab_g)
        qc_df = pd.DataFrame(rows)
        qc_df.to_csv(out / "qc.tsv", sep="\t", index=False, float_format="%.4f")
        n_flag = int((qc_df.doubtful_transfer == True).sum())  # noqa: E712
        n_eval = int(qc_df.doubtful_transfer.notna().sum())
        record["qc"] = {"identities_evaluated": n_eval, "doubtful": n_flag}
        _verbosity.say(f"[SCARLET] transfer check: {n_flag} of {n_eval} cell types flagged 'doubtful transfer' (share "
                       f"of the reference markers confirmed on the query below {T.QC_FLAG}; how reliable the flag is: "
                       "the documentation, Reading the output)", _verbosity.at_top("transfer"))
    n_u = int((final.scarlet_unique.astype(str) != T.UNASSIGNED).sum())
    n_b = int((final.scarlet_best.astype(str) != T.UNASSIGNED).sum())
    why_u = final.scarlet_unique_reason.value_counts()
    cc_ref = _run_summary.cell_cycle_text(resolved["blocklist"])
    record["cell_cycle"] = {"reference": cc_ref, "learns": [{"learn": t, "cell_cycle_genes": cc} for t, cc in cell_cycle]}
    others = sorted({f"{t}: {cc}" for t, cc in cell_cycle if cc != cc_ref})
    if quick:
        notes = [f"Quick mode: one draw (seed {seed_list[0]}); the calls are that draw's (no consensus, no per-cell "
                 "agreement); another seed can change some of them."]
    else:
        notes = [f"Draws: {n_draws} (seeds {_seed_text(seed_list)}); a cell's final call is the outcome most draws "
                 f"give, when at least {k_need} of the {n_draws} give it; otherwise it is left unassigned "
                 "(undetermined across draws)."]
    if refine != "no":
        n_ref = sum(len(r["refined"]) for r in plan)
        n_nr = sum(len(r["not_refined"]) for r in plan if "skipped" not in r)
        skipped = [r["origin"] for r in plan if "skipped" in r]
        notes.append(f"Refinement: {int((anchor != T.UNASSIGNED).sum()):,} stable anchors (cells called uniquely to the "
                     + ("same cell type in the draw" if quick else f"same cell type in each of the {n_draws} draws")
                     + f"); {n_ref} cell type(s) refined and {n_nr} not refined (not anchored stably), over the origins"
                     + (f"; origin(s) not refined (fewer than 2 cell types anchored stably): {', '.join(skipped)}"
                        if skipped else "") + ".")
    if others:
        notes.append(f"Cell-cycle genes in the other learns of this transfer, where they differ: {'; '.join(others)}.")
    notes.append(f"Result: {n_u:,} of {len(final):,} query cells received a unique call and {n_b:,} a best call; "
                 f"unique calls refused: {int(why_u.get(T.REASON_BACKGROUND, 0)):,} background, "
                 f"{int(why_u.get(T.REASON_AMBIGUITY, 0)):,} ambiguity"
                 + ("" if quick else f", {int(why_u.get(T.UNDETERMINED, 0)):,} undetermined across draws")
                 + (f", {int(why_u.get(T.REASON_FLIP, 0)):,} flip test" if flip == "refuse" else "") + ".")
    text = _run_summary.summary(
        "transfer",
        inputs=[("Reference", dict(path=params["reference"], n_cells=rd["n_cells"], n_genes=rd["n_genes"],
                                   mean_detected_genes=rd.get("mean_detected_genes"), labels=labels,
                                   n_labels=rd.get("n_labels"))),
                ("Query", dict(path=params["query"], n_cells=qd["n_cells"], n_genes=qd["n_genes"],
                               mean_detected_genes=qd.get("mean_detected_genes")))],
        resolved=resolved, versions=record["versions"],
        params=dict(params, technology_how=None if tech_note == "as given" else tech_note),
        not_default=not_default, notes=notes,
        outputs=[(str(out / "labels.tsv"), "the final calls per query cell"),
                 (str(out / "first_pass.tsv"), "the first-pass calls, before refinement"),
                 (str(out / "transfer_record.json"), "the full record of this run")]
        + ([(str(out / "qc.tsv"), "the transfer check per identity")] if qc else []))
    record["summary"] = text
    with open(out / "transfer_record.json", "w") as fh:
        json.dump(record, fh, indent=2, default=str)
        fh.write("\n")
    _run_summary.write(out / _run_summary.FILE_NAME, text)
    _verbosity.say(f"transfer finished: {out}. {notes[-1]}", _verbosity.at_top("transfer"))
    if report:   # after every output is written; read-only; in-memory inputs are passed to the report as given
        _report_or_warn(str(out), reference=None if isinstance(reference, (str, os.PathLike)) else reference,
                        query=None if isinstance(query, (str, os.PathLike)) else query, ncpu=ncpu)
    return TransferResult(ref_models, query_models, first, final, qc_df, record)


def _rows(A, pos):
    """A[pos].copy(), but sharing A's X (read only by its callers) when pos is every row in order: everything else
    is copied as AnnData.copy copies it."""
    if len(pos) != A.n_obs or not np.array_equal(pos, np.arange(A.n_obs)):
        return A[pos].copy()
    X = A.X
    try:
        A.X = None
        S = A.copy()
    finally:
        A.X = X
    S.X = X
    return S


def _seed_text(seeds):
    """1, 2, ..., 20 as '1-20'; other lists as written."""
    if seeds == list(range(seeds[0], seeds[0] + len(seeds))):
        return f"{seeds[0]}-{seeds[-1]}" if len(seeds) > 1 else str(seeds[0])
    return ", ".join(str(s) for s in seeds)


def _aligned(X, col):
    """X (cells x query genes, CSR) with its columns put in the reference's gene order `col` (query column per
    reference gene, -1 = absent: an empty column)."""
    import scipy.sparse as sp
    X = sp.csc_matrix(X)
    present = np.flatnonzero(col >= 0)
    sub = X[:, col[present]]
    P = sp.csr_matrix((np.ones(len(present)), (np.arange(len(present)), present)), shape=(len(present), len(col)))
    return sp.csr_matrix(sub @ P)


def _safe(name):
    return "".join(c if (c.isalnum() or c in "-_.") else "_" for c in str(name))


@_command
def report(run, *, reference=None, query=None, model=None, query_labels=None, output=None, self_assessment=True,
           n_folds=5, seed=0, umap_cells=50000, umap_all_genes=True, alpha=0.01, ncpu=1, verbosity="normal") -> str:
    """Write the HTML report of a finished run.

    The report reads a run's outputs and changes nothing: a learn report (from a model file) shows the reference,
    the markers and a cross-validated self-assessment; a transfer report (from a transfer folder) shows the first
    pass and the refined calls side by side and the transfer check; an assign report (from an annotated ``.h5ad``,
    with its model) shows the calls and the refusals. Each report also shows the run's summary. The guide:
    [Reports](../reports/index.md).

    Args:
        run (str): What to report on: a model file (learn), a transfer output folder, or an annotated ``.h5ad``
            (assign).
        reference (str | None): The reference a learn or transfer run used, when its recorded path is no longer
            readable.
        query (str | None): The query a transfer run used, when its recorded path is no longer readable.
        model (str | None): The model an assign run used (needed for an assign report).
        query_labels (str | None): Transfer and assign reports: an obs column of the query (transfer) or of the
            annotated data (assign) holding labels you already know for its cells, e.g. in a tutorial. The report
            then adds a section comparing SCARLET's calls with them (share refused, F1 on the called cells, per
            label), marked as labels supplied by you; SCARLET never reads them.
        output (str | None): The report folder; by default next to the run (``<stem>_report/``, or ``report/``
            in a transfer folder).
        self_assessment (bool): Learn and transfer reports: estimate how well the model calls its own reference,
            by cross-validation (it learns once per fold, so it takes longer).
        n_folds (int): Folds of the self-assessment.
        seed (int): Seed of the self-assessment's folds and of the UMAP layout.
        umap_cells (int): At most this many cells in the report's UMAP plots.
        umap_all_genes (bool): Also compute each UMAP on all genes, a second space the reader can switch to (the
            place to look for populations the model lacks); the default space is the model's own genes. No skips it
            (faster on very large data); the report then says so.
        alpha (float): Advanced. The threshold of a unique call in the learn report's self-assessment. An assign
            report uses the alpha the assign run recorded (a different value given here is an error).
        ncpu (int): Number of CPU cores to use; it does not change the result.
        verbosity (str): What the command prints: ``normal``, its pre-flight once, one progress bar and one line
            per fold or draw; ``quiet``, nothing; ``full``, everything each step prints. It does not change the
            result.

    Returns:
        str: The path of ``report.html``.

    Writes:
        The report folder: ``report.html`` (one self-contained page) and ``bundle/`` (its tables and data); a
        transfer report also writes ``refinement.html`` (what refinement changed, linked from its section 6) and
        ``supplement_A/`` (the report of the reference model).

    Raises:
        ValueError: ``run`` is not a finished SCARLET run; an input the report needs (the reference or the query)
            cannot be read from its recorded path and was not given; an assign report without its model, or with
            a model the data were not annotated with.

    Examples:
        >>> import scarlet
        >>> scarlet.report("reference.scarlet")

        The same from the command line:

            scarlet report reference.scarlet
    """
    _stop_if_unguarded_import()
    from scarlet.html_report.core import report as _report
    return _report(run, reference=reference, query=query, model=model, query_labels=query_labels, output=output,
                   self_assessment=self_assessment,
                   n_folds=n_folds, seed=seed, umap_cells=umap_cells, umap_all_genes=umap_all_genes, alpha=alpha,
                   ncpu=ncpu, log=_verbosity.say)


COMMANDS = {"learn": learn, "assign": assign, "transfer": transfer, "report": report}
