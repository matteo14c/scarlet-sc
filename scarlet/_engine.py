"""SCARLET internal engine (not a public interface).

The steps behind learn and transfer (scarlet.api): the technology presets (TECHNOLOGIES), the marker search's
argument list (build_marker_finder_argv). Nothing here parses a command line.
"""

from __future__ import annotations

import os
import tempfile

os.environ.setdefault("NUMBA_CACHE_DIR", tempfile.gettempdir())
os.environ.setdefault("MPLCONFIGDIR", tempfile.gettempdir())
os.environ.setdefault("HDF5_USE_FILE_LOCKING", "FALSE")

import argparse
import gc
import sys
from pathlib import Path



def ensure_parent(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)

def default_background_path(output_universe: str) -> Path:
    universe = Path(output_universe)
    stem = universe.with_suffix("")
    return stem.parent / f"{stem.name}_background_universe.tsv"

# The technology presets: marker sizes (50 markers per cell type, ranked on a shortlist of 150) and the marker
# score's background penalty (gamma). Both presets read detection (count > 0); they differ only in the penalty.
TECHNOLOGIES = {
    "sparse": {
        "top_final": 50,
        "background_genes_N": 150,
        "marker_gamma": 6.0,
    },
    "dense": {
        "top_final": 50,
        "background_genes_N": 150,
        "marker_gamma": 0.0,
    },
}
SIZE_KEYS = ("top_final", "background_genes_N")


def apply_technology(args: argparse.Namespace) -> None:
    tech_name = getattr(args, "technology", "sparse")
    args.technology = tech_name

    technology = TECHNOLOGIES[tech_name]
    if getattr(args, "marker_gamma", None) is None:
        args.marker_gamma = technology["marker_gamma"]


def apply_preset_sizes(args: argparse.Namespace) -> None:
    """The marker sizes not given explicitly are the technology preset's; `sizes_from` records where they come
    from."""
    for key in SIZE_KEYS:
        if getattr(args, key, None) is None:
            setattr(args, key, TECHNOLOGIES[args.technology][key])
    args.sizes_from = f"the technology preset ({args.technology})"


def configure_marker_defaults(args: argparse.Namespace) -> None:
    apply_technology(args)
    apply_preset_sizes(args)

# The method's fixed settings (scarlet.api._FIXED_SETTINGS): the marker finder's parser holds them as defaults, not
# options, so they are not forwarded; a namespace that sets another value is refused.
_FIXED_MARKER_SETTINGS = dict(bootstrap_order="keyed", min_genes_expressed=0, normalize_active_genes=False,
                              active_genes_mean=1000, active_genes_sd=250, marker_alpha=1.5, marker_beta=1.0,
                              marker_eps=0.0001, marker_alpha_rule=None, marker_alpha_rho=0.05, classifier_bg_mode=None,
                              universe_mode="natural", universe_temper=0.5, external_universe=None,
                              exclude_blocklist_categories=None, include_blocklist_categories=None, compute_universe=True)


def build_marker_finder_argv(
    args: argparse.Namespace,
    *,
    output_markers: Path,
    output_universe: Path,
    output_background: Path,
) -> list[str]:
    """The marker-finder argv for a resolved namespace (every option forwarded)."""
    for key, value in _FIXED_MARKER_SETTINGS.items():
        if getattr(args, key, value) != value:
            raise ValueError(f"{key} is a fixed setting of the method ({value!r}); got {getattr(args, key)!r}")
    cmd_args = [
        "--in", args.input,
        "--group_key", args.group_key,
        "--gene_col", args.gene_col,
        "--output_markers", str(output_markers),
        "--output_universe", str(output_universe),
        "--output_background_universe", str(output_background),
        "--n_iterations", str(args.n_iterations),
        "--ncpu", str(args.ncpu),
        "--n_cells", str(args.n_cells),
        "--background_genes_N", str(args.background_genes_N),
        "--top_final", str(args.top_final),
        "--active_genes_seed", str(args.active_genes_seed),
        "--gamma", str(args.marker_gamma),
    ]
    if args.blocklist:
        cmd_args.extend(["--blocklist", args.blocklist])
    # forwarded only when not the default, so default command lines are unchanged
    if getattr(args, "imbalanced", "auto") not in (None, "auto"):
        cmd_args.extend(["--imbalanced", args.imbalanced])
    if getattr(args, "blocklist_modules", "auto") not in (None, "auto"):
        cmd_args.extend(["--blocklist-modules", args.blocklist_modules])

    cmd_args.extend(args.marker_extra)
    return cmd_args

def run_marker_finder(
    args: argparse.Namespace,
    *,
    output_markers: Path,
    output_universe: Path,
    output_background: Path,
    dry_run: bool,
    adata=None,
) -> None:
    """adata: the input in memory (args.input then only names it in the records)."""
    cmd_args = build_marker_finder_argv(
        args,
        output_markers=output_markers,
        output_universe=output_universe,
        output_background=output_background,
    )

    print("+ scarlet.marker_finder: " + " ".join(cmd_args), flush=True)
    if dry_run:
        return

    from scarlet.marker_finder import build_parser, validate_args, run_marker_finder_core
    parser = build_parser()
    if adata is None:
        parsed_args = validate_args(parser, parser.parse_args(cmd_args))
    else:   # nothing to read: the file checks do not apply
        from scarlet.marker_finder import validate_marker_params
        parsed_args = validate_marker_params(parser, parser.parse_args(cmd_args))
    context = {
        "entry": getattr(args, "_entry", "scarlet"),
        "subcommand": getattr(args, "marker_command", None),
        "argv": sys.argv[1:],
        "marker_finder_argv": cmd_args,
        "technology": getattr(args, "technology", None),
        "sizes_from": getattr(args, "sizes_from", None),
    }
    # the full in-class table (the model's in_class member): a private setting, so the recorded parameters and
    # command line stay as they were
    parsed_args._output_in_class = getattr(args, "_output_in_class", None)
    run_marker_finder_core(parsed_args, context=context, adata=adata)


def run_find(args: argparse.Namespace, adata=None) -> None:
    configure_marker_defaults(args)
    output_markers = Path(args.output_markers)
    output_universe = Path(args.output_universe)
    output_background = Path(args.output_background_universe or default_background_path(args.output_universe))
    for path in (output_markers, output_universe, output_background):
        ensure_parent(path)
    run_marker_finder(
        args,
        output_markers=output_markers,
        output_universe=output_universe,
        output_background=output_background,
        dry_run=args.dry_run,
        adata=adata,
    )
    # memory: the marker search's AnnData objects (the blocklist's column view and its parent) form reference
    # cycles, which Python frees only when its cyclic collector runs, and large arrays do not trigger it; a transfer
    # runs about twenty marker searches, so without this they pile up (mouse brain, blocklist auto: 7.2 GB instead
    # of 3.1 GB with one CPU). Collecting once per search frees them; the results do not change.
    gc.collect()
