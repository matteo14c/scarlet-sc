import os
import sys
import argparse
import multiprocessing as mp
from pathlib import Path
import anndata as ad
import pandas as pd

_SCORING_X = None
_SCORING_BASE_LL = None
_SCORING_DELTA_LOG_P = None
_SCORING_BASE_LL_BG = None
_SCORING_DELTA_LOG_P_BG = None

class DefaultsAndRawTextFormatter(
    argparse.ArgumentDefaultsHelpFormatter,
    argparse.RawDescriptionHelpFormatter,
):
    pass

def init_scoring_worker(X, base_ll, delta_log_p, base_ll_bg, delta_log_p_bg):
    global _SCORING_X
    global _SCORING_BASE_LL
    global _SCORING_DELTA_LOG_P
    global _SCORING_BASE_LL_BG
    global _SCORING_DELTA_LOG_P_BG

    _SCORING_X = X
    _SCORING_BASE_LL = base_ll
    _SCORING_DELTA_LOG_P = delta_log_p
    _SCORING_BASE_LL_BG = base_ll_bg
    _SCORING_DELTA_LOG_P_BG = delta_log_p_bg

def score_cell_chunk(bounds):
    import numpy as np

    start, end = bounds
    X_chunk = _SCORING_X[start:end, :]
    ll_chunk = np.asarray(_SCORING_BASE_LL + X_chunk @ _SCORING_DELTA_LOG_P)

    if _SCORING_DELTA_LOG_P_BG is None:
        ll_bg_chunk = None
    else:
        ll_bg_chunk = np.asarray(_SCORING_BASE_LL_BG + X_chunk @ _SCORING_DELTA_LOG_P_BG).ravel()

    return start, end, ll_chunk, ll_bg_chunk

def make_chunks(n_items, chunk_size):
    return [(start, min(start + chunk_size, n_items)) for start in range(0, n_items, chunk_size)]

def log_likelihoods_to_probabilities(ll_matrix):
    import numpy as np

    shifted = ll_matrix - np.max(ll_matrix, axis=1, keepdims=True)
    probs = np.exp(shifted)
    probs /= probs.sum(axis=1, keepdims=True)
    return probs.astype(np.float32, copy=False)

def has_hdf5_signature(path):
    try:
        with open(path, "rb") as fh:
            return fh.read(8) == b"\x89HDF\r\n\x1a\n"
    except OSError:
        return False

def validate_input_files(parser, args):
    input_path = Path(args.input)
    universe_path = Path(args.universe)

    if not input_path.is_file():
        parser.error(f"input AnnData file not found: {input_path}")
    if not universe_path.is_file():
        parser.error(f"universe TSV file not found: {universe_path}")

    input_is_hdf5 = has_hdf5_signature(input_path)
    universe_is_hdf5 = has_hdf5_signature(universe_path)
    if not input_is_hdf5 and universe_is_hdf5:
        parser.error(
            "the input file is not an HDF5/.h5ad file, while the universe "
            "file is HDF5. Did you swap input.h5ad and universe.tsv?"
        )
    if not input_is_hdf5:
        parser.error(
            f"input file is not a valid .h5ad/HDF5 file: {input_path}. "
            "Use the .h5ad file with --in/--input, or regenerate it if it "
            "is truncated/corrupted."
        )

    try:
        with open(universe_path, "rb") as fh:
            sample = fh.read(4096)
    except OSError as exc:
        parser.error(f"cannot read universe TSV file {universe_path}: {exc}")
    if b"\t" not in sample:
        parser.error(
            f"universe file does not look like a tab-separated TSV: {universe_path}"
        )

def compute_likelihood_chunks(X, base_ll, delta_log_p, base_ll_bg, delta_log_p_bg, args):
    import numpy as np

    n_cells = X.shape[0]
    chunks = make_chunks(n_cells, args.chunk_size)
    ll_matrix = np.empty((n_cells, delta_log_p.shape[1]), dtype=float)
    ll_bg = np.empty(n_cells, dtype=float) if delta_log_p_bg is not None else None

    init_scoring_worker(X, base_ll, delta_log_p, base_ll_bg, delta_log_p_bg)
    try:
        if args.ncpu == 1:
            for bounds in chunks:
                start, end, ll_chunk, ll_bg_chunk = score_cell_chunk(bounds)
                ll_matrix[start:end, :] = ll_chunk
                if ll_bg is not None:
                    ll_bg[start:end] = ll_bg_chunk
            return ll_matrix, ll_bg

        ctx = mp.get_context("fork") if "fork" in mp.get_all_start_methods() else mp.get_context()
        chunksize = max(1, len(chunks) // (args.ncpu * 4))
        with ctx.Pool(
            processes=args.ncpu,
            initializer=init_scoring_worker,
            initargs=(X, base_ll, delta_log_p, base_ll_bg, delta_log_p_bg),
        ) as pool:
            for start, end, ll_chunk, ll_bg_chunk in pool.imap(score_cell_chunk, chunks, chunksize=chunksize):
                ll_matrix[start:end, :] = ll_chunk
                if ll_bg is not None:
                    ll_bg[start:end] = ll_bg_chunk

        return ll_matrix, ll_bg
    finally:
        # memory: the module's references to the scored cells' detection and the tables are released once the
        # scores are computed (they were kept until the next scoring replaced them)
        init_scoring_worker(None, None, None, None, None)

def build_parser():
    parser = argparse.ArgumentParser(
        description=(
            "SCARLET classifier.\n\n"
            "The tool annotates cells from an AnnData matrix (detection: count > 0) using a marker "
            "universe table produced by SCARLET marker discovery. If the universe contains "
            "a bg column, bg is used as a background rejection model and cells can be "
            "reported as unassigned."
        ),
        formatter_class=DefaultsAndRawTextFormatter,
        epilog=(
            "Examples:\n"
            "  python scarlet_classifier_alpha_gold_gc.py --in input.h5ad --universe universe.tsv\n"
            "  python scarlet_classifier_alpha_gold_gc.py --input input.h5ad --universe universe.tsv "
            "--output annotated.h5ad"
        ),
    )

    input_group = parser.add_argument_group("Required inputs")
    input_group.add_argument(
        "--in",
        "--input",
        dest="input",
        required=True,
        help="Input AnnData file (.h5ad).",
    )
    input_group.add_argument(
        "--universe",
        dest="universe",
        required=True,
        help="Input SCARLET universe scoring table (.tsv).",
    )

    io_group = parser.add_argument_group("General outputs")
    io_group.add_argument("--output", default="scarlet_annotated.h5ad", help="Output annotated AnnData file.")
    io_group.add_argument(
        "--inplace",
        action="store_true",
        help="Write annotations back to the input AnnData file instead of creating a separate output file.",
    )
    io_group.add_argument(
        "--no-save-probs",
        action="store_true",
        help="Do not store the SCARLET identity probability matrix in the output AnnData file.",
    )
    io_group.add_argument(
        "--prob-obsm-key",
        default="scarlet_probs",
        help="adata.obsm key used to store the identity probability matrix.",
    )
    io_group.add_argument(
        "--identity-names-key",
        default="scarlet_identity_names",
        help="adata.uns key used to store the identity names matching the probability matrix columns.",
    )
    io_group.add_argument(
        "--save-probs-to-obs",
        action="store_true",
        help="Also store each identity probability as an adata.obs column using --obs-prob-prefix.",
    )
    io_group.add_argument(
        "--obs-prob-prefix",
        default="scarlet_prob__",
        help="Prefix for optional per-identity probability columns when --save-probs-to-obs is used.",
    )

    data_group = parser.add_argument_group("Data columns")
    data_group.add_argument(
        "--gene_col",
        type=str,
        default="index",
        help="Column in adata.var containing gene names. Use 'index' to use adata.var_names.",
    )

    scoring_group = parser.add_argument_group("Scoring parameters")
    scoring_group.add_argument("--alpha", type=float, default=0.01, help="FDR threshold used for uniqueness and local multi-label inclusion.")
    scoring_group.add_argument("--eps", type=float, default=0.00001, help="Small constant used to stabilize Naive Bayes log probabilities.")

    performance_group = parser.add_argument_group("Scoring performance")
    performance_group.add_argument(
        "--ncpu",
        type=int,
        default=1,
        help=(
            "Number of worker processes used to parallelize likelihood scoring. "
            "Use 1 as the recommended baseline; parallel mode is not yet "
            "memory-optimized for very large .h5ad files."
        ),
    )
    performance_group.add_argument("--chunk_size", type=int, default=10000, help="Number of cells scored per worker task.")

    return parser

def validate_args(parser, args):
    if not 0 < args.alpha < 1:
        parser.error("--alpha must be > 0 and < 1.")
    if args.eps <= 0:
        parser.error("--eps must be > 0.")
    if args.ncpu < 1:
        parser.error("--ncpu must be >= 1.")
    if args.chunk_size < 1:
        parser.error("--chunk_size must be >= 1.")
    return args

UNASSIGNED_LABEL = "unassigned"


def score_cells(adata, universe, eps=0.00001, ncpu=1, chunk_size=10000, gene_col="index", floor=None, _cache=None):
    """The per-cell scoring of `classify`: the Bernoulli naive-Bayes
    log-likelihood of every cell under every identity of `universe` and under its `bg` column. Returns (ll_matrix cells
    x identities, ll_bg or None, identity names, has_bg).
    floor (scarlet.scoring.jeffreys): None (eps: log(p + eps), the scoring without the floor) or (identity -> f_k, f_bg):
    the rates clamped to [f, 1 - f] and scored as log p, log(1 - p) with p clipped to [1e-5, 1 - 1e-5].
    _cache (classify's flip test scores the same cells against the same universe twice): a dict that keeps the
    detection matrix of the scoring genes from the first call for the second."""
    import numpy as np
    import pandas as pd
    universe.index = universe.index.astype(str)

    all_categories = [c for c in universe.columns if c != 'bg']
    has_bg = 'bg' in universe.columns
    if len(all_categories) == 0:
        raise ValueError("Error: universe does not contain biological categories other than bg.")
    if len(all_categories) == 1 and not has_bg:
        raise ValueError("Error: LRT requires at least a second category or a bg column.")
    if universe.index.has_duplicates:
        raise ValueError("Error: universe contains duplicate genes in the index.")
    genes_in_universe = pd.Index(universe.index)

    # Gene synchronization: keep only genes present in the universe and align
    if gene_col == "index":
        adata_genes = pd.Index(adata.var_names)
    else:
        if gene_col not in adata.var.columns:
            raise KeyError(f"Error: gene_col '{gene_col}' was not found in adata.var.")
        adata_genes = pd.Index(adata.var[gene_col].astype(str))

    mask = adata_genes.isin(genes_in_universe)
    if not mask.any():
        raise ValueError("Error: AnnData and universe have no genes in common.")

    found_genes = adata_genes[mask]
    if found_genes.has_duplicates:
        duplicated = found_genes[found_genes.duplicated()].unique()[:5].tolist()
        raise ValueError(f"Error: duplicate genes in the scoring gene set: {duplicated}")

    adata_sub = adata[:, mask]
    universe_sub = universe.loc[found_genes, all_categories]
    bg_sub = universe.loc[found_genes, 'bg'] if has_bg else None
    print(f"Genes used for scoring: {len(found_genes)}/{adata.n_vars} present in universe.tsv")

    # Naive Bayes parameters
    P = universe_sub.values
    if floor is None:
        log_P = np.log(P + eps)
        log_1_minus_P = np.log(1 - P + eps)
    else:   # the Jeffreys floor (scarlet.scoring)
        from scarlet.scoring import clamp
        Pf = clamp(P, np.array([floor[0][k] for k in all_categories]))
        log_P, log_1_minus_P = np.log(Pf), np.log(1 - Pf)
    if has_bg:
        P_bg = bg_sub.values.reshape(-1, 1)
        if floor is None:
            log_P_bg = np.log(P_bg + eps)
            log_1_minus_P_bg = np.log(1 - P_bg + eps)
        else:
            from scarlet.scoring import clamp
            Pb = clamp(P_bg, floor[1])
            log_P_bg, log_1_minus_P_bg = np.log(Pb), np.log(1 - Pb)

    if _cache is not None and "X" in _cache:
        X = _cache["X"]
    else:
        X = adata_sub.X > 0
        if _cache is not None:
            _cache["X"] = X
    # --- 2. LOG-LIKELIHOOD CALCULATION ---
    print(f"Computing likelihoods for {adata.n_obs} cells (ncpu={ncpu}, chunk_size={chunk_size})...")
    base_ll = log_1_minus_P.sum(axis=0)
    delta_log_P = log_P - log_1_minus_P
    if has_bg:
        base_ll_bg = log_1_minus_P_bg.sum(axis=0)
        delta_log_P_bg = log_P_bg - log_1_minus_P_bg
    else:
        base_ll_bg = None
        delta_log_P_bg = None

    class MockArgs:
        def __init__(self, ncpu, chunk_size):
            self.ncpu = ncpu
            self.chunk_size = chunk_size
    mock_args = MockArgs(ncpu, chunk_size)

    ll_matrix, ll_bg = compute_likelihood_chunks(
        X,
        base_ll,
        delta_log_P,
        base_ll_bg,
        delta_log_P_bg,
        mock_args,
    )

    return ll_matrix, ll_bg, all_categories, has_bg


def classify(
    adata: ad.AnnData,
    universe: pd.DataFrame,
    alpha: float = 0.01,
    eps: float = 0.00001,
    ncpu: int = 1,
    chunk_size: int = 10000,
    gene_col: str = "index",
    no_save_probs: bool = False,
    prob_obsm_key: str = "scarlet_probs",
    identity_names_key: str = "scarlet_identity_names",
    save_probs_to_obs: bool = False,
    obs_prob_prefix: str = "scarlet_prob__",
    floor=None,
    flip=None,
    flip_mode: str = "refuse",
) -> ad.AnnData:
    """R0: per cell the best identity (unassigned when the background wins), its unique flag (BH-adjusted LRT of the
    best against the next best or the background, below alpha) and the set of identities within alpha of the best.
    floor: the scoring's Jeffreys floor (score_cells) or None (eps). flip: None, or the Jeffreys
    floor (identity -> f_k, f_bg) of the flip test: the cells are scored under eps and under j, and a cell whose
    highest-likelihood identity differs between the two is recorded (obs scarlet_flip, scarlet_top1_eps,
    scarlet_top1_j, scarlet_reason); flip_mode refuse: such a cell is also refused (best unassigned, not unique),
    flag: its call is kept. obs also holds per cell scarlet_top1 / scarlet_top2 and scarlet_ll_top1 / ll_top2 /
    ll_bg (R0's statistic ll_top1 - max(ll_top2, ll_bg))."""
    import numpy as np
    import pandas as pd
    from scipy.stats import chi2
    from statsmodels.stats.multitest import multipletests
    from scarlet.input_checks import check_non_negative

    check_non_negative(adata.X, "the classifier input")  # non-count input stops SCARLET
    # Memory: a copy of everything but X: the result shares the caller's X (read-only here) instead of
    # copying it; obs / var / uns / obsm etc. are copied as before and the caller's object is left unchanged.
    if adata.is_view or adata.X is None:
        adata = adata.copy()
    else:
        _X = adata.X
        try:
            adata.X = None
            _copy = adata.copy()
        finally:
            adata.X = _X
        _copy.X = _X
        adata = _copy
    det_cache = {}   # the flip test scores the same cells twice: their detection matrix is made once
    ll_matrix, ll_bg, all_categories, has_bg = score_cells(adata, universe, eps=eps, ncpu=ncpu, chunk_size=chunk_size,
                                                           gene_col=gene_col, floor=floor, _cache=det_cache)

    # --- 3. STATISTICAL TEST (LRT) ---
    print("Running LRT anchored to the best label...")
    sorted_indices = np.argsort(ll_matrix, axis=1)

    best_idx = sorted_indices[:, -1]
    best_ll = ll_matrix[np.arange(len(adata)), best_idx]
    if len(all_categories) > 1:
        sec_idx = sorted_indices[:, -2]
        sec_ll = ll_matrix[np.arange(len(adata)), sec_idx]
    else:
        sec_ll = np.full(len(adata), -np.inf)
    if has_bg:
        competitor_ll = np.maximum(sec_ll, ll_bg)
    else:
        competitor_ll = sec_ll
    lrt_1v2 = 2 * (best_ll - competitor_ll)
    lrt_1v2 = np.maximum(lrt_1v2, 0)
    raw_p_1v2 = chi2.sf(lrt_1v2, df=1)

    # Global FDR correction
    _, fdr_p_1v2, _, _ = multipletests(raw_p_1v2, alpha=alpha, method='fdr_bh')

    # --- 4. MULTI-LABEL ASSIGNMENT ---
    # per cell, in descending order of likelihood: the identities whose LRT against the best is not significant, up to
    # the first one that is (chi2.sf on all pairs at once: the same values as one call per pair)
    desc = sorted_indices[:, ::-1]
    S_desc = np.take_along_axis(ll_matrix, desc, axis=1)
    s_top = S_desc[:, 0]
    bg_wins = (ll_bg >= s_top) if has_bg else np.zeros(len(s_top), dtype=bool)
    if S_desc.shape[1] > 1:
        n_inc = np.cumprod(chi2.sf(2 * (s_top[:, None] - S_desc[:, 1:]), df=1) >= alpha, axis=1).sum(axis=1)
    else:
        n_inc = np.zeros(len(s_top), dtype=np.int64)
    cats_ = np.asarray(all_categories, dtype=object)
    top_names = cats_[desc[:, 0]]
    best_labels = np.where(bg_wins, "unassigned", top_names).astype(object).tolist()
    all_labels_list = list(best_labels)
    for i in np.flatnonzero(~bg_wins & (n_inc > 0)):
        all_labels_list[i] = ";".join(cats_[desc[i, :1 + n_inc[i]]].tolist())
    unique_flags = [False if b else u for b, u in zip(bg_wins.tolist(), list(fdr_p_1v2 < alpha))]

    flipped = None
    if flip is not None:   # the call must not depend on how never-seen genes are treated
        from scarlet.scoring import FLIP_REASON
        if flip_mode not in ("refuse", "flag"):
            raise ValueError(f"flip_mode must be refuse or flag (got {flip_mode!r})")
        cats_ = np.asarray(all_categories, dtype=object)
        if floor is None:   # scored under eps: the j scoring is the other one
            ll_e = ll_matrix
            ll_j = score_cells(adata, universe, eps=eps, ncpu=ncpu, chunk_size=chunk_size, gene_col=gene_col,
                               floor=flip, _cache=det_cache)[0]
        else:               # scored under j: the eps scoring is the other one
            ll_e = score_cells(adata, universe, eps=eps, ncpu=ncpu, chunk_size=chunk_size, gene_col=gene_col,
                               _cache=det_cache)[0]
            ll_j = ll_matrix
        top_e = cats_[np.argsort(ll_e, axis=1)[:, -1]].astype(str)
        top_j = cats_[np.argsort(ll_j, axis=1)[:, -1]].astype(str)
        flipped = top_e != top_j
        if flip_mode == "refuse":
            best_labels = [UNASSIGNED_LABEL if f_ else b_ for b_, f_ in zip(best_labels, flipped)]
            all_labels_list = [UNASSIGNED_LABEL if f_ else a_ for a_, f_ in zip(all_labels_list, flipped)]
            unique_flags = [bool(u_) and not f_ for u_, f_ in zip(unique_flags, flipped)]
        print(f"Flip test ({flip_mode}): {int(flipped.sum())} of {len(flipped)} cells flip ({FLIP_REASON})")

    # --- 5. UPDATE AND OUTPUT ---
    adata.obs['scarlet_best'] = best_labels
    adata.obs['scarlet_fdr'] = fdr_p_1v2
    adata.obs['is_scarlet_unique'] = unique_flags
    adata.obs['scarlet_all'] = all_labels_list
    # per cell the two highest-likelihood identities and the log-likelihoods R0 compares (R0 statistic = ll_top1 -
    # max(ll_top2, ll_bg)); top1 is scarlet_best wherever the background does not win (and no flip refuses it)
    cats = np.asarray(all_categories, dtype=object)
    adata.obs['scarlet_top1'] = cats[best_idx].astype(str)
    adata.obs['scarlet_top2'] = cats[sec_idx].astype(str) if len(all_categories) > 1 else ""
    adata.obs['scarlet_ll_top1'] = best_ll
    adata.obs['scarlet_ll_top2'] = sec_ll
    adata.obs['scarlet_ll_bg'] = ll_bg if has_bg else np.nan
    if flipped is not None:
        adata.obs['scarlet_flip'] = flipped
        adata.obs['scarlet_top1_eps'] = top_e
        adata.obs['scarlet_top1_j'] = top_j
        adata.obs['scarlet_reason'] = np.where(flipped, FLIP_REASON, "")
    if not no_save_probs:
        scarlet_probs = log_likelihoods_to_probabilities(ll_matrix)
        adata.obsm[prob_obsm_key] = scarlet_probs
        adata.uns[identity_names_key] = list(all_categories)
        if save_probs_to_obs:
            for idx, identity in enumerate(all_categories):
                adata.obs[f"{obs_prob_prefix}{identity}"] = scarlet_probs[:, idx]

    return adata

def run_classifier_core(args):
    import numpy as np
    import pandas as pd
    import anndata as ad
    from scarlet.input_checks import NonCountInputError

    # --- 1. PREPARATION ---
    try:
        adata = ad.read_h5ad(args.input)
        universe = pd.read_csv(args.universe, sep='\t', index_col=0)
    except Exception as e:
        print(f"Loading error: {e}"); sys.exit(1)

    try:
        adata_annotated = classify(
            adata,
            universe,
            alpha=args.alpha,
            eps=args.eps,
            ncpu=args.ncpu,
            chunk_size=args.chunk_size,
            gene_col=args.gene_col,
            no_save_probs=args.no_save_probs,
            prob_obsm_key=args.prob_obsm_key,
            identity_names_key=args.identity_names_key,
            save_probs_to_obs=args.save_probs_to_obs,
            obs_prob_prefix=args.obs_prob_prefix,
        )
    except NonCountInputError as e:
        print(f"Error: {e}"); sys.exit(1)
    except Exception as e:
        print(f"Scoring error: {e}"); sys.exit(1)

    output_path = args.input if args.inplace else args.output
    if args.inplace:
        print(f"Saving in place to {output_path}...")
    else:
        print(f"Saving to {output_path}...")
    adata_annotated.write_h5ad(output_path)

    unique_flags = adata_annotated.obs['is_scarlet_unique'].values
    n_unique = sum(unique_flags)
    print(f"Done. Unique calls: {n_unique}/{len(adata)} ({(n_unique/len(adata))*100:.1f}%)")

