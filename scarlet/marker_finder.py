import os
import argparse
import functools
import multiprocessing as mp
from pathlib import Path
import numpy as np
import pandas as pd
import anndata as ad
from tqdm import tqdm

from scarlet import _verbosity
import sys
from scipy.stats import mannwhitneyu
from scipy import sparse
import networkx as nx
import h5py

from scarlet.input_checks import NonCountInputError, check_non_negative

DEFAULT_BLOCKLIST_PATH = Path(__file__).resolve().with_name("scarlet_blocklist_universal.tsv")

class DefaultsAndRawTextFormatter(
    argparse.ArgumentDefaultsHelpFormatter,
    argparse.RawDescriptionHelpFormatter,
):
    pass

def has_hdf5_signature(path):
    try:
        with open(path, "rb") as fh:
            return fh.read(8) == b"\x89HDF\r\n\x1a\n"
    except OSError:
        return False

def get_gene_names(adata, gene_col="index"):
    """Safely extract gene names from the index or from an adata.var column."""
    if gene_col == "index":
        return adata.var_names.tolist()

    if gene_col in adata.var.columns:
        return adata.var[gene_col].tolist()
    else:
        raise KeyError(f"Column '{gene_col}' does not exist in adata.var.")

def normalize_gene_symbol(gene):
    """Normalize gene symbols for robust cross-dataset matching."""
    return str(gene).strip().upper()

def calculate_empirical_universe(adata):
    """Compute P_universe as positive-expression frequency across all cells."""
    X_positive = positive_pattern(adata.X)
    return np.asarray(matrix_mean_vector(X_positive), dtype=float)

def load_external_universe(path):
    """Load an external universe with gene_symbol and expression_probability columns."""
    if path.lower().endswith(".parquet"):
        universe_df = pd.read_parquet(path)
    else:
        universe_df = pd.read_csv(path, sep=None, engine="python")

    required_cols = {"gene_symbol", "expression_probability"}
    missing_cols = required_cols.difference(universe_df.columns)
    if missing_cols:
        raise ValueError(
            f"--external-universe must contain columns {sorted(required_cols)}. "
            f"Missing columns: {sorted(missing_cols)}"
        )

    universe_df = universe_df[["gene_symbol", "expression_probability"]].copy()
    universe_df["gene_symbol_norm"] = universe_df["gene_symbol"].map(normalize_gene_symbol)
    universe_df["expression_probability"] = pd.to_numeric(
        universe_df["expression_probability"],
        errors="coerce"
    )
    universe_df = universe_df.dropna(subset=["gene_symbol_norm", "expression_probability"])
    universe_df["expression_probability"] = universe_df["expression_probability"].clip(0.0, 1.0)

    if universe_df.empty:
        raise ValueError(f"--external-universe file '{path}' does not contain valid probabilities.")

    if universe_df["gene_symbol_norm"].duplicated().any():
        print("[!] WARNING: duplicate genes in the external universe after normalization; using mean probabilities.")

    return universe_df.groupby("gene_symbol_norm")["expression_probability"].mean()

def resolve_universe_probabilities(adata, gene_names, args, empirical=None):
    """
    Return the P_universe vector aligned to the current genes.
    Priority: --compute-universe, then --external-universe, then empirical fallback with warning.
    """
    if empirical is None:
        empirical = calculate_empirical_universe(adata)

    if args.compute_universe:
        if args.external_universe:
            print("[!] WARNING: both --compute-universe and --external-universe are set; using the empirical universe.")
        print("[SCARLET] P_universe: computing empirical probabilities from the current dataset.")
        args._universe_source = "empirical"
        return empirical

    if args.external_universe:
        print(f"[SCARLET] P_universe: loading external universe from {args.external_universe}")
        external_universe = load_external_universe(args.external_universe)
        normalized_genes = [normalize_gene_symbol(gene) for gene in gene_names]
        external_values = pd.Series(normalized_genes).map(external_universe).to_numpy(dtype=float)
        missing_mask = np.isnan(external_values)

        if missing_mask.any():
            print(
                f"[!] WARNING: {missing_mask.sum()} dataset genes are missing from "
                "the external universe; using local empirical frequencies for these genes."
            )
            external_values[missing_mask] = empirical[missing_mask]

        args._universe_source = "external"
        return external_values

    print("[!] WARNING: neither --compute-universe nor --external-universe was provided; using empirical P_universe from the current dataset.")
    args._universe_source = "empirical_fallback"
    return empirical

# The matching rules, recorded in every params.json. The Symbol match is case-insensitive, so a
# human list also removes the mouse genes whose symbol differs only by case (Fos = FOS, Rpl10 = RPL10).
# "symbol_ci_canonical" (the default) also undoes upstream decorations of the symbol the dataset
# uses: a '__chrN'-style suffix (Muraro 'ATF3__chr1') and R's '-' -> '.' mangling ('MT.CO1', 'Rpl10.ps1').
# "symbol_ci" is the earlier matching (kept for reproducing earlier runs).
BLOCKLIST_SYMBOL_MATCH_MODES = ("ci_canonical", "ci")
BLOCKLIST_MATCH_RULES = {
    "ci_canonical": {
        "id": "ensg_exact+symbol_ci_canonical",
        "ENSG_ID": "exact, case-sensitive, whitespace-stripped dataset gene name == ENSG_ID entry",
        "Symbol": ("canonical key on both sides: whitespace-stripped, a trailing '__...' decoration removed, "
                   "'.' treated as '-', upper-cased"),
        "combine": "a gene is removed if either column matches",
    },
    "ci": {
        "id": "ensg_exact+symbol_ci",
        "ENSG_ID": "exact, case-sensitive, whitespace-stripped dataset gene name == ENSG_ID entry",
        "Symbol": "case-insensitive: upper-cased, whitespace-stripped dataset gene name == upper-cased Symbol entry",
        "combine": "a gene is removed if either column matches",
    },
}


def canonical_symbol_key(gene):
    """Canonical Symbol key: 'ATF3__chr1' -> 'ATF3', 'MT.CO1' -> 'MT-CO1', 'Rpl10.ps1' -> 'RPL10-PS1'."""
    g = str(gene).strip()
    if "__" in g and not g.startswith("__"):
        g = g.split("__", 1)[0]
    return g.replace(".", "-").upper()


def _name_keys(gene_names):
    """(stripped names, their upper-case symbols, their canonical keys) as pd.Index, computed once per process for a
    list of names (a transfer learns many times on the same genes)."""
    from scarlet.blocklists import _memo
    names = tuple(str(g) for g in gene_names)

    def compute():
        stripped = pd.Index([str(gene).strip() for gene in names])
        return (stripped, pd.Index([normalize_gene_symbol(gene) for gene in stripped]),
                pd.Index([canonical_symbol_key(g) for g in stripped]))
    return _memo(("name_keys", names), compute)


def blocklist_gene_mask(gene_names, blocklist_df, symbol_match="ci_canonical"):
    """Boolean mask of the blocklisted genes among ``gene_names``, and match counts.

    A gene is blocklisted when its stripped name equals an ``ENSG_ID`` entry (exact, case-sensitive), or
    its Symbol key equals the key of a ``Symbol`` entry. ``symbol_match``: "ci" = upper-cased names (the
    earlier rule); "ci_canonical" = ``canonical_symbol_key`` on both sides (default). The dataset's
    own identifier is therefore matched whichever kind it is (Ensembl id or symbol).
    """
    if symbol_match not in BLOCKLIST_SYMBOL_MATCH_MODES:
        raise ValueError(f"symbol_match must be one of {BLOCKLIST_SYMBOL_MATCH_MODES}, not {symbol_match!r}")
    blocklist_columns = set(blocklist_df.columns)
    block_ids = set()
    if "ENSG_ID" in blocklist_columns:
        block_ids = set(blocklist_df["ENSG_ID"].dropna().astype(str).str.strip())

    symbols = blocklist_df["Symbol"].dropna() if "Symbol" in blocklist_columns else pd.Series([], dtype=object)
    block_symbols = set(symbols.map(normalize_gene_symbol))
    block_keys = set(symbols.map(canonical_symbol_key))

    gene_names, gene_symbols, canonical_keys = _name_keys(gene_names)
    by_id = np.asarray(gene_names.isin(block_ids))
    by_symbol_plain = np.asarray(gene_symbols.isin(block_symbols))
    if symbol_match == "ci_canonical":
        by_symbol = by_symbol_plain | np.asarray(canonical_keys.isin(block_keys))
    else:
        by_symbol = by_symbol_plain
    blocked = by_id | by_symbol
    plain_total = by_id | by_symbol_plain
    counts = {
        "n_genes_in": int(len(gene_names)),
        "n_removed": int(blocked.sum()),
        "n_matched_by_id": int(by_id.sum()),
        "n_matched_by_symbol": int(by_symbol.sum()),
        "n_removed_plain_matching": int(plain_total.sum()),
        "n_removed_canonical_only": int((blocked & ~plain_total).sum()),
        "canonical_only_genes": [str(g) for g in gene_names[blocked & ~plain_total]],
        "matched_columns": [c for c, n in (("ENSG_ID", by_id.sum()), ("Symbol", by_symbol.sum())) if n],
        "symbol_match": symbol_match,
        "match_rule": BLOCKLIST_MATCH_RULES[symbol_match]["id"],
        "match_rule_detail": BLOCKLIST_MATCH_RULES[symbol_match],
    }
    return blocked, counts


def _blocklist_record(args, blocklist_path, gene_names=None, counts=None, warning=None):
    """Provenance of the blocklist actually applied; stored on args for params.json."""
    from scarlet.blocklists import blocklist_identity, gene_id_type
    record = {"identity": blocklist_identity(blocklist_path), "gene_col": getattr(args, "gene_col", "index"),
              "exclude_categories": getattr(args, "exclude_blocklist_categories", None),
              "include_categories": getattr(args, "include_blocklist_categories", None),
              # the gene groups applied (scarlet.blocklist_modules)
              "modules": getattr(args, "_blocklist_modules_record", None)}
    if gene_names is not None:
        record["gene_id_type"] = gene_id_type(gene_names)
    record.update(counts or {"n_removed": 0})
    if warning:
        record["warning"] = warning
    args._blocklist_record = record
    return record


def apply_universal_blocklist(adata, args):
    """
    Remove technical and biological noise genes from the AnnData object.
    """
    if not getattr(args, "blocklist", None):
        args._blocklist_record = None
        return adata
    from scarlet.blocklists import choose_blocklist
    # auto / human / mouse / none or an existing file (a missing file is an error)
    names = get_gene_names(adata, getattr(args, "gene_col", "index")) if str(args.blocklist) == "auto" else None
    blocklist_path = Path(choose_blocklist(args.blocklist, names))

    print(f"\n[SCARLET] Loading the blocklist: {blocklist_path}")
    
    from scarlet.blocklists import read_table
    blocklist_df = read_table(blocklist_path)

    exclude_categories = getattr(args, "exclude_blocklist_categories", None)
    include_categories = getattr(args, "include_blocklist_categories", None)
    if "Category" in blocklist_df.columns:
        # on a list with gene groups (the bundled lists) the rule decides which groups apply, from the labels
        # being learned (decided once per marker-finder run, see _find_markers_from_args); explicit categories override
        from scarlet.blocklist_modules import resolve_modules
        labels = getattr(args, "_blocklist_module_labels", None)
        if labels is None:
            gk = getattr(args, "group_key", None)
            labels = adata.obs[gk].dropna().astype(str).unique().tolist() if gk and gk in adata.obs.columns else []
        res = resolve_modules(blocklist_df["Category"].dropna().astype(str).unique(), labels,
                              getattr(args, "blocklist_modules", "auto"), exclude_categories, include_categories)
        args._blocklist_modules_record = res["record"]
        exclude_categories, include_categories = res["exclude"], res["include"]
        if res["record"]["modular_list"] and not res["record"]["override"]:
            print(f"  -> Blocklist modules ({res['record']['mode']}): ON {res['record']['on']}, OFF {res['record']['off']}"
                  + (f" (cell_cycle OFF by: {res['record']['cell_cycle_off_by']})" if res['record']['cell_cycle_off_by'] else ""))
        if exclude_categories:
            if isinstance(exclude_categories, str):
                exclude_categories = [c.strip() for c in exclude_categories.split(",") if c.strip()]
            blocklist_df = blocklist_df[~blocklist_df["Category"].isin(exclude_categories)]
            print(f"  -> Excluded blocklist categories: {exclude_categories}")
        if include_categories:
            if isinstance(include_categories, str):
                include_categories = [c.strip() for c in include_categories.split(",") if c.strip()]
            blocklist_df = blocklist_df[blocklist_df["Category"].isin(include_categories)]
            print(f"  -> Kept only blocklist categories: {include_categories}")

    blocklist_columns = set(blocklist_df.columns)
    supported_columns = [col for col in ("ENSG_ID", "Symbol") if col in blocklist_columns]
    if not supported_columns:
        print(
            f"[!] WARNING: blocklist '{blocklist_path}' has neither ENSG_ID nor Symbol columns. "
            "No filtering applied."
        )
        _blocklist_record(args, blocklist_path,
                          warning="blocklist has neither ENSG_ID nor Symbol columns; no filtering applied")
        return adata

    raw_names = get_gene_names(adata, args.gene_col)
    blocked, counts = blocklist_gene_mask(raw_names, blocklist_df,
                                          symbol_match=getattr(args, "blocklist_symbol_match", None) or "ci_canonical")
    valid_genes_mask = ~blocked
    n_blocked = counts["n_removed"]
    counts["n_list_rows_used"] = int(len(blocklist_df))
    warning = None
    from scarlet.blocklists import gene_id_type, list_id_counts
    id_type = gene_id_type(raw_names)
    id_counts = list_id_counts(blocklist_df)
    counts["list_id_counts"] = id_counts
    if id_type in ("ensembl_human", "ensembl_mouse"):
        # on Ensembl-indexed data only the ENSG_ID column can match. A list (almost) without ids of the
        # data's species would silently remove (almost) nothing, so a list with ids of the data's species in fewer
        # than 10 % of its rows is refused.
        n_rows = int(len(blocklist_df))
        if n_rows and id_counts[id_type] < 0.1 * n_rows:
            kind = "mouse (ENSMUSG)" if id_type == "ensembl_mouse" else "human (ENSG)"
            msg = (f"blocklist '{blocklist_path}' has a {kind} id in its ENSG_ID column for only "
                   f"{id_counts[id_type]} of {n_rows} rows, but the dataset's gene names are {kind} Ensembl ids, so "
                   "(almost) no gene could be blocklisted. Use a list with ids of the data's species "
                   "(blocklist='mouse' or 'human'), or give gene_col, a column of the genes' symbols in var.")
            print(f"[SCARLET] ERROR: {msg}")
            raise SystemExit(f"ERROR: {msg}")
        # rows without an id of the data's species cannot match
        counts["n_list_rows_without_data_id"] = int(len(blocklist_df)) - id_counts[id_type]
        print(f"  -> {counts['n_list_rows_without_data_id']} of {len(blocklist_df)} list rows carry no "
              f"{id_type} id and cannot match these data.")
    _blocklist_record(args, blocklist_path, raw_names, counts, warning)
    
    if n_blocked:
        print(f"  -> Found {n_blocked} blocklisted genes in the dataset. Removing them...")
        # Filter AnnData using a lightweight view to avoid duplicating memory
        args._blocklist_gene_mask = np.asarray(valid_genes_mask, dtype=bool)  # for run_discovery's pattern
        adata = adata[:, valid_genes_mask]
    else:
        print("  -> No blocklisted genes found in the dataset.")
        
    return adata

def filter_cells_by_detected_genes(adata, min_genes_expressed):
    """
    Keep only cells with at least min_genes_expressed genes detected above zero.
    """
    if min_genes_expressed is None or min_genes_expressed <= 0:
        return adata

    detected = np.asarray((adata.X > 0).sum(axis=1)).ravel()
    keep_mask = detected >= min_genes_expressed
    kept = int(keep_mask.sum())
    total = int(len(keep_mask))
    print(
        f"[SCARLET] Keeping {kept}/{total} cells with at least "
        f"{min_genes_expressed} detected genes."
    )
    return adata[keep_mask].copy()

def collapse_duplicate_gene_names(adata, args):
    """
    Collapse duplicate selected gene names into one binary presence column.

    SCARLET scores detected/not-detected expression. If the selected gene
    column contains repeated symbols, the merged gene is considered expressed
    when any duplicate feature is expressed.
    """
    gene_names = pd.Index([str(gene).strip() for gene in get_gene_names(adata, args.gene_col)])
    duplicate_count = int(gene_names.duplicated().sum())
    if duplicate_count == 0:
        return adata

    print(f"[SCARLET] Collapsing {duplicate_count} duplicate gene-name columns with binary OR.")
    codes, unique_genes = pd.factorize(gene_names, sort=False)
    first_indices = np.full(len(unique_genes), -1, dtype=int)
    for position, code in enumerate(codes):
        if first_indices[code] == -1:
            first_indices[code] = position

    mapper = sparse.csr_matrix(
        (
            np.ones(adata.n_vars, dtype=np.uint8),
            (np.arange(adata.n_vars), codes),
        ),
        shape=(adata.n_vars, len(unique_genes)),
    )
    collapsed_X = ((adata.X > 0) @ mapper) > 0

    collapsed = adata[:, first_indices].copy()
    collapsed.X = collapsed_X.astype(np.uint8)
    collapsed.var_names = pd.Index(unique_genes.astype(str))
    if args.gene_col != "index":
        collapsed.var[args.gene_col] = collapsed.var_names
    return collapsed

def class_balanced_universe(X_bin, group_indices, cell_types):
    """Strategy C: per gene, the mean over classes of the within-class detection rate (each class weighted equally)."""
    rates = [np.asarray(matrix_mean_vector(X_bin[group_indices[ct], :]), dtype=float)
             for ct in cell_types if len(group_indices[ct])]
    return np.mean(np.vstack(rates), axis=0)


def tempered_class_weights(sizes, t):
    """Strategy E: class weights proportional to share**t, renormalised to sum to 1
    (t = 1: the natural proportions; t = 0: equal weights, as class_balanced)."""
    sizes = np.asarray(sizes, dtype=float)
    share = sizes / sizes.sum()
    w = np.power(share, float(t))
    return w / w.sum()


def tempered_universe(X_bin, group_indices, cell_types, t):
    """Strategy E: per gene, sum over classes of w_c * (within-class detection rate), w = tempered_class_weights."""
    cts = [ct for ct in cell_types if len(group_indices[ct])]
    rates = np.vstack([np.asarray(matrix_mean_vector(X_bin[group_indices[ct], :]), dtype=float) for ct in cts])
    w = tempered_class_weights([len(group_indices[ct]) for ct in cts], t)
    return w @ rates


def alpha_fmax(f_max, alpha_min=1.5, rho=0.05):
    """Strategy B: the marker-filter alpha that still lets a gene expressed only in the largest class (share f_max)
    pass p_in - p_bg >= p_in / alpha with margin rho: max(alpha_min, 1 / (1 - f_max - (1 - f_max) * rho)).
    With the natural universe such a gene has p_bg >= f_max * p_in, so under a fixed alpha it is zeroed whenever
    f_max > 1 - 1/alpha (1/3 at alpha 1.5); a larger alpha is more permissive."""
    denom = 1.0 - f_max - (1.0 - f_max) * rho
    if denom <= 0:
        return float(alpha_min)
    return float(max(alpha_min, 1.0 / denom))


IMBALANCED_THRESHOLD = 1 / 3


STANDARD, BALANCED = "standard", "balanced for unequal class sizes"


def resolve_imbalanced_profile(args, labels):
    """imbalanced auto | yes | no: how the marker search treats unequal class sizes, decided once per run.

    f_max is the largest class's share of the labels learned from. auto: balanced when f_max > 1/3, else standard;
    yes: balanced; no: standard.
      standard: the marker filter's alpha is fixed; the classifier's background is the gene universe;
      balanced: alpha = max(alpha, 1/(1 - f_max - (1 - f_max)*rho)), and the classifier's background uses the
         class-balanced rates (markers and the background file keep the natural universe).
    1/3 is where a fixed filter fails: at alpha 1.5 a gene expressed only in a class with more than a third of the cells
    has p_bg > p_in / 3 and is always removed. Sets args._imbalanced_record, args._eff_alpha_rule and
    args._eff_classifier_bg_mode."""
    mode = getattr(args, "imbalanced", "auto") or "auto"
    lab = pd.Series(np.asarray(labels)).dropna()
    counts = lab.value_counts()
    n = int(counts.sum())
    f_max = float(counts.max() / n) if n else None
    note = None
    if mode == "yes":
        profile = BALANCED
    elif mode == "no":
        profile = STANDARD
    else:
        profile = BALANCED if (f_max is not None and f_max > IMBALANCED_THRESHOLD) else STANDARD
        if profile == BALANCED and getattr(args, "external_universe", None):
            profile, note = STANDARD, "external universe given: auto resolves to standard (balanced needs the class rates of the input)"
    eff_alpha = "fmax" if profile == BALANCED else "fixed"
    eff_bg = "class_balanced" if profile == BALANCED else "natural"
    overrides = []
    if getattr(args, "alpha_rule", None) is not None:
        if args.alpha_rule != eff_alpha:
            overrides.append(f"alpha_rule={args.alpha_rule}")
        eff_alpha = args.alpha_rule
    if getattr(args, "classifier_bg_mode", None) is not None:
        if args.classifier_bg_mode != eff_bg:
            overrides.append(f"classifier_bg_mode={args.classifier_bg_mode}")
        eff_bg = args.classifier_bg_mode
    args._eff_alpha_rule, args._eff_classifier_bg_mode = eff_alpha, eff_bg
    args._imbalanced_record = {"mode": mode, "threshold": IMBALANCED_THRESHOLD, "f_max": f_max, "n_cells": n,
                               "n_classes": int(len(counts)), "profile": profile, "alpha_rule": eff_alpha,
                               "classifier_bg_mode": eff_bg, "overrides": overrides, "note": note}
    return profile


def resolve_alpha_rule(args, group_indices):
    """Set args._alpha_used (and the rule's inputs) for this run; alpha_rule 'fixed' keeps --alpha."""
    rule = getattr(args, "_eff_alpha_rule", None) or getattr(args, "alpha_rule", None) or "fixed"
    sizes = np.array([len(v) for v in group_indices.values()], dtype=float)
    f_max = float(sizes.max() / sizes.sum()) if sizes.sum() else None
    if rule == "fmax":
        rho = float(getattr(args, "alpha_rule_rho", 0.05))
        args._alpha_used = alpha_fmax(f_max, args.alpha, rho)
    else:
        args._alpha_used = float(args.alpha)
    args._alpha_rule_record = {"alpha_rule": rule, "alpha_min": float(args.alpha),
                               "rho": float(getattr(args, "alpha_rule_rho", 0.05)) if rule == "fmax" else None,
                               "f_max": f_max, "n_cells": int(sizes.sum()), "n_classes": int(len(sizes)),
                               "alpha_used": args._alpha_used}
    return args._alpha_used


def scoring_f(p_in, p_bg, alpha, beta, gamma, eps):
    log_ratio = np.log2((p_in + eps) / (p_bg + eps))
    bg_penalty = (1 - p_bg)**gamma
    score = bg_penalty * (p_in**beta) * log_ratio
    score[(p_in - p_bg) < (p_in / alpha)] = 0
    return score

_BOOTSTRAP_X_BIN = None
_BOOTSTRAP_GROUP_INDICES = None
_BOOTSTRAP_CELL_TYPES = None
_BOOTSTRAP_P_UNIVERSE = None
_BOOTSTRAP_PARAMS = None
_BOOTSTRAP_KEEP_HISTORY = False

def _bool_csr_column_mean(X):
    """Column means of a boolean CSR matrix, bit-identical to scipy's X.mean(axis=0), without its copies.

    scipy casts the matrix to float64, multiplies every stored value by c = 1/n_rows and sums
    each column by a product with a vector of ones: a column's mean is c added to 0.0 once per True entry, in row
    order. All addends are equal, so the result depends only on the column's count k: the k-th running sum of c.
    Those running sums are computed once (np.cumsum is sequential) and looked up per column; scipy's float64
    copies of the matrix (data and indices, twice) are not made."""
    n_rows, n_cols = X.shape
    c = 1.0 / n_rows
    nz = X.indices if X.data.all() else X.indices[X.data]
    # integer counts per column, in chunks (np.bincount would cast the whole index array to int64 at once); exact
    counts = np.zeros(n_cols, dtype=np.int64)
    for start in range(0, len(nz), 1 << 26):
        counts += np.bincount(nz[start:start + (1 << 26)], minlength=n_cols)
    table = np.zeros(int(counts.max()) + 1 if counts.size else 1, dtype=np.float64)
    if table.size > 1:
        table[1:] = np.cumsum(np.full(table.size - 1, c, dtype=np.float64))
    return table[counts]


def matrix_mean_vector(X):
    """Return per-gene means as a 1D ndarray without densifying sparse matrices."""
    if sparse.issparse(X) and X.ndim == 2 and X.format == "csr" and X.dtype == np.bool_ and X.shape[0] > 0:
        return _bool_csr_column_mean(X)
    return np.asarray(X.mean(axis=0)).ravel()


def binary_core_enabled():
    """The binary core (default on; SCARLET_BINARY_CORE=0 turns it off): once an input is loaded and
    checked, the steps that only read detection (X > 0) keep its detection pattern instead of the value matrix."""
    return os.environ.get("SCARLET_BINARY_CORE", "1").strip().lower() not in ("0", "false", "no", "off")


def positive_pattern(X):
    """X > 0 as a boolean sparse matrix. Memory: for a CSR matrix with sorted indices whose stored
    values are all positive (counts without stored zeros, or a pattern already), the pattern shares X's index arrays
    instead of copying them, and its values are a slice of the shared read-only all-True buffer (scarlet._patterns)
    instead of a new array; the matrix is the same as X > 0. Anything else: X > 0."""
    if sparse.issparse(X) and X.format == "csr" and X.has_sorted_indices and X.nnz and \
            bool(X.data.all() if X.dtype == np.bool_ else (X.data > 0).all()):
        from scarlet._patterns import all_true
        return sparse.csr_matrix((all_true(X.nnz), X.indices, X.indptr), shape=X.shape, copy=False)
    return X > 0

def normalize_active_genes_by_sampling(X, target_mean=1000, target_sd=250, seed=42):
    """
    Convert expression values to a binary matrix while sampling each cell toward
    a comparable number of active genes.

    Cells with an observed positive-gene count at or below target_mean are left
    untouched. Only cells above target_mean are sampled down toward a
    comparable active-gene count.

    Memory: same generator, same draws in the same order and the same result as before; the selected
    positions go into preallocated arrays instead of a Python list of Python ints.
    """
    rng = np.random.default_rng(seed)
    X_csr = X.tocsr() if sparse.issparse(X) else sparse.csr_matrix(X)
    n_cells, n_genes = X_csr.shape
    target_counts = rng.normal(target_mean, target_sd, n_cells)
    target_counts = np.clip(target_counts, 1, n_genes).astype(int)

    n_stored = np.diff(X_csr.indptr)
    n_pos = np.zeros(n_cells, dtype=np.int64)  # positive entries per row (the loop's n_positive)
    nonempty = n_stored > 0
    if nonempty.any():
        n_pos[nonempty] = np.add.reduceat(X_csr.data > 0, X_csr.indptr[:-1][nonempty], dtype=np.int64)
    bound = np.where(n_pos <= target_mean, n_pos, np.minimum(target_counts, n_pos))
    out_indices = np.empty(int(bound.sum()), dtype=np.int32)
    out_indptr = np.zeros(n_cells + 1, dtype=np.int64)
    pos = 0
    for i in tqdm(range(n_cells), desc="Active-gene normalization", disable=not _verbosity.step_bars()):
        start = X_csr.indptr[i]
        end = X_csr.indptr[i + 1]
        row_indices = X_csr.indices[start:end]
        row_values = X_csr.data[start:end]

        positive_mask = row_values > 0
        row_indices = row_indices[positive_mask]
        row_values = row_values[positive_mask]
        n_positive = len(row_indices)
        if n_positive == 0:
            out_indptr[i + 1] = pos
            continue

        if n_positive <= target_mean:
            selected = row_indices
        else:
            target_n = min(target_counts[i], n_positive)
            log_values = np.log1p(row_values)
            log_sum = log_values.sum()
            if log_sum > 0:
                probs = log_values / log_sum
            else:
                probs = None
            selected_positions = rng.choice(
                n_positive,
                size=target_n,
                replace=False,
                p=probs,
            )
            selected = row_indices[selected_positions]

        selected = np.sort(selected)
        out_indices[pos:pos + len(selected)] = selected
        pos += len(selected)
        out_indptr[i + 1] = pos

    out_indices = out_indices[:pos]
    out_data = np.ones(len(out_indices), dtype=bool)
    return sparse.csr_matrix(
        (out_data, out_indices, out_indptr),
        shape=(n_cells, n_genes),
        dtype=bool,
    )

def init_bootstrap_worker(X_bin, group_indices, cell_types, p_universe, params, keep_history):
    global _BOOTSTRAP_X_BIN
    global _BOOTSTRAP_GROUP_INDICES
    global _BOOTSTRAP_CELL_TYPES
    global _BOOTSTRAP_P_UNIVERSE
    global _BOOTSTRAP_PARAMS
    global _BOOTSTRAP_KEEP_HISTORY

    _BOOTSTRAP_X_BIN = X_bin
    _BOOTSTRAP_GROUP_INDICES = group_indices
    _BOOTSTRAP_CELL_TYPES = cell_types
    _BOOTSTRAP_P_UNIVERSE = p_universe
    _BOOTSTRAP_PARAMS = params
    _BOOTSTRAP_KEEP_HISTORY = keep_history

BOOTSTRAP_ORDERS = ("keyed", "legacy")


def class_stream_key(indices):
    """Stable 128-bit key of a class: sha256 of its sorted cell positions (never Python's hash()).

    It depends only on which cells form the class, not on the class label or on its position in the
    sorted label list."""
    import hashlib
    digest = hashlib.sha256(np.sort(np.asarray(indices, dtype="<i8")).tobytes()).digest()
    return [int.from_bytes(digest[i:i + 4], "little") for i in range(0, 16, 4)]


def class_rng(iteration_seed, class_key):
    """The keyed bootstrap stream of one class in one iteration: SeedSequence(iteration seed, class key)."""
    return np.random.Generator(np.random.PCG64(np.random.SeedSequence([int(iteration_seed), *class_key])))


_MEAN_TABLES = {}


def _mean_table(n_rows, max_count):
    """_bool_csr_column_mean's table for a sample of n_rows rows: entry k is the k-th running sum of 1 / n_rows
    (np.cumsum is sequential, so a longer table has the same first entries); kept per process and sample size."""
    t = _MEAN_TABLES.get(n_rows)
    if t is None or len(t) <= max_count:
        size = max(n_rows, max_count) + 1
        t = np.zeros(size, dtype=np.float64)
        t[1:] = np.cumsum(np.full(size - 1, 1.0 / n_rows, dtype=np.float64))
        _MEAN_TABLES[n_rows] = t
    return t


def _marker_score(p_in, p_bg, pbg_eps, penalty, alpha, beta, eps):
    """scoring_f with its parts that depend only on the background computed once per learn: pbg_eps = p_bg + eps,
    penalty = (1 - p_bg)**gamma (the same operations on the same array, so the same values)."""
    log_ratio = np.log2((p_in + eps) / pbg_eps)
    score = penalty * (p_in**beta) * log_ratio
    score[(p_in - p_bg) < (p_in / alpha)] = 0
    return score


def _top_positive(scores, k):
    """The genes bootstrap_iteration counts as stable in a round: those of argpartition(scores, -k)[-k:] with a
    score > 0 (None when no score is > 0). Found from the positive scores: all of them when there are at most k;
    otherwise those above the k-th largest value and those equal to it when they are exactly the ones needed. Only when
    ties straddle that value does the choice depend on argpartition, which is then called as before."""
    if np.isnan(scores).any():
        if not np.any(scores > 0):
            return None
        top_idx = np.argpartition(scores, -k)[-k:]
        return top_idx[scores[top_idx] > 0]
    pos = np.flatnonzero(scores > 0)
    if len(pos) == 0:
        return None
    if len(pos) <= k:
        return pos
    v = scores[pos]
    kth = len(pos) - k
    t = np.partition(v, kth)[kth]
    gt, eq = pos[v > t], pos[v == t]
    if len(gt) + len(eq) == k:
        return np.concatenate([gt, eq])
    top_idx = np.argpartition(scores, -k)[-k:]
    return top_idx[scores[top_idx] > 0]


def _gene_name_rank(names, n_genes):
    """Each gene's rank in the string order of its name as given (var_names); duplicate names keep their column order.
    Without names (internal callers), the column order."""
    if names is None:
        return np.arange(n_genes, dtype=np.int64)
    names = np.array([str(g) for g in names], dtype=str)
    if len(names) != n_genes:
        raise ValueError(f"{len(names)} gene names for {n_genes} genes")
    rank = np.empty(n_genes, dtype=np.int64)
    rank[np.argsort(names, kind="stable")] = np.arange(n_genes)
    return rank


def _top_positive_tied(scores_kg, c, k, name_rank):
    """The genes a bootstrap round counts as stable for class c (rows of scores_kg: the round's scores per class): those
    of _top_positive, except that genes exactly tied at the k-th largest positive score are taken in a fixed order:
    lower specificity first (the gene's maximum score over the other classes in the same round: the more specific gene
    first), then the gene name as a string (name_rank). Not argpartition, whose choice among tied values depends on
    numpy's CPU dispatch."""
    scores = scores_kg[c]
    if np.isnan(scores).any():
        return _top_positive(scores, k)
    pos = np.flatnonzero(scores > 0)
    if len(pos) <= k:
        return pos if len(pos) else None
    v = scores[pos]
    t = np.partition(v, len(pos) - k)[len(pos) - k]
    gt, eq = pos[v > t], pos[v == t]
    if len(gt) + len(eq) == k:
        return np.concatenate([gt, eq])
    others = np.delete(scores_kg[:, eq], c, axis=0)
    spec = others.max(axis=0) if len(others) else np.zeros(len(eq))
    order = np.lexsort((name_rank[eq], spec))
    return np.concatenate([gt, np.sort(eq[order[:k - len(gt)]])])


def bootstrap_iteration(seed):
    """One round of the marker search. The fast path (keyed order, a detection pattern without stored False): each
    class's per-gene counts from its sampled rows' positions, the column means from the cached table, the score's
    background parts once per learn, the stable genes by value; arrays filled per class (contiguous) and returned as
    (genes x classes) views. The same values as _bootstrap_iteration_reference."""
    P = _BOOTSTRAP_PARAMS
    if not P.get("fast"):
        return _bootstrap_iteration_reference(seed)
    X = _BOOTSTRAP_X_BIN
    ind, ip = X.indices, X.indptr
    n_genes, K = X.shape[1], len(_BOOTSTRAP_CELL_TYPES)
    scores_kg = np.zeros((K, n_genes), dtype=float)
    p_in_kg = np.zeros((K, n_genes), dtype=float)
    pres_kg = np.zeros((K, n_genes), dtype=np.uint8)
    alpha, beta, eps, n_cells = P["alpha"], P["beta"], P["eps"], P["n_cells"]
    p_bg, pbg_eps, penalty = _BOOTSTRAP_P_UNIVERSE, P["pbg_eps"], P["penalty"]
    top_k = min(P["top_k"], n_genes)
    keep = _BOOTSTRAP_KEEP_HISTORY or P.get("collect_drawn")
    samples = {} if keep else None
    for ct_idx, ct in enumerate(_BOOTSTRAP_CELL_TYPES):
        indices = _BOOTSTRAP_GROUP_INDICES[ct]
        if len(indices) == 0:
            continue
        size = min(n_cells, len(indices))
        sample_idx = indices[class_rng(seed, P["class_keys"][ct_idx]).integers(0, len(indices), size=size)]
        nz = np.concatenate([ind[ip[r]:ip[r + 1]] for r in sample_idx.tolist()])
        counts = np.bincount(nz, minlength=n_genes)
        p_in = _mean_table(size, int(counts.max()) if counts.size else 0)[counts]
        scores = _marker_score(p_in, p_bg, pbg_eps, penalty, alpha, beta, eps)
        scores_kg[ct_idx] = scores
        p_in_kg[ct_idx] = p_in
        if keep:
            samples[ct] = sample_idx.astype(np.int64, copy=False)
    # the stable genes, once every class's scores of the round are known (ties use specificity across classes)
    for c, ct in enumerate(_BOOTSTRAP_CELL_TYPES):
        if len(_BOOTSTRAP_GROUP_INDICES[ct]) == 0:
            continue
        sel = _top_positive_tied(scores_kg, c, top_k, P["name_rank"])
        if sel is not None:
            pres_kg[c, sel] = 1
    return scores_kg.T, p_in_kg.T, pres_kg.T, None, samples


def _bootstrap_iteration_reference(seed):
    keyed = _BOOTSTRAP_PARAMS.get("bootstrap_order", "keyed") == "keyed"
    class_keys = _BOOTSTRAP_PARAMS.get("class_keys")
    rng = None if keyed else np.random.RandomState(seed)
    n_genes = _BOOTSTRAP_X_BIN.shape[1]
    n_cell_types = len(_BOOTSTRAP_CELL_TYPES)
    sum_scores = np.zeros((n_genes, n_cell_types), dtype=float)
    sum_p_in = np.zeros((n_genes, n_cell_types), dtype=float)
    count_presence = np.zeros((n_genes, n_cell_types), dtype=np.int64)
    iter_p_in = None
    iter_sample_indices = {} if _BOOTSTRAP_KEEP_HISTORY else None

    alpha = _BOOTSTRAP_PARAMS["alpha"]
    beta = _BOOTSTRAP_PARAMS["beta"]
    gamma = _BOOTSTRAP_PARAMS["gamma"]
    eps = _BOOTSTRAP_PARAMS["eps"]
    n_cells = _BOOTSTRAP_PARAMS["n_cells"]
    top_k = min(_BOOTSTRAP_PARAMS["top_k"], n_genes)

    for ct_idx, ct in enumerate(_BOOTSTRAP_CELL_TYPES):
        indices = _BOOTSTRAP_GROUP_INDICES[ct]
        if len(indices) == 0:
            continue

        size = min(n_cells, len(indices))
        if keyed:
            # each class draws from its own stream keyed by (iteration seed, its cells), so the
            # sample does not depend on the label names or on the sorted-label order.
            sample_idx = indices[class_rng(seed, class_keys[ct_idx]).integers(0, len(indices), size=size)]
        else:
            # legacy: one stream per iteration, consumed in sorted-label order
            sample_idx = rng.choice(indices, size=size, replace=True)
        X_sample = _BOOTSTRAP_X_BIN[sample_idx, :]
        p_in = matrix_mean_vector(X_sample)
        scores = scoring_f(p_in, _BOOTSTRAP_P_UNIVERSE, alpha, beta, gamma, eps)

        sum_scores[:, ct_idx] = scores
        sum_p_in[:, ct_idx] = p_in

        if _BOOTSTRAP_KEEP_HISTORY:
            iter_sample_indices[ct] = sample_idx.astype(np.int64, copy=False)

    scores_kg = np.ascontiguousarray(sum_scores.T)
    for c, ct in enumerate(_BOOTSTRAP_CELL_TYPES):
        if len(_BOOTSTRAP_GROUP_INDICES[ct]) == 0:
            continue
        sel = _top_positive_tied(scores_kg, c, top_k, _BOOTSTRAP_PARAMS["name_rank"])
        if sel is not None:
            count_presence[sel, c] = 1
    return sum_scores, sum_p_in, count_presence, iter_p_in, iter_sample_indices

def run_bootstrap_tasks(X_bin, group_indices, cell_types, p_universe, args, keep_history):
    params = {
        "alpha": getattr(args, "_alpha_used", None) or args.alpha,
        "beta": args.beta,
        "gamma": args.gamma,
        "eps": args.eps,
        "n_cells": args.n_cells,
        "top_k": args.top_k,
        "bootstrap_order": getattr(args, "bootstrap_order", "keyed"),
        "class_keys": [class_stream_key(group_indices[ct]) for ct in cell_types],
        "name_rank": _gene_name_rank(getattr(args, "_gene_names", None), X_bin.shape[1]),
    }
    # the fast path of bootstrap_iteration: keyed order on a detection pattern without stored False values; the
    # score's background parts once per learn
    params["fast"] = (params["bootstrap_order"] == "keyed" and sparse.issparse(X_bin) and X_bin.format == "csr"
                      and X_bin.dtype == np.bool_ and bool(X_bin.data.all()))
    if params["fast"]:
        p_u = np.asarray(p_universe)
        params.update(pbg_eps=p_u + params["eps"], penalty=(1 - p_u)**params["gamma"],
                      collect_drawn=bool(getattr(args, "_collect_drawn", False)))
    active_seed = getattr(args, "active_genes_seed", 42)
    if active_seed is not None:
        np.random.seed(active_seed)
    seeds = np.random.randint(0, np.iinfo(np.int32).max, size=args.n_iterations)
    ncpu = max(1, int(args.ncpu))

    init_bootstrap_worker(X_bin, group_indices, cell_types, p_universe, params, keep_history)
    try:
        if ncpu == 1:
            for seed in tqdm(seeds, desc="Bootstrap", disable=not _verbosity.step_bars()):
                yield bootstrap_iteration(int(seed))
            return

        # the workers return each iteration's arrays, in iteration order, and the caller adds them one by one as
        # with one CPU: the sums, and so the model, are the same bits whatever ncpu (summing per batch first was not)
        ctx = mp.get_context("fork") if "fork" in mp.get_all_start_methods() else mp.get_context()
        with ctx.Pool(
            processes=ncpu,
            initializer=init_bootstrap_worker,
            initargs=(X_bin, group_indices, cell_types, p_universe, params, keep_history),
        ) as pool:
            iterator = pool.imap(bootstrap_iteration, [int(seed) for seed in seeds], chunksize=1)
            for result in tqdm(iterator, total=len(seeds), desc=f"Bootstrap ({ncpu} CPU)",
                               disable=not _verbosity.step_bars()):
                yield result
    finally:
        # memory: the module's references to this learn's detection matrix and tables are released once its
        # rounds are done (they were kept until the next learn replaced them); every round sets them first
        init_bootstrap_worker(None, None, None, None, None, False)

def scoring_stats(X_bin, group_indices, cell_types, args, drawn=None):
    """(learn records it, assign's Jeffreys floor reads it): per identity, the bootstrap's distinct drawn
    cells n_eff (the draws rebuilt with the same seeds and streams as run_bootstrap_tasks: keyed order,
    min(n_cells, n_k) per round with replacement), the mean number of detected genes over them (Kbar) and over all n_k
    cells; for the background, the classifier's mode, the reference size N and its mean K (natural: over all N cells;
    class_balanced: the mean over identities of their mean K over all cells). K = genes detected per cell in the
    marker finder's detection matrix. Not computed (None) with the legacy bootstrap order."""
    if getattr(args, "bootstrap_order", "keyed") != "keyed":
        return None
    seeds = np.random.RandomState(getattr(args, "active_genes_seed", 42)).randint(0, np.iinfo(np.int32).max,
                                                                                    size=args.n_iterations)
    K = np.asarray(X_bin.sum(axis=1)).ravel().astype(float)
    ids = {}
    for ct in cell_types:
        idx = np.asarray(group_indices[ct])
        if len(idx) == 0:
            continue
        if drawn is not None:   # the cells the bootstrap drew (the same streams), sorted as np.unique sorts them
            drawn_ct = np.flatnonzero(drawn[ct])
        else:
            key = class_stream_key(idx)
            size = min(int(args.n_cells), len(idx))
            drawn_ct = np.unique(np.concatenate([idx[class_rng(int(sd), key).integers(0, len(idx), size=size)]
                                                 for sd in seeds]))
        ids[str(ct)] = {"n_k": int(len(idx)), "n_eff": int(len(drawn_ct)), "Kbar": float(K[drawn_ct].mean()),
                        "Kbar_all": float(K[idx].mean())}
    mode = getattr(args, "_classifier_bg_mode", None) or "natural"
    kbg = float(K.mean()) if mode == "natural" else float(np.mean([v["Kbar_all"] for v in ids.values()]))
    return {"schema": "scarlet.scoring_stats/1", "identities": ids,
            "background": {"mode": mode, "N": int(X_bin.shape[0]), "Kbar": kbg},
            "n_iterations": int(args.n_iterations), "n_cells": int(args.n_cells)}


def run_discovery(adata, args, keep_history=False,skip_top=0, group_indices=None):
    """Run bootstrap discovery and return metrics and history when requested.

    group_indices (discovery v2, step 4): optional precomputed {class: cell positions in `adata`}; the
    classes may overlap (a cell in several classes is sampled in each), and the background universe is computed on
    the cells of `adata` once (unique cells). Without it the classes are the values of obs[args.group_key]."""
    n_obs_in = adata.n_obs
    adata_in = adata
    args._blocklist_gene_mask = None
    adata = apply_universal_blocklist(adata,args)
    blocklist_view = adata if (adata is not adata_in and args._blocklist_gene_mask is not None) else None
    adata = collapse_duplicate_gene_names(adata, args)
    adata = filter_cells_by_detected_genes(adata, args.min_genes_expressed)
    if group_indices is not None and adata.n_obs != n_obs_in:
        raise ValueError("group_indices are cell positions of the input: min_genes_expressed must not remove cells")
    gene_names = get_gene_names(adata, args.gene_col) 
    if group_indices is None:
        cell_types = sorted(adata.obs[args.group_key].unique())
    else:
        cell_types = sorted(group_indices)
    if args.normalize_active_genes:
        print(
            "[SCARLET] Active-gene normalization enabled "
            f"(mean={args.active_genes_mean}, sd={args.active_genes_sd}, seed={args.active_genes_seed})."
        )
        X_bin = normalize_active_genes_by_sampling(
            adata.X,
            target_mean=args.active_genes_mean,
            target_sd=args.active_genes_sd,
            seed=args.active_genes_seed,
        )
    elif blocklist_view is not None and adata is blocklist_view and not adata_in.is_view:
        # memory: the blocklist view's X would be materialised as a float copy of the kept columns just to
        # take its pattern; the pattern of the input's X, restricted to the same columns, is the same matrix
        # ((X > 0)[:, mask] == X[:, mask] > 0) without that copy.
        X_bin = positive_pattern(adata_in.X)[:, args._blocklist_gene_mask]
    else:
        X_bin = positive_pattern(adata.X)
    args._gene_names = np.asarray(adata.var_names).astype(str)   # the tie-break of the stable genes (_top_positive_tied)
    if group_indices is None:
        group_indices = {
            ct: np.where(adata.obs[args.group_key].to_numpy() == ct)[0]
            for ct in cell_types
        }
    else:
        group_indices = {ct: np.asarray(group_indices[ct], dtype=np.int64) for ct in cell_types}
    # Background strategies (opt-in; the defaults leave the outputs unchanged):
    # --universe_mode class_balanced: the empirical universe is the mean over classes of each class's
    # detection rate (classes weighted equally) instead of the detection rate over all cells. It is the single
    # p_universe of this run, so it reaches scoring_f (log ratio, alpha filter and the gamma penalty), the written
    # universe ('bg' column, the classifier's background model) and the background-universe file alike.
    if getattr(args, "universe_mode", "natural") == "class_balanced":
        empirical_universe = class_balanced_universe(X_bin, group_indices, cell_types)
    elif getattr(args, "universe_mode", "natural") == "tempered":
        # tempered: class weights proportional to share**t (t = --universe_temper), used everywhere the universe acts
        empirical_universe = tempered_universe(X_bin, group_indices, cell_types, args.universe_temper)
    else:
        empirical_universe = np.asarray(matrix_mean_vector(X_bin), dtype=float)
    p_universe = resolve_universe_probabilities(adata, gene_names, args, empirical=empirical_universe)
    args._universe_mode = getattr(args, "universe_mode", "natural")
    # --classifier_bg_mode class_balanced: the class-balanced rates go ONLY into the 'bg' column of the
    # written universe (the classifier's background rejection model); scoring, the gamma penalty and the
    # background-universe file keep p_universe.
    if getattr(args, "_imbalanced_record", None) is None:  # direct calls of run_discovery: decide on these labels
        resolve_imbalanced_profile(args, adata.obs[args.group_key].to_numpy())
    args._classifier_bg_mode = getattr(args, "_eff_classifier_bg_mode", None) or "natural"
    p_bg_classifier = None
    if args._classifier_bg_mode == "class_balanced":
        p_bg_classifier = (empirical_universe if args._universe_mode == "class_balanced"
                           else class_balanced_universe(X_bin, group_indices, cell_types))
    # --alpha_rule fmax: alpha = max(--alpha, 1 / (1 - f_max - (1 - f_max) * rho)), f_max = the largest
    # class's share of the cells this run learns from (the training fold in CV). A fixed rule, never tuned.
    resolve_alpha_rule(args, group_indices)
    n_genes = len(gene_names)
    n_cell_types = len(cell_types)
    
    # the sums per class (contiguous rows; transposed once at the end: the same additions in the same order)
    sum_scores_kg = np.zeros((n_cell_types, n_genes), dtype=float)
    sum_p_in_kg = np.zeros((n_cell_types, n_genes), dtype=float)
    count_presence_kg = np.zeros((n_cell_types, n_genes), dtype=np.int64)
    args._collect_drawn = (bool(getattr(args, "scoring_stats", False))
                           and getattr(args, "bootstrap_order", "keyed") == "keyed")
    drawn = {ct: np.zeros(X_bin.shape[0], dtype=bool) for ct in cell_types} if args._collect_drawn else None
    
    storage_sample_indices = [] if keep_history else None

    print(f"  -> Starting bootstrap discovery ({len(cell_types)} clusters, ncpu={args.ncpu})...")
    for scores_arr, p_in_arr, presence_arr, _iter_p_in, iter_sample_indices in run_bootstrap_tasks(
        X_bin,
        group_indices,
        cell_types,
        p_universe,
        args,
        keep_history,
    ):
        sum_scores_kg += scores_arr.T
        sum_p_in_kg += p_in_arr.T
        count_presence_kg += presence_arr.T
        if keep_history:
            storage_sample_indices.append(iter_sample_indices)
        if drawn is not None and iter_sample_indices:
            for ct, smp in iter_sample_indices.items():
                drawn[ct][smp] = True
    sum_scores_arr = np.ascontiguousarray(sum_scores_kg.T)
    sum_p_in_arr = np.ascontiguousarray(sum_p_in_kg.T)
    count_presence_arr = np.ascontiguousarray(count_presence_kg.T)
    del sum_scores_kg, sum_p_in_kg, count_presence_kg

    if getattr(args, "scoring_stats", False):   # what assign's Jeffreys floor (j) needs (learn only)
        args._scoring_stats = scoring_stats(X_bin, group_indices, cell_types, args,
                                            drawn=drawn if drawn is not None and all(m.any() for ct, m in drawn.items()
                                                                                     if len(group_indices[ct])) else None)
    sum_scores = pd.DataFrame(sum_scores_arr, index=gene_names, columns=cell_types)
    count_presence = pd.DataFrame(count_presence_arr, index=gene_names, columns=cell_types)
    avg_p_in = pd.DataFrame(sum_p_in_arr / args.n_iterations, index=gene_names, columns=cell_types)
    avg_p_bg = pd.DataFrame(
        np.tile(p_universe.reshape(-1, 1), (1, n_cell_types)),
        index=gene_names,
        columns=cell_types,
    )
   
    # Initialize output containers before the per-cell-type loop.
    final_rows = []
    final_marker_indices = {}

    for ct in cell_types:
        ct_data = pd.DataFrame({
            'idx': range(len(gene_names)), 
            'gene': gene_names, 
            'stab': count_presence[ct], 
            'score': sum_scores[ct]
        })
        
        sorted_df = ct_data.sort_values(['stab', 'score'], ascending=False)
        
        # the class's markers (top_final genes)
        top_df = sorted_df.iloc[skip_top : skip_top + args.top_final]
        final_rows.append(top_df.assign(cell_type=ct))
        final_marker_indices[ct] = top_df['idx'].values

    # Explicit return variables.
    final_df_concatenated = pd.concat(final_rows)

    return {
        'cell_types': cell_types,
        'final_rows': final_df_concatenated,
        'avg_p_in': avg_p_in,
        'avg_p_bg': avg_p_bg,
        'p_universe': pd.Series(p_universe, index=gene_names, name='expression_probability'),
        'p_bg_classifier': (pd.Series(p_bg_classifier, index=gene_names, name='bg')
                            if p_bg_classifier is not None else None),
        'final_marker_indices': final_marker_indices,
        'iter_p_in_all': None,
        'storage_sample_indices': storage_sample_indices,
        'X_bin': X_bin if keep_history else None
    }
def build_parser():
    """The marker search's own argument list, built by scarlet._engine from learn's and transfer's options (it is not a
    public command). The method's fixed settings are defaults, not options (scarlet.api._FIXED_SETTINGS)."""
    parser = argparse.ArgumentParser(
        description=(
            "SCARLET marker discovery.\n\n"
            "The tool discovers marker genes from an AnnData matrix (detection: count > 0) using bootstrap "
            "sampling and enrichment against the background universe of the input's cells."
        ),
        formatter_class=DefaultsAndRawTextFormatter,
    )

    input_group = parser.add_argument_group("Required inputs")
    input_group.add_argument(
        "--in",
        "--input",
        dest="input",
        required=True,
        help="Input AnnData file (.h5ad).",
    )

    io_group = parser.add_argument_group("General outputs")
    io_group.add_argument("--output_markers", default="top_markers.csv", help="Output marker table.")
    io_group.add_argument(
        "--output_universe",
        default="universe.tsv",
        help="Marker-restricted universe table containing only genes present in the final marker output.",
    )
    io_group.add_argument(
        "--output_background_universe",
        default="background_universe.tsv",
        help="Complete empirical background universe of the input's cells.",
    )

    data_group = parser.add_argument_group("Data columns")
    data_group.add_argument("--group_key", default="cell_type", help="Column in adata.obs containing cluster or cell-type labels.")
    data_group.add_argument(
        "--gene_col",
        type=str,
        default="index",
        help="Column in adata.var containing gene names. Use 'index' to use adata.var_names.",
    )
    data_group.add_argument(
        "--blocklist",
        default=str(DEFAULT_BLOCKLIST_PATH),
        help=(
            "TSV file containing genes to exclude before marker discovery. "
            "The default points to scarlet_blocklist_universal.tsv next to this script."
        ),
    )
    data_group.add_argument(
        "--blocklist-modules", "--blocklist_modules",
        dest="blocklist_modules",
        default="auto",
        help=("Gene groups of the bundled lists to apply: auto (default; the fixed rule on the labels "
              "being learned: core modules ON, stress_extended OFF, cell_cycle OFF when a label names a proliferative "
              "state), all, or a comma-separated list of the modules to keep ON. Ignored for a list without gene "
              "groups."),
    )

    bootstrap_group = parser.add_argument_group("Bootstrap discovery")
    bootstrap_group.add_argument("-n", "--n_iterations", type=int, default=100, help="Number of bootstrap iterations.")
    bootstrap_group.add_argument(
        "--ncpu",
        type=int,
        default=1,
        help="Number of worker processes used to parallelize bootstraps (the result is the same at any ncpu).",
    )
    bootstrap_group.add_argument("-c", "--n_cells", type=int, default=100, help="Cells sampled with replacement per group in each bootstrap.")
    bootstrap_group.add_argument(
        "-N",
        "--background_genes_N",
        dest="top_k",
        metavar="N",
        type=int,
        default=200,
        help=(
            "Number of top candidate genes retained per group as the internal background marker pool. "
            "This larger pool is used for bootstrap stability; it is not the "
            "final number of markers written to --output_markers."
        ),
    )
    bootstrap_group.add_argument(
        "--scoring_stats",
        action="store_true",
        default=argparse.SUPPRESS,
        help="(set by learn and transfer): record per identity the bootstrap's distinct drawn cells (n_eff) and "
             "the mean detected genes, which assign's Jeffreys floor reads (scarlet.scoring).",
    )
    bootstrap_group.add_argument(
        "--top_final",
        type=int,
        default=50,
        help="Final number of marker genes written per group to --output_markers.",
    )
    bootstrap_group.add_argument(
        "--active_genes_seed",
        type=int,
        default=42,
        help="Random seed of the run: it seeds the bootstrap seed stream. It is the only seed the marker finder reads.",
    )

    scoring_group = parser.add_argument_group("Scoring parameters")
    scoring_group.add_argument("--imbalanced", choices=["auto", "yes", "no"], default="auto",
                               help="Unequal class sizes (default auto): auto = balanced when the largest class of the labels "
                                    "learned from holds more than a third of the cells, else standard; yes = always "
                                    "balanced; no = always standard. balanced = the alpha rule fmax plus the "
                                    "class-balanced background of the classifier.")
    scoring_group.add_argument("--gamma", type=float, default=6.0, help="Exponent applied to the background penalty (1 - P_universe).")

    # The method's fixed settings: no public command changes them, so they are not options; they stay defaults so the
    # resolved-parameter record names them as before.
    parser.set_defaults(output_params=None, blocklist_symbol_match="ci_canonical", exclude_blocklist_categories=None,
                        include_blocklist_categories=None, compute_universe=True, external_universe=None,
                        normalize_active_genes=False, active_genes_mean=1000, active_genes_sd=250,
                        bootstrap_order="keyed", min_genes_expressed=0, alpha=1.5, alpha_rule=None, alpha_rule_rho=0.05,
                        classifier_bg_mode=None, universe_mode="natural", universe_temper=0.5, beta=1.0, eps=0.0001)
    return parser

def validate_input_file(parser, args):
    input_path = Path(args.input)
    if not input_path.is_file():
        parser.error(f"input AnnData file not found: {input_path}")
    if not has_hdf5_signature(input_path):
        parser.error(
            f"input file is not a valid .h5ad/HDF5 file: {input_path}. "
            "Use the .h5ad file with --in/--input, or regenerate it if it is truncated/corrupted."
        )


def validate_args(parser, args):
    validate_input_file(parser, args)
    return validate_marker_params(parser, args)


def validate_marker_params(parser, args):
    """Check and resolve every marker-finder parameter except the input file.

    Shared by the CLI (through validate_args) and the in-memory Python API, so both
    resolve defaults with exactly the same rules.
    """
    if args.n_iterations < 1:
        parser.error("--n_iterations must be >= 1.")
    if args.ncpu < 1:
        parser.error("--ncpu must be >= 1.")
    if args.n_cells < 1:
        parser.error("--n_cells must be >= 1.")
    if args.top_k < 1:
        parser.error("--background_genes_N must be >= 1.")
    if args.top_final < 1:
        parser.error("--top_final must be >= 1.")
    return args


# ---------------------------------------------------------------------------
# Shared entry: every SCARLET command that learns markers (scarlet.api learn and transfer, through scarlet._engine)
# ends in _find_markers_from_args with one fully resolved namespace.
# ---------------------------------------------------------------------------

SEED_POLICY = (
    "active_genes_seed is the single marker-finder seed: it seeds the bootstrap seed stream "
    "(np.random.seed before drawing the per-iteration seeds) and, with --normalize_active_genes, "
    "the active-gene sampling. The API/CLI parameter 'seed' (scarlet learn / transfer) sets it. "
    "With bootstrap_order=keyed (default) every class of every iteration samples from its own stream "
    "SeedSequence(iteration seed, sha256 of the class's sorted cell positions); bootstrap_order=legacy "
    "uses one RandomState(iteration seed) per iteration consumed in sorted-label order."
)




def _find_markers_from_args(adata, args):
    """Run marker discovery from a fully resolved namespace.

    This is the single code path behind every SCARLET command that learns markers (scarlet.api learn and transfer,
    through scarlet._engine). It uses `args` as given: nothing is re-parsed and nothing is reset to defaults. Sets
    `args._universe_source`.
    """
    check_non_negative(adata.X, "the marker-finder input")  # non-count input stops SCARLET
    args._universe_source = None
    # the --imbalanced profile is decided once per marker-finder run, on the input labels
    args._imbalanced_record = None
    resolve_imbalanced_profile(args, adata.obs[args.group_key].to_numpy())
    # the blocklist module rule reads the input labels once per run
    args._blocklist_module_labels = adata.obs[args.group_key].dropna().astype(str).unique().tolist()
    print("--- Bootstrap marker discovery ---")
    disc_res = run_discovery(adata, args, keep_history=False)
    final_rows = disc_res['final_rows']
    avg_p_in = disc_res['avg_p_in']
    avg_p_bg = disc_res['avg_p_bg']
    p_universe = disc_res['p_universe']
    p_bg_classifier = disc_res.get('p_bg_classifier')

    bg_col = avg_p_bg.iloc[:, 0] if isinstance(avg_p_bg, pd.DataFrame) else avg_p_bg
    if p_bg_classifier is not None:  # balanced: the classifier's background model only
        bg_col = p_bg_classifier
    universe_df = avg_p_in.loc[final_rows['gene'].unique()].assign(bg=bg_col)
    if getattr(args, "_output_in_class", None):   # the same table for every gene the search scored (the model's in_class)
        args._in_class_df = avg_p_in.assign(bg=bg_col)

    return adata, final_rows, universe_df, p_universe


def _to_jsonable(value):
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, dict):
        return {str(k): _to_jsonable(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [_to_jsonable(v) for v in value]
    if value is None or isinstance(value, (str, int, float, bool)):
        return value
    return str(value)


@functools.lru_cache(maxsize=None)
def _git_commit():
    import subprocess
    try:
        out = subprocess.run(
            ["git", "-C", str(Path(__file__).resolve().parent), "rev-parse", "HEAD"],
            capture_output=True, text=True, timeout=10,
        )
        return out.stdout.strip() or None
    except Exception:
        return None


def marker_params_record(args, context=None, outputs=None):
    """Every marker-finder parameter after default resolution, as actually used."""
    params = {key: _to_jsonable(value) for key, value in sorted(vars(args).items()) if not key.startswith("_")}
    params["background_genes_N"] = params.get("top_k")
    try:
        from scarlet import __version__ as scarlet_version
    except Exception:
        scarlet_version = None
    return {
        "schema": "scarlet.marker_finder.params/1",
        "scarlet_version": scarlet_version,
        "scarlet_file": str(Path(__file__).resolve()),
        "git_commit": _git_commit(),
        "written_at": pd.Timestamp.now().isoformat(timespec="seconds"),
        "cwd": os.getcwd(),
        "seed_policy": SEED_POLICY,
        "universe_source": getattr(args, "_universe_source", None),
        # background strategies: universe mode and the alpha rule as applied (alpha_used, f_max)
        "background_strategy": {"universe_mode": getattr(args, "_universe_mode", getattr(args, "universe_mode", "natural")),
                                "classifier_bg_mode": getattr(args, "_classifier_bg_mode", None) or getattr(args, "classifier_bg_mode", None) or "natural",
                                "imbalanced": getattr(args, "_imbalanced_record", None),
                                "universe_temper": (float(getattr(args, "universe_temper", 0.5))
                                                    if getattr(args, "universe_mode", "natural") == "tempered" else None),
                                **(getattr(args, "_alpha_rule_record", None) or {"alpha_rule": getattr(args, "alpha_rule", "fixed")})},
        # blocklist file (md5, version, species), matched column(s) and genes actually removed
        "blocklist": _to_jsonable(getattr(args, "_blocklist_record", None)),
        "marker_finder": params,
        **({"scoring_stats": _to_jsonable(args._scoring_stats)} if getattr(args, "_scoring_stats", None) else {}),
        "context": _to_jsonable(context or {}),
        "outputs": _to_jsonable(outputs or {}),
    }


def default_params_path(output_markers):
    output_markers = Path(output_markers)
    return output_markers.parent / f"{output_markers.stem}.params.json"


def write_params_json(path, record):
    """Write `record` to `path`; if a JSON object already exists there, extend it."""
    import json
    path = Path(path)
    payload = {}
    if path.exists():
        try:
            with open(path, "r", encoding="utf-8") as handle:
                existing = json.load(handle)
            if isinstance(existing, dict):
                payload = existing
        except (OSError, ValueError):
            payload = {}
    payload.update(record)
    with open(path, "w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2, sort_keys=False)
        handle.write("\n")
    return path


def run_marker_finder_core(args, context=None, adata=None):
    """adata: the input in memory (learn and transfer), read from args.input when None. The caller's object is not
    changed: the marker search works on a new AnnData sharing its obs and var (copied) and its X (or X's pattern)."""
    # binary core: the marker finder reads detection only (X > 0), except the opt-in active-gene
    # normalisation, which samples from the values. The values are checked and released here.
    binary_core = binary_core_enabled() and not args.normalize_active_genes
    if adata is None:
        adata = ad.read_h5ad(args.input)
        if binary_core:
            check_non_negative(adata.X, "the marker-finder input")  # non-count input stops SCARLET
            adata.X = positive_pattern(adata.X)
    else:
        if binary_core:
            check_non_negative(adata.X, "the marker-finder input")
        adata = ad.AnnData(X=positive_pattern(adata.X) if binary_core else adata.X, obs=adata.obs.copy(),
                           var=adata.var.copy())

    adata_res, final_rows, universe_df, p_universe = _find_markers_from_args(adata, args)

    print("\n--- SAVING FINAL OUTPUTS ---")
    final_rows.to_csv(args.output_markers, sep='\t', index=False)
    universe_df.to_csv(args.output_universe, sep='\t')
    if getattr(args, "_output_in_class", None):   # not a recorded output: the model's in_class member (report only)
        args._in_class_df.to_csv(args._output_in_class, sep='\t')
        del args._in_class_df
    outputs = {
        "markers": args.output_markers,
        "universe": args.output_universe,
        "background_universe": None,
    }

    if getattr(args, "_universe_source", None) in {"empirical", "empirical_fallback"}:
        complete_bg = p_universe.rename_axis("gene_symbol").reset_index()
        complete_bg.columns = ["gene_symbol", "expression_probability"]
        complete_bg.to_csv(args.output_background_universe, sep='\t', index=False)
        outputs["background_universe"] = args.output_background_universe
        print(f"[V] Complete empirical universe saved to: {args.output_background_universe}")

    params_path = Path(args.output_params) if getattr(args, "output_params", None) else default_params_path(args.output_markers)
    write_params_json(params_path, marker_params_record(args, context=context, outputs=outputs))
    print(f"[V] Resolved parameters saved to: {params_path}")

    print("\n[V] Process complete. Markers and universe saved.")


