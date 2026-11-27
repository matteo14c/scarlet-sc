"""The steps of transfer (scarlet.api.transfer).

A transfer runs in draws: each draw repeats the first pass and the refinement with its own seed (depth matching and
the marker search), and each cell's final call is the consensus of the draws.

  - stable_anchors: the cells called uniquely to the same cell type in every draw of the run;
  - refine_plan: per origin, the cell types refined (at least MIN_CELLS stable anchors) and the origins not refined
    (fewer than 2 refined cell types);
  - keep_first: in a draw, the cells that keep that draw's first-pass calls (first-pass best call a cell type that is
    not refined);
  - pooled_background: the refinement's background, the detection rate over all the origin's query cells;
  - consensus: per cell, the most frequent outcome across the draws ("unassigned" is an outcome), kept when it holds
    at least needed(N) of the N draws (8 of 10; the fixed share AGREE_SHARE, rounded up), else unassigned
    (undetermined across draws); with one draw (quick mode) the draw's own call;
  - qc_rows: the transfer check per cell type: the share of the reference markers that recur in markers learned on
    the query with the transferred labels; below QC_FLAG the cell type is flagged 'doubtful transfer'.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

MIN_CELLS = 3  # a cell type is refined with at least 3 anchors (cells the first pass called uniquely)
UNASSIGNED = "unassigned"
QC_FLAG = 0.2
QC_MIN_CELLS = 3
QC_PERFORMANCE = {  # measured when the check was designed (see the documentation's benchmarks); not tuned
    "flag": "share of the reference markers conserved < 0.2",
    "spearman_rho_with_true_f1": 0.58, "precision_for_true_f1_below_0.5": 0.61, "base_rate": 0.30, "recall": 0.49,
    "source": "the SCARLET documentation, Reading the output: qc.tsv",
}
NOT_ASSIGNED = {"unassigned", "Unassigned", "Unknown", "unknown", "rand", "Node", "ambiguous", "nan", "NA"}


NOT_STABLE_NOTE = "not refined: not anchored stably"
ORIGIN_NOT_REFINED_NOTE = "origin not refined: fewer than 2 cell types anchored stably"
UNDETERMINED = "undetermined across draws"
REASON_BACKGROUND = "background"
REASON_AMBIGUITY = "ambiguity"
REASON_FLIP = "flip test"
REASON_ORDER = (REASON_BACKGROUND, REASON_AMBIGUITY, REASON_FLIP)   # a tie between reasons goes to the first
AGREE_SHARE = 0.8   # a final call needs its outcome in at least 80 % of the draws (8 of 10: one-sided p 0.055)


def needed(n):
    """The draws a final call needs out of n: the share AGREE_SHARE of n, rounded up (8 of 10, 16 of 20, 1 of 1)."""
    return -(-4 * int(n) // 5)


def stable_anchors(unique_draws):
    """unique_draws: cells x N first-pass unique calls (UNASSIGNED where not unique), one column per draw of the run.
    Returns per cell the cell type it was called uniquely to in all N draws, UNASSIGNED elsewhere."""
    U_ = np.asarray(unique_draws, dtype=object).astype(str)
    first = U_[:, 0]
    same = (U_ == first[:, None]).all(axis=1) & (first != UNASSIGNED)
    return np.where(same, first, UNASSIGNED).astype(object)


def refine_plan(anchor, origin, first_best_draws, first_unique_draws, min_cells=MIN_CELLS):
    """Per origin (sorted): the stable anchors per cell type, the cell types refined (>= min_cells stable anchors),
    those not refined (any other cell type the first pass called in the origin, in any draw, or with stable anchors)
    and whether the origin is refined (>= 2 refined cell types). Returns the records (one per origin)."""
    anchor = np.asarray(anchor, dtype=object).astype(str)
    origin = np.asarray(origin).astype(str)
    fb = np.asarray(first_best_draws, dtype=object).astype(str)
    fu = np.asarray(first_unique_draws, dtype=object).astype(str)
    records = []
    for g in sorted(set(origin)):
        pos = np.flatnonzero(origin == g)
        a = anchor[pos]
        vc = pd.Series(a[a != UNASSIGNED]).value_counts()
        refined = sorted(vc[vc >= min_cells].index)
        seen = (set(np.unique(fb[pos])) | set(np.unique(fu[pos])) | set(vc.index)) - {UNASSIGNED}
        rec = dict(origin=str(g), n=int(len(pos)), stable_anchors={str(k): int(v) for k, v in vc.sort_index().items()},
                   refined=[str(k) for k in refined], not_refined=sorted(str(k) for k in seen - set(refined)),
                   n_train=int(np.isin(a, refined).sum()))
        if len(refined) < 2:
            rec["skipped"] = "fewer than 2 cell types anchored stably: first-pass calls kept"
        records.append(rec)
    return records


def keep_first(first_best, refined):
    """In one draw: the cells whose first-pass best call is a cell type that is not refined (they keep that draw's
    first-pass calls)."""
    fb = np.asarray(first_best, dtype=object).astype(str)
    return (fb != UNASSIGNED) & ~np.isin(fb, list(refined))


def pooled_background(model, data):
    """A copy of `model` whose classifier background is the detection rate over all the cells of `data` (the origin's
    query cells as the model scores them): the universe's bg column, and a natural background of N = the cells in the
    scoring statistics (Jeffreys floor 0.5 / (N + 1)). Returns (model, record)."""
    import dataclasses
    import io
    import json
    import scipy.sparse as sp
    uni = pd.read_csv(io.BytesIO(model.files["universe.tsv"]), sep="\t", index_col=0)
    genes = pd.Index(data.var_names.astype(str))
    col = genes.get_indexer(uni.index.astype(str))
    if (col < 0).any():
        raise ValueError(f"{int((col < 0).sum())} programme genes are missing from the query")
    n = int(data.n_obs)
    if n == 0:
        raise ValueError("no query cell to build the background from")
    X = data.X if sp.issparse(data.X) else sp.csr_matrix(data.X)
    X = X.tocsr()
    det = np.asarray((X[:, col] > 0).sum(axis=0)).ravel() / n
    K = np.asarray((X > 0).sum(axis=1)).ravel()
    uni["bg"] = det
    mfp = json.loads(model.files["marker_finder_params.json"])
    st = mfp.get("scoring_stats")
    if not st or "background" not in st:
        raise ValueError("the refinement model has no scoring statistics")
    st["background"] = {"mode": "natural", "N": float(n), "Kbar": float(K.mean())}
    files = dict(model.files)
    buf = io.StringIO()
    uni.to_csv(buf, sep="\t")
    files["universe.tsv"] = buf.getvalue().encode()
    files["marker_finder_params.json"] = json.dumps(mfp, indent=1).encode()
    if "in_class.tsv" in model.files:   # the full in-class table: the same background, for every gene
        inc = pd.read_csv(io.BytesIO(model.files["in_class.tsv"]), sep="\t", index_col=0)
        ci = genes.get_indexer(inc.index.astype(str))
        if (ci < 0).any():
            raise ValueError(f"{int((ci < 0).sum())} genes of the in-class table are missing from the query")
        inc["bg"] = np.asarray((X[:, ci] > 0).sum(axis=0)).ravel() / n
        buf = io.StringIO()
        inc.to_csv(buf, sep="\t")
        files["in_class.tsv"] = buf.getvalue().encode()
    return dataclasses.replace(model, files=files), dict(cells=n, mode="pooled", f_bg=0.5 / (n + 1.0))


def consensus(calls, reasons=None):
    """Per cell (row), the consensus of the draws (columns) of one output.

    calls: cells x N outcomes (a cell type, or UNASSIGNED). reasons: cells x N, each draw's refusal reason where its
    outcome is UNASSIGNED ('' elsewhere), or None. The winner is the most frequent outcome; it is kept when it holds
    at least needed(N) of the N draws (more than half, so a kept winner is never tied), else the cell is UNASSIGNED
    with the reason UNDETERMINED. A kept UNASSIGNED takes the most frequent reason among the refusing draws (ties:
    REASON_ORDER). Returns (label, agreement = the winner's count, reason)."""
    C = np.asarray(calls, dtype=object).astype(str)
    n_cells, N = C.shape
    vals, codes = np.unique(C, return_inverse=True)
    codes = codes.reshape(n_cells, N)
    counts = np.zeros((n_cells, len(vals)), dtype=np.int64)
    for j in range(len(vals)):
        counts[:, j] = (codes == j).sum(axis=1)
    order = np.argsort(-counts, axis=1, kind="stable")
    rows = np.arange(n_cells)
    kw_ = counts[rows, order[:, 0]]
    ok = kw_ >= needed(N)
    label = np.where(ok, vals[order[:, 0]], UNASSIGNED).astype(object)
    reason = np.where(ok, "", UNDETERMINED).astype(object)
    refused = ok & (label == UNASSIGNED)
    if refused.any():   # the most frequent reason among the refusing draws; ties: REASON_ORDER, then alphabetical
        R = np.full(C.shape, "", dtype=object) if reasons is None else np.asarray(reasons, dtype=object).astype(str)
        idx = np.flatnonzero(refused)
        Rr = R[idx].astype(str)
        ok = (C[idx] == UNASSIGNED) & (Rr != "")
        rv, rc = np.unique(Rr, return_inverse=True)
        rc = rc.reshape(Rr.shape)
        n_r, V = len(idx), len(rv)
        rows = np.repeat(np.arange(n_r), N).reshape(n_r, N)
        cnt = np.bincount((rows * V + rc)[ok], minlength=n_r * V).reshape(n_r, V)
        rank = np.array([REASON_ORDER.index(x) if x in REASON_ORDER else len(REASON_ORDER) + j
                         for j, x in enumerate(rv.tolist())])   # rv is sorted: alphabetical after REASON_ORDER
        mx = cnt.max(axis=1)
        cand = np.where((cnt == mx[:, None]) & (mx[:, None] > 0), rank[None, :], len(rank) + len(REASON_ORDER) + 1)
        pick = rv[cand.argmin(axis=1)].astype(object)
        reason[idx] = np.where(mx > 0, pick, REASON_BACKGROUND)
    return label, kw_.astype(int), reason


def _mode_rows(M):
    """Per row of M (cells x draws, strings; '' = none): the most frequent non-empty value (ties: alphabetical) and
    its count; ('', 0) where a row has none."""
    M = np.asarray(M, dtype=object)
    n = M.shape[0]
    out_v = np.full(n, "", dtype=object)
    out_n = np.zeros(n, dtype=np.int64)
    if n == 0 or M.size == 0:
        return out_v, out_n
    vals, codes = np.unique(M.astype(str), return_inverse=True)   # sorted: argmax takes the first of a tie
    codes = codes.reshape(M.shape)
    V = len(vals)
    cnt = np.bincount((np.arange(n)[:, None] * V + codes).ravel(), minlength=n * V).reshape(n, V)
    cnt[:, vals == ""] = 0
    j = cnt.argmax(axis=1)
    k = cnt[np.arange(n), j]
    has = k > 0
    out_v[has] = vals[j[has]].astype(object)
    out_n[has] = k[has]
    return out_v, out_n


def runner_up(best, top1, top2):
    """Per cell, from the draws' ranked identities (cells x N; scarlet_top1 / scarlet_top2 of each draw's assign):
    top1 = the most frequent top-ranked identity across the draws (ties: alphabetical); the runner-up = the most
    frequent, across the draws, of the highest-ranked identity other than the cell's reference identity (the
    output's best call where it holds one, else top1): the draw's top1 if it differs from it, else its top2.
    Returns (top1, runner_up, runner_up_draws); '' and 0 where there is none (a single identity)."""
    T1 = np.asarray(top1, dtype=object).astype(str)
    T2 = np.asarray(top2, dtype=object).astype(str)
    T1[np.isin(T1, ["None", "nan", UNASSIGNED])] = ""
    T2[np.isin(T2, ["None", "nan", UNASSIGNED])] = ""
    t1, _ = _mode_rows(T1)
    b = np.asarray(best, dtype=object).astype(str)
    ref = np.where(b != UNASSIGNED, b, t1).astype(str)
    alt = np.where(T1 != ref[:, None], T1, T2)
    alt[alt == ref[:, None]] = ""
    ru, n = _mode_rows(alt)
    return t1, ru, n


def draw_reasons(best, unique, flipped=None):
    """One draw's refusal reasons for its best and unique outputs: flip test where a flipped call was refused, else
    background where best is unassigned, else (unique only) ambiguity; '' where the output holds a call."""
    b = np.asarray(best, dtype=object).astype(str)
    u = np.asarray(unique, dtype=object).astype(str)
    fl = np.zeros(len(b), bool) if flipped is None else np.asarray(flipped, bool)
    rb = np.where(b == UNASSIGNED, np.where(fl, REASON_FLIP, REASON_BACKGROUND), "").astype(object)
    ru = np.where(u == UNASSIGNED, np.where(fl, REASON_FLIP, np.where(b == UNASSIGNED, REASON_BACKGROUND,
                                                                    REASON_AMBIGUITY)), "").astype(object)
    return rb, ru


def check_counts(X, what):
    """Depth matching reads detection and expression values: non-negative values are required
    (the classifier's own rule; a negative value stops it as non-count input)."""
    import scipy.sparse as sp
    from scarlet.input_checks import NonCountInputError
    v = X.data if sp.issparse(X) else np.asarray(X).ravel()
    if v.size and np.any(v < 0):
        raise NonCountInputError(f"depth matching needs non-negative expression values, but the {what} has negative "
                                 "values (corrected or scaled data?)")


def qc_classes(labels, min_cells=QC_MIN_CELLS):
    """The identities a query programme is learned for: transferred labels with >= min_cells cells (not
    unassigned)."""
    lab = np.asarray(labels).astype(str)
    vc = pd.Series(lab[~np.isin(lab, list(NOT_ASSIGNED))]).value_counts()
    return vc, vc[vc >= min_cells].index


def qc_rows(origin, ref_markers, query_markers, labels, min_cells=QC_MIN_CELLS):
    """Per identity: markers of the reference programme, of the query programme, the share conserved and the flag.

    ref_markers / query_markers: marker tables (columns cell_type, gene); query_markers None when fewer than 2
    identities had >= min_cells transferred cells. labels: the transferred labels of the query cells.
    """
    vc, ok = qc_classes(labels, min_cells)
    lab = np.asarray(labels).astype(str)
    idents = sorted(set(vc.index) | set(ref_markers["cell_type"].astype(str)))
    rows = []
    for c in idents:
        sa = set(ref_markers.gene[ref_markers.cell_type.astype(str) == c].astype(str))
        sb = set(query_markers.gene[query_markers.cell_type.astype(str) == c].astype(str)) \
            if query_markers is not None else set()
        inb = query_markers is not None and c in ok
        share = (len(sa & sb) / len(sa)) if (sa and inb) else np.nan
        rows.append(dict(origin=origin, identity=c, n_cells_transferred=int((lab == c).sum()), in_query_programme=inb,
                         reference_markers=len(sa), query_markers=len(sb), conserved=len(sa & sb),
                         share_conserved=share,
                         doubtful_transfer=(bool(share < QC_FLAG) if not np.isnan(share) else None)))
    return rows
