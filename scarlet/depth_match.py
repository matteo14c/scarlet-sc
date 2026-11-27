"""Depth matching in transfer, between a reference side and a query side of the same genes.

Depth of a side = the median, over its cells, of the genes detected per cell on the shared genes (genes detected in
at least one cell of each side), on the matrix the cells are scored on. The deeper side is thinned towards the
shallower one; no threshold (equal medians: nothing is thinned). The thinning (capture), seeded: across the shared
genes, an isotonic (non-decreasing) curve f from the gene's mean log expression on the deep side (mean over the deep
cells of log1p(value), zeros included) to its detection rate on the shallow side (pool-adjacent-violators,
unweighted; linear between fitted points, constant outside); on the deep side each detection of gene g is kept with
probability min(1, f(m_g) / d_g), d_g = g's detection rate on the deep side. Label-free.

Kept detections keep their values, so counts stay counts.
"""
from __future__ import annotations

import numpy as np
from scipy import sparse

METHODS = ("capture",)


def _csr(X):
    return sparse.csr_matrix(X, dtype=np.float64) if not sparse.issparse(X) else X.tocsr()


def shared_genes(X_ref, X_query):
    """Genes detected (> 0) in at least one cell of each side."""
    det_r = np.asarray((_csr(X_ref) > 0).sum(axis=0)).ravel() > 0
    det_q = np.asarray((_csr(X_query) > 0).sum(axis=0)).ravel() > 0
    return det_r & det_q


def genes_per_cell(X, mask):
    """Genes detected per cell on the genes in mask."""
    X = _csr(X)
    return np.asarray((X[:, np.flatnonzero(mask)] > 0).sum(axis=1)).ravel().astype(np.int64)


def direction(depth_ref, depth_query):
    if depth_ref == depth_query:
        return "same"
    return "deep -> shallow" if depth_ref > depth_query else "shallow -> deep"


def isotonic_fit(x, y):
    """Non-decreasing least-squares fit of y on x (pool-adjacent-violators, equal weights; ties in x pooled first).
    Returns (knots, values): the fitted function is linear between knots and constant outside."""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    o = np.argsort(x, kind="mergesort")
    xs, ys = x[o], y[o]
    ux, start = np.unique(xs, return_index=True)
    cnt = np.diff(np.r_[start, len(xs)])
    sums = np.add.reduceat(ys, start)
    # PAV over the tied-x blocks: (sum, weight) per block, merged while decreasing
    bs, bw, bx0 = [], [], []
    for i in range(len(ux)):
        s, w = float(sums[i]), float(cnt[i])
        bs.append(s)
        bw.append(w)
        bx0.append(i)
        while len(bs) > 1 and bs[-2] / bw[-2] > bs[-1] / bw[-1]:
            s1, w1 = bs.pop(), bw.pop()
            bx0.pop()
            bs[-1] += s1
            bw[-1] += w1
    vals = np.empty(len(ux))
    ends = bx0[1:] + [len(ux)]
    for s, w, a, b in zip(bs, bw, bx0, ends):
        vals[a:b] = s / w
    return ux, vals


def _apply_fit(knots, vals, x):
    return np.interp(np.asarray(x, dtype=float), knots, vals)


def _keep(X, keep_mask_data):
    """X (CSR) with only the stored entries where keep_mask_data is True."""
    X = X.tocsr()
    rows = np.repeat(np.arange(X.shape[0]), np.diff(X.indptr))
    k = np.asarray(keep_mask_data, dtype=bool)
    out = sparse.csr_matrix((X.data[k], (rows[k], X.indices[k])), shape=X.shape)
    out.sort_indices()
    return out


def _detected_entries(X):
    X = _csr(X).copy()
    X.eliminate_zeros()
    X.sort_indices()
    return X


def thin_capture(X_deep, X_shallow, shared, seed):
    """Capture-informed thinning. Returns (thinned X, record)."""
    X = _detected_entries(X_deep)
    S = _detected_entries(X_shallow)
    rng = np.random.default_rng(seed)
    n_d, n_s = X.shape[0], S.shape[0]
    d = np.bincount(X.indices, minlength=X.shape[1]) / n_d
    s = np.bincount(S.indices, minlength=S.shape[1]) / n_s
    lv = X.copy()
    lv.data = np.log1p(np.maximum(lv.data, 0))
    m = np.asarray(lv.sum(axis=0)).ravel() / n_d
    knots, vals = isotonic_fit(m[shared], s[shared])
    f = _apply_fit(knots, vals, m)
    p = np.where(d > 0, np.minimum(1.0, f / np.where(d > 0, d, 1.0)), 1.0)
    keep = rng.random(X.nnz) < p[X.indices]
    Y = _keep(X, keep)
    sh = shared & (d > 0)
    rec = dict(method="capture", genes_fitted=int(shared.sum()), knots=int(len(knots)),
               keep_prob_median_shared=float(np.median(p[sh])) if sh.any() else None,
               keep_prob_mean_detection=float(p[X.indices].mean()) if X.nnz else None,
               genes_thinned=int((p[d > 0] < 1).sum()), detections_before=int(X.nnz), detections_after=int(Y.nnz))
    return Y, rec


def match_depth(X_ref, X_query, method, seed):
    """Thin the deeper of the two sides towards the shallower one. Returns (X_ref', X_query', record); the side
    that is not thinned is returned as given."""
    if method not in METHODS:
        raise ValueError(f"depth matching: unknown method {method!r}")
    sh = shared_genes(X_ref, X_query)
    kr, kq = genes_per_cell(X_ref, sh), genes_per_cell(X_query, sh)
    dr, dq = float(np.median(kr)), float(np.median(kq))
    rec = dict(shared_genes=int(sh.sum()), depth_ref=dr, depth_query=dq,
               ratio=(dr / dq) if dq > 0 else None, direction=direction(dr, dq), thinned="none")
    if rec["direction"] == "same":
        return X_ref, X_query, rec
    if dr > dq:
        deep, shallow, side = X_ref, X_query, "reference"
    else:
        deep, shallow, side = X_query, X_ref, "query"
    Y, r = thin_capture(deep, shallow, sh, seed)
    rec.update(thinned=side, **r, depth_after=float(np.median(genes_per_cell(Y, sh))))
    return (Y, X_query, rec) if side == "reference" else (X_ref, Y, rec)


# ---------------------------------------------------------------------------------------------------------------------
# The same depth matching, with its seed-independent part computed once per run (prepare) and only the keep draw per
# draw (thin_draw): prepare + thin_draw(seed) give the matrix and the record that match_depth(X_ref, X_query, method,
# seed) gives, bit for bit. Counts are taken from the stored values without copying the matrices, the deep side's
# detected entries are its own matrix when it is already in canonical form without stored zeros (else the copy
# match_depth makes), and the kept entries are written straight into CSR (no COO, no row index of every entry).
# ---------------------------------------------------------------------------------------------------------------------
_CHUNK = 8192 * 512   # entries per chunk: a multiple of NumPy's 8192-element buffer, so casts see the same blocks


def _col_counts(X, nonzero=False):
    """Per column, the stored values > 0 (nonzero=True: != 0), counted in chunks of entries."""
    X = _csr(X)
    out = np.zeros(X.shape[1], dtype=np.int64)
    for s in range(0, X.nnz, _CHUNK):
        v = X.data[s:s + _CHUNK]
        ok = (v != 0) if nonzero else (v > 0)
        out += np.bincount(X.indices[s:s + _CHUNK][ok], minlength=X.shape[1])
    return out


def _row_counts_in(X, mask, chunk_rows=20000):
    """genes_per_cell(X, mask) (per row, the stored values > 0 in the columns of mask), in chunks of whole rows capped
    by rows and by stored entries (settings._row_chunks: integer counts per row, whatever the chunks)."""
    from scarlet.settings import _row_chunks
    X = _csr(X)
    n = X.shape[0]
    out = np.zeros(n, dtype=np.int64)
    for i, j in _row_chunks(X.indptr, chunk_rows):
        a, b = X.indptr[i], X.indptr[j]
        ok = (X.data[a:b] > 0) & mask[X.indices[a:b]]
        rows = np.repeat(np.arange(j - i), np.diff(X.indptr[i:j + 1]))
        out[i:j] = np.bincount(rows[ok], minlength=j - i)
    return out


def _detected(X):
    """_detected_entries(X): X itself when it is CSR in canonical form (sorted, no duplicates) without stored zeros,
    else the copy _detected_entries makes."""
    X = _csr(X)
    if X.has_canonical_format and (X.nnz == 0 or bool(np.all(X.data != 0))):
        return X
    return _detected_entries(X)


def prepare_depth(X_ref, X_query, method, cache=None):
    """The seed-independent part of match_depth(X_ref, X_query, method, seed). cache (a dict, one per run): the
    reference's detected entries and column statistics, shared by the origins it is matched against."""
    if method not in METHODS:
        raise ValueError(f"depth matching: unknown method {method!r}")
    cache = {} if cache is None else cache
    kr_key = ("ref", id(X_ref))
    if kr_key not in cache:
        cache[kr_key] = {"pos": _col_counts(X_ref) > 0}
    cr = cache[kr_key]
    sh = cr["pos"] & (_col_counts(X_query) > 0)
    kr, kq = _row_counts_in(X_ref, sh), _row_counts_in(X_query, sh)
    dr, dq = float(np.median(kr)), float(np.median(kq))
    rec = dict(shared_genes=int(sh.sum()), depth_ref=dr, depth_query=dq,
               ratio=(dr / dq) if dq > 0 else None, direction=direction(dr, dq), thinned="none")
    prep = {"rec": rec, "sh": sh}
    if rec["direction"] == "same":
        return prep
    if dr > dq:
        deep, shallow, side = X_ref, X_query, "reference"
    else:
        deep, shallow, side = X_query, X_ref, "query"
    if side == "reference" and "det" in cr:
        X, d_cnt, m = cr["det"], cr["d_cnt"], cr["m"]
    else:
        X = _detected(deep)
        d_cnt = np.bincount(X.indices, minlength=X.shape[1])
        m = _log_means(X)
        if side == "reference":
            cr.update(det=X, d_cnt=d_cnt, m=m)
    n_d, n_s = X.shape[0], _csr(shallow).shape[0]
    d = d_cnt / n_d
    s = _col_counts(shallow, nonzero=True) / n_s
    knots, vals = isotonic_fit(m[sh], s[sh])
    f = _apply_fit(knots, vals, m)
    p = np.where(d > 0, np.minimum(1.0, f / np.where(d > 0, d, 1.0)), 1.0)
    shd = sh & (d > 0)
    r = dict(method="capture", genes_fitted=int(sh.sum()), knots=int(len(knots)),
             keep_prob_median_shared=float(np.median(p[shd])) if shd.any() else None,
             keep_prob_mean_detection=float(p[X.indices].mean()) if X.nnz else None,
             genes_thinned=int((p[d > 0] < 1).sum()), detections_before=int(X.nnz))
    prep.update(side=side, X=X, p=p, r=r)
    return prep


def _log_means(X):
    """thin_capture's m: the column sums of log1p(max(value, 0)) over the detected entries, divided by the rows; the
    same SciPy sum on a matrix sharing X's index arrays (the log values made in chunks of whole 8192-element blocks,
    in the dtype thin_capture gets)."""
    lv = np.empty(X.nnz, dtype=np.log1p(np.maximum(X.data[:0], 0)).dtype)   # thin_capture's dtype
    for s in range(0, X.nnz, _CHUNK):
        np.log1p(np.maximum(X.data[s:s + _CHUNK], 0), out=lv[s:s + _CHUNK])
    L = sparse.csr_matrix((lv, X.indices, X.indptr), shape=X.shape, copy=False)
    return np.asarray(L.sum(axis=0)).ravel() / X.shape[0]


def thin_draw(prep, seed):
    """One draw of match_depth from prepare_depth's result: (the thinned matrix or None, the record)."""
    rec = dict(prep["rec"])
    if rec["direction"] == "same":
        return None, rec
    X, p = prep["X"], prep["p"]
    rng = np.random.default_rng(seed)
    n = X.shape[0]
    counts = np.zeros(n, dtype=np.int64)
    idx_parts, dat_parts = [], []
    # a detection pattern (all values True): the kept values are all True, taken from the shared read-only buffer
    # (scarlet._patterns) instead of copied
    pattern = X.dtype == np.bool_ and bool(X.data.all())
    rows_per = max(1, int(_CHUNK // max(1, X.nnz / max(1, n))))
    for i in range(0, n, rows_per):   # whole rows per chunk: the uniforms are drawn in entry order, as one array
        j = min(n, i + rows_per)
        a, b = X.indptr[i], X.indptr[j]
        keep = rng.random(b - a) < p[X.indices[a:b]]
        rows = np.repeat(np.arange(j - i), np.diff(X.indptr[i:j + 1]))
        counts[i:j] = np.bincount(rows[keep], minlength=j - i)
        idx_parts.append(X.indices[a:b][keep])
        if not pattern:
            dat_parts.append(X.data[a:b][keep])
    indptr = np.concatenate([[0], np.cumsum(counts)])
    if pattern:
        from scarlet._patterns import all_true
        data = all_true(int(indptr[-1]))
    else:
        data = np.concatenate(dat_parts) if dat_parts else X.data[:0]
    Y = sparse.csr_matrix((data, np.concatenate(idx_parts) if idx_parts else X.indices[:0], indptr), shape=X.shape)
    rec.update(thinned=prep["side"], **prep["r"], detections_after=int(Y.nnz),
               depth_after=float(np.median(_row_counts_in(Y, prep["sh"]))))
    return Y, rec
