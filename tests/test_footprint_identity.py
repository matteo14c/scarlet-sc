"""Memory footprint (the verbs): each change keeps the same values, byte for byte, as the computation it replaces
(the replaced code is copied here as the reference)."""
from __future__ import annotations

import numpy as np
import pytest
from scipy import sparse


def _matrix(seed, n=60, g=40, density=0.4, dtype=np.int64, zeros=True):
    rng = np.random.default_rng(seed)
    X = sparse.random(n, g, density=density, format="csr", random_state=rng,
                      data_rvs=lambda k: rng.integers(1, 5000, k)).astype(dtype)
    if dtype != np.int64:
        X.data = X.data * rng.random(X.nnz).astype(dtype)
    X[3] = 0                      # an empty row
    X[7, :] = rng.integers(1, 9, g)   # a full (long) row
    X = sparse.csr_matrix(X)
    if zeros:                     # stored zeros
        X.data[::11] = 0
    else:                         # canonical counts: no stored zeros, sorted indices
        X.eliminate_zeros()
        X.sort_indices()
    return X


# ---- P1: row passes in chunks of whole rows capped by entries --------------------------------------------------
def _row_stats_rows(X, chunk=20000):
    """settings._row_stats before P1 (chunks of rows)"""
    X = sparse.csr_matrix(X)
    n = X.shape[0]
    detected = np.zeros(n)
    totals = np.zeros(n)
    for i in range(0, n, chunk):
        j = min(n, i + chunk)
        v = X.data[X.indptr[i]:X.indptr[j]]
        rows = np.repeat(np.arange(j - i), np.diff(X.indptr[i:j + 1]))
        detected[i:j] = np.bincount(rows, weights=(v > 0), minlength=j - i)
        totals[i:j] = np.bincount(rows, weights=v.astype(np.float64), minlength=j - i)
    return np.asarray(detected, dtype=float), np.asarray(totals, dtype=float)


def _row_counts_in_rows(X, mask, chunk_rows=20000):
    """depth_match._row_counts_in before P1"""
    X = sparse.csr_matrix(X)
    n = X.shape[0]
    out = np.zeros(n, dtype=np.int64)
    for i in range(0, n, chunk_rows):
        j = min(n, i + chunk_rows)
        a, b = X.indptr[i], X.indptr[j]
        ok = (X.data[a:b] > 0) & mask[X.indices[a:b]]
        rows = np.repeat(np.arange(j - i), np.diff(X.indptr[i:j + 1]))
        out[i:j] = np.bincount(rows[ok], minlength=j - i)
    return out


def _row_stats_backed_rows(path, n_obs, chunk=20000):
    """settings._row_stats_backed before P1 (csr files)"""
    import h5py
    with h5py.File(path, "r") as f:
        X = f["X"]
        indptr = X["indptr"][:]
        data = X["data"]
        detected = np.zeros(n_obs)
        totals = np.zeros(n_obs)
        for i in range(0, n_obs, chunk):
            j = min(n_obs, i + chunk)
            v = data[indptr[i]:indptr[j]].astype(np.float64)
            rows = np.repeat(np.arange(j - i), np.diff(indptr[i:j + 1]))
            detected[i:j] = np.bincount(rows, weights=(v > 0), minlength=j - i)
            totals[i:j] = np.bincount(rows, weights=v, minlength=j - i)
        return detected, totals


def _same(a, b):
    a, b = np.asarray(a), np.asarray(b)
    return a.dtype == b.dtype and a.shape == b.shape and a.tobytes() == b.tobytes()


@pytest.mark.parametrize("cap", [1, 7, 40, 100, 1 << 22])
@pytest.mark.parametrize("dtype", [np.int64, np.float32, np.float64])
def test_row_chunks_cover_whole_rows(cap, dtype):
    from scarlet import settings as S
    X = _matrix(1, dtype=dtype)
    ranges = list(S._row_chunks(X.indptr, 20000, cap))
    assert ranges[0][0] == 0 and ranges[-1][1] == X.shape[0]
    assert all(a[1] == b[0] for a, b in zip(ranges, ranges[1:]))
    for i, j in ranges:   # within the cap, or one row alone
        assert X.indptr[j] - X.indptr[i] <= cap or j == i + 1


@pytest.mark.parametrize("cap", [1, 7, 40, 100, 1 << 22])
@pytest.mark.parametrize("dtype", [np.int64, np.float32, np.float64])
def test_row_stats_same_as_row_chunks(monkeypatch, cap, dtype):
    from scarlet import depth_match as D, settings as S
    monkeypatch.setattr(S._row_chunks, "__defaults__", (20000, cap))
    for seed in range(3):
        X = _matrix(seed, dtype=dtype)
        for chunk in (5, 20000):
            new, old = S._row_stats(X, chunk), _row_stats_rows(X, chunk)
            assert _same(new[0], old[0]) and _same(new[1], old[1])
        mask = np.random.default_rng(seed).random(X.shape[1]) < 0.6
        assert _same(D._row_counts_in(X, mask), _row_counts_in_rows(X, mask))
        assert _same(D._row_counts_in(X, mask, 4), _row_counts_in_rows(X, mask, 4))


def test_row_stats_backed_same(tmp_path, monkeypatch):
    import anndata as ad
    from scarlet import settings as S
    monkeypatch.setattr(S._row_chunks, "__defaults__", (20000, 13))
    X = _matrix(5)
    A = ad.AnnData(X)
    p = tmp_path / "a.h5ad"
    A.write_h5ad(p)
    new = S._row_stats_backed(str(p), A.n_obs)
    old = _row_stats_backed_rows(str(p), A.n_obs)
    assert _same(new[0], old[0]) and _same(new[1], old[1])


# ---- P6: the module globals of the marker search and of the scoring released after their step -----------------
def test_globals_released_after_learn_and_assign(tmp_path):
    from conftest import make_synthetic_adata
    from scarlet import api, classifier, marker_finder
    a = make_synthetic_adata(seed=2)
    m = tmp_path / "m.scarlet"
    api.learn(a, "cell_type", str(m), blocklist=None)
    assert marker_finder._BOOTSTRAP_X_BIN is None and marker_finder._BOOTSTRAP_PARAMS is None
    q = a.copy()
    q.obs = q.obs[[]]
    api.assign(str(m), q, str(tmp_path / "q.h5ad"))
    assert classifier._SCORING_X is None and classifier._SCORING_DELTA_LOG_P is None


# ---- P7: detection patterns take their values from one shared read-only all-True buffer ------------------------
def _positive_pattern_own(X):
    """marker_finder.positive_pattern before P7"""
    if sparse.issparse(X) and X.format == "csr" and X.has_sorted_indices and X.nnz and bool((X.data > 0).all()):
        return sparse.csr_matrix((np.ones(X.nnz, dtype=bool), X.indices, X.indptr), shape=X.shape, copy=False)
    return X > 0


def _thin_draw_own(prep, seed):
    """depth_match.thin_draw before P7"""
    from scarlet import depth_match as D
    rec = dict(prep["rec"])
    if rec["direction"] == "same":
        return None, rec
    X, p = prep["X"], prep["p"]
    rng = np.random.default_rng(seed)
    n = X.shape[0]
    counts = np.zeros(n, dtype=np.int64)
    idx_parts, dat_parts = [], []
    rows_per = max(1, int(D._CHUNK // max(1, X.nnz / max(1, n))))
    for i in range(0, n, rows_per):
        j = min(n, i + rows_per)
        a, b = X.indptr[i], X.indptr[j]
        keep = rng.random(b - a) < p[X.indices[a:b]]
        rows = np.repeat(np.arange(j - i), np.diff(X.indptr[i:j + 1]))
        counts[i:j] = np.bincount(rows[keep], minlength=j - i)
        idx_parts.append(X.indices[a:b][keep])
        dat_parts.append(X.data[a:b][keep])
    indptr = np.concatenate([[0], np.cumsum(counts)])
    Y = sparse.csr_matrix((np.concatenate(dat_parts) if dat_parts else X.data[:0],
                           np.concatenate(idx_parts) if idx_parts else X.indices[:0], indptr), shape=X.shape)
    rec.update(thinned=prep["side"], **prep["r"], detections_after=int(Y.nnz),
               depth_after=float(np.median(D._row_counts_in(Y, prep["sh"]))))
    return Y, rec


def _same_csr(A, B):
    return (A.shape == B.shape and A.format == B.format and _same(A.data, B.data) and _same(A.indices, B.indices)
            and _same(A.indptr, B.indptr))


@pytest.mark.parametrize("dtype", [np.int64, np.float32, bool])
def test_positive_pattern_shared_buffer(dtype):
    from scarlet import _patterns, marker_finder as M
    for seed in range(3):
        X = _matrix(seed, dtype=np.int64 if dtype is bool else dtype, zeros=False)
        if dtype is bool:
            X = X > 0
        new, old = M.positive_pattern(X), _positive_pattern_own(X)
        assert _same_csr(new, old)
        assert not new.data.flags.writeable and _patterns.shares_buffer(new.data)
        assert _same_csr(M.positive_pattern(new), old)   # a pattern of a pattern: the same matrix
    Z = _matrix(7)   # stored zeros: X > 0, as before
    assert _same_csr(M.positive_pattern(Z), _positive_pattern_own(Z))


@pytest.mark.parametrize("values", [True, False])
def test_thin_draw_shared_buffer(values):
    from scarlet import _patterns, depth_match as D, marker_finder as M
    R = _matrix(11, n=80, g=50, density=0.6, zeros=False)
    Q = _matrix(12, n=70, g=50, density=0.2, zeros=False)
    prep = D.prepare_depth(R, Q, "capture")
    assert prep["rec"]["direction"] != "same"
    if not values:
        prep["X"] = M.positive_pattern(prep["X"])
    for seed in (1, 2, 99):
        (Y, r), (Y0, r0) = D.thin_draw(prep, seed), _thin_draw_own(prep, seed)
        assert _same_csr(Y, Y0) and r == r0
        if not values:
            assert not Y.data.flags.writeable and _patterns.shares_buffer(Y.data)


def test_shared_buffer_kept_by_scipy():
    """SciPy keeps every slice of the shared buffer as the pattern's value array (no copy), whatever the lengths."""
    from scarlet import _patterns
    for n in (1, 2, 3, 5, 100, 129, 1000, 5442, 18, 65, 64, 4097, 3000):
        d = _patterns.all_true(n)
        Z = sparse.csr_matrix((d, np.arange(n, dtype=np.int32) % 7, np.array([0, n], dtype=np.int32)), shape=(1, 7),
                              copy=False)
        assert not Z.data.flags.writeable and _patterns.shares_buffer(Z.data), n


# ---- P3: int64 counts read as int32 (scarlet._read) -------------------------------------------------------------
def _write(tmp_path, X, name, **kw):
    import anndata as ad
    import pandas as pd
    n, g = X.shape
    A = ad.AnnData(X, obs=pd.DataFrame({"cell_type": pd.Categorical(np.array(["a", "b", "c"])[np.arange(n) % 3]),
                                        "depth": np.arange(n, dtype=float)}, index=[f"c{i}" for i in range(n)]),
                   var=pd.DataFrame(index=[f"G{i}" for i in range(g)]), uns={"note": "x", "nested": {"k": [1, 2]}},
                   **kw)
    p = tmp_path / f"{name}.h5ad"
    A.write_h5ad(p)
    return p


def _same_adata_but_x(A, B):
    import anndata as ad
    assert A.obs.equals(B.obs) and A.var.equals(B.var) and A.obs_names.equals(B.obs_names)
    assert list(A.obs.dtypes) == list(B.obs.dtypes)
    C, D = A.copy(), B.copy()
    C.X = D.X = None
    ad.tests.helpers.assert_equal(C, D) if hasattr(ad, "tests") else None
    assert repr(A.uns) == repr(B.uns) and list(A.obsm) == list(B.obsm) and list(A.layers) == list(B.layers)


@pytest.mark.parametrize("fmt", ["csr", "csc"])
def test_read_int32_same_values(tmp_path, fmt):
    import anndata as ad
    from scarlet._read import read_h5ad
    X = _matrix(3, zeros=False).asformat(fmt)
    X.data[0] = 2 ** 31 - 1   # the int32 limit still fits
    p = _write(tmp_path, X, fmt, obsm={"emb": np.ones((X.shape[0], 2))})
    A, B = ad.read_h5ad(p), read_h5ad(p)
    assert A.X.dtype == np.int64 and B.X.dtype == np.int32 and type(A.X) is type(B.X) and A.X.format == B.X.format
    assert np.array_equal(A.X.data, B.X.data) and _same(A.X.indices, B.X.indices) and _same(A.X.indptr, B.X.indptr)
    assert A.X.has_sorted_indices == B.X.has_sorted_indices
    _same_adata_but_x(A, B)


@pytest.mark.parametrize("fmt", ["csr", "csc"])
def test_read_int32_every_field_as_anndata(tmp_path, fmt):
    """Every element but X's values as anndata.read_h5ad reads it: obs, var, obsm, varm, obsp, varp, layers, uns."""
    import anndata as ad
    import scipy.sparse as sp
    from scarlet._read import read_h5ad
    X = _matrix(5, zeros=False).asformat(fmt)
    n, g = X.shape
    p = _write(tmp_path, X, f"all_{fmt}", obsm={"emb": np.ones((n, 2))}, varm={"load": np.arange(g * 2.0).reshape(g, 2)},
               obsp={"conn": sp.identity(n, format="csr")}, varp={"gg": sp.identity(g, format="csr")},
               layers={"counts": X.copy()})
    A, B = ad.read_h5ad(p), read_h5ad(p)
    assert B.X.dtype == np.int32 and type(A.X) is type(B.X) and np.array_equal(A.X.data, B.X.data)
    _same_adata_but_x(A, B)
    for k in ("obsm", "varm", "obsp", "varp", "layers"):
        assert list(getattr(A, k)) == list(getattr(B, k)), k   # in order (anndata 0.13 lists X as layers[None])


def test_read_uses_only_public_anndata():
    """The reader follows anndata's releases: no private anndata module, no experimental dispatch."""
    import inspect
    import scarlet._read as R
    src = inspect.getsource(R)
    for bad in ("anndata._", "read_dispatched", "get_memory_class"):
        assert bad not in src, bad


def test_read_int32_falls_back(tmp_path):
    import anndata as ad
    from scarlet._read import read_h5ad
    big = _matrix(4, zeros=False)
    big.data[5] = 2 ** 31   # outside int32: as stored
    cases = {"big": big, "float": _matrix(4, dtype=np.float32, zeros=False), "dense": _matrix(4, zeros=False).toarray(),
             "int16": _matrix(4, zeros=False).astype(np.int16)}
    for name, X in cases.items():
        p = _write(tmp_path, X, name)
        A, B = ad.read_h5ad(p), read_h5ad(p)
        assert type(A.X) is type(B.X) and A.X.dtype == B.X.dtype, name
        _same_adata_but_x(A, B)
    A0 = ad.AnnData(_matrix(4, zeros=False), obs=ad.read_h5ad(_write(tmp_path, big, "o")).obs)   # with raw: as stored
    A0.raw = A0.copy()
    A0.write_h5ad(tmp_path / "raw.h5ad")
    assert read_h5ad(tmp_path / "raw.h5ad").X.dtype == np.int64


def test_int32_counts_same_summaries():
    """the steps that compute with the values: the same floats from int32 as from int64"""
    from scarlet import depth_match as D, settings as S
    X = _matrix(9, n=90, g=60, density=0.5, zeros=False)
    X.data = X.data * 40   # read-count sized values
    X32 = X.astype(np.int32)
    for a, b in zip(S._row_stats(X), S._row_stats(X32)):
        assert _same(a, b)
    Q = _matrix(10, n=70, g=60, density=0.2, zeros=False)
    p64, p32 = D.prepare_depth(X, Q, "capture"), D.prepare_depth(X32, Q.astype(np.int32), "capture")
    assert p64["rec"] == p32["rec"] and p64["r"] == p32["r"] and _same(p64["p"], p32["p"])
    assert _same(D._log_means(p64["X"]), D._log_means(p32["X"]))
    for seed in (1, 5):
        (Y, r), (Y2, r2) = D.thin_draw(p64, seed), D.thin_draw(p32, seed)
        assert r == r2 and _same(Y.indices, Y2.indices) and _same(Y.indptr, Y2.indptr)
        assert np.array_equal(Y.data, Y2.data)


# ---- P5: learn keeps the detection only of an input it read --------------------------------------------------------
def test_learn_from_file_same_model_and_caller_unchanged(tmp_path):
    import json
    import zipfile
    from conftest import make_synthetic_adata
    from scarlet import api
    a = make_synthetic_adata(seed=4)
    a.X = sparse.csr_matrix(a.X).astype(np.int64)
    before = a.X.copy()
    p = tmp_path / "ref.h5ad"
    a.write_h5ad(p)
    api.learn(str(p), "cell_type", str(tmp_path / "f.scarlet"), blocklist=None)
    api.learn(a, "cell_type", str(tmp_path / "m.scarlet"), blocklist=None)
    assert a.X.dtype == np.int64 and (a.X != before).nnz == 0   # the caller's object keeps its values
    zf, zm = zipfile.ZipFile(tmp_path / "f.scarlet"), zipfile.ZipFile(tmp_path / "m.scarlet")
    assert zf.namelist() == zm.namelist()
    differ = [n for n in zf.namelist() if zf.read(n) != zm.read(n)]
    for n in differ:   # only what names the input (its path, or "in-memory AnnData"; the model path)
        assert n in ("marker_finder_params.json", "provenance.json", "settings.json", "summary.txt",
                     "manifest.json"), n   # (the manifest holds the members' checksums)
    mp = [json.loads(z.read("marker_finder_params.json")) for z in (zf, zm)]
    for d in mp:
        d["marker_finder"].pop("input")
        d["context"].pop("marker_finder_argv")
        for k in ("output_markers", "output_universe", "output_background_universe"):   # temporary folders
            d["marker_finder"].pop(k, None)
    assert mp[0]["marker_finder"] == mp[1]["marker_finder"]
    for n in ("markers.tsv", "universe.tsv", "background_universe.tsv", "programmes.json", "positive_evidence.json"):
        if n in zf.namelist():
            assert zf.read(n) == zm.read(n), n
