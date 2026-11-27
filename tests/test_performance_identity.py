"""The faster paths give the same bits as the computations they replace."""
import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp


def _counts(n, g, density, seed, zeros=False):
    rng = np.random.default_rng(seed)
    X = sp.random(n, g, density=density, format="csr", random_state=seed,
                  data_rvs=lambda k: rng.integers(1, 50, size=k)).astype(np.int64)
    if zeros:   # stored zeros
        X.data[::7] = 0
    return X


@pytest.mark.parametrize("zeros", [False, True])
def test_detected_per_cell_equals_row_stats(zeros):
    from scarlet import settings
    X = _counts(300, 80, 0.2, 1, zeros)
    a, b = settings.detected_per_cell(X), settings._row_stats(X)[0]
    assert a.dtype == b.dtype and a.tobytes() == b.tobytes()
    P = X > 0
    assert settings.detected_per_cell(P).tobytes() == settings._row_stats(P)[0].tobytes()
    D = X.toarray()
    assert settings.detected_per_cell(D).tobytes() == settings._row_stats(D)[0].tobytes()


@pytest.mark.parametrize("zeros", [False, True])
def test_marker_counts_equal_detection_then_columns(zeros):
    from scarlet import evidence
    import anndata as ad
    X = _counts(200, 60, 0.3, 2, zeros)
    var = pd.Index([f"g{i}" for i in range(60)])
    lists = {"a": ["g1", "g5", "g7"], "b": ["g5", "g30", "g59"], "c": []}
    ids = ["a", "b", "c"]
    old = evidence._marker_matrix(evidence._detection(ad.AnnData(X)), var, lists, ids)
    new = evidence._marker_counts(X, var, lists, ids)
    assert old.dtype == new.dtype and old.tobytes() == new.tobytes()


def test_rows_shares_X_and_copies_the_rest():
    import anndata as ad
    from scarlet import api
    X = _counts(50, 20, 0.3, 3)
    A = ad.AnnData(X, obs=pd.DataFrame({"o": ["x"] * 50}, index=[f"c{i}" for i in range(50)]))
    S = api._rows(A, np.arange(50))
    assert S.X is A.X and A.X is X
    S.obs["new"] = 1
    assert "new" not in A.obs
    T = api._rows(A, np.arange(10))
    assert T.X is not A.X and (T.X != X[:10]).nnz == 0


def test_memo_returns_what_the_functions_compute():
    from scarlet import blocklists as B
    names = ["CD3E", "Cd3e", "ENSG00000000003", "MT-CO1"] * 10
    assert B.symbol_case_species(names) == B._symbol_case_species(names)
    assert B.symbol_case_species(names) == B.symbol_case_species(list(names))
    assert B.gene_id_type(names) == B._gene_id_type(names)
    p = B.bundled_path("human")
    assert B.file_md5(p) == B._file_md5(p)
    assert B.blocklist_identity(p) == B._blocklist_identity(p)
    assert B.read_table(p).equals(pd.read_csv(p, sep="\t"))


def _pair(deep_ref, dtype, unsorted=False):
    rng = np.random.default_rng(5)
    def m(n, dens, seed):
        X = sp.random(n, 120, density=dens, format="csr", random_state=seed,
                      data_rvs=lambda k: rng.integers(1, 30, size=k)).astype(dtype)
        return X
    R, Q = (m(150, 0.5, 6), m(90, 0.15, 7)) if deep_ref else (m(150, 0.15, 6), m(90, 0.5, 7))
    if unsorted:   # a stored zero and reversed rows: match_depth copies and sorts, the fast path must too
        R = R.tolil().tocsr()
        R.data[3] = 0
        R.indices[:R.indptr[1]] = R.indices[:R.indptr[1]][::-1].copy()
        R.data[:R.indptr[1]] = R.data[:R.indptr[1]][::-1].copy()
        R.has_sorted_indices = False
    return R, Q


@pytest.mark.parametrize("deep_ref", [True, False])
@pytest.mark.parametrize("dtype", [np.int64, np.float32, np.float64])
@pytest.mark.parametrize("unsorted", [False, True])
def test_prepared_depth_matching_equals_match_depth(deep_ref, dtype, unsorted):
    from scarlet import depth_match as DM
    R, Q = _pair(deep_ref, dtype, unsorted)
    prep = DM.prepare_depth(R, Q, "capture", {})
    for seed in (1, 2, 12345):
        Xr, Xq, rec0 = DM.match_depth(R.copy(), Q.copy(), "capture", seed)
        Y, rec1 = DM.thin_draw(prep, seed)
        assert list(rec0) == list(rec1) and rec0 == rec1
        old = Xr if rec0["thinned"] == "reference" else Xq
        assert old.dtype == Y.dtype and old.shape == Y.shape
        for a in ("data", "indices", "indptr"):
            assert np.array_equal(getattr(old, a), getattr(Y, a)), a


def test_prepared_depth_matching_same_depth():
    from scarlet import depth_match as DM
    R, Q = _pair(True, np.int64)
    prep = DM.prepare_depth(R, R, "capture", {})
    Y, rec = DM.thin_draw(prep, 3)
    assert Y is None and rec == DM.match_depth(R, R, "capture", 3)[2]


def test_uniforms_in_chunks_are_the_same_stream():
    a = np.random.default_rng(9).random(100003)
    g = np.random.default_rng(9)
    b = np.concatenate([g.random(n) for n in (7, 8192, 50000, 41804)])
    assert a.tobytes() == b.tobytes()


def test_learn_in_memory_equals_learn_from_file(synthetic_h5ad, tmp_path):
    import anndata as ad
    import zipfile
    from scarlet import api
    A = ad.read_h5ad(synthetic_h5ad)
    lab = "cell_type" if "cell_type" in A.obs else A.obs.columns[0]
    api.learn(str(synthetic_h5ad), lab, str(tmp_path / "f.scarlet"), blocklist="none")
    api.learn(A, lab, str(tmp_path / "m.scarlet"), blocklist="none")
    z0, z1 = zipfile.ZipFile(tmp_path / "f.scarlet"), zipfile.ZipFile(tmp_path / "m.scarlet")
    for n in ("markers.tsv", "universe.tsv", "background_universe.tsv", "programmes.json", "positive_evidence.json"):
        if n in z0.namelist():
            assert z0.read(n) == z1.read(n), n


def _bootstrap_sums(X, groups, types, p_u, fast, top_k=30):
    from types import SimpleNamespace
    from scarlet import marker_finder as M
    args = SimpleNamespace(alpha=1.5, beta=1.0, gamma=6.0, eps=1e-4, n_cells=40, top_k=top_k, n_iterations=17,
                           active_genes_seed=42, bootstrap_order="keyed", ncpu=1)
    out = [np.zeros((X.shape[1], len(types))), np.zeros((X.shape[1], len(types))),
           np.zeros((X.shape[1], len(types)), dtype=np.int64)]
    real = M.bootstrap_iteration
    if not fast:
        M.bootstrap_iteration = M._bootstrap_iteration_reference
    try:
        for scores, p_in, presence, _, _ in M.run_bootstrap_tasks(X, groups, types, p_u, args, False):
            out[0] += scores
            out[1] += p_in
            out[2] += presence
    finally:
        M.bootstrap_iteration = real
    return out


@pytest.mark.parametrize("ties", [False, True])
@pytest.mark.parametrize("top_k", [5, 30, 400])
def test_fast_bootstrap_equals_the_reference(ties, top_k):
    rng = np.random.default_rng(1)
    dens = rng.uniform(0.02, 0.6, 400)
    if ties:   # blocks of identical genes: equal scores straddle the top-k boundary
        dens = np.repeat(dens[:40], 10)
    D = rng.random((300, 400)) < dens
    if ties:
        D = np.repeat(D[:, ::10], 10, axis=1)
    X = sp.csr_matrix(D)
    types = ["a", "b", "c", "d"]
    groups = {t: np.arange(i, 300, 4) for i, t in enumerate(types)}
    p_u = np.asarray(X.mean(axis=0)).ravel()
    for a, b in zip(_bootstrap_sums(X, groups, types, p_u, True, top_k),
                    _bootstrap_sums(X, groups, types, p_u, False, top_k)):
        assert a.tobytes() == b.tobytes()


def _old_multilabel(ll, ll_bg, cats, fdr, alpha):
    from scipy.stats import chi2
    srt = np.argsort(ll, axis=1)
    best, alls, uniq = [], [], []
    for i in range(ll.shape[0]):
        ci = srt[i, :][::-1]
        cs = ll[i, ci]
        if ll_bg is not None and ll_bg[i] >= cs[0]:
            best.append("unassigned"); alls.append("unassigned"); uniq.append(False)
            continue
        best.append(cats[ci[0]])
        inc = [ci[0]]
        for j in range(1, len(ci)):
            if chi2.sf(2 * (cs[0] - cs[j]), df=1) >= alpha:
                inc.append(ci[j])
            else:
                break
        alls.append(";".join([cats[x] for x in inc]))
        uniq.append(fdr[i] < alpha)
    return best, alls, uniq


def test_classify_multilabel_equals_the_loop(synthetic_h5ad, tmp_path):
    import anndata as ad
    from scarlet import api, classifier
    A = ad.read_h5ad(synthetic_h5ad)
    lab = "cell_type" if "cell_type" in A.obs else A.obs.columns[0]
    m = api.learn(A, lab, None, blocklist="none")
    U = m.universe.copy()
    out = classifier.classify(A, U, alpha=0.2, no_save_probs=True)
    ll, cats = classifier.score_cells(A, U)[0], [c for c in U.columns if c != "bg"]
    ll_bg = classifier.score_cells(A, U)[1]
    b, al, u = _old_multilabel(ll, ll_bg, cats, out.obs["scarlet_fdr"].to_numpy(), 0.2)
    assert list(out.obs["scarlet_best"]) == b and list(out.obs["scarlet_all"]) == al
    assert list(out.obs["is_scarlet_unique"]) == [bool(x) for x in u]


def _old_mode_rows(M):
    out_v = np.full(M.shape[0], "", dtype=object)
    out_n = np.zeros(M.shape[0], dtype=np.int64)
    for i, row in enumerate(M):
        r = row[row != ""]
        if len(r):
            vals, cnt = np.unique(r, return_counts=True)
            j = int(np.argmax(cnt))
            out_v[i], out_n[i] = vals[j], int(cnt[j])
    return out_v, out_n


def test_mode_rows_and_consensus_reasons_equal_the_loops():
    from scarlet import transfer_steps as T
    rng = np.random.default_rng(3)
    ids = np.array(["B", "a", "T cell", "", T.UNASSIGNED], dtype=object)
    M = ids[rng.integers(0, 4, size=(500, 10))]
    v0, n0 = _old_mode_rows(M)
    v1, n1 = T._mode_rows(M)
    assert list(v0) == list(v1) and n0.tobytes() == n1.tobytes()
    calls = np.where(rng.random((500, 10)) < 0.85, T.UNASSIGNED, "x").astype(object)
    reasons = np.array([T.REASON_AMBIGUITY, T.REASON_BACKGROUND, T.REASON_FLIP, "zz other", ""],
                       dtype=object)[rng.integers(0, 5, size=(500, 10))]
    reasons[calls != T.UNASSIGNED] = ""
    lab, k, why = T.consensus(calls, reasons)
    for i in np.flatnonzero((lab == T.UNASSIGNED) & (k >= T.needed(10))):
        r = reasons[i][calls[i] == T.UNASSIGNED]
        r = r[r != ""]
        if len(r):
            vc = pd.Series(r).value_counts()
            top = vc[vc == vc.max()].index
            exp = next(x for x in T.REASON_ORDER if x in top) if any(x in top for x in T.REASON_ORDER) \
                else sorted(top)[0]
        else:
            exp = T.REASON_BACKGROUND
        assert why[i] == exp


def test_table_payload_rebuilds_the_rows():
    from scarlet.html_report.certainty import table_payload, table_rows
    rng = np.random.default_rng(4)
    n = 300
    df = pd.DataFrame({"cell": [f"c{i}" for i in range(n)],
                       "label": rng.choice(["T cell", "B cell", "unassigned"], n),
                       "maybe": rng.choice(["x", None], n),
                       "agree": rng.integers(0, 11, n), "fdr": np.where(rng.random(n) < 0.1, np.nan, rng.random(n)),
                       "flag": rng.random(n) < 0.3, "few": ["only"] * n})
    old = df.astype(object).where(pd.notna(df), None).to_numpy().tolist()
    T = table_payload(df, "t.tsv")
    assert table_rows(T) == old
    assert [type(v) for v in table_rows(T)[0]] == [type(v) for v in old[0]]
    assert any("codes" in c for c in T["columns"]) and any("values" in c for c in T["columns"])
