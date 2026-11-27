"""Lean memory: chunked row statistics and the copy-free classify give the same results as before."""
import numpy as np
import pandas as pd
import anndata as ad
from scipy import sparse

from scarlet import settings
from scarlet.classifier import classify


def ref_row_stats(X):  # 0.2.0's _row_stats, verbatim apart from the name
    import scipy.sparse as sp
    if sp.issparse(X):
        X = sp.csr_matrix(X)
        rows = np.repeat(np.arange(X.shape[0]), np.diff(X.indptr))
        detected = np.bincount(rows, weights=(X.data > 0), minlength=X.shape[0])
        totals = np.bincount(rows, weights=X.data.astype(np.float64), minlength=X.shape[0])
    else:
        X = np.asarray(X)
        detected = (X > 0).sum(axis=1)
        totals = X.sum(axis=1)
    return np.asarray(detected, dtype=float), np.asarray(totals, dtype=float)


def _X(seed, n=5000, g=300, d=0.2):
    X = sparse.random(n, g, density=d, format="csr", random_state=seed, dtype=np.float32)
    X.data = np.round(X.data * 50, 3).astype(np.float32)
    return X


def test_row_stats_identical():
    for seed, chunk in ((1, 7), (2, 20000), (3, 333)):
        X = _X(seed)
        X.data[::11] = 0  # stored zeros
        d1, t1 = settings._row_stats(X, chunk=chunk)
        d0, t0 = ref_row_stats(X)
        assert np.array_equal(d1, d0) and np.array_equal(t1, t0)
    Z = sparse.csr_matrix((3, 4))
    assert all(np.array_equal(a, b) for a, b in zip(settings._row_stats(Z), ref_row_stats(Z)))


def _adata(seed):
    X = _X(seed, n=400, g=60, d=0.3)
    obs = pd.DataFrame({"x": np.arange(400)}, index=[f"c{i}" for i in range(400)])
    A = ad.AnnData(X=X, obs=obs, var=pd.DataFrame(index=[f"g{j}" for j in range(60)]))
    A.uns["k"] = {"a": 1}
    A.obsm["e"] = np.ones((400, 2))
    return A


def _universe():
    rng = np.random.default_rng(0)
    U = pd.DataFrame(rng.uniform(0.01, 0.9, size=(40, 3)), index=[f"g{j}" for j in range(40)], columns=["A", "B", "bg"])
    return U


def test_classify_copy_free_equals_full_copy():
    A = _adata(5)
    before_obs, before_X = A.obs.copy(), A.X.copy()
    new = classify(A, _universe())
    old = classify(_adata(5)[:, :], _universe())  # a view takes the former full-copy path
    assert (new.X != old.X).nnz == 0
    assert new.obs.equals(old.obs) and new.var.equals(old.var)
    assert np.array_equal(new.obsm["scarlet_probs"], old.obsm["scarlet_probs"])
    assert set(new.uns) == set(old.uns)
    assert A.obs.equals(before_obs) and (A.X != before_X).nnz == 0 and "scarlet_best" not in A.obs
    assert new.X is A.X  # shared, not copied
