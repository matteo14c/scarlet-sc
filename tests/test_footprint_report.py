"""Memory footprint (the report): the report's UMAP input built in float32 (P2) and the self-assessment's folds taken
from the detection pattern (P4) give the same outputs, byte for byte, as the code they replace (swapped back in)."""
from __future__ import annotations

import json

import numpy as np
import pytest
from scipy import sparse

pytest.importorskip("scanpy")


def _counts_adata(seed=0, n=240, g=300):
    import anndata as ad
    import pandas as pd
    rng = np.random.default_rng(seed)
    lab = np.array(["a", "b", "c", "d"])[np.arange(n) % 4]
    lam = rng.gamma(0.4, 2.0, size=(4, g)) * 3
    X = rng.poisson(lam[np.arange(n) % 4]).astype(np.int64)
    X[:, :5] *= 300   # read-count sized values
    A = ad.AnnData(sparse.csr_matrix(X), obs=pd.DataFrame({"cell_type": lab}, index=[f"c{i}" for i in range(n)]),
                   var=pd.DataFrame(index=[f"G{i}" for i in range(g)]), uns={"note": "x"})
    A.obsm["emb"] = rng.random((n, 2))
    return A


@pytest.mark.parametrize("max_cells", [10 ** 6, 150])
@pytest.mark.parametrize("dtype", [np.int64, np.int32, np.float32])
def test_umap_input_same_coordinates(monkeypatch, max_cells, dtype):
    from scarlet.html_report import learn_report as L
    A = _counts_adata()
    A.X = A.X.astype(dtype)
    X0, d0 = A.X, A.X.data.copy()
    new = L.umap_coords(A, max_cells, 0, log=lambda *a: None)
    assert A.X is X0 and np.array_equal(A.X.data, d0) and A.X.dtype == dtype   # A unchanged
    monkeypatch.setattr(L, "_umap_input", lambda A, idx: A[idx].copy())          # the replaced code
    old = L.umap_coords(A, max_cells, 0, log=lambda *a: None)
    assert np.array_equal(new[0], old[0]) and new[1].dtype == old[1].dtype and new[1].tobytes() == old[1].tobytes()
    assert new[2] == old[2]


def test_fold_data_same_self_assessment(tmp_path, monkeypatch):
    from scarlet import api
    from scarlet.html_report import selfassess as S
    A = _counts_adata(seed=3)
    m = tmp_path / "m.scarlet"
    api.learn(A, "cell_type", str(m), blocklist=None)
    before = A.X.copy()
    C1, r1 = S.run(str(m), A, n_folds=3, seed=0, log=lambda *a: None)
    assert (A.X != before).nnz == 0 and A.X.dtype == np.int64
    monkeypatch.setattr(S, "_fold_data", lambda A: A)   # the replaced code: folds copied with their values
    C0, r0 = S.run(str(m), A, n_folds=3, seed=0, log=lambda *a: None)
    assert C1.to_csv() == C0.to_csv()
    for r in (r0, r1):   # the folds' learn and assign seconds: run times
        assert r["times"]
        for t in r["times"]:
            t.pop("learn_s")
            t.pop("assign_s")
    assert json.dumps(r1, sort_keys=True, default=str) == json.dumps(r0, sort_keys=True, default=str)
