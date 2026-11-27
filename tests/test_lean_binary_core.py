"""Lean memory: the binary core (default) gives the same results and outputs as the value matrix
(SCARLET_BINARY_CORE=0): learn's model, assign's labels and annotated output (X = the input's X)."""
import os
import zipfile

import anndata as ad
import numpy as np
import pandas as pd
import pytest
from scipy import sparse


def _data(path):
    rng = np.random.default_rng(0)
    n, g = 240, 120
    labels = np.repeat(["A", "B", "C", "D"], n // 4)
    X = rng.poisson(0.3, size=(n, g)).astype(np.float32)
    for k, lab in enumerate(["A", "B", "C"]):  # D shares C's markers (a merge candidate)
        X[labels == lab, k * 10:(k + 1) * 10] += rng.poisson(4, size=((labels == lab).sum(), 10))
    X[labels == "D", 20:30] += rng.poisson(4, size=((labels == "D").sum(), 10))
    A = ad.AnnData(X=sparse.csr_matrix(X), obs=pd.DataFrame({"label": labels}, index=[f"c{i}" for i in range(n)]),
                   var=pd.DataFrame(index=[f"G{j}" for j in range(g)]))
    A.write_h5ad(path)
    return A


def _run(tmp, mode, monkeypatch):
    from scarlet import api
    monkeypatch.setenv("SCARLET_BINARY_CORE", mode)
    d = tmp / f"core{mode}"
    d.mkdir()
    data = str(tmp / "data.h5ad")
    kw = dict(n_iterations=6, n_cells=30)
    api.learn(data, "label", str(d / "m.scarlet"), **kw)
    lab = api.assign(str(d / "m.scarlet"), data, str(d / "assigned.h5ad"))
    return d, lab


def _model_members(p):
    with zipfile.ZipFile(p) as z:
        return {n: z.read(n) for n in z.namelist() if n.endswith((".tsv",))}


def test_binary_core_identical(tmp_path, monkeypatch):
    A = _data(tmp_path / "data.h5ad")
    d1, lab1 = _run(tmp_path, "1", monkeypatch)
    d0, lab0 = _run(tmp_path, "0", monkeypatch)
    assert _model_members(d1 / "m.scarlet") == _model_members(d0 / "m.scarlet")
    assert lab1.equals(lab0)
    a1, a0 = ad.read_h5ad(d1 / "assigned.h5ad"), ad.read_h5ad(d0 / "assigned.h5ad")
    assert (a1.X != A.X).nnz == 0 and (a0.X != A.X).nnz == 0 and a1.X.dtype == a0.X.dtype
    assert a1.obs.equals(a0.obs) and np.array_equal(a1.obsm["scarlet_probs"], a0.obsm["scarlet_probs"])
