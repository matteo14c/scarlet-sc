"""technology auto (the default): learn picks its preset from the reference with the documented rule (dense when
the reference detects at least 5,000 genes per cell on average, else sparse); the choice and its reason are printed
and recorded; an explicit technology is used as given; other values are refused."""
import inspect
import zipfile

import anndata as ad
import numpy as np
import pandas as pd
import pytest
from scipy import sparse

from scarlet import api, cli
from scarlet.settings import HIGH_DENS_MIN_GENES
from synthetic_dense import _data

KW = dict(n_iterations=6, n_cells=20)


def _members(p):
    with zipfile.ZipFile(p) as z:
        return {n: z.read(n) for n in z.namelist()
                if n.endswith(".tsv") and any(k in n for k in ("marker", "programme", "universe", "background"))}


def _depth(A):
    return float(np.mean(np.diff(sparse.csr_matrix(A.X).indptr)))


def _rows(path, k, n_genes=5200, n=60):
    """n cells with exactly k detected genes each (3 labels, each with its own block of 200 genes)."""
    rng = np.random.default_rng(3)
    labels = np.repeat(["A", "B", "C"], n // 3)
    X = np.zeros((n, n_genes), dtype=np.float32)
    for i, lab in enumerate(labels):
        own = np.arange(200) + 200 * "ABC".index(lab)
        rest = rng.choice(np.setdiff1d(np.arange(n_genes), own), size=k - 200, replace=False)
        X[i, own] = rng.poisson(5, 200) + 1
        X[i, rest] = rng.poisson(1, k - 200) + 1
    A = ad.AnnData(X=sparse.csr_matrix(X), obs=pd.DataFrame({"label": labels}, index=[f"c{i}" for i in range(n)]),
                   var=pd.DataFrame(index=[f"G{j}" for j in range(n_genes)]))
    A.write_h5ad(path)
    return A


@pytest.fixture(scope="module")
def dense(tmp_path_factory):
    d = tmp_path_factory.mktemp("r062auto")
    A = _data(d / "dense.h5ad")
    assert _depth(A) >= HIGH_DENS_MIN_GENES
    S = A.copy()
    S.X = sparse.csr_matrix(np.minimum(np.asarray(A.X.todense()), 1) * (np.random.default_rng(2).random(A.shape) < .1))
    S.write_h5ad(d / "sparse.h5ad")
    assert _depth(S) < HIGH_DENS_MIN_GENES
    return d


def test_auto_resolves_high_dens_and_equals_the_explicit_run(dense, capsys):
    m = api.learn(str(dense / "dense.h5ad"), "label", str(dense / "auto.scarlet"), **KW)
    out = capsys.readouterr().out
    api.learn(str(dense / "dense.h5ad"), "label", str(dense / "hd.scarlet"), technology="dense", **KW)
    assert m.settings["resolved"]["technology"]["name"] == "dense"
    rec = m.settings["parameters"]["technology_resolved"]
    assert rec["technology"] == "dense" and rec["how"].startswith("resolved from the reference (")
    assert "resolved     technology dense (resolved from the reference" in out
    assert _members(dense / "auto.scarlet") == _members(dense / "hd.scarlet")


def test_auto_resolves_low_dens_and_equals_the_explicit_run(dense):
    m = api.learn(str(dense / "sparse.h5ad"), "label", str(dense / "auto_s.scarlet"), **KW)
    api.learn(str(dense / "sparse.h5ad"), "label", str(dense / "ld_s.scarlet"), technology="sparse", **KW)
    assert m.settings["resolved"]["technology"]["name"] == "sparse"
    assert m.settings["parameters"]["technology_resolved"]["technology"] == "sparse"
    assert _members(dense / "auto_s.scarlet") == _members(dense / "ld_s.scarlet")


@pytest.mark.parametrize("k,want", [(HIGH_DENS_MIN_GENES - 1, "sparse"), (HIGH_DENS_MIN_GENES, "dense")])
def test_the_boundary(tmp_path, k, want):
    A = _rows(tmp_path / "x.h5ad", k)
    assert _depth(A) == k
    m = api.learn(str(tmp_path / "x.h5ad"), "label", None, n_iterations=2, n_cells=10)
    assert m.settings["resolved"]["technology"]["name"] == want


def test_explicit_technology(dense):
    m = api.learn(str(dense / "dense.h5ad"), "label", None, technology="sparse", **KW)
    assert m.settings["resolved"]["technology"]["name"] == "sparse"
    assert "technology_resolved" not in m.settings["parameters"]   # as given: recorded as before
    for bad in ("other", None):
        with pytest.raises(ValueError, match="technology"):
            api.learn(str(dense / "dense.h5ad"), "label", None, technology=bad, **KW)


def test_cli_default_is_auto():
    ns = cli.build_parser().parse_args(["learn", "--reference", "x.h5ad", "--labels", "g", "--model", "m"])
    assert ns.technology == "auto"
    for name in ("learn", "transfer"):
        assert inspect.signature(getattr(api, name)).parameters["technology"].default == "auto"
