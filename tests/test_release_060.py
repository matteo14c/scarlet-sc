"""The scoring defaults: learn records each identity's effective number of cells and the positive-evidence reference;
assign scores with the Jeffreys floor, applies the flip test (refused with a sparse model, a flag with a dense one) and
flags the positive evidence; transfer matches depths, resolves the technology from the shallower side, refines where
anchored and records the flip as a flag."""
from __future__ import annotations

import json

import numpy as np
import pandas as pd
import pytest
from scipy import sparse

from conftest import make_synthetic_adata

U = "unassigned"


def _split(a, seed=0):
    rng = np.random.default_rng(seed)
    idx = rng.permutation(a.n_obs)
    return a[np.sort(idx[: a.n_obs // 2])].copy(), a[np.sort(idx[a.n_obs // 2:])].copy()


@pytest.fixture(scope="module")
def models(tmp_path_factory):
    from scarlet import api
    tmp = tmp_path_factory.mktemp("r060")
    a = make_synthetic_adata(seed=3)
    ref, qry = _split(a)
    rp, qp = tmp / "ref.h5ad", tmp / "query.h5ad"
    ref.write_h5ad(rp)
    qry.write_h5ad(qp)
    m_def = api.learn(str(rp), "cell_type", str(tmp / "default.scarlet"), blocklist=None)
    return dict(tmp=tmp, ref=rp, query=qp, default=m_def)


def test_learn_records_n_eff(models):
    from scarlet.scoring import jeffreys
    for m in (models["default"],):
        st = m.marker_finder_params["scoring_stats"]
        for k, v in st["identities"].items():
            assert 1 <= v["n_eff"] <= v["n_k"]
        fk, fbg = jeffreys(m)
        assert all(abs(fk[k] - 0.5 / (st["identities"][k]["n_eff"] + 1)) < 1e-15 for k in fk)
    assert "positive_evidence.json" in models["default"].files


def test_default_assign_j_flip_refuse_and_pe(models, tmp_path):
    import anndata as ad
    from scarlet import api
    from scarlet.scoring import FLIP_REASON
    q = str(models["query"])
    out = tmp_path / "ann.h5ad"
    lab = api.assign(models["default"], q, str(out))
    A = ad.read_h5ad(out)
    assert "f_k" in A.uns["scarlet_floor"] and A.uns["scarlet_assign"]["flip"] == "refuse"  # sparse
    fl = A.obs.scarlet_flip.astype(bool).to_numpy()
    assert (A.obs.scarlet_top1_eps.astype(str).to_numpy() != A.obs.scarlet_top1_j.astype(str).to_numpy()).tolist() \
        == fl.tolist()
    assert (lab.scarlet_best.to_numpy()[fl] == U).all() and set(A.obs.scarlet_reason[fl]) <= {FLIP_REASON}
    # flag: the same scoring, the calls kept
    fg = api.assign(models["default"], q, None, flip="flag")
    nf = api._assign(models["default"], q, None, flip="off")
    assert (fg.scarlet_best == nf.scarlet_best).all() and (fg.scarlet_flip.astype(bool).to_numpy() == fl).all()
    assert (lab.scarlet_best.to_numpy()[~fl] == nf.scarlet_best.to_numpy()[~fl]).all()
    # positive evidence (the ratio rule): label only, flagged <=> detected < floor(E_k / 5)
    for c in ("scarlet_pe_detected", "scarlet_pe_expected", "scarlet_pe_required", "scarlet_pe_typical",
              "scarlet_pe_low"):
        assert c in A.obs.columns
    o = A.obs
    assert (o.scarlet_pe_required == np.floor(o.scarlet_pe_expected / 5)).all()
    assert (o.scarlet_pe_low.astype(bool) == (o.scarlet_pe_detected < o.scarlet_pe_required)).all()
    rec = json.loads(models["default"].files["positive_evidence.json"])
    assert all(v["typical_median"] is not None and v["required"] == int(np.floor(v["expected"] / 5))
               for v in rec["identities"].values())


def test_j_equals_floor_arithmetic(models):
    """j = the rates clamped to [f_k, 1 - f_k] then the 1e-5 guard (scarlet.scoring.clamp), by hand."""
    import anndata as ad
    from scarlet.classifier import score_cells
    from scarlet.scoring import GUARD, jeffreys
    m = models["default"]
    A = ad.read_h5ad(models["query"])
    U_ = m.universe
    fk, fbg = jeffreys(m)
    ll, llbg, cats, _ = score_cells(A, U_.copy(), floor=(fk, fbg))
    genes = [g for g in A.var_names if g in U_.index]
    X = (A[:, genes].X > 0).toarray().astype(float)
    for j, k in enumerate(cats):
        p = np.clip(np.minimum(np.maximum(U_.loc[genes, k].to_numpy(), fk[k]), 1 - fk[k]), GUARD, 1 - GUARD)
        assert np.allclose(ll[:, j], X @ np.log(p) + (1 - X) @ np.log(1 - p))


def test_flip_modes():
    """assign refuses flipped calls by default on every technology; flag and off (the refinement) are the others;
    auto is gone."""
    import inspect
    from scarlet import api
    from scarlet.scoring import resolve_flip
    assert inspect.signature(api.assign).parameters["flip"].default == "refuse"
    assert inspect.signature(api.transfer).parameters["flip"].default == "flag"
    assert [resolve_flip(x) for x in ("refuse", "flag", "off")] == ["refuse", "flag", "off"]
    for bad in ("auto", "yes"):
        with pytest.raises(ValueError):
            resolve_flip(bad)


def _two_origin_query(a, seed=1):
    rng = np.random.default_rng(seed)
    X = a.X.toarray()
    o = np.where(np.arange(a.n_obs) % 2 == 0, "shallow", "deep").astype(object)
    sh = o == "shallow"
    X[sh] = X[sh] * (rng.random(X[sh].shape) < 0.5)
    a = a.copy()
    a.X = sparse.csr_matrix(X.astype(np.float32))
    a.obs["origin"] = o
    return a


def test_transfer_defaults(models, tmp_path):
    import anndata as ad
    from scarlet import api
    q = _two_origin_query(ad.read_h5ad(models["query"]))
    qp = tmp_path / "q2.h5ad"
    q.write_h5ad(qp)
    r = api.transfer(str(models["ref"]), "cell_type", str(qp), str(tmp_path / "def"), origin_col="origin",
                     blocklist=None, qc=False)
    rec = json.loads((tmp_path / "def" / "transfer_record.json").read_text())
    assert rec["parameters"]["depth_match"] == "capture" and rec["parameters"]["refine"] == "anchored"
    assert rec["parameters"]["flip"] == "flag" and rec["parameters"]["technology"] == "sparse"
    assert "shallower side" in rec["preflight"]["technology_resolved"]["how"]
    assert {d["origin"] for d in rec["depth_match"][0]["origins"]} == {"deep", "shallow"}
    assert rec["refine"]["anchors"] == "stable" and "anchor_draws" not in rec["refine"]
    assert rec["draws"]["n"] == 10 and rec["draws"]["seeds"] == list(range(1, 11)) and rec["draws"]["mode"] == "accurate"
    L = pd.read_csv(tmp_path / "def" / "labels.tsv", sep="\t", index_col=0)
    for c in ("scarlet_flip", "scarlet_refine_note", "scarlet_pe_low", "scarlet_stable_anchor",
              "scarlet_best_agreement", "scarlet_unique_reason"):
        assert c in L.columns
    F = pd.read_csv(tmp_path / "def" / "first_pass.tsv", sep="\t", index_col=0)
    assert not F.scarlet_flip.isna().any()
    rn = api.transfer(str(models["ref"]), "cell_type", str(qp), str(tmp_path / "nor"), origin_col="origin",
                      blocklist=None, qc=False, refine=False, depth_match=False, draws=7)
    recn = json.loads((tmp_path / "nor" / "transfer_record.json").read_text())
    assert recn["depth_match"] == [] and recn["parameters"]["refine"] == "no"
    # refuse: the first pass's flipped cells stay refused
    rr = api.transfer(str(models["ref"]), "cell_type", str(qp), str(tmp_path / "ref"), origin_col="origin",
                      blocklist=None, qc=False, flip="refuse", draws=7)
    Lr = pd.read_csv(tmp_path / "ref" / "labels.tsv", sep="\t", index_col=0)
    fl = Lr.scarlet_flip.astype(str).str.lower().eq("true").to_numpy()
    assert (Lr.scarlet_best.to_numpy()[fl] == U).all()


def test_transfer_dense_reads_the_values(models, tmp_path, monkeypatch):
    """With the dense preset every learn of transfer reads detection as count > 0 of the cells' input values (no
    binarisation by a threshold): a learn receives the values or their detection pattern, which is exactly the input's
    count > 0 for the refinement (the query's own cells) and for an unthinned reference, and exactly the thinned
    matrix's detection for a depth-matched reference, whose thinning only removes detections of the input."""
    import anndata as ad
    from scarlet import api, depth_match
    seen, thinned = [], []
    orig, orig_thin = api._learn_into, depth_match.thin_draw

    def spy(model_path, data_path, labels, kw, resolved, command, params, inputs, **k):
        A = k.get("data")
        A = ad.read_h5ad(data_path) if A is None else A
        seen.append((str(params.get("step", "")), A.obs_names.astype(str).to_numpy(), sparse.csr_matrix(A.X) != 0))
        return orig(model_path, data_path, labels, kw, resolved, command, params, inputs, **k)

    def spy_thin(prep, seed):
        Y, r = orig_thin(prep, seed)
        if Y is not None:
            thinned.append((r["thinned"], sparse.csr_matrix(prep["X"]) != 0, sparse.csr_matrix(Y) != 0))
        return Y, r
    monkeypatch.setattr(api, "_learn_into", spy)
    monkeypatch.setattr(depth_match, "thin_draw", spy_thin)
    api.transfer(str(models["ref"]), "cell_type", str(models["query"]), str(tmp_path / "t"), blocklist=None,
                 qc=False, technology="dense", n_markers=10, marker_shortlist=30, draws=7)
    R, Q = ad.read_h5ad(models["ref"]), ad.read_h5ad(models["query"])
    R.obs_names, Q.obs_names = R.obs_names.astype(str), Q.obs_names.astype(str)
    steps = {s.split(" ")[0] for s, _, _ in seen}
    assert {"reference", "refine"} <= steps, steps
    for side, before, after in thinned:
        assert (after > before).nnz == 0, side   # thinning removes detections, never adds one
    exact = thin = 0
    for step, cells, got in seen:
        src = Q if step.startswith("refine") else R
        det = sparse.csr_matrix(src[cells].X) > 0
        if (det != got).nnz == 0:
            exact += 1
            continue
        assert not step.startswith("refine"), step
        assert (got > det).nnz == 0, step
        assert any(s == "reference" and Y.shape == got.shape and (Y != got).nnz == 0 and (b != det).nnz == 0
                   for s, b, Y in thinned), step   # the draw's thinned reference, itself thinned from count > 0
        thin += 1
    print(f"dense learns: {exact} equal to count > 0 of the input, {thin} equal to a thinned reference's detection")


def test_depth_match_units():
    from scarlet import depth_match as D
    rng = np.random.default_rng(6)
    lam = np.exp(rng.normal(0, 1.5, 300))
    deep = sparse.csr_matrix(rng.poisson(lam * 2.0, size=(400, 300)).astype(np.float32))
    shallow = sparse.csr_matrix(rng.poisson(lam * 0.3, size=(400, 300)).astype(np.float32))
    R, Q, rec = D.match_depth(deep, shallow, "capture", seed=2)
    assert rec["thinned"] == "reference" and Q is shallow
    d_after = np.asarray((R > 0).mean(axis=0)).ravel()
    s = np.asarray((shallow > 0).mean(axis=0)).ravel()
    assert np.abs(d_after - s).mean() < 0.05
    R2, _, _ = D.match_depth(deep, shallow, "capture", seed=2)
    assert (R2 != R).nnz == 0
    with pytest.raises(ValueError):
        D.match_depth(deep, shallow, "trim", seed=2)        # capture is the one thinning
    _, _, same = D.match_depth(deep, deep[::-1].copy(), "capture", seed=1)
    assert same["direction"] == "same"
