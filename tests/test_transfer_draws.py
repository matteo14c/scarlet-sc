"""transfer in draws: the stable anchors, the refinement's pooled background, the consensus of the draws (a final
call needs its outcome in at least 80 % of the draws, rounded up: 8 of 10), quick mode (one draw), the seed list, and
reproducibility (the same seeds give byte-identical outputs, whatever ncpu)."""
from __future__ import annotations

import json

import numpy as np
import pandas as pd
import pytest
from scipy import sparse

from scarlet import transfer_steps as T

U = "unassigned"


# ----------------------------------------------------------------------------------------------- the consensus rule
@pytest.mark.parametrize("n, needed", [(1, 1), (2, 2), (3, 3), (5, 4), (7, 6), (10, 8), (20, 16), (50, 40)])
def test_consensus_needs_80_percent(n, needed):
    """The winner needs `needed` of n draws: 80 % of n, rounded up."""
    assert T.needed(n) == needed
    rows = [["a"] * k + ["b"] * (n - k) for k in range(n // 2 + 1, n + 1)]
    lab, agree, why = T.consensus(np.array(rows, dtype=object))
    for (k, l_, w) in zip(range(n // 2 + 1, n + 1), lab, why):
        if k >= needed:
            assert l_ == "a" and w == "", (n, k)
        else:
            assert l_ == U and w == T.UNDETERMINED, (n, k)
    assert list(agree) == list(range(n // 2 + 1, n + 1))


def test_consensus_ties_runner_up_and_refusals():
    calls = np.array([
        ["a"] * 5 + ["b"] * 5,                   # a tie for first: undetermined
        ["a"] * 8 + ["b"] * 1 + [U] * 1,         # 8 of 10: kept
        [U] * 9 + ["a"] * 1,                     # refused by consensus
        [U] * 8 + ["a"] * 2,                     # refused, reasons tied: background first
        ["a"] * 7 + ["b"] * 3,                   # 7 of 10: undetermined
    ], dtype=object)
    reasons = np.full(calls.shape, "", dtype=object)
    reasons[2, :9] = [T.REASON_AMBIGUITY] * 6 + [T.REASON_BACKGROUND] * 3
    reasons[3, :8] = [T.REASON_AMBIGUITY] * 4 + [T.REASON_BACKGROUND] * 4
    lab, agree, why = T.consensus(calls, reasons)
    assert list(lab) == [U, "a", U, U, U]
    assert list(why) == [T.UNDETERMINED, "", T.REASON_AMBIGUITY, T.REASON_BACKGROUND, T.UNDETERMINED]
    assert list(agree) == [5, 8, 9, 8, 7]


def test_one_draw_is_its_own_call():
    calls = np.array([["a"], [U], [U]], dtype=object)
    reasons = np.array([[""], [T.REASON_AMBIGUITY], [""]], dtype=object)
    lab, agree, why = T.consensus(calls, reasons)
    assert list(lab) == ["a", U, U] and list(agree) == [1, 1, 1]
    assert list(why) == ["", T.REASON_AMBIGUITY, T.REASON_BACKGROUND]


def test_draw_reasons():
    rb, ru = T.draw_reasons(["a", U, "b", U], ["a", U, U, U], [False, False, False, True])
    assert list(rb) == ["", T.REASON_BACKGROUND, "", T.REASON_FLIP]
    assert list(ru) == ["", T.REASON_BACKGROUND, T.REASON_AMBIGUITY, T.REASON_FLIP]


# ------------------------------------------------------------------------------------- stable anchors and the plan
def test_stable_anchors_need_every_draw():
    uq = np.array([["a", "a", "a"], ["a", "a", U], ["b", "b", "b"], [U, U, U], ["a", "b", "a"]], dtype=object)
    assert list(T.stable_anchors(uq)) == ["a", U, "b", U, U]


def test_refine_plan_and_keep_first():
    anchor = np.array(["a"] * 3 + ["b"] * 3 + ["c"] * 2 + [U] * 2 + ["a"] * 3 + [U] * 2, dtype=object)
    origin = np.array(["o1"] * 10 + ["o2"] * 5)
    fb = np.array([["a"]] * 3 + [["b"]] * 3 + [["c"]] * 2 + [["d"], [U]] + [["a"]] * 3 + [["b"], [U]], dtype=object)
    plan = T.refine_plan(anchor, origin, fb, fb)
    o1, o2 = plan
    assert o1["refined"] == ["a", "b"] and o1["not_refined"] == ["c", "d"] and o1["n_train"] == 6
    assert o1["stable_anchors"] == {"a": 3, "b": 3, "c": 2} and "skipped" not in o1
    assert o2["refined"] == ["a"] and "skipped" in o2                       # fewer than 2 refined cell types
    keep = T.keep_first(fb[:10, 0], o1["refined"])
    assert list(keep) == [False] * 6 + [True] * 3 + [False]                 # c, c, d keep; unassigned does not


def test_pooled_background_is_the_detection_over_all_cells(synthetic_adata, tmp_path):
    from scarlet import api
    from scarlet.scoring import jeffreys
    m = api.learn(synthetic_adata, "cell_type", None, n_iterations=5, n_cells=30)
    q = synthetic_adata.copy()
    m2, rec = T.pooled_background(m, q)
    uni = m2.universe
    X = q.X if sparse.issparse(q.X) else sparse.csr_matrix(q.X)
    col = pd.Index(q.var_names.astype(str)).get_indexer(uni.index.astype(str))
    det = np.asarray((X[:, col] > 0).sum(axis=0)).ravel() / q.n_obs
    np.testing.assert_allclose(uni["bg"].to_numpy(), det, rtol=0, atol=1e-12)
    assert rec == dict(cells=q.n_obs, mode="pooled", f_bg=0.5 / (q.n_obs + 1.0))
    assert jeffreys(m2)[1] == pytest.approx(0.5 / (q.n_obs + 1.0))
    np.testing.assert_allclose(m2.universe.drop(columns="bg").to_numpy(), m.universe.drop(columns="bg").to_numpy(),
                               rtol=1e-12, atol=0)   # the identities' columns unchanged


# ------------------------------------------------------------------------------------------- options and their checks
@pytest.mark.parametrize("kw, msg", [(dict(draws=0), "at least 1"), (dict(seeds=[1, 2, 3, 3]), "repeat"),
                                     (dict(draws=8, seeds=list(range(7))), "differ"),
                                     (dict(seeds="1,2,x"), "integers"), (dict(ncpu=0), "ncpu")])
def test_draw_options_are_checked(synthetic_h5ad, tmp_path, kw, msg):
    from scarlet import api
    with pytest.raises(ValueError, match=msg):
        api.transfer(str(synthetic_h5ad), "cell_type", str(synthetic_h5ad), str(tmp_path / "t"), **kw)


def test_draw_plan_defaults():
    from scarlet import api
    assert api._draw_plan(None, None) == (10, list(range(1, 11)))
    assert api._draw_plan(None, "5, 9,11,2") == (4, [5, 9, 11, 2])
    assert api._draw_plan(1, None) == (1, [1])


# ------------------------------------------------------------------------------------------------- whole transfers
@pytest.fixture(scope="module")
def two_origin(tmp_path_factory):
    from conftest import make_synthetic_adata
    d = tmp_path_factory.mktemp("draws")
    ref = make_synthetic_adata(seed=2)
    q = make_synthetic_adata(seed=3)
    q.obs = pd.DataFrame({"batch": np.where(np.arange(q.n_obs) % 2, "b1", "b2")}, index=q.obs_names)
    ref.write_h5ad(d / "ref.h5ad")
    q.write_h5ad(d / "q.h5ad")
    return d


def _run(d, name, **kw):
    from scarlet import api
    args = dict(origin_col="batch", n_iterations=10, n_cells=30, draws=7)
    args.update(kw)
    return api.transfer(str(d / "ref.h5ad"), "cell_type", str(d / "q.h5ad"), str(d / name), **args)


def test_outputs_record_and_anchors(two_origin):
    res = _run(two_origin, "a")
    L, F = res.labels, res.first_pass
    for c in ("origin", "scarlet_best", "scarlet_unique", "scarlet_best_agreement", "scarlet_unique_agreement",
              "scarlet_best_reason", "scarlet_unique_reason", "scarlet_refine_note", "scarlet_stable_anchor"):
        assert c in L.columns, c
    assert "scarlet_stable_anchor" in F.columns and "scarlet_refine_note" not in F.columns
    assert "scarlet_best_agreement" in F.columns and "scarlet_runner_up_draws" in L.columns
    for X in (L, F):
        for o in ("best", "unique"):
            refused = X[f"scarlet_{o}"] == U
            assert (X.loc[refused, f"scarlet_{o}_reason"] != "").all()
            assert (X.loc[~refused, f"scarlet_{o}_reason"] == "").all()
            assert X[f"scarlet_{o}_agreement"].between(1, 7).all()
    rec = res.record
    assert rec["schema"] == "scarlet.transfer.record/4"
    assert rec["draws"]["n"] == 7 and rec["draws"]["seeds"] == list(range(1, 8)) and rec["draws"]["mode"] == "accurate"
    assert rec["draws"]["consensus"]["needed"] == 6 and "anchor_draws" not in rec["draws"]
    assert len(rec["draws"]["per_draw"]) == 7 and len(rec["depth_match"]) == 7
    assert "seed" not in rec["parameters"] and rec["parameters"]["seeds"] == list(range(1, 8))
    # every stable anchor keeps its anchor label in every draw: unanimous, so the consensus is the anchor
    sel = L.scarlet_stable_anchor.astype(bool)
    assert sel.any() and (F.scarlet_stable_anchor.astype(bool) == sel).all()
    a = L.scarlet_unique.where(sel, "")
    assert (L.scarlet_best[sel] == a[sel]).all() and (F.scarlet_unique[sel] == a[sel]).all()
    assert (L.scarlet_unique_agreement[sel] == 7).all() and (F.scarlet_unique_agreement[sel] == 7).all()
    for g in rec["refine"]["groups"]:
        om = L.origin.astype(str) == g["origin"]
        for k, v in g["stable_anchors"].items():
            assert int(((a == k) & om).sum()) == v
    assert all(b["mode"] == "pooled" and b["cells"] == int((L.origin == b["origin"]).sum())
               for b in rec["refine"]["backgrounds"])
    s = (two_origin / "a" / "summary.txt").read_text()
    assert "at least 6 of the 7" in s and "stable anchors" in s and "undetermined across draws" in s


def test_same_seeds_same_bytes_and_ncpu_free(two_origin):
    """The draws run side by side with ncpu > 1: the outputs are the same bytes."""
    _run(two_origin, "r1")
    _run(two_origin, "r2", ncpu=3)
    for f in ("labels.tsv", "first_pass.tsv", "qc.tsv"):
        assert (two_origin / "r1" / f).read_bytes() == (two_origin / "r2" / f).read_bytes(), f
    r1 = json.loads((two_origin / "r1" / "transfer_record.json").read_text())
    r2 = json.loads((two_origin / "r2" / "transfer_record.json").read_text())
    for k in ("draws", "refine", "depth_match", "qc"):
        assert r1[k] == r2[k], k


def test_other_seeds_and_the_draws_file(two_origin):
    from scarlet import api
    kw = dict(origin_col="batch", n_iterations=10, n_cells=30, refine=True, depth_match=True, qc=False,
              technology="auto", blocklist="auto", blocklist_modules="auto", imbalanced="auto", gene_col="index",
              n_markers=None, marker_shortlist=None, draws=None, seeds=[3, 1, 4, 5, 9, 2, 6],
              alpha=0.01, flip="flag", chunk_size=10000, ncpu=1, report=False, not_default={})
    res = api._transfer(str(two_origin / "ref.h5ad"), "cell_type", str(two_origin / "q.h5ad"),
                        str(two_origin / "d"), keep_draws=True, **kw)
    D = pd.read_csv(two_origin / "d" / "draws.tsv.gz", sep="\t", keep_default_na=False, dtype={"cell": str})
    assert sorted(D.seed.unique()) == [1, 2, 3, 4, 5, 6, 9]
    for out in ("best", "unique"):
        W = D.pivot(index="cell", columns="draw", values=out).loc[res.labels.index]
        lab, k, _ = T.consensus(W.to_numpy(dtype=object))
        assert list(lab) == list(res.labels[f"scarlet_{out}"]) and list(k) == list(res.labels[f"scarlet_{out}_agreement"])
        W = D.pivot(index="cell", columns="draw", values=f"first_{out}").loc[res.first_pass.index]
        assert list(T.consensus(W.to_numpy(dtype=object))[0]) == list(res.first_pass[f"scarlet_{out}"])


def test_given_anchors_and_no_refinement(two_origin):
    from scarlet import api
    base = dict(origin_col="batch", n_iterations=10, n_cells=30, refine=True, depth_match=True, qc=False,
                technology="auto", blocklist="auto", blocklist_modules="auto", imbalanced="auto", gene_col="index",
                n_markers=None, marker_shortlist=None, draws=7, seeds=None, alpha=0.01, flip="flag",
                chunk_size=10000, ncpu=1, report=False, not_default={})
    import anndata as ad
    q = ad.read_h5ad(two_origin / "q.h5ad")
    labs = sorted(ad.read_h5ad(two_origin / "ref.h5ad").obs["cell_type"].astype(str).unique())
    b1 = list(q.obs_names[q.obs.batch.to_numpy() == "b1"].astype(str))
    given = {c: labs[0] for c in b1[:6]} | {c: labs[1] for c in b1[6:12]}
    res = api._transfer(str(two_origin / "ref.h5ad"), "cell_type", str(two_origin / "q.h5ad"),
                        str(two_origin / "g"), anchors=given, **base)
    assert res.record["internal"]["anchors_given"] is True
    L = res.labels
    assert L.loc[list(given), "scarlet_stable_anchor"].astype(bool).all()
    assert int(L.scarlet_stable_anchor.astype(bool).sum()) == len(given)
    assert (L.loc[list(given), "scarlet_best"] == pd.Series(given)).all()
    rn = api.transfer(str(two_origin / "ref.h5ad"), "cell_type", str(two_origin / "q.h5ad"), str(two_origin / "n"),
                      origin_col="batch", n_iterations=10, n_cells=30, draws=7, refine=False, qc=False)
    assert "scarlet_stable_anchor" not in rn.labels.columns and rn.record["refine"] is None
    assert (rn.labels.scarlet_best == rn.first_pass.scarlet_best).all()
    assert "a first pass only" in (two_origin / "n" / "summary.txt").read_text()


def test_quick_mode(two_origin):
    """draws=1: the one draw's calls, its own stable anchors (its unique calls), no agreement columns."""
    from scarlet import api
    kw = dict(origin_col="batch", n_iterations=10, n_cells=30, refine=True, depth_match=True, qc=False,
              technology="auto", blocklist="auto", blocklist_modules="auto", imbalanced="auto", gene_col="index",
              n_markers=None, marker_shortlist=None, draws=None, seeds=[4], alpha=0.01, flip="flag",
              chunk_size=10000, ncpu=2, report=False, not_default={})
    res = api._transfer(str(two_origin / "ref.h5ad"), "cell_type", str(two_origin / "q.h5ad"),
                        str(two_origin / "quick"), keep_draws=True, **kw)
    L, F, rec = res.labels, res.first_pass, res.record
    assert rec["draws"]["mode"] == "quick" and rec["draws"]["consensus"] is None and rec["draws"]["n"] == 1
    for X in (L, F):
        assert not [c for c in X.columns if c.endswith("_agreement") or c.endswith("_draws")]
        assert X.scarlet_best_reason.isin(["", T.REASON_BACKGROUND]).all()
    D = pd.read_csv(two_origin / "quick" / "draws.tsv.gz", sep="\t", keep_default_na=False, dtype={"cell": str})
    D = D.set_index("cell").loc[L.index]
    assert (D.best == L.scarlet_best).all() and (D.unique == L.scarlet_unique).all()
    assert (D.first_unique == F.scarlet_unique).all()
    anchors = (D.first_unique != U).to_numpy()
    plan_refined = any(g.get("refined") and "skipped" not in g for g in rec["refine"]["groups"])
    assert (L.scarlet_stable_anchor.astype(bool).to_numpy() == anchors).all() and plan_refined
    assert "Quick mode" in (two_origin / "quick" / "summary.txt").read_text()


# ------------------------------------------------------------------------------------------- documentation coverage
def test_documentation_covers_the_draws():
    """Every option, column, note, reason and record field of the draws is named in the documentation."""
    from pathlib import Path
    root = Path(__file__).resolve().parents[1]
    guide = (root / "docs" / "guide" / "transfer.md").read_text()
    output = (root / "docs" / "output.md").read_text()
    adv = (root / "docs" / "advanced.md").read_text()
    gloss = (root / "docs" / "glossary.md").read_text()
    report = "\n".join(f.read_text() for f in sorted((root / "docs" / "reports").glob("*.md")))   # the report pages
    readme = (root / "README.md").read_text()
    for opt in ("--draws", "--seeds"):
        assert opt in guide and opt in adv, opt
    assert "--anchor-draws" not in guide + adv + output
    assert "quick mode" in guide.lower() and "quick mode" in adv.lower() and "8 of the 10" in guide
    for col in ("scarlet_best_agreement", "scarlet_unique_agreement", "scarlet_best_reason", "scarlet_unique_reason",
                "scarlet_flip_draws", "scarlet_pe_low_draws", "scarlet_refine_note", "scarlet_stable_anchor"):
        assert col in output, col
    for text in (T.NOT_STABLE_NOTE, T.ORIGIN_NOT_REFINED_NOTE, T.UNDETERMINED, T.REASON_BACKGROUND,
                 T.REASON_AMBIGUITY, T.REASON_FLIP):
        assert text in output, text
    for term in ("### Draw", "### Consensus", "### Agreement", "### Stable anchors", "### Undetermined across draws",
                 "### Quick mode"):
        assert term in gloss, term
    assert "7. Agreement across the draws" in report and "`per_draw`" in report
    assert "origin not refined" in guide.lower() or "origin with fewer than 2 refined" in guide
    limits = (root / "docs" / "limitations.md").read_text()   # the known limitations (linked from the README)
    assert "Stable anchors keep their label" in limits and "Zheng 68K" in limits and "docs/limitations.md" in readme
    assert not (root / "API.md").exists() and not (root / "DEFAULTS.md").exists()
