"""The report: refusal reasons, the CLI, and end to end that the report changes no call."""
from __future__ import annotations

import json
import re
import zipfile

import numpy as np
import pandas as pd
import pytest

from conftest import make_synthetic_adata

U = "unassigned"


# --- refusal reasons -------------------------------------------------------------------------------
def test_assign_reasons_follow_scarlets_steps():
    from scarlet.html_report import reasons as R
    best = np.array(["a", U, U, "b", "c"])
    uniq = np.array(["a", U, U, U, "c"])
    ll1 = np.array([-1.0, -5.0, -2.0, -1.0, -1.0])
    llbg = np.array([-9.0, -4.0, -9.0, -9.0, -9.0])      # cell 1: background wins
    fdr = np.array([0.001, 1.0, 0.001, 0.2, 0.001])      # cell 3: ambiguous
    flip = np.array([False, False, True, False, False])  # cell 2: flip, refused under 'refuse'
    rb, ru = R.assign_reasons(best, uniq, ll1, llbg, fdr, 0.01, flip=flip, flip_mode="refuse")
    assert list(rb) == ["", R.BACKGROUND, R.FLIP, "", ""]
    assert list(ru) == ["", R.BACKGROUND, R.FLIP, R.AMBIGUOUS, ""]
    assert R.check_exact(best, uniq, rb, ru)


def test_unexplained_refusal_fails_the_build():
    from scarlet.html_report import reasons as R
    best, uniq = np.array([U]), np.array([U])
    rb, ru = R.assign_reasons(best, uniq, [-1.0], [-9.0], [0.5], 0.01, flip=[True], flip_mode="flag")
    with pytest.raises(AssertionError):
        R.check_exact(best, uniq, rb, ru)


def test_reason_names_are_the_outputs():
    """The report names the reasons as transfer's outputs and the documentation do."""
    from scarlet import transfer_steps as T
    from scarlet.html_report import reasons as R
    assert (R.BACKGROUND, R.AMBIGUOUS, R.FLIP, R.UNDETERMINED) == (T.REASON_BACKGROUND, T.REASON_AMBIGUITY,
                                                                 T.REASON_FLIP, T.UNDETERMINED)


def test_counts_have_no_residual_bucket():
    from scarlet.html_report import reasons as R
    best = np.array(["a", U, U, "b"])
    uniq = np.array(["a", U, U, U])
    rb, ru = R.assign_reasons(best, uniq, [-1, -5, -2, -1], [-9, -4, -9, -9], [0, 1, 0, 0.5], 0.01,
                              flip=[False, False, True, False], flip_mode="refuse")
    c = R.counts(best, uniq, rb, ru, flip=[False, False, True, False])
    assert c["best_other"] == 0 and c["unique_other"] == 0
    assert (c["best_background"], c["best_flip_refused"], c["unique_ambiguity"]) == (1, 1, 1)


# --- the CLI ---------------------------------------------------------------------------------------------------
def test_cli_report_verb_and_flags():
    from scarlet.cli import build_parser
    p = build_parser()
    a = p.parse_args(["report", "run.scarlet", "--no-self-assessment"])
    assert a.run == "run.scarlet" and a.self_assessment is False
    assert p.parse_args(["report", "run.scarlet"]).self_assessment is True
    assert p.parse_args(["learn", "--reference", "r", "--labels", "l", "--model", "m", "--report"]).report is True
    assert p.parse_args(["learn", "--reference", "r", "--labels", "l", "--model", "m"]).report is False
    assert p.parse_args(["transfer", "--reference", "r", "--labels", "l", "--query", "q", "--output", "o",
                         "--report"]).report is True


def _members(path):
    with zipfile.ZipFile(path) as z:
        return {n: z.read(n) for n in z.namelist()}


def _strict_json(html_path):
    h = open(html_path, encoding="utf-8").read()
    m = re.search(r"<script type='application/json' id='report-data'>(.*?)</script>", h, re.S)
    json.loads(m.group(1), parse_constant=lambda c: (_ for _ in ()).throw(ValueError(c)))
    # self-contained: no external script or stylesheet (plotly.js is embedded; its own default URLs for geo maps and
    # MathJax are not used by the report)
    assert not re.search(r"<script[^>]+src=", h) and not re.search(r"<link[^>]+href=", h)
    return h


@pytest.fixture(scope="module")
def ref_path(tmp_path_factory):
    pytest.importorskip("plotly")
    p = tmp_path_factory.mktemp("report") / "ref.h5ad"
    make_synthetic_adata(seed=5).write_h5ad(p)
    return p


VOLATILE = {"provenance.json", "settings.json", "marker_finder_params.json", "manifest.json", "summary.txt"}


def test_learn_report_changes_no_call(ref_path, tmp_path):
    from scarlet import api
    m0 = api.learn(str(ref_path), "cell_type", str(tmp_path / "a" / "m.scarlet"), n_iterations=20, report=False)
    (tmp_path / "b").mkdir()
    m1 = api.learn(str(ref_path), "cell_type", str(tmp_path / "b" / "m.scarlet"), n_iterations=20, report=True)
    a, b = _members(tmp_path / "a" / "m.scarlet"), _members(tmp_path / "b" / "m.scarlet")
    assert set(a) == set(b)
    assert {n for n in a if a[n] != b[n]} <= VOLATILE
    q = make_synthetic_adata(seed=6)
    q.obs = pd.DataFrame(index=q.obs_names)
    l0 = api.assign(str(tmp_path / "a" / "m.scarlet"), q.copy(), None)
    l1 = api.assign(str(tmp_path / "b" / "m.scarlet"), q.copy(), None)
    pd.testing.assert_frame_equal(l0, l1)
    out = tmp_path / "b" / "m_report"
    h = _strict_json(out / "report.html")
    text = re.sub(r"<script.*?</script>", "", h, flags=re.S)   # the page without the embedded plotly code
    assert "learn at package defaults" in text and "auto: " in text and not re.search(r"\bU4\b", text)
    man = json.load(open(out / "bundle" / "manifest.json"))
    for t in ("depth", "self_assessment", "analytics", "umap", "weak_classes", "markers", "universe"):
        assert t in man["tables"], t
    sa = pd.read_csv(out / "bundle" / "tables" / "self_assessment.tsv", sep="\t")
    assert len(sa) == 240 and sorted(sa.fold.unique()) == [1, 2, 3, 4, 5]   # numbered as the log's 'fold k of 5'
    assert m0.identities == m1.identities


def test_transfer_and_assign_reports_change_no_call(ref_path, tmp_path):
    import anndata as ad
    from scarlet import api
    q = make_synthetic_adata(seed=7)
    q.obs = pd.DataFrame(index=q.obs_names)
    qp = tmp_path / "q.h5ad"
    q.write_h5ad(qp)
    api.transfer(str(ref_path), "cell_type", str(qp), str(tmp_path / "t0"), n_iterations=20, report=False, draws=7)
    api.transfer(str(ref_path), "cell_type", str(qp), str(tmp_path / "t1"), n_iterations=20, report=True, draws=7)
    for f in ("first_pass.tsv", "labels.tsv", "qc.tsv"):
        assert open(tmp_path / "t0" / f, "rb").read() == open(tmp_path / "t1" / f, "rb").read(), f
    h = _strict_json(tmp_path / "t1" / "report" / "report.html")
    _strict_json(tmp_path / "t1" / "report" / "supplement_A" / "report.html")
    # first pass and refined, best and unique, side by side, from first_pass.tsv and labels.tsv (on the report's own
    # page on refinement, linked from section 6)
    assert 'href="refinement.html"' in h
    assert "First pass and refined, side by side" in _strict_json(tmp_path / "t1" / "report" / "refinement.html")
    SB = pd.read_csv(tmp_path / "t1" / "report" / "bundle" / "tables" / "first_pass_vs_refined.tsv", sep="\t")
    FP = pd.read_csv(tmp_path / "t1" / "first_pass.tsv", sep="\t")
    LB = pd.read_csv(tmp_path / "t1" / "labels.tsv", sep="\t")
    for e in ("best", "unique"):
        for side, X in (("first_pass", FP), ("refined", LB)):
            assert SB[f"{side}_{e}"].sum() == len(X)
            got = SB.set_index("identity")[f"{side}_{e}"]
            assert (got[got > 0] == X[f"scarlet_{e}"].astype(str).value_counts().reindex(got[got > 0].index)).all()
    api.assign(str(tmp_path / "t0" / "reference.scarlet"), str(qp), str(tmp_path / "an.h5ad"))
    before = ad.read_h5ad(tmp_path / "an.h5ad").obs.copy()
    api.report(str(tmp_path / "an.h5ad"), model=str(tmp_path / "t0" / "reference.scarlet"))
    pd.testing.assert_frame_equal(before, ad.read_h5ad(tmp_path / "an.h5ad").obs)   # the run's output unchanged
    _strict_json(tmp_path / "an_report" / "report.html")


def test_report_finds_relative_recorded_inputs(ref_path, tmp_path, monkeypatch):
    """a transfer run with relative paths, reported later from another working directory"""
    import shutil
    from scarlet import api
    run = tmp_path / "run"
    run.mkdir()
    shutil.copy(ref_path, run / "ref.h5ad")
    q = make_synthetic_adata(seed=8)
    q.obs = pd.DataFrame(index=q.obs_names)
    q.write_h5ad(run / "q.h5ad")
    monkeypatch.chdir(run)
    api.transfer("ref.h5ad", "cell_type", "q.h5ad", "./t", n_iterations=20, draws=7)
    monkeypatch.chdir(tmp_path)
    p = api.report(str(run / "t"), self_assessment=False)
    _strict_json(p)


def test_learn_report_defaults_with_resolved_technology(ref_path, tmp_path):
    """0.6.2: learn at defaults resolves and records its technology; the report still says 'package defaults'.
    An explicit technology is reported as a non-default parameter."""
    from scarlet import api
    for d in ("a", "b"):
        (tmp_path / d).mkdir()
    api.learn(str(ref_path), "cell_type", str(tmp_path / "a" / "m.scarlet"), report=True)
    p = json.load(open(tmp_path / "a" / "m_report" / "bundle" / "manifest.json"))["meta"]["provenance"]
    assert p["model"]["at_package_defaults"] is True, p["model"]["non_default_parameters"]
    api.learn(str(ref_path), "cell_type", str(tmp_path / "b" / "m.scarlet"), technology="sparse", report=True)
    p = json.load(open(tmp_path / "b" / "m_report" / "bundle" / "manifest.json"))["meta"]["provenance"]
    assert p["model"]["non_default_parameters"] == {"technology": {"value": "sparse", "default": "auto"}}


def test_the_assign_view_reads_the_recorded_alpha(ref_path, tmp_path):
    """0.6.3: assign records its alpha; the report reproduces the stored calls with it and states the option."""
    from scarlet import api
    api.learn(str(ref_path), "cell_type", str(tmp_path / "m.scarlet"), n_iterations=20)
    q = make_synthetic_adata(seed=9)
    q.obs = pd.DataFrame(index=q.obs_names)
    api.assign(str(tmp_path / "m.scarlet"), q, str(tmp_path / "an.h5ad"), alpha=0.2)
    api.report(str(tmp_path / "an.h5ad"), model=str(tmp_path / "m.scarlet"))
    p = json.load(open(tmp_path / "an_report" / "bundle" / "manifest.json"))["meta"]["provenance"]
    assert p["assign"]["alpha"] == 0.2 and p["assign"]["alpha_note"] == "recorded by the assign run"
    assert p["assign"]["non_default_parameters"] == {"alpha": {"value": 0.2, "default": 0.01}}
    h = re.sub(r"<script.*?</script>", "", open(tmp_path / "an_report" / "report.html").read(), flags=re.S)
    assert "the published benchmarks used the defaults of the advanced options" in h
