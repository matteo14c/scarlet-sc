"""The report's certainty, boundary, anchor and refusal views (transfer and learn), the transfer's runner-up columns,
the report with in-memory inputs, the plotly check and the palette."""
from __future__ import annotations

import builtins
import json
import re

import numpy as np
import pandas as pd
import pytest

from conftest import make_synthetic_adata
from scarlet.html_report.certainty import table_rows

U = "unassigned"


def _page_data(html_path):
    h = open(html_path, encoding="utf-8").read()
    m = re.search(r"<script type='application/json' id='report-data'>(.*?)</script>", h, re.S)
    data = json.loads(m.group(1).replace("<\\/", "</"))
    # self-contained: no external script, stylesheet, image or frame (plotly.js is embedded)
    assert not re.search(r"<(script|img|iframe)[^>]+src=", h) and not re.search(r"<link[^>]+href=", h)
    return h, data


def test_runner_up_rule():
    from scarlet.transfer_steps import runner_up
    best = np.array(["A", U, "A", "B"], dtype=object)
    t1 = np.array([["A", "A", "B"], ["C", "C", "B"], ["A", "A", "A"], ["B", "B", "B"]], dtype=object)
    t2 = np.array([["B", "C", "A"], ["B", "B", "C"], ["", "", ""], ["C", "A", "A"]], dtype=object)
    top1, ru, n = runner_up(best, t1, t2)
    assert list(top1) == ["A", "C", "A", "B"]
    # cell 0: A called; per draw the best other identity: B, C, B -> B (2); cell 1: refused, reference = top1 C:
    # B, B, B -> B (3); cell 2: a single identity -> none; cell 3: B called: C, A, A -> A (2)
    assert list(ru) == ["B", "B", "", "A"] and list(n) == [2, 3, 0, 2]
    # ties go to the alphabetically first
    _, ru2, n2 = runner_up(np.array(["A"], dtype=object), np.array([["A", "A"]], dtype=object),
                           np.array([["C", "B"]], dtype=object))
    assert ru2[0] == "B" and n2[0] == 1


@pytest.fixture(scope="module")
def runs(tmp_path_factory):
    pytest.importorskip("plotly")
    from scarlet import api
    d = tmp_path_factory.mktemp("views")
    ref = make_synthetic_adata(seed=5)
    q = make_synthetic_adata(seed=7)
    q.obs = pd.DataFrame(index=q.obs_names)
    # in-memory inputs, report=True: the report reads the objects given
    res = api.transfer(ref, "cell_type", q, str(d / "t"), n_iterations=20, draws=7, report=True)
    api.learn(ref, "cell_type", str(d / "m.scarlet"), n_iterations=20, report=True)
    return d, res


def test_transfer_writes_top1_and_runner_up(runs):
    d, res = runs
    for X in (res.first_pass, res.labels):
        for c in ("scarlet_top1", "scarlet_runner_up", "scarlet_runner_up_draws"):
            assert c in X.columns, c
        b, r = X.scarlet_best.astype(str), X.scarlet_runner_up.astype(str)
        assert ((r == "") | (r != b)).all()                       # never the cell's own best call
        assert X.scarlet_runner_up_draws.between(0, 7).all()
        assert ((X.scarlet_runner_up_draws == 0) == (r == "")).all()
        assert (X.scarlet_top1.astype(str) != U).all()


def test_transfer_report_views(runs):
    d, res = runs
    out = d / "t" / "report"
    h, data = _page_data(out / "report.html")
    _page_data(out / "supplement_A" / "report.html")
    tab = out / "bundle" / "tables"
    LB = pd.read_csv(d / "t" / "labels.tsv", sep="\t", keep_default_na=False)
    FP = pd.read_csv(d / "t" / "first_pass.tsv", sep="\t", keep_default_na=False)
    # boundary map: the counts are the outputs' cells with a best call and no unique call
    BR = pd.read_csv(tab / "boundary.tsv", sep="\t", keep_default_na=False)
    for st, X in (("first pass", FP), ("refined", LB)):
        n = int(((X.scarlet_best != U) & (X.scarlet_unique == U)).sum())
        assert int(BR[BR.stage == st].cells.sum()) == n, st
    # per identity: calls and stable anchors as in the outputs
    PI = pd.read_csv(tab / "per_identity.tsv", sep="\t", keep_default_na=False)
    assert PI.called_best.sum() == int((LB.scarlet_best != U).sum())
    assert PI.called_unique.sum() == int((LB.scarlet_unique != U).sum())
    assert PI.stable_anchors.sum() == int(LB.scarlet_stable_anchor.astype(str).str.lower().eq("true").sum())
    assert PI.ranked_first.sum() == len(LB)
    assert PI.refined.map(lambda v: v == "yes" or v.startswith("no: ") or v.startswith("origin not refined")).all()
    # the per-cell table: every cell, with agreement, reason, runner-up and anchor flag; all embedded for download
    QC = pd.read_csv(tab / "query_cells.tsv", sep="\t", keep_default_na=False)
    for c in ("best_agreement", "reason_unique", "runner_up", "anchor", "first_runner_up", "stable_anchor"):
        assert c in QC.columns, c
    assert len(QC) == len(LB) and len(table_rows(data["tables"]["cells"])) == len(LB)
    # the switchable query UMAP: four outputs, four views, every query cell (fewer than umap_cells)
    # section 2: A's labels and B's switchable UMAP side by side, each with its own highlight
    assert [q["div"] for q in data["qumaps"]] == ["umapA", "qumap"]
    S = data["qumaps"][1]
    assert 'id="umapA_hl"' in h and 'id="qumap_hl"' in h
    assert [o["name"] for o in S["outputs"]] == ["final best", "final unique", "first pass best", "first pass unique"]
    assert S["views"] == ["label", "agreement", "ambiguity", "anchors", "unlabelled"]
    assert len(S["x"]) == len(LB) and sum(S["anchor"]) == int(LB.scarlet_stable_anchor.astype(str).str.lower().eq("true").sum())
    assert all(0 < v <= 1 for o in S["outputs"] for a in o["shade"].values() for v in a)
    assert all(set(o["shade"]) == {"agreement", "ambiguity"} for o in S["outputs"])
    for k in ("Boundary cells", "Per identity", "Per-cell data", "Download all", "colour-blind-safe"):
        assert k in h, k
    prov = json.load(open(out / "bundle" / "manifest.json"))["meta"]["provenance"]
    assert prov["reference"]["file"]["path"] == "<AnnData>" and "not hashed" in prov["reference"]["file"]["note"]


def test_learn_report_views(runs):
    d, _ = runs
    out = d / "m_report"
    h, data = _page_data(out / "report.html")
    tab = out / "bundle" / "tables"
    SA = pd.read_csv(tab / "self_assessment.tsv", sep="\t", keep_default_na=False)
    for c in ("top2", "ll_top2", "margin", "runner_up"):
        assert c in SA.columns, c
    called = SA.best != U
    assert (SA.runner_up[called] != SA.best[called]).all()
    BR = pd.read_csv(tab / "boundary.tsv", sep="\t", keep_default_na=False)
    assert int(BR.cells.sum()) == int((called & (SA.unique == U)).sum())
    CPI = pd.read_csv(tab / "certainty_per_identity.tsv", sep="\t")
    assert CPI.cells.sum() == len(SA)
    # section 2: its own switchable UMAP (reference labels and the cross-validated calls; highlight off by default)
    assert [q["div"] for q in data["qumaps"]] == ["umap2", "qumap"]
    S2 = data["qumaps"][0]
    assert [o["name"] for o in S2["outputs"]][:3] == ["reference labels", "cross-validated best",
                                                      "cross-validated unique"]
    assert S2["views"] == ["label"] and '<option value="-1">(none)</option>' in h
    S = data["qumaps"][1]
    assert S["views"] == ["label", "margin", "ambiguity", "unlabelled"] and S["truth"] is not None
    assert all(set(o["shade"]) == {"margin", "ambiguity"} for o in S["outputs"])
    assert len(table_rows(data["tables"]["cells"])) == len(SA)
    assert "Certainty and boundary cells (cross-validation)" in h


def test_without_plotly_the_run_goes_on_and_the_report_is_skipped(tmp_path, monkeypatch, capsys):
    """A missing plotly never stops a run (report=True): a warning, the run's outputs written, no report."""
    from scarlet import api
    real = builtins.__import__

    def no_plotly(name, *a, **kw):
        if name == "plotly" or name.startswith("plotly."):
            raise ImportError("No module named 'plotly'")
        return real(name, *a, **kw)

    monkeypatch.setattr(builtins, "__import__", no_plotly)
    ref = make_synthetic_adata(seed=5)
    q = make_synthetic_adata(seed=7)
    api.transfer(ref, "cell_type", q, str(tmp_path / "t"), draws=7, report=True)
    assert (tmp_path / "t" / "labels.tsv").exists() and not (tmp_path / "t" / "report").exists()
    err = capsys.readouterr().err
    assert "the report will not be written" in err and "scarlet-sc[report]" in err
    api.learn(ref, "cell_type", str(tmp_path / "m.scarlet"), report=True)
    assert (tmp_path / "m.scarlet").exists() and not (tmp_path / "m_report").exists()
    assert "the report will not be written" in capsys.readouterr().err


def test_palette_is_colour_blind_safe_first():
    from scarlet.html_report import bundle as B
    c = B.colours(["b", "a", "c"])
    assert [c["a"], c["b"], c["c"]] == B.OKABE_ITO[:3]
    assert len(set(B.BASE)) == B.CVD_SAFE_MAX == 24 and B.REFUSED_COLOUR not in B.BASE
    assert "not colour-blind-safe" in B.palette_note(30) and "not" not in B.palette_note(10)


def test_documentation_covers_the_views():
    from pathlib import Path
    root = Path(__file__).resolve().parents[1]
    report = "\n".join(f.read_text() for f in sorted((root / "docs" / "reports").glob("*.md")))   # the report pages
    output = (root / "docs" / "output.md").read_text()
    gloss = (root / "docs" / "glossary.md").read_text()
    start = (root / "docs" / "getting-started.md").read_text()
    for k in ("Boundary cells", "certainty", "anchors", "unlabelled", "runner-up", "Per-cell data", "Download"):
        assert k in report, k
    for c in ("scarlet_top1", "scarlet_runner_up", "scarlet_runner_up_draws"):
        assert c in output, c
    for term in ("### Runner-up", "### Boundary cell"):
        assert term in gloss, term
    assert "scarlet-sc[report]" in start and "scarlet-sc[report]" in report


def test_pages_open_on_a_summary_and_end_on_provenance(runs):
    d, _ = runs
    for page in (d / "t" / "report" / "report.html", d / "m_report" / "report.html"):
        h, data = _page_data(page)
        assert h.index("<h2>Summary</h2>") < h.index("<h2>1. ") < h.index("Provenance")
        last = h.rindex("<h2>")
        assert "Provenance" in h[last:last + 40]   # the last section
        assert "colour-blind-safe" in h[h.rindex("Provenance"):]   # the palette note is a footnote
    h, data = _page_data(d / "t" / "report" / "report.html")
    assert "SCARLET transfer report: A &rarr; B" in h   # in-memory inputs: A and B; files: their stems
    assert "Depth matching, draw by draw" in h[h.index("Provenance"):]
    T = data["tables"]["cells"]
    assert len(T["groups"]) == T["n"] and T["group_names"][0].startswith("worth a look")
    assert sorted(T["groups"]) == list(T["groups"])   # the cells worth a look first
    assert "agree_all" in data["figs"] and "agree_id" in data["figs"]


def test_tables_read_well():
    from scarlet.html_report.html import cell_text, header, point_size, table_html
    assert cell_text("pct_boundary", 1.6148) == "1.6%" and cell_text("cv_recall_best", 0.9054) == "90.5%"
    assert cell_text("median_genes_detected", 1787.0) == "1,787" and cell_text("cells", 2539) == "2,539"
    assert cell_text("fdr", 5.155e-09) == "5.2e-09" and cell_text("median_fdr_best", 0.0) == "0"
    assert cell_text("seed", 1835504127) == "1835504127" and cell_text("flip", False) == "no"
    assert header("pct_unlabelled_best: background") == "unlabelled best: background"
    h = table_html(pd.DataFrame({"identity": ["a"], "runner_up": ["b"], "refused_flip test": [0]}))
    assert 'title="runner_up: for a called cell' in h and "flip test" not in h   # all-zero refusal column left out
    assert point_size(253) > point_size(6444) > point_size(68000) == 2


def test_svg_fallback_samples_within_each_identity():
    from scarlet.html_report.certainty import svg_subsample
    lab = np.array(["a"] * 20000 + ["b"] * 50 + ["c"] * 300, dtype=object)
    idx, t = svg_subsample(lab, cap=10000)
    kept = pd.Series(lab[idx]).value_counts()
    assert len(idx) <= 10000 and kept["b"] == 50 and kept["c"] == 300 and kept["a"] == t == 10000 - 350
    assert svg_subsample(lab[:5000], cap=10000) == (None, None)   # small data: every cell


def test_certainty_is_relative_to_the_run():
    from scarlet.html_report.certainty import MIN_INTENSITY, intensity_from_fdr
    fdr = 10.0 ** -np.arange(1, 101)    # -log10 FDR 1 .. 100
    called = np.ones(100, bool)
    v, lo, hi = intensity_from_fdr(fdr, called)
    assert abs(lo - np.percentile(np.arange(1, 101), 5)) < 1e-9 and abs(hi - np.percentile(np.arange(1, 101), 95)) < 1e-9
    assert v[0] == MIN_INTENSITY and v[-1] == 1.0 and len(set(np.round(v, 6))) > 50   # a spread, not saturated
    assert (np.diff(v) >= 0).all()


def test_margin_agreement_and_ambiguity_views():
    """The views that replace 'certainty': margin (learn report), agreement (transfer), ambiguity (both)."""
    from scarlet.html_report.certainty import (MIN_INTENSITY, intensity_ambiguous, intensity_from_margin,
                                               intensity_from_need)
    m = np.arange(1, 101, dtype=float)
    v, lo, hi = intensity_from_margin(m, np.ones(100, bool))
    assert v[0] == MIN_INTENSITY and v[-1] == 1.0 and (np.diff(v) >= 0).all() and len(set(np.round(v, 6))) > 50
    k = intensity_from_need([8, 9, 10, 5, None], 10, 8)      # 8 of 10 faintest, 10 of 10 full
    assert k[0] == MIN_INTENSITY and 0.5 < k[1] < 0.6 and k[2] == 1.0 and k[3] == MIN_INTENSITY and k[4] == MIN_INTENSITY
    assert (intensity_from_need([3], 1, 1) == 1.0).all()     # quick mode: one draw, no agreement
    a = intensity_ambiguous(["B", "B", U, "T"], ["B", U, U, "T"])
    assert list(a) == [MIN_INTENSITY, 1.0, MIN_INTENSITY, MIN_INTENSITY]


def test_known_labels_in_the_transfer_and_assign_reports(tmp_path):
    pytest.importorskip("plotly")
    from scarlet import api
    ref, q = make_synthetic_adata(seed=5), make_synthetic_adata(seed=7)
    q.obs["known"] = q.obs["cell_type"].astype(str)
    q.write_h5ad(tmp_path / "q.h5ad")
    ref.write_h5ad(tmp_path / "r.h5ad")
    api.transfer(str(tmp_path / "r.h5ad"), "cell_type", str(tmp_path / "q.h5ad"), str(tmp_path / "t"), n_iterations=20,
                 draws=3)
    with pytest.raises(ValueError, match="no obs column"):
        api.report(str(tmp_path / "t"), query_labels="nope", self_assessment=False, output=str(tmp_path / "bad"))
    page = api.report(str(tmp_path / "t"), query_labels="known", self_assessment=False, output=str(tmp_path / "rk"))
    h, data = _page_data(page)
    assert "Against labels supplied by the user" in h and "SCARLET did not read them" in h and "known_conf" in data["figs"]
    S = pd.read_csv(tmp_path / "rk" / "bundle" / "tables" / "known_summary.tsv", sep="\t")
    assert list(S.output) == ["final best", "final unique", "first pass best", "first pass unique"]
    assert S.median_f1_called.between(0, 1).all() and S.cells.eq(q.n_obs).all()
    plain = api.report(str(tmp_path / "t"), self_assessment=False, output=str(tmp_path / "r0"))
    assert "Against labels supplied by the user" not in open(plain).read()   # absent by default
    api.learn(str(tmp_path / "r.h5ad"), "cell_type", str(tmp_path / "m.scarlet"), n_iterations=20)
    with pytest.raises(ValueError, match="transfer or an assign report"):
        api.report(str(tmp_path / "m.scarlet"), query_labels="known", self_assessment=False)
    api.assign(str(tmp_path / "m.scarlet"), str(tmp_path / "q.h5ad"), str(tmp_path / "a.h5ad"))
    page = api.report(str(tmp_path / "a.h5ad"), model=str(tmp_path / "m.scarlet"), query_labels="known",
                      output=str(tmp_path / "ra"))
    assert "Against labels supplied by the user" in open(page).read()


def test_umaps_in_scarlets_space_and_on_all_genes(runs, tmp_path):
    """Default space: the model's universe genes as detection (Jaccard); all genes as the second space, unless
    umap_all_genes is off (then the report says so). Every UMAP names its space."""
    from scarlet import api
    d, _ = runs
    for page, keys in ((d / "m_report" / "report.html", ["umap2", "qumap"]),
                       (d / "t" / "report" / "report.html", ["umapA", "qumap"])):
        h, data = _page_data(page)
        for S in data["qumaps"]:
            names = [s["name"] for s in S["spaces"]]
            assert names[0].startswith("universe genes of") and "Jaccard" in names[0] and names[1].startswith("all genes")
            assert len(S["spaces"][1]["x"]) == len(S["x"])
        assert "the place to look for populations the model lacks" in h and 'id="qumap_space"' in h
    U = pd.read_csv(d / "m_report" / "bundle" / "tables" / "umap.tsv", sep="\t")
    assert {"x", "y", "x_all", "y_all"} <= set(U.columns)
    page = api.report(str(d / "m.scarlet"), reference=make_synthetic_adata(seed=5), self_assessment=False,
                      umap_all_genes=False, output=str(tmp_path / "r1"))   # the fixture learned from memory
    h, data = _page_data(page)
    assert [len(S["spaces"]) for S in data["qumaps"]] == [1] and "umap_all_genes = no" in h
