"""The transfer report's pages after the third review: section 3's bars without numbers (the assign report's too),
the drill-down in three parts read from the models' in-class tables, the ranks in both models, 'what refinement
changed' on a page of its own, and the per-cell data counted first with the table folded and paged."""
from __future__ import annotations

import json
import re

import numpy as np
import pandas as pd
import pytest

from conftest import make_synthetic_adata

U = "unassigned"
PARTS = ["common to both models", "A only (dropped by refinement)", "B only (added by refinement)"]


def _page(path):
    h = open(path, encoding="utf-8").read()
    m = re.search(r"<script type='application/json' id='report-data'>(.*?)</script>", h, re.S)
    return h, json.loads(m.group(1).replace("<\\/", "</"))


@pytest.fixture(scope="module")
def run(tmp_path_factory):
    pytest.importorskip("plotly")
    from scarlet import api
    d = tmp_path_factory.mktemp("pages")
    ref = make_synthetic_adata(seed=5)
    q = make_synthetic_adata(seed=7)
    q.obs = pd.DataFrame(index=q.obs_names)
    api.transfer(ref, "cell_type", q, str(d / "t"), n_iterations=20, draws=7, report=True)
    api.learn(ref, "cell_type", str(d / "m.scarlet"), n_iterations=20)
    api.assign(str(d / "m.scarlet"), q, str(d / "q_annotated.h5ad"))
    api.report(str(d / "q_annotated.h5ad"), model=str(d / "m.scarlet"))
    return d, ref, q


def test_bars_without_numbers(run):
    d, _, _ = run
    for page in (d / "t" / "report" / "report.html", d / "q_annotated_report" / "report.html"):
        _, data = _page(page)
        bars = [t for t in data["figs"]["analytics"]["data"] if t["type"] == "bar"]
        assert bars and all(not t.get("text") for t in bars), page
        assert all("of the cells" in t["hovertemplate"] and len(t["customdata"]) == len(t["x"]) for t in bars)
        assert "hover for the count" in data["figs"]["analytics"]["layout"]["title"]["text"]


def test_models_keep_the_in_class_table(run):
    """learn's model and transfer's saved models hold the full in-class table (every gene the marker search scored);
    its marker rows are the universe's rows, as text; the refined model's bg is its pooled background in both."""
    from scarlet.model import load_model
    d, ref, _ = run
    for p in (d / "m.scarlet", d / "t" / "reference.scarlet", d / "t" / "refine_all.scarlet"):
        m = load_model(p)
        U = m.files["universe.tsv"].decode().splitlines()
        rows = {r.split("\t")[0]: r for r in m.files["in_class.tsv"].decode().splitlines()[1:]}
        assert m.files["in_class.tsv"].decode().splitlines()[0] == U[0]
        assert all(rows[r.split("\t")[0]] == r for r in U[1:]), p.name
        assert len(rows) > len(U) - 1 and list(m.in_class.columns) == list(m.universe.columns)


def test_drill_down_in_three_parts_from_the_models(run):
    """every drill-down value is the model's own in-class estimate; no blank when every gene was scored"""
    from scarlet.model import load_model
    d, _, _ = run
    out = d / "t" / "report"
    _, data = _page(out / "report.html")
    mk = pd.read_csv(out / "bundle" / "tables" / "markers.tsv", sep="\t", keep_default_na=False)
    inc = {n: load_model(d / "t" / n).in_class for n in ("reference.scarlet", "refine_all.scarlet")}
    refined = 0
    for k, dr in data["drill"].items():
        A_, B_ = dr["heat"]
        if not B_["genes"]:
            continue
        refined += 1
        ra = set(mk[(mk.model == "reference.scarlet") & (mk.cell_type == k)].gene.astype(str))
        rq = set(mk[mk.model.str.startswith("refine_") & (mk.cell_type == k)].gene.astype(str))
        for h, model in ((A_, "reference.scarlet"), (B_, "refine_all.scarlet")):
            assert h["split"] and h["genes"] == A_["genes"]
            names = [g[0] for g in h["groups"]]
            assert names == [p for p in PARTS if p in names] and len(names) >= 1
            bounds = {g[0]: h["genes"][g[1]:g[2]] for g in h["groups"]}
            assert set(bounds.get(PARTS[0], [])) == ra & rq
            assert set(bounds.get(PARTS[1], [])) == ra - rq and set(bounds.get(PARTS[2], [])) == rq - ra
            assert all(v is not None for row in h["z"] for v in row) and not h["gaps"]   # no blank
            F = inc[model]
            for i, g in enumerate(h["genes"]):
                for j, c in enumerate(h["cols"]):
                    assert abs(h["z"][i][j] - float(F.loc[g, c])) <= 5e-5 + 1e-12   # the value, to 4 decimals
        r = dr["rank"][0]
        assert {x[0] for x in r["common"]} == ra & rq and {x[0] for x in r["a_only"]} == ra - rq
        assert {x[0] for x in r["b_only"]} == rq - ra
        assert sorted([x[1] for x in r["common"]] + [x[1] for x in r["a_only"]]) == list(range(1, len(ra) + 1))
    assert refined >= 1


def test_models_without_the_table_leave_blanks(run, tmp_path):
    """models written before 0.6.3 (no in_class member): the other model's markers are blank, 'not a marker of this
    model' in the hover; the report computes nothing from the data"""
    import shutil
    from scarlet import api
    from scarlet.model import load_model
    d, ref, q = run
    t = tmp_path / "t"
    shutil.copytree(d / "t", t, ignore=shutil.ignore_patterns("report"))
    for n in ("reference.scarlet", "refine_all.scarlet"):
        m = load_model(t / n)
        del m.files["in_class.tsv"]
        m.save(t / n)
    api.report(str(t), reference=ref, query=q)
    h, data = _page(t / "report" / "report.html")
    assert not (t / "report" / "bundle" / "tables" / "in_class.tsv").exists()
    assert "not a marker of this model" in h
    blanks = 0
    for k, dr in data["drill"].items():
        A_, B_ = dr["heat"]
        if not B_["genes"]:
            continue
        for h_ in (A_, B_):
            for i, row in enumerate(h_["z"]):
                for j, v in enumerate(row):
                    if v is None:
                        blanks += 1
                        assert h_["gaps"][i][j] == 1
    assert blanks > 0


def test_refinement_on_its_own_page(run):
    d, _, _ = run
    out = d / "t" / "report"
    h, _ = _page(out / "report.html")
    six = h[h.index("<h2>6. "):h.index("<h2>7. ")]
    assert 'href="refinement.html"' in six and "<table>" not in six
    hr, dr = _page(out / "refinement.html")
    assert 'href="report.html"' in hr and "First pass and refined, side by side" in hr
    X6 = pd.read_csv(out / "bundle" / "tables" / "refinement_changes.tsv", sep="\t", keep_default_na=False)
    lst = hr[hr.index("<details>"):hr.index("</details>")]
    assert lst.count("<tr>") == len(X6) + 1   # header + every change
    for e in ("best", "unique"):
        if (X6.entry == e).any():
            assert f"chg_{e}" in dr["figs"]


def test_per_cell_data_counted_first_then_folded(run):
    d, _, _ = run
    h, data = _page(d / "t" / "report" / "report.html")
    eight = h[h.index("<h2>8. "):h.index("<h2>9. ")]
    T = data["tables"]["cells"]
    assert T["page"] == 25 and 'id="cells_more"' in eight
    assert eight.index("class='tiles'") < eight.index("Download all") < eight.index("<details><summary>Browse the cells")
    nums = [int(x.replace(",", "")) for x in re.findall(r'<div class="tile"><b>([\d,]+)</b>', eight)]
    assert sum(nums) == T["n"]


def test_documentation_names_the_pages():
    from pathlib import Path
    doc = (Path(__file__).resolve().parents[1] / "docs" / "reports" / "transfer.md").read_text()
    for k in ("refinement.html", "three parts", "in_class", "Browse the cells", "25 rows at a time"):
        assert k in doc, k
