"""The public interface: exactly learn, assign, transfer and report; the command line and Python take the same
options with the same defaults and give the same outputs; values outside the documented ones are refused."""
from __future__ import annotations

import inspect
import json
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pytest

import scarlet
import scarlet.cli as cli
import scarlet.marker_finder as mf
from scarlet import api
from scarlet.model import load_model

FAST = ["--n-iterations", "5", "--n-cells", "30"]
FAST_API = dict(n_iterations=5, n_cells=30)
POSITIONAL = {("report", "run")}
OFF_SWITCHES = {("report", "self_assessment"): "--no-self-assessment"}


def run_cli(argv):
    old = sys.argv
    sys.argv = ["scarlet"] + list(argv)
    try:
        cli.main(list(argv))
    finally:
        sys.argv = old


def members(path, names=("markers.tsv", "universe.tsv", "background_universe.tsv", "programmes.json")):
    m = load_model(path)
    return {n: m.files[n] for n in names}


# ---------------------------------------------------------------------------------------------- the public surface
def test_the_public_api_is_learn_assign_transfer_report():
    assert list(api.COMMANDS) == ["learn", "assign", "transfer", "report"]
    assert set(scarlet.__all__) == {"learn", "assign", "transfer", "report", "load_model", "Model", "datasets"}


def test_the_package_never_imports_igraph_or_leidenalg():
    root = Path(scarlet.__file__).parent
    for p in root.rglob("*.py"):
        text = p.read_text()
        assert "igraph" not in text and "leidenalg" not in text, p


@pytest.mark.parametrize("name", list(api.COMMANDS))
def test_cli_and_api_are_symmetric(name):
    """Every command-line option is the Python parameter of the same name, with the same default; each has help."""
    sig = inspect.signature(api.COMMANDS[name])
    sub = cli.build_parser()._subparsers._group_actions[0].choices[name]
    actions = {a.dest: a for a in sub._actions if a.dest != "help"}
    assert set(actions) == set(sig.parameters), set(actions) ^ set(sig.parameters)
    for pname, p in sig.parameters.items():
        act = actions[pname]
        assert act.help and len(act.help) > 20, (name, pname)
        if (name, pname) in POSITIONAL:
            assert act.option_strings == [] and p.default is inspect.Parameter.empty
            continue
        assert act.option_strings == [OFF_SWITCHES.get((name, pname), "--" + pname.replace("_", "-"))]
        if p.default is inspect.Parameter.empty:
            assert act.required
        else:
            assert not act.required and act.default == p.default, (pname, act.default, p.default)


@pytest.mark.parametrize("name", [None] + list(api.COMMANDS))
def test_cli_help_prints(name, capsys):
    """scarlet --help and every command's --help print and exit 0; a percent sign in a description is shown as
    written (transfer's --draws: "80 % of the draws")."""
    with pytest.raises(SystemExit) as e:
        cli.main(([name] if name else []) + ["--help"])
    assert e.value.code == 0
    out = capsys.readouterr().out
    assert out.startswith("usage: scarlet")
    if name == "transfer":
        assert "80 % of the draws" in " ".join(out.split())


def test_cli_parses_to_the_api_call(monkeypatch):
    import functools
    seen = {}

    @functools.wraps(api.transfer)
    def fake(**kw):
        seen.update(kw)

    monkeypatch.setitem(api.COMMANDS, "transfer", fake)
    cli.main(["transfer", "--reference", "r.h5ad", "--labels", "cell_type", "--query", "q.h5ad", "--output", "o",
              "--refine", "no", "--depth-match", "no", "--n-markers", "9"])
    defaults = {k: p.default for k, p in inspect.signature(api.transfer).parameters.items()}
    assert seen["refine"] is False and seen["depth_match"] is False and seen["n_markers"] == 9
    rest = ("reference", "labels", "query", "output", "refine", "depth_match", "n_markers")
    assert {k: v for k, v in seen.items() if k not in rest} == {k: v for k, v in defaults.items() if k not in rest}


@pytest.mark.parametrize("kw", [dict(technology="other"), dict(imbalanced="maybe"), dict(n_iterations=0),
                                dict(blocklist="no-such-list")])
def test_learn_refuses_values_outside_the_documented_ones(synthetic_h5ad, tmp_path, kw):
    with pytest.raises((ValueError, FileNotFoundError)):
        api.learn(str(synthetic_h5ad), "cell_type", str(tmp_path / "m.scarlet"), **kw)


@pytest.mark.parametrize("verb", ["learn", "assign", "transfer", "report"])
def test_unknown_options_are_refused(verb):
    with pytest.raises(TypeError):
        getattr(api, verb)("a", "b", "c", no_such_option=1) if verb != "report" else api.report("a", no_such_option=1)


def test_marker_shortlist_replaces_n_background():
    for name in ("learn", "transfer"):
        params = inspect.signature(api.COMMANDS[name]).parameters
        assert "marker_shortlist" in params and "n_background" not in params
        assert params["marker_shortlist"].default is None             # the technology preset's (150)
    with pytest.raises(SystemExit):
        cli.build_parser().parse_args(["learn", "--reference", "r", "--labels", "l", "--model", "m",
                                       "--n-background", "20"])


def test_assign_and_transfer_refuse_flip_values_outside_the_documented_ones(synthetic_h5ad, tmp_path):
    with pytest.raises(ValueError):
        api.assign("m.scarlet", str(synthetic_h5ad), str(tmp_path / "a.h5ad"), flip="off")
    with pytest.raises(ValueError):
        api.transfer(str(synthetic_h5ad), "cell_type", str(synthetic_h5ad), str(tmp_path / "t"), flip="off")


# ------------------------------------------------------------------------------------- the same outputs, CLI = API
def test_learn_cli_equals_api(synthetic_h5ad, tmp_path):
    run_cli(["learn", "--reference", str(synthetic_h5ad), "--labels", "cell_type", "--model",
             str(tmp_path / "cli.scarlet"), "--n-markers", "7", "--marker-shortlist", "20", "--seed", "3"] + FAST)
    api.learn(str(synthetic_h5ad), "cell_type", str(tmp_path / "api.scarlet"), n_markers=7, marker_shortlist=20, seed=3,
              **FAST_API)
    assert members(tmp_path / "cli.scarlet") == members(tmp_path / "api.scarlet")


def test_learn_anndata_input_equals_file_input(synthetic_h5ad, synthetic_adata, tmp_path):
    api.learn(str(synthetic_h5ad), "cell_type", str(tmp_path / "f.scarlet"), **FAST_API)
    api.learn(synthetic_adata, "cell_type", str(tmp_path / "a.scarlet"), **FAST_API)
    assert members(tmp_path / "f.scarlet") == members(tmp_path / "a.scarlet")


def test_learn_dense_cli_equals_api(synthetic_h5ad, tmp_path):
    run_cli(["learn", "--reference", str(synthetic_h5ad), "--labels", "cell_type", "--model",
             str(tmp_path / "cli.scarlet"), "--technology", "dense"] + FAST)
    api.learn(str(synthetic_h5ad), "cell_type", str(tmp_path / "api.scarlet"), technology="dense", **FAST_API)
    assert members(tmp_path / "cli.scarlet") == members(tmp_path / "api.scarlet")


def test_assign_cli_equals_api(synthetic_h5ad, tmp_path):
    api.learn(str(synthetic_h5ad), "cell_type", str(tmp_path / "m.scarlet"), **FAST_API)
    run_cli(["assign", "--model", str(tmp_path / "m.scarlet"), "--data", str(synthetic_h5ad), "--output",
             str(tmp_path / "cli.h5ad"), "--predictions", str(tmp_path / "cli.tsv")])
    lab = api.assign(str(tmp_path / "m.scarlet"), str(synthetic_h5ad), str(tmp_path / "api.h5ad"),
                     predictions=str(tmp_path / "api.tsv"))
    assert (tmp_path / "cli.tsv").read_bytes() == (tmp_path / "api.tsv").read_bytes()
    a, b = ad.read_h5ad(tmp_path / "cli.h5ad"), ad.read_h5ad(tmp_path / "api.h5ad")
    pd.testing.assert_frame_equal(a.obs, b.obs)
    assert np.array_equal(a.obsm["scarlet_probs"], b.obsm["scarlet_probs"])
    assert list(lab.columns) == ["scarlet_best", "scarlet_unique", "is_scarlet_unique", "scarlet_fdr", "scarlet_flip",
                                 "scarlet_reason", "scarlet_pe_low", "scarlet_pe_note"]
    assert (lab.scarlet_unique[~lab.is_scarlet_unique] == "unassigned").all()
    assert (tmp_path / "api.summary.txt").read_text().startswith("SCARLET ")


def test_transfer_cli_equals_api(synthetic_h5ad, synthetic_adata, tmp_path):
    q = synthetic_adata.copy()
    q.obs["batch"] = np.where(np.arange(q.n_obs) % 2, "b1", "b2")
    qp = tmp_path / "q.h5ad"
    q.write_h5ad(qp)
    run_cli(["transfer", "--reference", str(synthetic_h5ad), "--labels", "cell_type", "--query", str(qp),
             "--output", str(tmp_path / "cli"), "--origin-col", "batch", "--draws", "7"]
            + FAST)
    res = api.transfer(str(synthetic_h5ad), "cell_type", str(qp), str(tmp_path / "api"), origin_col="batch",
                       draws=7, **FAST_API)
    for name in ("labels.tsv", "first_pass.tsv", "qc.tsv"):
        assert (tmp_path / "cli" / name).read_bytes() == (tmp_path / "api" / name).read_bytes(), name
    for name in ("reference.scarlet", "refine_b1.scarlet", "refine_b2.scarlet"):
        assert members(tmp_path / "cli" / name) == members(tmp_path / "api" / name), name
    assert set(res.query_models) == {"b1", "b2"} and res.qc is not None
    assert (tmp_path / "api" / "summary.txt").exists()
    rec = json.load(open(tmp_path / "api" / "transfer_record.json"))
    assert rec["summary"] == (tmp_path / "api" / "summary.txt").read_text()


# --------------------------------------------------------------------------------------------- presets and records
@pytest.mark.parametrize("technology, gamma, n_markers, marker_shortlist",
                         [("dense", 0.0, 50, 150), ("sparse", 6.0, 50, 150)])
def test_technology_presets(synthetic_h5ad, tmp_path, technology, gamma, n_markers, marker_shortlist):
    m = api.learn(str(synthetic_h5ad), "cell_type", str(tmp_path / "m.scarlet"), technology=technology, **FAST_API)
    s = m.settings["resolved"]
    assert (s["technology"]["name"], s["technology"]["marker_gamma"], s["sizes"]["n_markers"],
            s["sizes"]["marker_shortlist"]) == (technology, gamma, n_markers, marker_shortlist)
    p = json.loads(m.files["marker_finder_params.json"])
    assert p["context"]["technology"] == technology
    assert "auto_binarize" not in p["context"] and s["technology"]["detection"] == "count > 0"


def test_explicit_sizes_win(synthetic_h5ad, tmp_path):
    m = api.learn(str(synthetic_h5ad), "cell_type", str(tmp_path / "m.scarlet"), n_markers=25, marker_shortlist=75,
                  **FAST_API)
    assert (m.settings["resolved"]["sizes"]["n_markers"], m.settings["resolved"]["sizes"]["marker_shortlist"]) == (25, 75)


def test_model_members_and_plain_records(synthetic_h5ad, tmp_path):
    import zipfile
    m = api.learn(str(synthetic_h5ad), "cell_type", str(tmp_path / "m.scarlet"), **FAST_API)
    with zipfile.ZipFile(tmp_path / "m.scarlet") as z:
        names = set(z.namelist())
    assert names == {"manifest.json", "markers.tsv", "programmes.json", "universe.tsv", "background_universe.tsv",
                     "marker_finder_params.json", "settings.json", "provenance.json", "positive_evidence.json",
                     "summary.txt", "in_class.tsv"}
    assert m.settings["resolved"]["imbalanced"]["handling"] in ("standard", "balanced for unequal class sizes")
    assert m.settings["parameters"]["n_iterations"] == 5
    assert (tmp_path / "m.summary.txt").read_text() == m.files["summary.txt"].decode()
    bad = tmp_path / "bad.scarlet"
    with zipfile.ZipFile(tmp_path / "m.scarlet") as zi, zipfile.ZipFile(bad, "w") as zo:
        for n in zi.namelist():
            zo.writestr(n, zi.read(n) + b"x" if n == "universe.tsv" else zi.read(n))
    with pytest.raises(ValueError):
        load_model(bad)


# ------------------------------------------------------------------------------- options not at their defaults
def test_learn_records_and_states_its_options_not_at_their_defaults(synthetic_h5ad, tmp_path):
    m = api.learn(str(synthetic_h5ad), "cell_type", str(tmp_path / "m.scarlet"), **FAST_API)
    nd = m.settings["parameters"]["options_not_default"]
    assert nd == {"n_iterations": {"value": 5, "default": 100}, "n_cells": {"value": 30, "default": 100}}
    text = m.files["summary.txt"].decode()
    assert "Options not at their defaults: n_iterations 5 (default: 100); n_cells 30 (default: 100)." in text
    m2 = api.learn(str(synthetic_h5ad), "cell_type", None, n_markers=7, marker_shortlist=20, seed=3, ncpu=2,
                   **FAST_API)
    nd2 = m2.settings["parameters"]["options_not_default"]
    assert set(nd2) == {"n_markers", "marker_shortlist", "seed", "n_iterations", "n_cells"}    # ncpu changes no call
    assert nd2["n_markers"] == {"value": 7, "default": None}
    assert "n_markers 7 (default: the technology preset's)" in m2.files["summary.txt"].decode()
    assert (m2.settings["resolved"]["sizes"]["n_markers"], m2.settings["resolved"]["sizes"]["marker_shortlist"]) == \
        (7, 20)


def test_assign_records_alpha_and_its_options(synthetic_h5ad, tmp_path):
    api.learn(str(synthetic_h5ad), "cell_type", str(tmp_path / "m.scarlet"), **FAST_API)
    a = api.assign(str(tmp_path / "m.scarlet"), str(synthetic_h5ad), str(tmp_path / "a.h5ad"))
    b = api.assign(str(tmp_path / "m.scarlet"), str(synthetic_h5ad), str(tmp_path / "b.h5ad"), alpha=0.2)
    A, Bz = ad.read_h5ad(tmp_path / "a.h5ad"), ad.read_h5ad(tmp_path / "b.h5ad")
    assert float(A.uns["scarlet_options"]["alpha"]) == 0.01
    assert json.loads(A.uns["scarlet_options"]["options_not_default"]) == {}
    assert float(Bz.uns["scarlet_options"]["alpha"]) == 0.2
    assert json.loads(Bz.uns["scarlet_options"]["options_not_default"]) == {"alpha": {"value": 0.2, "default": 0.01}}
    assert "All options at their defaults." in str(A.uns["scarlet_summary"])
    assert "Options not at their defaults: alpha 0.2 (default: 0.01)." in str(Bz.uns["scarlet_summary"])
    pd.testing.assert_series_equal(a.scarlet_best, b.scarlet_best)             # alpha changes no best call
    assert b.is_scarlet_unique.sum() >= a.is_scarlet_unique.sum()


def test_transfer_records_its_options_and_the_cell_cycle_decision(synthetic_h5ad, tmp_path):
    api.transfer(str(synthetic_h5ad), "cell_type", str(synthetic_h5ad), str(tmp_path / "t"), alpha=0.05,
                 refine=False, draws=7, **FAST_API)
    rec = json.load(open(tmp_path / "t" / "transfer_record.json"))
    assert rec["options_not_default"]["alpha"] == {"value": 0.05, "default": 0.01}
    assert rec["options_not_default"]["refine"] == {"value": False, "default": True}
    assert rec["cell_cycle"]["reference"] == "blocked"                        # labels A, B, C: none names dividing cells
    assert all(e["cell_cycle_genes"] == "blocked" for e in rec["cell_cycle"]["learns"])
    text = (tmp_path / "t" / "summary.txt").read_text()
    assert "Cell-cycle genes: blocked." in text and "alpha 0.05 (default: 0.01)" in text


def test_the_cell_cycle_statement(synthetic_adata, tmp_path):
    from scarlet.run_summary import cell_cycle_text
    A = synthetic_adata.copy()
    A.obs["cell_type"] = A.obs["cell_type"].astype(str).replace({"A": "cycling T cells", "B": "proliferating B"})
    m = api.learn(A, "cell_type", None, **FAST_API)
    assert cell_cycle_text(m.settings["resolved"]["blocklist"]) == \
        "not blocked, because the labels 'cycling T cells', 'proliferating B' name dividing cells"
    assert "Cell-cycle genes: not blocked, because the labels" in m.files["summary.txt"].decode()
    m = api.learn(A, "cell_type", None, blocklist_modules="all", **FAST_API)
    assert cell_cycle_text(m.settings["resolved"]["blocklist"]) == "blocked"
    m = api.learn(A, "cell_type", None, blocklist="none", **FAST_API)
    assert cell_cycle_text(m.settings["resolved"]["blocklist"]) == "not blocked (no blocklist applied)"
