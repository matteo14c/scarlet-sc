"""learn and transfer record the sha256 of their input files: learn in the model's provenance (inputs.reference),
transfer in transfer_record.json (inputs.reference and inputs.query, with the cells and genes as read); an
in-memory AnnData has none."""
import hashlib
import json

import pytest

from scarlet import api, load_model
from conftest import make_synthetic_adata

KW = dict(n_iterations=6, n_cells=20)


def _sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def test_learn_records_the_reference_sha256(tmp_path):
    A = make_synthetic_adata(seed=5)
    p = tmp_path / "ref.h5ad"
    A.write_h5ad(p)
    api.learn(str(p), "cell_type", str(tmp_path / "m.scarlet"), **KW)
    ref = load_model(str(tmp_path / "m.scarlet")).provenance["inputs"]["reference"]
    assert ref["sha256"] == _sha(p) and ref["path"] == str(p)
    api.learn(A, "cell_type", str(tmp_path / "mem.scarlet"), **KW)
    assert load_model(str(tmp_path / "mem.scarlet")).provenance["inputs"]["reference"]["sha256"] is None


def test_transfer_records_both_sha256(tmp_path):
    R, Q = make_synthetic_adata(seed=5), make_synthetic_adata(seed=6)
    pr, pq = tmp_path / "ref.h5ad", tmp_path / "qry.h5ad"
    R.write_h5ad(pr)
    Q.write_h5ad(pq)
    api.transfer(str(pr), "cell_type", str(pq), str(tmp_path / "out"), draws=2, **KW)
    rec = json.loads((tmp_path / "out" / "transfer_record.json").read_text())["inputs"]
    assert rec["reference"] == {"path": str(pr), "sha256": _sha(pr), "n_cells": R.n_obs, "n_genes": R.n_vars}
    assert rec["query"] == {"path": str(pq), "sha256": _sha(pq), "n_cells": Q.n_obs, "n_genes": Q.n_vars}


def test_the_learn_report_checks_the_recorded_hash(tmp_path):
    """The report states the hash the run recorded and whether the file it reads is the same one."""
    pytest.importorskip("plotly")
    A = make_synthetic_adata(seed=5)
    p = tmp_path / "ref.h5ad"
    A.write_h5ad(p)
    api.learn(str(p), "cell_type", str(tmp_path / "m.scarlet"), **KW)
    api.report(str(tmp_path / "m.scarlet"), self_assessment=False, output=str(tmp_path / "r1"))
    ref = json.load(open(tmp_path / "r1" / "bundle" / "manifest.json"))["meta"]["provenance"]["input"]["reference"]
    assert ref["sha256"] == _sha(p) and ref["check"] == "recorded by the run; the file the report read is the same"
    A.obs["extra"] = 1   # another file with the same cells and genes
    A.write_h5ad(p)
    api.report(str(tmp_path / "m.scarlet"), self_assessment=False, output=str(tmp_path / "r2"))
    ref = json.load(open(tmp_path / "r2" / "bundle" / "manifest.json"))["meta"]["provenance"]["input"]["reference"]
    assert "differs: it is not the run's input" in ref["check"]


def test_an_empty_git_commit_is_not_shown():
    from scarlet.html_report import html
    assert html._rows([("git commit", ""), ("git commit", None), ("seed", 42)]) == [("seed", 42)]
    assert html._rows([("git commit", "abc123")]) == [("git commit", "abc123")]
