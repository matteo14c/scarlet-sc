"""Report defects found on rc3 by a run queue on atlases (gene_col, repeated gene names, a report stopping a run)."""
import numpy as np
import pandas as pd
import pytest

from conftest import make_synthetic_adata


def test_self_assessment_reads_the_genes_as_learn_did(tmp_path):
    """A model learned with gene_col: the folds' assign must read the same column (before: the var index, so no
    programme gene was found and every cell was refused)."""
    pytest.importorskip("plotly")
    from scarlet import api
    A = make_synthetic_adata(seed=5)
    A.var["symbol"] = A.var_names.astype(str)
    A.var_names = [f"ENSG{i:011d}" for i in range(A.n_vars)]
    A.write_h5ad(tmp_path / "r.h5ad")
    api.learn(str(tmp_path / "r.h5ad"), "cell_type", str(tmp_path / "m.scarlet"), gene_col="symbol", n_iterations=20)
    api.report(str(tmp_path / "m.scarlet"), output=str(tmp_path / "rep"), umap_all_genes=False)
    sa = pd.read_csv(tmp_path / "rep" / "bundle" / "tables" / "self_assessment.tsv", sep="\t")
    assert (sa.unique != "unassigned").mean() > 0.5


def test_repeated_gene_names_do_not_stop_the_report(tmp_path):
    pytest.importorskip("plotly")
    from scarlet import api
    A = make_synthetic_adata(seed=5)
    names = A.var_names.astype(str).tolist()
    names[1] = names[0]                     # one repeated name, as in some atlases
    A.var_names = names
    A.write_h5ad(tmp_path / "r.h5ad")
    api.learn(str(tmp_path / "r.h5ad"), "cell_type", str(tmp_path / "m.scarlet"), n_iterations=20)
    out = api.report(str(tmp_path / "m.scarlet"), output=str(tmp_path / "rep"), umap_all_genes=False)
    assert (tmp_path / "rep" / "report.html").exists() and out


def test_a_failing_report_is_a_warning_not_a_stop(tmp_path, monkeypatch, capsys):
    from scarlet import api
    A = make_synthetic_adata(seed=5)
    A.write_h5ad(tmp_path / "r.h5ad")
    monkeypatch.setattr(api, "_report_available", lambda: True)
    def broken(*a, **k):
        raise RuntimeError("broken report")
    monkeypatch.setattr(api, "report", broken)
    m = api.learn(str(tmp_path / "r.h5ad"), "cell_type", str(tmp_path / "m.scarlet"), n_iterations=20, report=True)
    assert m is not None and (tmp_path / "m.scarlet").exists()
    assert "the report of" in capsys.readouterr().err
