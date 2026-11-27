"""The verbosity option of learn, assign, transfer and report: normal (the default) prints the command's pre-flight
once, one progress bar and one line per fold or draw; quiet prints nothing; full prints everything each step prints.
It changes no output file and is not recorded as a setting. Captured at the file-descriptor level, so the draw
workers' output (new Python processes) is included."""
import json

import pytest

from scarlet import api, load_model
from conftest import make_synthetic_adata

KW = dict(n_iterations=6, n_cells=20)
FILES = ("labels.tsv", "first_pass.tsv", "qc.tsv")


@pytest.fixture(scope="module")
def pair(tmp_path_factory):
    d = tmp_path_factory.mktemp("verbosity")
    make_synthetic_adata(seed=5).write_h5ad(d / "ref.h5ad")
    make_synthetic_adata(seed=6).write_h5ad(d / "qry.h5ad")
    return d


def _transfer(pair, name, capfd, **kw):
    out = pair / name
    api.transfer(str(pair / "ref.h5ad"), "cell_type", str(pair / "qry.h5ad"), str(out), draws=3, **KW, **kw)
    return out, capfd.readouterr()


@pytest.mark.parametrize("ncpu", [1, 2])
def test_transfer_levels_print_what_they_say_and_change_no_output(pair, capfd, ncpu):
    runs = {v: _transfer(pair, f"{v}{ncpu}", capfd, verbosity=v, ncpu=ncpu) for v in ("quiet", "normal", "full")}
    quiet, normal, full = (runs[v][1].out for v in ("quiet", "normal", "full"))
    assert quiet == ""
    assert normal.count("[SCARLET pre-flight]") == 1 and "[SCARLET pre-flight] transfer" in normal
    for step in ("first pass", "refinement"):
        assert all(f"{step}, draw {d} of 3:" in normal for d in (1, 2, 3)), normal
    assert "transfer finished:" in normal and "Result:" in normal
    assert "+ scarlet.marker_finder" not in normal and "Bootstrap" not in normal
    assert "[SCARLET pre-flight] assign" in full and "+ scarlet.marker_finder" in full
    assert all(f"first pass, draw {d} of 3:" in full for d in (1, 2, 3))
    for f in FILES:
        assert len({(runs[v][0] / f).read_bytes() for v in runs}) == 1, f
    recs = [json.loads((runs[v][0] / "transfer_record.json").read_text()) for v in runs]
    for k in ("options_not_default", "settings", "draws", "inputs"):
        assert recs[0][k] == recs[1][k] == recs[2][k], k
    assert all("verbosity" not in r["parameters"] for r in recs)
    assert "Bootstrap" not in runs["normal"][1].err and "Bootstrap" not in runs["quiet"][1].err


def test_learn_and_assign(pair, tmp_path, capfd):
    m = {}
    for v in ("quiet", "normal", "full"):
        m[v] = tmp_path / f"{v}.scarlet"
        api.learn(str(pair / "ref.h5ad"), "cell_type", str(m[v]), verbosity=v, **KW)
        out, err = capfd.readouterr()
        if v == "quiet":
            assert out == "" and "Bootstrap" not in err
        elif v == "normal":
            assert out.count("[SCARLET pre-flight] learn") == 1 and "learn finished:" in out
            assert "+ scarlet.marker_finder" not in out and "Bootstrap" in err   # learn's one progress bar
        else:
            assert "+ scarlet.marker_finder" in out
    loaded = [load_model(str(p)) for p in m.values()]   # the same model (the files differ by path and time only)
    for part in ("markers", "universe", "background"):
        assert all(getattr(x, part).equals(getattr(loaded[0], part)) for x in loaded), part
    assert all("verbosity" not in x.settings["parameters"] for x in loaded)
    lab = {}
    for v in ("quiet", "normal", "full"):
        lab[v] = api.assign(str(m["normal"]), str(pair / "qry.h5ad"), None, verbosity=v)
        out = capfd.readouterr().out
        if v == "quiet":
            assert out == ""
        elif v == "normal":
            assert out.count("[SCARLET pre-flight] assign") == 1 and "Done. Unique calls" in out
            assert "Computing likelihoods" not in out
        else:
            assert "Computing likelihoods" in out
    assert lab["quiet"].equals(lab["normal"]) and lab["normal"].equals(lab["full"])


def test_a_command_inside_another_keeps_its_level(pair, tmp_path, capfd):
    """learn(report=True): the report's fold learns and assigns print nothing of their own at normal."""
    pytest.importorskip("plotly")
    api.learn(str(pair / "ref.h5ad"), "cell_type", str(tmp_path / "m.scarlet"), report=True, **KW)
    out, err = capfd.readouterr()
    assert out.count("[SCARLET pre-flight]") == 1
    assert all(f"self-assessment fold {k} of 5" in out for k in range(1, 6)) and "learn report:" in out
    assert err.count("Bootstrap:   0%") == 1   # learn's own bar only (each bar starts at 0 % once), not the folds'


def test_an_unknown_level_stops(pair, tmp_path):
    with pytest.raises(ValueError, match="verbosity"):
        api.learn(str(pair / "ref.h5ad"), "cell_type", str(tmp_path / "m.scarlet"), verbosity="loud", **KW)
