"""Release 0.6.3: the blocklist option is auto (the default), human, mouse, none or a file path. auto picks the bundled
human or mouse list from the gene names and records how; a path that does not exist is an error (before 0.6.3 it
meant no blocklist, with a warning); the list's version is recorded but not printed."""
import argparse
import zipfile

import anndata as ad
import numpy as np
import pandas as pd
import pytest

from scarlet import api, settings
from scarlet.blocklists import (BLOCKLIST_CHOICES, KNOWN_BLOCKLISTS, bundled_path, check_blocklist_option,
                                choose_blocklist, file_md5)
from synthetic_dense import _data

KW = dict(n_iterations=6, n_cells=20)
HUMAN = [f"GENE{i}" for i in range(30)]
MOUSE = [f"Gene{i}" for i in range(30)]


def _members(p):
    with zipfile.ZipFile(p) as z:
        return {n: z.read(n) for n in z.namelist()
                if n.endswith(".tsv") and any(k in n for k in ("marker", "programme", "universe", "background"))}


@pytest.fixture(scope="module")
def data(tmp_path_factory):
    d = tmp_path_factory.mktemp("r063bl")
    A = _data(d / "human.h5ad")
    M = A.copy()
    M.var_names = [f"Gm{j}" for j in range(M.n_vars)]        # capitalised: mouse-style
    M.write_h5ad(d / "mouse.h5ad")
    return d


def test_bundled_lists_and_their_versions():
    assert BLOCKLIST_CHOICES == ("auto", "human", "mouse", "none")
    versions = {k: KNOWN_BLOCKLISTS[file_md5(bundled_path(k))][0] for k in ("human", "mouse", "none")}
    assert versions == {"human": "U4", "mouse": "M4", "none": "E0"}
    assert open(bundled_path("none")).read() == "ENSG_ID\tSymbol\tCategory\n"   # header only: nothing removed


@pytest.mark.parametrize("names, choice", [(HUMAN, "human"), (MOUSE, "mouse"),
                                           ([f"ENSG{i:011d}" for i in range(30)], "human"),
                                           ([f"ENSMUSG{i:011d}" for i in range(30)], "mouse")])
def test_auto_picks_the_species_list(names, choice):
    c = choose_blocklist("auto", names)
    assert str(c) == bundled_path(choice)
    assert c.record["option"] == "auto" and c.record["choice"] == choice
    assert c.record["basis"].startswith("auto: 100.0% of 30 gene names are " + choice)
    assert choose_blocklist(None, names) == c                     # None is auto


@pytest.mark.parametrize("names", [HUMAN[:15] + MOUSE[:15], HUMAN[:10]])
def test_auto_refuses_when_the_species_is_unclear(names):
    with pytest.raises(ValueError, match="cannot tell whether the data are human or mouse"):
        choose_blocklist("auto", names)


def test_explicit_choices_and_files(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)                                   # no copy of a bundled file in the working folder
    for k in ("human", "mouse", "none"):
        c = choose_blocklist(k)
        assert str(c) == bundled_path(k) and c.record == dict(option=k, choice=k, basis=c.record["basis"])
    f = tmp_path / "mine.tsv"
    f.write_text("Symbol\nGENE1\n")
    c = choose_blocklist(str(f))
    assert str(c) == str(f) and c.record["choice"] == "file"     # an existing file is used as given
    for bad in (str(tmp_path / "missing.tsv"), "scarlet_blocklist_mouse.tsv", "U4"):
        with pytest.raises(FileNotFoundError, match="not found"):
            choose_blocklist(bad)
        with pytest.raises(FileNotFoundError):
            check_blocklist_option(bad)


def test_settings_record_and_preflight_name_no_version():
    r = settings.resolve(gene_names=pd.Index(HUMAN))
    b = r["blocklist"]
    assert (b["option"], b["choice"], b["given"], b["version"]) == ("auto", "human", False, "U4")
    text = settings.preflight_text("learn", r)
    assert "auto: 100.0% of 30 gene names are human-style" in text
    for code in ("U4", "M4", "E0", "md5"):
        assert code not in text


def test_learn_auto_equals_the_explicit_list(data):
    for sp in ("human", "mouse"):
        ma = api.learn(str(data / f"{sp}.h5ad"), "label", str(data / f"{sp}_auto.scarlet"), **KW)
        api.learn(str(data / f"{sp}.h5ad"), "label", str(data / f"{sp}_given.scarlet"), blocklist=sp, **KW)
        b = ma.settings["resolved"]["blocklist"]
        assert (b["option"], b["choice"]) == ("auto", sp)
        assert _members(data / f"{sp}_auto.scarlet") == _members(data / f"{sp}_given.scarlet")


def test_a_missing_file_stops_before_any_work(data):
    out = data / "never.scarlet"
    with pytest.raises(FileNotFoundError, match="blocklist file"):
        api.learn(str(data / "human.h5ad"), "label", str(out), blocklist=str(data / "no_such_list.tsv"), **KW)
    assert not out.exists()


def test_marker_finder_refuses_a_missing_file(tmp_path):
    from scarlet import marker_finder
    A = ad.AnnData(np.ones((3, 3), dtype=np.float32), var=pd.DataFrame(index=["A1", "B1", "C1"]))
    with pytest.raises(FileNotFoundError):
        marker_finder.apply_universal_blocklist(A.copy(), argparse.Namespace(blocklist=str(tmp_path / "x.tsv")))
