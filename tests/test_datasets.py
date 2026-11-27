"""scarlet.datasets.fetch: downloads the tutorial files, checks each SHA-256, keeps a file already present with the
right checksum, deletes a download whose checksum differs and stops; python -m scarlet.datasets does the same.
Served here from a local folder (file://) with stand-in files, so no network is needed."""
import hashlib
import os
import sys

import pytest

import scarlet
from scarlet import datasets


@pytest.fixture
def served(tmp_path, monkeypatch):
    src = tmp_path / "served"
    src.mkdir()
    files = {}
    for name, body in (("a.h5ad", b"first file\n" * 1000), ("b.h5ad", b"second file\n")):
        (src / name).write_bytes(body)
        files[name] = (hashlib.sha256(body).hexdigest(), len(body))
    monkeypatch.setattr(datasets, "FILES", files)
    return src.as_uri() + "/"


def test_a_download_shows_its_progress(served, tmp_path, capsys):
    """A long download must not look stuck: the file's line shows the MB downloaded and
    ends at 100 %, then the checksum line."""
    datasets.fetch(["a.h5ad"], dest=str(tmp_path / "p"), base_url=served)
    out = capsys.readouterr().out
    assert "a.h5ad: 0.0 of 0.0 MB (100%)" in out and out.index("(100%)") < out.index("a.h5ad: checksum ok")


def test_the_package_exposes_it():
    assert scarlet.datasets is datasets and len(datasets.FILES) == 7
    assert all(len(sha) == 64 and size > 0 for sha, size in datasets.FILES.values())


def test_fetch_downloads_checks_and_reuses(served, tmp_path, capsys):
    dest = tmp_path / "data"
    got = datasets.fetch(dest=str(dest), base_url=served)
    assert got == {n: str(dest / n) for n in ("a.h5ad", "b.h5ad")}
    assert (dest / "a.h5ad").read_bytes() == b"first file\n" * 1000
    assert "checksum ok" in capsys.readouterr().out
    got = datasets.fetch(["b.h5ad"], dest=str(dest), base_url=served)
    assert got == {"b.h5ad": str(dest / "b.h5ad")} and "present, checksum ok" in capsys.readouterr().out


def test_the_environment_variable_sets_the_source(served, tmp_path, monkeypatch):
    monkeypatch.setenv("SCARLET_EXAMPLE_DATA_URL", served)
    assert set(datasets.fetch(dest=str(tmp_path / "d"))) == {"a.h5ad", "b.h5ad"}


def test_a_wrong_checksum_stops_and_deletes(served, tmp_path, monkeypatch):
    monkeypatch.setitem(datasets.FILES, "b.h5ad", ("0" * 64, 12))
    dest = tmp_path / "data"
    with pytest.raises(RuntimeError, match="differs"):
        datasets.fetch(["b.h5ad"], dest=str(dest), base_url=served)
    assert list(dest.iterdir()) == []


def test_an_unknown_name_stops(served, tmp_path):
    with pytest.raises(ValueError, match="unknown file"):
        datasets.fetch(["c.h5ad"], dest=str(tmp_path), base_url=served)


def test_the_command_line(served, tmp_path):
    datasets.main(["--dest", str(tmp_path / "cli"), "--base-url", served, "a.h5ad"])
    assert [p.name for p in (tmp_path / "cli").iterdir()] == ["a.h5ad"]
    with pytest.raises(SystemExit, match="unknown file"):
        datasets.main(["--dest", str(tmp_path / "cli"), "--base-url", served, "c.h5ad"])


def _fake_gh(tmp_path, monkeypatch, src, releases):
    """A stand-in for the GitHub CLI on PATH: `gh api .../releases` lists `releases`; `gh release download TAG -p NAME
    -D DIR` copies NAME from src/TAG/ (exit 1 when that tag does not have it)."""
    bin_ = tmp_path / "bin"
    bin_.mkdir()
    gh = bin_ / "gh"
    gh.write_text(f"""#!{sys.executable}
import json, os, shutil, sys
a = sys.argv[1:]
if a[0] == "api":
    print(json.dumps([{{"tag_name": t}} for t in {releases!r}]))
    sys.exit(0)
tag, name, dest = a[2], a[a.index("-p") + 1], a[a.index("-D") + 1]
f = os.path.join({str(src)!r}, tag, name)
if not os.path.exists(f):
    print("release not found", file=sys.stderr)
    sys.exit(1)
shutil.copy(f, os.path.join(dest, name))
""")
    gh.chmod(0o755)
    monkeypatch.setenv("PATH", f"{bin_}{os.pathsep}{os.environ['PATH']}")


def test_a_private_repository_falls_back_to_the_github_cli(tmp_path, monkeypatch, capsys):
    """HTTPS fails (as for a private repository's assets): the newest v0.6.3 pre-release is fetched with gh."""
    body = b"tutorial data\n"
    monkeypatch.setattr(datasets, "FILES", {"a.h5ad": (hashlib.sha256(body).hexdigest(), len(body))})
    monkeypatch.setattr(datasets, "BASE_URL", (tmp_path / "nowhere").as_uri() + "/")
    monkeypatch.delenv("SCARLET_EXAMPLE_DATA_URL", raising=False)
    src = tmp_path / "releases"
    (src / "v0.6.3rc2").mkdir(parents=True)
    (src / "v0.6.3rc2" / "a.h5ad").write_bytes(body)
    _fake_gh(tmp_path, monkeypatch, src, ["v0.6.3rc2", "v0.6.3rc1", "v0.6.2"])
    got = datasets.fetch(dest=str(tmp_path / "data"))
    assert open(got["a.h5ad"], "rb").read() == body
    out = capsys.readouterr().out
    assert "downloaded with gh from release v0.6.3rc2" in out
    assert "\ra.h5ad: 0.0 of 0.0 MB (100%)\n" in out   # the gh route's line ends at 100 %, as the HTTPS one
    assert sorted(os.listdir(tmp_path / "data")) == ["a.h5ad"]   # no leftovers


def test_neither_https_nor_gh_gives_a_clear_message(tmp_path, monkeypatch):
    monkeypatch.setattr(datasets, "FILES", {"a.h5ad": ("0" * 64, 1)})
    monkeypatch.setattr(datasets, "BASE_URL", (tmp_path / "nowhere").as_uri() + "/")
    monkeypatch.delenv("SCARLET_EXAMPLE_DATA_URL", raising=False)
    monkeypatch.setenv("PATH", str(tmp_path))   # no gh
    with pytest.raises(RuntimeError, match=r"gh auth login.*by hand"):
        datasets.fetch(dest=str(tmp_path / "data"))
    with pytest.raises(RuntimeError, match="could not download"):   # another source: no gh attempt
        datasets.fetch(dest=str(tmp_path / "data"), base_url=(tmp_path / "elsewhere").as_uri() + "/")


def test_the_notebook_script_is_scarlet_datasets():
    import importlib.util
    path = os.path.join(os.path.dirname(__file__), os.pardir, "notebooks", "fetch_example_data.py")
    spec = importlib.util.spec_from_file_location("fetch_example_data", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    assert mod.fetch is datasets.fetch and mod.FILES is datasets.FILES and mod.BASE_URL == datasets.BASE_URL
