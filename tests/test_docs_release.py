"""The documentation and notebooks as the third review left them: benchmark links pinned to the release tag,
troubleshooting as message / cause / fix tables, the tutorials page's memory (what each number measures, on how many
CPUs) and laptop time, and notebook cells with ids."""
import json
import re
from pathlib import Path

import scarlet

ROOT = Path(__file__).resolve().parents[1]
DOCS = ROOT / "docs"


def test_benchmark_links_are_pinned_to_the_release_tag():
    tag = "v" + re.match(r"\d+\.\d+\.\d+", scarlet.__version__).group(0)
    text = "\n".join(f.read_text() for f in DOCS.rglob("*.md"))
    links = re.findall(r"github\.com/matteo14c/scarlet-sc/blob/([^/]+)/benchmark/", text)
    assert links and set(links) == {tag}


def test_troubleshooting_messages_are_tables():
    t = (DOCS / "troubleshooting.md").read_text()
    sections = re.split(r"^## ", t, flags=re.M)[1:]
    with_messages = [s for s in sections if not s.startswith("Memory and time")]
    assert len(with_messages) == 6
    for s in with_messages:
        assert "| Message | Cause | Fix |" in s, s.splitlines()[0]
    assert not re.search(r"^\*\*`", t, flags=re.M)   # no message as a bold paragraph title any more
    for row in (r for r in t.splitlines() if r.startswith("| ") and not r.startswith("| :")):
        assert row.count(" | ") >= 2 and row.rstrip().endswith("|"), row


def test_tutorials_state_what_the_memory_numbers_measure():
    t = (DOCS / "tutorials.md").read_text()
    assert "12 GB" not in t and "about 7 GB" not in t
    for k in ("measured with 2 CPUs, SCARLET alone", "Jupyter and\n    the browser come on top", "about 5 GB",
              "about 7.5 GB", "about 33 minutes", "Close and Shut Down Notebook", "activate it again"):
        assert k in t, k
    nb = json.loads((ROOT / "notebooks" / "03_transfer.ipynb").read_text())
    setup = "".join(nb["cells"][1]["source"])
    assert "12 GB" not in setup and "about 5 GB" in setup and "about 7.5 GB" in setup


def test_notebook_cells_have_ids():
    for f in sorted((ROOT / "notebooks").glob("*.ipynb")):
        nb = json.loads(f.read_text())
        assert (nb["nbformat"], nb["nbformat_minor"]) >= (4, 5), f.name
        ids = [c.get("id") for c in nb["cells"]]
        assert all(ids) and len(set(ids)) == len(ids), f.name
        assert all(re.fullmatch(r"[a-zA-Z0-9_-]{1,64}", i) for i in ids), f.name


def test_every_release_asset_has_its_source_and_licence():
    from scarlet import datasets
    t = (DOCS / "tutorials.md").read_text()
    sec = t[t.index("## Data sources and licences"):t.index("## Run them on your computer")]
    rows = [r for r in sec.splitlines() if r.startswith("| `")]
    names = [r.split("`")[1] for r in rows]
    assert sorted(names) == sorted(list(datasets.FILES) + ["MANIFEST.tsv", "SHA256SUMS", "ATTRIBUTION.md"])
    for r in rows:
        cells = [c.strip() for c in r.strip("|").split("|")]
        assert len(cells) == 4 and all(cells), r   # file, study, accession, licence
        assert "CC BY 4.0" in cells[3] or "BSD 3-Clause" in cells[3], r
