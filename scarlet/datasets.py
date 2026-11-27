"""The tutorial datasets: download them and check each file's SHA-256.

The files are the release assets of SCARLET v0.6.3: the inputs of SCARLET's formal benchmark (Abdelaal et al. 2019,
Zenodo 3357167, CC BY 4.0) with each benchmark task's own cell selection, compressed; their content is proven identical
to the benchmark's (MANIFEST.tsv in the release assets). Baron human pancreas (notebook 01), Segerstolpe human
pancreas (02), Allen mouse brain VISp and ALM and PBMC pbmc1 Smart-seq2 and 10x v3 (03). Sources and licences:
``ATTRIBUTION.md`` in the release assets.

In Python::

    from scarlet import datasets
    paths = datasets.fetch(["baron_human.h5ad"])          # {name: path}, in ./data

From the command line: ``python -m scarlet.datasets [--dest DIR] [--base-url URL] [FILE ...]``.

How the files are fetched: by HTTPS from the release (standard library only; no login). While the repository is
private, its release assets need a login: when HTTPS fails, the GitHub CLI (``gh``, logged in with ``gh auth login``)
downloads them from the v0.6.3 release, or from the newest v0.6.3 pre-release while there is no release yet. When
neither works, the function stops with a message saying what to do.
"""
from __future__ import annotations

import hashlib
import json
import os
import shutil
import subprocess
import sys
import tempfile
import urllib.error
import urllib.request

REPOSITORY = "matteo14c/scarlet-sc"
TAG = "v0.6.3"
BASE_URL = f"https://github.com/{REPOSITORY}/releases/download/{TAG}/"

# name: (SHA-256, size in bytes)
FILES = {
    "baron_human.h5ad": ("922416d0d8a78cc11142147cee63e7467256d2893921dcab0b06dfa64069fa8c", 37096717),
    "segerstolpe.h5ad": ("4ee0c405d2d8518f0c38bacfc5a7020261d48d9b236e9b99fbd2ed61173c3a2b", 38139993),
    "pbmc1_smartseq2.h5ad": ("64d20248c26f7444bec798a81bed94f61f0592e528aa8902493d673123ac15b9", 2757292),
    "pbmc1_10xv3_query.h5ad": ("6b389343549ecb0643ad29bcfb5f37eb6512365e6599b44aa98ce930d929300c", 11022072),
    "pbmc1_10xv3.h5ad": ("b527f436b5bad6d457040ec99aaad23071f61fc9d6170e0f18954ac26cc34a83", 12110345),
    "amb34_visp.h5ad": ("b8d9f00d73e3ffbb3108ac601987dc06d3af612ec4b1aa9364274fe8f0963169", 399930484),
    "amb34_alm.h5ad": ("7d303787875e67556e323f5ec416903b0f7a862b6fd9ec893366b7553c351b89", 277253460),
}


def _sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def _gh_tags(gh):
    """The release itself, then its pre-releases (newest first), as the GitHub CLI lists them."""
    r = subprocess.run([gh, "api", f"repos/{REPOSITORY}/releases?per_page=100"], capture_output=True, text=True)
    if r.returncode != 0:
        return [TAG]
    tags = [x.get("tag_name", "") for x in json.loads(r.stdout or "[]")]
    return [TAG] + [t for t in tags if t.startswith(TAG) and t != TAG]


class _Progress:
    """One line per file, rewritten in place while it downloads: 'name: 12.3 of 399.9 MB (3%)', at most twice a second;
    end() finishes the line."""

    def __init__(self, name, size):
        self.name, self.size, self.last, self.shown = name, size, 0.0, False

    def update(self, done):
        import time
        now = time.monotonic()
        if now - self.last < 0.5 and done < self.size:
            return
        self.last = now
        pct = f" ({100 * done / self.size:.0f}%)" if self.size else ""
        print(f"\r{self.name}: {done / 1e6:.1f} of {self.size / 1e6:.1f} MB{pct}", end="", flush=True)
        self.shown = True

    def end(self):
        if self.shown:
            print()
            self.shown = False


def _gh_download(name, part, tags_cache, progress=None):
    """name into part with the GitHub CLI; returns (the tag it came from, None) or (None, the reason)."""
    gh = shutil.which("gh")
    if gh is None:
        return None, "the GitHub CLI (gh) is not installed"
    if "tags" not in tags_cache:
        tags_cache["tags"] = _gh_tags(gh)
    err = ""
    for tag in tags_cache["tags"]:
        with tempfile.TemporaryDirectory(dir=os.path.dirname(part) or ".") as tmp:
            cmd = [gh, "release", "download", tag, "-R", REPOSITORY, "-p", name, "-D", tmp]
            with tempfile.TemporaryFile("w+") as out, tempfile.TemporaryFile("w+") as err:
                proc = subprocess.Popen(cmd, stdout=out, stderr=err, text=True)
                while True:   # the file's size in gh's folder, as it grows
                    try:
                        proc.wait(timeout=0.5)
                        break
                    except subprocess.TimeoutExpired:
                        if progress is not None:
                            progress.update(sum(os.path.getsize(os.path.join(tmp, f)) for f in os.listdir(tmp)))
                if progress is not None:
                    if proc.returncode == 0 and os.path.exists(os.path.join(tmp, name)):   # the last line: 100 %
                        progress.update(os.path.getsize(os.path.join(tmp, name)))
                    progress.end()
                out.seek(0)
                err.seek(0)
                r = subprocess.CompletedProcess(cmd, proc.returncode, out.read(), err.read())
            if r.returncode == 0 and os.path.exists(os.path.join(tmp, name)):
                os.replace(os.path.join(tmp, name), part)
                return tag, None
            lines = (r.stderr or r.stdout).strip().splitlines()
            err = lines[-1] if lines else ""
    return None, (f"gh release download failed: {err}" if err else "no release holds this file")


def fetch(files=None, dest="data", base_url=None) -> dict:
    """Download tutorial datasets and check their SHA-256.

    A file already in ``dest`` with the right checksum is not downloaded again; a download whose checksum differs is
    deleted and the function stops. Each file comes by HTTPS; when that fails and ``base_url`` is the release's, the
    GitHub CLI (``gh``, logged in) fetches it, which a private repository needs.

    Args:
        files (list[str] | None): The files to fetch, by name (see ``FILES``); ``None`` fetches all of them.
        dest (str): The folder to download into (created when missing).
        base_url (str | None): Where the files are; ``None`` uses the environment variable
            ``SCARLET_EXAMPLE_DATA_URL`` or else the v0.6.3 release assets. Another folder with the same files works
            too, e.g. a local copy (``file:///path/to/folder/``).

    Returns:
        dict: ``{name: path}`` of the checked files.

    Raises:
        ValueError: A name is not one of the tutorial files.
        RuntimeError: A downloaded file's SHA-256 differs from the expected one (the file is deleted); or a file
            could be downloaded neither by HTTPS nor with the GitHub CLI.
    """
    base_url = base_url or os.environ.get("SCARLET_EXAMPLE_DATA_URL") or BASE_URL
    tags_cache = {}
    os.makedirs(dest, exist_ok=True)
    out = {}
    for name in files or FILES:
        if name not in FILES:
            raise ValueError(f"unknown file {name!r}; known: {', '.join(FILES)}")
        want, size = FILES[name]
        path = os.path.join(dest, name)
        if os.path.exists(path) and _sha256(path) == want:
            print(f"{name}: present, checksum ok")
        else:
            url = base_url.rstrip("/") + "/" + name
            print(f"{name}: downloading {url} ({size / 1e6:.1f} MB)")
            try:
                bar = _Progress(name, size)
                try:
                    urllib.request.urlretrieve(url, path + ".part",
                                               reporthook=lambda n, bs, _t: bar.update(min(n * bs, size)))
                finally:
                    bar.end()
            except (urllib.error.URLError, OSError) as e:
                if os.path.exists(path + ".part"):
                    os.remove(path + ".part")
                if base_url.rstrip("/") != BASE_URL.rstrip("/"):
                    raise RuntimeError(f"{name}: could not download {url} ({e})") from None
                print(f"{name}: HTTPS failed ({e}); trying the GitHub CLI (a private repository needs a login)")
                tag, why = _gh_download(name, path + ".part", tags_cache, _Progress(name, size))
                if tag is None:
                    raise RuntimeError(
                        f"{name}: could not download it by HTTPS ({e}) nor with the GitHub CLI ({why}). While the "
                        f"repository is private: install the GitHub CLI and log in (gh auth login), then run this "
                        f"again; or download the release assets by hand into {dest}") from None
                print(f"{name}: downloaded with gh from release {tag}")
            got = _sha256(path + ".part")
            if got != want:
                os.remove(path + ".part")
                raise RuntimeError(f"{name}: SHA-256 {got} differs from {want}; the download was deleted")
            os.replace(path + ".part", path)
            print(f"{name}: checksum ok")
        out[name] = path
    return out


def main(argv=None):
    import argparse
    ap = argparse.ArgumentParser(prog="python -m scarlet.datasets", description=__doc__.splitlines()[0])
    ap.add_argument("--base-url", default=None, help=f"default: $SCARLET_EXAMPLE_DATA_URL or {BASE_URL}")
    ap.add_argument("--dest", default="data")
    ap.add_argument("files", nargs="*", help=f"default: all ({', '.join(FILES)})")
    a = ap.parse_args(argv)
    try:
        fetch(a.files or None, a.dest, a.base_url)
    except (RuntimeError, ValueError) as e:
        sys.exit(str(e))


if __name__ == "__main__":
    main()
