#!/usr/bin/env python3
"""Download SCARLET's example data and verify each file's SHA-256: scarlet.datasets, which this script calls.

    python fetch_example_data.py [--dest DIR] [--base-url URL] [FILE ...]

The same as ``python -m scarlet.datasets``. The files are the release assets of SCARLET v0.6.3
(https://github.com/matteo14c/scarlet-sc/releases/tag/v0.6.3), downloaded by HTTPS, or with the GitHub CLI (gh,
logged in) while the repository is private. --base-url (or the environment variable SCARLET_EXAMPLE_DATA_URL) points
to another folder with the same files, e.g. a local copy (file:///path/to/folder/). Files already present with the
right checksum are not downloaded again; a download whose checksum differs is deleted and the script stops.

In Python: ``from scarlet import datasets; paths = datasets.fetch(["baron_human.h5ad"])``.

Attribution: ATTRIBUTION.md (Abdelaal et al. 2019, Zenodo 3357167, CC BY 4.0, and the original studies).
"""
from scarlet.datasets import BASE_URL, FILES, fetch, main  # noqa: F401  (the notebooks and CI import these)

if __name__ == "__main__":
    main()
