"""Shared synthetic fixtures for the fast SCARLET unit tests (no real data)."""
from __future__ import annotations

import os
import tempfile

os.environ.setdefault("NUMBA_CACHE_DIR", tempfile.gettempdir())
os.environ.setdefault("MPLCONFIGDIR", tempfile.gettempdir())
os.environ.setdefault("HDF5_USE_FILE_LOCKING", "FALSE")

import numpy as np
import pandas as pd
import pytest
from scipy import sparse

# Genes from the universal blocklist, one module each at least (planted as markers of group A
# so that they would surface if they were not filtered). Categories are the U4 modules (the default list); the
# synthetic labels (A, B, C) leave every one of these modules ON under the module rule.
BLOCKLIST_GENES = {
    "FOS": "core_stress",
    "JUN": "core_stress",
    "MKI67": "cell_cycle",
    "TOP2A": "cell_cycle",
    "MT-CO1": "core_mitochondrial",
    "MT-ND1": "core_mitochondrial",
    "RPL3": "core_ribosomal",
    "RPS6": "core_ribosomal",
    "XIST": "core_sex",
}
GROUPS = ("A", "B", "C")
N_PER_GROUP = 80
N_PLAIN_GENES = 141
MARKERS_PER_GROUP = 10


def make_synthetic_adata(seed: int = 0):
    """240 cells x 150 genes, 3 groups with 10 planted markers each, integer counts."""
    import anndata as ad

    rng = np.random.default_rng(seed)
    plain = [f"G{i:03d}" for i in range(N_PLAIN_GENES)]
    genes = list(BLOCKLIST_GENES) + plain
    n_cells = N_PER_GROUP * len(GROUPS)
    labels = np.repeat(GROUPS, N_PER_GROUP)

    p = np.full((n_cells, len(genes)), 0.06)
    for g_idx, group in enumerate(GROUPS):
        rows = labels == group
        start = len(BLOCKLIST_GENES) + g_idx * MARKERS_PER_GROUP
        p[np.ix_(rows, range(start, start + MARKERS_PER_GROUP))] = 0.85
        # a few shared, weaker genes so the merge has something to look at
        p[np.ix_(rows, range(len(genes) - 10, len(genes)))] = 0.3 + 0.1 * g_idx
    p[np.ix_(labels == "A", range(len(BLOCKLIST_GENES)))] = 0.9
    detected = rng.random(p.shape) < p
    counts = detected * rng.integers(1, 6, size=p.shape)
    adata = ad.AnnData(
        X=sparse.csr_matrix(counts.astype(np.float32)),
        obs=pd.DataFrame({"cell_type": pd.Categorical(labels)}, index=[f"cell{i}" for i in range(n_cells)]),
        var=pd.DataFrame(index=genes),
    )
    return adata


@pytest.fixture(scope="session")
def synthetic_adata():
    return make_synthetic_adata()


@pytest.fixture(scope="session")
def synthetic_h5ad(tmp_path_factory, synthetic_adata):
    path = tmp_path_factory.mktemp("synthetic") / "synthetic.h5ad"
    synthetic_adata.write_h5ad(path)
    return path


