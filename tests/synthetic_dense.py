"""Synthetic dense data for the tests (about 5,200 detected genes per cell, three labels)."""
import zipfile

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse



def _data(path):
    rng = np.random.default_rng(1)
    n, g = 120, 6400
    labels = np.repeat(["A", "B", "C"], n // 3)
    X = rng.poisson(1.6, size=(n, g)).astype(np.float32)  # ~5,200 detected genes per cell: many cells above 5,000
    for k, lab in enumerate(["A", "B", "C"]):
        X[labels == lab, k * 300:(k + 1) * 300] += rng.poisson(6, size=((labels == lab).sum(), 300))
    A = ad.AnnData(X=sparse.csr_matrix(X), obs=pd.DataFrame({"label": labels}, index=[f"c{i}" for i in range(n)]),
                   var=pd.DataFrame(index=[f"G{j}" for j in range(g)]))
    A.write_h5ad(path)
    return A


def _members(p):
    with zipfile.ZipFile(p) as z:
        return {n: z.read(n) for n in z.namelist()
                if n.endswith(".tsv") and any(k in n for k in ("marker", "programme", "universe", "background"))}
