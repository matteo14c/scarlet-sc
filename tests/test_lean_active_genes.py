"""Lean memory: the rewritten active-gene sampling of the marker search is bit-identical to 0.2.0's.

The reference implementation below is 0.2.0's code, copied verbatim (scarlet/marker_finder.py); only the function
name is changed."""
import numpy as np
import anndata as ad
import pandas as pd
import pytest
from scipy import sparse
from scipy.sparse import csr_matrix
from tqdm import tqdm

from scarlet.marker_finder import normalize_active_genes_by_sampling


def ref_normalize_active_genes_by_sampling(X, target_mean=1000, target_sd=250, seed=42):
    """
    Convert expression values to a binary matrix while sampling each cell toward
    a comparable number of active genes.

    Cells with an observed positive-gene count at or below target_mean are left
    untouched. Only cells above target_mean are sampled down toward a
    comparable active-gene count.
    """
    rng = np.random.default_rng(seed)
    X_csr = X.tocsr() if sparse.issparse(X) else sparse.csr_matrix(X)
    n_cells, n_genes = X_csr.shape
    target_counts = rng.normal(target_mean, target_sd, n_cells)
    target_counts = np.clip(target_counts, 1, n_genes).astype(int)

    out_indptr = [0]
    out_indices = []
    for i in tqdm(range(n_cells), desc="Active-gene normalization"):
        start = X_csr.indptr[i]
        end = X_csr.indptr[i + 1]
        row_indices = X_csr.indices[start:end]
        row_values = X_csr.data[start:end]

        positive_mask = row_values > 0
        row_indices = row_indices[positive_mask]
        row_values = row_values[positive_mask]
        n_positive = len(row_indices)
        if n_positive == 0:
            out_indptr.append(len(out_indices))
            continue

        if n_positive <= target_mean:
            selected = row_indices
        else:
            target_n = min(target_counts[i], n_positive)
            log_values = np.log1p(row_values)
            log_sum = log_values.sum()
            if log_sum > 0:
                probs = log_values / log_sum
            else:
                probs = None
            selected_positions = rng.choice(
                n_positive,
                size=target_n,
                replace=False,
                p=probs,
            )
            selected = row_indices[selected_positions]

        selected = np.sort(selected)
        out_indices.extend(selected.tolist())
        out_indptr.append(len(out_indices))

    out_data = np.ones(len(out_indices), dtype=bool)
    return sparse.csr_matrix(
        (out_data, np.asarray(out_indices, dtype=np.int32), np.asarray(out_indptr, dtype=np.int64)),
        shape=(n_cells, n_genes),
        dtype=bool,
    )



def _matrix(seed, n_cells=300, n_genes=400, density=0.3, unsorted=False, zeros=False):
    rng = np.random.default_rng(seed)
    X = sparse.random(n_cells, n_genes, density=density, format="csr", random_state=seed, dtype=np.float32)
    X.data = np.floor(X.data * 20).astype(np.float32)  # counts, some explicit zeros
    if not zeros:
        X.eliminate_zeros()
    X = X.tocsr()
    if unsorted:  # rows with unsorted column order (non-canonical CSR)
        for i in range(n_cells):
            s, e = X.indptr[i], X.indptr[i + 1]
            p = rng.permutation(e - s)
            X.indices[s:e] = X.indices[s:e][p]
            X.data[s:e] = X.data[s:e][p]
        X.has_sorted_indices = False
    return X


def _adata(X):
    obs = pd.DataFrame({"label": [f"c{i % 5}" for i in range(X.shape[0])]}, index=[f"cell{i}" for i in range(X.shape[0])])
    var = pd.DataFrame(index=[f"g{j}" for j in range(X.shape[1])])
    A = ad.AnnData(X=X, obs=obs, var=var)
    A.obsm["emb"] = np.arange(X.shape[0] * 2, dtype=float).reshape(-1, 2)
    A.uns["note"] = "x"
    return A


def _both(f_new, f_ref):
    """Run both; they must both succeed, or both raise the same error (0.2.0's behaviour kept)."""
    try:
        r_new, e_new = f_new(), None
    except Exception as exc:  # noqa: BLE001
        r_new, e_new = None, exc
    try:
        r_ref, e_ref = f_ref(), None
    except Exception as exc:  # noqa: BLE001
        r_ref, e_ref = None, exc
    assert (e_new is None) == (e_ref is None), (e_new, e_ref)
    if e_new is not None:
        assert type(e_new) is type(e_ref) and str(e_new) == str(e_ref)
        return None, None
    return r_new, r_ref


def _same(a, b):
    a, b = a.tocsr(), b.tocsr()
    assert a.shape == b.shape and a.dtype == b.dtype
    assert np.array_equal(a.indptr, b.indptr) and np.array_equal(a.indices, b.indices) and np.array_equal(a.data, b.data)


@pytest.mark.parametrize("mean,sd", [(60, 15), (200, 50), (5, 2)])
@pytest.mark.parametrize("unsorted,zeros", [(False, False), (True, False), (False, True)])
def test_normalize_active_genes_identical_to_main(mean, sd, unsorted, zeros):
    X = _matrix(13, unsorted=unsorted, zeros=zeros)
    new, ref = _both(lambda: normalize_active_genes_by_sampling(X.copy(), target_mean=mean, target_sd=sd, seed=9),
                     lambda: ref_normalize_active_genes_by_sampling(X.copy(), target_mean=mean, target_sd=sd, seed=9))
    if new is not None:
        _same(new, ref)
