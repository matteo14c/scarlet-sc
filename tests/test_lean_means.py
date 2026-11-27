"""Lean memory: the boolean column mean and the shared-index positive pattern equal scipy's results bit for bit."""
import numpy as np
import pytest
from scipy import sparse

from scarlet.marker_finder import _bool_csr_column_mean, matrix_mean_vector, positive_pattern


def _counts(seed, n_rows, n_cols, density):
    X = sparse.random(n_rows, n_cols, density=density, format="csr", random_state=seed, dtype=np.float32)
    X.data = np.floor(X.data * 9) .astype(np.float32)
    return X


@pytest.mark.parametrize("n_rows,n_cols,density", [(1, 5, 1.0), (7, 3, 0.5), (1000, 50, 0.3), (100, 400, 0.9),
                                                    (65943, 40, 0.6), (329762, 12, 0.97), (12345, 30, 0.0)])
def test_bool_column_mean_bit_identical_to_scipy(n_rows, n_cols, density):
    B = _counts(n_rows + n_cols, n_rows, n_cols, density) > 0
    if B.nnz:
        B[0, :] = True  # one full row, so some columns reach every cell when density is high
    B = sparse.csr_matrix(B)
    ref = np.asarray(B.mean(axis=0)).ravel()
    new = matrix_mean_vector(B)
    assert new.dtype == ref.dtype and np.array_equal(new, ref)
    assert np.array_equal(_bool_csr_column_mean(B), ref)


def test_full_columns_long_rounding_chain():
    # every column fully populated: the running sum of 1/n over n cells (the case where rounding accumulates)
    for n in (3, 10, 4999, 117459, 329762):
        B = sparse.csr_matrix(np.ones((n, 2), dtype=bool))
        assert np.array_equal(matrix_mean_vector(B), np.asarray(B.mean(axis=0)).ravel()), n


def test_bool_with_explicit_false_entries():
    B = sparse.csr_matrix((np.array([True, False, True, True]), np.array([0, 1, 1, 2]), np.array([0, 2, 4])), shape=(2, 3))
    assert np.array_equal(matrix_mean_vector(B), np.asarray(B.mean(axis=0)).ravel())


def test_row_subsets_as_the_bootstrap_uses_them():
    B = (_counts(5, 5000, 80, 0.4) > 0).tocsr()
    rng = np.random.default_rng(1)
    for _ in range(20):
        idx = rng.integers(0, B.shape[0], size=100)
        S = B[idx, :]
        assert np.array_equal(matrix_mean_vector(S), np.asarray(S.mean(axis=0)).ravel())


def test_positive_pattern_equals_comparison():
    X = _counts(3, 400, 60, 0.3)
    X.eliminate_zeros()
    P = positive_pattern(X)
    R = X > 0
    assert (P != R).nnz == 0 and P.dtype == R.dtype and P.shape == R.shape
    assert np.shares_memory(P.indices, X.indices)  # no copy of the indices
    Z = _counts(4, 400, 60, 0.3)  # stored zeros: falls back to X > 0
    Z.data[::7] = 0
    P2 = positive_pattern(Z)
    assert (P2 != (Z > 0)).nnz == 0 and not np.shares_memory(P2.indices, Z.indices)
    assert np.array_equal(matrix_mean_vector(P), np.asarray(R.mean(axis=0)).ravel())


def test_pattern_then_column_subset_equals_subset_then_pattern():
    X = _counts(8, 500, 90, 0.4)
    X.eliminate_zeros()
    mask = np.random.default_rng(2).random(90) > 0.2
    A = positive_pattern(X)[:, mask]
    B = X[:, mask] > 0
    assert (A != B).nnz == 0 and A.dtype == B.dtype and A.shape == B.shape
    assert np.array_equal(matrix_mean_vector(sparse.csr_matrix(A)), np.asarray(B.mean(axis=0)).ravel())
