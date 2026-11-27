"""Input checks shared by the marker finder and the classifier.

SCARLET reads detection (value > 0) from counts. Negative values mean the matrix is not counts (e.g. batch-corrected
or scaled data); reading detection from it would silently turn corrected values into detection calls, so SCARLET
stops.
"""
import numpy as np
from scipy import sparse

NON_COUNT_MESSAGE = "non-count input, e.g. batch-corrected, is outside SCARLET's scope"


class NonCountInputError(ValueError):
    """The expression matrix holds negative values."""


def check_non_negative(X, where="the input matrix"):
    """Raise NonCountInputError if X (dense or sparse) has any negative value."""
    if sparse.issparse(X):
        values = X.data if hasattr(X, "data") else X.tocsr().data
    else:
        values = np.asarray(X)
    if values.size == 0:
        return
    vmin = values.min()
    if vmin < 0:
        n_neg = int(np.count_nonzero(values < 0))
        raise NonCountInputError(
            f"{where} has {n_neg} negative values (min {float(vmin):.4g}): {NON_COUNT_MESSAGE}. "
            "SCARLET needs raw or normalised counts (all values >= 0).")
