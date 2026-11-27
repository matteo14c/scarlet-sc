"""Reading an input .h5ad with its count values held in the narrowest exact type.

anndata reads X's values in the type the file stores; many count files store integers as int64 (8 bytes per value).
read_h5ad(path) reads such a file with X's values as int32 (4 bytes) when every value fits, read in chunks into the
int32 array (no int64 copy of the whole matrix is held), and everything else exactly as anndata.read_h5ad does (its
own reading of each element; X built as anndata builds a sparse matrix). The values are the same integers; every step
that computes with them converts them to the same floating type as before (float64 sums and logarithms; scanpy's
float32), so the results are the same. Anything else (dense X, values stored as float or another integer type, a
value outside int32, raw or legacy layouts): anndata.read_h5ad as it is. Not used where the values are written back
(assign's annotated output keeps the input's type). Only anndata's public API and h5py are used, so the reading
follows anndata's releases."""
from __future__ import annotations

import numpy as np

CHUNK = 1 << 24
_I32 = np.iinfo(np.int32)
# the AnnData fields stored as top-level elements of an .h5ad file (X is read separately)
_FIELDS = ("obs", "var", "obsm", "varm", "obsp", "varp", "layers", "uns")


def _read_elem():
    """anndata's public element reader: anndata.io.read_elem (anndata >= 0.11), anndata.experimental.read_elem before."""
    try:
        from anndata.io import read_elem
    except ImportError:
        from anndata.experimental import read_elem
    return read_elem


def _int32_sparse(X, enc):
    """X's matrix (an h5py group of a csr / csc matrix with int64 values) with int32 values, or None when a value does
    not fit int32. Read with h5py; the matrix is scipy's sparse array or matrix class as anndata itself returns it
    (anndata.settings.use_sparse_array_on_read), its index arrays as stored."""
    import scipy.sparse as sp
    from anndata import settings
    ds = X["data"]
    n = int(ds.shape[0])
    data = np.empty(n, dtype=np.int32)
    for s in range(0, n, CHUNK):
        c = ds[s:s + CHUNK]
        if c.size and (c.min() < _I32.min or c.max() > _I32.max):
            return None
        data[s:s + len(c)] = c
    shape = tuple(int(x) for x in X.attrs["shape"])
    fmt = "csr" if enc == "csr_matrix" else "csc"
    array = bool(getattr(settings, "use_sparse_array_on_read", False))
    cls = {("csr", False): sp.csr_matrix, ("csc", False): sp.csc_matrix,
           ("csr", True): sp.csr_array, ("csc", True): sp.csc_array}[(fmt, array)]
    mtx = cls(shape, dtype=np.int32)
    mtx.data = data
    mtx.indices = X["indices"][...]
    mtx.indptr = X["indptr"][...]
    return mtx


def read_h5ad(path):
    """anndata.read_h5ad(path), with X's int64 count values held as int32 when they all fit (see the module). Uses
    only anndata's public API (read_elem, settings) and h5py."""
    import anndata as ad
    import h5py
    with h5py.File(path, "r") as f:
        X = f.get("X")
        enc = X.attrs.get("encoding-type", "") if isinstance(X, h5py.Group) else ""
        enc = enc.decode() if isinstance(enc, bytes) else enc
        plain = not (enc in ("csr_matrix", "csc_matrix") and X["data"].dtype == np.int64
                     and isinstance(f.get("obs"), h5py.Group) and "encoding-type" in f["obs"].attrs
                     and "raw" not in f
                     and not any(k.startswith("raw.") for k in f)
                     and set(f) - {"X"} <= set(_FIELDS))
        if not plain:
            X32 = _int32_sparse(X, enc)
            if X32 is not None:
                read_elem = _read_elem()
                return ad.AnnData(X=X32, **{k: read_elem(f[k]) for k in _FIELDS if k in f})
    return ad.read_h5ad(path)
