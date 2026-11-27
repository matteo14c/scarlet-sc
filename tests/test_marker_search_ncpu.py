"""The marker search gives the same bits whatever ncpu: the bootstrap's sums are added in iteration order."""
from __future__ import annotations

from types import SimpleNamespace

import numpy as np
from scipy import sparse

from scarlet import marker_finder as M


def _sums(X, groups, types, p_u, ncpu):
    args = SimpleNamespace(alpha=1.5, beta=1.0, gamma=6.0, eps=1e-4, n_cells=40, top_k=30, n_iterations=23,
                           active_genes_seed=42, bootstrap_order="keyed", ncpu=ncpu)
    out = [np.zeros((X.shape[1], len(types))), np.zeros((X.shape[1], len(types))),
           np.zeros((X.shape[1], len(types)), dtype=np.int64)]
    for scores, p_in, presence, _, _ in M.run_bootstrap_tasks(X, groups, types, p_u, args, False):
        out[0] += scores
        out[1] += p_in
        out[2] += presence
    return out


def test_bootstrap_sums_same_bits_whatever_ncpu():
    rng = np.random.default_rng(0)
    X = sparse.csr_matrix(rng.random((300, 400)) < rng.uniform(0.02, 0.6, 400))
    types = ["a", "b", "c", "d"]
    groups = {t: np.arange(i, 300, 4) for i, t in enumerate(types)}
    p_u = np.asarray(X.mean(axis=0)).ravel()
    ref = _sums(X, groups, types, p_u, 1)
    for n in (2, 4):
        for a, b in zip(ref, _sums(X, groups, types, p_u, n)):
            assert a.tobytes() == b.tobytes()
