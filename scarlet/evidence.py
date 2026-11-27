"""The positive-evidence flag (the ratio rule; a label only).

For a call of cell c to identity k (the forced call: the highest-likelihood identity, whether or not R0 keeps it):
  m_c,k = k's markers (its programme, in the model) that c detects;
  E_k   = the sum of k's detection rates over those markers (the model's universe; the markers the data carry): the
          number of k's markers a typical k cell is expected to detect;
  the call is flagged 'low positive evidence' when m_c,k < floor(E_k / 5).
The typical value reported beside it is k's in-sample median of markers detected (over k's training cells, with the
model's own markers), stored by learn in the model member positive_evidence.json; no extra learning. A label only: the
call and its refusal are unchanged. assign writes obs scarlet_pe_detected, scarlet_pe_markers, scarlet_pe_expected,
scarlet_pe_required, scarlet_pe_typical, scarlet_pe_low and scarlet_pe_note.
"""
from __future__ import annotations

import json
import os

import numpy as np
import pandas as pd

RATIO = 5
MEMBER = "positive_evidence.json"
NOTE = ("low positive evidence for {k}: {m} of its {n} markers detected, fewer than floor({e:.1f} / 5) = {req} "
        "(E_k = {e:.1f} expected; a typical {k} training cell detects {t})")


def _marker_matrix(X, var, lists, ids):
    """cells x identities: the identity's markers each cell detects (X: detection CSR, var: gene names)."""
    pos = {}
    for i, g in enumerate(var):
        pos.setdefault(g, i)
    G = sorted({g for v in lists.values() for g in v})
    gi = {g: i for i, g in enumerate(G)}
    Ind = np.zeros((len(G), len(ids)))
    for j, k in enumerate(ids):
        for g in lists[k]:
            Ind[gi[g], j] = 1
    return np.asarray(X[:, [pos[g] for g in G]] @ Ind) if G else np.zeros((X.shape[0], len(ids)))


def _marker_counts(X, var, lists, ids):
    """_marker_matrix(_detection(A), var, lists, ids), reading only the marker columns of X (detection then column
    subset, or column subset then detection: the same 0/1 entries; their sums per cell are exact integers)."""
    from scipy import sparse
    pos = {}
    for i, g in enumerate(var):
        pos.setdefault(g, i)
    G = sorted({g for v in lists.values() for g in v})
    if not G:
        return np.zeros((X.shape[0], len(ids)))
    gi = {g: i for i, g in enumerate(G)}
    Ind = np.zeros((len(G), len(ids)))
    for j, k in enumerate(ids):
        for g in lists[k]:
            Ind[gi[g], j] = 1
    Xs = sparse.csr_matrix(X)[:, [pos[g] for g in G]]
    return np.asarray(sparse.csr_matrix(Xs > 0, dtype=np.float64) @ Ind)


def _detection(A):
    from scipy import sparse
    return sparse.csr_matrix(sparse.csr_matrix(A.X) > 0, dtype=np.float64)


def reference(data, labels, model, gene_col="index"):
    """The model's positive-evidence record: per identity its markers, E_k, the required floor(E_k / 5) and its
    in-sample median of markers detected over its training cells (`data`: an AnnData, or the path of one; obs
    `labels`). No learning; `data` is not changed."""
    import anndata as ad
    A = ad.read_h5ad(data) if isinstance(data, (str, os.PathLike)) else data
    lab = A.obs[labels].astype(str).to_numpy()
    var = (pd.Index(A.var[gene_col].astype(str).to_numpy()) if gene_col != "index" else A.var_names).astype(str)
    U = model.universe
    U.index = U.index.astype(str)
    ids = list(model.identities)
    have = set(var)
    lists = {k: [g for g in model.programmes[k] if g in U.index and g in have] for k in ids}
    M = _marker_counts(A.X, var, lists, ids)
    out = {}
    for j, k in enumerate(ids):
        E = float(U.loc[lists[k], k].sum()) if lists[k] else 0.0
        own = lab == k
        out[k] = {"markers": len(lists[k]), "expected": E, "required": int(np.floor(E / RATIO)),
                  "typical_median": float(np.median(M[own, j])) if own.any() else None, "n_cells": int(own.sum())}
    return {"schema": "scarlet.positive_evidence/2", "decision": "ratio rule: low when markers detected < floor(expected / ratio); label only", "ratio": RATIO,
            "identities": out}


def of_model(model):
    b = model.files.get(MEMBER)
    return None if b is None else json.loads(b)


def flags(A, model, forced, gene_col="index"):
    """Per cell the positive evidence of its forced call (`forced`: identity names, cells in A's order); A holds the
    detection the call was scored on. Returns a DataFrame (cells x the obs columns)."""
    rec = of_model(model) or {"identities": {}}
    var = pd.Index(A.var_names if gene_col == "index" else A.var[gene_col]).astype(str)
    U = model.universe
    U.index = U.index.astype(str)
    ids = list(model.identities)
    have = set(var)
    lists = {k: [g for g in model.programmes[k] if g in U.index and g in have] for k in ids}
    M = _marker_counts(A.X, var, lists, ids)
    col = {k: j for j, k in enumerate(ids)}
    forced = np.asarray(forced).astype(str)
    jj = np.array([col.get(k, -1) for k in forced])
    ok = jj >= 0
    m = np.where(ok, M[np.arange(len(forced)), np.maximum(jj, 0)], np.nan)
    n = np.array([len(lists.get(k, [])) for k in forced], dtype=float)
    Ek = {k: float(U.loc[lists[k], k].sum()) if lists[k] else np.nan for k in ids}
    E = np.array([Ek.get(k, np.nan) for k in forced])
    req = np.floor(E / RATIO)
    typ = np.array([(rec["identities"].get(k) or {}).get("typical_median") for k in forced], dtype=float)
    low = np.isfinite(m) & np.isfinite(req) & (m < req)
    note = np.full(len(forced), "", dtype=object)   # formatted for the flagged calls only
    m0, n0, r0 = np.nan_to_num(m), np.nan_to_num(n), np.nan_to_num(req)
    for i in np.flatnonzero(low):
        note[i] = NOTE.format(k=forced[i], m=int(m0[i]), n=int(n0[i]), e=E[i], req=int(r0[i]),
                              t="?" if not np.isfinite(typ[i]) else f"{typ[i]:g}")
    return pd.DataFrame({"scarlet_pe_detected": m, "scarlet_pe_markers": n, "scarlet_pe_expected": E,
                         "scarlet_pe_required": req, "scarlet_pe_typical": typ, "scarlet_pe_low": low,
                         "scarlet_pe_note": note}, index=A.obs_names)
