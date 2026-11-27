"""Assign's scoring: the Jeffreys floor and the flip test.

The floor (j, the scoring of every call): each identity's detection rates clamped to [f_k, 1 - f_k], f_k = 0.5 / (n_eff,k + 1), n_eff,k =
the distinct training cells of k drawn across learn's bootstrap rounds (recorded by learn in the model's
marker_finder_params.json, `scoring_stats`); the background's to [f_bg, 1 - f_bg] (natural: 0.5 / (N + 1), N the
reference cells; class_balanced: the mean over identities of 0.5 / (n_k + 1)); then log p, log(1 - p) with p clipped to
[1e-5, 1 - 1e-5] as a numerical guard. A rate measured on n cells cannot honestly be closer to 0 or 1 than about
1 / (n + 1): the floor stops a gene that an identity never showed from vetoing it with log(1e-5).
Without the floor (eps; used only by the flip test): log(p + eps) with eps = 1e-5.

The flip test: every cell scored under eps and under j; a cell whose highest-likelihood (forced) identity differs
between the two depends on how never-seen genes are treated. Recorded per cell (obs scarlet_flip, scarlet_top1_eps,
scarlet_top1_j, scarlet_reason); refused (best unassigned, not unique) when the flip is applied (assign's default on every
model; an option in transfer, whose default is to flag).
"""
from __future__ import annotations

import numpy as np

GUARD = 1e-5
FLIP_REASON = "the call depends on genes the candidate class never showed"
FLIP_MODES = ("refuse", "flag", "off")   # off: the refinement inside transfer


def clamp(P, f):
    """Rates P (genes x identities, or genes) clamped to [f, 1 - f] (f per column, or a scalar), then the guard."""
    P = np.asarray(P, dtype=float)
    f = np.asarray(f, dtype=float)
    return np.clip(np.minimum(np.maximum(P, f), 1.0 - f), GUARD, 1.0 - GUARD)


def scoring_stats_of(model):
    return model.marker_finder_params.get("scoring_stats")


def jeffreys(model):
    """(identity -> f_k, f_bg) of a model whose learn recorded the scoring statistics, or None (a model without
    them)."""
    st = scoring_stats_of(model)
    if not st:
        return None
    ids = list(model.identities)
    I = st["identities"]
    missing = [k for k in ids if k not in I]
    if missing:
        raise ValueError(f"scoring_stats lack identities {missing[:5]}")
    fk = {k: 0.5 / (I[k]["n_eff"] + 1.0) for k in ids}
    bg = st["background"]
    if bg["mode"] == "natural":
        fbg = 0.5 / (bg["N"] + 1.0)
    else:
        fbg = float(np.mean([0.5 / (I[k]["n_k"] + 1.0) for k in ids]))
    return fk, fbg


def resolve_flip(flip):
    """The flip mode, checked: refuse (assign's default, every technology), flag (transfer's default), off."""
    if flip not in FLIP_MODES:
        raise ValueError(f"flip must be one of {FLIP_MODES} (got {flip!r})")
    return flip
