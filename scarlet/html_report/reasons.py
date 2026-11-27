"""Refusal reasons per cell, as the outputs name them (background, ambiguity, flip test, undetermined across draws).

Nothing here decides a call. Transfer records each cell's reason (first_pass.tsv / labels.tsv); for assign and the
self-assessment the reasons are read back from the stored scores, in the order SCARLET applies its steps
(scarlet/classifier.py classify):
  1. background: the background likelihood is at least the best identity's (ll_bg >= ll_top1) -> best and unique
     refused;
  2. ambiguity: best called, but the stored FDR-adjusted p of the top identity against its competitor is >= alpha
     -> unique refused;
  3. flip test (flip mode 'refuse'): the call depends on genes the candidate class never showed (stored scarlet_flip)
     -> best and unique refused.
One primary reason per refused cell; flags beside it (flip when the run's flip mode is 'flag'; low positive evidence;
not refined in transfer). `check_exact` is the build requirement: every refusal and every non-refusal explained, no
residual bucket.
"""
import numpy as np
import pandas as pd

U = "unassigned"
LABEL = "derived from SCARLET's rules on stored scores"
BACKGROUND = "background"
AMBIGUOUS = "ambiguity"
FLIP = "flip test"
UNDETERMINED = "undetermined across draws"
RECORDED = (BACKGROUND, AMBIGUOUS, FLIP, UNDETERMINED)   # transfer's recorded reasons (labels.tsv / first_pass.tsv)
CALLED = ""


def _bool(x):
    return pd.Series(x).astype(str).str.lower().isin(["true", "1"]).to_numpy()


def assign_reasons(best, unique, ll_top1, ll_bg, fdr, alpha, flip=None, flip_mode="off"):
    """Per cell (arrays of equal length): reason_best, reason_unique (CALLED where the call was made).
    best / unique: the stored calls ('unassigned' = refused); ll_top1, ll_bg, fdr: the stored scores; flip: the
    stored flip column (None when the flip test was off); flip_mode: 'refuse', 'flag' or 'off' (the run's)."""
    best, unique = np.asarray(best, str), np.asarray(unique, str)
    rb, ru = best == U, unique == U
    bg = np.asarray(ll_bg, float) >= np.asarray(ll_top1, float)   # NaN (no background) -> False
    fl = _bool(flip) if flip is not None else np.zeros(len(best), bool)
    amb = np.asarray(fdr, float) >= float(alpha)
    r_best = np.full(len(best), CALLED, dtype=object)
    r_best[rb & bg] = BACKGROUND
    if flip_mode == "refuse":
        r_best[rb & ~bg & fl] = FLIP
    r_uni = np.where(rb, r_best, CALLED).astype(object)
    r_uni[ru & ~rb & amb] = AMBIGUOUS
    return r_best, r_uni


def check_exact(best, unique, r_best, r_uni, *, where=""):
    """The build requirement: a refused cell has a reason, a called cell none. Raises AssertionError naming the first
    offending cells."""
    best, unique = np.asarray(best, str), np.asarray(unique, str)
    r_best, r_uni = np.asarray(r_best, object), np.asarray(r_uni, object)
    for name, call, r in (("best", best, r_best), ("unique", unique, r_uni)):
        refused = call == U
        bad = np.where(refused != (r != CALLED))[0]
        if len(bad):
            raise AssertionError(f"{where} {name}: {len(bad)} cells whose refusal is not explained or whose reason "
                                 f"does not match the call (first rows {bad[:5].tolist()})")
    return True


def counts(best, unique, r_best, r_uni, flip=None, pe_low=None, not_refined=None):
    """The refusal-reason counts (best / unique refused, background, flip, ambiguity, other), from the per-cell
    reasons, for the consistency check; there is no residual bucket, so the 'other' keys must come out 0."""
    best, unique = np.asarray(best, str), np.asarray(unique, str)
    rb, ru = best == U, unique == U
    r_best, r_uni = np.asarray(r_best, object), np.asarray(r_uni, object)
    bg = r_best == BACKGROUND
    fl = r_best == FLIP
    amb = (r_uni == AMBIGUOUS) & ~rb
    n = len(best)
    r = dict(best_refused=int(rb.sum()), best_background=int(bg.sum()), best_flip_refused=int(fl.sum()),
             best_other=int((rb & ~bg & ~fl).sum()), unique_refused=int(ru.sum()),
             unique_background=int((ru & bg).sum()), unique_flip_refused=int((ru & fl).sum()),
             unique_ambiguity=int(amb.sum()), unique_other=int((ru & ~bg & ~fl & ~amb).sum()),
             flip_flagged=int(_bool(flip).sum()) if flip is not None else 0,
             pe_low_flagged=int(_bool(pe_low).sum()) if pe_low is not None else 0)
    if not_refined is not None:
        r["not_refined_noted"] = int(np.asarray(not_refined, bool).sum())
    assert r["best_refused"] <= n
    return r
