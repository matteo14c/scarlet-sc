"""The learn report's self-assessment: cross-validated learn + assign on the labelled reference.

Five folds by default (sklearn StratifiedKFold by label, shuffle, seeded): each fold's model is learned with the
model's own recorded learn parameters (settings.json 'parameters') on the other folds' cells, and assigns the fold's
cells, which it never saw. Nothing here touches the model being reported; the fold models are temporary.
"""
import inspect
import json
import os
import tempfile
import time

import numpy as np
import pandas as pd

from scarlet.html_report import reasons as R

#: learn parameters that are not settings (inputs, outputs and resources)
NOT_SETTINGS = ("reference", "labels", "model", "ncpu", "verbosity")


def learn_params(model):
    """The model's recorded learn parameters, restricted to api.learn's settings."""
    from scarlet import api
    p = dict(model.settings["parameters"])
    names = set(inspect.signature(api.learn).parameters) - set(NOT_SETTINGS)
    return {k: v for k, v in p.items() if k in names}


def make_folds(obs_names, y, *, n_folds=5, seed=0):
    """fold (0-based) per cell: StratifiedKFold by label, shuffled, seeded"""
    from sklearn.model_selection import StratifiedKFold
    n = len(y)
    f = np.full(n, -1)
    skf = StratifiedKFold(n_splits=n_folds, shuffle=True, random_state=seed)
    for k, (_, te) in enumerate(skf.split(np.zeros(n), y)):
        f[te] = k
    return f, dict(source="StratifiedKFold by label, shuffle", n_folds=n_folds, seed=seed)


def _fold_data(A):
    """The data the folds are taken from. Memory: the folds' learn and assign read detection only (the binary core:
    X > 0), so they take A's detection pattern (A's index arrays, the shared all-True values; rows copied per fold:
    5 bytes per entry) instead of copies of its values; the obs and var as A's. A is not changed."""
    import anndata as ad
    import scipy.sparse as sp
    from scarlet.marker_finder import binary_core_enabled, positive_pattern
    if not binary_core_enabled() or not sp.issparse(A.X):
        return A
    return ad.AnnData(X=positive_pattern(sp.csr_matrix(A.X)), obs=A.obs, var=A.var)


def run(model, reference, *, n_folds=5, seed=0, alpha=0.01, ncpu=1, workdir=None, log=print):
    """Per-cell out-of-fold calls and stored scores; returns (cells DataFrame, record dict)."""
    import anndata as ad
    from scarlet import api
    from scarlet.model import Model, load_model
    m = model if isinstance(model, Model) else load_model(model)
    params = learn_params(m)
    labels = m.settings["parameters"]["labels"]
    from scarlet._read import read_h5ad   # int64 counts held as int32 (the same values; scarlet._read)
    A = read_h5ad(reference) if isinstance(reference, (str, os.PathLike)) else reference
    y = A.obs[labels].astype(str).to_numpy()
    fold, src = make_folds(A.obs_names, y, n_folds=n_folds, seed=seed)
    D = _fold_data(A)
    tmp = tempfile.mkdtemp(prefix="scarlet_selfassess_", dir=workdir)
    rows, times, flip_modes = [], [], set()
    from scarlet import _verbosity
    progress = _verbosity.bar(len(set(fold)), "self-assessment")
    try:
        for k in sorted(set(fold)):
            te = fold == k
            Rf, Qf = D[~te].copy(), D[te].copy()
            Qf.obs = pd.DataFrame(index=Qf.obs_names)
            t0 = time.time()
            api.learn(Rf, labels, f"{tmp}/fold{k}.scarlet", ncpu=ncpu, **params)
            t1 = time.time()
            lab = api.assign(f"{tmp}/fold{k}.scarlet", Qf, f"{tmp}/fold{k}.h5ad", ncpu=ncpu, alpha=alpha,
                             gene_col=params.get("gene_col", "index"))
            t2 = time.time()
            Z = ad.read_h5ad(f"{tmp}/fold{k}.h5ad", backed="r")
            obs = Z.obs.copy()
            flip_modes.add(str(Z.uns.get("scarlet_assign", {}).get("flip", "off")))
            Z.file.close()
            lab.index = lab.index.astype(str)
            obs.index = obs.index.astype(str)
            idx = Qf.obs_names.astype(str)
            P = pd.DataFrame({"cell": idx, "fold": k + 1, "truth": y[te],   # folds numbered 1 to n, as the log
                              "forced": obs["scarlet_top1"].astype(str).reindex(idx).to_numpy(),
                              "best": lab["scarlet_best"].astype(str).reindex(idx).to_numpy(),
                              "unique": lab["scarlet_unique"].astype(str).reindex(idx).to_numpy(),
                              "fdr": obs["scarlet_fdr"].reindex(idx).to_numpy(float),
                              "ll_top1": obs["scarlet_ll_top1"].reindex(idx).to_numpy(float),
                              "ll_bg": obs["scarlet_ll_bg"].reindex(idx).to_numpy(float),
                              "top2": (obs["scarlet_top2"].astype(str).reindex(idx).to_numpy() if "scarlet_top2" in obs
                                       else np.full(len(idx), "", dtype=object)),
                              "ll_top2": (obs["scarlet_ll_top2"].reindex(idx).to_numpy(float)
                                          if "scarlet_ll_top2" in obs else np.full(len(idx), np.nan))})
            for c, name in (("scarlet_flip", "flip"), ("scarlet_pe_low", "pe_low"),
                            ("scarlet_pe_detected", "pe_detected"), ("scarlet_pe_expected", "pe_expected"),
                            ("scarlet_pe_required", "pe_required")):
                if c in obs:
                    P[name] = obs[c].reindex(idx).to_numpy()
            rows.append(P)
            times.append(dict(fold=int(k) + 1, train_cells=int((~te).sum()), test_cells=int(te.sum()),
                              learn_s=round(t1 - t0, 1), assign_s=round(t2 - t1, 1)))
            os.remove(f"{tmp}/fold{k}.h5ad")
            progress.update(1)
            log(f"[report] self-assessment fold {k + 1} of {len(set(fold))}: learn {t1 - t0:.0f} s, assign "
                f"{t2 - t1:.0f} s")
    finally:
        progress.close()
        for f in os.listdir(tmp):
            os.remove(os.path.join(tmp, f))
        os.rmdir(tmp)
    C = (pd.concat(rows, ignore_index=True).set_index("cell")
         .reindex(pd.Index(A.obs_names.astype(str), name="cell")).reset_index())
    assert C.fold.notna().all(), "self-assessment: a reference cell was not assigned in any fold"
    if len(flip_modes) != 1:
        raise RuntimeError(f"self-assessment: the folds were assigned with different flip modes {flip_modes}")
    flip_mode = flip_modes.pop()
    rb, ru = R.assign_reasons(C.best, C.unique, C.ll_top1, C.ll_bg, C.fdr, alpha,
                              flip=C["flip"] if "flip" in C else None, flip_mode=flip_mode)
    R.check_exact(C.best, C.unique, rb, ru, where="self-assessment")
    C["reason_best"], C["reason_unique"] = rb, ru
    # certainty and boundary: the margin of the top identity over its competitor (R0's statistic) and the
    # runner-up identity (the second-ranked where the best call is the top identity, else the top identity)
    C["margin"] = C.ll_top1 - np.fmax(C.ll_top2, C.ll_bg)
    C["runner_up"] = np.where(C.best.astype(str) == C.forced.astype(str), C.top2.astype(str),
                              C.forced.astype(str))
    C.loc[C.runner_up.isin(["", "nan", "None"]), "runner_up"] = ""
    rec = dict(folds=src, alpha=float(alpha), flip_mode=flip_mode,
               learn_parameters="the model's recorded learn parameters (settings.json 'parameters')",
               times=times, cells=int(len(C)))
    return C, rec
