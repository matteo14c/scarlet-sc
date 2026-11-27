"""The learn report: sections 0-5, computed into the bundle.

0 provenance; 1 depth (genes per cell per identity); 2 UMAPs (original labels; cross-validated best; cross-validated
unique, refusals in grey); 3 analytics per identity (correct, wrong and to which identity, refused by reason);
4 identity drill-down (the identity's markers x every identity of the universe, detection rates; the identity on the
UMAPs); 5 weak-class flags ranked weakest first with their evidence. Only values SCARLET computed are shown; the
under-50-cells marker is a reporting convention, labelled as such, not a SCARLET threshold.
"""
import datetime
import io
import os

import numpy as np
import pandas as pd

from scarlet import run_summary as _rs
from scarlet import settings as _settings
from scarlet.html_report import bundle as B
from scarlet.html_report import reasons as R

ENTRIES = ("forced", "best", "unique")
UNDER = 50


def cv_per_identity(C, ids):
    """Per reference label, the out-of-fold certainty: cells; % best correct; % with a unique call; %
    boundary (best call, no unique call); % unlabelled by reason (best, unique); median FDR and margin of its cells'
    best calls."""
    y, b, u = C.truth.astype(str).to_numpy(), C.best.astype(str).to_numpy(), C.unique.astype(str).to_numpy()
    rb, ru = C.reason_best.astype(str).to_numpy(), C.reason_unique.astype(str).to_numpy()
    fdr = pd.to_numeric(C.fdr, errors="coerce").to_numpy()
    mg = pd.to_numeric(C.margin, errors="coerce").to_numpy()
    rows = []
    for k in ids:
        m = y == k
        if not m.any():
            continue
        called = m & (b != R.U)
        r = dict(identity=k, cells=int(m.sum()), pct_best_correct=100 * float((b[m] == k).mean()),
                 pct_unique_called=100 * float((u[m] != R.U).mean()),
                 pct_boundary=100 * float(((b[m] != R.U) & (u[m] == R.U)).mean()))
        for out, rr, whys in (("best", rb, (R.BACKGROUND, R.FLIP)), ("unique", ru, (R.BACKGROUND, R.AMBIGUOUS, R.FLIP))):
            for w in whys:
                r[f"pct_unlabelled_{out}: {w}"] = 100 * float((rr[m] == w).mean())
        r["median_fdr_best"] = float(np.nanmedian(fdr[called])) if called.any() else np.nan
        r["median_margin_best"] = float(np.nanmedian(mg[called])) if called.any() else np.nan
        rows.append(r)
    return pd.DataFrame(rows)


def genes_per_cell(X):
    import scipy.sparse as sp
    if sp.issparse(X):
        X = X.tocsr()
        return np.diff(X.indptr) if not (X.data <= 0).any() else np.asarray((X > 0).sum(axis=1)).ravel()
    return (np.asarray(X) > 0).sum(axis=1)


def total_counts(X):
    return np.asarray(X.sum(axis=1)).ravel()


def _umap_input(A, idx):
    """The cells idx of A as the UMAP's AnnData (everything copied but the obs columns, as A[idx].copy()). Memory: for
    integer counts in CSR, X is built directly in float32 (scanpy's normalize_total converts integer X to float32 before
    anything else: the same float32 values) instead of a copy of the integer values; the index arrays are shared with
    A's when they are A's own sorted arrays (nothing writes them), copied otherwise."""
    import scipy.sparse as sp
    X = A.X
    if A.is_view or not (sp.issparse(X) and X.format == "csr" and np.issubdtype(X.dtype, np.integer)):
        return A[idx].copy()
    whole = len(idx) == A.n_obs
    Xs = X if whole else X[idx]
    keep = not whole or X.has_sorted_indices
    Ux = sp.csr_matrix((Xs.data.astype(np.float32), Xs.indices if keep else Xs.indices.copy(),
                        Xs.indptr if keep else Xs.indptr.copy()), shape=Xs.shape, copy=False)
    del Xs
    try:   # everything but X, as A[idx].copy()
        A.X = None
        U = A[idx].copy()
    finally:
        A.X = X
    U.X = Ux
    return U


def umap_coords(A, max_cells, seed, log=print):
    """display only: scanpy defaults (normalize_total, log1p, highly_variable_genes, pca, neighbors, umap), seeded;
    a uniform random subsample of max_cells cells (seeded) on larger data"""
    import scanpy as sc
    n = A.n_obs
    rng = np.random.default_rng(seed)
    idx = np.sort(rng.choice(n, size=max_cells, replace=False)) if n > max_cells else np.arange(n)
    U = _umap_input(A, idx)
    U.obs = pd.DataFrame(index=U.obs_names)
    sc.pp.normalize_total(U)
    sc.pp.log1p(U)
    sc.pp.highly_variable_genes(U)
    U = U[:, U.var["highly_variable"].to_numpy()].copy()
    sc.pp.pca(U, random_state=seed)
    sc.pp.neighbors(U, random_state=seed)
    sc.tl.umap(U, random_state=seed)
    log(f"[report] UMAP on {len(idx)} of {n} cells")
    return idx, U.obsm["X_umap"], dict(cells_shown=int(len(idx)), cells_total=int(n), seed=int(seed),
                                        subsampled=bool(n > max_cells),
                                        method="scanpy defaults: normalize_total, log1p, highly_variable_genes, pca, "
                                               "neighbors, umap (random_state = seed); display only, never used for "
                                               "a call")


ALL_GENES = "all genes (scanpy defaults: normalize_total, log1p, highly_variable_genes, pca, neighbors, umap)"


def universe_genes(models):
    """The genes SCARLET decides on: the union of the models' universe genes (their marker genes), sorted."""
    out = set()
    for m in models:
        U = m.universe
        out |= set(map(str, U["gene"] if "gene" in U.columns else U.index))
    return sorted(out)


def umap_spaces(A, max_cells, seed, genes, genes_from, *, gene_col="index", all_genes=True, log=print):
    """Display only: the UMAP of the same cells in each space, seeded; a uniform random subsample of max_cells cells
    (seeded) on larger data. Default space: the universe genes (the marker genes SCARLET decides on) as detection
    (X > 0, as SCARLET reads them), neighbours by Jaccard distance; second space (all_genes): umap_coords on every
    gene. Returns (idx, {'universe': xy, 'all': xy or None}, record)."""
    import scanpy as sc
    n = A.n_obs
    rng = np.random.default_rng(seed)
    idx = np.sort(rng.choice(n, size=max_cells, replace=False)) if n > max_cells else np.arange(n)
    names = A.var_names.astype(str) if gene_col in (None, "index") else A.var[gene_col].astype(str)
    # a name repeated in the data (learn joins such columns: detected in any) maps to every column that carries it
    where = pd.Series(np.arange(len(names))).groupby(np.asarray(names)).apply(list).to_dict()
    groups = [where[g] for g in genes if g in where]
    pos = np.array([g[0] for g in groups], dtype=int) if groups else np.array([], dtype=int)
    cols = np.unique(pos)
    xy = {"universe": None, "all": None}
    spaces = []
    if len(cols) >= 3:
        allc = np.unique(np.concatenate([np.asarray(g) for g in groups]))
        X = A.X[idx][:, allc] if hasattr(A.X, "tocsr") else np.asarray(A.X)[idx][:, allc]
        d_all = (X.toarray() if hasattr(X, "toarray") else np.asarray(X)) > 0
        at = {c: i for i, c in enumerate(allc)}
        # one column per gene, in the order of cols (its first column), detected in any column of its name
        det = np.column_stack([d_all[:, [at[c] for c in g]].any(axis=1) for g in groups])[:, np.argsort(pos, kind="stable")]
        import anndata as ad
        U = ad.AnnData(np.zeros((len(idx), 1), dtype=np.float32))
        U.obsm["X_det"] = det
        sc.pp.neighbors(U, use_rep="X_det", metric="jaccard", random_state=seed)
        sc.tl.umap(U, random_state=seed)
        xy["universe"] = U.obsm["X_umap"]
        spaces.append(dict(key="universe", genes=int(len(cols)),
                           name=f"universe genes of {genes_from} ({len(cols):,} genes: detection, Jaccard)"))
    if all_genes or xy["universe"] is None:
        _, xy["all"], _ = umap_coords(A, max_cells, seed, log=lambda *_: None)
        spaces.append(dict(key="all", name=ALL_GENES))
    log(f"[report] UMAP on {len(idx)} of {n} cells ({', '.join(s['key'] for s in spaces)})")
    rec = dict(cells_shown=int(len(idx)), cells_total=int(n), seed=int(seed), subsampled=bool(n > max_cells),
               spaces=spaces, all_genes_skipped=not all_genes and xy["universe"] is not None,
               method=("; ".join(s["name"] for s in spaces) + "; random_state = seed; display only, never used for "
                       "a call"))
    if xy["universe"] is None:
        xy["universe"], xy["all"] = xy["all"], None
    return idx, xy, rec


def umap_frame(cells, xy):
    """cell, x, y (the default space) and x_all, y_all (all genes) when computed"""
    F = pd.DataFrame({"cell": cells, "x": xy["universe"][:, 0], "y": xy["universe"][:, 1]})
    if xy.get("all") is not None:
        F["x_all"], F["y_all"] = xy["all"][:, 0], xy["all"][:, 1]
    return F


def analytics(C, ids):
    """per identity and entry: cells, correct, wrong, refused by reason (counts and shares)"""
    rows, wrong = [], []
    for e in ENTRIES:
        reason = C[f"reason_{e}"] if e != "forced" else pd.Series(R.CALLED, index=C.index)
        for k in ids:
            m = (C.truth == k).to_numpy()
            p = C[e].to_numpy()[m]
            r = np.asarray(reason)[m]
            n = int(m.sum())
            ok = int((p == k).sum())
            ref = p == R.U
            wr = int((~ref & (p != k)).sum())
            row = dict(entry=e, identity=k, cells=n, correct=ok, wrong=wr, refused=int(ref.sum()))
            for why in (R.BACKGROUND, R.AMBIGUOUS, R.FLIP):
                row[f"refused_{why}"] = int((r == why).sum())
            row["recall"] = ok / n if n else np.nan
            rows.append(row)
            vc = pd.Series(p[~ref & (p != k)]).value_counts()
            for tgt, c in vc.items():
                wrong.append(dict(entry=e, identity=k, called=tgt, cells=int(c), share=c / n))
    return pd.DataFrame(rows), pd.DataFrame(wrong, columns=["entry", "identity", "called", "cells", "share"])


def build(model_path, reference, out, *, self_assessment=True, n_folds=5, seed=0, umap_cells=50000,
          alpha=0.01, ncpu=1, log=print, reference_path=None, note=None, umap=None, umap_all_genes=True):
    """reference: a path or an AnnData (then reference_path names its file, for the hash); note: one line shown
    under the provenance (the transfer report's supplement says which model it reports)"""
    import scarlet
    from scarlet.model import load_model
    m = load_model(model_path)
    S, P = m.settings, m.provenance
    labels = S["parameters"]["labels"]
    rec_in = (P.get("inputs") or {}).get("reference")
    rec_in = rec_in if isinstance(rec_in, dict) else {"path": rec_in}   # a transfer's model records the path only
    if reference is None:
        from scarlet.html_report.transfer_report import resolve_recorded
        reference = resolve_recorded(rec_in.get("path"), S["parameters"].get("model"), model_path)
        if not reference or reference == "<AnnData>" or not os.path.exists(reference):
            raise ValueError("the learn report needs the reference data: the model does not record a readable path "
                             f"({reference!r}); pass reference=")
    if isinstance(reference, (str, os.PathLike)):
        from scarlet._read import read_h5ad   # int64 counts held as int32 (the same values; scarlet._read)
        A, reference_path = read_h5ad(reference), str(reference)
    else:
        A = reference
    if labels not in A.obs:
        raise ValueError(f"the reference has no obs column {labels!r} (the model's labels)")
    if rec_in.get("n_cells") not in (None, A.n_obs) or rec_in.get("n_genes") not in (None, A.n_vars):
        raise ValueError(f"the reference ({A.n_obs} x {A.n_vars}) is not the model's ({rec_in.get('n_cells')} x "
                         f"{rec_in.get('n_genes')})")
    y = A.obs[labels].astype(str).to_numpy()
    ids = list(m.identities)
    col = B.colours(ids)
    bd = B.Bundle(out, "learn")
    # 0 provenance
    res = S["resolved"]
    bl = res.get("blocklist", {})
    prov = dict(
        report=dict(scarlet=scarlet.__version__, written_at=datetime.datetime.now().astimezone().isoformat(
            timespec="seconds")),
        model=dict(path=str(model_path), created_at=P.get("created_at"), host=P.get("host"),
                   scarlet=P.get("versions", {}).get("scarlet"), git_commit=P.get("versions", {}).get("git_commit"),
                   technology=res.get("technology", {}).get("name"), sizes=res.get("sizes"),
                   seed=S["parameters"].get("seed"),
                   blocklist=dict(version=bl.get("version"), md5=bl.get("md5"),
                                  given=bl.get("given"), note=_settings.blocklist_text(bl),
                                  modules_on=(bl.get("modules") or {}).get("on"),
                                  modules_off=(bl.get("modules") or {}).get("off")),
                   imbalanced=res.get("imbalanced", {}).get("handling"),
                   positive_evidence=S["parameters"].get("positive_evidence"),
                   at_package_defaults=None),
        note=note,
        input=dict(reference=B.checked_hash(reference_path, rec_in.get("sha256")),
                   cells=int(A.n_obs), genes=int(A.n_vars), labels=labels, identities=len(ids)))
    # learn at package defaults? every recorded learn parameter equal to api.learn's default
    import inspect
    from scarlet import api
    from scarlet.html_report.selfassess import NOT_SETTINGS
    sig = inspect.signature(api.learn).parameters
    given = {k: v for k, v in S["parameters"].items() if k in sig and k not in NOT_SETTINGS
             and sig[k].default is not inspect.Parameter.empty and v != sig[k].default
             # technology auto records the preset it resolved (technology_resolved)
             and not (k == "technology" and sig[k].default == "auto" and "technology_resolved" in S["parameters"])}
    nd = S["parameters"].get("options_not_default")      # recorded by the run
    if nd is None:
        nd = {k: {"value": v, "default": sig[k].default} for k, v in given.items()}
    prov["model"]["at_package_defaults"] = not nd
    prov["model"]["non_default_parameters"] = nd
    prov["model"]["cell_cycle"] = _rs.cell_cycle_text(res.get("blocklist"))
    prov["summary"] = m.files["summary.txt"].decode() if "summary.txt" in m.files else None
    bd.meta("provenance", prov)
    bd.table("colours", pd.DataFrame({"identity": list(col), "colour": list(col.values())}),
             "one fixed colour per identity, used in every panel; refusals grey " + B.REFUSED_COLOUR)
    # 1 depth
    g = genes_per_cell(A.X)
    t = total_counts(A.X)
    D = pd.DataFrame({"cell": A.obs_names.astype(str), "identity": y, "genes_detected": g, "total_counts": t})
    bd.table("depth", D, "per cell: genes detected (X > 0) and total counts, by the reference label")
    # self-assessment
    if self_assessment:
        from scarlet.html_report import selfassess
        C, srec = selfassess.run(m, A, n_folds=n_folds, seed=seed, alpha=alpha, ncpu=ncpu, log=log)
        bd.table("self_assessment", C, "per cell: the fold, the label, the out-of-fold calls (forced = top1, best, "
                                       "unique), the stored scores and flags, and the refusal reason "
                                       f"({R.LABEL})")
        bd.meta("self_assessment", srec)
        AN, WR = analytics(C, ids)
        bd.table("analytics", AN, "per identity and entry: cells, correct, wrong, refused by reason (out-of-fold)")
        bd.table("wrong_to", WR, "per identity and entry: the identities its wrong calls went to")
        from scarlet.html_report import certainty as C_
        if "runner_up" in C:
            bd.table("boundary", C_.boundary_rows(C.best, C.unique, C.runner_up, "cross-validation"),
                     "out-of-fold boundary cells (a best call, no unique call) by best identity x runner-up")
            bd.table("certainty_per_identity", cv_per_identity(C, ids),
                     "per reference label (out-of-fold): % best correct, % with a unique call, % boundary, % "
                     "unlabelled by reason, median FDR and margin of the best calls")
    else:
        C = None
        bd.meta("self_assessment", dict(off=True, note="self-assessment switched off (--no-self-assessment)"))
    # 2 UMAP
    # umap: umap_spaces(A, ...) already computed by the caller (the transfer report, on the same A and model)
    idx, xy, urec = umap if umap is not None else umap_spaces(
        A, umap_cells, seed, universe_genes([m]), "the model", gene_col=S["parameters"].get("gene_col", "index"),
        all_genes=umap_all_genes, log=log)
    Uc = umap_frame(A.obs_names.astype(str)[idx], xy).assign(label=y[idx])
    if C is not None:
        Ci = C.set_index("cell").loc[Uc.cell]
        for e in ENTRIES:
            Uc[f"cv_{e}"] = Ci[e].to_numpy()
    bd.table("umap", Uc, "UMAP coordinates (display only) with the label and the cross-validated calls")
    bd.meta("umap", urec)
    # 4 drill-down: markers and the universe's detection rates
    mk = m.markers.copy()
    uni = m.universe.copy()
    uni.index = uni.index.astype(str)
    bd.table("markers", mk, "the model's markers per identity (markers.tsv: rank, gene, stability, score)")
    UT = uni.reset_index().rename(columns={uni.index.name or "index": "gene"})
    bd.table("universe", UT, "the model's universe: per marker gene, each identity's detection rate (bootstrap mean) "
                             "and the background rate (bg)")
    # 5 weak-class flags
    pe = (m.files.get("positive_evidence.json") and __import__("json").loads(m.files["positive_evidence.json"])) or {}
    pei = pe.get("identities", {}) if isinstance(pe, dict) else {}
    ss = m.marker_finder_params.get("scoring_stats", {}).get("identities", {})
    W = []
    for k in ids:
        mm = y == k
        r = dict(identity=k, cells=int(mm.sum()), under_50_cells=bool(mm.sum() < UNDER),
                 median_genes_detected=float(np.median(g[mm])) if mm.any() else np.nan,
                 n_eff=ss.get(k, {}).get("n_eff"), markers=int((mk.cell_type.astype(str) == k).sum()),
                 pe_expected=pei.get(k, {}).get("expected"), pe_required=pei.get(k, {}).get("required"),
                 pe_typical_median=pei.get(k, {}).get("typical_median"))
        if C is not None:
            cm = (C.truth == k).to_numpy()
            for e in ENTRIES:
                r[f"cv_recall_{e}"] = float((C[e].to_numpy()[cm] == k).mean()) if cm.any() else np.nan
            wb = WR[(WR.entry == "best") & (WR.identity == k)].sort_values("cells", ascending=False)
            r["most_confused_with"] = wb.called.iloc[0] if len(wb) else ""
            r["most_confused_share"] = float(wb.share.iloc[0]) if len(wb) else 0.0
            r["refused_best_share"] = float((C.best.to_numpy()[cm] == R.U).mean()) if cm.any() else np.nan
            if "pe_detected" in C:
                r["pe_median_detected_cv"] = float(np.nanmedian(C.pe_detected.to_numpy(float)[cm])) if cm.any() \
                    else np.nan
            if "pe_low" in C:
                r["pe_low_share_cv"] = float(R._bool(C.pe_low.to_numpy()[cm]).mean()) if cm.any() else np.nan
        W.append(r)
    W = pd.DataFrame(W)
    key = "cv_recall_best" if C is not None else "cells"
    W = W.sort_values([key, "cells"], ascending=[True, True], kind="mergesort").reset_index(drop=True)
    W.insert(0, "rank", np.arange(1, len(W) + 1))
    bd.table("weak_classes", W, f"identities ranked weakest first by {key} (ties: fewer cells first); evidence only, "
                                "no threshold; under_50_cells is a reporting convention, not a SCARLET "
                                "threshold")
    bd.meta("weak_rank_key", key)
    bd.meta("identities", ids)
    _figures(bd, D, Uc, C, AN if C is not None else None, col, ids)
    bd.save()
    return bd


def _figures(bd, D, Uc, C, AN, col, ids):
    """the bundle's static figures (matplotlib, from the same tables)"""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    h = max(3.0, 0.28 * len(ids) + 1)
    fig, ax = plt.subplots(figsize=(7, h))
    data = [D.genes_detected[D.identity == k].to_numpy() for k in ids]
    bp = ax.boxplot(data, vert=False, patch_artist=True, showfliers=False)
    for p, k in zip(bp["boxes"], ids):
        p.set_facecolor(col[k])
        p.set_alpha(0.7)
    ax.set_yticks(range(1, len(ids) + 1), ids, fontsize=7)
    ax.set_xlabel("genes detected per cell")
    ax.set_title("Depth per identity")
    bd.figure("depth", fig, "genes detected per cell, per identity")
    plt.close(fig)
    panels = [("label", "Labels")] + ([("cv_best", "Cross-validated best"), ("cv_unique", "Cross-validated unique")]
                                      if C is not None else [])
    fig, axes = plt.subplots(1, len(panels), figsize=(5 * len(panels), 5), squeeze=False)
    for ax, (c, title) in zip(axes[0], panels):
        v = Uc[c].astype(str).to_numpy()
        cc = np.array([col.get(x, B.REFUSED_COLOUR) for x in v])
        order = np.argsort(v == R.U)[::-1]
        ax.scatter(Uc.x.to_numpy()[order], Uc.y.to_numpy()[order], s=2, c=cc[order], linewidths=0)
        ax.set_title(title)
        ax.set_xticks([])
        ax.set_yticks([])
    bd.figure("umap", fig, "UMAPs (display only)")
    plt.close(fig)
    if AN is not None:
        fig, axes = plt.subplots(1, 3, figsize=(15, h), sharey=True)
        for ax, e in zip(axes, ENTRIES):
            a = AN[AN.entry == e].set_index("identity").reindex(ids)
            n = a.cells.replace(0, np.nan)
            left = np.zeros(len(ids))
            for c, lab, cl in (("correct", "correct", "#4CAF50"), ("wrong", "wrong", "#E53935"),
                               (f"refused_{R.BACKGROUND}", "refused: background", "#9E9E9E"),
                               (f"refused_{R.AMBIGUOUS}", "refused: ambiguity", "#BDBDBD"),
                               (f"refused_{R.FLIP}", "refused: flip test", "#616161")):
                w = (a[c] / n).fillna(0).to_numpy()
                ax.barh(range(len(ids)), w, left=left, color=cl, label=lab)
                left += w
            ax.set_title(e)
            ax.set_xlim(0, 1)
        axes[0].set_yticks(range(len(ids)), ids, fontsize=7)
        axes[0].invert_yaxis()
        axes[1].legend(fontsize=7, loc="upper center", bbox_to_anchor=(0.5, -0.12), ncol=5, frameon=False)
        bd.figure("analytics", fig, "out-of-fold outcome shares per identity")
        plt.close(fig)
