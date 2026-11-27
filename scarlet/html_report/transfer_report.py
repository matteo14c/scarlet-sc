"""The transfer report: main page = transfer A -> B, supplement = the learn
report on A (its self-assessment).

Main sections: 0 provenance for both sides; 1 depth of A and B, before and after depth matching (per cell before; the
record's medians after); 2 UMAPs: A by its labels, B with the anchors highlighted, B with the final calls (refusals
grey); 3 analytics on the final calls (best / unique / refused by reason, overall and per identity; identities of A
with no calls in B); 4 identity drill-down (A and B marker heatmaps; A-to-B marker comparison; markers before vs after
refinement; stable anchors); 5 weak-class flags from the reference and from the transfer; 6 what refinement changed
(labels changed, cells refused in the first pass but called after refinement); 7 agreement across the draws (per-cell
agreement, the draws' own results). The refusal reasons are the ones transfer recorded per cell (background, ambiguity,
undetermined across draws, flip test).
"""
import json
import os

import numpy as np
import pandas as pd

from scarlet import run_summary as _rs
from scarlet import settings as _settings
from scarlet.html_report import bundle as B
from scarlet.html_report import reasons as R
from scarlet.html_report.learn_report import genes_per_cell, total_counts, umap_frame, umap_spaces, universe_genes


def resolve_recorded(recorded, recorded_anchor, actual_anchor):
    """a path the run recorded relative to its working directory, found from where the run's output (or model) is
    now: anchor = the output path as recorded and as it is"""
    if not recorded or recorded == "<AnnData>" or os.path.isabs(str(recorded)) or not recorded_anchor:
        return recorded
    rel = os.path.relpath(str(recorded), start=str(recorded_anchor))
    return os.path.normpath(os.path.join(str(actual_anchor), rel))


IN_MEMORY = dict(path="<AnnData>", note="in memory, not hashed")


def _read(path_or_none, recorded, what):
    """The data as given (an AnnData in memory: returned with no path) or read from the path given or recorded."""
    import anndata as ad
    if isinstance(path_or_none, ad.AnnData):
        return path_or_none, None
    p = path_or_none or recorded
    if not p or p == "<AnnData>" or not os.path.exists(str(p)):
        raise ValueError(f"the transfer report needs the {what} data: the record has no readable path ({recorded!r}); "
                         f"pass {what}=")
    from scarlet._read import read_h5ad   # int64 counts held as int32 (the same values; scarlet._read)
    return read_h5ad(p), p


def _model_files(run):
    fs = sorted(f for f in os.listdir(run) if f.endswith(".scarlet"))
    ref = [f for f in fs if f.startswith("reference")]
    ref_by = {("all" if f == "reference.scarlet" else f[len("reference_"):-len(".scarlet")]): f for f in ref}
    qry_by = {f[len("refine_"):-len(".scarlet")]: f for f in fs if f.startswith("refine_")}
    return ref_by, qry_by


def detection_in(X, cells_mask, gene_idx):
    """share of the given cells detecting each gene (X > 0): a description of the query data"""
    import scipy.sparse as sp
    sub = X[np.where(cells_mask)[0]][:, gene_idx]
    if sp.issparse(sub):
        return np.asarray((sub > 0).mean(axis=0)).ravel()
    return (np.asarray(sub) > 0).mean(axis=0)


def build(run, reference, query, out, *, self_assessment=True, n_folds=5, seed=0, umap_cells=50000,
          ncpu=1, log=print, query_labels=None, umap_all_genes=True):
    import scarlet
    from scarlet.html_report import learn_report
    from scarlet.model import load_model
    rec = json.load(open(os.path.join(run, "transfer_record.json")))
    par = rec["parameters"]
    out_rec = par.get("output")
    A, ref_path = _read(reference, resolve_recorded(par.get("reference"), out_rec, run), "reference")
    Q, qry_path = _read(query, resolve_recorded(par.get("query"), out_rec, run), "query")
    A.obs_names, Q.obs_names = A.obs_names.astype(str), Q.obs_names.astype(str)
    labels = par["labels"]
    FP = pd.read_csv(os.path.join(run, "first_pass.tsv"), sep="\t", dtype={"cell": str}, keep_default_na=False)
    LB = pd.read_csv(os.path.join(run, "labels.tsv"), sep="\t", dtype={"cell": str}, keep_default_na=False)
    if not (FP.cell.tolist() == LB.cell.tolist() and set(FP.cell) == set(Q.obs_names)):
        raise ValueError("the query data are not the run's query (cells differ from labels.tsv)")
    order = pd.Index(Q.obs_names, name="cell")
    FP = FP.set_index("cell").loc[order].reset_index()
    LB = LB.set_index("cell").loc[order].reset_index()
    ref_by, qry_by = _model_files(run)
    refm = {o: load_model(os.path.join(run, f)) for o, f in ref_by.items()}
    qm = {o: load_model(os.path.join(run, f)) for o, f in qry_by.items()}
    m0 = next(iter(refm.values()))
    ids = list(m0.identities)
    col = B.colours(ids)
    truth = None
    if query_labels is not None:   # the user's labels: read before any work, so a missing column stops early
        from scarlet.html_report import known as K
        truth = K.labels_of(Q.obs, query_labels)
    bd = B.Bundle(out, "transfer")
    if truth is not None:
        K.write(bd, truth, {"final best": LB.scarlet_best, "final unique": LB.scarlet_unique,
                            "first pass best": FP.scarlet_best, "first pass unique": FP.scarlet_unique},
                ids, query_labels)
    tflip = par.get("flip")
    # reasons: as transfer recorded them per cell
    for stage, T in (("first", FP), ("final", LB)):
        rb, ru = T.scarlet_best_reason.astype(str).to_numpy(), T.scarlet_unique_reason.astype(str).to_numpy()
        unknown = sorted((set(rb) | set(ru)) - set(R.RECORDED) - {R.CALLED})
        if unknown:
            raise ValueError(f"transfer {stage}: reasons this report does not know: {unknown}")
        R.check_exact(T.scarlet_best, T.scarlet_unique, rb, ru, where=f"transfer {stage}")
        T["reason_best"], T["reason_unique"] = rb, ru
    # stable anchors: labels.tsv's scarlet_stable_anchor flag (a stable anchor keeps its label, so its cell type is its
    # final call), checked against the record per origin and cell type
    groups = (rec.get("refine") or {}).get("groups") or []
    anchors = LB["scarlet_stable_anchor"].astype(str).str.lower().eq("true").to_numpy() \
        if "scarlet_stable_anchor" in LB.columns else np.zeros(len(LB), dtype=bool)
    sa = np.where(anchors, LB.scarlet_unique.astype(str).to_numpy(), "")
    anc_rows = []
    for g in groups:
        o = g["origin"]
        om = (LB.origin.astype(str) == str(o)).to_numpy()
        refined = set(g.get("refined", []))
        for k, v in (g.get("stable_anchors") or {}).items():
            n_k = int((om & (sa == k)).sum())
            if n_k != int(v):
                raise AssertionError(f"stable anchors of {k} ({o}): record {v}, labels.tsv {n_k}")
            anc_rows.append(dict(origin=o, identity=k, anchors=int(v), refined=k in refined and "skipped" not in g))
    ANC = pd.DataFrame(anc_rows, columns=["origin", "identity", "anchors", "refined"])
    bd.table("anchors", ANC, "per origin and identity: stable anchors (cells called uniquely to it in every draw of the "
                             "run) and whether it was refined "
                             f"(transfer_record.json; SCARLET's minimum {(rec.get('refine') or {}).get('min_cells')} "
                             "stable anchors, and at least 2 refined cell types per origin)")
    # 0 provenance
    rin, qin = ((rec.get("inputs") or {}).get(k) or {} for k in ("reference", "query"))   # recorded by the run
    v = rec.get("versions", {})
    st = rec.get("settings", {})
    bl = st.get("blocklist", {})
    prov = dict(
        report=dict(scarlet=scarlet.__version__),
        transfer=dict(folder=str(run), scarlet=v.get("scarlet"), git_commit=v.get("git_commit"),
                      refine=par.get("refine"), depth_match=par.get("depth_match"),
                      draws=(rec.get("draws") or {}).get("n"), seeds=(rec.get("draws") or {}).get("seeds"),
                      mode=(rec.get("draws") or {}).get("mode"), flip=tflip, alpha=par.get("alpha"),
                      positive_evidence=par.get("positive_evidence"),
                      technology=(st.get("technology") or {}).get("name"), sizes=st.get("sizes"),
                      blocklist=dict(version=bl.get("version"), md5=bl.get("md5"), given=bl.get("given"),
                                     note=_settings.blocklist_text(bl),
                                     modules_on=(bl.get("modules") or {}).get("on"),
                                     modules_off=(bl.get("modules") or {}).get("off")),
                      origins=rec.get("origins"), qc=rec.get("qc")),
        reference=dict(file=B.checked_hash(ref_path, rin.get("sha256")) or IN_MEMORY, cells=int(A.n_obs),
                       genes=int(A.n_vars), labels=labels,
                       identities=len(ids)),
        query=dict(file=B.checked_hash(qry_path, qin.get("sha256")) or IN_MEMORY, cells=int(Q.n_obs),
                   genes=int(Q.n_vars)))
    import inspect
    from scarlet import api
    sig = inspect.signature(api.transfer).parameters
    skip = {"reference", "labels", "query", "output", "ncpu", "report", "verbosity"}
    pub = dict(par)   # the recorded values as the public options name them
    pub["refine"] = {"anchored": True, "no": False}.get(par.get("refine"), par.get("refine"))
    pub["depth_match"] = {"capture": True, "none": False}.get(par.get("depth_match"), par.get("depth_match"))
    if ((rec.get("preflight") or {}).get("technology_resolved") or {}).get("how", "as given") != "as given":
        pub["technology"] = "auto"
    given = {k: x for k, x in pub.items() if k in sig and k not in skip and sig[k].default is not inspect.Parameter.empty
             and x != sig[k].default and not (sig[k].default is None)}
    nd = rec.get("options_not_default")                  # recorded by the run
    if nd is None:
        nd = {k: {"value": x, "default": sig[k].default} for k, x in given.items()}
    prov["transfer"]["at_package_defaults"] = not nd
    prov["transfer"]["non_default_parameters"] = nd
    cc = rec.get("cell_cycle") or {}
    ref_cc = cc.get("reference") or _rs.cell_cycle_text(bl)
    others = [f"{e['learn']}: {e['cell_cycle_genes']}" for e in cc.get("learns", [])
              if e["cell_cycle_genes"] != ref_cc]
    prov["transfer"]["cell_cycle"] = (ref_cc + (f"; where the other learns differ: {'; '.join(others)}" if others
                                                else "")) if ref_cc else None
    prov["summary"] = rec.get("summary")
    bd.meta("provenance", prov)
    bd.table("colours", pd.DataFrame({"identity": list(col), "colour": list(col.values())}),
             "one fixed colour per identity, used in every panel; refusals grey " + B.REFUSED_COLOUR)
    # 1 depth
    ga, gq = genes_per_cell(A.X), genes_per_cell(Q.X)
    D = pd.concat([pd.DataFrame({"side": "A (reference)", "cell": A.obs_names, "group": A.obs[labels].astype(str),
                                 "genes_detected": ga, "total_counts": total_counts(A.X)}),
                   pd.DataFrame({"side": "B (query)", "cell": Q.obs_names, "group": LB.scarlet_best.to_numpy(),
                                 "genes_detected": gq, "total_counts": total_counts(Q.X)})], ignore_index=True)
    bd.table("depth", D, "per cell, before depth matching: genes detected (X > 0); A grouped by its labels, B by the "
                         "final best call")
    dm_rows = _depth_rows(rec)
    DM = pd.DataFrame(dm_rows)
    bd.table("depth_match", DM, "depth matching per draw and origin, as recorded by transfer (median genes per cell on "
                                "the shared genes before; depth_after on the thinned side)")
    # 2 UMAPs (display only)
    # in SCARLET's own space: A on its model's universe genes, B on the refined model's (the final calls' model;
    # the reference model's when nothing was refined); all genes as the second space
    gcol = par.get("gene_col", "index")
    ia, xa, ua = umap_spaces(A, umap_cells, seed, universe_genes(refm.values()), "A's model", gene_col=gcol,
                             all_genes=umap_all_genes, log=log)
    iq, xq, uq = umap_spaces(Q, umap_cells, seed, universe_genes((qm or refm).values()),
                             "the refined model" if qm else "A's model (nothing refined)", gene_col=gcol,
                             all_genes=umap_all_genes, log=log)
    UA = umap_frame(A.obs_names[ia], xa).assign(label=A.obs[labels].astype(str).to_numpy()[ia])
    UQ = umap_frame(Q.obs_names[iq], xq).assign(anchor=np.where(anchors[iq], sa[iq], R.U),
                                                final_best=LB.scarlet_best.to_numpy()[iq],
                                                final_unique=LB.scarlet_unique.to_numpy()[iq])
    bd.table("umap_A", UA, "UMAP of A (display only), by its labels")
    bd.table("umap_B", UQ, "UMAP of B (display only): stable anchors, final calls")
    bd.meta("umap", dict(A=ua, B=uq))
    # 3 analytics on the final calls
    rows = []
    for e in ("best", "unique"):
        c = LB[f"scarlet_{e}"].to_numpy()
        for k in ids:
            rows.append(dict(entry=e, identity=k, called=int((c == k).sum()), share=float((c == k).mean())))
        r = LB[f"reason_{e}"].to_numpy()
        for why in sorted(set(r) - {R.CALLED}):
            rows.append(dict(entry=e, identity=f"refused: {why}", called=int((r == why).sum()),
                             share=float((r == why).mean())))
    AN = pd.DataFrame(rows)
    bd.table("analytics", AN, "final calls per identity and refusals by reason (as transfer recorded them), best and "
                              "unique")
    none = [dict(entry=e, identity=k) for e in ("best", "unique") for k in ids
            if not (LB[f"scarlet_{e}"].to_numpy() == k).any()]
    bd.table("no_calls", pd.DataFrame(none, columns=["entry", "identity"]), "identities of A with no call in B")
    # 4 drill-down
    genes_q = pd.Index(Q.var_names.astype(str))
    MK, SC = [], []
    for o, m in refm.items():
        mr = m.markers.assign(origin=o, model=ref_by[o])
        MK.append(mr)
    for o, m in qm.items():
        MK.append(m.markers.assign(origin=o, model=qry_by[o]))
    MK = pd.concat(MK, ignore_index=True)
    bd.table("markers", MK, "markers per identity of the reference model(s) (before refinement) and the refined "
                            "model(s) (after); model file named")
    UNI = []
    for name, m in [(ref_by[o], refm[o]) for o in refm] + [(qry_by[o], qm[o]) for o in qm]:
        u = m.universe.copy()
        u.index = u.index.astype(str)
        UNI.append(u.reset_index().rename(columns={u.index.name or "index": "gene"}).assign(model=name))
    bd.table("universe", pd.concat(UNI, ignore_index=True),
             "each model's universe: per marker gene, each identity's detection rate and the background (bg)")
    # the drill-down's values: each model's in-class table (learn's bootstrap estimate of the share of each cell type's
    # cells detecting a gene, for every gene its marker search scored; the universe is its marker rows), on the genes
    # of the drill-down (every model's markers); models written before 0.6.3 have none
    drill_genes = set(MK.gene.astype(str))
    INC = []
    for name, m in [(ref_by[o], refm[o]) for o in refm] + [(qry_by[o], qm[o]) for o in qm]:
        inc = m.in_class
        if inc is not None:
            inc.index = inc.index.astype(str)
            inc = inc.loc[[g for g in inc.index if g in drill_genes]]
            INC.append(inc.reset_index().rename(columns={inc.index.name or "index": "gene"}).assign(model=name))
    if INC:
        bd.table("in_class", pd.concat(INC, ignore_index=True),
                 "each model's in-class table on the drill-down's genes (every model's markers): per gene, the share "
                 "of each identity's cells detecting it as learn estimated it (the bootstrap average; the universe "
                 "holds the same values for the model's own markers) and the background (bg)")
    fb = LB.scarlet_best.to_numpy()
    for o, m in refm.items():
        u = m.universe
        u.index = u.index.astype(str)
        for k in ids:
            g = m.markers[m.markers.cell_type.astype(str) == k].sort_values("idx").gene.astype(str).tolist()
            g = [x for x in g if x in u.index]
            called = fb == k
            gi = genes_q[~genes_q.duplicated()].get_indexer(g)   # a repeated name: its first column
            gi = np.where(gi >= 0, np.flatnonzero(~genes_q.duplicated())[np.maximum(gi, 0)], -1)
            ok = gi >= 0
            det = np.full(len(g), np.nan)
            if called.any() and ok.any():
                det[ok] = detection_in(Q.X, called, gi[ok])
            for x, a_rate, b_rate, present in zip(g, u.loc[g, k].to_numpy(), det, ok):
                SC.append(dict(origin=o, identity=k, gene=x, detection_A=float(a_rate), detection_B_called=b_rate,
                               in_query=bool(present), B_cells_called=int(called.sum())))
    bd.table("marker_A_vs_B", pd.DataFrame(SC), "per reference marker: its detection rate in A (the reference "
                                                "model's universe) and the share of B cells called that identity "
                                                "(final best) that detect it (computed by the report from the query "
                                                "data, X > 0)")
    # 5 weak-class flags: the reference (from the supplement) and the transfer
    qc = pd.read_csv(os.path.join(run, "qc.tsv"), sep="\t") if os.path.exists(os.path.join(run, "qc.tsv")) else None
    if qc is not None:
        bd.table("qc", qc, "transfer's own QC per identity (qc.tsv: conserved markers, SCARLET's 'doubtful transfer' "
                           "flag)")
    ch = []
    for k in ids:
        f = FP.scarlet_best.to_numpy() == k
        ch.append(dict(identity=k, first_pass_best=int(f.sum()),
                       changed_by_refinement=int((f & (fb != k) & (fb != R.U)).sum()),
                       refused_after=int((f & (fb == R.U)).sum())))
    bd.table("changed_per_identity", pd.DataFrame(ch), "per first-pass identity: cells, labels changed by refinement, "
                                                       "cells refused after refinement")
    # 6 what refinement changed
    X6 = []
    for e in ("best", "unique"):
        a, b = FP[f"scarlet_{e}"].to_numpy(), LB[f"scarlet_{e}"].to_numpy()
        ct = pd.crosstab(pd.Series(a, name="first_pass"), pd.Series(b, name="final"))
        for i in ct.index:
            for j in ct.columns:
                if ct.loc[i, j] and i != j:
                    X6.append(dict(entry=e, first_pass=i, final=j, cells=int(ct.loc[i, j])))
    bd.table("refinement_changes", pd.DataFrame(X6, columns=["entry", "first_pass", "final", "cells"]),
             "cells whose call changed between the first pass and the final labels (refused = 'unassigned')")
    six = {e: dict(refused_then_called=int(((FP[f"scarlet_{e}"] == R.U) & (LB[f"scarlet_{e}"] != R.U)).sum()),
                   called_then_refused=int(((FP[f"scarlet_{e}"] != R.U) & (LB[f"scarlet_{e}"] == R.U)).sum()),
                   changed_identity=int(((FP[f"scarlet_{e}"] != R.U) & (LB[f"scarlet_{e}"] != R.U) &
                                         (FP[f"scarlet_{e}"] != LB[f"scarlet_{e}"])).sum()))
           for e in ("best", "unique")}
    bd.meta("refinement_summary", six)
    # first pass and refined, best and unique, side by side (both documented outputs: first_pass.tsv and
    # labels.tsv); cells per call and share of B
    n = len(FP)
    SB = []
    for k in list(ids) + [R.U]:
        r = dict(identity=k)
        for e in ("best", "unique"):
            for side, X in (("first_pass", FP), ("refined", LB)):
                c = int((X[f"scarlet_{e}"].to_numpy() == k).sum())
                r[f"{side}_{e}"] = c
                r[f"{side}_{e}_pct"] = round(100 * c / n, 2) if n else float("nan")
        SB.append(r)
    bd.table("first_pass_vs_refined", pd.DataFrame(SB), "cells per call in B, first pass (first_pass.tsv) and refined "
                                                         "(labels.tsv), best and unique; 'unassigned' = refused")
    # per query cell: calls, agreement, reasons, runner-up, stable anchor and flags, first pass and final
    n_dr = int((rec.get("draws") or {}).get("n") or 0)

    def cols(X, pre):
        g = lambda c: X[c].to_numpy() if c in X.columns else np.full(len(X), "", dtype=object)  # noqa: E731
        return {f"{pre}best": X.scarlet_best.to_numpy(), f"{pre}unique": X.scarlet_unique.to_numpy(),
                f"{pre}best_agreement": g("scarlet_best_agreement"),
                f"{pre}unique_agreement": g("scarlet_unique_agreement"),
                f"{pre}reason_best": X.reason_best.to_numpy(), f"{pre}reason_unique": X.reason_unique.to_numpy(),
                f"{pre}top1": g("scarlet_top1"), f"{pre}runner_up": g("scarlet_runner_up"),
                f"{pre}runner_up_draws": g("scarlet_runner_up_draws")}

    g_ = lambda X, c: X[c].to_numpy() if c in X.columns else np.full(len(X), "", dtype=object)  # noqa: E731
    QC = pd.DataFrame(dict(cell=LB.cell.to_numpy(), origin=LB.origin.to_numpy(), **cols(LB, ""), **cols(FP, "first_"),
                           stable_anchor=np.where(anchors, sa, ""), anchor=anchors, flip=g_(LB, "scarlet_flip"),
                           first_flip=g_(FP, "scarlet_flip"), pe_low=g_(LB, "scarlet_pe_low"),
                           refine_note=g_(LB, "scarlet_refine_note")))
    bd.table("query_cells", QC, "per query cell: final and first-pass calls with their agreement (draws giving the "
                                "consensus outcome, of N), reasons, top-ranked identity and runner-up (with its "
                                "draws), stable anchor, flags")
    QU = QC.iloc[iq].reset_index(drop=True)
    for c in [c for c in ("x", "y", "x_all", "y_all") if c in UQ]:
        QU[c] = UQ[c].to_numpy()
    bd.table("query_umap", QU,
             "the query UMAP's cells (display only) with their per-cell values")
    from scarlet.html_report import certainty as C_
    BR = pd.concat([C_.boundary_rows(FP.scarlet_best, FP.scarlet_unique, g_(FP, "scarlet_runner_up"), "first pass"),
                    C_.boundary_rows(LB.scarlet_best, LB.scarlet_unique, g_(LB, "scarlet_runner_up"), "refined")],
                   ignore_index=True)
    bd.table("boundary", BR, "boundary cells (a best call, no unique call) by best identity x runner-up, first pass "
                             "and refined")
    bd.table("per_identity", per_identity(LB, ids, groups, ANC, n_dr, (rec.get("refine") or {}).get("min_cells")),
             "per identity: calls, stable anchors, refined or not (why), % of its best calls at full agreement and "
             "boundary, % unlabelled by reason among the cells it ranks first")
    bd.meta("n_draws", n_dr)
    # 7 agreement across the draws
    n_dr = (rec.get("draws") or {}).get("n")
    AG = []
    if n_dr:
        for side, X in (("first_pass", FP), ("refined", LB)):
            for e in ("best", "unique"):
                c = f"scarlet_{e}_agreement"
                if c in X.columns:
                    vc = X[c].astype(int).value_counts().sort_index()
                    for k, v in vc.items():
                        AG.append(dict(output=f"{side} {e}", agreement=int(k), draws=int(n_dr), cells=int(v),
                                       share=round(float(v) / len(X), 4)))
    bd.table("agreement", pd.DataFrame(AG, columns=["output", "agreement", "draws", "cells", "share"]),
             "per output: cells by agreement (the number of draws giving the cell's most frequent outcome)")
    bd.table("per_draw", pd.DataFrame((rec.get("draws") or {}).get("per_draw") or []),
             "each draw's own result: cells labelled, first pass and refined, best and unique (the spread across the "
             "draws)")
    bd.meta("identities", ids)
    bd.meta("supplement", "supplement_A/report.html")
    _figures(bd, D, UA, UQ, col)
    bd.save()
    # supplement: the learn report on A (the reference model of the transfer)
    o0 = next(iter(refm))
    sup = os.path.join(out, "supplement_A")
    learn_report.build(os.path.join(run, ref_by[o0]), A, sup, self_assessment=self_assessment,
                       n_folds=n_folds, seed=seed, umap_cells=umap_cells, alpha=par.get("alpha", 0.01), ncpu=ncpu,
                       log=log, reference_path=ref_path, umap=(ia, xa, ua),
                       note=("the transfer's reference model (" + ref_by[o0] + ")" +
                             ("; the reference was depth-matched (thinned) before learning, the self-assessment "
                              "learns on A as given" if any(d.get("thinned") == "reference"
                                                             for d in dm_rows) else "")))
    return bd


def per_identity(LB, ids, groups, ANC, n_draws, min_cells):
    """The per-identity table on the final labels: cells called best / unique; stable anchors; refined or not
    (with the reason); of its best calls, % at full agreement (all N draws) and % boundary (no unique call); among
    the cells whose most frequent top-ranked identity (scarlet_top1) it is, % unlabelled by reason, best and unique."""
    b, u = LB.scarlet_best.astype(str).to_numpy(), LB.scarlet_unique.astype(str).to_numpy()
    kb = pd.to_numeric(LB.get("scarlet_best_agreement"), errors="coerce").to_numpy() \
        if "scarlet_best_agreement" in LB.columns else np.full(len(LB), np.nan)
    t1 = LB["scarlet_top1"].astype(str).to_numpy() if "scarlet_top1" in LB.columns else b
    rb, ru = LB.reason_best.astype(str).to_numpy(), LB.reason_unique.astype(str).to_numpy()
    rows = []
    for k in ids:
        st = []
        for g in groups:
            if "skipped" in g:
                why = "origin not refined (fewer than 2 cell types anchored stably)"
            elif k in set(g.get("refined", [])):
                why = "yes"
            else:
                why = f"no: not anchored stably (fewer than {min_cells} stable anchors)"
            st.append(why if len(groups) == 1 else f"{g['origin']}: {why}")
        mb, mt = b == k, t1 == k
        r = dict(identity=k, called_best=int(mb.sum()), called_unique=int((u == k).sum()),
                 stable_anchors=int(ANC[ANC.identity.astype(str) == k].anchors.sum()) if len(ANC) else 0,
                 refined="; ".join(st) if st else "refinement off",
                 pct_full_agreement=100 * float((kb[mb] == n_draws).mean()) if mb.any() and n_draws else np.nan,
                 pct_boundary=100 * float((u[mb] == R.U).mean()) if mb.any() else np.nan,
                 ranked_first=int(mt.sum()))
        for out, rr, whys in (("best", rb, (R.BACKGROUND, R.UNDETERMINED, R.FLIP)),
                              ("unique", ru, (R.BACKGROUND, R.AMBIGUOUS, R.UNDETERMINED, R.FLIP))):
            for w in whys:
                r[f"pct_unlabelled_{out}: {w}"] = 100 * float((rr[mt] == w).mean()) if mt.any() else np.nan
        rows.append(r)
    return pd.DataFrame(rows)


def _depth_rows(rec):
    """The record's depth-matching plans as rows (one per draw and origin)."""
    return [dict(draw=x.get("draw"), **o) for x in rec.get("depth_match") or [] for o in x["origins"]]


def _figures(bd, D, UA, UQ, col):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(7, 3))
    sides = list(dict.fromkeys(D.side))
    ax.boxplot([D.genes_detected[D.side == s].to_numpy() for s in sides], vert=False, showfliers=False)
    ax.set_yticks(range(1, len(sides) + 1), sides)
    ax.set_xlabel("genes detected per cell (before depth matching)")
    bd.figure("depth", fig, "genes detected per cell, A and B")
    plt.close(fig)
    panels = [(UA, "label", "A: labels"), (UQ, "anchor", "B: stable anchors"), (UQ, "final_best", "B: final best"),
              (UQ, "final_unique", "B: final unique")]
    fig, axes = plt.subplots(1, 4, figsize=(20, 5))
    for ax, (U, c, t) in zip(axes, panels):
        v = U[c].astype(str).to_numpy()
        cc = np.array([col.get(x, B.REFUSED_COLOUR) for x in v])
        order = np.argsort(v == R.U)[::-1]
        ax.scatter(U.x.to_numpy()[order], U.y.to_numpy()[order], s=3, c=cc[order], linewidths=0)
        ax.set_title(t)
        ax.set_xticks([])
        ax.set_yticks([])
    bd.figure("umap", fig, "UMAPs of A and B (display only)")
    plt.close(fig)
