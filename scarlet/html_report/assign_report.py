"""The assign view: a query without labels, the learn report's layout with calls and refusals
only. Input: the annotated .h5ad that assign wrote, and the model it used. Refusal reasons from the stored scores
(reasons.assign_reasons) with the alpha the assign run recorded (uns['scarlet_options']); the exactness check fails,
saying so, if they do not reproduce the stored calls."""
import json
import os

import numpy as np
import pandas as pd

from scarlet import run_summary as _rs
from scarlet import settings as _settings
from scarlet.html_report import bundle as B
from scarlet.html_report import reasons as R
from scarlet.html_report.learn_report import genes_per_cell, total_counts, umap_frame, umap_spaces, universe_genes
from scarlet.html_report.transfer_report import detection_in


def build(annotated, model, out, *, alpha=0.01, seed=0, umap_cells=50000, log=print, query_labels=None,
          umap_all_genes=True):
    import scarlet
    from scarlet.model import load_model
    if model is None:
        raise ValueError("the assign view needs the model the run used: pass model=")
    m = load_model(model)
    from scarlet._read import read_h5ad   # int64 counts held as int32 (the same values; scarlet._read)
    Z = read_h5ad(annotated)
    Z.obs_names = Z.obs_names.astype(str)
    obs = Z.obs
    for c in ("scarlet_best", "is_scarlet_unique", "scarlet_fdr", "scarlet_top1", "scarlet_ll_top1", "scarlet_ll_bg"):
        if c not in obs:
            raise ValueError(f"{annotated}: obs {c!r} missing (not an annotated output of scarlet assign)")
    names = list(Z.uns.get("scarlet_identity_names", m.identities))
    if set(map(str, names)) != set(map(str, m.identities)):
        raise ValueError("the annotated data were not assigned with this model (identities differ)")
    ids = list(m.identities)
    col = B.colours(ids)
    best = obs.scarlet_best.astype(str).to_numpy()
    uniq = np.where(R._bool(obs.is_scarlet_unique), best, R.U)
    truth = None
    if query_labels is not None:   # the user's labels (an obs column of the annotated data); display only
        from scarlet.html_report import known as K
        truth = K.labels_of(obs, query_labels)
    prof = dict(Z.uns.get("scarlet_assign", {}))
    rec = dict(Z.uns.get("scarlet_options", {}))      # the run's alpha and its options not at their defaults
    nd = json.loads(str(rec["options_not_default"])) if "options_not_default" in rec else None
    if "alpha" in rec:
        if alpha not in (0.01, float(rec["alpha"])):
            raise ValueError(f"the assign run recorded alpha {float(rec['alpha'])}; the report uses it (got {alpha})")
        alpha, alpha_note = float(rec["alpha"]), "recorded by the assign run"
    else:
        alpha_note = "the run records no alpha: the value given to the report; the stored calls are reproduced with it"
    fmode = str(prof.get("flip", "off"))
    rb, ru = R.assign_reasons(best, uniq, obs.scarlet_ll_top1, obs.scarlet_ll_bg, obs.scarlet_fdr, alpha,
                              flip=obs["scarlet_flip"] if "scarlet_flip" in obs else None, flip_mode=fmode)
    try:
        R.check_exact(best, uniq, rb, ru, where="assign view")
    except AssertionError as e:
        raise ValueError(f"{e}. The refusal reasons do not reproduce the stored calls with alpha {alpha}: pass the "
                         "alpha the assign run used") from None
    bd = B.Bundle(out, "assign")
    if truth is not None:
        K.write(bd, truth, {"best": best, "unique": uniq}, ids, query_labels)
    sp = m.settings
    bl = sp["resolved"].get("blocklist", {})
    bd.meta("provenance", dict(
        report=dict(scarlet=scarlet.__version__),
        model=dict(path=str(model), scarlet=m.provenance.get("versions", {}).get("scarlet"),
                   git_commit=m.provenance.get("versions", {}).get("git_commit"),
                   technology=sp["resolved"]["technology"]["name"],
                   blocklist=dict(version=bl.get("version"), md5=bl.get("md5"),
                                  note=_settings.blocklist_text(bl))),
        assign=dict(annotated=B.file_hash(annotated), steps=prof, alpha=alpha, alpha_note=alpha_note,
                    at_package_defaults=(not nd) if nd is not None else None, non_default_parameters=nd,
                    cell_cycle=_rs.cell_cycle_text(bl)),
        data=dict(cells=int(Z.n_obs), genes=int(Z.n_vars)), hash_note="input hash computed by the report",
        summary=str(Z.uns["scarlet_summary"]) if "scarlet_summary" in Z.uns else None))
    bd.table("colours", pd.DataFrame({"identity": list(col), "colour": list(col.values())}),
             "one fixed colour per identity; refusals grey " + B.REFUSED_COLOUR)
    g = genes_per_cell(Z.X)
    C = pd.DataFrame({"cell": Z.obs_names, "best": best, "unique": uniq, "reason_best": rb, "reason_unique": ru,
                      "genes_detected": g, "total_counts": total_counts(Z.X)})
    for c, n in (("scarlet_flip", "flip"), ("scarlet_pe_low", "pe_low")):
        if c in obs:
            C[n] = obs[c].to_numpy()
    bd.table("cells", C, f"per cell: calls, refusal reasons ({R.LABEL}), depth, flags")
    idx, xy, urec = umap_spaces(Z, umap_cells, seed, universe_genes([m]), "the model",
                                gene_col=m.settings["parameters"].get("gene_col", "index"), all_genes=umap_all_genes,
                                log=log)
    bd.table("umap", umap_frame(Z.obs_names[idx], xy).assign(best=best[idx], unique=uniq[idx]),
             "UMAP (display only; x, y: the model's universe genes; x_all, y_all: all genes) with the calls")
    bd.meta("umap", urec)
    rows = []
    for e, v, r in (("best", best, rb), ("unique", uniq, ru)):
        for k in ids:
            mk = v == k
            row = dict(entry=e, identity=k, called=int(mk.sum()), share=float(mk.mean()))
            if "flip" in C:
                row["flip_flagged"] = int(R._bool(C.flip.to_numpy()[mk]).sum())
            if "pe_low" in C:
                row["pe_low_flagged"] = int(R._bool(C.pe_low.to_numpy()[mk]).sum())
            rows.append(row)
        for why in sorted(set(r) - {R.CALLED}):
            rows.append(dict(entry=e, identity=f"refused: {why}", called=int((r == why).sum()),
                             share=float((r == why).mean())))
    bd.table("analytics", pd.DataFrame(rows), "calls per identity and refusals by reason, best and unique")
    bd.table("no_calls", pd.DataFrame([dict(entry=e, identity=k) for e, v in (("best", best), ("unique", uniq))
                                       for k in ids if not (v == k).any()], columns=["entry", "identity"]),
             "identities of the model with no call")
    mkr = m.markers
    u = m.universe
    u.index = u.index.astype(str)
    gz = pd.Index(Z.var_names.astype(str))
    SC = []
    for k in ids:
        gl = [x for x in mkr[mkr.cell_type.astype(str) == k].sort_values("idx").gene.astype(str) if x in u.index]
        called = best == k
        gi = gz[~gz.duplicated()].get_indexer(gl)   # a repeated name: its first column
        gi = np.where(gi >= 0, np.flatnonzero(~gz.duplicated())[np.maximum(gi, 0)], -1)
        det = np.full(len(gl), np.nan)
        if called.any() and (gi >= 0).any():
            det[gi >= 0] = detection_in(Z.X, called, gi[gi >= 0])
        for x, a, b in zip(gl, u.loc[gl, k].to_numpy(), det):
            SC.append(dict(identity=k, gene=x, detection_model=float(a), detection_called=b,
                           cells_called=int(called.sum())))
    bd.table("marker_detection", pd.DataFrame(SC), "per marker: the model's detection rate and the share of cells "
                                                   "called that identity (best) detecting it (computed by the report)")
    bd.table("markers", mkr, "the model's markers")
    bd.table("universe", u.reset_index().rename(columns={u.index.name or "index": "gene"}), "the model's universe")
    bd.meta("identities", ids)
    bd.save()
    return bd
