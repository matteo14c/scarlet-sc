"""One self-contained HTML per run, rendered from the bundle only: Plotly (plotly.js embedded once; no server,
no network at view time). Each page opens on a short summary and ends with the provenance; its UMAPs have their own
switches (output, view, highlight), and the identity drill-down has its own dropdown."""
import html as _h
import json
import os

import numpy as np
import pandas as pd

from scarlet import run_summary as _rs
from scarlet.html_report import bundle as B
from scarlet.html_report import reasons as R

CSS = """
:root{--fg:#1b1b1b;--muted:#5f6368;--line:#e0e0e0;--bg:#ffffff;--card:#fafafa;--accent:#1a5fb4}
@media (prefers-color-scheme: dark){:root{--fg:#e8e8e8;--muted:#a0a0a0;--line:#3a3a3a;--bg:#161616;--card:#1f1f1f;
--accent:#78aeed}}
body{font:14px/1.45 system-ui,-apple-system,Segoe UI,Roboto,sans-serif;color:var(--fg);background:var(--bg);
margin:0;padding:0 16px 48px;max-width:1280px;margin-inline:auto}
h1{font-size:22px;margin:24px 0 4px}h2{font-size:17px;margin:28px 0 8px;border-bottom:1px solid var(--line);
padding-bottom:4px}p.note,.note{color:var(--muted);font-size:12.5px}
table{border-collapse:collapse;font-size:12.5px;margin:6px 0}td,th{border-bottom:1px solid var(--line);
padding:3px 8px;text-align:left;vertical-align:top}th{font-weight:600}th[title]{cursor:help}
.wrap{overflow-x:auto;max-width:100%;border-right:1px solid var(--line)}
.card{background:var(--card);border:1px solid var(--line);border-radius:6px;padding:8px 12px}
.grid3{display:grid;grid-template-columns:repeat(auto-fit,minmax(320px,1fr));gap:8px}
.grid2{display:grid;grid-template-columns:repeat(auto-fit,minmax(420px,1fr));gap:12px}
.scrollbox{overflow:auto;max-height:660px;max-width:100%;border:1px solid var(--line);border-radius:4px}
ul.summary{margin:4px 0 8px 18px;padding:0}ul.summary li{margin:2px 0}
select{font-size:14px;padding:3px 6px}.sw{display:inline-block;width:10px;height:10px;border-radius:2px;
margin-right:5px;vertical-align:middle}.flag{font-weight:600}.footnote{margin-top:32px}
details{margin:8px 0}details>summary{cursor:pointer;font-weight:600;color:var(--accent);margin:4px 0}
.tiles{display:flex;flex-wrap:wrap;gap:8px;margin:6px 0}.tile{background:var(--card);border:1px solid var(--line);
border-radius:6px;padding:6px 10px;min-width:150px}.tile b{display:block;font-size:17px}
"""


def _esc(x):
    return _h.escape("" if x is None or (isinstance(x, float) and np.isnan(x)) else str(x))


#: readable table headers (the column's own name, as in the bundle's TSV, is the header's tooltip)
HEADERS = {
    "entry": "output", "cells": "cells", "refused_flip test": "refused: flip test", "best": "best call",
    "unique": "unique call", "forced": "forced (top identity)", "runner_up": "runner-up", "truth": "reference label",
    "fdr": "FDR", "reason_best": "reason (best)", "reason_unique": "reason (unique)", "flip": "flip test",
    "pe_low": "low positive evidence", "pct_best_correct": "best correct", "pct_unique_called": "with a unique call",
    "pct_boundary": "boundary", "median_fdr_best": "median FDR (best calls)",
    "median_margin_best": "median margin (best calls)", "under_50_cells": "under 50 cells",
    "median_genes_detected": "median genes detected", "n_eff": "effective cells", "pe_expected": "expected markers",
    "pe_required": "required markers", "pe_typical_median": "typical markers detected",
    "pe_median_detected_cv": "median markers detected (CV)", "pe_low_share_cv": "low evidence (CV)",
    "cv_recall_forced": "recall forced (CV)", "cv_recall_best": "recall best (CV)",
    "cv_recall_unique": "recall unique (CV)", "most_confused_share": "share confused",
    "refused_best_share": "refused best (CV)", "agreement": "agreement (draws)", "anchors": "stable anchors",
    "first_pass_best": "first pass best", "first_pass_unique": "first pass unique", "refined_best": "refined best",
    "refined_unique": "refined unique", "first_pass_best_pct": "first pass best (%)", "refined_best_pct": "refined best (%)",
    "first_pass_unique_pct": "first pass unique (%)", "refined_unique_pct": "refined unique (%)",
    "changed_by_refinement": "changed by refinement", "refused_after": "refused after refinement",
    "called_best": "called best", "called_unique": "called unique", "stable_anchors": "stable anchors",
    "pct_full_agreement": "at full agreement", "ranked_first": "ranked first", "n_cells_transferred": "cells called",
    "in_query_programme": "in the query's programme", "reference_markers": "reference markers",
    "query_markers": "query markers", "conserved": "conserved markers", "share_conserved": "share conserved",
    "doubtful_transfer": "doubtful transfer", "best_agreement": "agreement best", "unique_agreement": "agreement unique",
    "top1": "top identity", "runner_up_draws": "runner-up (draws)", "first_best": "first pass best",
    "first_unique": "first pass unique", "first_best_agreement": "first pass agreement best",
    "first_unique_agreement": "first pass agreement unique", "first_reason_best": "first pass reason (best)",
    "first_reason_unique": "first pass reason (unique)", "first_top1": "first pass top identity",
    "first_runner_up": "first pass runner-up", "first_runner_up_draws": "first pass runner-up (draws)",
    "stable_anchor": "stable anchor of", "anchor": "stable anchor", "first_flip": "first pass flip test",
    "refine_note": "refinement note", "first_pass": "first pass", "flip_flagged": "flip test flagged",
    "pe_low_flagged": "low evidence flagged", "first_best_labelled": "first pass best called",
    "first_unique_labelled": "first pass unique called", "best_labelled": "final best called",
    "unique_labelled": "final unique called", "depth_ref": "depth A", "depth_query": "depth B",
    "depth_after": "depth after", "shared_genes": "shared genes", "median_f1_called": "median F1 on called cells",
    "pct_refused": "refused", "cells_label_unknown_to_model": "cells of labels unknown to the model",
    "pct_refused_label_unknown_to_model": "of those, refused", "known_to_model": "known to the model",
    "label": "label (supplied)",
}
#: tooltips beyond the column's name
NOTES = {
    "fdr": "Benjamini-Hochberg FDR of the call against the background",
    "margin": "log-likelihood of the top identity minus the larger of the runner-up's and the background's",
    "runner_up": "for a called cell, the identity ranked second; for a cell refused for background, its top identity, "
                 "which the background beat (its margin is then negative)",
}
#: columns holding a share between 0 and 1 (shown as percent)
SHARES = {"recall", "share", "most_confused_share", "refused_best_share", "pe_low_share_cv", "share_conserved",
          "cv_recall_forced", "cv_recall_best", "cv_recall_unique", "keep_prob_median_shared",
          "keep_prob_mean_detection"}
#: integer columns written without a thousands separator
PLAIN_INT = {"seed", "draw", "fold", "rank", "idx", "agreement", "draws", "best_agreement", "unique_agreement",
             "first_best_agreement", "first_unique_agreement", "runner_up_draws", "first_runner_up_draws"}


def header(c):
    c = str(c)
    if c in HEADERS:
        return HEADERS[c]
    h = c[4:] if c.startswith("pct_") else c
    return h.replace("_", " ")


def cell_text(c, v):
    """A table value as the reader wants it: percent with one decimal, counts as integers, FDR in scientific
    notation, other numbers to three significant digits; flags as yes / no."""
    c = str(c)
    if v is None or (isinstance(v, float) and np.isnan(v)):
        return ""
    if isinstance(v, (bool, np.bool_)) or (isinstance(v, str) and v in ("True", "False")):
        return "yes" if str(v) == "True" else "no"
    if isinstance(v, (int, np.integer)) or (isinstance(v, (float, np.floating)) and float(v).is_integer()
                                             and c not in SHARES and "fdr" not in c):
        if c.startswith("pct_") or c.endswith("_pct"):
            return f"{float(v):.1f}%"
        return f"{int(v)}" if c in PLAIN_INT else f"{int(v):,}"
    if isinstance(v, (float, np.floating)):
        v = float(v)
        if c.startswith("pct_") or c.endswith("_pct"):
            return f"{v:.1f}%"
        if c in SHARES:
            return f"{100 * v:.1f}%"
        if "fdr" in c:
            return f"{v:.1e}" if v else "0"
        if "margin" in c or c.startswith("ll_") or c.startswith("depth"):
            return f"{v:,.1f}"
        return f"{v:.3g}"
    return _esc(v)


def table_html(df, cols=None, fmt=None, max_rows=None, notes=None, drop_zero=True):
    """A table with readable headers (the column's name in the tooltip) and readable values (cell_text); wide tables
    scroll sideways. drop_zero: columns of refusal shares or counts that are zero in every row are left out."""
    cols = list(cols or df.columns)
    if drop_zero and len(df):
        def zero(c):
            if not (str(c).startswith(("pct_unlabelled", "refused_")) or str(c).endswith("_flagged")):
                return False
            v = pd.to_numeric(df[c], errors="coerce")
            return bool(((v.fillna(0)) == 0).all())
        cols = [c for c in cols if not zero(c)]
    fmt = fmt or {}
    nt = dict(NOTES, **(notes or {}))
    d = df if max_rows is None else df.head(max_rows)
    head = "".join(f'<th title="{_esc(c)}{(": " + _esc(nt[c])) if c in nt else ""}">{_esc(header(c))}</th>'
                   for c in cols)
    body = "".join("<tr>" + "".join(f"<td>{fmt[c](r[c]) if c in fmt else cell_text(c, r[c])}</td>" for c in cols)
                   + "</tr>" for _, r in d.iterrows())
    return f'<div class="wrap"><table><thead><tr>{head}</tr></thead><tbody>{body}</tbody></table></div>'


def _kv(d, prefix=""):
    rows = []
    for k, v in d.items():
        if isinstance(v, dict):
            rows += _kv(v, f"{prefix}{k}.")
        else:
            rows.append((f"{prefix}{k}", json.dumps(v) if isinstance(v, (list, tuple)) else v))
    return rows


def _plotlyjs():
    from plotly.offline import get_plotlyjs
    return get_plotlyjs()


def point_size(n):
    """A UMAP's marker size from the number of cells shown: larger with fewer cells (8 at 100 cells or fewer,
    about 5 at 1,000, 3 at 10,000, 2 from 30,000)."""
    return round(float(np.clip(12 - 2.2 * np.log10(max(int(n), 1)), 2, 8)), 1)


def _layout_title(text):
    """A plot title at the top left, clear of the plot's toolbar (top right, shown on hover)."""
    return dict(text=text, x=0.01, xanchor="left", y=0.98, yanchor="top", font=dict(size=13))


def _depth_fig(D, ids, col):
    import plotly.graph_objects as go
    fig = go.Figure()
    for k in ids:
        fig.add_trace(go.Box(x=D.genes_detected[D.identity == k], name=k, marker_color=col[k], boxpoints=False,
                             orientation="h"))
    fig.update_layout(height=max(300, 22 * len(ids) + 80), showlegend=False, xaxis_title="genes detected per cell",
                      margin=dict(l=10, r=10, t=10, b=40), template="plotly_white",
                      yaxis=dict(autorange="reversed"))
    return fig


def _analytics_fig(AN, ids):
    import plotly.graph_objects as go
    from plotly.subplots import make_subplots
    ents = ["forced", "best", "unique"]
    fig = make_subplots(rows=1, cols=3, shared_yaxes=True, subplot_titles=ents, horizontal_spacing=0.03)
    parts = [("correct", "correct", "#4CAF50"), ("wrong", "wrong", "#E53935"),
             (f"refused_{R.BACKGROUND}", "refused: background", "#9E9E9E"),
             (f"refused_{R.AMBIGUOUS}", "refused: ambiguity", "#BDBDBD"),
             (f"refused_{R.FLIP}", "refused: flip test", "#616161")]
    for j, e in enumerate(ents, 1):
        a = AN[AN.entry == e].set_index("identity").reindex(ids)
        n = a.cells.replace(0, np.nan)
        for c, lab, cl in parts:
            fig.add_trace(go.Bar(y=ids, x=(a[c] / n).fillna(0), orientation="h", name=lab, marker_color=cl,
                                 legendgroup=lab, showlegend=j == 1, customdata=a[c],
                                 hovertemplate="%{y}<br>" + lab + ": %{customdata} cells (%{x:.1%})<extra></extra>"),
                          row=1, col=j)
    fig.update_layout(barmode="stack", height=max(320, 22 * len(ids) + 100), template="plotly_white",
                      margin=dict(l=10, r=10, t=30, b=30), legend=dict(orientation="h", y=-0.05, font=dict(size=11)))
    fig.update_yaxes(autorange="reversed")
    fig.update_xaxes(range=[0, 1], tickformat=".0%")
    return fig


def _heat(div, title, genes, cols, z, groups=None):
    return dict(div=div, title=title, genes=genes, cols=cols, z=z, groups=groups or [])


def _by_detection(uni, genes, k):
    """genes ordered by their detection rate in identity k (highest first; the marker rank breaks ties)"""
    if k not in uni.columns:
        return list(genes)
    d = pd.to_numeric(uni[k].reindex(genes), errors="coerce").to_numpy(float)
    order = sorted(range(len(genes)), key=lambda i: (-(d[i] if np.isfinite(d[i]) else -1.0), i))
    return [genes[i] for i in order]


def _drill_heat(mk, uni, ids, k, div, title, genes=None, groups=None, inc=None, split=False):
    """the identity's markers (or `genes`, in that order) x every identity of the universe (+ bg); markers ordered
    by their detection in k (item order otherwise: rank). inc: the model's in-class table (every value read from it;
    its marker rows are the universe's); without it (a model written before 0.6.3) only the model's own markers
    have values. gaps: per cell, why it is blank (1 not a marker of this model, 2 not among the genes the model's
    marker search scored). split: the groups drawn as separate blocks."""
    if genes is None:
        g = [x for x in mk[mk.cell_type.astype(str) == k].sort_values("idx").gene.astype(str).tolist()
             if x in uni.index]
        g = _by_detection(uni, g, k)
    else:
        g = list(genes)
    cols = [c for c in ids if c in uni.columns and uni[c].notna().any()] + (["bg"] if "bg" in uni.columns else [])
    z = (inc if inc is not None else uni).reindex(index=g, columns=cols).round(4)
    h = _heat(div, title, g, cols, None, groups)
    if split:
        gap = 2 if inc is not None else 1
        h["gaps"] = np.where(z.isna().to_numpy(), gap, 0).tolist() if z.isna().any().any() else []
    h["z"] = z.astype(object).where(z.notna(), None).to_numpy().tolist()
    h["split"] = bool(split)
    return h


def _rank_plot(div, k, mA, mB):
    """each marker's rank in A's model and in B's refined model (by score, highest first): common markers joined by
    a line, A-only markers on A's side, B-only on B's"""
    def ranks(m):
        s = m[m.cell_type.astype(str) == k].reset_index(drop=True)
        s = s.iloc[np.lexsort((np.arange(len(s)), -pd.to_numeric(s.score, errors="coerce").to_numpy()))]
        return {g: i + 1 for i, g in enumerate(s.gene.astype(str))}
    a, b = ranks(mA), ranks(mB)
    return dict(div=div, title=f"{k}: each marker's rank in A's model and in B's refined model (1 = highest score)",
                common=[[g, a[g], b[g]] for g in a if g in b], a_only=[[g, a[g]] for g in a if g not in b],
                b_only=[[g, b[g]] for g in b if g not in a])


def _wrong_html(wrong, k):
    if wrong is None:
        return ""
    w = wrong[wrong.identity.astype(str) == k]
    if not len(w):
        return "<span class='note'>No wrong calls (out-of-fold).</span>"
    return ("<b>Wrong calls of " + _esc(k) + " went to</b> (out-of-fold):" +
            table_html(w, ["entry", "called", "cells", "share"]))


def _prov_table(rows):
    return "<table>" + "".join(f"<tr><th>{_esc(a)}</th><td>{_esc(b)}</td></tr>" for a, b in _rows(rows)) + "</table>"


def _bl_text(blk):
    s = f"{blk.get('note')}"   # how the list was chosen, in plain words (its version stays in provenance)
    if blk.get("modules_on") is not None:
        s += f"; modules on: {', '.join(blk.get('modules_on') or [])}; off: {', '.join(blk.get('modules_off') or [])}"
    return s


def _summary_block(prov):
    """The run's plain-language summary (summary.txt), as the run wrote it."""
    t = prov.get("summary")
    if not t:
        return "<p><em>This run has no summary.</em></p>"
    return f'<h3>What this run did</h3><pre style="white-space:pre-wrap">{_esc(t)}</pre>'


def _hash_text(f):
    f = f or {}
    text = f.get("sha256") or f.get("note") or ""
    return text + (f"  ({f['check']})" if f.get("check") else "")


def _rows(head):
    """Provenance rows; an empty git commit (an installed package has none) is left out."""
    return [(k, v) for k, v in head if not (k == "git commit" and not v)]


def _num_text(v):
    """a margin for the page: three significant digits; '-' when there is none"""
    return "-" if v is None or not np.isfinite(v) else f"{v:.3g}"


def _fdr_text(lg):
    """an FDR given as -log10: 1e-3.9"""
    return "-" if lg is None else f"1e-{lg:.1f}"


def _pct(n, d):
    return f"{100 * n / d:.1f}%" if d else "-"


def _ul(items):
    return "<ul class='summary'>" + "".join(f"<li>{x}</li>" for x in items if x) + "</ul>"


def _reasons_text(reasons, n):
    """'background 2.1%, ambiguity 1.6%' for the refusal reasons present"""
    r = pd.Series(reasons).astype(str)
    r = r[(r != "") & (r != "nan")]
    return ", ".join(f"{_esc(k)} {_pct(v, n)}" for k, v in r.value_counts().items())


def _footnote(n_ids):
    return (f"<p class='note footnote'>Colours: one fixed colour per identity in every panel, refusals grey "
            f"({_esc(B.palette_note(n_ids))}).</p>")


def _label_spec(div, title, U, outputs, col, truth=None, spaces=None):
    """A switchable UMAP whose outputs are labels only (no certainty, reasons or anchors)."""
    from scarlet.html_report import certainty as C_
    outs = [dict(name=n, label=lab, reason=None, runner_up=None, hover=[])
            for n, lab in outputs]
    return C_.umap_spec(div, title, U, outs, col, truth=truth, views=["label"], spaces=spaces)


def _space_names(urec):
    """The UMAP's gene spaces as the report recorded them (an older bundle: its one method)."""
    sp = (urec or {}).get("spaces")
    return [s["name"] for s in sp] if sp else [(urec or {}).get("method", "scanpy defaults on all genes")]


def _space_note(urec):
    """The caption on the UMAPs' gene spaces."""
    names = _space_names(urec)
    s = (f"Space: {_esc(names[0])}, the genes SCARLET decides on, by default"
         + (f"; or {_esc(names[1])}, the place to look for populations the model lacks (on the model's genes such "
            "cells sit with their nearest relative, as SCARLET calls them; on all genes they may stand apart)"
            if len(names) > 1 else "") + ".")
    if (urec or {}).get("all_genes_skipped"):
        s += " The all-genes UMAP was not computed (umap_all_genes = no)."
    return s


def _cell_table(df, key, file, groups, group_names, default, notes=None, page=None):
    """A per-cell table drawn in the page: a 'show' filter (groups, the cells worth a look first), headers that sort
    the rows, the first 200 matching rows (or `page` rows, with a button for more); the button downloads every row
    (the payload's columns, unchanged)."""
    from scarlet.html_report import certainty as C_
    T = C_.table_payload(df, file)
    nt = dict(NOTES, **(notes or {}))
    kinds = []
    for c in df.columns:
        v = df[c]
        if pd.api.types.is_bool_dtype(v) or set(v.dropna().astype(str)) <= {"True", "False"} and len(v.dropna()):
            kinds.append("bool")
        elif "fdr" in c:
            kinds.append("sci")
        elif "margin" in c or c.startswith("ll_"):
            kinds.append("f1")
        elif pd.api.types.is_integer_dtype(v):
            kinds.append("int")
        elif pd.api.types.is_float_dtype(v):
            kinds.append("g3")
        else:
            kinds.append("text")
    T.update(headers=[header(c) for c in df.columns], tips=[f"{c}" + (f": {nt[c]}" if c in nt else "")
                                                             for c in df.columns],
             kinds=kinds, groups=[int(g) for g in groups], group_names=group_names, default=default)
    if page:
        T["page"] = int(page)
    opts = "".join(f'<option value="{i}">{_esc(n)}</option>' for i, n in enumerate(group_names))
    html = (f'<p>Show: <select id="{key}_show">{opts}</select> <span class="note" id="{key}_count"></span> '
            '<span class="note">(click a header to sort)</span></p><div class="wrap" id="' + key + '_table"></div>'
            + (f'<p><button type="button" id="{key}_more">Show {int(page)} more</button></p>' if page else ""))
    return T, html


def _calls_fig(AN, ids, col, title):
    """Calls per identity, best and unique side by side, in each identity's colour (unique hatched), with the cells
    left unassigned as a grey bar; the count and its share of the cells in the hover (no numbers above the bars)."""
    import plotly.graph_objects as go
    fig = go.Figure()
    called = AN[~AN.identity.astype(str).str.startswith("refused")]
    refused = AN[AN.identity.astype(str).str.startswith("refused")]
    xs = list(ids) + [R.U]
    for e, pattern in (("best", ""), ("unique", "/")):
        a = called[called.entry == e].set_index(called[called.entry == e].identity.astype(str))
        r = refused[refused.entry == e]
        n = [int(a.called.get(k, 0)) for k in ids] + [int(r.called.sum())]
        sh = [float(a.share.get(k, 0.0)) for k in ids] + [float(r.share.sum())]
        fig.add_trace(go.Bar(x=xs, y=n, name=e, customdata=[f"{c:,} ({100 * s:.1f}%)" for c, s in zip(n, sh)],
                             marker=dict(color=[col.get(k, B.REFUSED_COLOUR) for k in ids] + [B.REFUSED_COLOUR],
                                         pattern=dict(shape=pattern, fillmode="overlay", fgcolor="#ffffff")),
                             hovertemplate="%{x}<br>%{customdata} of the cells<extra>" + e + "</extra>"))
    fig.update_layout(barmode="group", title=_layout_title(title + " (best solid, unique hatched; hover for the "
                                                                   "count and share of the cells)"),
                      height=460, template="plotly_white", showlegend=False,
                      margin=dict(l=10, r=10, t=60, b=140), xaxis=dict(tickangle=-45))
    return fig


def _id_fmt(col):
    return {"identity": lambda v: f'<span class="sw" style="background:{col.get(str(v), B.REFUSED_COLOUR)}"></span>'
                                  f'{_esc(v)}',
            "under_50_cells": lambda v: '<span class="flag">yes</span>' if str(v) == "True" else ""}


def _picker(ids):
    opts = "".join(f'<option value="{_esc(k)}">{_esc(k)}</option>' for k in ids)
    return f'<p>Identity: <select id="pick">{opts}</select> <span class="note">(drives this section)</span></p>'


def _data(figs, drill):
    return dict(figs={k: json.loads(f.to_json()) for k, f in figs.items()}, drill=drill)


def _scatter(div, title, x, y, text, xt, yt):
    return dict(div=div, title=title, x=x, y=y, text=text, xt=xt, yt=yt)


BG_NOTE = ("bg: the background, the gene's detection rate across the reference's cells (with unequal class sizes "
           "balanced, the mean over the cell types); a cell whose best cell type explains it no better than the "
           "background is refused (reason background)")


# --------------------------------------------------------------------------------------------------------------
# learn report
# --------------------------------------------------------------------------------------------------------------
def _learn_summary(T, meta, ids, sa):
    items = [f"{len(ids)} identities; {meta['umap']['cells_total']:,} reference cells."]
    if sa.get("off") or "self_assessment" not in T:
        return _ul(items + ["Self-assessment off: no cross-validated numbers."])
    SA = T["self_assessment"]
    n = len(SA)
    AN = T["analytics"].assign(identity=T["analytics"].identity.astype(str))
    for e in ("best", "unique"):
        a = AN[AN.entry == e].set_index("identity").recall.dropna()
        if len(a):
            items.append(f"Cross-validated recall, {e} call: {100 * a.min():.1f}% ({_esc(a.idxmin())}) to "
                         f"{100 * a.max():.1f}% ({_esc(a.idxmax())}) over the identities.")
    for e in ("best", "unique"):
        ref = SA[e].astype(str) == R.U
        items.append(f"Refused, {e} call: {_pct(int(ref.sum()), n)} of the cells "
                     f"({_reasons_text(SA[f'reason_{e}'][ref], n) or 'none'}).")
    W = T.get("weak_classes")
    if W is not None and len(W):
        w = W.sort_values("rank").head(3)
        items.append("Weakest identities: " + "; ".join(
            f"{_esc(r.identity)} (recall best {100 * float(r.cv_recall_best):.1f}%)" if "cv_recall_best" in W
            and pd.notna(r.cv_recall_best) else _esc(r.identity) for r in w.itertuples()) + ".")
    return _ul(items)


def render_learn(out, back_link=None):
    man, T = B.load(out)
    meta = man["meta"]
    ids = meta["identities"]
    col = dict(zip(T["colours"].identity.astype(str), T["colours"].colour))
    prov = meta["provenance"]
    sa = meta.get("self_assessment", {})
    figs = {"depth": _depth_fig(T["depth"].assign(identity=T["depth"].identity.astype(str)), ids, col)}
    U = T["umap"]
    outs = [("reference labels", U["label"])]
    if "cv_best" in U:
        outs += [("cross-validated best", U["cv_best"]), ("cross-validated unique", U["cv_unique"])]
        if "cv_forced" in U:
            outs.append(("cross-validated forced (top identity)", U["cv_forced"]))
    spn = _space_names(meta["umap"])
    specs = [_label_spec("umap2", "Reference", U, outs, col, truth=U["label"], spaces=spn)]
    if "analytics" in T:
        figs["analytics"] = _analytics_fig(T["analytics"].assign(identity=T["analytics"].identity.astype(str)), ids)
    mk = T["markers"]
    uni = T["universe"].set_index("gene")
    uni.index = uni.index.astype(str)
    wrong = T.get("wrong_to")
    drill = {k: dict(heat=[_drill_heat(mk, uni, ids, k, "heat", f"{k}: its markers (highest detection in {k} "
                                                                f"first) x every identity, detection rate")],
                     html={"drillhtml": _wrong_html(wrong, k)}) for k in ids}
    pm = prov["model"]
    sz = pm.get("sizes") or {}
    head = [("SCARLET (model)", pm.get("scarlet")), ("git commit", pm.get("git_commit")),
            ("learned", pm.get("created_at")), ("technology preset", pm.get("technology")),
            ("markers / shortlist", f"{sz.get('n_markers')} / {sz.get('marker_shortlist')} ({sz.get('from')})"),
            ("seed", pm.get("seed")), ("blocklist", _bl_text(pm["blocklist"])),
            ("unequal class sizes", pm.get("imbalanced")),
            ("cell-cycle genes", pm.get("cell_cycle")),
            ("learn at package defaults", _rs.options_report_text(pm.get("non_default_parameters"))),
            ("reference", (prov["input"].get("reference") or {}).get("path")),
            ("reference sha256", _hash_text(prov["input"].get("reference"))),
            ("cells x genes", f"{prov['input']['cells']} x {prov['input']['genes']}"),
            ("labels", f"{prov['input']['labels']} ({prov['input']['identities']} identities)")]
    if prov.get("note"):
        head.insert(0, ("this report", prov["note"]))
    if not sa.get("off"):
        f = sa["folds"]
        head.append(("self-assessment", f"{f['n_folds']} folds, {f['source']}" +
                     (f", seed {f['seed']}" if f.get("seed") is not None else "") +
                     f"; assign alpha {sa['alpha']}, flip {sa['flip_mode']}; fold "
                     "models learned with the model's recorded parameters"))
    else:
        head.append(("self-assessment", "off"))
    head += [("UMAP", f"{meta['umap']['cells_shown']} of {meta['umap']['cells_total']} cells; {meta['umap']['method']}"),
             ("report", f"SCARLET {prov['report']['scarlet']}, {prov['report'].get('written_at', '')}")]
    from scarlet.html_report import certainty as C_
    cert_html, cells_html, tables = "", "", {}
    if not sa.get("off") and "self_assessment" in T and "margin" in T["self_assessment"]:   # certainty and boundary views
        SA = T["self_assessment"].set_index("cell")
        UU = U.set_index("cell").join(SA[[c for c in SA.columns if c not in U.columns]], how="left").reset_index()
        outs3 = []
        mgn, m_lo, m_hi = C_.intensity_from_margin(UU.margin, UU["best"].astype(str) != R.U)
        amb = C_.intensity_ambiguous(UU["best"], UU["unique"])
        for name, e, rcol in (("cross-validated best", "best", "reason_best"),
                              ("cross-validated unique", "unique", "reason_unique"),
                              ("cross-validated forced (top identity)", "forced", None)):
            outs3.append(dict(name=name, label=UU[e], shade={"margin": mgn, "ambiguity": amb},
                              reason=UU[rcol] if rcol else None, runner_up=UU["runner_up"],
                              hover=[("FDR", pd.to_numeric(UU.fdr, errors="coerce")),
                                     ("margin", pd.to_numeric(UU.margin, errors="coerce"))]))
        specs.append(C_.umap_spec("qumap", "Reference (out-of-fold calls)", UU, outs3, col, truth=UU["label"],
                                  views=["label", "margin", "ambiguity", "unlabelled"], spaces=spn))
        BRl = T.get("boundary", pd.DataFrame(columns=["stage", "best", "runner_up", "cells"]))
        nbl = int(BRl.cells.sum()) if len(BRl) else 0
        f = C_.boundary_fig(BRl, "cross-validation", f"Boundary cells (out-of-fold): best call x runner-up ({nbl} cells)",
                            ids)
        if f is not None:
            figs["bnd_cv"] = f
        CPI = T.get("certainty_per_identity")
        cert_html = (f"<h2>3b. Certainty and boundary cells (cross-validation)</h2>"
                     f"<p class='note'>Out-of-fold. A boundary cell has a best call but no unique call (its best "
                     f"identity is not separated from the runner-up at alpha {sa['alpha']}); {nbl} cells. Views of the "
                     f"UMAP: <i>margin</i>, the identity's colour, darker the larger the call's margin, relative to this "
                     f"run: faintest at a margin of {_num_text(m_lo)} or below, full at {_num_text(m_hi)} or above (the "
                     f"5th and 95th percentiles of this run's best calls; the values are in the hover and the per-cell "
                     f"table); <i>ambiguity</i>, the boundary cells in full colour, the others faint (seen on the best "
                     f"call; in the unique call they are the unassigned ones). "
                     f"Margin = ll_top1 &minus; max(ll_top2, ll_bg).</p>" +
                     C_.umap_controls(specs[-1]) +
                     ('<div id="bnd_cv"></div>' if "bnd_cv" in figs else "<p class='note'>no boundary cell</p>") +
                     ("<p><b>Per reference label</b> (out-of-fold): % best correct, % with a unique call, % "
                      "boundary, % unlabelled by reason, median FDR and margin of its best calls.</p>" +
                      table_html(CPI, fmt=_id_fmt(col)) if CPI is not None else ""))
        keep = [c for c in ["cell", "fold", "truth", "forced", "best", "unique", "fdr", "margin", "runner_up",
                            "reason_best", "reason_unique", "flip", "pe_low"] if c in T["self_assessment"]]
        SAt = T["self_assessment"][keep]
        tr, b, u = (SAt[c].astype(str).to_numpy() for c in ("truth", "best", "unique"))
        grp = np.where((b != R.U) & (b != tr), 1, np.where(b == R.U, 2, np.where(u == R.U, 3, 4)))
        order = np.argsort(grp, kind="stable")
        SAt = SAt.iloc[order].reset_index(drop=True)
        tables["cells"], tab = _cell_table(
            SAt, "cells", "scarlet_learn_selfassessment_cells.tsv", grp[order],
            ["worth a look: wrong, refused or boundary", "wrong (best call not the label)", "refused (no best call)",
             "boundary (best call, no unique call)", "correct", "every cell"], 0)
        cells_html = ("<h2>6. Per-cell data (cross-validation)</h2><p class='note'>Every reference cell: its fold, "
                      "label, out-of-fold calls, FDR, margin, runner-up, reasons and flags; the cells worth a look "
                      "first. The runner-up of a cell refused for background is its top identity, the one the "
                      "background beat. The button downloads every cell.</p>" +
                      C_.download_button("cells", len(SAt), "cells") + tab)
    W = T["weak_classes"]
    wcols = [c for c in ["rank", "identity", "cv_recall_best", "cv_recall_unique", "cv_recall_forced",
                         "most_confused_with", "most_confused_share", "refused_best_share", "cells", "under_50_cells",
                         "median_genes_detected", "n_eff", "markers", "pe_expected", "pe_required",
                         "pe_median_detected_cv", "pe_low_share_cv"] if c in W]
    sec5 = table_html(W, wcols, fmt=_id_fmt(col))
    sa_note = ("" if sa.get("off") else
               f"<p class='note'>Out-of-fold: each cell is called by a model learned without it ({sa['cells']} cells). "
               f"forced = the top identity (scarlet_top1). Refusal reasons are {R.LABEL}.</p>")
    back = f"<p><a href='{_esc(back_link)}'>Back to the transfer report</a></p>" if back_link else ""
    body = f"""{back}
<h1>SCARLET learn report</h1>
<p class="note">Only values SCARLET computed are shown; the report reads SCARLET's outputs and changes no call. UMAPs are
display only.</p>
<h2>Summary</h2>{_learn_summary(T, meta, ids, sa)}
<h2>1. Depth: genes per cell per identity</h2><div id="depth"></div>
<h2>2. UMAP</h2><p class="note">Choose what colours the cells; highlight one identity if you wish (all in full colour
by default). {_space_note(meta["umap"])}</p>{C_.umap_controls(specs[0], highlight=True)}
<h2>3. Analytics per identity</h2>{sa_note}<div id="analytics"></div>
{cert_html}
<h2>4. Identity drill-down</h2>{_picker(ids)}
<p class="note">{_esc(BG_NOTE)}.</p><div class="scrollbox"><div id="heat"></div></div><div id="drillhtml" class="card"></div>
<h2>5. Weak-class flags</h2>
<p class="note">Ranked weakest first by {_esc(meta.get('weak_rank_key'))} (ties: fewer cells first). Evidence only:
no threshold is applied. <b>under 50 cells</b> marks identities with fewer than 50 reference cells, a reporting
convention, not a SCARLET threshold. The markers' columns are SCARLET's positive-evidence values (expected markers
E_k, required = floor(E_k / 5), the in-sample typical median; out-of-fold median detected and share flagged low).
Hover a header for the column's name in the bundle.</p>
{sec5}
{cells_html}
<h2>7. Provenance</h2>{_summary_block(prov)}{_prov_table(head)}
{_footnote(len(ids))}"""
    data = _data(figs, drill)
    data["qumaps"] = specs
    data["tables"] = tables
    return _page("SCARLET learn report", body, data)


# --------------------------------------------------------------------------------------------------------------
# transfer report
# --------------------------------------------------------------------------------------------------------------
def _names(prov):
    """The two datasets' names: their file stems (A and B for data passed in memory)."""
    out = []
    for side, dflt in (("reference", "A"), ("query", "B")):
        p = (prov[side].get("file") or {}).get("path") or ""
        out.append(os.path.splitext(os.path.basename(p))[0] if p and p != "<AnnData>" else dflt)
    return out


def _depth_line(DM, names):
    """'smartseq2_reference thinned from 2,357 to ~814 detected genes per cell (median over 10 draws), to match
    10xv2_query (845)', one clause per origin"""
    if DM is None or not len(DM):
        return "no depth matching"
    parts = []
    for g, d in DM.groupby("origin", sort=False):
        th = str(d.thinned.iloc[0])
        if th not in ("reference", "query"):
            parts.append(f"origin {_esc(g)}: not thinned")
            continue
        src, other = ("depth_ref", "depth_query") if th == "reference" else ("depth_query", "depth_ref")
        who, whom = (names[0], names[1]) if th == "reference" else (names[1], names[0])
        after = pd.to_numeric(d.depth_after, errors="coerce")
        parts.append(f"{_esc(who)}{'' if g == 'all' else f' (origin {_esc(g)})'} thinned from "
                     f"{float(d[src].iloc[0]):,.0f} to ~{after.median():,.0f} detected genes per cell (median over "
                     f"{len(d)} draws), to match {_esc(whom)} ({float(d[other].iloc[0]):,.0f})")
    return "; ".join(parts[:4]) + (f"; and {len(parts) - 4} more origins" if len(parts) > 4 else "")


def _transfer_depth_fig(D, DM, names):
    """Per cell before matching (boxes, every gene); the record's medians on the shared genes before and after
    matching (after: median over the draws, with their range)."""
    import plotly.graph_objects as go
    fd = go.Figure()
    sides = list(dict.fromkeys(D.side))
    for sd in sides:
        fd.add_trace(go.Box(x=D.genes_detected[D.side == sd], name=sd, boxpoints=False, orientation="h",
                            legendgroup="cells", showlegend=False))
    if DM is not None and len(DM):
        for g, d in DM.groupby("origin", sort=False):
            th = str(d.thinned.iloc[0])
            if th not in ("reference", "query"):
                continue
            side = next((s for s in sides if s.startswith("A" if th == "reference" else "B")), sides[0])
            src = "depth_ref" if th == "reference" else "depth_query"
            after = pd.to_numeric(d.depth_after, errors="coerce")
            fd.add_trace(go.Scatter(x=[float(d[src].iloc[0])], y=[side], mode="markers", name="before (record)",
                                    marker=dict(symbol="circle-open", size=12, color="#333333"),
                                    legendgroup="before", showlegend=g == DM.origin.iloc[0],
                                    hovertemplate=f"origin {g}: before matching, median %{{x:,.0f}} genes on the "
                                                  "shared genes<extra></extra>"))
            fd.add_trace(go.Scatter(x=[float(after.median())], y=[side], mode="markers", name="after (median over "
                                                                                             "draws, range)",
                                    marker=dict(symbol="diamond", size=12, color="#1a5fb4"),
                                    error_x=dict(type="data", symmetric=False,
                                                 array=[float(after.max() - after.median())],
                                                 arrayminus=[float(after.median() - after.min())]),
                                    legendgroup="after", showlegend=g == DM.origin.iloc[0],
                                    hovertemplate=f"origin {g}: after matching, median %{{x:,.0f}} (draws "
                                                  f"{after.min():,.0f} to {after.max():,.0f})<extra></extra>"))
    fd.update_layout(height=300, xaxis_title="genes detected per cell", template="plotly_white",
                     margin=dict(l=10, r=10, t=10, b=40), legend=dict(orientation="h", y=-0.3, font=dict(size=11)))
    return fd


def _agreement_figs(QC, ids, col, n_dr, k_need):
    """The agreement (k of n draws) of the final calls: overall (best and unique), and by final best call."""
    import plotly.graph_objects as go
    ks = list(range(1, int(n_dr) + 1))
    fo = go.Figure()
    for e, pattern in (("best", ""), ("unique", "/")):
        a = pd.to_numeric(QC[f"{e}_agreement"], errors="coerce").dropna().astype(int)
        cnt = [int((a == k).sum()) for k in ks]
        fo.add_trace(go.Bar(x=ks, y=cnt, name=e, marker=dict(color=["#1a5fb4" if k >= k_need else "#9E9E9E"
                                                                    for k in ks],
                                                             pattern=dict(shape=pattern, fgcolor="#ffffff")),
                            hovertemplate="%{x} of " + str(n_dr) + " draws: %{y} cells<extra>" + e + "</extra>"))
    fo.update_layout(barmode="group", height=320, template="plotly_white", showlegend=False,
                     title=_layout_title(f"Agreement of the final calls (best solid, unique hatched); blue: at least "
                                         f"{k_need} of {n_dr}, the final call; grey: below, unassigned"),
                     xaxis=dict(title="draws giving the cell's most frequent outcome", dtick=1),
                     yaxis_title="cells", margin=dict(l=10, r=10, t=50, b=40))
    lab = QC["best"].astype(str)
    rows = [k for k in ids if (lab == k).any()] + ([R.U] if (lab == R.U).any() else [])
    a = pd.to_numeric(QC["best_agreement"], errors="coerce")
    z = [[float(((lab == r) & (a == k)).sum() / max(1, (lab == r).sum())) for k in ks] for r in rows]
    fi = go.Figure(go.Heatmap(z=z, x=ks, y=rows, colorscale="Blues", zmin=0, zmax=1,
                              hovertemplate="%{y}<br>%{x} of " + str(n_dr) + " draws: %{z:.1%} of its cells"
                                            "<extra></extra>", colorbar=dict(title=dict(text="share"))))
    fi.update_layout(height=max(280, 24 * len(rows) + 120), template="plotly_white",
                     title=_layout_title("By final best call: share of its cells at each agreement"),
                     xaxis=dict(title="agreement (draws)", dtick=1), yaxis=dict(autorange="reversed"),
                     margin=dict(l=10, r=10, t=50, b=40))
    return fo, fi


def render_transfer(out):
    man, T = B.load(out)
    meta = man["meta"]
    ids = meta["identities"]
    col = dict(zip(T["colours"].identity.astype(str), T["colours"].colour))
    prov = meta["provenance"]
    pt = prov["transfer"]
    nA, nB = _names(prov)
    DM = T.get("depth_match")
    figs = {"depth": _transfer_depth_fig(T["depth"], DM, (nA, nB))}
    UA = T["umap_A"]
    from scarlet.html_report import certainty as C_
    specs = [_label_spec("umapA", f"A ({nA})", UA, [("labels of A", UA["label"])], col,
                         spaces=_space_names(meta["umap"]["A"]))]
    n_dr = meta.get("n_draws") or pt.get("draws") or 0
    k_need = -(-4 * int(n_dr) // 5) if n_dr else 0
    cert_html, umap_html = "", ""
    if "query_umap" in T:   # the switchable query UMAP, boundary maps, per-identity table, per-cell data
        QU = T["query_umap"]
        outs = []
        for name, pre, e in (("final best", "", "best"), ("final unique", "", "unique"),
                             ("first pass best", "first_", "best"), ("first pass unique", "first_", "unique")):
            k = pd.to_numeric(QU[f"{pre}{e}_agreement"], errors="coerce")
            hover = [] if n_dr <= 1 else [("agreement", [f"{int(v)} of {n_dr}" if np.isfinite(v) else "" for v in k])]
            amb = C_.intensity_ambiguous(QU[f"{pre}best"], QU[f"{pre}unique"])
            outs.append(dict(name=name, label=QU[f"{pre}{e}"],
                             shade={"agreement": C_.intensity_from_need(k, n_dr, k_need), "ambiguity": amb},
                             reason=QU[f"{pre}reason_{e}"], runner_up=QU[f"{pre}runner_up"], hover=hover))
        anc_flag = QU["anchor"].astype(str).str.lower().eq("true").to_numpy()
        specs.append(C_.umap_spec("qumap", f"B ({nB})", QU, outs, col, anchor=anc_flag,
                                  views=["label"] + (["agreement"] if n_dr > 1 else []) +
                                        ["ambiguity", "anchors", "unlabelled"],
                                  spaces=_space_names(meta["umap"]["B"])))
        umap_html = (f"<div class='grid2'><div>{C_.umap_controls(specs[0], highlight=True)}</div>"
                     f"<div>{C_.umap_controls(specs[1], highlight=True)}</div></div>"
                     f"<p class='note'>Each UMAP is computed on its dataset's own cells (display only). B's views: "
                     f"<i>label</i>; <i>agreement</i> (the label's colour by how many of the {n_dr} draws agree: faintest "
                     f"at {k_need} of {n_dr}, the least a final call needs, full at {n_dr} of {n_dr}); <i>ambiguity</i> "
                     f"(the boundary cells, a best call and no unique call, in full colour, the others faint; seen on "
                     f"the best calls); <i>anchors</i> (stable anchors outlined); <i>unlabelled</i> "
                     "(labelled cells faint; unlabelled cells grey, one marker per reason: circle background, diamond "
                     "ambiguity, x undetermined across draws, triangle flip test). Hover: label, agreement, "
                     "runner-up, reason, stable anchor. Highlight: one identity in full colour, the others faint. "
                     f"A: {_space_note(meta['umap']['A'])} B: {_space_note(meta['umap']['B'])}</p>")
        BR = T.get("boundary", pd.DataFrame(columns=["stage", "best", "runner_up", "cells"]))
        nb = {st: int(BR[BR.stage.astype(str) == st].cells.sum()) if len(BR) else 0 for st in ("first pass", "refined")}
        for st, div in (("first pass", "bnd_first"), ("refined", "bnd_refined")):
            f = C_.boundary_fig(BR, st, f"Boundary cells, {st}: best call x runner-up ({nb[st]} cells)", ids)
            if f is not None:
                figs[div] = f
        PI = T["per_identity"]
        cert_html = f"""
<h2>3b. Certainty, boundary cells and refusals</h2>
<p class="note">A boundary cell has a best call but no unique call: SCARLET could not separate its best identity from
the runner-up (the identity ranked next across the draws). Counts per best identity x runner-up; first pass
{nb['first pass']} cells, refined {nb['refined']}.</p>
<div class="grid3">{''.join(f'<div id="{d}"></div>' if d in figs else
                           f"<p class='note'>{_esc(t)}: no boundary cell</p>"
                           for d, t in (("bnd_first", "first pass"), ("bnd_refined", "refined")))}</div>
<p><b>Per identity</b> (final labels): cells called best and unique; stable anchors; refined or not, and why; of its
best calls, % at full agreement ({n_dr} of {n_dr} draws) and % boundary; among the cells it ranks first
(scarlet_top1), % unlabelled by reason, best and unique.</p>{table_html(PI, fmt=_id_fmt(col))}"""
    else:   # an older bundle without the query table
        UB = T["umap_B"]
        specs.append(_label_spec("qumap", f"B ({nB})", UB, [("final best", UB["final_best"]),
                                                            ("final unique", UB["final_unique"]),
                                                            ("stable anchors", UB["anchor"])], col,
                                 spaces=_space_names(meta["umap"]["B"])))
        umap_html = (f"<div class='grid2'><div>{C_.umap_controls(specs[0], highlight=True)}</div>"
                     f"<div>{C_.umap_controls(specs[1], highlight=True)}</div></div>")
    AN = T["analytics"]
    figs["analytics"] = _calls_fig(AN, ids, col, f"Final calls per identity in B ({nB})")
    mk, un = T["markers"], T["universe"]
    sc = T["marker_A_vs_B"]
    anc = T["anchors"]
    ref_models = sorted(set(mk.model) - {x for x in mk.model if str(x).startswith("refine_")})
    q_models = sorted(x for x in set(mk.model) if str(x).startswith("refine_"))
    drill = {}
    uA = un[un.model == ref_models[0]].drop(columns=["model"]).set_index("gene")
    uA.index = uA.index.astype(str)
    uB, fA, fB = None, None, None
    INC = T.get("in_class")   # each model's in-class table on the drill-down's genes (models since 0.6.3)

    def inc_of(m_):
        if INC is None or not (INC.model == m_).any():
            return None
        f = INC[INC.model == m_].drop(columns=["model"]).set_index("gene")
        f.index = f.index.astype(str)
        return f
    if q_models:
        uB = un[un.model == q_models[0]].drop(columns=["model"]).set_index("gene")
        uB.index = uB.index.astype(str)
        fA, fB = inc_of(ref_models[0]), inc_of(q_models[0])
    for k in ids:
        ra_l = mk[(mk.model == ref_models[0]) & (mk.cell_type.astype(str) == k)].sort_values("idx").gene.astype(str)
        rq_l = mk[(mk.model.isin(q_models)) & (mk.cell_type.astype(str) == k)].sort_values("idx").gene.astype(str)
        ra, rq = set(ra_l), set(rq_l)
        heats, ranks = [], []
        if rq:   # three parts, one gene order in both maps: common, A only (dropped), B only (added)
            common = _by_detection(uA, [g for g in ra_l if g in rq], k)
            a_only = _by_detection(uA, [g for g in ra_l if g not in rq], k)
            b_only = _by_detection(uB, [g for g in rq_l if g not in ra], k)
            genes = common + a_only + b_only
            groups, i = [], 0
            for nm, gs in (("common to both models", common), ("A only (dropped by refinement)", a_only),
                           ("B only (added by refinement)", b_only)):
                if gs:
                    groups.append([nm, i, i + len(gs)])
                i += len(gs)
            heats.append(_drill_heat(None, uA, ids, k, "heatA", f"{k}: detection in A (A's model, "
                                                               f"{ref_models[0]})", genes=genes, groups=groups,
                                     inc=fA, split=True))
            heats.append(_drill_heat(None, uB, ids, k, "heatB", f"{k}: detection in B (B's refined model, "
                                                               f"{q_models[0]})", genes=genes, groups=groups,
                                     inc=fB, split=True))
            ranks.append(_rank_plot("rankplot", k, mk[mk.model == ref_models[0]], mk[mk.model.isin(q_models)]))
        else:
            heats.append(_drill_heat(mk[mk.model == ref_models[0]], uA, ids, k, "heatA",
                                     f"{k}: A's model, {ref_models[0]} (highest detection in {k} first)"))
            heats.append(_heat("heatB", f"{k}: not refined (not anchored stably)", [], [], []))
            if q_models:
                ranks.append(dict(div="rankplot", title=f"{k}: not refined, no ranks to compare", common=[],
                                  a_only=[], b_only=[]))
        s = sc[sc.identity.astype(str) == k]
        scat = [_scatter("scat", f"{k}: each reference marker's detection in A vs in the B cells called {k} "
                                 f"({int(s.B_cells_called.max()) if len(s) else 0} cells)",
                         s.detection_A.round(4).tolist(), s.detection_B_called.round(4).fillna(-0.05).tolist(),
                         s.gene.astype(str).tolist(), "detection in A (reference model)",
                         f"detection in B cells called {k}")]
        a_ = anc[anc.identity.astype(str) == k]
        info = (f"<b>Stable anchors</b>: {', '.join(f'{r.origin}: {r.anchors}' for r in a_.itertuples()) or 'none'}; "
                f"refined: {'yes' if a_.refined.astype(str).eq('True').any() else 'no (not anchored stably)'}"
                f"<br><b>Markers before vs after refinement</b>: common {len(ra & rq)}, A only (dropped) "
                f"{len(ra - rq)}, B only (added) {len(rq - ra)}"
                if rq else f"<b>Stable anchors</b>: {', '.join(f'{r.origin}: {r.anchors}' for r in a_.itertuples()) or 'none'}; "
                           "not refined")
        drill[k] = dict(heat=heats, scatter=scat, rank=ranks, html={"drillhtml": info})
    head = [("transfer folder", pt.get("folder")), ("SCARLET", pt.get("scarlet")), ("git commit", pt.get("git_commit")),
            ("refine", pt.get("refine")),
            ("depth matching", pt.get("depth_match")),
            ("draws", _draws_text(pt)),
            ("technology", pt.get("technology")), ("flip", pt.get("flip")), ("alpha", pt.get("alpha")),
            ("blocklist", _bl_text(pt["blocklist"])),
            ("cell-cycle genes", pt.get("cell_cycle")),
            ("transfer at package defaults", _rs.options_report_text(pt.get("non_default_parameters"))),
            (f"A ({nA}, reference)", f"{prov['reference']['file'].get('path')}; {prov['reference']['cells']} x "
                                     f"{prov['reference']['genes']}; labels {prov['reference']['labels']} "
                                     f"({prov['reference']['identities']} identities)"),
            ("A sha256", _hash_text(prov["reference"]["file"])),
            (f"B ({nB}, query)", f"{prov['query']['file'].get('path')}; {prov['query']['cells']} x "
                                 f"{prov['query']['genes']}"),
            ("B sha256", _hash_text(prov["query"]["file"])),
            ("origins", ", ".join(map(str, pt.get("origins") or []))),
            ("UMAPs", f"A {meta['umap']['A']['cells_shown']} of {meta['umap']['A']['cells_total']} cells, B "
                      f"{meta['umap']['B']['cells_shown']} of {meta['umap']['B']['cells_total']}; "
                      f"{meta['umap']['A']['method']}"),
            ("report", f"SCARLET {prov['report']['scarlet']}")]
    dmcols = [c for c in ["draw", "origin", "cells", "direction", "thinned", "method", "depth_ref", "depth_query",
                          "depth_after", "shared_genes", "seed"] if DM is not None and c in DM]
    refused = AN[AN.identity.astype(str).str.startswith("refused")]
    nc = T["no_calls"]
    six = meta["refinement_summary"]
    qc = T.get("qc")
    ch = T["changed_per_identity"]
    sup = T.get("weak_classes")
    weak_t = ch.merge(anc.groupby("identity", as_index=False).agg(anchors=("anchors", "sum"),
                                                                   refined=("refined", lambda v: any(map(
                                                                       lambda x: str(x) == "True", v)))),
                      on="identity", how="left")
    if qc is not None:
        weak_t = weak_t.merge(qc.groupby("identity", as_index=False).agg(
            query_markers=("query_markers", "sum"), conserved=("conserved", "sum"),
            share_conserved=("share_conserved", "mean"),
            doubtful_transfer=("doubtful_transfer", lambda v: any(map(lambda x: str(x) == "True", v)))),
            on="identity", how="left")
        weak_t = weak_t.sort_values(["share_conserved", "anchors"], ascending=[True, True], na_position="first")
        rank_note = ("ranked by SCARLET's transfer QC share of conserved markers, lowest first (ties: fewer anchors); "
                     "doubtful transfer is SCARLET's own QC flag")
    else:
        weak_t = weak_t.sort_values("anchors")
        rank_note = "ranked by anchors, fewest first"
    sup_html = ""
    if sup is not None:
        sup_html = table_html(sup, [c for c in ["rank", "identity", "cv_recall_best", "most_confused_with", "cells",
                                               "under_50_cells", "median_genes_detected"] if c in sup],
                              fmt=_id_fmt(col))
    from scarlet.html_report import bundle as _b
    sup_dir = os.path.join(out, "supplement_A")
    if os.path.exists(os.path.join(sup_dir, "bundle", "manifest.json")):
        _, ST = _b.load(sup_dir)
        if "weak_classes" in ST:
            s2 = ST["weak_classes"]
            sup_html = table_html(s2, [c for c in ["rank", "identity", "cv_recall_best", "most_confused_with",
                                                   "most_confused_share", "cells", "under_50_cells",
                                                   "median_genes_detected"] if c in s2], fmt=_id_fmt(col))
    tables, cells_html, agree_figs = {}, "", ""
    if n_dr > 1:
        if "query_cells" in T:
            fo, fi = _agreement_figs(T["query_cells"], ids, col, n_dr, k_need)
            figs["agree_all"], figs["agree_id"] = fo, fi
            agree_figs = '<div id="agree_all"></div><div id="agree_id"></div>'
        agree_html = (f"""<p class="note">Each cell's final call is the outcome most draws give, when at least 80 % of
the draws give it (here {k_need} of {n_dr}); otherwise the cell is unassigned (undetermined across draws).
Agreement: the number of draws giving the cell's most frequent outcome ('unassigned' counts as an outcome). It
measures how stable a call is under the random depth thinning and seeds of the draws, not whether it is correct.</p>"""
                      + agree_figs
                      + (table_html(T["agreement"]) if "agreement" in T else "")
                      + '<p class="note">Each draw\'s own result (cells called): the spread across the draws.</p>'
                      + (table_html(T["per_draw"]) if "per_draw" in T else ""))
    else:
        agree_html = ('<p class="note">Quick mode: one draw, so no consensus and no per-cell agreement; another seed '
                      'can change some calls. The draw\'s result:</p>'
                      + (table_html(T["per_draw"]) if "per_draw" in T else ""))
    from scarlet.html_report import known as K
    known_html, kfig = K.section(T, meta, col, ids, "3c", table_html, _esc, _id_fmt(col))
    if kfig is not None:
        figs["known_conf"] = kfig
    if "query_cells" in T and "query_umap" in T:
        QC = T["query_cells"]
        u, b = QC.unique.astype(str).to_numpy(), QC.best.astype(str).to_numpy()
        ag = pd.to_numeric(QC.unique_agreement, errors="coerce").fillna(0).to_numpy()
        grp = np.where((b == R.U), 1, np.where(u == R.U, 2, np.where((n_dr > 1) & (ag < n_dr), 3, 4)))
        order = np.argsort(grp, kind="stable")
        QCs = QC.iloc[order].reset_index(drop=True)
        gnames = ["worth a look: refused, boundary or below full agreement", "refused (no best call)",
                  "boundary (best call, no unique call)", f"called, below full agreement ({n_dr} of {n_dr})",
                  "called at full agreement", "every cell"]
        tables["cells"], tab = _cell_table(
            QCs, "cells", "scarlet_transfer_cells.tsv", grp[order], gnames, 0,
            notes={"runner_up": "the identity ranked after the top identity (in the most draws); for a refused cell "
                                "too, unlike the learn report's runner-up"}, page=25)
        what = ((f"agreement (draws of {n_dr} giving the consensus outcome), reasons, top-ranked identity, runner-up "
                 "(and in how many draws)") if n_dr > 1 else "reasons, top-ranked identity, runner-up (quick mode: one "
                "draw, no agreement)")
        ng = [int((grp == g).sum()) for g in (1, 2, 3, 4)]
        tiles = "".join(f'<div class="tile"><b>{c:,}</b>{_esc(nm)} ({_pct(c, len(QC))})</div>'
                        for c, nm in zip(ng, gnames[1:5]) if not (nm.startswith("called, below") and n_dr <= 1))
        cells_html = (f"<h2>8. Per-cell data</h2><p class='note'>B's {len(QC):,} cells by how their final call came "
                      "out:</p><div class='tiles'>" + tiles + "</div>" +
                      C_.download_button("cells", len(QC), "cells") +
                      f"<details><summary>Browse the cells</summary><p class='note'>Per cell: final and first-pass "
                      f"calls, {what}, stable anchor, flags. <i>Show</i> picks a group (the cells worth a look "
                      "first); 25 rows at a time.</p>" + tab + "</details>")
    LB_cells = int(prov["query"]["cells"])
    best_called = AN[(AN.entry == "best") & ~AN.identity.astype(str).str.startswith("refused")].called.sum()
    uniq_called = AN[(AN.entry == "unique") & ~AN.identity.astype(str).str.startswith("refused")].called.sum()
    ur = refused[refused.entry == "unique"]
    summary = _ul([
        f"{_esc(nA)} (A, {prov['reference']['cells']:,} cells, {prov['reference']['identities']} identities) &rarr; "
        f"{_esc(nB)} (B, {LB_cells:,} cells).",
        f"Final calls: best {_pct(int(best_called), LB_cells)} of B's cells, unique {_pct(int(uniq_called), LB_cells)}; "
        "unique calls refused: " + (", ".join(f"{_esc(str(r.identity).replace('refused: ', ''))} "
                                             f"{_pct(int(r.called), LB_cells)}" for r in ur.itertuples()) or "none")
        + ".",
        f"Depth matching: {_depth_line(DM, (nA, nB))}." if DM is not None and len(DM) else "No depth matching.",
        (f"Draws: {n_dr}; a final call needs {k_need} of {n_dr}." if n_dr > 1 else "Quick mode: one draw."),
        (f"Transfer check: {int(qc.doubtful_transfer.astype(str).eq('True').sum())} of {len(qc)} identities flagged "
         "doubtful transfer." if qc is not None and len(qc) else ""),
        (f"Identities of A with no call in B: {len(nc)} (section 3)." if len(nc) else "")])
    body = f"""
<h1>SCARLET transfer report: {_esc(nA)} &rarr; {_esc(nB)}</h1>
<p class="note">A = {_esc(nA)} (the reference), B = {_esc(nB)} (the query). Only values SCARLET computed are shown; the
report reads the transfer's outputs and changes no call. UMAPs are display only. Supplement:
<a href="{_esc(meta['supplement'])}">the learn report on A</a> (its self-assessment).</p>
<h2>Summary</h2>{summary}
<h2>1. Depth of A ({_esc(nA)}) and B ({_esc(nB)}), before and after depth matching</h2><div id="depth"></div>
<p class="note">Boxes: genes detected per cell before matching (every gene). Markers: the record's medians on the
genes both share, before matching (open circle) and after (diamond: median over the draws, bar: their range). The
per-draw table is in the provenance (section 9).</p>
<h2>2. UMAPs: A ({_esc(nA)}) and B ({_esc(nB)})</h2>{umap_html}
<h2>3. Analytics on the final calls</h2><div id="analytics"></div>
<p class="note">Refusals by reason (as transfer recorded them per cell; 'undetermined across draws': no outcome
won clearly across the draws):</p>{table_html(refused, ["entry", "identity", "called", "share"])}
<p class="note">Identities of A with no call in B:</p>{table_html(nc) if len(nc) else "<p class='note'>none</p>"}
{cert_html}
{known_html}
<h2>4. Identity drill-down</h2>{_picker(ids)}
<p class="note">{_esc(BG_NOTE)}. With refinement, both maps list the same genes in three parts: the markers common
to A's model and B's refined model; those only in A's model (dropped by refinement); those only in B's refined model
(added). Left, their detection in A's model (A's cells per label); right, in B's refined model (B's stable anchors per
cell type, the cells it learned from). Every value is the model's own estimate, read from its in-class table: the
share of the cell type's cells detecting the gene, averaged over the bootstrap of learn, the same estimate the model
uses for its markers. A model written before SCARLET 0.6.3 has no such table: the other model's markers are then blank
(the hover says 'not a marker of this model'). Below, each marker's rank in both models.</p>
<div class="grid2"><div class="scrollbox"><div id="heatA"></div></div><div class="scrollbox"><div id="heatB"></div></div></div>
<div class="grid2"><div id="rankplot"></div><div id="scat"></div></div><div id="drillhtml" class="card"></div>
<h2>5. Weak-class flags</h2>
<p><b>From the reference</b> (the supplement's self-assessment; ranked by cross-validated best recall, weakest
first; under 50 cells is a reporting convention, not a SCARLET threshold):</p>{sup_html}
<p><b>From the transfer</b> ({rank_note}; evidence only, no threshold):</p>
{table_html(weak_t, fmt=_id_fmt(col))}
<h2>6. What refinement changed</h2>
{_six_text(six)}
<p>The details are on a page of their own: <a href="{REFINEMENT_PAGE}">what refinement changed</a> (first pass and
refined calls side by side, a map of every change of call, and the list).</p>
<h2>7. Agreement across the draws</h2>
{agree_html}
{cells_html}
<h2>9. Provenance (both sides)</h2>{_summary_block(prov)}{_prov_table(head)}
<p class="note">Depth matching, draw by draw (the record's medians on the shared genes):</p>
{table_html(DM, dmcols) if DM is not None and len(DM) else "<p class='note'>no depth matching</p>"}
{_footnote(len(ids))}"""
    data = _data(figs, drill)
    data["qumaps"] = specs
    data["tables"] = tables
    return _page(f"SCARLET transfer report: {nA} to {nB}", body, data)


REFINEMENT_PAGE = "refinement.html"
SB_COLS = ["identity", "first_pass_best", "refined_best", "first_pass_unique", "refined_unique", "first_pass_best_pct",
           "refined_best_pct", "first_pass_unique_pct", "refined_unique_pct"]


def _six_text(six):
    def n(v):
        return f"{v:,} cell" + ("" if v == 1 else "s")
    return "<ul>" + "".join(
        f"<li><b>{e}</b>: {n(six[e]['refused_then_called'])} refused in the first pass and called after "
        f"refinement; {n(six[e]['changed_identity'])} changed cell type; {n(six[e]['called_then_refused'])} called "
        f"in the first pass and refused after.</li>" for e in ("best", "unique")) + "</ul>"


def _changes_fig(X6, e, ids, title):
    """first-pass call (rows) x final call (columns): the cells whose call changed; the diagonal is empty (no
    change) and so are the pairs with no cell"""
    import plotly.graph_objects as go
    X = X6[X6.entry.astype(str) == e]
    if not len(X):
        return None
    order = [k for k in list(ids) + [R.U] if k in set(X.first_pass.astype(str)) | set(X.final.astype(str))]
    M = X.pivot_table(index="first_pass", columns="final", values="cells", aggfunc="sum").reindex(index=order,
                                                                                                 columns=order)
    z = M.to_numpy(float)
    fig = go.Figure(go.Heatmap(z=np.where(np.isfinite(z), z, None).tolist(), x=order, y=order, colorscale="Blues",
                               colorbar=dict(title=dict(text="cells")), hoverongaps=False,
                               hovertemplate="first pass %{y}<br>final %{x}<br>%{z} cells<extra></extra>"))
    fig.update_layout(title=_layout_title(title), height=max(360, 18 * len(order) + 200),
                      width=max(480, 18 * len(order) + 330), template="plotly_white",
                      margin=dict(l=10, r=10, t=50, b=10), xaxis=dict(title="final call", tickangle=-45, side="top"),
                      yaxis=dict(title="first-pass call", autorange="reversed", automargin=True))
    return fig


def render_refinement(out):
    """the transfer report's page on refinement: first pass and refined side by side, a map of every change of call
    (first pass x final, best and unique), and the list of changes, largest first"""
    man, T = B.load(out)
    meta = man["meta"]
    ids = meta["identities"]
    col = dict(zip(T["colours"].identity.astype(str), T["colours"].colour))
    nA, nB = _names(meta["provenance"])
    X6, SB, six = T["refinement_changes"], T["first_pass_vs_refined"], meta["refinement_summary"]
    figs, maps = {}, []
    for e in ("best", "unique"):
        f = _changes_fig(X6, e, ids, f"{e}: cells whose call changed, first pass (rows) to final (columns)")
        if f is not None:
            figs[f"chg_{e}"] = f
            maps.append(f'<div class="scrollbox"><div id="chg_{e}"></div></div>')
        else:
            maps.append(f"<p class='note'>{e}: no call changed</p>")
    X = X6.assign(_o=X6.entry.astype(str).map({"best": 0, "unique": 1})).sort_values(
        ["_o", "cells", "first_pass", "final"], ascending=[True, False, True, True]).drop(columns="_o")
    body = f"""
<h1>What refinement changed: {_esc(nA)} &rarr; {_esc(nB)}</h1>
<p class="note">A page of the <a href="report.html">transfer report</a> (section 6). The first pass calls B with A's
model; refinement learns again from B's stable anchors and gives the final calls. Both are SCARLET outputs
(first_pass.tsv and labels.tsv); refused = 'unassigned'.</p>
<h2>Summary</h2>{_six_text(six)}
<h2>First pass and refined, side by side</h2>
<p class="note">Cells per call in B, and their share of B's cells.</p>{table_html(SB, SB_COLS, fmt=_id_fmt(col))}
<h2>Every change of call</h2>
<p class="note">Rows: the first-pass call; columns: the final call; colour and hover: the cells that went from one to
the other. 'unassigned' as a row: cells refused in the first pass and called after refinement; as a column: called,
then refused.</p>
<div class="grid2">{''.join(maps)}</div>
<details><summary>The list of changes ({len(X6):,} rows, largest first)</summary>{table_html(X)}</details>
{_footnote(len(ids))}"""
    return _page(f"What refinement changed: {nA} to {nB}", body, _data(figs, {}))


def _draws_text(pt):
    n = pt.get("draws")
    if not n:
        return None
    seeds = pt.get("seeds") or []
    s = ", ".join(map(str, seeds[:5])) + (f", ... {seeds[-1]}" if len(seeds) > 5 else "")
    if n == 1:
        return f"1, quick mode (seed {s}; no consensus, no per-cell agreement)"
    return f"{n} (seeds {s}); stable anchors over all {n} draws"


# --------------------------------------------------------------------------------------------------------------
# assign report
# --------------------------------------------------------------------------------------------------------------
def render_assign(out):
    man, T = B.load(out)
    meta = man["meta"]
    ids = meta["identities"]
    col = dict(zip(T["colours"].identity.astype(str), T["colours"].colour))
    prov = meta["provenance"]
    C = T["cells"]
    import plotly.graph_objects as go
    fd = go.Figure()
    for k in ids:
        g = C.genes_detected[C.best.astype(str) == k]
        if len(g):
            fd.add_trace(go.Box(x=g, name=k, marker_color=col[k], boxpoints=False, orientation="h"))
    fd.update_layout(height=max(300, 22 * len(ids) + 80), showlegend=False, xaxis_title="genes detected per cell",
                     margin=dict(l=10, r=10, t=10, b=40), template="plotly_white", yaxis=dict(autorange="reversed"))
    U = T["umap"]
    from scarlet.html_report import certainty as C_
    specs = [_label_spec("umap2", "Query", U, [("best", U["best"]), ("unique", U["unique"])], col,
                         spaces=_space_names(meta["umap"]))]
    figs = {"depth": fd}
    AN = T["analytics"]
    figs["analytics"] = _calls_fig(AN, ids, col, "Calls per identity")
    mk = T["markers"]
    uni = T["universe"].set_index("gene")
    uni.index = uni.index.astype(str)
    sc = T["marker_detection"]
    drill = {}
    for k in ids:
        s = sc[sc.identity.astype(str) == k]
        drill[k] = dict(heat=[_drill_heat(mk, uni, ids, k, "heat", f"{k}: its markers (highest detection in {k} "
                                                                    "first) x every identity (model)")],
                        scatter=[_scatter("scat", f"{k}: each marker's detection in the model vs in cells called {k}",
                                          s.detection_model.round(4).tolist(),
                                          s.detection_called.round(4).fillna(-0.05).tolist(), s.gene.astype(str).tolist(),
                                          "detection rate in the model", f"detection in cells called {k}")],
                        html={})
    pm, pa = prov["model"], prov["assign"]
    head = [("model", pm.get("path")), ("SCARLET (model)", pm.get("scarlet")), ("git commit", pm.get("git_commit")),
            ("technology", pm.get("technology")),
            ("blocklist", _bl_text(pm["blocklist"])), ("annotated data", pa["annotated"].get("path")),
            ("sha256", _hash_text(pa["annotated"]) + f"  ({prov['hash_note']})"),
            ("assign steps", json.dumps(pa.get("steps"))), ("alpha", f"{pa['alpha']} ({pa['alpha_note']})"),
            ("cell-cycle genes (model)", pa.get("cell_cycle")),
            ("assign at package defaults", _rs.options_report_text(pa.get("non_default_parameters"))),
            ("cells x genes", f"{prov['data']['cells']} x {prov['data']['genes']}"),
            ("UMAP", f"{meta['umap']['cells_shown']} of {meta['umap']['cells_total']} cells; {meta['umap']['method']}")]
    refused = AN[AN.identity.astype(str).str.startswith("refused")]
    called = AN[~AN.identity.astype(str).str.startswith("refused")].sort_values(["entry", "called"])
    nc = T["no_calls"]
    from scarlet.html_report import known as K
    known_html, kfig = K.section(T, meta, col, ids, "3b", table_html, _esc, _id_fmt(col))
    if kfig is not None:
        figs["known_conf"] = kfig
    n = len(C)
    items = [f"{n:,} cells annotated with a model of {len(ids)} identities."]
    for e in ("best", "unique"):
        ref = C[e].astype(str) == R.U
        items.append(f"{e.capitalize()} call: {_pct(int((~ref).sum()), n)} called; refused "
                     f"{_pct(int(ref.sum()), n)} ({_reasons_text(C[f'reason_{e}'][ref], n) or 'none'}).")
    if len(nc):
        items.append(f"Identities of the model with no call: {len(nc)} (section 5).")
    body = f"""
<h1>SCARLET assign report (query without labels)</h1>
<p class="note">Calls and refusals only; only values SCARLET computed are shown. UMAPs are display only.</p>
<h2>Summary</h2>{_ul(items)}
<h2>1. Depth: genes per cell per called identity</h2><div id="depth"></div>
<h2>2. UMAP</h2><p class="note">Choose the output; highlight one identity if you wish (all in full colour by
default). {_space_note(meta["umap"])}</p>{C_.umap_controls(specs[0], highlight=True)}
<h2>3. Analytics</h2><div id="analytics"></div>
<p class="note">Refusals by reason ({R.LABEL}):</p>{table_html(refused, ["entry", "identity", "called", "share"])}
{known_html}
<h2>4. Identity drill-down</h2>{_picker(ids)}<p class="note">{_esc(BG_NOTE)}.</p>
<div class="scrollbox"><div id="heat"></div></div><div id="scat"></div>
<h2>5. Flags</h2><p class="note">Identities with fewest calls first; flip and low positive evidence are SCARLET's
per-cell flags. Identities of the model with no call:</p>{table_html(nc) if len(nc) else "<p class='note'>none</p>"}
{table_html(called, fmt=_id_fmt(col))}
<h2>6. Provenance</h2>{_summary_block(prov)}{_prov_table(head)}
{_footnote(len(ids))}
"""
    data = _data(figs, drill)
    data["qumaps"] = specs
    return _page("SCARLET assign report", body, data)


JS = r"""
const D = JSON.parse(document.getElementById('report-data').textContent);
const cfg = {responsive: true, displaylogo: false};
for (const [k, f] of Object.entries(D.figs)) {
  const el = document.getElementById(k); if (el) Plotly.newPlot(el, f.data, f.layout, cfg);
}
// the identity drill-down: heatmaps in scrolling boxes (sized by genes and identities), scatter, notes
function pick(k) {
  const d = D.drill[k]; if (!d) return;
  for (const h of (d.heat || [])) {
    if (!document.getElementById(h.div)) continue;
    if (h.split && (h.groups || []).length && h.genes.length) { splitHeat(h); continue; }
    const shapes = [], notes = [];
    for (const [nm, a, b] of (h.groups || [])) {
      if (a > 0) shapes.push({type: 'line', xref: 'paper', x0: 0, x1: 1, yref: 'y', y0: a - 0.5, y1: a - 0.5,
                              line: {color: '#000000', width: 1.5}});
      notes.push({xref: 'paper', x: 1.0, xanchor: 'left', yref: 'y', y: (a + b - 1) / 2, text: nm + ' (' + (b - a) + ')',
                  showarrow: false, font: {size: 12}});
    }
    Plotly.react(h.div, h.genes.length ? [{type: 'heatmap', z: h.z, x: h.cols, y: h.genes, zmin: 0, zmax: 1,
      colorscale: 'Viridis', hovertemplate: 'gene %{y}<br>identity %{x}<br>detection %{z:.3f}<extra></extra>',
      colorbar: {title: {text: 'detection'}, x: 1.12}}] : [],
      {title: {text: h.title, x: 0.01, xanchor: 'left', font: {size: 13}},
       height: Math.max(320, 17 * h.genes.length + 190), width: Math.max(480, 36 * h.cols.length + 330),
       margin: {l: 120, r: 150, t: 50, b: 140}, shapes: shapes, annotations: notes,
       yaxis: {autorange: 'reversed', tickfont: {size: 12}, automargin: true},
       xaxis: {tickangle: -45, tickfont: {size: 12}, automargin: true}, template: 'plotly_white'},
      {displaylogo: false, responsive: false});
  }
  for (const r of (d.rank || [])) { if (document.getElementById(r.div)) rankPlot(r); }
  for (const s of (d.scatter || [])) {
    if (!document.getElementById(s.div)) continue;
    Plotly.react(s.div, [{type: 'scatter', mode: 'markers', x: s.x, y: s.y, text: s.text,
      hovertemplate: '%{text}<br>x %{x:.3f}<br>y %{y:.3f}<extra></extra>', marker: {size: 6}},
      {type: 'scatter', mode: 'lines', x: [0, 1], y: [0, 1], line: {dash: 'dot', color: '#999'}, hoverinfo: 'skip',
       showlegend: false}],
      {title: {text: s.title, x: 0.01, xanchor: 'left', font: {size: 13}}, height: 440,
       xaxis: {title: {text: s.xt}, range: [-0.07, 1.02]},
       yaxis: {title: {text: s.yt}, range: [-0.07, 1.02], automargin: true}, showlegend: false,
       template: 'plotly_white', margin: {l: 80, r: 20, t: 50, b: 50}}, cfg);
  }
  for (const [div, h] of Object.entries(d.html || {})) {
    const el = document.getElementById(div); if (el) el.innerHTML = h;
  }
}
// a heatmap in separate blocks (its groups: common, A only, B only), one shared colour scale; values the model
// has none for, computed by the report from the data, say so in the hover
const GAPS = {1: 'not a marker of this model', 2: 'not among the genes this model\'s marker search scored'};
function splitHeat(h) {
  const G = h.groups, n = h.genes.length, gap = 40, rowPx = 17, top = 50, bottom = 140;
  const plotPx = rowPx * n + gap * (G.length - 1), H = top + bottom + plotPx;
  const tr = [], lay = {title: {text: h.title, x: 0.01, xanchor: 'left', font: {size: 13}}, height: H,
    width: Math.max(480, 36 * h.cols.length + 330), margin: {l: 120, r: 150, t: top, b: bottom},
    template: 'plotly_white', annotations: []};
  let y0 = plotPx;
  G.forEach(([nm, a, b], gi) => {
    const ax = gi === 0 ? '' : String(gi + 1), hi = y0 / plotPx, lo = (y0 - rowPx * (b - a)) / plotPx;
    y0 -= rowPx * (b - a) + gap;
    const z = h.z.slice(a, b), gp = (h.gaps || []).slice(a, b), gs = h.genes.slice(a, b);
    tr.push({type: 'heatmap', z: z, x: h.cols, y: gs, yaxis: 'y' + ax, xaxis: 'x', zmin: 0, zmax: 1,
      colorscale: 'Viridis', showscale: gi === 0, colorbar: {title: {text: 'detection'}, x: 1.12, len: 0.5, y: 1,
      yanchor: 'top'}, hoverongaps: true,
      customdata: z.map((row, i) => row.map((v, j) => 'gene ' + gs[i] + '<br>identity ' + h.cols[j] + '<br>' +
        (v === null ? (GAPS[gp.length ? gp[i][j] : 1] || 'no value') : 'detection ' + v.toFixed(3)))),
      hovertemplate: '%{customdata}<extra></extra>'});
    lay['yaxis' + ax] = {domain: [Math.max(0, lo), Math.min(1, hi)], autorange: 'reversed', tickfont: {size: 12},
                         automargin: true};
    lay.annotations.push({xref: 'paper', x: 1.0, xanchor: 'left', yref: 'paper', y: (lo + hi) / 2, yanchor: 'middle',
                          text: nm.replace(' (', '<br>(') + ' (' + (b - a) + ')', showarrow: false, font: {size: 12}});
  });
  lay.xaxis = {tickangle: -45, tickfont: {size: 12}, automargin: true, anchor: 'y' + (G.length > 1 ? G.length : '')};
  Plotly.react(h.div, tr, lay, {displaylogo: false, responsive: false});
}
// each marker's rank in both models: common markers joined by a line, the others on their own side
function rankPlot(r) {
  const lx = [], ly = [], lt = [];
  for (const [g, a, b] of r.common) { lx.push('A', 'B', null); ly.push(a, b, null); lt.push(g, g, null); }
  const tr = [
    {type: 'scatter', mode: 'lines+markers', name: 'common (' + r.common.length + ')', x: lx, y: ly, text: lt,
     line: {color: '#8c8c8c', width: 1}, marker: {color: '#555555', size: 6},
     hovertemplate: '%{text}<br>rank %{y} in %{x}<extra>common</extra>'},
    {type: 'scatter', mode: 'markers', name: 'A only (' + r.a_only.length + ')', x: r.a_only.map(() => 'A'),
     y: r.a_only.map(v => v[1]), text: r.a_only.map(v => v[0]), marker: {color: '#E69F00', size: 9, symbol: 'triangle-right'},
     hovertemplate: '%{text}<br>rank %{y} in A; not in B\'s model<extra>A only</extra>'},
    {type: 'scatter', mode: 'markers', name: 'B only (' + r.b_only.length + ')', x: r.b_only.map(() => 'B'),
     y: r.b_only.map(v => v[1]), text: r.b_only.map(v => v[0]), marker: {color: '#0072B2', size: 9, symbol: 'triangle-left'},
     hovertemplate: '%{text}<br>rank %{y} in B; not in A\'s model<extra>B only</extra>'}];
  Plotly.react(r.div, tr, {title: {text: r.title, x: 0.01, xanchor: 'left', font: {size: 13}}, height: 440,
    template: 'plotly_white', margin: {l: 60, r: 20, t: 50, b: 50},
    xaxis: {type: 'category', categoryorder: 'array', categoryarray: ['A', 'B'],
            ticktext: ["A's model", "B's refined model"], tickvals: ['A', 'B']},
    yaxis: {title: {text: 'rank'}, autorange: 'reversed', zeroline: false},
    legend: {orientation: 'h', y: -0.12}}, cfg);
}
const sel = document.getElementById('pick');
if (sel) { sel.addEventListener('change', e => pick(e.target.value)); pick(sel.value); }
// the switchable UMAPs (output x view, and a highlight): certainty, anchors, unlabelled reasons
function blend(hex, f) {
  const n = parseInt(hex.slice(1), 16), r = n >> 16, g = (n >> 8) & 255, b = n & 255;
  const m = c => Math.round(255 - (255 - c) * f);
  return 'rgb(' + m(r) + ',' + m(g) + ',' + m(b) + ')';
}
function qHover(S, O, i) {
  let h = S.cell[i] + (S.truth ? '<br>reference label: ' + S.cats[S.truth[i]] : '') + '<br>' + O.name + ': ' +
    S.cats[O.lab[i]];
  for (const [nm, arr] of O.hover) h += '<br>' + nm + ': ' + (arr[i] === null ? '' : arr[i]);
  if (O.ru) h += '<br>runner-up: ' + (O.ru[i] < 0 ? '' : S.cats[O.ru[i]]);
  if (O.rs) h += '<br>reason: ' + (O.rs[i] < 0 ? '' : S.reasons[O.rs[i]]);
  if (S.anchor) h += '<br>stable anchor: ' + (S.anchor[i] ? 'yes' : 'no');
  return h;
}
// without WebGL (switched off, or a browser without it) the UMAPs are drawn as SVG, on a sample taken within each
// identity when there are many cells (S.svg; every cell of the small identities kept)
const WEBGL = (() => { try { const c = document.createElement('canvas');
  return !!(window.WebGLRenderingContext && (c.getContext('webgl') || c.getContext('experimental-webgl'))); }
  catch (e) { return false; } })();
function qTrace(S, O, idx, marker, name) {
  if (!WEBGL && S.svg) { const keep = new Set(S.svg); const j = [];
    idx.forEach((i, p) => { if (keep.has(i)) j.push(p); });
    const pick = a => (Array.isArray(a) ? j.map(p => a[p]) : a);
    marker = Object.assign({}, marker, {color: pick(marker.color), size: pick(marker.size)});
    if (marker.line) marker.line = Object.assign({}, marker.line, {width: pick(marker.line.width)});
    idx = j.map(p => idx[p]); }
  return {type: WEBGL ? 'scattergl' : 'scatter', mode: 'markers', name: name, x: idx.map(i => S.cx[i]),
          y: idx.map(i => S.cy[i]), marker: marker, text: idx.map(i => qHover(S, O, i)),
          hovertemplate: '%{text}<extra></extra>'};
}
function qDraw(S) {
  const oi = +document.getElementById(S.div + '_out').value, view = document.getElementById(S.div + '_view').value;
  const hlEl = document.getElementById(S.div + '_hl'), hl = hlEl ? +hlEl.value : -1;
  const spEl = document.getElementById(S.div + '_space'), si = spEl ? +spEl.value : 0;
  const SP = (S.spaces || [{}])[si] || {};
  S.cx = si > 0 && SP.x ? SP.x : S.x; S.cy = si > 0 && SP.y ? SP.y : S.y;
  const O = S.outputs[oi], U = S.cats.indexOf('unassigned'), n = S.x.length, tr = [];
  const ps = (!WEBGL && S.svg) ? S.svg_size : S.size;
  if (view === 'unlabelled') {
    const lab = [], byR = new Map();
    for (let i = 0; i < n; i++) {
      if (O.lab[i] === U) { const r = O.rs ? O.rs[i] : -1; if (!byR.has(r)) byR.set(r, []); byR.get(r).push(i); }
      else lab.push(i);
    }
    tr.push(qTrace(S, O, lab, {size: ps, opacity: 0.15, color: lab.map(i => S.col[O.lab[i]])}, 'labelled'));
    for (const [r, idx] of [...byR.entries()].sort((a, b) => a[0] - b[0])) {
      const why = r < 0 ? 'reason not recorded' : S.reasons[r];
      tr.push(qTrace(S, O, idx, {size: ps + 3, opacity: 0.9, color: S.unlabelled, symbol: S.symbols[why] || 'circle'},
                     'unlabelled: ' + why));
    }
  } else {
    const by = new Map();
    for (let i = 0; i < n; i++) { const c = O.lab[i]; if (!by.has(c)) by.set(c, []); by.get(c).push(i); }
    const keys = [...by.keys()].sort((a, b) => ((a === U) - (b === U)) || (a - b));
    for (const c of keys) {
      const idx = by.get(c), base = S.col[c], mk = {size: ps, opacity: 0.85, color: base};
      if (O.shade && O.shade[view] && c !== U) mk.color = idx.map(i => blend(base, O.shade[view][i]));
      if (view === 'anchors' && S.anchor) {
        mk.size = idx.map(i => (S.anchor[i] ? ps + 3 : ps));
        mk.line = {color: '#000000', width: idx.map(i => (S.anchor[i] ? 1.2 : 0))};
      }
      if (hl >= 0) { mk.opacity = (c === hl ? 1.0 : 0.07); if (c === hl && !Array.isArray(mk.size)) mk.size = ps + 2; }
      tr.push(qTrace(S, O, idx, mk, S.cats[c]));
    }
  }
  Plotly.react(S.div, tr, {height: 560, margin: {l: 10, r: 10, t: 50, b: 10}, template: 'plotly_white',
    title: {text: S.title + ': ' + O.name + (S.views.length > 1 ? ', ' + view : '') +
           (hl >= 0 ? ' (highlight: ' + S.cats[hl] + ')' : '') +
           (SP.name ? '<br><span style="font-size:11px">space: ' + SP.name + '</span>' : '') +
           ((!WEBGL && S.svg) ? '<br><span style="font-size:11px">WebGL unavailable: drawn as SVG on ' +
             S.svg.length.toLocaleString() + ' of ' + n.toLocaleString() + ' cells, sampled within each identity ' +
             '(identities of up to ' + S.svg_limit.toLocaleString() + ' cells kept whole)</span>' : ''),
           x: 0.01, xanchor: 'left', font: {size: 13}},
    xaxis: {visible: false}, yaxis: {visible: false}, legend: {font: {size: 11}, itemsizing: 'constant'}}, cfg);
}
for (const S of (D.qumaps || [])) {
  for (const w of ['_out', '_view', '_hl', '_space']) {
    const el = document.getElementById(S.div + w); if (el) el.addEventListener('change', () => qDraw(S));
  }
  if (document.getElementById(S.div)) qDraw(S);
}
function tableRows(T) {
  if (T.rows) return T.rows;
  const cols = T.columns.map(c => (c.codes ? c.codes.map(k => (k === null ? null : c.dict[k])) : c.values));
  const rows = [];
  for (let i = 0; i < T.n; i++) rows.push(cols.map(col => col[i]));
  return rows;
}
function tsvDownload(key) {
  const T = D.tables[key]; if (!T) return;
  const lines = [T.cols.join('\t')];
  for (const r of tableRows(T)) lines.push(r.map(v => (v === null ? '' : String(v))).join('\t'));
  const blob = new Blob([lines.join('\n') + '\n'], {type: 'text/tab-separated-values'});
  const a = document.createElement('a'); a.href = URL.createObjectURL(blob); a.download = T.file;
  document.body.appendChild(a); a.click(); a.remove();
}
// the per-cell tables: a 'show' filter (cells worth a look first), headers that sort, the first 200 rows
function esc(s) { return String(s).replace(/[&<>"]/g, c => ({'&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;'})[c]); }
function fmtCell(v, kind) {
  if (v === null || v === undefined || v === '') return '';
  if (kind === 'bool') return (v === true || v === 'True') ? 'yes' : 'no';
  if (typeof v !== 'number') return esc(v);
  if (kind === 'sci') return v === 0 ? '0' : v.toExponential(1);
  if (kind === 'f1') return v.toFixed(1);
  if (kind === 'int') return String(v);
  return Number.isInteger(v) ? String(v) : String(Number(v.toPrecision(3)));
}
function cellTable(key) {
  const T = D.tables[key], showEl = document.getElementById(key + '_show'), box = document.getElementById(key + '_table');
  if (!T || !showEl || !box) return;
  const rows = tableRows(T), state = {col: -1, dir: 1}, step = T.page || 200;
  let shown = step;
  const more = document.getElementById(key + '_more');
  function draw() {
    const g = +showEl.value, last = T.group_names.length - 1;
    let idx = [];
    for (let i = 0; i < rows.length; i++) {
      const gi = T.groups[i];
      if (g === last || (g === 0 ? gi < last - 1 : gi === g)) idx.push(i);
    }
    if (state.col >= 0) {
      const c = state.col;
      idx.sort((a, b) => { const x = rows[a][c], y = rows[b][c];
        if (x === y) return 0; if (x === null) return 1; if (y === null) return -1;
        return (x < y ? -1 : 1) * state.dir; });
    }
    document.getElementById(key + '_count').textContent = idx.length ? ('showing ' + Math.min(shown, idx.length) +
      ' of ' + idx.length.toLocaleString() + ' cells') : 'no cell';
    if (more) more.style.display = idx.length > shown ? '' : 'none';
    let h = '<table><thead><tr>' + T.headers.map((c, j) => '<th title="' + esc(T.tips[j]) + '" data-j="' + j +
      '" style="cursor:pointer">' + esc(c) + (state.col === j ? (state.dir > 0 ? ' &#9650;' : ' &#9660;') : '') +
      '</th>').join('') + '</tr></thead><tbody>';
    for (const i of idx.slice(0, shown)) h += '<tr>' + rows[i].map((v, j) => '<td>' + fmtCell(v, T.kinds[j]) +
      '</td>').join('') + '</tr>';
    box.innerHTML = h + '</tbody></table>';
    for (const th of box.querySelectorAll('th')) th.addEventListener('click', () => {
      const j = +th.dataset.j; state.dir = (state.col === j ? -state.dir : 1); state.col = j; draw(); });
  }
  showEl.value = String(T.default || 0);
  showEl.addEventListener('change', () => { shown = step; draw(); });
  if (more) more.addEventListener('click', () => { shown += step; draw(); });
  draw();
}
for (const key of Object.keys(D.tables || {})) cellTable(key);
"""


def _strict(o):
    """JSON for the browser: NaN / inf become null"""
    if isinstance(o, dict):
        return {k: _strict(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [_strict(v) for v in o]
    if isinstance(o, float) and not np.isfinite(o):
        return None
    return o


def _page(title, body, data):
    payload = json.dumps(_strict(data), default=B._jsonable, allow_nan=False).replace("</", "<\\/")
    return (f"<!doctype html><html lang='en'><head><meta charset='utf-8'><meta name='viewport' "
            f"content='width=device-width, initial-scale=1'><title>{_esc(title)}</title><style>{CSS}</style>"
            f"<script>{_plotlyjs()}</script></head><body>{body}"
            f"<script type='application/json' id='report-data'>{payload}</script><script>{JS}</script>"
            f"</body></html>")
