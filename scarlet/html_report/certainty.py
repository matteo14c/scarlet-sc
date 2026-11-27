"""Certainty, boundary cells, anchors and unlabelled reasons: the views shared by the learn report (from the
cross-validation) and the transfer report (from the draws).

Only values SCARLET wrote are shown: the transfer's per-cell agreement, reasons, runner-up and stable anchors
(first_pass.tsv, labels.tsv); the self-assessment's FDR, margin and runner-up. A boundary cell has a best call but no
unique call: SCARLET could not separate its best identity from the runner-up. The query UMAP is display only, drawn
on a WebGL canvas, with two switches (output and view) and the same colour per identity as every other panel.
"""
import numpy as np
import pandas as pd

from scarlet.html_report import bundle as B
from scarlet.html_report import reasons as R

NONE = "(none)"
#: marker symbol per refusal reason in the 'unlabelled' view
REASON_SYMBOLS = {R.BACKGROUND: "circle", R.AMBIGUOUS: "diamond", R.UNDETERMINED: "x", R.FLIP: "triangle-up"}
UNLABELLED_COLOUR = "#6E6E6E"
#: the learn report's certainty: the colour's intensity from -log10 FDR, between the run's 5th and 95th
#: percentiles over its called cells (relative to the run)
FDR_PERCENTILES = (5, 95)
MIN_INTENSITY = 0.12


def text(v):
    """strings of a column read back from the bundle ('' for missing)"""
    return pd.Series(v, dtype=object).fillna("").astype(str).replace({"nan": "", "None": ""}).to_numpy(dtype=object)


def boundary_rows(best, unique, runner_up, stage):
    """Cells with a best call and no unique call, counted by best identity x runner-up identity."""
    b, u, r = text(best), text(unique), text(runner_up)
    m = (b != R.U) & (b != "") & (u == R.U)
    cols = ["stage", "best", "runner_up", "cells"]
    if not m.any():
        return pd.DataFrame(columns=cols)
    df = pd.DataFrame({"best": b[m], "runner_up": np.where(r[m] == "", NONE, r[m])})
    g = df.groupby(["best", "runner_up"]).size().reset_index(name="cells")
    g.insert(0, "stage", stage)
    return g.sort_values(["cells", "best", "runner_up"], ascending=[False, True, True])[cols]


def boundary_fig(rows, stage, title, ids):
    """Heatmap best x runner-up (cells) of one stage; None when the stage has no boundary cell."""
    import plotly.graph_objects as go
    d = rows[rows.stage.astype(str) == stage] if len(rows) else rows
    if not len(d):
        return None
    d = d.assign(best=text(d.best), runner_up=text(d.runner_up))
    order = list(map(str, ids))
    rws = [k for k in order if k in set(d.best)] + sorted(set(d.best) - set(order))
    cls = [k for k in order if k in set(d.runner_up)] + sorted(set(d.runner_up) - set(order) - {NONE}) + \
        ([NONE] if NONE in set(d.runner_up) else [])
    Z = d.pivot_table(index="best", columns="runner_up", values="cells", aggfunc="sum").reindex(index=rws, columns=cls)
    z = Z.fillna(0).astype(int).to_numpy()
    fig = go.Figure(go.Heatmap(z=z, x=cls, y=rws, colorscale="Viridis", zmin=0,
                               hovertemplate="best %{y}<br>runner-up %{x}<br>%{z} cells<extra></extra>",
                               colorbar=dict(title=dict(text="cells"))))
    fig.update_layout(title=dict(text=title, font=dict(size=13)), height=max(320, 22 * len(rws) + 170),
                      margin=dict(l=10, r=10, t=40, b=130), template="plotly_white",
                      xaxis=dict(title=dict(text="runner-up identity"), tickangle=-45, tickfont=dict(size=10)),
                      yaxis=dict(title=dict(text="best call"), autorange="reversed", tickfont=dict(size=10)))
    return fig


def intensity_from_agreement(k, n):
    """The colour intensity of a call agreed by k of n draws; full intensity with one draw (quick mode: no
    agreement)."""
    k = pd.to_numeric(pd.Series(k), errors="coerce").to_numpy(float)
    return np.clip(np.where(np.isfinite(k), k / float(n), 0.0), MIN_INTENSITY, 1.0) if n and n > 1 else np.ones(len(k))


def intensity_from_need(k, n, k_need):
    """The 'agreement' view (transfer): a call's colour from how many of the n draws agreed, spread over the range a
    final call can have: faintest at k_need (the least a final call needs, 8 of 10), full at n of n; below k_need
    (first-pass calls) faintest too. Full intensity with one draw (quick mode: no agreement)."""
    k = pd.to_numeric(pd.Series(k), errors="coerce").to_numpy(float)
    if not n or n <= 1:
        return np.ones(len(k))
    span = float(n - k_need) if n > k_need else 1.0
    v = np.where(np.isfinite(k), (k - k_need) / span, 0.0)
    return MIN_INTENSITY + (1 - MIN_INTENSITY) * np.clip(v, 0.0, 1.0)


def intensity_ambiguous(best, unique):
    """The 'ambiguity' view: full colour for the boundary cells (a best call and no unique call), the others faint."""
    b, u = text(best), text(unique)
    m = (b != R.U) & (b != "") & (u == R.U)
    return np.where(m, 1.0, MIN_INTENSITY)


def intensity_from_margin(margin, called):
    """(intensity per cell, the margin at its faintest and at full intensity), for the 'margin' view: the identity's
    colour, darker the larger the margin, scaled between the run's 5th and 95th percentiles over its called cells
    (relative to the run, the rule of the certainty view it replaces); below the 5th percentile faintest."""
    m = pd.to_numeric(pd.Series(margin), errors="coerce").to_numpy(float)
    ref = m[np.asarray(called, bool) & np.isfinite(m)]
    if not len(ref):
        return np.full(len(m), 1.0), None, None
    lo, hi = (float(np.percentile(ref, q)) for q in FDR_PERCENTILES)
    v = (m - lo) / (hi - lo) if hi > lo else np.ones(len(m))
    v = MIN_INTENSITY + (1 - MIN_INTENSITY) * np.clip(np.where(np.isfinite(v), v, 0.0), 0.0, 1.0)
    return v, lo, hi


def intensity_from_fdr(fdr, called):
    """(intensity per cell, the -log10 FDR at its faintest and at full intensity): -log10 FDR scaled between the
    run's 5th and 95th percentiles over its called cells, so the colours show how the calls of this run differ
    (relative to the run, not an absolute scale); below the 5th percentile faintest, above the 95th full."""
    f = pd.to_numeric(pd.Series(fdr), errors="coerce").to_numpy(float)
    with np.errstate(divide="ignore", invalid="ignore"):
        lg = -np.log10(np.clip(f, 1e-300, 1.0))
    ref = lg[np.asarray(called, bool) & np.isfinite(lg)]
    if not len(ref):
        return np.full(len(lg), 1.0), None, None
    lo, hi = (float(np.percentile(ref, q)) for q in FDR_PERCENTILES)
    v = (lg - lo) / (hi - lo) if hi > lo else np.ones(len(lg))
    v = MIN_INTENSITY + (1 - MIN_INTENSITY) * np.clip(np.where(np.isfinite(v), v, 0.0), 0.0, 1.0)
    return v, lo, hi


#: without WebGL the UMAPs are drawn as SVG on at most this many cells
SVG_CELLS = 10000


def svg_subsample(labels, cap=SVG_CELLS, seed=0):
    """The cells drawn when the browser has no WebGL (SVG is slow with many points): every cell when there are at most
    `cap`; otherwise a seeded sample taken within each identity, keeping every cell of an identity smaller than the
    per-identity limit t (the largest t with sum(min(n_k, t)) <= cap), so small populations stay visible. Returns
    (sorted indices or None, t)."""
    lab = text(labels)
    n = len(lab)
    if n <= cap:
        return None, None
    sizes = pd.Series(lab).value_counts()
    lo, hi = 0, int(sizes.max())
    while lo < hi:   # the largest per-identity limit that fits
        mid = (lo + hi + 1) // 2
        if int(np.minimum(sizes.to_numpy(), mid).sum()) <= cap:
            lo = mid
        else:
            hi = mid - 1
    rng = np.random.default_rng(seed)
    keep = []
    for k in sizes.index:
        idx = np.flatnonzero(lab == k)
        keep.append(idx if len(idx) <= lo else rng.choice(idx, size=lo, replace=False))
    return sorted(int(i) for i in np.concatenate(keep)), int(lo)


def umap_spec(div, title, U, outputs, colours, *, anchor=None, truth=None, views=None, spaces=None):
    """The data of one switchable UMAP for the page's script, labels coded against one category list.

    U: the UMAP table (x, y, cell); outputs: [dict(name, label (array), shade ({view: intensity array, 0-1}, the views
    drawn in the label's colour at that intensity), reason (array or None), runner_up (array or None), hover ([(title,
    array)]))]; colours: {identity: colour}; anchor: bool array or None; truth: array or None (the reference label,
    learn report)."""
    cats = set()
    for o in outputs:
        cats |= set(text(o["label"]))
        if o.get("runner_up") is not None:
            cats |= set(text(o["runner_up"]))
    if truth is not None:
        cats |= set(text(truth))
    cats.discard("")
    cats = sorted(cats - {R.U}) + [R.U]
    code = {c: i for i, c in enumerate(cats)}
    reasons = [R.BACKGROUND, R.AMBIGUOUS, R.UNDETERMINED, R.FLIP]

    def enc(v, table):
        return [table.get(x, -1) for x in text(v)]

    outs = []
    for o in outputs:
        rs = None
        if o.get("reason") is not None:
            r = text(o["reason"])
            extra = sorted(set(r) - set(reasons) - {""})
            reasons += extra
            rs = enc(r, {x: i for i, x in enumerate(reasons)})
        hov = []
        for t, arr in o.get("hover", []):
            a = pd.Series(arr)
            if pd.api.types.is_float_dtype(a):
                a = a.map(lambda v: None if not np.isfinite(v) else float(f"{v:.3g}"))
            hov.append([t, a.astype(object).where(a.notna(), None).tolist()])
        outs.append(dict(name=o["name"], lab=enc(o["label"], code),
                         shade={v: [round(float(x), 3) for x in a] for v, a in (o.get("shade") or {}).items()}, rs=rs,
                         ru=enc(o["runner_up"], code) if o.get("runner_up") is not None else None, hover=hov))
    from scarlet.html_report.html import point_size
    svg, svg_t = svg_subsample(truth if truth is not None else outputs[0]["label"])
    # the gene spaces: the first is U.x, U.y; any other carries its own coordinates (U.x_all, U.y_all)
    sp = [dict(name=(spaces or [None])[0])]
    if spaces and len(spaces) > 1 and "x_all" in U:
        sp.append(dict(name=spaces[1], x=[round(float(v), 3) for v in U.x_all],
                       y=[round(float(v), 3) for v in U.y_all]))
    return dict(div=div, title=title, x=[round(float(v), 3) for v in U.x], y=[round(float(v), 3) for v in U.y],
                cell=list(map(str, U.cell)), cats=cats, size=point_size(len(U)), svg=svg, svg_limit=svg_t, spaces=sp,
                svg_size=point_size(len(svg)) if svg is not None else point_size(len(U)),
                col=[colours.get(c, B.REFUSED_COLOUR) if c != R.U else B.REFUSED_COLOUR for c in cats],
                refused=B.REFUSED_COLOUR, unlabelled=UNLABELLED_COLOUR, symbols=REASON_SYMBOLS, reasons=reasons,
                outputs=outs, anchor=None if anchor is None else [bool(a) for a in anchor],
                truth=None if truth is None else enc(truth, code),
                views=views or (["label"] + (["anchors"] if anchor is not None else []) + ["unlabelled"]))


def umap_controls(spec, highlight=False):
    """The switches of a switchable UMAP (output; view, when there are several; highlight: one identity in full
    colour, the others faint, none by default) and its canvas."""
    from scarlet.html_report.html import _esc
    o = "".join(f'<option value="{i}">{_esc(x["name"])}</option>' for i, x in enumerate(spec["outputs"]))
    v = "".join(f'<option value="{x}">{_esc(x)}</option>' for x in spec["views"])
    d = spec["div"]
    one = len(spec["views"]) == 1
    view = (f'<select id="{d}_view" style="display:none">{v}</select>' if one else
            f' &nbsp; View: <select id="{d}_view">{v}</select>')
    out = (f'Output: <select id="{d}_out">{o}</select>' if len(spec["outputs"]) > 1 else
           f'<select id="{d}_out" style="display:none">{o}</select>')
    sp = spec.get("spaces") or []
    space = ""
    if len(sp) > 1:
        so = "".join(f'<option value="{i}">{_esc(x["name"])}</option>' for i, x in enumerate(sp))
        space = f'<br>Space: <select id="{d}_space">{so}</select>'
    hl = ""
    if highlight:
        h = "".join(f'<option value="{i}">{_esc(c)}</option>' for i, c in enumerate(spec["cats"]))
        hl = f' &nbsp; Highlight: <select id="{d}_hl"><option value="-1">(none)</option>{h}</select>'
    return f'<p>{out}{view}{hl}{space}</p><div id="{d}"></div>'


def table_payload(df, file):
    """A table embedded in the page for its download button (built in the browser; no network), stored column by
    column: a column of texts with few distinct values as codes into its list of texts, any other column as its
    values. The page rebuilds the rows from it (table_rows; tsvDownload in the page's script)."""
    d = df.astype(object).where(pd.notna(df), None)
    cols = []
    for c in d.columns:
        v = d[c].tolist()
        texts = [x for x in v if x is not None]
        if texts and all(isinstance(x, str) for x in texts):
            uniq = list(dict.fromkeys(texts))
            if 2 * len(uniq) < len(v):
                ix = {t: i for i, t in enumerate(uniq)}
                cols.append(dict(dict=uniq, codes=[None if x is None else ix[x] for x in v]))
                continue
        cols.append(dict(values=v))
    return dict(cols=list(map(str, df.columns)), n=int(len(df)), columns=cols, file=file)


def table_rows(T):
    """The rows of a table_payload, as the page rebuilds them."""
    if "rows" in T:
        return T["rows"]
    cols = [[None if k is None else c["dict"][k] for k in c["codes"]] if "codes" in c else c["values"]
            for c in T["columns"]]
    return [[col[i] for col in cols] for i in range(T["n"])]


def download_button(key, n, what):
    from scarlet.html_report.html import _esc
    return (f'<p><button type="button" onclick="tsvDownload(\'{_esc(key)}\')">Download all {n:,} {_esc(what)} '
            f'(TSV)</button> <span class="note">(built in your browser from the data in this page)</span></p>')
