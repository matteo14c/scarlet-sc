"""Known labels supplied by the user (report's query_labels): the calls of a transfer or an assign run against labels
the user already has for the query (the tutorials and benchmarks). Display only: SCARLET never reads them, and the
report shows them in their own section, marked as supplied by the user.

Scores as the tutorials compute them: over the query cells whose label the model knows, the share refused and the
median over labels of the F1 on the called cells (refusals are not counted as errors: they are read beside the
share refused). Per label: the share called correctly, called as another identity, and refused; for a label the
model does not know, what its cells were called."""
import numpy as np
import pandas as pd

from scarlet.html_report import reasons as R


def labels_of(obs, column):
    """The user's labels, as text, from the query's obs column."""
    if column not in obs.columns:
        raise ValueError(f"query_labels: the query has no obs column {column!r} (its columns: "
                         f"{', '.join(map(str, obs.columns[:20]))})")
    return obs[column].astype(str).to_numpy(dtype=object)


def tables(truth, outputs, known, column):
    """(summary, per_label, confusion, meta) for the bundle. outputs: {name: calls per cell, in the cells' order}."""
    from sklearn.metrics import f1_score
    truth = np.asarray(truth, dtype=object)
    known = set(map(str, known))
    in_model = np.isin(truth, list(known))
    summ, per, conf = [], {}, []
    labels = sorted(set(truth), key=lambda x: (x not in known, x))
    for name, pred in outputs.items():
        pred = np.asarray(pred, dtype=object)
        p, t = pred[in_model], truth[in_model]
        called = p != R.U
        labs = sorted(set(t[called]))
        f1 = f1_score(t[called], p[called], labels=labs, average=None, zero_division=0) if len(labs) else []
        summ.append({"output": name, "cells": int(in_model.sum()), "pct_refused": 100 * float((~called).mean())
                     if len(p) else np.nan, "median_f1_called": float(np.median(f1)) if len(f1) else np.nan,
                     "cells_label_unknown_to_model": int((~in_model).sum()),
                     "pct_refused_label_unknown_to_model": 100 * float((pred[~in_model] == R.U).mean())
                     if (~in_model).any() else np.nan})
        for k in labels:
            m = truth == k
            pk = pred[m]
            row = per.setdefault(k, {"label": k, "known_to_model": k in known, "cells": int(m.sum())})
            row[f"pct_correct, {name}"] = 100 * float((pk == k).mean()) if k in known else np.nan
            row[f"pct_other identity, {name}"] = 100 * float(((pk != k) & (pk != R.U)).mean())
            row[f"pct_refused, {name}"] = 100 * float((pk == R.U).mean())
            vc = pd.Series(pk).value_counts()
            for c, n in vc.items():
                conf.append({"output": name, "label": k, "call": c, "cells": int(n)})
    S = pd.DataFrame(summ)
    P = pd.DataFrame(list(per.values()))
    first = next(iter(outputs))
    P[f"most frequent call, {first}"] = [pd.Series(np.asarray(outputs[first], dtype=object)[truth == k])
                                          .value_counts().index[0] for k in P.label]
    return S, P, pd.DataFrame(conf), {"column": column, "outputs": list(outputs)}


def write(bd, truth, outputs, known, column):
    S, P, C, meta = tables(truth, outputs, known, column)
    bd.table("known_summary", S, f"labels supplied by the user (obs column {column!r}): per output, over the cells "
                                 "whose label the model knows, the share refused and the median F1 per label on the "
                                 "called cells")
    bd.table("known_per_label", P, "labels supplied by the user: per label, the share called correctly, called as "
                                   "another identity, and refused")
    bd.table("known_confusion", C, "labels supplied by the user: cells per label and call")
    bd.meta("known_labels", meta)


def section(T, meta, col, ids, number, table_html, esc, id_fmt):
    """(html, figure) of the section, or ("", None) without labels."""
    if "known_summary" not in T:
        return "", None
    import plotly.graph_objects as go
    km = meta["known_labels"]
    first = km["outputs"][0]
    C = T["known_confusion"]
    C = C[C.output.astype(str) == first]
    P = T["known_per_label"]
    rows = list(map(str, P.label))
    cols = [k for k in ids if k in set(C.call.astype(str))] + ([R.U] if R.U in set(C.call.astype(str)) else [])
    M = C.assign(label=C.label.astype(str), call=C.call.astype(str)).pivot_table(
        index="label", columns="call", values="cells", aggfunc="sum").reindex(index=rows, columns=cols).fillna(0)
    share = M.div(M.sum(axis=1).replace(0, np.nan), axis=0).fillna(0)
    fig = go.Figure(go.Heatmap(z=share.to_numpy().tolist(), x=cols, y=rows, colorscale="Blues", zmin=0, zmax=1,
                               customdata=M.to_numpy().astype(int).tolist(),
                               hovertemplate="label %{y}<br>called %{x}<br>%{customdata} cells (%{z:.1%} of the "
                                             "label)<extra></extra>", colorbar=dict(title=dict(text="share"))))
    fig.update_layout(height=max(320, 26 * len(rows) + 170), template="plotly_white",
                      title=dict(text=f"Labels supplied by the user x {first} (share of each label's cells)", x=0.01,
                                 xanchor="left", font=dict(size=13)),
                      xaxis=dict(title="SCARLET's call", tickangle=-45), yaxis=dict(title="label supplied by the user",
                                                                                  autorange="reversed"),
                      margin=dict(l=10, r=10, t=50, b=140))
    html = (f"<h2>{number}. Against labels supplied by the user</h2>"
            f"<p class='card'><b>Labels supplied by the user</b> (obs column {esc(km['column'])!s}), not by SCARLET: "
            "SCARLET did not read them; this section only compares its calls with them.</p>"
            "<p class='note'>Over the cells whose label the model knows: the share refused, and the median over labels "
            "of the F1 on the called cells (a refusal is not counted as an error; read the two together). Labels "
            "the model does not know are listed last, with what their cells were called.</p>"
            + table_html(T["known_summary"])
            + table_html(P, fmt={"label": id_fmt.get("identity")} if id_fmt else None)
            + '<div id="known_conf"></div>')
    return html, fig
