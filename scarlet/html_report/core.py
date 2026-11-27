"""scarlet.report: the report of a finished run. Read-only over the run's outputs; it changes no call."""
import os

from scarlet.html_report import bundle as B


REPORT_EXTRA = "pip install 'scarlet-sc[report]'"


def require_plotly():
    """The HTML reports need plotly (the package's `report` extra): checked before any work, so a run with
    report=True does not stop at its very end."""
    try:
        import plotly  # noqa: F401
    except ImportError as e:
        raise ImportError(f"the HTML report needs plotly, the 'report' extra of SCARLET: {REPORT_EXTRA}") from e


def kind_of(run):
    """learn (a model file), transfer (a transfer output folder) or assign (an annotated .h5ad)"""
    p = str(run)
    if os.path.isdir(p) and os.path.exists(os.path.join(p, "transfer_record.json")):
        return "transfer"
    if os.path.isfile(p) and p.endswith(".scarlet"):
        return "learn"
    if os.path.isfile(p) and p.endswith(".h5ad"):
        return "assign"
    raise ValueError(f"not a finished SCARLET run: {p!r} (expected a model file .scarlet, a transfer output folder "
                     "with transfer_record.json, or an annotated .h5ad from assign)")


def default_output(run, kind):
    p = str(run).rstrip("/")
    if kind == "transfer":
        return os.path.join(p, "report")
    stem = p[:-len(".scarlet")] if p.endswith(".scarlet") else p[:-len(".h5ad")] if p.endswith(".h5ad") else p
    return stem + "_report"


def report(run, *, reference=None, query=None, model=None, query_labels=None, output=None, self_assessment=True,
           n_folds=5, seed=0, umap_cells=50000, umap_all_genes=True, alpha=0.01, ncpu=1, log=print):
    require_plotly()
    kind = kind_of(run)
    if query_labels is not None and kind == "learn":
        raise ValueError("query_labels applies to a transfer or an assign report (a learn report has the reference's "
                         "own labels)")
    out = output or default_output(run, kind)
    os.makedirs(out, exist_ok=True)
    if kind == "learn":
        from scarlet.html_report import html, learn_report
        learn_report.build(run, reference, out, self_assessment=self_assessment, n_folds=n_folds,
                           seed=seed, umap_cells=umap_cells, alpha=alpha, ncpu=ncpu, log=log,
                           umap_all_genes=umap_all_genes)
        page = html.render_learn(out)
    elif kind == "transfer":
        from scarlet.html_report import html, transfer_report
        transfer_report.build(run, reference, query, out, self_assessment=self_assessment,
                              n_folds=n_folds, seed=seed, umap_cells=umap_cells, ncpu=ncpu, log=log,
                              query_labels=query_labels, umap_all_genes=umap_all_genes)
        sup = os.path.join(out, "supplement_A")
        with open(os.path.join(sup, "report.html"), "w", encoding="utf-8") as fh:
            fh.write(html.render_learn(sup, back_link="../report.html"))
        page = html.render_transfer(out)
        with open(os.path.join(out, html.REFINEMENT_PAGE), "w", encoding="utf-8") as fh:   # section 6's own page
            fh.write(html.render_refinement(out))
    else:
        from scarlet.html_report import assign_report, html
        assign_report.build(run, model, out, alpha=alpha, seed=seed, umap_cells=umap_cells, log=log,
                            query_labels=query_labels, umap_all_genes=umap_all_genes)
        page = html.render_assign(out)
    path = os.path.join(out, "report.html")
    with open(path, "w", encoding="utf-8") as fh:
        fh.write(page)
    man, _ = B.load(out)
    log(f"[report] {kind} report: {path} ({os.path.getsize(path) / 1e6:.1f} MB); bundle {os.path.join(out, 'bundle')}")
    return path
