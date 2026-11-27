"""The settings a learn uses, resolved once and recorded.

`resolve()` turns the options (technology preset, blocklist and its gene groups, imbalanced and the marker sizes)
into the settings a run uses, with the same rules the marker search applies:
  - technology preset and marker sizes: scarlet._engine.configure_marker_defaults (TECHNOLOGIES);
  - blocklist gene groups: scarlet.blocklist_modules.resolve_modules on the labels being learned;
  - imbalanced: scarlet.marker_finder.resolve_imbalanced_profile on the labels being learned.
It only reports: the marker search receives the options as given and decides again on the same inputs, so a run is
the same with or without this step. `preflight_text()` is the short summary every command prints before it runs;
the resolved dict is recorded in the model file (settings.json) and in provenance.
"""

from __future__ import annotations

import argparse
import os

import numpy as np
import pandas as pd

TECHNOLOGY_CHOICES = ("sparse", "dense")
IMBALANCED_CHOICES = ("auto", "yes", "no")
HIGH_DENS_MIN_GENES = 5000  # the dense rule of technology auto (mean detected genes per cell)


def default_blocklist():
    from scarlet.marker_finder import DEFAULT_BLOCKLIST_PATH
    return str(DEFAULT_BLOCKLIST_PATH)


def resolve(*, technology="sparse", blocklist=None, blocklist_modules="auto",
            exclude_blocklist_categories=None, include_blocklist_categories=None, imbalanced="auto",
            n_markers=None, marker_shortlist=None, marker_gamma=None, universe="compute",
            marker_alpha=1.5, marker_alpha_rule=None, classifier_bg_mode=None, labels=None, gene_names=None):
    """The settings a run uses. `labels` (the labels being learned) decide the gene groups and the class-size
    handling; without them those two are reported as decided per run."""
    from scarlet import _engine
    from scarlet.blocklists import blocklist_identity, symbol_case_species

    if technology not in TECHNOLOGY_CHOICES:
        raise ValueError(f"technology must be one of {TECHNOLOGY_CHOICES} (got {technology!r})")
    if imbalanced not in IMBALANCED_CHOICES:
        raise ValueError(f"imbalanced must be one of {IMBALANCED_CHOICES} (got {imbalanced!r})")
    ns = argparse.Namespace(technology=technology, marker_gamma=marker_gamma, top_final=n_markers,
                            background_genes_N=marker_shortlist)
    _engine.configure_marker_defaults(ns)
    out = {
        "schema": "scarlet.settings/3",
        "technology": {"name": ns.technology, "detection": "count > 0", "marker_gamma": float(ns.marker_gamma)},
        "sizes": {"from": ns.sizes_from, "n_markers": int(ns.top_final), "marker_shortlist": int(ns.background_genes_N),
                  "explicit": {k: v for k, v in (("n_markers", n_markers), ("marker_shortlist", marker_shortlist))
                               if v is not None}},
        "universe": "computed from the input cells" if universe == "compute" else f"external: {universe}",
        "marker_alpha": float(marker_alpha),
    }
    # auto / human / mouse / none / a file (scarlet.blocklists.choose_blocklist); auto without the gene
    # names (a record made without data) is reported as decided from the data
    from scarlet.blocklists import choose_blocklist
    if (blocklist is None or str(blocklist) == "auto") and gene_names is None:
        chosen, how = None, dict(option="auto", choice=None, basis="auto: chosen from the data's gene names")
    else:
        chosen = choose_blocklist(blocklist, gene_names)
        how = chosen.record
    path = str(chosen) if chosen is not None else None
    ident = blocklist_identity(path)
    bl = {"path": path, "given": how["option"] != "auto", "option": how["option"], "choice": how["choice"],
          "basis": how["basis"], "version": ident.get("version"),
          "species": ident.get("species"), "md5": ident.get("md5"), "n_rows": ident.get("n_rows"),
          "modules_option": blocklist_modules, "exclude_categories": exclude_blocklist_categories,
          "include_categories": include_blocklist_categories, "modules": None}
    if labels is not None and ident.get("md5"):
        from scarlet.blocklist_modules import resolve_modules
        from scarlet.blocklists import read_table
        cats = read_table(path)
        cats = cats["Category"].dropna().astype(str).unique() if "Category" in cats.columns else []
        names = pd.Series(np.asarray(labels)).dropna().astype(str).unique().tolist()
        rec = resolve_modules(cats, names, blocklist_modules, exclude_blocklist_categories,
                              include_blocklist_categories)["record"]
        bl["modules"] = {k: rec.get(k) for k in ("modular_list", "mode", "on", "off", "cell_cycle_off_by", "override")}
    if gene_names is not None:
        sp, frac, n = symbol_case_species(gene_names)
        bl["data_species"] = sp
        bl["species_check"] = ("not checked" if sp is None or ident.get("species") is None
                               else ("ok" if sp == ident.get("species") else "mismatch"))
    out["blocklist"] = bl
    imb = {"option": imbalanced, "handling": None, "f_max": None, "threshold": None,
           "marker_alpha_rule": marker_alpha_rule, "classifier_bg_mode": classifier_bg_mode}
    if labels is not None:
        from scarlet.marker_finder import resolve_imbalanced_profile
        probe = argparse.Namespace(imbalanced=imbalanced, alpha_rule=marker_alpha_rule,
                                   classifier_bg_mode=classifier_bg_mode,
                                   external_universe=None if universe == "compute" else universe)
        resolve_imbalanced_profile(probe, labels)
        r = probe._imbalanced_record
        imb.update(handling=r["profile"], f_max=r["f_max"], threshold=r["threshold"], n_cells=r["n_cells"],
                   n_classes=r["n_classes"], alpha_rule=r["alpha_rule"], classifier_bg_mode=r["classifier_bg_mode"])
    out["imbalanced"] = imb
    return out


def data_summary(path_or_adata, labels=None, gene_col="index"):
    """Cells, genes, labels and density (mean detected genes per cell; median total counts) of an input, read
    without loading the matrix into memory when a path is given."""
    import anndata as ad
    if isinstance(path_or_adata, (str, os.PathLike)):
        A = ad.read_h5ad(path_or_adata, backed="r")
        try:
            obs, var = A.obs, A.var
            detected, totals = _row_stats_backed(str(path_or_adata), A.n_obs)
        finally:
            A.file.close()
    else:
        A = path_or_adata
        obs, var = A.obs, A.var
        detected, totals = _row_stats(A.X)
    genes = var.index if gene_col == "index" else var[gene_col]
    out = {"n_cells": int(len(obs)), "n_genes": int(len(var)),
           "mean_detected_genes": float(np.mean(detected)) if len(detected) else None,
           "median_total_counts": float(np.median(totals)) if len(totals) else None,
           "gene_names": pd.Index(np.asarray(genes).astype(str))}
    if labels is not None:
        if labels not in obs.columns:
            raise KeyError(f"label column {labels!r} is not an obs column of the input")
        lab = obs[labels]
        out["labels"] = lab.to_numpy()
        out["n_labels"] = int(lab.dropna().nunique())
    return out


# Row passes in chunks of whole rows, capped by rows and by stored entries: a chunk never splits a row, so each row's
# values are summed by the same bincount in the same order whatever the chunk boundaries; the cap on entries bounds
# the row index (int64) and the float64 weights built per chunk (8 bytes per stored value each), which a cap on rows
# alone does not on deep data (20,000 rows of 10,000 genes: 200 M entries).
ROW_CHUNK_ENTRIES = 1 << 22


def _row_chunks(indptr, max_rows=20000, max_entries=ROW_CHUNK_ENTRIES):
    """(i, j) row ranges covering all rows: at most max_rows rows and max_entries stored entries each (a longer row
    alone)."""
    indptr = np.asarray(indptr)
    n = len(indptr) - 1
    i = 0
    while i < n:
        j = int(np.searchsorted(indptr, indptr[i] + max_entries, side="right")) - 1
        j = max(i + 1, min(j, i + max_rows, n))
        yield i, j
        i = j


def _row_stats(X, chunk=20000):
    import scipy.sparse as sp
    if sp.issparse(X):
        X = sp.csr_matrix(X)
        # Memory: in chunks of whole rows (_row_chunks; as _row_stats_backed): each row lies in one chunk and its
        # values are summed in the same order, so the result is identical; the whole-matrix row index (int64) and
        # float64 weights (8 bytes per stored value each) are not built at once.
        n = X.shape[0]
        detected = np.zeros(n)
        totals = np.zeros(n)
        for i, j in _row_chunks(X.indptr, chunk):
            v = X.data[X.indptr[i]:X.indptr[j]]
            rows = np.repeat(np.arange(j - i), np.diff(X.indptr[i:j + 1]))
            detected[i:j] = np.bincount(rows, weights=(v > 0), minlength=j - i)
            totals[i:j] = np.bincount(rows, weights=v.astype(np.float64), minlength=j - i)
    else:
        X = np.asarray(X)
        detected = (X > 0).sum(axis=1)
        totals = X.sum(axis=1)
    return np.asarray(detected, dtype=float), np.asarray(totals, dtype=float)


def detected_per_cell(X, chunk=20000):
    """Genes detected per cell: _row_stats(X)[0], the same floats, without the totals. For a CSR matrix whose stored
    values are all positive (counts without stored zeros, or a detection pattern) it is the number of stored values
    per row (1.0 added once per stored value > 0 in _row_stats: the same integers)."""
    import scipy.sparse as sp
    if sp.issparse(X) and X.format == "csr" and (X.nnz == 0 or X.data.min() > 0):
        return np.diff(X.indptr).astype(float)
    return _row_stats(X, chunk)[0]


def _row_stats_backed(path, n_obs, chunk=20000):
    import h5py
    with h5py.File(path, "r") as f:
        X = f["X"]
        enc = X.attrs.get("encoding-type", "")
        enc = enc.decode() if isinstance(enc, bytes) else enc
        if isinstance(X, h5py.Group) and enc == "csr_matrix":
            indptr = X["indptr"][:]
            data = X["data"]
            detected = np.zeros(n_obs)
            totals = np.zeros(n_obs)
            for i, j in _row_chunks(indptr, chunk):
                v = data[indptr[i]:indptr[j]].astype(np.float64)
                rows = np.repeat(np.arange(j - i), np.diff(indptr[i:j + 1]))
                detected[i:j] = np.bincount(rows, weights=(v > 0), minlength=j - i)
                totals[i:j] = np.bincount(rows, weights=v, minlength=j - i)
            return detected, totals
    import anndata as ad
    A = ad.read_h5ad(path)
    return _row_stats(A.X)


def blocklist_text(b):
    """The blocklist in plain words: how it was chosen."""
    return b.get("basis") or "not recorded"


def preflight_text(command, settings, data=None, extra=None):
    """The short pre-flight summary printed before a command runs."""
    t, b, i = settings["technology"], settings["blocklist"], settings["imbalanced"]
    p = settings["sizes"]
    lines = [f"[SCARLET pre-flight] {command}"]
    if data is not None:
        d = f"{data['n_cells']:,} cells x {data['n_genes']:,} genes"
        if data.get("n_labels") is not None:
            d += f", {data['n_labels']} labels"
        if data.get("mean_detected_genes") is not None:
            d += f"; mean {data['mean_detected_genes']:,.0f} detected genes per cell"
        lines.append(f"  data         {d}")
    tech = (f"{t['name']} (a gene is detected when its count is above 0; background penalty "
            f"{t['marker_gamma']:g})")
    lines.append(f"  technology   {tech}")
    pr = f"{p['n_markers']} markers per cell type, shortlist {p['marker_shortlist']}"
    pr += f" ({p['from']})"
    pr += f"; given: {', '.join(f'{k} {v}' for k, v in p['explicit'].items())}" if p["explicit"] else ""
    lines.append(f"  sizes        {pr}")
    blk = blocklist_text(b)
    m = b.get("modules")
    if m and m.get("modular_list"):
        blk += f"; gene groups ({m['mode']}) applied: {', '.join(m['on'] or []) or 'none'}; not applied: " \
               f"{', '.join(m['off'] or []) or 'none'}"
        if m.get("cell_cycle_off_by"):
            blk += f" (cell-cycle genes usable: {', '.join(m['cell_cycle_off_by'][:3])} name dividing cells)"
    elif m is not None:
        blk += "; the list has no gene groups"
    else:
        blk += f"; gene groups {b['modules_option']} (decided on each learn's labels)"
    if b.get("species_check") == "mismatch":
        blk += f"; WARNING: the data look {b['data_species']}, the list is {b['species']}"
    lines.append(f"  blocklist    {blk}")
    if i["handling"] is None:
        im = f"{i['option']} (decided on each learn's labels: balanced when the largest class holds more than a third " \
             "of the cells)"
    else:
        im = f"{i['option']} -> {i['handling']} (the largest class holds {i['f_max']:.0%} of {i['n_cells']:,} cells; " \
             "balanced above a third)"
    lines.append(f"  imbalanced   {im}")
    for k, v in (extra or {}).items():
        lines.append(f"  {k:<12} {v}")
    return "\n".join(lines)


def transfer_preflight(ref, query, technology, depth_match="none"):
    """Density and gene overlap of reference and query, with a plain recommendation (transfer pre-flight).
    depth_match: the transfer's depth matching (capture or none), for the depth advice."""
    shared = ref["gene_names"].intersection(query["gene_names"])
    o_ref = len(shared) / max(1, ref["n_genes"])
    o_q = len(shared) / max(1, query["n_genes"])
    dr, dq = ref["mean_detected_genes"], query["mean_detected_genes"]
    tr, tq = ref.get("median_total_counts"), query.get("median_total_counts")
    rec = []
    if len(shared) == 0:
        rec.append("no shared genes: check the gene identifiers (symbols vs Ensembl ids) or --gene-col")
    elif min(o_ref, o_q) < 0.5:
        rec.append(f"low gene overlap ({len(shared):,} shared genes): check the gene identifiers; markers outside the "
                   "shared genes cannot be scored")
    ratio = (max(dr, dq) / max(1.0, min(dr, dq))) if (dr and dq) else 1.0
    if ratio >= 1.5:
        deeper = "reference" if dr > dq else "query"
        if depth_match == "none":
            rec.append(f"the {deeper} detects {ratio:.1f}x more genes per cell: consider depth matching "
                       "(--depth-match yes, the default)")
        else:
            rec.append(f"the {deeper} detects {ratio:.1f}x more genes per cell: depth matching thins the {deeper} "
                       "in the first pass")
    if technology == "sparse" and dr and dq and min(dr, dq) >= HIGH_DENS_MIN_GENES:
        rec.append(f"both datasets detect >= {HIGH_DENS_MIN_GENES} genes per cell (plate-based): consider "
                   "--technology dense")
    if technology == "dense" and dr and dq and min(dr, dq) < HIGH_DENS_MIN_GENES:
        rec.append(f"the dense preset expects >= {HIGH_DENS_MIN_GENES} detected genes per cell; with fewer, "
                   "--technology sparse is the usual choice")
    if not rec:
        rec.append("no concerns: densities and gene overlap are compatible")
    return {"shared_genes": int(len(shared)), "share_of_reference_genes": round(o_ref, 4),
            "share_of_query_genes": round(o_q, 4), "reference_mean_detected_genes": dr,
            "query_mean_detected_genes": dq, "reference_median_total_counts": tr, "query_median_total_counts": tq,
            "density_ratio": round(ratio, 3), "recommendation": rec}
