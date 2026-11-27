"""The modular blocklists U4 (human) and M4 (mouse) and the module rule (scarlet/blocklist_modules.py).

1. symbol / id consistency of U4 (CELLxGENE census, HPA and HGNC gene tables) and M4 (MGI MRK_List2 / MRK_ENSEMBL);
2. U4/M4 are the bundled lists (blocklist human / mouse);
3. the rule reproduces its reference table (tests/data/blocklist_u4_module_rule_per_dataset.tsv) exactly;
4. the existing category flags switch exactly each module's genes in the marker finder.
5. the rule is applied per marker-finder run on its own labels, with a user override (--blocklist-modules), and
   recorded with the --imbalanced preset.
"""
import json
import os
import re
from argparse import Namespace

import anndata as ad
import numpy as np
import pandas as pd
import pytest

from scarlet import blocklist_modules as bm
from scarlet import blocklists as bl
from scarlet import marker_finder as mf

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
DATA = os.path.join(ROOT, "tests", "data")
U4 = os.path.join(ROOT, "scarlet", "scarlet_blocklist_universal.tsv")      # the bundled human list
M4 = os.path.join(ROOT, "scarlet", "scarlet_blocklist_mouse.tsv")          # the bundled mouse list
PRIMARY_HUMAN = ("cellxgene_census_2025-11-08", "hpa_rna_single_cell_cluster_2025-12-12")


def read(path):
    return pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)


# ------------------------------------------------------------------------------------------------ 1. consistency
U4_COUNTS = {"core_stress": 21, "core_ribosomal": 265, "core_mitochondrial": 49, "core_sex": 8, "cell_cycle": 97, "stress_extended": 14}


def test_u4_structure_and_counts():
    u = read(U4)
    assert list(u.columns) == ["ENSG_ID", "Symbol", "Category"]
    assert tuple(dict.fromkeys(u.Category)) == bm.MODULES
    assert u.Category.value_counts().to_dict() == U4_COUNTS
    assert not u.Symbol.str.upper().duplicated().any()
    ids = u.ENSG_ID[u.ENSG_ID != ""]
    assert not ids.duplicated().any() and ids.str.match(r"^ENSG\d{11}$").all()
    assert sorted(u.loc[u.ENSG_ID == "", "Symbol"]) == sorted(f"MTRNR2L{i}" for i in range(1, 13))
    assert bl.file_md5(U4) == "3b7fddd9cb2dd8667dce9f9c63a6d938"   # U4 as the package identifies it (blocklists.py)
    # the category decisions of 2026-10-07
    cat = dict(zip(u.Symbol, u.Category))
    assert cat["RHOB"] == "stress_extended" and cat["JUND"] == "core_stress" and cat["BAG3"] == "core_stress"
    assert cat["MT-TS1"] == "core_mitochondrial"
    for g in ("HMOX1", "HSPA5", "HSPA8", "HSP90AB1", "FOSL1", "RPL3L", "RPS6KA1", "RPS6KL1"):
        assert g not in cat


PKG_U = os.path.join(ROOT, "scarlet", "scarlet_blocklist_universal.tsv")
PKG_M = os.path.join(ROOT, "scarlet", "scarlet_blocklist_mouse.tsv")


@pytest.mark.parametrize("path", [U4, PKG_U], ids=["U4", "default-package"])
def test_u4_ids_belong_to_their_symbols(path):
    ref = pd.read_csv(os.path.join(DATA, "blocklist_u4_id_reference.tsv"), sep="\t", comment="#", dtype=str, keep_default_na=False)
    pairs = set(zip(ref.gene_id, ref.symbol))
    prim = ref[ref.source.isin(PRIMARY_HUMAN)]
    primary_ids = prim.groupby("symbol").gene_id.apply(set).to_dict()
    problems = []
    for sym, gid in read(path)[["Symbol", "ENSG_ID"]].itertuples(index=False):
        if gid and (gid, sym) not in pairs:
            problems.append(f"{sym}: {gid} is not {sym} in the gene tables")
        elif sym in primary_ids and gid not in primary_ids[sym]:
            problems.append(f"{sym}: id {gid or '(empty)'}, gene tables give {sorted(primary_ids[sym])}")
    assert problems == []
    # six ribosomal pseudogenes (e.g. RPL18P3) have their id from org.Hs.eg.db only; that source is accepted (as for U3.1)
    src = ref.groupby(["gene_id", "symbol"]).source.apply(set)
    only_orghs = [sym for sym, gid in read(path)[["Symbol", "ENSG_ID"]].itertuples(index=False)
                  if gid and src.get((gid, sym)) == {"org.Hs.eg.db_3.18.0"}]
    assert len(only_orghs) == 6 and "RPL18P3" in only_orghs


@pytest.mark.parametrize("path", [M4, PKG_M], ids=["M4", "default-package"])
def test_m4_ids_and_symbols_are_mgi(path):
    ref = pd.read_csv(os.path.join(DATA, "blocklist_m4_mgi_reference.tsv"), sep="\t", comment="#", dtype=str, keep_default_na=False)
    off = ref[ref.status == "official"]
    ens_of = {s: set(e.split(";")) - {""} for s, e in zip(off.symbol, off.ensmusg)}
    first_of = {s: (e.split(";")[0] if e else "") for s, e in zip(off.symbol, off.ensmusg)}
    withdrawn = set(ref.loc[ref.status == "withdrawn", "symbol"])
    m = read(path)
    assert list(m.columns) == ["ENSG_ID", "Symbol", "Category"]
    assert tuple(dict.fromkeys(m.Category)) == bm.MODULES
    assert not m.Symbol.map(mf.canonical_symbol_key).duplicated().any()
    ids = m.ENSG_ID[m.ENSG_ID != ""]
    assert not ids.duplicated().any() and ids.str.match(r"^ENSMUSG\d{11}$").all()
    problems = []
    for gid, sym in zip(m.ENSG_ID, m.Symbol):
        if sym in ens_of:
            if gid and gid not in ens_of[sym]:
                problems.append(f"{sym}: {gid} is not an MGI id of {sym}")
            if not gid and first_of[sym]:
                problems.append(f"{sym}: MGI gives {first_of[sym]}, the row has none")
        elif sym in withdrawn:
            if gid:
                problems.append(f"{sym}: alias rows carry no id")
        else:
            problems.append(f"{sym}: not an MGI mouse symbol")
    assert problems == []


def test_m4_mirrors_the_u4_modules():
    m = read(M4)
    mod_of = dict(zip(m.Symbol, m.Category))
    mp = pd.read_csv(os.path.join(DATA, "blocklist_u4_to_m4.tsv"), sep="\t", comment="#", dtype=str, keep_default_na=False)
    wrong = [(h, c, ms) for h, c, ms in zip(mp.human_symbol, mp.Category, mp.mouse_symbols)
             for x in ms.split(";") if x and mod_of.get(x) != c]
    assert wrong == []
    # the special-care rows (U4_mouse_special_care.tsv)
    row = mp.set_index("human_symbol")
    assert row.loc["HSPA6", "mouse_symbols"] == ""                      # no mouse orthologue
    assert row.loc["HSP90AA1", "mouse_symbols"] == "Hsp90aa1"           # Hsp90ab1 stays out (HSP90AB1 dropped)
    for h, mo in (("CENPU", "Cenpu"), ("PIMREG", "Pimreg"), ("JPT1", "Jpt1"), ("MT-TS1", "mt-Ts1")):
        assert row.loc[h, "mouse_symbols"] == mo and row.loc[h, "mapping_source"] == "MGI HOM"
    assert "Hsp90ab1" not in mod_of and "Hmox1" not in mod_of and "Rpl3l" not in mod_of and "Rps6ka1" not in mod_of
    assert mod_of["Eif2s3y"] == "core_sex" and mod_of["mt-Ts1"] == "core_mitochondrial"


# ------------------------------------------------------------------------------------------------ 2. defaults


def test_defaults_are_u4_and_m4():
    """The bundled human and mouse lists, recorded by version in provenance."""
    for p in (mf.DEFAULT_BLOCKLIST_PATH, PKG_U):
        assert bl.blocklist_identity(str(p))["version"] == "U4" and bl.file_md5(str(p)) == "3b7fddd9cb2dd8667dce9f9c63a6d938"
    for p in (PKG_M,):
        assert bl.blocklist_identity(p)["version"] == "M4" and bl.file_md5(p) == "d15c550f4fd8031325c5c97b7a01b745"
    # the marker finder's parser defaults to it; the commands default to 'auto' (0.6.3: U4 for human-style gene
    # names, M4 for mouse-style ones)
    assert mf.build_parser().get_default("blocklist") == str(mf.DEFAULT_BLOCKLIST_PATH)
    assert mf.build_parser().get_default("blocklist_modules") == "auto"
    import inspect
    from scarlet import api, settings
    assert settings.default_blocklist() == str(mf.DEFAULT_BLOCKLIST_PATH)
    for fn in (api.learn, api.transfer):
        ps = inspect.signature(fn).parameters
        assert ps["blocklist"].default == "auto" and ps["blocklist_modules"].default == "auto"
    assert bl.blocklist_identity(bl.choose_blocklist("auto", ["CD3E"] * 30))["version"] == "U4"
    assert bl.blocklist_identity(bl.choose_blocklist("auto", ["Cd3e"] * 30))["version"] == "M4"


# ------------------------------------------------------------------------------------------------ 3. the rule
def test_rule_reproduces_the_design_table():
    labels = pd.read_csv(os.path.join(DATA, "blocklist_u4_rule_labels.tsv"), sep="\t", dtype=str, keep_default_na=False)
    table = pd.read_csv(os.path.join(DATA, "blocklist_u4_module_rule_per_dataset.tsv"), sep="\t", dtype=str, keep_default_na=False)
    assert len(labels) == len(table) == 110
    for (_, lab), (_, row) in zip(labels.iterrows(), table.iterrows()):
        assert (lab.dataset, lab.label_set) == (row.dataset, row.label_set)
        res = bm.module_set(json.loads(lab.labels_json))
        st = res["states"]
        got = dict(cell_cycle=st["cell_cycle"], matched_classes="; ".join(res["cell_cycle_off_by"]), stress_core=st["core_stress"],
                   ribosomal=st["core_ribosomal"], mitochondrial=st["core_mitochondrial"], sex=st["core_sex"], stress_extended=st["stress_extended"])
        want = {k: row[k] for k in got}
        assert got == want, (row.dataset, row.label_set)


def test_rule_basic_cases():
    assert bm.module_set(["alpha", "beta"])["off"] == ["stress_extended"]
    r = bm.module_set(["T cells", "Cycling T&NK", "fibro-adipogenic_progenitors"])
    assert r["states"]["cell_cycle"] == "OFF" and r["cell_cycle_off_by"] == ["Cycling T&NK"]
    assert bm.module_set(["oligodendrocyte_progenitor_cells", "mesenchymal stem cell"])["states"]["cell_cycle"] == "ON"
    assert bm.excluded_categories(["Plasmablasts"]) == "cell_cycle,stress_extended"
    assert bm.excluded_categories(["B cells"]) == "stress_extended"


@pytest.mark.parametrize("opc", ["oligodendrocyte_progenitor", "oligodendrocyte progenitor cell",
                                 "oligodendrocyte precursor cell"])
def test_rule_opc_spellings_keep_cell_cycle_on(opc):
    # the rule names oligodendrocyte progenitors as a quiescent exception, under every spelling
    # (the Cell Ontology name "oligodendrocyte precursor cell" was matched by the bare term "precursor")
    for name in (opc, opc.capitalize(), opc.upper()):
        r = bm.module_set(["astrocyte", name, "neuron"])
        assert r["states"]["cell_cycle"] == "ON" and r["cell_cycle_off_by"] == [], name
    # a proliferative label next to it still turns cell_cycle off, and only that label is the cause
    r = bm.module_set(["astrocyte", opc, "cycling cell"])
    assert r["states"]["cell_cycle"] == "OFF" and r["cell_cycle_off_by"] == ["cycling cell"]
    # other precursors are unaffected
    assert bm.module_set(["neural precursor cell"])["cell_cycle_off_by"] == ["neural precursor cell"]


@pytest.mark.parametrize("label", ["cycling oligodendrocyte precursor cell", "cycling oligodendrocyte_progenitor",
                                   "proliferating oligodendrocyte progenitor cell", "Dividing OPC-like oligodendrocyte_progenitor",
                                   "mitotic oligodendrocyte precursor cell", "cycling mesenchymal stem cell",
                                   "proliferating fibro-adipogenic_progenitors"])
def test_rule_cycling_terms_override_the_excludes(label):
    # an exclude neutralises only the generic stage terms; an explicit cycling term still turns cell_cycle off
    r = bm.module_set(["astrocyte", label])
    assert r["states"]["cell_cycle"] == "OFF" and r["cell_cycle_off_by"] == [label], label


def test_rule_state_terms_are_rule_terms():
    assert set(bm.CC_OFF_STATE_TERMS) <= set(bm.CC_OFF_TERMS)
    # the excluded lineage names alone keep cell_cycle on
    names = ["fibro-adipogenic_progenitors", "oligodendrocyte_progenitor_cells", "oligodendrocyte precursor cell",
             "mesenchymal stem cell", "pituitary_stem_cells", "fallopian_basal_stem_cells", "hepatic_stellate_cells"]
    assert bm.module_set(names)["states"]["cell_cycle"] == "ON"

# ------------------------------------------------------------------------------------------------ 4. category flags
def _adata(names):
    return ad.AnnData(X=np.ones((2, len(names)), dtype=np.float32), var=pd.DataFrame(index=pd.Index(names, dtype=str)))


DECOYS = ["GAPDH", "ACTB", "CD3E", "Gapdh", "Actb", "Cd3e", "IL32", "TLR4", "Hspa5", "RPS6KA1"]


def _cases():
    u, m = read(U4), read(M4)
    return [("U4 symbols", U4, u, list(u.Symbol), u.Symbol), ("U4 Ensembl ids", U4, u, list(u.ENSG_ID[u.ENSG_ID != ""]), u.ENSG_ID),
            ("M4 symbols", M4, m, list(m.Symbol), m.Symbol), ("M4 Ensembl ids", M4, m, list(m.ENSG_ID[m.ENSG_ID != ""]), m.ENSG_ID)]


def _removed_marker_finder(names, path, exclude=None, include=None, modules="all"):
    # modules="all": every module ON when no category flag is given (the flags themselves override the rule)
    args = Namespace(blocklist=path, gene_col="index", exclude_blocklist_categories=exclude, include_blocklist_categories=include,
                     blocklist_symbol_match="ci_canonical", blocklist_modules=modules)
    out = mf.apply_universal_blocklist(_adata(names), args)
    return set(names) - set(out.var_names)


@pytest.mark.parametrize("remover", [_removed_marker_finder], ids=["marker_finder"])
@pytest.mark.parametrize("case", range(4), ids=["U4-symbols", "U4-ids", "M4-symbols", "M4-ids"])
def test_category_flags_switch_exactly_each_module(remover, case):
    label, path, df, list_names, key = _cases()[case]
    names = list_names + [d for d in DECOYS if d not in set(list_names)]
    genes_of = {mod: set(key[(df.Category == mod) & (key != "")]) for mod in bm.MODULES}
    assert remover(names, path) == set().union(*genes_of.values())
    for mod in bm.MODULES:
        others = set().union(*(g for k, g in genes_of.items() if k != mod))
        assert remover(names, path, exclude=mod) == others, (label, "exclude", mod)
        assert remover(names, path, include=mod) == genes_of[mod], (label, "include", mod)


@pytest.mark.parametrize("labels,on", [(["alpha", "beta"], ("core_stress", "core_ribosomal", "core_mitochondrial", "core_sex", "cell_cycle")),
                                       (["Cycling", "B cells"], ("core_stress", "core_ribosomal", "core_mitochondrial", "core_sex"))])
def test_rule_output_drives_the_flags(labels, on):
    u = read(U4)
    names = list(u.Symbol) + DECOYS
    removed = _removed_marker_finder(names, U4, exclude=bm.excluded_categories(labels))
    assert removed == set(u.Symbol[u.Category.isin(on)])


# ------------------------------------------------------------------------------------------------ 5. wiring
def test_resolve_modules_modes():
    cats = list(bm.MODULES)
    r = bm.resolve_modules(cats, ["T cells", "Cycling T"])
    assert r["exclude"] == ["cell_cycle", "stress_extended"] and r["include"] is None
    assert r["record"]["on"] == ["core_stress", "core_ribosomal", "core_mitochondrial", "core_sex"]
    assert r["record"]["cell_cycle_off_by"] == ["Cycling T"] and r["record"]["override"] is None
    r = bm.resolve_modules(cats, ["Cycling T"], "all")
    assert r["exclude"] is None and r["record"]["off"] == []
    r = bm.resolve_modules(cats, [], "core_stress, core_sex")
    assert r["record"]["on"] == ["core_stress", "core_sex"] and set(r["exclude"]) == set(cats) - {"core_stress", "core_sex"}
    with pytest.raises(ValueError, match="unknown blocklist module"):
        bm.resolve_modules(cats, [], "core_stress,Cell_Cycle")
    # explicit category flags are a user override and win over the mode
    r = bm.resolve_modules(cats, ["Cycling T"], "auto", exclude="core_sex")
    assert r["exclude"] == ["core_sex"] and r["record"]["override"] == {"exclude_categories": ["core_sex"], "include_categories": None}
    # a list without gene groups is left alone, and a non-auto mode is reported as ignored
    u31 = read(os.path.join(ROOT, "tests", "data", "blocklist_human_ungrouped.tsv")).Category.unique()
    r = bm.resolve_modules(u31, ["Cycling T"], "auto")
    assert r["exclude"] is None and r["include"] is None and r["record"]["modular_list"] is False
    assert "ignored" in bm.resolve_modules(u31, [], "all")["record"]["note"]


def _cycling(synthetic_adata, tmp_path, name):
    a = synthetic_adata.copy()
    a.obs["cell_type"] = a.obs["cell_type"].astype(str).replace({"A": "Cycling A"}).astype("category")
    path = tmp_path / name
    a.write_h5ad(path)
    return path


def _find(h5, tmp_path, tag, *extra):
    from scarlet.marker_finder import build_parser, validate_args, run_marker_finder_core
    d = tmp_path / tag
    d.mkdir()
    argv = ["--input", str(h5), "--group_key", "cell_type", "--output_markers", str(d / "m.tsv"),
            "--output_universe", str(d / "u.tsv"), "--n_iterations", "5", "--n_cells", "30", "--ncpu", "1",
            "--output_background_universe", str(d / "bg.tsv"), *extra]
    p = build_parser()
    run_marker_finder_core(validate_args(p, p.parse_args(argv)))
    params = json.load(open(d / "m.params.json"))
    genes = set(pd.read_csv(d / "bg.tsv", sep="\t")["gene_symbol"])
    return params, genes


def test_rule_is_applied_per_dataset_and_recorded(synthetic_h5ad, synthetic_adata, tmp_path):
    """The default list (U4) with the default --blocklist-modules auto: each run decides on its own labels."""
    plain, g_plain = _find(synthetic_h5ad, tmp_path, "plain")
    cyc, g_cyc = _find(_cycling(synthetic_adata, tmp_path, "cyc.h5ad"), tmp_path, "cycling")
    for params in (plain, cyc):
        rec = params["blocklist"]
        assert rec["identity"]["version"] == "U4" and rec["modules"]["modular_list"] is True
        assert rec["modules"]["mode"] == "auto" and rec["modules"]["states"]["stress_extended"] == "OFF"
        # the --imbalanced preset and the blocklist settings are both in params.json
        imb = params["background_strategy"]["imbalanced"]
        assert imb["mode"] == "auto" and imb["profile"] in ("standard", "balanced for unequal class sizes") and imb["threshold"] == pytest.approx(1 / 3)
    assert plain["blocklist"]["modules"]["states"]["cell_cycle"] == "ON"
    assert cyc["blocklist"]["modules"]["states"]["cell_cycle"] == "OFF"
    assert cyc["blocklist"]["modules"]["cell_cycle_off_by"] == ["Cycling A"]
    assert {"MKI67", "TOP2A"}.isdisjoint(g_plain) and {"MKI67", "TOP2A"} <= g_cyc
    for g in ("FOS", "JUN", "MT-CO1", "RPL3", "XIST"):
        assert g not in g_plain and g not in g_cyc


def test_module_override_is_applied_and_recorded(synthetic_h5ad, synthetic_adata, tmp_path):
    cyc = _cycling(synthetic_adata, tmp_path, "cyc.h5ad")
    p_all, g_all = _find(cyc, tmp_path, "all", "--blocklist-modules", "all")
    assert p_all["blocklist"]["modules"]["mode"] == "all" and p_all["blocklist"]["modules"]["off"] == []
    assert {"MKI67", "TOP2A"}.isdisjoint(g_all)
    p_one, g_one = _find(synthetic_h5ad, tmp_path, "one", "--blocklist_modules", "core_mitochondrial")
    assert p_one["blocklist"]["modules"]["on"] == ["core_mitochondrial"]
    assert {"MT-CO1", "MT-ND1"}.isdisjoint(g_one) and {"FOS", "MKI67", "RPL3", "XIST"} <= g_one
