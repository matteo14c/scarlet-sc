"""Blocklist identity and matching (synthetic names; the repo blocklist files are read, no data)."""
import os
import re
from argparse import Namespace

import anndata as ad
import numpy as np
import pandas as pd
import pytest

from scarlet import blocklists as bl
from scarlet import marker_finder as mf

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
# the bundled lists (tests/unit/test_blocklist_u4.py) and two lists without gene groups (test data)
DEFAULT_UNIVERSAL = os.path.join(ROOT, "scarlet", "scarlet_blocklist_universal.tsv")
DEFAULT_MOUSE = os.path.join(ROOT, "scarlet", "scarlet_blocklist_mouse.tsv")
UNIVERSAL = os.path.join(ROOT, "tests", "data", "blocklist_human_ungrouped.tsv")
MOUSE = os.path.join(ROOT, "tests", "data", "blocklist_mouse_ungrouped.tsv")


OLD_U3 = os.path.join(ROOT, "tests", "data", "blocklist_U3_original.tsv")
OLD_M2 = os.path.join(ROOT, "tests", "data", "blocklist_M2_original.tsv")


def test_the_package_copies_are_the_bundled_lists():
    for name in ("scarlet_blocklist_universal.tsv", "scarlet_blocklist_mouse.tsv"):
        # a copy at the repository root, where one exists, is the same file; the package holds the one it uses
        root = os.path.join(ROOT, name)
        assert not os.path.exists(root) or bl.file_md5(root) == bl.file_md5(os.path.join(ROOT, "scarlet", name))
    assert bl.blocklist_identity(DEFAULT_UNIVERSAL)["species"] == "human"
    assert bl.blocklist_identity(DEFAULT_MOUSE)["species"] == "mouse"


def test_unknown_and_missing_lists(tmp_path):
    p = tmp_path / "x.tsv"
    p.write_text("ENSG_ID\tSymbol\tCategory\nENSG00000170345\tFOS\tStress\n")
    rec = bl.blocklist_identity(str(p))
    assert rec["version"] is None and rec["species"] is None and rec["n_rows"] == 1
    assert bl.blocklist_identity(str(tmp_path / "nope.tsv"))["error"] == "file not found"
    assert bl.blocklist_identity(None)["md5"] is None


def test_mouse_list_ensembl_column_holds_mouse_ids_only():
    """M3 holds ENSMUSG ids (org.Mm.eg.db 3.18.0) or nothing; M2 held human ids (466) and 26 malformed."""
    df = pd.read_csv(MOUSE, sep="\t")
    ids = df["ENSG_ID"].dropna().astype(str)
    assert ids.str.match(r"^ENSMUSG\d{11}$").all() and len(ids) == 459
    assert dict(zip(df["Symbol"], df["ENSG_ID"]))["Hsp90aa1"] == "ENSMUSG00000021270"
    assert dict(zip(df["Symbol"], df["ENSG_ID"]))["Fos"] == "ENSMUSG00000021250"
    old = pd.read_csv(OLD_M2, sep="\t")["ENSG_ID"].dropna().astype(str)
    assert old.str.match(r"^ENSG\d{11}$").sum() == 466 and old.str.startswith("ENSMUSG").sum() == 26


def test_gene_id_type_and_species():
    assert bl.gene_id_type([f"ENSG{i:011d}" for i in range(20)]) == "ensembl_human"
    assert bl.gene_id_type([f"ENSMUSG{i:011d}" for i in range(20)]) == "ensembl_mouse"
    assert bl.gene_id_type(["FOS", "JUN", "CD3E"] * 10) == "symbol"
    assert bl.gene_id_type([f"ENSG{i:011d}" for i in range(15)] + ["FOS"] * 5) == "mixed"
    human = ["TSPAN6", "DPM1", "SCYL3", "C1orf112", "RP11-34P13.7", "AL627309.1", "MT-CO1", "A1BG__chr19"] * 3
    mouse = ["Xkr4", "Rp1", "Sox17", "0610005C13Rik", "X0610007P14Rik", "Gm12345", "mt-Co1", "Fos"] * 3
    assert bl.symbol_case_species(human)[0] == "human"
    assert bl.symbol_case_species(mouse)[0] == "mouse"
    assert bl.symbol_case_species(["FOS"] * 3)[0] is None             # too few names to decide


def _args(path, **kw):
    a = dict(blocklist=path, gene_col="index", exclude_blocklist_categories=None, include_blocklist_categories=None)
    a.update(kw)
    return Namespace(**a)


def _adata(names):
    return ad.AnnData(X=np.ones((3, len(names)), dtype=np.float32), var=pd.DataFrame(index=names))


@pytest.mark.parametrize("names,lst,removed,by_id,by_sym", [
    (["FOS", "JUN", "MT-CO1", "G1", "G2"], UNIVERSAL, 3, 0, 3),                 # human symbols
    (["Fos", "Jun", "mt-Co1", "G1", "G2"], UNIVERSAL, 3, 0, 3),                 # mouse case, human list
    (["Fos", "Jun", "Rpl10-ps1", "G1"], MOUSE, 3, 0, 3),                        # mouse-only symbol
    (["Fos", "Jun", "Rpl10-ps1", "G1"], UNIVERSAL, 2, 0, 2),
    (["ENSG00000170345", "ENSG00000177606", "ENSG00000000003"], UNIVERSAL, 2, 2, 0),  # human Ensembl
    (["ENSMUSG00000021250", "ENSMUSG00000052684", "ENSMUSG00000000001"], MOUSE, 2, 2, 0),  # mouse Ensembl (M3)
])
def test_marker_finder_matches_the_identifier_the_data_use(names, lst, removed, by_id, by_sym):
    args = _args(lst)
    out = mf.apply_universal_blocklist(_adata(names), args)
    rec = args._blocklist_record
    assert out.n_vars == len(names) - removed
    assert (rec["n_removed"], rec["n_matched_by_id"], rec["n_matched_by_symbol"]) == (removed, by_id, by_sym)
    assert rec["identity"]["md5"] == bl.file_md5(lst) and rec["n_genes_in"] == len(names)
    assert rec["matched_columns"] == [c for c, n in (("ENSG_ID", by_id), ("Symbol", by_sym)) if n]
    assert "warning" not in rec
    if names[0].startswith("ENS"):
        kind = "ensembl_mouse" if names[0].startswith("ENSMUSG") else "ensembl_human"
        assert rec["gene_id_type"] == kind and rec["n_list_rows_without_data_id"] == {
            "ensembl_mouse": 780 - 459, "ensembl_human": 12}[kind]


@pytest.mark.parametrize("kind", ["mouse data, M2 (human ids only)", "human data, a list with mouse ids only"])
def test_ensembl_data_with_a_list_without_ids_of_that_species_is_refused(kind, tmp_path):
    if kind.startswith("mouse"):
        names, lst = ["ENSMUSG00000021250", "ENSMUSG00000052684", "ENSMUSG00000000001"] * 7, OLD_M2
    else:
        names, lst = [f"ENSG{i:011d}" for i in range(20)], str(tmp_path / "mouse_ids_only.tsv")
        pd.DataFrame({"ENSG_ID": ["ENSMUSG00000021250"], "Symbol": ["Fos"], "Category": ["Stress"]}).to_csv(
            lst, sep="\t", index=False)
    with pytest.raises(SystemExit, match=r"has a (mouse|human) .* id in its ENSG_ID column for only [01] of"):
        mf.apply_universal_blocklist(_adata(names), _args(lst))


def test_mask_equals_the_old_inline_rule():
    """blocklist_gene_mask is the inline rule of apply_universal_blocklist before the refactor."""
    df = pd.read_csv(UNIVERSAL, sep="\t")
    rng = np.random.default_rng(0)
    names = list(rng.choice(df["Symbol"].tolist(), 40)) + [s.lower() for s in df["Symbol"][:10]] + \
        list(df["ENSG_ID"][:10]) + [f"G{i}" for i in range(40)] + [" FOS ", "fos"]
    block_ids = set(df["ENSG_ID"].dropna().astype(str).str.strip())
    block_symbols = set(df["Symbol"].dropna().map(mf.normalize_gene_symbol))
    gn = pd.Index([str(g).strip() for g in names])
    old = ~(gn.isin(block_ids) | pd.Index([mf.normalize_gene_symbol(g) for g in gn]).isin(block_symbols))
    mask, counts = mf.blocklist_gene_mask(names, df)
    assert np.array_equal(~mask, np.asarray(old)) and counts["n_removed"] == int((~old).sum())


def test_no_blocklist_and_missing_file_records(tmp_path):
    a = _args(None)
    mf.apply_universal_blocklist(_adata(["FOS"]), a)
    assert a._blocklist_record is None
    a = _args(str(tmp_path / "missing.tsv"))
    with pytest.raises(FileNotFoundError, match="not found"):    # 0.6.3: an error, no longer 'no filtering'
        mf.apply_universal_blocklist(_adata(["FOS"]), a)


def test_params_json_records_the_blocklist(synthetic_h5ad, tmp_path):
    import json
    from scarlet.marker_finder import build_parser, validate_args, run_marker_finder_core
    m = tmp_path / "m.tsv"
    argv = ["--input", str(synthetic_h5ad), "--group_key", "cell_type", "--output_markers", str(m),
            "--output_universe", str(tmp_path / "u.tsv"), "--n_iterations", "5", "--ncpu", "1",
            "--blocklist", UNIVERSAL, "--output_background_universe", str(tmp_path / "bg.tsv")]
    p = build_parser()
    run_marker_finder_core(validate_args(p, p.parse_args(argv)))
    rec = json.load(open(tmp_path / "m.params.json"))["blocklist"]
    assert rec["identity"]["version"] == "U3.1" and rec["identity"]["species"] == "human"
    assert rec["n_removed"] == 9 and rec["matched_columns"] == ["Symbol"] and rec["gene_id_type"] == "symbol"


# ------------------------------------------------------------------ canonical Symbol key
@pytest.mark.parametrize("raw,key", [
    ("ATF3__chr1", "ATF3"), ("A2M-AS1__chr12", "A2M-AS1"), ("MT.CO1", "MT-CO1"), ("MT-CO1", "MT-CO1"),
    ("Rpl10.ps1", "RPL10-PS1"), ("Rpl10-ps1", "RPL10-PS1"), (" fos ", "FOS"), ("RPS10.NUDT3", "RPS10-NUDT3"),
    ("__odd", "__ODD"), ("ENSG00000170345", "ENSG00000170345"),
])
def test_canonical_symbol_key(raw, key):
    assert mf.canonical_symbol_key(raw) == key


def test_canonical_matching_finds_decorated_symbols_only():
    df = pd.read_csv(UNIVERSAL, sep="\t")
    names = ["ATF3__chr1", "FOS__chr14", "MT.CO1", "MT.ND1", "RPS10.NUDT3",            # decorated blocklist genes
             "FOS", "jun",                                                                # plain matches
             "RP11-34P13.7", "AC012345.1", "MTCO1", "FOSX__chr1", "A1BG__chr19", "CD3E.1"]  # must stay
    plain, pc = mf.blocklist_gene_mask(names, df, symbol_match="ci")
    canon, cc = mf.blocklist_gene_mask(names, df)                   # default = ci_canonical
    assert pc["n_removed"] == 2 and pc["match_rule"] == "ensg_exact+symbol_ci" and pc["n_removed_canonical_only"] == 0
    assert cc["n_removed"] == 7 and cc["match_rule"] == "ensg_exact+symbol_ci_canonical"
    assert cc["n_removed_plain_matching"] == 2 and cc["n_removed_canonical_only"] == 5
    assert cc["canonical_only_genes"] == ["ATF3__chr1", "FOS__chr14", "MT.CO1", "MT.ND1", "RPS10.NUDT3"]
    assert not canon[7:].any()
    with pytest.raises(ValueError):
        mf.blocklist_gene_mask(names, df, symbol_match="exact")


def test_canonical_matching_mouse_pseudogenes_and_ensembl_unchanged():
    df = pd.read_csv(MOUSE, sep="\t")
    names = ["Rpl10.ps1", "Rpl10-ps1", "Fos", "Gm12345"]
    assert mf.blocklist_gene_mask(names, df, symbol_match="ci")[1]["n_removed"] == 2
    assert mf.blocklist_gene_mask(names, df)[1]["n_removed"] == 3
    u = pd.read_csv(UNIVERSAL, sep="\t")
    ens = ["ENSG00000170345", "ENSG00000170345.5", "ensg00000170345"]   # ENSG_ID stays exact
    assert mf.blocklist_gene_mask(ens, u)[1]["n_matched_by_id"] == 1
    assert mf.blocklist_gene_mask(ens, u)[1]["n_removed"] == 1


@pytest.mark.parametrize("mode,removed", [(None, 3), ("ci_canonical", 3), ("ci", 1)])
def test_apply_blocklist_symbol_match_mode(mode, removed):
    args = _args(UNIVERSAL, **({} if mode is None else {"blocklist_symbol_match": mode}))
    out = mf.apply_universal_blocklist(_adata(["FOS", "ATF3__chr1", "MT.CO1", "G1"]), args)
    assert out.n_vars == 4 - removed and args._blocklist_record["symbol_match"] == (mode or "ci_canonical")


def test_marker_finder_symbol_match_is_fixed():
    """The canonical match is the method's (not an option of the marker finder's parser)."""
    from scarlet.marker_finder import build_parser
    assert build_parser().parse_args(["--input", "x"]).blocklist_symbol_match == "ci_canonical"
    with pytest.raises(SystemExit):
        build_parser().parse_args(["--input", "x", "--blocklist_symbol_match", "ci"])
