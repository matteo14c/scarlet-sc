"""The gene groups of the bundled blocklists and the fixed rule that decides which apply.

The groups (modules) are the ``Category`` values of the bundled human and mouse lists (``scarlet_blocklist_universal
.tsv``, ``scarlet_blocklist_mouse.tsv``). The marker search applies the rule through ``resolve_modules``
(``blocklist_modules='auto'``, the default); ``blocklist_modules='all'`` or a list of groups overrides it. The groups
applied are recorded in the model's settings and provenance.

The rule (decided from the dataset and its reference labels, never from a result):
  * core_stress, core_ribosomal, core_mitochondrial, core_sex: always ON;
  * stress_extended: always OFF;
  * cell_cycle: OFF when at least one reference label names a proliferative state or a proliferative developmental
    stage (``CC_OFF_TERMS``, case-insensitive regular expressions on the stripped label name); otherwise ON.
    A label matching one of ``CC_OFF_EXCLUDE`` (quiescent lineage names) is not matched by the generic stage terms
    (progenitor, precursor, stem cell, ...), but an explicit cycling term (``CC_OFF_STATE_TERMS``: cycling,
    proliferating, dividing, mitotic) still switches cell_cycle off: "cycling oligodendrocyte precursor
    cell" and "cycling oligodendrocyte_progenitor" turn it OFF. Oligodendrocyte progenitors are excluded under each
    spelling: ``oligodendrocyte_progenitor`` and the Cell Ontology name "oligodendrocyte
    precursor cell" (which the bare term ``precursor`` matched); "oligodendrocyte progenitor cell" matches no
    term.
The rule reads the labels being learned. Known residual collateral: a tissue whose labels hide cycling sub-clusters
(no label names cycling) keeps cell_cycle ON.
"""
from __future__ import annotations

import re
from typing import Iterable

MODULES = ("core_stress", "core_ribosomal", "core_mitochondrial", "core_sex", "cell_cycle", "stress_extended")
CORE_MODULES = ("core_stress", "core_ribosomal", "core_mitochondrial", "core_sex")
ALWAYS_OFF = ("stress_extended",)

CC_OFF_TERMS = (
    r"cycl", r"prolif", r"dividing", r"mitotic", r"transit[ _-]?amplif", r"\bTA\b", r"germinal", r"GC_B", r"plasmablast",
    r"\b(pro|pre)-?B\b", r"pro-T", r"erythroblast", r"myeloblast", r"monoblast", r"promonocyte", r"myelocyte", r"promyelocyte",
    r"granulocytopoietic", r"spermatogon", r"spermatocyt", r"stem[ _]cell", r"^progenitor$", r"hematopoietic[ _](stem|precursor|progenitor)",
    r"(neutrophil|monocyte|erythrocyte|erythroid|megakaryocyte|myeloid|multipotent|fetal)[ _-]*(erythroid[ _-]*)?progenitor", r"precursor",
)
# explicit cycling terms (a subset of CC_OFF_TERMS): an exclude does not neutralise them
CC_OFF_STATE_TERMS = (r"cycl", r"prolif", r"dividing", r"mitotic")
CC_OFF_EXCLUDE = (
    r"fibro-adipogenic", r"oligodendrocyte_progenitor", r"mesenchymal stem cell", r"pituitary_stem", r"fallopian_basal_stem",
    r"hepatic_stellate", r"oligodendrocyte[ _]precursor",
)
_TERMS = tuple(re.compile(t, re.I) for t in CC_OFF_TERMS)
_STATE_TERMS = tuple(re.compile(t, re.I) for t in CC_OFF_STATE_TERMS)
_EXCLUDE = tuple(re.compile(t, re.I) for t in CC_OFF_EXCLUDE)


def proliferative_labels(label_names: Iterable[str]) -> list[str]:
    """The reference labels that switch cell_cycle off, in input order (stripped names)."""
    hits = []
    for lab in label_names:
        name = str(lab).strip()
        if any(x.search(name) for x in _STATE_TERMS):
            hits.append(name)
            continue
        if any(x.search(name) for x in _EXCLUDE):
            continue
        if any(x.search(name) for x in _TERMS):
            hits.append(name)
    return hits


def module_set(label_names: Iterable[str]) -> dict:
    """Module states for a dataset with these reference label names.

    Returns {"on": [...], "off": [...], "states": {module: "ON"|"OFF"}, "cell_cycle_off_by": [labels]}.
    """
    hits = proliferative_labels(label_names)
    states = {m: "ON" for m in CORE_MODULES}
    states["cell_cycle"] = "OFF" if hits else "ON"
    for m in ALWAYS_OFF:
        states[m] = "OFF"
    states = {m: states[m] for m in MODULES}
    return {"on": [m for m in MODULES if states[m] == "ON"], "off": [m for m in MODULES if states[m] == "OFF"],
            "states": states, "cell_cycle_off_by": hits}


def excluded_categories(label_names: Iterable[str]) -> str:
    """The value for ``--exclude-blocklist-categories`` (comma-separated module names) for these labels."""
    return ",".join(module_set(label_names)["off"])


MODES_HELP = ("auto (default): the rule above on the labels being learned; all: every module ON (stress_extended "
              "included); or a comma-separated list of the modules to keep ON")


def resolve_modules(categories: Iterable[str], label_names: Iterable[str], mode: str | None = "auto",
                    exclude=None, include=None) -> dict:
    """Which categories of a blocklist to apply, and the record of why.

    categories: the Category values of the list. The rule acts only on a list with gene groups (every category is one
    of MODULES, as in the bundled lists); on any other list nothing is decided here.
    mode: "auto" (default) = module_set(label_names); "all" = every module ON; otherwise a comma-separated list of the
    modules to keep ON. Explicit exclude / include categories (--exclude/--include-blocklist-categories) are a user
    override and take precedence over the mode.
    Returns {"exclude": [...] | None, "include": [...] | None, "record": {...}} for apply_universal_blocklist.
    """
    cats = sorted({str(c) for c in categories if str(c)})
    modular = bool(cats) and set(cats) <= set(MODULES)
    mode = (mode or "auto").strip()
    split = lambda v: [c.strip() for c in v.split(",") if c.strip()] if isinstance(v, str) else (list(v) if v else [])  # noqa: E731
    rec = {"modular_list": modular, "mode": mode, "rule": "scarlet.blocklist_modules.module_set (the fixed rule on the labels)",
           "override": None, "states": None, "on": None, "off": None, "cell_cycle_off_by": None}
    if exclude or include:
        rec["override"] = {"exclude_categories": split(exclude) or None, "include_categories": split(include) or None}
        return {"exclude": split(exclude) or None, "include": split(include) or None, "record": rec}
    if not modular:
        if mode != "auto":
            rec["note"] = "the list has no gene groups (its categories are not the bundled lists' groups): blocklist_modules ignored"
        return {"exclude": None, "include": None, "record": rec}
    if mode == "auto":
        ms = module_set(label_names)
        rec.update(states=ms["states"], on=ms["on"], off=ms["off"], cell_cycle_off_by=ms["cell_cycle_off_by"])
        return {"exclude": ms["off"] or None, "include": None, "record": rec}
    on = list(MODULES) if mode == "all" else split(mode)
    unknown = [m for m in on if m not in MODULES]
    if unknown:
        raise ValueError(f"unknown blocklist module(s) {unknown}; modules are {list(MODULES)}")
    states = {m: ("ON" if m in on else "OFF") for m in MODULES}
    rec.update(states=states, on=[m for m in MODULES if states[m] == "ON"], off=[m for m in MODULES if states[m] == "OFF"])
    return {"exclude": rec["off"] or None, "include": None, "record": rec}
