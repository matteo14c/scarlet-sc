"""The documentation standard: every public function, option and model member is documented, the command line shows
the same descriptions, and the command-line reference in docs/ is the one generated from the code.

Regenerate docs/reference/cli.md after changing an option or its description:

    SCARLET_WRITE_DOCS=1 python -m pytest -q tests/test_docs.py
"""
from __future__ import annotations

import inspect
import os
from pathlib import Path

import pytest

import scarlet.cli as cli
from scarlet import _docstrings, api
from scarlet.model import Model, load_model

ROOT = Path(__file__).resolve().parents[1]
CLI_MD = ROOT / "docs" / "reference" / "cli.md"
COMMANDS = list(api.COMMANDS.items())
SECTIONS = {"learn": ("Args", "Returns", "Writes", "Raises", "Examples"),
            "assign": ("Args", "Returns", "Writes", "Raises", "Examples"),
            "transfer": ("Args", "Returns", "Writes", "Raises", "Examples"),
            "report": ("Args", "Returns", "Writes", "Raises", "Examples"),
            "load_model": ("Args", "Returns", "Raises", "Examples")}
FUNCTIONS = COMMANDS + [("load_model", load_model)]


@pytest.mark.parametrize("name, fn", FUNCTIONS)
def test_every_public_function_has_the_required_sections(name, fn):
    sec = _docstrings.sections(fn)
    assert len(sec["summary"]) > 40, name
    for title in SECTIONS[name]:
        assert sec.get(title, "").strip(), (name, title)


@pytest.mark.parametrize("name, fn", FUNCTIONS)
def test_every_option_is_described_in_signature_order(name, fn):
    described = _docstrings.args(fn)
    assert list(described) == list(inspect.signature(fn).parameters), name
    for p, (typ, text) in described.items():
        assert typ and len(text) > 20, (name, p)


@pytest.mark.parametrize("name, fn", COMMANDS)
def test_the_command_line_shows_the_docstring(name, fn):
    sub = cli.build_parser()._subparsers._group_actions[0].choices[name]
    described = _docstrings.args(fn)
    for a in sub._actions:
        if a.dest in described and not any(o.startswith("--no-") for o in a.option_strings):
            assert a.help == described[a.dest][1], (name, a.dest)


def test_model_and_transfer_result_members_are_documented():
    for cls in (Model, api.TransferResult):
        assert inspect.getdoc(cls) and len(inspect.getdoc(cls)) > 40
        for k, v in vars(cls).items():
            if k.startswith("_") or not (callable(v) or isinstance(v, property)):
                continue
            assert inspect.getdoc(v), f"{cls.__name__}.{k}"


def coverage():
    """The documentation's coverage numbers (reported with the documentation review)."""
    fns = [(n, f) for n, f in FUNCTIONS]
    params = [(n, p) for n, f in fns for p in inspect.signature(f).parameters]
    described = sum(p in _docstrings.args(f) for n, f in fns for p in inspect.signature(f).parameters)
    sections = sum(bool(_docstrings.sections(f).get(t, "").strip()) for n, f in fns for t in SECTIONS[n])
    members = [(c, k) for c in (Model, api.TransferResult) for k, v in vars(c).items()
               if not k.startswith("_") and (callable(v) or isinstance(v, property))]
    documented = sum(bool(inspect.getdoc(getattr(c, k))) for c, k in members)
    opts = [a for n, _ in COMMANDS for a in cli.build_parser()._subparsers._group_actions[0].choices[n]._actions
            if a.dest != "help" and not any(o.startswith("--no-") for o in a.option_strings)]
    return {"functions": (len(fns), len(fns)), "parameters": (described, len(params)),
            "sections": (sections, sum(len(SECTIONS[n]) for n, _ in fns)),
            "model members": (documented, len(members)),
            "command-line options with help": (sum(bool(a.help) for a in opts), len(opts))}


def cli_reference() -> str:
    """docs/reference/cli.md, generated from the commands' signatures, option types and docstrings."""
    out = ["# Command-line reference", "",
           "Generated from the code (tests/test_docs.py); every option is the Python parameter of the same name with",
           "the same default ([Python reference](python.md)). yes/no options take `yes` or `no`.", ""]
    for name, fn in COMMANDS:
        sec = _docstrings.sections(fn)
        described = _docstrings.args(fn)
        own = api.PARAMS_BY_COMMAND.get(name, {})
        out += [f"## scarlet {name}", "", _docstrings.plain(sec["summary"].split("\n\n")[0].replace("\n", " ")), "",
                "| option | values | default | description |", "| :--- | :--- | :--- | :--- |"]
        for p, prm in inspect.signature(fn).parameters.items():
            meta = own.get(p, api.PARAMS[p])
            opt = p.upper() if (name, p) in cli.POSITIONAL else cli.option_name(p)
            if meta[0] is bool:
                vals = "yes, no"
            elif len(meta) > 1:
                vals = ", ".join(meta[1])
            else:
                vals = {str: "text", int: "integer", float: "number"}[meta[0]]
            if prm.default is inspect.Parameter.empty:
                default = "required"
            elif prm.default is None:
                default = "none"
            elif isinstance(prm.default, bool):
                default = "yes" if prm.default else "no"
            else:
                default = str(prm.default)
            out.append(f"| `{opt}` | {vals} | {default} | {described[p][1].replace('|', '/')} |")
        out.append("")
    return "\n".join(out)


def test_the_command_line_reference_is_current():
    text = cli_reference()
    if os.environ.get("SCARLET_WRITE_DOCS") == "1":
        CLI_MD.parent.mkdir(parents=True, exist_ok=True)
        CLI_MD.write_text(text)
    assert CLI_MD.read_text() == text, "docs/reference/cli.md is out of date: regenerate it (see this file's docstring)"
