#!/usr/bin/env python3
"""SCARLET command line: one command per function of the Python interface (scarlet.api).

Each command's options are built from its function's signature, so the command line and Python take the same
options with the same defaults (option --some-name = parameter some_name; yes/no options take yes or no).

    scarlet learn    --reference REF.h5ad --labels COLUMN --model MODEL.scarlet
    scarlet assign   --model MODEL.scarlet --data DATA.h5ad --output ANNOTATED.h5ad
    scarlet transfer --reference REF.h5ad --labels COLUMN --query QUERY.h5ad --output FOLDER
    scarlet report   RUN      (RUN: a model file, a transfer folder or an annotated .h5ad)
"""

from __future__ import annotations

import os
import tempfile

os.environ.setdefault("NUMBA_CACHE_DIR", tempfile.gettempdir())
os.environ.setdefault("MPLCONFIGDIR", tempfile.gettempdir())
os.environ.setdefault("HDF5_USE_FILE_LOCKING", "FALSE")

import argparse
import inspect
import sys


def _yes_no(value):
    v = str(value).strip().lower()
    if v in ("yes", "true", "1"):
        return True
    if v in ("no", "false", "0"):
        return False
    raise argparse.ArgumentTypeError(f"expected yes or no (got {value!r})")


def option_name(param: str) -> str:
    return "--" + param.replace("_", "-")


class _HelpFormatter(argparse.ArgumentDefaultsHelpFormatter):
    """The options' help texts are their plain docstring descriptions, never format templates: a percent sign in
    them ("80 % of the draws") is shown as written; the default is added as ArgumentDefaultsHelpFormatter adds it."""

    def _get_help_string(self, action):
        import copy
        plain = copy.copy(action)
        plain.help = (action.help or "").replace("%", "%%")
        return super()._get_help_string(plain)


POSITIONAL = {("report", "run")}           # scarlet report RUN
FLAGS = {"report"}                          # scarlet learn ... --report
NEGATABLE = {"self_assessment"}             # scarlet report RUN --no-self-assessment


def add_command(subparsers, name, fn):
    """One subcommand whose options are the parameters of `fn` (scarlet.api). An option's help text is its
    description in fn's docstring (Args); its type and allowed values are scarlet.api.PARAMS_BY_COMMAND[name]'s,
    then PARAMS'."""
    from scarlet import _docstrings
    from scarlet.api import PARAMS, PARAMS_BY_COMMAND
    sec = _docstrings.sections(fn)
    described = _docstrings.args(fn)
    own = PARAMS_BY_COMMAND.get(name, {})
    summary = sec.get("summary", "").split("\n\n")
    sp = subparsers.add_parser(name, help=summary[0].replace("\n", " ").replace("%", "%%"),
                               description=_docstrings.plain(" ".join(x.replace("\n", " ") for x in summary)),
                               formatter_class=_HelpFormatter)
    for pname, p in inspect.signature(fn).parameters.items():
        meta = own.get(pname, PARAMS[pname])
        help_text = described[pname][1]
        if (name, pname) in POSITIONAL:   # scarlet report RUN
            sp.add_argument(pname, metavar=pname.upper(), help=help_text)
            continue
        typ = meta[0]
        choices = meta[1] if len(meta) > 1 else None
        kwargs = dict(dest=pname, help=help_text)
        if p.default is inspect.Parameter.empty:
            kwargs["required"] = True
        else:
            kwargs["default"] = p.default
        if typ is bool and pname in FLAGS:   # --report (a flag)
            kwargs.update(type=_yes_no, nargs="?", const=True, metavar="{yes,no}")
        elif typ is bool:
            kwargs.update(type=_yes_no, metavar="{yes,no}")
        elif typ is list:
            kwargs.update(nargs=argparse.REMAINDER)
        else:
            kwargs["type"] = typ
            if choices:
                kwargs["choices"] = choices
        sp.add_argument(option_name(pname), **kwargs)
        if pname in NEGATABLE:   # --no-self-assessment
            sp.add_argument("--no-" + pname.replace("_", "-"), dest=pname, action="store_false",
                            help=f"Same as {option_name(pname)} no.")
    return sp


def build_parser() -> argparse.ArgumentParser:
    from scarlet import __version__
    from scarlet.api import COMMANDS
    parser = argparse.ArgumentParser(prog="scarlet", description=__doc__.split("\n\n")[0],
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--version", action="version", version=f"SCARLET {__version__}")
    sub = parser.add_subparsers(dest="_command", required=True, metavar="{" + ",".join(COMMANDS) + "}")
    for name, fn in COMMANDS.items():
        add_command(sub, name, fn)
    return parser


def main(argv=None) -> None:
    from scarlet.api import COMMANDS
    from scarlet.input_checks import NonCountInputError
    args = build_parser().parse_args(argv)
    kwargs = {k: v for k, v in vars(args).items() if not k.startswith("_")}
    try:
        COMMANDS[args._command](**kwargs)
    except (NonCountInputError, ValueError, KeyError, FileNotFoundError) as e:
        sys.exit(f"Error: {e}")


if __name__ == "__main__":
    main()
