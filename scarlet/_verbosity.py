"""How much learn, assign, transfer and report print on the standard output (their ``verbosity`` option).

- ``quiet``: nothing (errors still stop the command; their messages are unchanged).
- ``normal`` (the default): the command's pre-flight once, one progress bar, one line per fold or draw, and the
  command's closing lines.
- ``full``: everything every step prints (each marker search, each learn and assign inside a transfer or a report,
  their own progress bars), as before this option existed.

The level is the outermost command's: a command called inside another (the learns and assigns of a transfer's draws
or of a report's folds) keeps it. Draw workers started as new Python processes read it from the environment variable
SCARLET_VERBOSITY, which the outermost command sets while it runs. Code that calls the steps outside a command
(``python -m scarlet.marker_finder``) prints everything, as before. The level never changes a result or a record.
"""
from __future__ import annotations

import contextlib
import io
import os
import sys

LEVELS = ("quiet", "normal", "full")
ENV = "SCARLET_VERBOSITY"
_state = {"level": None, "top": None, "out": None, "nested": 0}


class _Discard(io.TextIOBase):
    """The standard output of the steps while a command runs below ``full``."""

    def write(self, s):
        return len(s)

    def writable(self):
        return True


def check(verbosity):
    if verbosity not in LEVELS:
        raise ValueError(f"verbosity must be one of {', '.join(LEVELS)} (got {verbosity!r})")


def level():
    return _state["level"] or os.environ.get(ENV) or "full"


def top():
    """The outermost command running in this process (None outside a command and in a draw worker)."""
    return _state["top"]


def shows(at):
    return LEVELS.index(level()) >= LEVELS.index(at)


def at_top(command):
    """The level at which `command`'s own lines (pre-flight, closing lines) show: normal when it is the outermost
    command (not one of the same name inside it, e.g. a report's fold learns inside learn(report=True)), otherwise
    full only."""
    return "normal" if _state["top"] == command and not _state["nested"] else "full"


@contextlib.contextmanager
def command(name, verbosity):
    """Entered by learn, assign, transfer and report. The outermost sets the level for everything it runs; below
    full, what the steps print is discarded and only say()'s lines reach the caller's standard output."""
    if _state["top"] is not None:   # inside another command: its level holds
        _state["nested"] += 1
        try:
            yield False
        finally:
            _state["nested"] -= 1
        return
    check(verbosity)
    env = os.environ.get(ENV)
    _state.update(level=verbosity, top=name, out=sys.stdout)
    os.environ[ENV] = verbosity
    try:
        if verbosity == "full":
            yield True
        else:
            with contextlib.redirect_stdout(_Discard()):
                yield True
    finally:
        _state.update(level=None, top=None, out=None)
        if env is None:
            os.environ.pop(ENV, None)
        else:
            os.environ[ENV] = env


def _out():
    return _state["out"] or sys.stdout


def say(msg, at="normal"):
    """msg on the caller's standard output when the level shows `at` (above an open progress bar)."""
    if not shows(at):
        return
    from tqdm import tqdm
    tqdm.write(msg, file=_out())
    _out().flush()


def bar(total, desc):
    """The command's one progress bar (normal only; at full the steps' own bars show instead)."""
    from tqdm import tqdm
    return tqdm(total=total, desc=desc, file=_out(), disable=level() != "normal", leave=True)


def step_bars():
    """Whether the steps' own progress bars show: at full; at normal only learn's own marker search (its one bar; not
    the learns of a report it writes)."""
    return level() == "full" or (level() == "normal" and _state["top"] == "learn" and not _state["nested"])


def enter_worker():
    """In a draw worker (a new Python process): below full, its standard output is discarded."""
    if level() != "full":
        sys.stdout = _Discard()
