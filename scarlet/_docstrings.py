"""The sections of the public functions' docstrings (Google style), read by the command line for its help texts.

The docstrings of scarlet.api.learn, assign, transfer and report are the one description of every option: the
command line shows the same text (scarlet/cli.py), the documentation's reference renders it, and
tests/test_docs.py checks that every option is described.
"""
from __future__ import annotations

import inspect
import re

_ENTRY = re.compile(r"^(\w+) \(([^)]*)\):\s*(.*)$")
_SECTION = re.compile(r"^([A-Z][A-Za-z ]*):\s*$")


def sections(fn) -> dict:
    """{section title: its text, dedented} of fn's docstring; the text before the first section is 'summary'."""
    doc = inspect.getdoc(fn) or ""
    out, title, buf = {}, "summary", []
    for line in doc.splitlines():
        m = _SECTION.match(line)
        if m and not line.startswith(" "):
            out[title] = "\n".join(buf).strip("\n")
            title, buf = m.group(1), []
            continue
        buf.append(line[4:] if title != "summary" and line.startswith("    ") else line)
    out[title] = "\n".join(buf).strip("\n")
    return out


def args(fn) -> dict:
    """{parameter: (type text, description in one line)} from the Args section of fn's docstring."""
    text = sections(fn).get("Args", "")
    out, name = {}, None
    for line in text.splitlines():
        m = _ENTRY.match(line) if not line.startswith(" ") else None
        if m:
            name = m.group(1)
            out[name] = [m.group(2), m.group(3).strip()]
        elif name is not None and line.strip():
            out[name][1] += " " + line.strip()
    return {k: (t, plain(d)) for k, (t, d) in out.items()}


def plain(text: str) -> str:
    """Docstring markup removed, for a terminal: ``x`` -> x, [text](link) -> text."""
    text = re.sub(r"``([^`]*)``", r"\1", text)
    text = re.sub(r"\[([^\]]*)\]\([^)]*\)", r"\1", text)
    return text
