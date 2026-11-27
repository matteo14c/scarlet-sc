# Contributing to SCARLET

Thank you for testing SCARLET during the beta. The most useful contributions now are reports of what you tried.

## Reporting

Open an [issue](../../issues) and pick a template:

- **A call looks wrong:** please include the call (`scarlet_best`, `scarlet_unique`, and after `transfer` the
  first-pass calls), its **reason and flags** from the report's bundle (`bundle/tables/cells.tsv` after `assign`,
  `bundle/tables/query_cells.tsv` after `transfer`), the pre-flight summary SCARLET printed, and what you expected.
- **Bug report:** the error, the command or code, the SCARLET version (`scarlet --version`) and your environment.
- **Beta feedback:** what worked, what did not, what was unclear.

Please do not attach data you are not allowed to share. A small example that reproduces the problem (for instance
with one of the notebooks' example files) helps most.

## Code changes

The package holds `learn`, `assign`, `transfer` and `report`. Before writing code, please open an issue to discuss
the change. Pull requests should:

- keep the command line and the Python API identical (the command line is built from the function signatures and
  their docstrings; a test checks it);
- document every option they add or change (purpose, values, default; tests check every option is described);
- not change a default, a call or a score without a stated reason and a test;
- pass the tests: `python -m pytest -q tests` (with `PYTHONHASHSEED=0`, as CI does).

CI runs the tests, the documentation build and the notebooks. The notebooks read the example data published with the
v0.6.3 release; until those files are published, the notebooks job reports them as pending (a notice, not a failure)
and does not run them.

## Building the documentation locally

Users read the documentation online (and, until the public release, as the Markdown files in `docs/` on GitHub).
To preview the site while you edit it, from a clone, in a fresh environment:

```bash
pip install -e ".[docs]"        # the pinned mkdocs packages, the versions CI builds with
mkdocs serve                    # then open http://127.0.0.1:8001
```

`mkdocs serve` keeps its terminal busy and rebuilds the site when a file changes (Ctrl+C stops it). It uses port
8001 (`dev_addr` in `mkdocs.yml`), so it does not clash with Jupyter. CI builds the site with `mkdocs build --strict`:
a broken link or a page missing from the navigation fails the build. The command-line reference
(`docs/reference/cli.md`) is generated from the code: after changing an option, regenerate it with
`SCARLET_WRITE_DOCS=1 python -m pytest tests/test_docs.py -k command_line_reference`.

## Code of conduct

Be respectful and constructive. Report unacceptable behaviour to the maintainers at matteo.chiara@unimi.it.
