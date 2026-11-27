# Getting started

This page installs SCARLET, downloads a small example, runs `learn`, `assign` and `report` at their defaults and
shows the files they write. It takes a few minutes on a laptop.

## Install

Python 3.11 to 3.13, in a fresh environment (a virtual environment or a conda environment):

```bash
pip install "scarlet-sc[report] @ git+https://github.com/matteo14c/scarlet-sc@v0.6.3"
```

This is all you need to use SCARLET. The `report` extra adds the HTML reports (plotly); `pip install scarlet-sc`
alone gives SCARLET without them. Two more extras exist, for other purposes: `tutorials` (to run the notebooks:
[Tutorials](tutorials.md)) and `docs` (only to build this documentation locally: contributors, see
`CONTRIBUTING.md`).

### Check the install

```bash
scarlet --version
python -c "import scarlet; print(scarlet.__version__)"
```

Both print SCARLET's version. `scarlet --help` lists the four commands, and `scarlet learn --help` (or `assign`,
`transfer`, `report`) lists a command's options with their defaults.

## What SCARLET needs

- **AnnData files** (`.h5ad`), cells in rows and genes in columns, with **counts** in `X`: raw counts are the
  safest input. SCARLET reads whether each gene is detected in each cell (a value above 0). Batch-corrected,
  scaled or centred data have negative values and are refused (see [Troubleshooting](troubleshooting.md)).
- **A reference with labels**: one `obs` column holds the cell type of each reference cell.
- **The same kind of gene names** in the reference and the query: gene symbols in both, or Ensembl ids in both.
  By default the gene names are the `var` names; a `var` column can be used instead (`--gene-col`).

## Example data

The tutorials' data are the inputs of SCARLET's benchmark, downloaded with their checksums checked:

```bash
python -m scarlet.datasets --dest data pbmc1_smartseq2.h5ad pbmc1_10xv3_query.h5ad
```

(in Python: `from scarlet import datasets; datasets.fetch(["pbmc1_smartseq2.h5ad", "pbmc1_10xv3_query.h5ad"],
dest="data")`). These are peripheral blood cells of one donor profiled with two technologies: Smart-seq2 (253 cells)
and 10x Chromium v3 (2,990 cells of the same cell types), each with its labels in `obs['cell_type']`.

## A first run

Learn a model from the Smart-seq2 cells:

```bash
scarlet learn --reference data/pbmc1_smartseq2.h5ad --labels cell_type --model pbmc.scarlet
```

Before it starts, SCARLET prints a **pre-flight summary**: the data, the technology preset it chose and why, the
marker sizes, the blocklist it chose from the gene names, and how it handles unequal class sizes. Then one progress
bar, and one line when the model is written (`--verbosity full` prints every step; `quiet` nothing).

Annotate the 10x v3 cells with the model, and also write the calls as a table:

```bash
scarlet assign --model pbmc.scarlet --data data/pbmc1_10xv3_query.h5ad --output 10xv3_annotated.h5ad \
  --predictions 10xv3_calls.tsv
```

`assign` ignores the query's own labels. The query here comes from another technology than the reference: for a
query of another depth or technology, [`transfer`](guide/transfer.md) is the command made for the job (it matches
the depths first and then refines on the query's own cells):

```bash
scarlet transfer --reference data/pbmc1_smartseq2.h5ad --labels cell_type --query data/pbmc1_10xv3_query.h5ad \
  --output pbmc_transfer/
```

Write the reports:

```bash
scarlet report pbmc.scarlet                                     # the model, with its self-assessment
scarlet report 10xv3_annotated.h5ad --model pbmc.scarlet        # the annotated cells
scarlet report pbmc_transfer/ --query-labels cell_type          # the transfer, scored against the query's labels
```

The same in Python:

```python
import scarlet

if __name__ == "__main__":
    model = scarlet.learn("data/pbmc1_smartseq2.h5ad", "cell_type", "pbmc.scarlet")
    calls = scarlet.assign("pbmc.scarlet", "data/pbmc1_10xv3_query.h5ad", "10xv3_annotated.h5ad",
                           predictions="10xv3_calls.tsv")
    res = scarlet.transfer("data/pbmc1_smartseq2.h5ad", "cell_type", "data/pbmc1_10xv3_query.h5ad",
                           "pbmc_transfer/")
    scarlet.report("pbmc_transfer/", query_labels="cell_type")
```

A Python script that calls `transfer` with `ncpu` above 1 must run it under `if __name__ == "__main__":` (its draws
run in new Python processes, which import the script); notebooks need nothing.

## The files

| file | written by | what it is |
| :--- | :--- | :--- |
| `pbmc.scarlet` | learn | the model: one file, used by `assign` ([what it holds](output.md#the-model-file)) |
| `pbmc.summary.txt` | learn | the run's summary in plain words |
| `10xv3_annotated.h5ad` | assign | the query, unchanged, with the calls added to `obs` ([columns](output.md#assign-the-annotated-file)) |
| `10xv3_annotated.summary.txt` | assign | the run's summary |
| `10xv3_calls.tsv` | assign | the calls per cell as a table ([columns](output.md#assign-the-calls-table)) |
| `pbmc_transfer/` | transfer | the calls, the models, the transfer check and the run's record and summary ([files](output.md#transfer)) |
| `pbmc_report/`, `10xv3_annotated_report/`, `pbmc_transfer/report/` | report | `report.html` and a bundle of tables ([reports](reports/index.md)) |

### The summary

Every run writes a few lines in plain words next to its outputs: what the run did, on which data, which technology
preset and blocklist it used and why, whether the cell-cycle genes were blocked, the marker sizes, how unequal class
sizes were handled, the options not at their defaults, the outputs and the versions. For example `pbmc.summary.txt`
after the `learn` above. The same text is stored with the output (inside the model file, in the annotated file's
`uns['scarlet_summary']`, in `transfer_record.json`) and shown by the report under "What this run did". Each line is
explained in [Reading the output](output.md#the-run-summary).

## Reading the calls

The two calls to look at first:

- `scarlet_best`: the most likely cell type, or `unassigned` when the background explains the cell at least as well
  (the cell looks like none of the reference's cell types);
- `scarlet_unique`: the best cell type only when it is clearly better than its competitor, else `unassigned`.

`unique` is the stricter call: it leaves more cells unlabelled, and the cells it labels are more reliable. Every
column is explained in [Reading the output](output.md).

## Next

- One guide per command: [learn](guide/learn.md), [assign](guide/assign.md), [transfer](guide/transfer.md),
  [report](reports/index.md).
- [Tutorials](tutorials.md): three notebooks run the same steps on benchmark datasets and read the results.
