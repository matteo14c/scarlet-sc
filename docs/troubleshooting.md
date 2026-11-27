# Troubleshooting

On the command line an error is printed as `Error: <message>` and SCARLET exits with status 1; in Python the same
message comes with an exception (`ValueError`, `KeyError` or `FileNotFoundError`). The tables below give SCARLET's
own messages (`...` stands for the part that depends on your data), what causes each, and what to do. A row in
italics is a situation rather than a message.

## The input

| Message | Cause | Fix |
| :--- | :--- | :--- |
| `... has ... negative values (min ...): non-count input, e.g. batch-corrected, is outside SCARLET's scope` | The matrix holds negative values: batch-corrected, scaled or centred data. SCARLET reads detection (a value above 0) and needs raw or normalised counts. | Use the counts (often in `adata.layers['counts']` or `adata.raw`): put them in `X` of the file you give SCARLET. |
| `label column '...' is not an obs column of the input` (or `of the reference`) | The `--labels` name is not a column of the reference's `obs`. | Check the spelling: `adata.obs.columns` lists them. |
| `origin column '...' is not an obs column of the query` | The `--origin-col` name is not a column of the query's `obs`. | Check the spelling, as above. |
| `reference: ... does not exist` | The file path is wrong, or relative to another folder than the one you run from. | Give the right path (or an absolute one). |
| *Few marker genes found in the query*: `assign` prints how many of the model's genes it found; `transfer`'s pre-flight advises "no shared genes: check the gene identifiers (symbols vs Ensembl ids) or --gene-col" or "low gene overlap (... shared genes): ..." | Reference and query name their genes differently (for example symbols in one, Ensembl ids in the other). | Use the same kind of identifiers in both, or point `--gene-col` at a `var` column that holds them. |

## The blocklist

What the blocklist is and how to change it: [The blocklist](blocklist.md).

| Message | Cause | Fix |
| :--- | :--- | :--- |
| `blocklist='auto' cannot tell whether the data are human or mouse: ... Pass blocklist='human', 'mouse', 'none' or a file.` | `auto` chooses the list from the gene names: at least 60 % mouse-style names (capitalised symbols such as `Cd3e`, or `ENSMUSG` ids) mean mouse, at most 40 % mean human (upper-case symbols such as `CD3E`, or `ENSG` ids), and at least 20 usable names are needed. Between the two, or with too few names, it stops. | Give the species yourself: `--blocklist human` or `--blocklist mouse` (or `none`, or your own file). |
| `blocklist file '...' not found. Give an existing file, or one of: auto, human, mouse, none ...` | The value is neither one of the four names nor an existing file. SCARLET never runs without the blocklist you asked for. | Give one of the four names or the path of an existing file. |
| `blocklist '...' has a mouse (ENSMUSG) id in its ENSG_ID column for only ... rows ...` | The data's genes are Ensembl ids of one species, but the list has almost no ids of that species, so (almost) no gene could be removed. | Use the list of the data's species (`--blocklist mouse` or `human`), or give `--gene-col`, a `var` column with the genes' symbols. |

## The settings

| Message | Cause | Fix |
| :--- | :--- | :--- |
| `technology must be auto, sparse or dense (got ...)`, `flip must be ...`, `imbalanced must be ...` | An option has a value outside the documented ones. | `scarlet <command> --help` lists each option's values. |
| *The technology preset looks wrong*: `transfer`'s pre-flight advises "both datasets detect >= 5000 genes per cell (plate-based): consider --technology dense", or "the dense preset expects >= 5000 detected genes per cell; with fewer, --technology sparse is the usual choice" | `auto` chose the preset from the detected genes per cell; the pre-flight prints the preset and the number behind the choice. | Give `--technology sparse` or `dense` to override `auto`; the summary and the report then list it among the options not at their defaults. |

## Models

| Message | Cause | Fix |
| :--- | :--- | :--- |
| `this model has no scoring statistics ...: it was not learned by this version of SCARLET; learn it again` | The model was written by another version of SCARLET. | Learn it again with this version. |
| `... is not a SCARLET model file` | The file is not a model. | Give the `.scarlet` file `learn` wrote. |
| `model format version ... is newer than this SCARLET` | The model was written by a newer SCARLET. | Update SCARLET. |
| `member ... does not match its checksum` | The model file is damaged. | Copy it again. |

## The report

| Message | Cause | Fix |
| :--- | :--- | :--- |
| `not a finished SCARLET run: ...` | `report` was given something other than a model file (`.scarlet`), a transfer output folder (with its `transfer_record.json`) or an annotated `.h5ad` written by `assign`. | Give one of those three. |
| `the assign view needs the model the run used: pass model=` | An annotated file was given without its model. | Give `--model` with the model the `assign` run used. |
| `the annotated data were not assigned with this model (identities differ)` | The model given is not the one the `assign` run used. | Give the run's model. |
| `the learn report needs the reference data ...` / `the transfer report needs the ... data ...` | The report reads the reference (and, for a transfer, the query) from the paths the run recorded; they are no longer readable (moved, or relative to another folder). | Give `--reference` (and `--query`). |
| `the assign run recorded alpha ...; the report uses it (got ...)` | An annotated file records the `alpha` its run used, and the assign view reads it from there. | Leave `--alpha` at its default for an assign view. |
| *An error or a warning naming plotly* | The report needs the `report` extra. `scarlet report` stops with this error; `learn` and `transfer` with `--report yes` print a warning, run as usual and skip the report. | `pip install "scarlet-sc[report]"` (or `pip install ".[report]"` from a clone). |
| *A warning that the report failed* | `learn` or `transfer` with `--report yes` finished its run, but the report could not be written. The run's outputs are complete; the warning names the error. | Write the report later with `scarlet report` on the same output. |
| *A reference with one cell type* | `learn`, `assign` and `transfer` run, but with one cell type SCARLET can only tell that cell type from the background, and the refinement of `transfer` does not run (it needs 2 cell types per origin). | Give a reference with at least two cell types. |

## Running

| Message | Cause | Fix |
| :--- | :--- | :--- |
| `transfer with ncpu above 1 runs its draws in new Python processes, which import your script: call transfer under if __name__ == "__main__": (or use ncpu=1, which gives the same result)` | A Python script called `transfer` with `ncpu` above 1 at its top level. The draws run in new Python processes, which import the script again; without the guard each would start the transfer again. | Put the call under the guard (below), or use `ncpu=1` (the same calls, slower). Notebooks and the command line need nothing. |

```python
import scarlet

if __name__ == "__main__":
    scarlet.transfer("reference.h5ad", "cell_type", "query.h5ad", "out", ncpu=4)
```

## Memory and time

- `--chunk-size` (assign and transfer, default 10,000) sets how many cells are scored at a time: lower it to use less
  memory. It does not change the result.
- `--ncpu` (default 1): in `transfer` the draws run side by side, up to `--ncpu` at a time; they share the input
  data and each adds its own working memory: faster, with more memory (measured on two benchmark tasks: 5 CPUs 3.1–3.2 times faster with
  3.3–4.1 times the memory). Choose a number of CPUs that divides the number of draws (with 10 draws: 1, 2, 5 or
  10): values between 5 and 10 bring no further speed-up; memory decides how high to go. A deep dataset of brain
  size (about 11,000 genes detected per cell) needs about 3 GB with one CPU; the second CPU adds at most the 1-CPU
  memory plus 0.5 GB, each further CPU at most half the 1-CPU memory plus 1 GB; the measured time and memory
  per mode and CPU count are in [transfer](guide/transfer.md#run-time-memory-and-cpus), with
  [how to choose `--ncpu`](guide/transfer.md#choosing-ncpu). It does not change the result. In a Python script,
  call `transfer` with `ncpu` above 1 under `if __name__ == "__main__":` (see [Running](#running)).
- **Quick mode.** `transfer --draws 1` runs once, about five times faster than the default 10 draws with one CPU,
  with about the same memory whatever `--ncpu`, seed-dependent calls and no consensus (see
  [transfer](guide/transfer.md#accurate-and-quick-mode)).
- The learn report's self-assessment learns once per fold (5 by default): `--no-self-assessment` skips it for very
  large references.
- A transfer learns several times in each of its 10 draws (the reference, the refinement per origin) and once more
  for the transfer check: it takes many times as long as one `learn`. Measured run times and memory:
  [transfer](guide/transfer.md#run-time-memory-and-cpus).

See also: [Reading the output](output.md), [Glossary](glossary.md).
