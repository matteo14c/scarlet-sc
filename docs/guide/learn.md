# learn

`learn` builds a model from a labelled reference. For each cell type of the reference it finds the marker genes
that the type's cells consistently detect and that are rarely detected elsewhere, and stores, for every marker, how
often each cell type and the background detect it. The model is one file; [`assign`](assign.md) uses it to annotate
new cells.

```bash
scarlet learn --reference reference.h5ad --labels cell_type --model reference.scarlet
```
```python
import scarlet
model = scarlet.learn("reference.h5ad", "cell_type", "reference.scarlet")
```

## Inputs

- `reference`: an `.h5ad` file (or, in Python, an AnnData object) of **counts**, cells in rows, genes in columns.
  SCARLET reads detection: a gene counts as detected in a cell when its value is above 0. Raw counts are the safest
  input; negative values (batch-corrected, scaled or centred data) stop the run.
- `labels`: the `obs` column with each cell's cell type. Any labels can be used.
- `model`: where to write the model file (`.scarlet` by convention). In Python, `None` returns the model without
  writing it.

## How learn works

In each of 100 rounds, learn draws 100 cells of every cell type (with replacement; a smaller cell type contributes as
many cells as it has, drawn with replacement) and ranks the genes for each cell type by how much more often its
cells detect them than the cells as a whole. A gene earns a point of **stability** in a round when it is on the cell
type's shortlist of the 150 best-ranked genes. A cell type's **markers** are its 50 most stable genes. The model
stores, for every marker, the detection rate of every cell type and of the background. These numbers are
[advanced options](../advanced.md); the defaults are what the benchmarks used.

Before the marker search, genes on the [blocklist](#blocklist) are removed. The ranking penalises genes that many
cells of the whole reference detect, more strongly with the sparse [technology preset](#technology-preset) than with
the dense one.

## Technology preset

`--technology` tells SCARLET how dense the data are. Both presets read the same thing from every cell: a gene is
**detected** when its count is above 0. They differ in one setting of the marker search, the **background
penalty**:

| preset | typical data | background penalty |
| :--- | :--- | :--- |
| `sparse` | droplet data (e.g. 10x), a few hundred to a few thousand genes detected per cell | 6 |
| `dense` | plate data (e.g. Smart-seq2, CEL-Seq2), 5,000 or more genes detected per cell | 0 |

**The background penalty in plain words.** A good marker of a cell type is detected in most of its cells and in few
other cells. When SCARLET scores a candidate marker, it multiplies the score by (1 − p)<sup>penalty</sup>, where p is
the share of all the reference's cells that detect the gene. With the sparse preset's 6, a gene detected in 30 % of
all cells keeps about an eighth of its score, so markers detected almost only in their own cell type win. On dense
data most genes are detected in a large share of all cells, so the same penalty would leave few usable markers; the
dense preset's 0 ranks markers on how much more often their own cell type detects them. Neither value is best on
every dataset of the published benchmarks, which is why both presets exist.

With `auto` (the default) learn chooses `dense` when the reference detects at least 5,000 genes per cell on average,
and `sparse` otherwise. The pre-flight prints the choice and the number it was based on, and the model records it.
`--technology sparse` or `--technology dense` is used as given. Both presets use 50 markers per cell type ranked on
a shortlist of 150.

## Blocklist

What is on the lists, how `auto` chooses and how to change them: [The blocklist](../blocklist.md).

Before the marker search, SCARLET keeps known technical and broadly expressed genes (stress, ribosomal,
mitochondrial and sex-linked genes, and, unless a label names dividing cells, cell-cycle genes) out of the markers.
`--blocklist` takes:

| value | what it does |
| :--- | :--- |
| `auto` (default) | chooses the bundled human or mouse list from the gene names: mostly upper-case symbols (`CD3E`) or `ENSG` ids mean human, mostly capitalised symbols (`Cd3e`) or `ENSMUSG` ids mean mouse; it stops with an error when it cannot tell |
| `human`, `mouse` | the bundled list for that species |
| `none` | no gene is kept out |
| a file path | your own tab-separated list with a `Symbol` and/or `ENSG_ID` column; a path that does not exist is an error |

The pre-flight, the summary and the report say which list was used and why. Genes are matched by Ensembl id or by
symbol (case-insensitive). The bundled lists are divided into gene groups; which groups apply is decided from the
labels ([advanced option `blocklist_modules`](../advanced.md#blocklist_modules)). The summary and the report state
whether the cell-cycle genes were blocked and, if not, which labels kept them usable.

## Unequal class sizes

`--imbalanced` (`auto`, `yes`, `no`) handles references where one cell type holds most of the cells. With `auto`
(the default) SCARLET corrects when the largest cell type holds more than a third of the cells. The pre-flight, the
summary and the report say whether it corrected. How the correction works:
[Advanced options](../advanced.md#unequal-class-sizes-imbalanced).

## Other options

| option | default | what it does |
| :--- | :--- | :--- |
| `--gene-col` | `index` | the `var` column with the gene names; `index` uses the `var` names |
| `--seed` | `42` | random seed; the same seed and inputs give the same model |
| `--ncpu` | `1` | number of CPU cores; it does not change the result |
| `--report` | off | also write the [learn report](../reports/learn.md) when the model is written |
| `--verbosity` | `normal` | what the command prints: `normal` (the pre-flight, one progress bar, a closing line), `quiet` (nothing), `full` (every step); it does not change the result |

The [advanced options](../advanced.md) `--n-markers`, `--marker-shortlist`, `--n-iterations`, `--n-cells` and
`--blocklist-modules` are also available; the published benchmarks hold at their defaults only.

## What learn writes

| file | what it is |
| :--- | :--- |
| `reference.scarlet` (the `--model` path) | the model ([members](../output.md#the-model-file)) |
| `reference.summary.txt` | the run's [summary](../output.md#the-run-summary), next to the model; also stored inside the model |
| `reference_report/` | with `--report`: the learn report |

In Python, `learn` returns the model (a `Model`): `model.identities` lists the cell types, `model.markers` the
markers, `model.settings` what the run resolved. A model file is read back with `scarlet.load_model(path)`.

## What can go wrong

- The label column is not in `obs`, the reference does not exist, or the data hold negative values: see
  [Troubleshooting](../troubleshooting.md).
- `blocklist auto` cannot tell human from mouse: pass `--blocklist human`, `mouse`, `none` or a file.
- A cell type with few cells (tens or fewer) gets markers from very few cells; the
  [learn report](../reports/learn.md)'s weak-class table shows how well each cell type is separated.

See also: [Reading the output](../output.md), [Advanced options](../advanced.md),
[Troubleshooting](../troubleshooting.md).
