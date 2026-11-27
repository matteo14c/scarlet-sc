# assign

`assign` annotates new cells with a model written by [`learn`](learn.md). It scores each cell against every cell
type of the model and against the background, gives the most likely cell type, tests whether it is clearly better
than its competitor, and flags calls that deserve a second look.

```bash
scarlet assign --model reference.scarlet --data query.h5ad --output query_annotated.h5ad
```
```python
import scarlet
calls = scarlet.assign("reference.scarlet", "query.h5ad", "query_annotated.h5ad")   # one row per cell
```

## Inputs

- `model`: a model file (or, in Python, a `Model` from `learn` or `scarlet.load_model`).
- `data`: the cells to annotate, an `.h5ad` file (or an AnnData) of counts, with the same kind of gene names as the
  reference. Labels in its `obs`, if any, are ignored. Marker genes the data do not contain are not scored.
- `output`: the annotated `.h5ad` to write: the input as given, with the calls added. In Python, `None` writes
  nothing and only returns the calls.

With a model learned with the **dense** technology preset, assign first reduces the query cells to a comparable
number of detected genes exactly as learn did for the reference (the model's targets and its learn seed). The
annotated file keeps your input's values. A sparse model reads the data as given.

`assign` does not match depths: it scores the query at its own depth. For a query of another depth or technology
than the reference, use [`transfer`](transfer.md).

## How a cell is scored

For each cell, SCARLET computes how likely its pattern of detected and undetected marker genes is under every cell
type's detection rates and under the background's. Then:

1. the most likely cell type is the **forced** call (`scarlet_top1`);
2. if the background explains the cell at least as well as that cell type, the cell is **refused**:
   `scarlet_best` is `unassigned`; otherwise `scarlet_best` is the forced call;
3. the best cell type is tested against its competitor, the stronger of the second cell type and the background.
   The p-values are adjusted over the cells of the run (Benjamini-Hochberg) into `scarlet_fdr`. When `scarlet_fdr`
   is below 0.01 (the [advanced option `alpha`](../advanced.md#alpha)) the call is **unique**: `scarlet_unique` holds
   the cell type; otherwise it is `unassigned`.

A detection rate measured on few cells cannot be trusted close to 0 or 1, so each cell type's rates are kept away
from 0 and 1 by a margin that shrinks as the cell type has more reference cells (the floor). A single marker that a
small cell type never showed therefore cannot rule it out on its own.

## Flags

Flags label a call; they do not change it, except the flip test when it refuses.

- **Flip test** (`scarlet_flip`): every cell is also scored without the floor. When the most likely cell type
  differs between the two scorings, the call depends on genes the candidate cell type never showed in the
  reference. In plain words: a cell type that never showed a gene in the reference may still have it, SCARLET just
  saw too few of its cells to know; the floor gives such a gene a small, honest probability. If the cell's top call
  changes when such genes are treated as almost impossible instead, the call rests on genes the cell type never
  showed. `--flip refuse` (the default, for every model) **refuses** such a cell (`scarlet_reason`: "the call
  depends on genes the candidate class never showed"); on the published benchmarks' annotation within a dataset,
  refusing them never lowered the median F1 of any dataset. `--flip flag` keeps the call and only marks it in
  `scarlet_flip`. ([transfer](transfer.md) flags by default: across datasets the accuracy of flipped calls varies
  widely, from 14 % to 90 % right depending on the task.)
- **Positive evidence** (`scarlet_pe_low`): for the forced call, SCARLET counts how many of the cell type's markers
  the cell detects. The call is flagged when that is fewer than a fifth (rounded down) of the number a cell of that
  type is expected to detect (the sum of the type's detection rates over its markers). `scarlet_pe_note` gives the
  numbers and the typical value among the reference's cells of that type. The call is kept.

## Options

| option | default | what it does |
| :--- | :--- | :--- |
| `--predictions` | none | also write the calls per cell as a tab-separated table to this file |
| `--gene-col` | `index` | the `var` column with the gene names; `index` uses the `var` names |
| `--flip` | `refuse` | the flip test: `refuse` (a flipped call is left unassigned) or `flag` (kept, marked) |
| `--save-probs` | `yes` | store each cell's probability for every cell type in `obsm['scarlet_probs']` |
| `--chunk-size` | `10000` | cells scored at a time; lower it to use less memory; it does not change the result |
| `--ncpu` | `1` | number of CPU cores; it does not change the result |
| `--verbosity` | `normal` | what the command prints: `normal` (the pre-flight and the result), `quiet` (nothing), `full` (every step); it does not change the result |

The [advanced option](../advanced.md#alpha) `--alpha` (default 0.01) is also available; the published benchmarks
hold at its default only.

## What assign writes

| file | what it is |
| :--- | :--- |
| `query_annotated.h5ad` (the `--output` path) | the input with the calls in `obs`, the probabilities in `obsm`, and the run's records in `uns` ([columns](../output.md#assign-the-annotated-file)) |
| `query_annotated.summary.txt` | the run's [summary](../output.md#the-run-summary), next to the output; also in `uns['scarlet_summary']` |
| the `--predictions` file | with `--predictions`: the calls per cell ([columns](../output.md#assign-the-calls-table)) |

In Python, `assign` returns the calls table (one row per cell).

## Reading the result

`unique` is the stricter call: it leaves more cells unlabelled, and the cells it labels are more reliable. `best`
labels more cells. For every refused cell, the [assign view of the report](../reports/assign.md) gives the
reason: background, ambiguity or flip test. Every column: [Reading the output](../output.md#assign).

## What can go wrong

- Negative values in the data, a model file that does not load, or a model without scoring statistics (written by
  another version): see [Troubleshooting](../troubleshooting.md).
- Few of the model's markers in the data (different gene identifiers, symbols against Ensembl ids): fewer genes
  support each call; check that reference and query use the same kind of gene names.

See also: [Reading the output](../output.md), [Advanced options](../advanced.md),
[Troubleshooting](../troubleshooting.md).
