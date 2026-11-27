# The assign report

[Reports](index.md) · [learn report](learn.md) · [transfer report](transfer.md)

The report of an annotated `.h5ad` written by `assign`, for a query without labels: its calls and its refusals.

```bash
scarlet report query_annotated.h5ad --model reference.scarlet
scarlet report query_annotated.h5ad --model reference.scarlet --query-labels my_labels   # labels you already know
```

## Sections

| section | content |
| :--- | :--- |
| Summary | cells annotated; the share called and refused, by reason, for the best and the unique call; cell types of the model with no call |
| 1. Depth | genes detected per cell, per called cell type |
| 2. UMAP | the query by its best or unique calls; highlight and space switches ([the UMAPs](index.md#the-umaps)) |
| 3. Analytics | calls per cell type (best solid, unique hatched, the cells left unassigned as a grey bar; the count and share in the hover); refusals by reason |
| 3b. Against labels supplied by the user | only with `query_labels`: the share refused and the median F1 per label on the called cells, per label, and a map of labels against calls |
| 4. Cell type drill-down | a cell type's markers and their detection rate in every cell type of the model and in the background; each marker's detection in the model against its detection in the cells called that cell type |
| 5. Flags | cell types with no call; flip test and low positive evidence per cell type |
| 6. Provenance | the run's summary; the model, the annotated data and its SHA-256, the assign steps, `alpha` (as the run recorded it), package defaults or not, the UMAP spaces |

The refusal reasons are read with the `alpha` the assign run recorded; another value given to the report is an
error.

## Options

`model` (required), `query_labels`, `output`, `seed`, `umap_cells`, `umap_all_genes`, `ncpu`: see the
[learn report](learn.md#options) and [Reports](index.md).

Tables of the bundle: `colours`, `cells` (per cell: calls, refusal reasons, depth, flags), `umap`, `analytics`,
`no_calls`, `marker_detection`, `markers`, `universe`; with `query_labels` also `known_summary`, `known_per_label`
and `known_confusion`.
