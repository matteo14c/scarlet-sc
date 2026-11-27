# The learn report

[Reports](index.md) · [assign report](assign.md) · [transfer report](transfer.md)

The report of a model file: the reference it was learned on, its markers, and how well the model calls its own
reference.

```bash
scarlet report reference.scarlet                      # the reference is read from the path the model records
scarlet report reference.scarlet --reference ref.h5ad # when that path is no longer readable
```

## Sections

| section | content |
| :--- | :--- |
| Summary | cell types and reference cells; the cross-validated recall range of the best and unique calls (lowest and highest cell type); the share refused, by reason; the three weakest cell types |
| 1. Depth | genes detected per cell, per cell type |
| 2. UMAP | the reference: its labels, or the cross-validated best, unique or forced calls (the *output* switch); highlight and space switches ([the UMAPs](index.md#the-umaps)) |
| 3. Analytics per cell type | from the self-assessment: correct, wrong (and to which cell type) and refused by reason, for the forced, best and unique calls |
| 3b. Certainty and boundary cells | a UMAP of the cross-validated calls with views label, margin, ambiguity and unlabelled; the boundary map; per reference label the % best correct, % with a unique call, % boundary, % unlabelled by reason, median FDR and margin |
| 4. Cell type drill-down | a cell type's markers (highest detection in it first) and their detection rate in every cell type and in the background (`bg`), in a scrolling box; where its wrong calls went |
| 5. Weak-class flags | cell types weakest first by cross-validated best recall, with the evidence: most confused cell type, refused share, cells, a mark for fewer than 50 cells (a reporting convention, not a SCARLET threshold), median genes detected, markers, positive-evidence values; no threshold is applied |
| 6. Per-cell data | every reference cell: fold (1 to 5), label, out-of-fold calls, FDR, margin, runner-up, reasons, flags. A *Show* switch lists the cells worth a look first (wrong, refused, boundary); headers sort; 200 rows shown; the button *Download all cells (TSV)* builds the whole table in your browser |
| 7. Provenance | the run's summary; the model, technology preset, marker sizes, seed, blocklist and gene groups; package defaults or not; the reference and its SHA-256 (as the model recorded it, checked against the file read); the self-assessment settings; the UMAP spaces |

## The self-assessment

"How well does this model separate its own cell types?" The reference is split into five folds (`n_folds`;
stratified by label, shuffled with `seed`). Each fold is assigned by a model learned on the other four, with the
model's own recorded options; every reference cell is thus called by a model that never saw it (out-of-fold). It
costs one extra learn and assign per fold; `self_assessment=False` (`--no-self-assessment`) switches it off.

## Certainty and the runner-up

- **Margin** (view *margin*): the cell type's colour, darker the larger the call's margin, **relative to this run**:
  faintest at the 5th percentile of the margin over the run's best calls, full colour at the 95th; the caption gives
  the two margins. The margins themselves are in the hover and the per-cell table.
- **Ambiguity** (view *ambiguity*): the boundary cells (a best call and no unique call) in full colour, the others
  faint; seen on the best call (in the unique call they are the unassigned cells).
- **Margin** = the log-likelihood of the top cell type minus the larger of the runner-up's and the background's.
- **Runner-up** in this report: for a called cell, the cell type ranked second; for a cell refused for background,
  its **top cell type, the one the background beat** (its margin is then negative). (In a transfer, the runner-up is
  the cell type ranked after the top one for every cell: [Reading the output](../output.md#transfer).)
- **Boundary cells** have a best call but no unique call: SCARLET could not separate the best cell type from the
  runner-up at `alpha`. The boundary map counts them by best call × runner-up.

## Options

| option | default | what it does |
| :--- | :--- | :--- |
| `reference` | none | the reference, when the path the model records is no longer readable (in Python also an AnnData) |
| `self_assessment` | yes | the cross-validated self-assessment |
| `n_folds` | 5 | its folds |
| `seed` | 0 | the folds' and the UMAP's seed |
| `umap_cells` | 50000 | at most this many cells on the UMAPs |
| `umap_all_genes` | yes | also the UMAP on all genes (the second space) |
| `alpha` | 0.01 | advanced: the threshold of a unique call in the self-assessment |
| `ncpu` | 1 | CPU cores; it does not change the result |

Tables of the bundle: `colours`, `depth`, `self_assessment` (per cell), `analytics`, `wrong_to`, `boundary`,
`certainty_per_identity`, `umap` (`x`, `y` on the model's genes; `x_all`, `y_all` on all genes), `markers`,
`universe`, `weak_classes`.
