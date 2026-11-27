# The transfer report

[Reports](index.md) · [learn report](learn.md) · [assign report](assign.md)

The report of a transfer output folder. A is the reference, B the query; the title names both by their file stems
(for example *smartseq2_reference → 10xv2_query*).

```bash
scarlet report transfer_out/                              # reads the reference and query from the recorded paths
scarlet report transfer_out/ --query-labels cell_type     # also compare the calls with labels you know
```

Three pages: the main page (`report.html`); *what refinement changed* (`refinement.html`), linked from section 6;
and the supplement (`supplement_A/report.html`), the [learn report](learn.md) of the transfer's reference model with
its self-assessment. If the reference was thinned by depth matching, the supplement says so; its self-assessment
learns on the reference as given.

## Sections

| section | content |
| :--- | :--- |
| Summary | A → B with their cells and cell types; the share of B called (best, unique) and refused by reason; depth matching in one line (which side was thinned, from and to how many detected genes per cell, median over the draws); the draws and the 8-of-10 rule; the transfer check; cell types of A with no call in B |
| 1. Depth | genes detected per cell in A and B before depth matching (boxes), and the record's medians on the shared genes before and after matching (after: the median over the draws with their range) |
| 2. UMAPs | A by its labels and B side by side. B: output (final best, final unique, first pass best, first pass unique) and view (label, agreement, ambiguity, anchors, unlabelled); both: highlight and space ([the UMAPs](index.md#the-umaps)) |
| 3. Analytics | final calls per cell type (best solid, unique hatched, unassigned grey; the count and share in the hover); refusals by reason (background, ambiguity, undetermined across draws, flip test); cell types of A with no call in B |
| 3b. Certainty, boundary cells and refusals | the boundary maps (first pass and refined) and the per-identity table |
| 3c. Against labels supplied by the user | only with `query_labels`: for each output, over the cells whose label the model knows, the share refused and the median F1 on the called cells; per label (labels the model does not know listed last, with what their cells were called); a map of labels against the final best calls |
| 4. Cell type drill-down | two heatmaps in one gene order and three parts: the markers common to A's model and B's refined model, those only in A's (dropped by refinement), those only in B's (added); left their detection in A's model (A's cells per label), right in B's refined model (B's stable anchors per cell type), with the background column `bg`. Every value is the model's own estimate, read from its in-class table (`in_class.tsv` in the model: learn's bootstrap estimate of the share of each cell type's cells detecting each gene, the same estimate the model uses for its markers). A model written before SCARLET 0.6.3 has no such table: the other model's markers are then blank, *not a marker of this model* in the hover. Each marker's rank in both models (common markers joined by a line); each reference marker's detection in A against the B cells called that cell type; stable anchors |
| 5. Weak-class flags | from the reference (the supplement's self-assessment) and from the transfer (ranked by the share of reference markers confirmed, lowest first, with stable anchors and *doubtful transfer*) |
| 6. What refinement changed | cells refused then called, changed cell type, and called then refused (best and unique); the details on their own page, `refinement.html`: first-pass and refined calls side by side, a map of every change of call (first-pass call against final call) and the list of changes, largest first |
| 7. Agreement across the draws | the rule, a histogram of the final calls' agreement (best and unique; at least 8 of 10 in blue, the final call; below in grey, unassigned) and its share by final best call; each draw's own result |
| 8. Per-cell data | B's cells counted by how their final call came out (refused, boundary, called below full agreement, called at full agreement); the button *Download all cells (TSV)* builds the whole table in your browser; *Browse the cells* opens the table: final and first-pass calls, agreement, reasons, top identity, runner-up (and in how many draws), stable anchor, flags, with a *Show* switch (the cells worth a look first), headers that sort, 25 rows at a time |
| 9. Provenance | the run's summary; both sides and their SHA-256 (as transfer recorded them, checked against the files read); refinement, depth matching, draws and seeds, technology preset, flip test, alpha, blocklist; package defaults or not; the UMAP spaces; depth matching draw by draw |

## Agreement, certainty and anchors

- **Agreement** of a cell: how many of the draws gave its most frequent outcome (`unassigned` counts as an outcome).
  The final call needs at least 8 of 10 draws; below that the cell is unassigned, *undetermined across draws*. It
  measures how stable a call is under the random depth thinning and seeds of the draws, not whether it is correct.
- **Agreement** (view *agreement* of B): the call's colour by how many draws agree: faintest at 8 of 10 (the least a
  final call needs), full at 10 of 10. Quick mode has one draw and no agreement, and no such view.
- **Ambiguity** (view *ambiguity* of B): the boundary cells (a best call and no unique call) in full colour, the
  others faint; seen on the best calls.
- **Stable anchors** (view *anchors*): the query cells called uniquely to the same cell type in every draw; a cell
  type with at least 3 is learned again from them (refined).
- **Runner-up:** the cell type ranked after the top one, in the most draws (`scarlet_runner_up`), for every cell.
- **Per identity:** cells called best and unique; stable anchors; refined or not, and why; of its best calls, % at
  full agreement and % boundary; among the cells it ranks first, % unlabelled by reason.

## Options

`reference`, `query` (when the recorded paths are no longer readable), `query_labels`, `output`,
`self_assessment` (the supplement's), `n_folds`, `seed`, `umap_cells`, `umap_all_genes`, `ncpu`: see
[the learn report](learn.md#options) and [Reports](index.md).

Tables of the bundle: `anchors`, `colours`, `depth`, `depth_match`, `umap_A`, `umap_B`, `query_umap`, `analytics`,
`no_calls`, `markers`, `universe`, `in_class` (the models' in-class tables on the drill-down's genes), `marker_A_vs_B`, `qc`, `changed_per_identity`, `refinement_changes`,
`first_pass_vs_refined`, `query_cells`, `boundary`, `per_identity`, `agreement`, `per_draw`; with `query_labels` also
`known_summary`, `known_per_label` and `known_confusion`.
