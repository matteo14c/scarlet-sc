# Known limitations

- **Small classes.** A cell type with few reference cells (tens or fewer) gets markers from very few cells; check
  the learn report's weak-class table before trusting its calls.
- **A population the reference lacks.** SCARLET refuses a cell when the background explains it as well as any cell
  type of the reference. A missing population *distant* from every cell type of the reference is refused; one
  *close* to a cell type it has is mostly called as that relative (in tutorial 3, natural killer cells, absent from
  the reference, are called cytotoxic T cells). Look for such cells in the report's UMAP on all genes, or use a reference that has
  them ([Benchmarks](https://github.com/matteo14c/scarlet-sc/blob/v0.6.3/benchmark/fig6.md), tutorial 3).
- **Depth mismatch.** A reference and a query of very different depth are what `transfer` (depth matching) is for;
  `assign` does not match depths: it scores the query at its own depth.
- **Detection only.** SCARLET scores whether genes are detected, not how much they are expressed, so cell types
  that differ mainly in expression level are hard for it to separate.
- **Refinement and refusals.** In `transfer`, the refinement can turn cells refused in the first pass into calls, and
  calls into refusals. The report shows both; the first-pass calls stay in `first_pass.tsv`.
- **A query dominated by one population, with few stable anchors.** The refinement learns each cell type again from
  the query's own stable anchors and scores against a background pooled over all the query's cells. When the query is
  small and one population dominates it, and that population has few stable anchors, its cells can lose to the
  background and be refused. Example (the benchmark's 10x → Smart-seq2 pair, 253 query cells): the 10x reference has
  natural killer, dendritic and plasmacytoid dendritic cells, which the Smart-seq2 annotation folds into other labels;
  so in the first pass many T cells are undetermined across draws; so few become stable anchors (7 CD4+ T, 13
  cytotoxic T); the pooled background is 69 % T cells, and after refinement about half of the cells are refused for
  background, while the calls made are correct. Planned for 0.6.4: a check on the stable anchors before a cell type is
  refined, and a background balanced across the first pass's cell types.
- **Dominant cell types refused more after refinement.** For the same reason, the cells of the largest cell types
  must fit their type clearly better than a background that resembles them; more of them are refused than in the first
  pass (in the pancreas benchmark, alpha and beta cells); the cells called stay correct.
- **Closely related cell types.** Cells that fit two cell types almost equally (for example CD4 and CD8 T cells in
  some protocols) are often refused, in the first pass and after refinement.
- **Stable anchors keep their label.** A cell the first pass called the same cell type in every draw keeps that
  label, even when the query holds no cell type of the reference: in a query of foreign cells such anchors are wrong
  calls by construction.
- **Draws, run time and memory.** A transfer repeats its random steps in 10 draws and keeps a cell's call when at
  least 8 of them agree; with one CPU it takes about five times as long as one draw (82 s and 0.5 GB on a PBMC
  transfer, 3.7 min and 3.1 GB on a deep mouse-brain one, without the report). `ncpu` runs draws side by side and
  needs more memory: see [Run time, memory and CPUs](guide/transfer.md#run-time-memory-and-cpus). Quick mode
  (`draws=1`) runs once, with seed-dependent calls and no consensus. Origins with fewer than 2 cell types anchored
  stably keep their first-pass calls.
- **Hard transfers.** Across technologies, between datasets of very different sizes, or among closely related
  populations, SCARLET refuses more cells rather than guess.
- **`learn` and `assign` run once.** Their marker search is random too (seeded), and on some datasets a different
  seed changes the call of a sizeable share of cells (about a fifth of the cells of the Zheng 68K PBMC dataset);
  they are not repeated in draws in this release.
- **Cross-species** annotation is out of scope.
