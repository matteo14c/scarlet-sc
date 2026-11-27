# How the benchmark was run

[Back to the benchmarks](index.md)

## How to read the tables

- **Every score: median F1 on assigned cells, with % unlabelled in brackets, both in %.** F1 is computed per class
  on the cells the output assigns; the median is taken over classes. The share of cells the output leaves
  unlabelled is always shown beside it (0.0 where nothing is refused), because an output can raise its F1 on
  assigned cells by refusing more cells.
- **Two views, side by side.** *Without abstention:* the outputs that label every cell (SCARLET best; Seurat predicted.id, SingleR labels, CellTypist best match, SVM, Harmony kNN). SCARLET best can refuse a
  cell that its background explains better than any cell type; its % unlabelled is shown. *With abstention:* the
  outputs that can refuse a cell (SCARLET unique; Seurat GTP, SingleR pruned, CellTypist prob match, SVMrejection,
  scmap-cluster, scmap-cell, scType). On a wide screen the two views sit side by side; on a narrow one, one under
  the other.
- **SCARLET 0.6.3 at its defaults, without a gene blocklist** (`blocklist none`), as every other tool, so that all
  tools saw the same genes. In normal use keep the blocklist on: SCARLET's default is `blocklist auto`.
- **SCARLET's outputs.** *forced:* the top cell type for every cell; never refuses. *best:* the top cell type, but
  refuses cells the background explains better than any cell type, and cells the flip test refuses (`assign`).
  On transfer, *forced* is the first pass's top label (`scarlet_top1`, the cell type ranked first most often
  across the draws, whatever the call): never refused, shown as the top-scoring label without refusal.
  *unique:* as best, and also refuses ambiguous cells (the top cell type is not significantly better than every
  other one).
- **First pass and refined (B and Part 1).** `transfer` first annotates the query with the reference's model (the
  first pass, `first_pass.tsv`), then re-learns on the query cells the first pass called uniquely and stably and
  re-annotates the query (refined, `labels.tsv`). Both are documented outputs, so both are shown.
- **Accurate and quick (transfer: Figs. 3a, 4, 5a, S8 and Fig. 6 Part 1).** *Accurate* is `transfer`'s default: 10 draws
  (seeds 1 to 10), and a final call needs the same outcome in at least 8 of the 10 draws, otherwise the cell is
  unassigned ("undetermined across draws"). The figures' tables show accurate. *Quick* is `--draws 1`: one draw, whose
  calls depend on its seed; each figure has one compact quick-mode table from 10 separate runs with seeds 1 to 10.
- **One section per figure panel**, each with its own rows, medians and better / equal / worse; the medians of F1
  and of % unlabelled are each taken over the panel's experiments, separately. The brain, pancreas and CellBench
  panels show each training / test direction, as the benchmark's figures do.
- **Better / equal / worse** is counted within each panel and view only: per experiment, SCARLET's F1 against each other output of
  the same view; F1 values equal to within 1e-9 count as equal. In brackets: the median % unlabelled over the
  experiments compared, SCARLET vs the tool. The count compares F1 only; it does not weigh one output's extra refusals against
  the other's.
- **Each tool only in the output modes its authors document**, as a user receives them. A multi-label output
  (CellTypist prob match "A|B") is a call, scored as that label string. No abstention mode is constructed for any
  tool.
- **Harmony kNN is transductive:** PCA and Harmony are fitted on reference and query together, so it sees the query
  before it labels it. It is flagged in every table.
- "not run": the output was not run on that task (see SVMrejection under [How each tool was run](#how-each-tool-was-run)).

## Provenance

### Code versions

Every SCARLET value on this page comes from SCARLET 0.6.3: its benchmark runs used a release candidate of 0.6.3 whose calls the release does not change. The package's `learn` and `assign` (A; abstention Parts
2–3) and `transfer` (B; abstention Part 1), every option at its default except `blocklist none`. The other tools'
versions are under [Tool versions](#tool-versions); their values are those of the benchmark runs described below.

### Gene blocklist

- **Every tool, SCARLET included, was run without a gene blocklist** (no genes removed). None of the other tools
  removes genes, so all tools saw the same genes.
- **In normal use keep it on.** SCARLET's default is `blocklist auto`: the bundled human or mouse list, chosen from
  the gene names.

### Folds

- **A tasks.** The benchmark's class filter: classes with 10 cells or fewer are removed. The 5 folds were regenerated
  with the rules and seed of its `Cross_Validation.R` script; the exact 2019 folds may differ, because R 3.6 changed
  `sample()`. Mapping to the benchmark: derived (its data, with this documented change).
- **B tasks.** The benchmark's own fold files, as distributed. Mapping: identical (its data and its folds).
- **Abstention.** The benchmark's own rejection files and folds; Part 3's seen-cell sample is ours (a stratified 20 %
  of each seen cell type, seeded).

### Transfer settings (B and abstention Part 1)

- Depth matching `capture` (the default): the deeper side is thinned towards the shallower one, in every draw with
  the draw's seed.
- Refinement: on (the default); stable anchors are the cells called uniquely to the same cell type in every draw.
- The flip test flags cells in `transfer` (the default) and refuses them in `assign` (the default).
- Accurate: 10 draws, seeds 1–10 (the default). Quick: `--draws 1`, 10 separate runs with `--seeds 1` to `--seeds 10`.

### SCARLET settings per task

Read from every run's own records: the technology preset SCARLET chose (sparse or dense, from the reference's
detected genes per cell; dense changes only the background penalty of the marker search) and the marker sizes.

| block | task | runs | technology | markers / shortlist |
| :--- | :--- | :--- | :--- | :--- |
| A | A/amb | 5 | dense | 50/150 |
| A | A/amb16 | 5 | dense | 50/150 |
| A | A/celseq2_5cl | 5 | dense | 50/150 |
| A | A/pancreas_muraro | 5 | dense | 50/150 |
| A | A/pancreas_segerstolpe | 5 | dense | 50/150 |
| A | A/pancreas_xin | 5 | dense | 50/150 |
| A | A/10x_5cl | 5 | sparse | 50/150 |
| A | A/pancreas_baron_human | 5 | sparse | 50/150 |
| A | A/pancreas_baron_mouse | 5 | sparse | 50/150 |
| A | A/zheng_sorted | 5 | sparse | 50/150 |
| A | A/tm | 5 | sparse | 50/150 |
| A | A/zheng_68k | 5 | sparse | 50/150 |
| B | B/Brain_Labels34_fold1 | 1 | dense | 50/150 |
| B | B/Brain_Labels34_fold2 | 1 | dense | 50/150 |
| B | B/Brain_Labels3_fold1 | 1 | dense | 50/150 |
| B | B/Brain_Labels3_fold2 | 1 | dense | 50/150 |
| B | B/CellBench | 2 | sparse | 50/150 |
| B | B/PbmcBench_10Xv2 | 1 | sparse | 50/150 |
| B | B/PbmcBench_10Xv3 | 1 | sparse | 50/150 |
| B | B/PbmcBench_CEL-Seq | 1 | sparse | 50/150 |
| B | B/PbmcBench_Drop-Seq | 1 | sparse | 50/150 |
| B | B/PbmcBench_Seq-Well | 1 | sparse | 50/150 |
| B | B/PbmcBench_Smart-Seq2 | 1 | sparse | 50/150 |
| B | B/PbmcBench_inDrop | 1 | sparse | 50/150 |
| B | NEW/Pancreatic_unaligned | 4 | sparse | 50/150 |
| abstention | REJ/part1_a | 1 | sparse | 50/150 |
| abstention | REJ/part1_b | 1 | sparse | 50/150 |
| abstention | REJ/part1_c | 1 | sparse | 50/150 |
| abstention | REJ/part1_d | 1 | sparse | 50/150 |
| abstention | part 2 | 3 | sparse | 50/150 |
| abstention | part 3 | 3 | sparse | 50/150 |

### Tool versions

As recorded in each tool's run provenance. `svm` covers SVM and SVMrejection; `scmap` covers scmap-cluster and
scmap-cell.

| tool | packages |
| :--- | :--- |
| seurat | Seurat 5.3.0; SeuratObject 5.2.0 |
| celltypist | anndata 0.11.4; celltypist 1.7.1; python 3.10.20; scanpy 1.11.5; scikit-learn 1.7.2 |
| singler | Seurat 5.3.0; SeuratObject 5.2.0; SingleCellExperiment 1.24.0; SingleR 2.4.0; python 3.14.6 |
| scmap | SingleCellExperiment 1.24.0; python 3.14.6; scmap 1.24.0 |
| sctype | Seurat 5.3.0; SeuratObject 5.2.0 |
| svm | anndata 0.12.14; harmonypy 2.0.0; python 3.11.15; scanpy 1.11.5; scikit-learn 1.8.0 |
| harmony_knn | anndata 0.12.14; harmonypy 2.0.0; python 3.11.15; scanpy 1.11.5; scikit-learn 1.8.0 |

### Input rule (all tools)

- Where a tool's documentation specifies its input, it was followed.
- Where it does not, the benchmark's published wrapper was followed.
- Any other input or reading of a tool's output is not shown in these tables.

### How each tool was run

| tool | input | settings | outputs shown |
| :--- | :--- | :--- | :--- |
| SCARLET 0.6.3 | counts as given; detection (count > 0) | the package's `learn` + `assign` (A, abstention Parts 2–3) and `transfer` (B, abstention Part 1), every option at its default except `blocklist none`; `transfer` accurate (default) and quick (`--draws 1`, seeds 1–10) | forced, best, unique; on transfer, forced = the first pass's top label, and first pass and refined best and unique |
| Seurat | Seurat's mapping vignette: `NormalizeData` (LogNormalize), `FindVariableFeatures`, `ScaleData`, `RunPCA` | `FindTransferAnchors`, `TransferData`, as in the vignette | `predicted.id` (never refuses); GTP (see below) |
| CellTypist | log1p of counts per 10,000, as its documentation requires; all genes | `celltypist.train` at its defaults; `annotate` | best match (never refuses); prob match (default threshold 0.5) |
| SingleR | log2(counts per 10,000 + 1) for reference and query (no scale is prescribed by its documentation) | `de.method = "wilcox"`, as its vignette recommends for single-cell references | `labels` (never refuses); `pruned.labels` (default pruning) |
| scmap-cluster | the benchmark's wrapper: raw counts as `normcounts`, `logcounts = log2(normcounts + 1)` (the scmap vignette prescribes no normalisation) | defaults (`scmapCluster`, similarity threshold 0.7) | its only output (refuses as "unassigned") |
| scmap-cell | as scmap-cluster | `indexCell` defaults, seeded per fold; `scmapCell`; `scmapCell2Cluster` defaults (w = 3, threshold 0.5) | its only output (refuses as "unassigned") |
| scType | its README workflow on each test fold: Seurat `NormalizeData`, variable features, `ScaleData`, PCA, clustering | `sctype_score` on the scaled data, one type per cluster, scType's "Unknown" rule (cluster score below a quarter of its cells); markers derived from the training data (no scType database entry exists for these labels) | its only output |
| SVM | the benchmark's wrapper: log1p of the values; all genes | the benchmark's SVM script (`LinearSVC`) | never refuses |
| SVMrejection | as SVM | the benchmark's SVM-rejection script: `LinearSVC` calibrated with `CalibratedClassifierCV` (3 folds, as its 2019 runs); "Unknown" below a calibrated probability of 0.7 | refuses as "Unknown" |
| Harmony kNN | log-normalised counts | not a published classifier: PCA and Harmony fitted on reference and query together, then a kNN classifier | never refuses; **transductive** (flagged in every table) |

- **SVMrejection was not run** on `B/Brain_MouseV1_MouseALM_Labels34`. There, a class has 2 training cells, and
  current scikit-learn versions refuse 3-fold calibration; the benchmark's 2019 runs used scikit-learn 0.21, which
  did not stop. Plain SVM was run.
- **Seurat GTP** (`GetTransferPredictions`, `score.filter` 0.75) is defined by Seurat's own documentation, so it is
  shown in the tables. Seurat 5.3.0's reference manual, help page "Get the predicted identity": "Utility function to
  easily pull out the name of the class with the maximum prediction. This is useful if you've set
  prediction.assay = TRUE in TransferData and want to have a vector with the predicted class". Its argument
  `score.filter`: "Return "Unassigned" for any cell with a score less than this value" (default 0.75). It is an
  exported function. It is not in Seurat's mapping vignette (which uses `TransferData`'s `predicted.id` only, with no
  threshold) and not in Seurat's NEWS.md. No published benchmark is its source. Choosing it as Seurat's abstaining
  output was our choice, under the rule that a tool abstains only through a mode the tool itself provides, at that
  mode's defaults. The refusal rule and its threshold (0.75) are Seurat's. Seurat `predicted.id` (the vignette's
  output) is shown beside it.

On the figure pages SCARLET appears once in each view: **SCARLET best** (without abstention) and **SCARLET unique** (with abstention); for `transfer`, these are the final calls, after refinement, in accurate mode (the default). Its other outputs (forced, first pass, quick mode) are on [a separate page](further.md).
