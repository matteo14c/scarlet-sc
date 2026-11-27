# Supp. Fig. S8: CellBench

[Back to the benchmarks](index.md) · [how to read the tables](methods.md#how-to-read-the-tables) · [SCARLET's other outputs](further.md)

## Supp. Fig. S8: CellBench

> Figure S8. Classification performance across the CellBench datasets. Heatmaps show the (A) median F1-score and (B) percentage of unlabeled cells across the CellBench datasets. The training set is indicated above the heatmap, the test set below.
>
> — Abdelaal et al. (2019), the legend of Additional file 1, Fig. S8

One row per experiment: the two directions between the CellBench datasets (10x Chromium, CEL-Seq2), training → test. Each value: median F1 (% unlabelled). Code: SCARLET 0.6.3, `transfer`, accurate (see [Transfer settings](methods.md#transfer-settings-b-and-abstention-part-1)); the other tools as in [How each tool was run](methods.md#how-each-tool-was-run).

<div class="views" markdown>

<div markdown>

**Without abstention: outputs that label every cell**

| training → test | SCARLET best | Seurat predicted.id | SingleR labels | CellTypist best match | SVM | Harmony kNN (transductive) |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| 10X → CEL-Seq2 | 100.0 (1.6) | 100.0 (0.0) | 99.4 (0.0) | 100.0 (0.0) | 100.0 (0.0) | 100.0 (0.0) |
| CEL-Seq2 → 10X | 100.0 (1.6) | 100.0 (0.0) | 100.0 (0.0) | 100.0 (0.0) | 100.0 (0.0) | 100.0 (0.0) |

</div>

<div markdown>

**With abstention: outputs that can refuse a cell**

| training → test | SCARLET unique | Seurat GTP | SingleR pruned | CellTypist prob match | SVMrejection | scmap-cluster | scmap-cell | scType |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| 10X → CEL-Seq2 | 100.0 (3.5) | 100.0 (0.0) | 100.0 (2.6) | 99.3 (0.0) | 100.0 (0.0) | 100.0 (0.0) | 100.0 (0.0) | 100.0 (0.0) |
| CEL-Seq2 → 10X | 100.0 (2.1) | 100.0 (0.0) | 100.0 (1.3) | 99.9 (0.1) | 100.0 (0.8) | 100.0 (5.3) | 100.0 (0.0) | 100.0 (0.0) |

</div>

</div>

### Supp. Fig. S8: medians over the 2 experiments

<div class="views" markdown>

<div markdown>

**Without abstention**

| output | experiments | median F1 (% unlabelled) |
| :--- | :--- | :--- |
| SCARLET best | 2 | 100.0 (1.6) |
| Seurat predicted.id | 2 | 100.0 (0.0) |
| SingleR labels | 2 | 99.7 (0.0) |
| CellTypist best match | 2 | 100.0 (0.0) |
| SVM | 2 | 100.0 (0.0) |
| Harmony kNN (transductive) | 2 | 100.0 (0.0) |

</div>

<div markdown>

**With abstention**

| output | experiments | median F1 (% unlabelled) |
| :--- | :--- | :--- |
| SCARLET unique | 2 | 100.0 (2.8) |
| Seurat GTP | 2 | 100.0 (0.0) |
| SingleR pruned | 2 | 100.0 (2.0) |
| CellTypist prob match | 2 | 99.6 (0.1) |
| SVMrejection | 2 | 100.0 (0.4) |
| scmap-cluster | 2 | 100.0 (2.6) |
| scmap-cell | 2 | 100.0 (0.0) |
| scType | 2 | 100.0 (0.0) |

</div>

</div>

### Supp. Fig. S8: better / equal / worse

<div class="views" markdown>

<div markdown>

**Without abstention**

| SCARLET output | Seurat predicted.id | SingleR labels | CellTypist best match | SVM | Harmony kNN (transductive) |
| :--- | :--- | :--- | :--- | :--- | :--- |
| SCARLET best | 0 / 2 / 0 (1.6 vs 0.0) | 1 / 1 / 0 (1.6 vs 0.0) | 0 / 2 / 0 (1.6 vs 0.0) | 0 / 2 / 0 (1.6 vs 0.0) | 0 / 2 / 0 (1.6 vs 0.0) |

</div>

<div markdown>

**With abstention**

| SCARLET output | Seurat GTP | SingleR pruned | CellTypist prob match | SVMrejection | scmap-cluster | scmap-cell | scType |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| SCARLET unique | 0 / 2 / 0 (2.8 vs 0.0) | 0 / 2 / 0 (2.8 vs 2.0) | 2 / 0 / 0 (2.8 vs 0.1) | 0 / 2 / 0 (2.8 vs 0.4) | 0 / 2 / 0 (2.8 vs 2.6) | 0 / 2 / 0 (2.8 vs 0.0) | 0 / 2 / 0 (2.8 vs 0.0) |

</div>

</div>
