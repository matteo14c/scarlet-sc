# Fig. 4a and 4b: mouse brain, 3 and 34 cell populations

[Back to the benchmarks](index.md) · [how to read the tables](methods.md#how-to-read-the-tables) · [SCARLET's other outputs](further.md)

## Fig. 4a: brain, 3 cell populations

> Classification performance across brain datasets. Heatmaps show the median F1-scores of the supervised classifiers when tested on a major lineage annotation with three cell populations and b deeper level of annotation with 34 cell populations. The training sets are indicated using the column labels on top of the heatmap. The test set is indicated in the gray box.
>
> — Abdelaal et al. (2019), the legend of Fig. 4

One row per experiment: the two directions between the mouse brain datasets (VISp, ALM), training → test, 3 populations. Each value: median F1 (% unlabelled). Code: SCARLET 0.6.3, `transfer`, accurate (see [Transfer settings](methods.md#transfer-settings-b-and-abstention-part-1)); the other tools as in [How each tool was run](methods.md#how-each-tool-was-run).

<div class="views" markdown>

<div markdown>

**Without abstention: outputs that label every cell**

| training → test | SCARLET best | Seurat predicted.id | SingleR labels | CellTypist best match | SVM | Harmony kNN (transductive) |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| VISp → ALM (Mouse) | 100.0 (0.8) | 100.0 (0.0) | 100.0 (0.0) | 99.8 (0.0) | 100.0 (0.0) | 100.0 (0.0) |
| ALM → VISp (Mouse) | 100.0 (1.7) | 100.0 (0.0) | 100.0 (0.0) | 100.0 (0.0) | 100.0 (0.0) | 100.0 (0.0) |

</div>

<div markdown>

**With abstention: outputs that can refuse a cell**

| training → test | SCARLET unique | Seurat GTP | SingleR pruned | CellTypist prob match | SVMrejection | scmap-cluster | scmap-cell | scType |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| VISp → ALM (Mouse) | 100.0 (1.1) | 100.0 (0.0) | 100.0 (4.5) | 99.5 (0.0) | 100.0 (0.0) | 100.0 (0.3) | 99.9 (0.2) | 99.5 (0.0) |
| ALM → VISp (Mouse) | 100.0 (2.5) | 100.0 (0.0) | 100.0 (1.6) | 100.0 (0.0) | 100.0 (0.0) | 100.0 (0.3) | 100.0 (0.0) | 99.9 (8.0) |

</div>

</div>

### Fig. 4a: medians over the 2 experiments

<div class="views" markdown>

<div markdown>

**Without abstention**

| output | experiments | median F1 (% unlabelled) |
| :--- | :--- | :--- |
| SCARLET best | 2 | 100.0 (1.3) |
| Seurat predicted.id | 2 | 100.0 (0.0) |
| SingleR labels | 2 | 100.0 (0.0) |
| CellTypist best match | 2 | 99.9 (0.0) |
| SVM | 2 | 100.0 (0.0) |
| Harmony kNN (transductive) | 2 | 100.0 (0.0) |

</div>

<div markdown>

**With abstention**

| output | experiments | median F1 (% unlabelled) |
| :--- | :--- | :--- |
| SCARLET unique | 2 | 100.0 (1.8) |
| Seurat GTP | 2 | 100.0 (0.0) |
| SingleR pruned | 2 | 100.0 (3.1) |
| CellTypist prob match | 2 | 99.7 (0.0) |
| SVMrejection | 2 | 100.0 (0.0) |
| scmap-cluster | 2 | 100.0 (0.3) |
| scmap-cell | 2 | 100.0 (0.1) |
| scType | 2 | 99.7 (4.0) |

</div>

</div>

### Fig. 4a: better / equal / worse

<div class="views" markdown>

<div markdown>

**Without abstention**

| SCARLET output | Seurat predicted.id | SingleR labels | CellTypist best match | SVM | Harmony kNN (transductive) |
| :--- | :--- | :--- | :--- | :--- | :--- |
| SCARLET best | 0 / 2 / 0 (1.3 vs 0.0) | 0 / 2 / 0 (1.3 vs 0.0) | 2 / 0 / 0 (1.3 vs 0.0) | 0 / 2 / 0 (1.3 vs 0.0) | 0 / 2 / 0 (1.3 vs 0.0) |

</div>

<div markdown>

**With abstention**

| SCARLET output | Seurat GTP | SingleR pruned | CellTypist prob match | SVMrejection | scmap-cluster | scmap-cell | scType |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| SCARLET unique | 0 / 2 / 0 (1.8 vs 0.0) | 0 / 2 / 0 (1.8 vs 3.1) | 2 / 0 / 0 (1.8 vs 0.0) | 0 / 2 / 0 (1.8 vs 0.0) | 0 / 2 / 0 (1.8 vs 0.3) | 2 / 0 / 0 (1.8 vs 0.1) | 2 / 0 / 0 (1.8 vs 4.0) |

</div>

</div>

## Fig. 4b: brain, 34 cell populations

> Classification performance across brain datasets. Heatmaps show the median F1-scores of the supervised classifiers when tested on a major lineage annotation with three cell populations and b deeper level of annotation with 34 cell populations. The training sets are indicated using the column labels on top of the heatmap. The test set is indicated in the gray box.
>
> — Abdelaal et al. (2019), the legend of Fig. 4

One row per experiment: the two directions between the mouse brain datasets (VISp, ALM), training → test, 34 populations. Each value: median F1 (% unlabelled). Code: SCARLET 0.6.3, `transfer`, accurate (see [Transfer settings](methods.md#transfer-settings-b-and-abstention-part-1)); the other tools as in [How each tool was run](methods.md#how-each-tool-was-run). SVMrejection was not run on these experiments (a cell type has 2 training cells; see there).

<div class="views" markdown>

<div markdown>

**Without abstention: outputs that label every cell**

| training → test | SCARLET best | Seurat predicted.id | SingleR labels | CellTypist best match | SVM | Harmony kNN (transductive) |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| VISp → ALM (Mouse) | 91.3 (0.8) | 93.9 (0.0) | 91.8 (0.0) | 93.2 (0.0) | 94.6 (0.0) | 92.6 (0.0) |
| ALM → VISp (Mouse) | 92.4 (2.6) | 94.6 (0.0) | 94.1 (0.0) | 93.1 (0.0) | 95.0 (0.0) | 93.7 (0.0) |

</div>

<div markdown>

**With abstention: outputs that can refuse a cell**

| training → test | SCARLET unique | Seurat GTP | SingleR pruned | CellTypist prob match | SVMrejection | scmap-cluster | scmap-cell | scType |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| VISp → ALM (Mouse) | 92.1 (2.7) | 98.2 (20.3) | 92.8 (1.1) | 93.4 (8.6) | not run | 91.5 (5.2) | 94.9 (38.5) | 61.2 (0.0) |
| ALM → VISp (Mouse) | 93.2 (3.7) | 98.6 (10.7) | 94.3 (1.5) | 94.4 (18.3) | not run | 93.1 (1.9) | 97.0 (20.5) | 51.9 (0.0) |

</div>

</div>

### Fig. 4b: medians over the 2 experiments

<div class="views" markdown>

<div markdown>

**Without abstention**

| output | experiments | median F1 (% unlabelled) |
| :--- | :--- | :--- |
| SCARLET best | 2 | 91.9 (1.7) |
| Seurat predicted.id | 2 | 94.3 (0.0) |
| SingleR labels | 2 | 93.0 (0.0) |
| CellTypist best match | 2 | 93.2 (0.0) |
| SVM | 2 | 94.8 (0.0) |
| Harmony kNN (transductive) | 2 | 93.1 (0.0) |

</div>

<div markdown>

**With abstention**

| output | experiments | median F1 (% unlabelled) |
| :--- | :--- | :--- |
| SCARLET unique | 2 | 92.6 (3.2) |
| Seurat GTP | 2 | 98.4 (15.5) |
| SingleR pruned | 2 | 93.5 (1.3) |
| CellTypist prob match | 2 | 93.9 (13.4) |
| SVMrejection | 0 | not run |
| scmap-cluster | 2 | 92.3 (3.6) |
| scmap-cell | 2 | 95.9 (29.5) |
| scType | 2 | 56.6 (0.0) |

</div>

</div>

### Fig. 4b: better / equal / worse

<div class="views" markdown>

<div markdown>

**Without abstention**

| SCARLET output | Seurat predicted.id | SingleR labels | CellTypist best match | SVM | Harmony kNN (transductive) |
| :--- | :--- | :--- | :--- | :--- | :--- |
| SCARLET best | 0 / 0 / 2 (1.7 vs 0.0) | 0 / 0 / 2 (1.7 vs 0.0) | 0 / 0 / 2 (1.7 vs 0.0) | 0 / 0 / 2 (1.7 vs 0.0) | 0 / 0 / 2 (1.7 vs 0.0) |

</div>

<div markdown>

**With abstention**

| SCARLET output | Seurat GTP | SingleR pruned | CellTypist prob match | SVMrejection | scmap-cluster | scmap-cell | scType |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| SCARLET unique | 0 / 0 / 2 (3.2 vs 15.5) | 0 / 0 / 2 (3.2 vs 1.3) | 0 / 0 / 2 (3.2 vs 13.4) | — | 2 / 0 / 0 (3.2 vs 3.6) | 0 / 0 / 2 (3.2 vs 29.5) | 2 / 0 / 0 (3.2 vs 0.0) |

</div>

</div>
