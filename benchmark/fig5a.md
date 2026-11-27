# Fig. 5a: pancreas, four datasets

[Back to the benchmarks](index.md) · [how to read the tables](methods.md#how-to-read-the-tables) · [SCARLET's other outputs](further.md)

## Fig. 5a: pancreas, unaligned

> Classification performance across pancreatic datasets. Heatmaps showing the median F1-score for each classifier for the a unaligned and b aligned datasets. The column labels indicate which of the four datasets was used as a test set, in which case the other three datasets were used as training.
>
> — Abdelaal et al. (2019), the legend of Fig. 5

One row per experiment: the four pancreas datasets, each the test set once, the other three the training set. Each value: median F1 (% unlabelled). Code: SCARLET 0.6.3, `transfer`, accurate (see [Transfer settings](methods.md#transfer-settings-b-and-abstention-part-1)); the other tools as in [How each tool was run](methods.md#how-each-tool-was-run).

<div class="views" markdown>

<div markdown>

**Without abstention: outputs that label every cell**

| training → test | SCARLET best | Seurat predicted.id | SingleR labels | CellTypist best match | SVM | Harmony kNN (transductive) |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| the other three → Xin | 100.0 (12.1) | 99.5 (0.0) | 99.4 (0.0) | 99.3 (0.0) | 98.9 (0.0) | 98.3 (0.0) |
| the other three → Segerstolpe | 99.8 (13.6) | 99.1 (0.0) | 98.9 (0.0) | 99.4 (0.0) | 99.8 (0.0) | 99.5 (0.0) |
| the other three → Muraro | 98.5 (6.2) | 97.0 (0.0) | 93.6 (0.0) | 94.3 (0.0) | 96.8 (0.0) | 95.5 (0.0) |
| the other three → Baron human | 99.4 (11.8) | 96.8 (0.0) | 94.7 (0.0) | 98.1 (0.0) | 98.6 (0.0) | 97.2 (0.0) |

</div>

<div markdown>

**With abstention: outputs that can refuse a cell**

| training → test | SCARLET unique | Seurat GTP | SingleR pruned | CellTypist prob match | SVMrejection | scmap-cluster | scmap-cell | scType |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| the other three → Xin | 100.0 (18.8) | 99.6 (11.0) | 99.9 (1.0) | 99.2 (0.6) | 99.3 (0.3) | 48.8 (38.5) | 99.6 (1.0) | 98.1 (0.0) |
| the other three → Segerstolpe | 99.8 (18.8) | 99.4 (0.8) | 100.0 (2.6) | 98.2 (0.1) | 99.9 (1.2) | 100.0 (36.7) | 96.4 (15.8) | 99.8 (2.6) |
| the other three → Muraro | 98.9 (8.8) | 97.9 (2.1) | 95.9 (2.9) | 80.5 (0.0) | 97.9 (1.2) | 31.7 (32.7) | 97.5 (11.0) | 95.7 (0.0) |
| the other three → Baron human | 99.8 (21.2) | 99.4 (6.1) | 95.2 (0.4) | 98.8 (1.6) | 100.0 (69.7) | 99.9 (88.5) | 99.3 (12.6) | 98.7 (1.9) |

</div>

</div>

### Fig. 5a: medians over the 4 experiments

<div class="views" markdown>

<div markdown>

**Without abstention**

| output | experiments | median F1 (% unlabelled) |
| :--- | :--- | :--- |
| SCARLET best | 4 | 99.6 (11.9) |
| Seurat predicted.id | 4 | 98.0 (0.0) |
| SingleR labels | 4 | 96.8 (0.0) |
| CellTypist best match | 4 | 98.7 (0.0) |
| SVM | 4 | 98.7 (0.0) |
| Harmony kNN (transductive) | 4 | 97.8 (0.0) |

</div>

<div markdown>

**With abstention**

| output | experiments | median F1 (% unlabelled) |
| :--- | :--- | :--- |
| SCARLET unique | 4 | 99.8 (18.8) |
| Seurat GTP | 4 | 99.4 (4.1) |
| SingleR pruned | 4 | 97.9 (1.8) |
| CellTypist prob match | 4 | 98.5 (0.3) |
| SVMrejection | 4 | 99.6 (1.2) |
| scmap-cluster | 4 | 74.3 (37.6) |
| scmap-cell | 4 | 98.4 (11.8) |
| scType | 4 | 98.4 (0.9) |

</div>

</div>

### Fig. 5a: better / equal / worse

<div class="views" markdown>

<div markdown>

**Without abstention**

| SCARLET output | Seurat predicted.id | SingleR labels | CellTypist best match | SVM | Harmony kNN (transductive) |
| :--- | :--- | :--- | :--- | :--- | :--- |
| SCARLET best | 4 / 0 / 0 (11.9 vs 0.0) | 4 / 0 / 0 (11.9 vs 0.0) | 4 / 0 / 0 (11.9 vs 0.0) | 4 / 0 / 0 (11.9 vs 0.0) | 4 / 0 / 0 (11.9 vs 0.0) |

</div>

<div markdown>

**With abstention**

| SCARLET output | Seurat GTP | SingleR pruned | CellTypist prob match | SVMrejection | scmap-cluster | scmap-cell | scType |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| SCARLET unique | 4 / 0 / 0 (18.8 vs 4.1) | 3 / 0 / 1 (18.8 vs 1.8) | 4 / 0 / 0 (18.8 vs 0.3) | 2 / 0 / 2 (18.8 vs 1.2) | 2 / 0 / 2 (18.8 vs 37.6) | 4 / 0 / 0 (18.8 vs 11.8) | 3 / 0 / 1 (18.8 vs 0.9) |

</div>

</div>
