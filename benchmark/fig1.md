# Fig. 1: annotation within a dataset

[Back to the benchmarks](index.md) · [how to read the tables](methods.md#how-to-read-the-tables) · [SCARLET's other outputs](further.md)

## Fig. 1: annotation within a dataset

The benchmark's Fig. 1 (A): 5-fold cross-validation within each of 12 datasets. Code: SCARLET 0.6.3 (see
[Code versions](methods.md#code-versions)).

<div class="views" markdown>

<div markdown>

**Without abstention: outputs that label every cell**

| task | SCARLET best | Seurat predicted.id | SingleR labels | CellTypist best match | SVM | Harmony kNN (transductive) |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| A/10x_5cl | 99.9 (1.1) | 100.0 (0.0) | 100.0 (0.0) | 100.0 (0.0) | 100.0 (0.0) | 100.0 (0.0) |
| A/amb | 100.0 (0.3) | 100.0 (0.0) | 100.0 (0.0) | 99.9 (0.0) | 100.0 (0.0) | 100.0 (0.0) |
| A/amb16 | 98.8 (0.8) | 99.3 (0.0) | 98.6 (0.0) | 99.1 (0.0) | 99.4 (0.0) | 99.2 (0.0) |
| A/celseq2_5cl | 100.0 (1.2) | 100.0 (0.0) | 100.0 (0.0) | 100.0 (0.0) | 100.0 (0.0) | 100.0 (0.0) |
| A/pancreas_baron_human | 96.9 (2.2) | 97.5 (0.0) | 96.9 (0.0) | 98.0 (0.0) | 98.2 (0.0) | 96.9 (0.0) |
| A/pancreas_baron_mouse | 100.0 (2.0) | 96.8 (0.0) | 96.8 (0.0) | 98.6 (0.0) | 98.2 (0.0) | 98.6 (0.0) |
| A/pancreas_muraro | 97.0 (1.7) | 96.3 (0.0) | 96.0 (0.0) | 96.6 (0.0) | 96.7 (0.0) | 95.9 (0.0) |
| A/pancreas_segerstolpe | 99.0 (2.1) | 98.2 (0.0) | 99.2 (0.0) | 99.0 (0.0) | 99.7 (0.0) | 99.2 (0.0) |
| A/pancreas_xin | 100.0 (3.1) | 98.0 (0.0) | 98.5 (0.0) | 95.3 (0.0) | 98.9 (0.0) | 99.6 (0.0) |
| A/tm | 95.2 (1.9) | 94.3 (0.0) | 91.9 (0.0) | 96.3 (0.0) | 98.0 (0.0) | 95.6 (0.0) |
| A/zheng_68k | 63.2 (2.5) | 66.2 (0.0) | 49.6 (0.0) | 68.3 (0.0) | 69.9 (0.0) | 68.0 (0.0) |
| A/zheng_sorted | 71.1 (1.2) | 82.1 (0.0) | 72.9 (0.0) | 94.0 (0.0) | 95.1 (0.0) | 89.0 (0.0) |

</div>

<div markdown>

**With abstention: outputs that can refuse a cell**

| task | SCARLET unique | Seurat GTP | SingleR pruned | CellTypist prob match | SVMrejection | scmap-cluster | scmap-cell | scType |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| A/10x_5cl | 100.0 (1.7) | 100.0 (0.0) | 100.0 (1.6) | 100.0 (0.0) | 100.0 (0.0) | 100.0 (0.0) | 100.0 (0.0) | 100.0 (0.0) |
| A/amb | 100.0 (0.5) | 100.0 (0.0) | 100.0 (1.4) | 99.9 (0.0) | 100.0 (0.0) | 100.0 (0.2) | 100.0 (0.1) | 99.7 (0.0) |
| A/amb16 | 99.0 (1.8) | 99.7 (2.0) | 98.7 (3.2) | 97.5 (0.1) | 99.8 (1.1) | 98.2 (1.3) | 99.5 (4.6) | 93.5 (0.0) |
| A/celseq2_5cl | 100.0 (2.8) | 100.0 (0.0) | 100.0 (1.6) | 100.0 (0.0) | 100.0 (0.2) | 100.0 (0.2) | 100.0 (0.0) | 99.7 (0.0) |
| A/pancreas_baron_human | 97.3 (3.9) | 98.8 (1.7) | 98.0 (4.1) | 96.7 (1.1) | 99.1 (1.6) | 94.9 (7.8) | 98.1 (4.2) | 71.8 (0.0) |
| A/pancreas_baron_mouse | 100.0 (3.5) | 99.6 (4.6) | 98.3 (3.1) | 99.3 (4.4) | 99.4 (2.3) | 96.6 (15.0) | 98.2 (6.4) | 74.3 (0.0) |
| A/pancreas_muraro | 97.0 (2.2) | 96.6 (0.9) | 97.8 (5.7) | 96.8 (1.0) | 97.5 (1.2) | 96.6 (1.1) | 97.4 (4.1) | 86.8 (0.0) |
| A/pancreas_segerstolpe | 99.0 (3.3) | 99.7 (1.2) | 99.7 (4.5) | 99.2 (1.1) | 100.0 (1.8) | 99.4 (3.3) | 99.9 (5.2) | 75.6 (0.0) |
| A/pancreas_xin | 100.0 (5.6) | 98.7 (2.1) | 99.9 (1.9) | 95.6 (3.1) | 99.4 (0.3) | 99.9 (4.1) | 75.9 (9.6) | 63.9 (0.0) |
| A/tm | 95.7 (3.2) | 97.6 (7.6) | 94.7 (3.1) | 94.7 (2.8) | 98.9 (2.9) | 87.2 (1.8) | 97.6 (9.6) | 8.8 (0.0) |
| A/zheng_68k | 84.9 (69.3) | 85.3 (55.4) | 49.6 (0.9) | 66.7 (16.0) | 91.9 (61.7) | 43.8 (20.4) | 66.7 (70.0) | 54.0 (0.0) |
| A/zheng_sorted | 84.2 (48.1) | 91.4 (39.3) | 73.1 (0.7) | 93.9 (13.1) | 99.2 (23.5) | 73.1 (7.2) | 72.8 (58.1) | 70.2 (0.0) |

</div>

</div>

### Fig. 1: medians over the 12 tasks

<div class="views" markdown>

<div markdown>

**Without abstention**

| output | tasks | median F1 (% unlabelled) |
| :--- | :--- | :--- |
| SCARLET best | 12 | 98.9 (1.8) |
| Seurat predicted.id | 12 | 97.7 (0.0) |
| SingleR labels | 12 | 97.7 (0.0) |
| CellTypist best match | 12 | 98.3 (0.0) |
| SVM | 12 | 98.5 (0.0) |
| Harmony kNN (transductive) | 12 | 98.9 (0.0) |

</div>

<div markdown>

**With abstention**

| output | tasks | median F1 (% unlabelled) |
| :--- | :--- | :--- |
| SCARLET unique | 12 | 99.0 (3.2) |
| Seurat GTP | 12 | 99.2 (1.9) |
| SingleR pruned | 12 | 98.5 (2.5) |
| CellTypist prob match | 12 | 97.1 (1.1) |
| SVMrejection | 12 | 99.4 (1.4) |
| scmap-cluster | 12 | 97.4 (2.6) |
| scmap-cell | 12 | 98.2 (4.9) |
| scType | 12 | 75.0 (0.0) |

</div>

</div>

### Fig. 1: better / equal / worse per task

<div class="views" markdown>

<div markdown>

**Without abstention**

| SCARLET output | Seurat predicted.id | SingleR labels | CellTypist best match | SVM | Harmony kNN (transductive) |
| :--- | :--- | :--- | :--- | :--- | :--- |
| SCARLET best | 6 / 1 / 5 (1.8 vs 0.0) | 7 / 1 / 4 (1.8 vs 0.0) | 5 / 1 / 6 (1.8 vs 0.0) | 4 / 1 / 7 (1.8 vs 0.0) | 5 / 1 / 6 (1.8 vs 0.0) |

</div>

<div markdown>

**With abstention**

| SCARLET output | Seurat GTP | SingleR pruned | CellTypist prob match | SVMrejection | scmap-cluster | scmap-cell | scType |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| SCARLET unique | 4 / 2 / 6 (3.2 vs 1.9) | 6 / 3 / 3 (3.2 vs 2.5) | 8 / 2 / 2 (3.2 vs 1.1) | 3 / 2 / 7 (3.2 vs 1.4) | 9 / 2 / 1 (3.2 vs 2.6) | 5 / 2 / 5 (3.2 vs 4.9) | 11 / 1 / 0 (3.2 vs 0.0) |

</div>

</div>
