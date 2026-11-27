# Benchmarks: SCARLET's other outputs

The figure pages show one SCARLET: its best and unique calls (for `transfer`, the final calls, after refinement), at the defaults the benchmark used. SCARLET writes more outputs; this page gives their values on the same experiments, from the same runs:

- **forced**: the top cell type for every cell, never refused (`scarlet_top1`; on transfer, the first pass's top label, the cell type ranked first most often across the draws);
- **first pass**: `transfer`'s calls before refinement (`first_pass.tsv`);
- **quick mode**: `transfer --draws 1`, one draw, 10 separate runs with seeds 1 to 10.

[Back to the benchmarks](index.md). How to read the values: [the tables](methods.md#how-to-read-the-tables).

## Fig. 1: annotation within a dataset

### Fig. 1: annotation within a dataset: forced and first pass

The rows of SCARLET's other outputs, from the tables of the figure page (same columns).

| output | tasks | median F1 (% unlabelled) |
| :--- | :--- | :--- |
| SCARLET forced | 12 | 98.4 (0.0) |

| SCARLET output | Seurat predicted.id | SingleR labels | CellTypist best match | SVM | Harmony kNN (transductive) |
| :--- | :--- | :--- | :--- | :--- | :--- |
| SCARLET forced | 3 / 1 / 8 (0.0 vs 0.0) | 6 / 1 / 5 (0.0 vs 0.0) | 3 / 1 / 8 (0.0 vs 0.0) | 3 / 1 / 8 (0.0 vs 0.0) | 2 / 2 / 8 (0.0 vs 0.0) |

## Fig. 3a: PBMC, between protocols and between samples

### Fig. 3a, left of the red line: protocols, sample pbmc1: forced and first pass

The rows of SCARLET's other outputs, from the tables of the figure page (same columns).

| output | experiments | median F1 (% unlabelled) |
| :--- | :--- | :--- |
| SCARLET forced (first pass top label) | 42 | 82.3 (0.0) |
| SCARLET first pass best | 42 | 92.2 (23.5) |

| output | experiments | median F1 (% unlabelled) |
| :--- | :--- | :--- |
| SCARLET first pass unique | 42 | 97.1 (57.3) |

| SCARLET output | Seurat predicted.id | SingleR labels | CellTypist best match | SVM | Harmony kNN (transductive) |
| :--- | :--- | :--- | :--- | :--- | :--- |
| SCARLET forced (first pass top label) | 9 / 0 / 33 (0.0 vs 0.0) | 18 / 0 / 24 (0.0 vs 0.0) | 26 / 0 / 16 (0.0 vs 0.0) | 23 / 0 / 19 (0.0 vs 0.0) | 24 / 1 / 17 (0.0 vs 0.0) |
| SCARLET first pass best | 29 / 2 / 11 (23.5 vs 0.0) | 35 / 1 / 6 (23.5 vs 0.0) | 37 / 0 / 5 (23.5 vs 0.0) | 38 / 0 / 4 (23.5 vs 0.0) | 35 / 0 / 7 (23.5 vs 0.0) |

| SCARLET output | Seurat GTP | SingleR pruned | CellTypist prob match | SVMrejection | scmap-cluster | scmap-cell | scType |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| SCARLET first pass unique | 27 / 0 / 15 (57.3 vs 22.7) | 35 / 1 / 6 (57.3 vs 0.6) | 37 / 1 / 4 (57.3 vs 9.3) | 32 / 0 / 10 (57.3 vs 53.7) | 39 / 2 / 1 (57.3 vs 69.2) | 36 / 1 / 5 (57.3 vs 43.0) | 35 / 0 / 7 (57.3 vs 1.0) |

### Fig. 3a, left of the red line: quick mode

SCARLET with `--draws 1`, 10 separate runs (seeds 1 to 10). For each seed, the median over the 42 experiments; shown: the median of the 10 seeds' values, then the lowest and highest seed. Accurate (the default) beside it.

| SCARLET output | accurate: median F1 (% unlabelled) | quick: median F1 (% unlabelled), median of 10 seeds | quick: median F1, lowest–highest seed |
| :--- | :--- | :--- | :--- |
| SCARLET forced (first pass top label) | 82.3 (0.0) | 78.6 (0.0) | 77.8–79.6 |
| SCARLET first pass best | 92.2 (23.5) | 83.1 (7.9) | 82.4–83.6 |
| SCARLET refined best | 94.0 (33.4) | 89.8 (8.2) | 88.0–90.4 |
| SCARLET first pass unique | 97.1 (57.3) | 92.9 (40.3) | 91.8–93.7 |
| SCARLET refined unique | 96.1 (50.7) | 93.0 (23.2) | 92.1–93.8 |

### Fig. 3a, right of the red line: pbmc1 → pbmc2, same protocol: forced and first pass

The rows of SCARLET's other outputs, from the tables of the figure page (same columns).

| output | experiments | median F1 (% unlabelled) |
| :--- | :--- | :--- |
| SCARLET forced (first pass top label) | 6 | 84.5 (0.0) |
| SCARLET first pass best | 6 | 92.5 (10.4) |

| output | experiments | median F1 (% unlabelled) |
| :--- | :--- | :--- |
| SCARLET first pass unique | 6 | 96.6 (30.6) |

| SCARLET output | Seurat predicted.id | SingleR labels | CellTypist best match | SVM | Harmony kNN (transductive) |
| :--- | :--- | :--- | :--- | :--- | :--- |
| SCARLET forced (first pass top label) | 1 / 0 / 5 (0.0 vs 0.0) | 5 / 0 / 1 (0.0 vs 0.0) | 2 / 0 / 4 (0.0 vs 0.0) | 2 / 0 / 4 (0.0 vs 0.0) | 2 / 0 / 4 (0.0 vs 0.0) |
| SCARLET first pass best | 4 / 0 / 2 (10.4 vs 0.0) | 5 / 0 / 1 (10.4 vs 0.0) | 4 / 0 / 2 (10.4 vs 0.0) | 3 / 0 / 3 (10.4 vs 0.0) | 3 / 0 / 3 (10.4 vs 0.0) |

| SCARLET output | Seurat GTP | SingleR pruned | CellTypist prob match | SVMrejection | scmap-cluster | scmap-cell | scType |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| SCARLET first pass unique | 4 / 0 / 2 (30.6 vs 16.7) | 5 / 0 / 1 (30.6 vs 0.2) | 4 / 0 / 2 (30.6 vs 9.1) | 4 / 0 / 2 (30.6 vs 20.7) | 5 / 0 / 1 (30.6 vs 48.9) | 4 / 0 / 2 (30.6 vs 37.6) | 6 / 0 / 0 (30.6 vs 0.3) |

### Fig. 3a, right of the red line: quick mode

SCARLET with `--draws 1`, 10 separate runs (seeds 1 to 10). For each seed, the median over the 6 experiments; shown: the median of the 10 seeds' values, then the lowest and highest seed. Accurate (the default) beside it.

| SCARLET output | accurate: median F1 (% unlabelled) | quick: median F1 (% unlabelled), median of 10 seeds | quick: median F1, lowest–highest seed |
| :--- | :--- | :--- | :--- |
| SCARLET forced (first pass top label) | 84.5 (0.0) | 84.8 (0.0) | 83.1–86.5 |
| SCARLET first pass best | 92.5 (10.4) | 86.2 (3.2) | 85.7–86.8 |
| SCARLET refined best | 94.3 (11.2) | 91.3 (4.1) | 90.9–92.0 |
| SCARLET first pass unique | 96.6 (30.6) | 93.8 (24.2) | 93.1–94.8 |
| SCARLET refined unique | 95.4 (20.9) | 92.8 (11.0) | 92.5–93.6 |

## Fig. 4a and 4b: mouse brain, 3 and 34 cell populations

### Fig. 4a: brain, 3 cell populations: forced and first pass

The rows of SCARLET's other outputs, from the tables of the figure page (same columns).

| output | experiments | median F1 (% unlabelled) |
| :--- | :--- | :--- |
| SCARLET forced (first pass top label) | 2 | 100.0 (0.0) |
| SCARLET first pass best | 2 | 100.0 (1.1) |

| output | experiments | median F1 (% unlabelled) |
| :--- | :--- | :--- |
| SCARLET first pass unique | 2 | 100.0 (2.0) |

| SCARLET output | Seurat predicted.id | SingleR labels | CellTypist best match | SVM | Harmony kNN (transductive) |
| :--- | :--- | :--- | :--- | :--- | :--- |
| SCARLET forced (first pass top label) | 0 / 1 / 1 (0.0 vs 0.0) | 0 / 1 / 1 (0.0 vs 0.0) | 1 / 0 / 1 (0.0 vs 0.0) | 0 / 1 / 1 (0.0 vs 0.0) | 0 / 1 / 1 (0.0 vs 0.0) |
| SCARLET first pass best | 0 / 2 / 0 (1.1 vs 0.0) | 0 / 2 / 0 (1.1 vs 0.0) | 2 / 0 / 0 (1.1 vs 0.0) | 0 / 2 / 0 (1.1 vs 0.0) | 0 / 2 / 0 (1.1 vs 0.0) |

| SCARLET output | Seurat GTP | SingleR pruned | CellTypist prob match | SVMrejection | scmap-cluster | scmap-cell | scType |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| SCARLET first pass unique | 0 / 2 / 0 (2.0 vs 0.0) | 0 / 2 / 0 (2.0 vs 3.1) | 2 / 0 / 0 (2.0 vs 0.0) | 0 / 2 / 0 (2.0 vs 0.0) | 0 / 2 / 0 (2.0 vs 0.3) | 2 / 0 / 0 (2.0 vs 0.1) | 2 / 0 / 0 (2.0 vs 4.0) |

### Fig. 4a: quick mode

SCARLET with `--draws 1`, 10 separate runs (seeds 1 to 10). For each seed, the median over the 2 experiments; shown: the median of the 10 seeds' values, then the lowest and highest seed. Accurate (the default) beside it.

| SCARLET output | accurate: median F1 (% unlabelled) | quick: median F1 (% unlabelled), median of 10 seeds | quick: median F1, lowest–highest seed |
| :--- | :--- | :--- | :--- |
| SCARLET forced (first pass top label) | 100.0 (0.0) | 100.0 (0.0) | 100.0–100.0 |
| SCARLET first pass best | 100.0 (1.1) | 100.0 (0.9) | 100.0–100.0 |
| SCARLET refined best | 100.0 (1.3) | 100.0 (0.6) | 100.0–100.0 |
| SCARLET first pass unique | 100.0 (2.0) | 100.0 (1.7) | 100.0–100.0 |
| SCARLET refined unique | 100.0 (1.8) | 100.0 (0.9) | 100.0–100.0 |

### Fig. 4b: brain, 34 cell populations: forced and first pass

The rows of SCARLET's other outputs, from the tables of the figure page (same columns).

| output | experiments | median F1 (% unlabelled) |
| :--- | :--- | :--- |
| SCARLET forced (first pass top label) | 2 | 91.3 (0.0) |
| SCARLET first pass best | 2 | 92.6 (9.4) |

| output | experiments | median F1 (% unlabelled) |
| :--- | :--- | :--- |
| SCARLET first pass unique | 2 | 94.4 (14.3) |

| SCARLET output | Seurat predicted.id | SingleR labels | CellTypist best match | SVM | Harmony kNN (transductive) |
| :--- | :--- | :--- | :--- | :--- | :--- |
| SCARLET forced (first pass top label) | 0 / 0 / 2 (0.0 vs 0.0) | 0 / 0 / 2 (0.0 vs 0.0) | 0 / 0 / 2 (0.0 vs 0.0) | 0 / 0 / 2 (0.0 vs 0.0) | 0 / 0 / 2 (0.0 vs 0.0) |
| SCARLET first pass best | 0 / 0 / 2 (9.4 vs 0.0) | 1 / 0 / 1 (9.4 vs 0.0) | 1 / 0 / 1 (9.4 vs 0.0) | 0 / 0 / 2 (9.4 vs 0.0) | 0 / 0 / 2 (9.4 vs 0.0) |

| SCARLET output | Seurat GTP | SingleR pruned | CellTypist prob match | SVMrejection | scmap-cluster | scmap-cell | scType |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| SCARLET first pass unique | 0 / 0 / 2 (14.3 vs 15.5) | 2 / 0 / 0 (14.3 vs 1.3) | 2 / 0 / 0 (14.3 vs 13.4) | — | 2 / 0 / 0 (14.3 vs 3.6) | 0 / 0 / 2 (14.3 vs 29.5) | 2 / 0 / 0 (14.3 vs 0.0) |

### Fig. 4b: quick mode

SCARLET with `--draws 1`, 10 separate runs (seeds 1 to 10). For each seed, the median over the 2 experiments; shown: the median of the 10 seeds' values, then the lowest and highest seed. Accurate (the default) beside it.

| SCARLET output | accurate: median F1 (% unlabelled) | quick: median F1 (% unlabelled), median of 10 seeds | quick: median F1, lowest–highest seed |
| :--- | :--- | :--- | :--- |
| SCARLET forced (first pass top label) | 91.3 (0.0) | 91.2 (0.0) | 91.1–91.6 |
| SCARLET first pass best | 92.6 (9.4) | 91.6 (7.2) | 91.4–91.9 |
| SCARLET refined best | 91.9 (1.7) | 91.8 (0.7) | 91.3–92.2 |
| SCARLET first pass unique | 94.4 (14.3) | 93.7 (12.0) | 93.4–93.8 |
| SCARLET refined unique | 92.6 (3.2) | 92.3 (1.7) | 92.0–92.5 |

## Fig. 5a: pancreas, four datasets

### Fig. 5a: pancreas, unaligned: forced and first pass

The rows of SCARLET's other outputs, from the tables of the figure page (same columns).

| output | experiments | median F1 (% unlabelled) |
| :--- | :--- | :--- |
| SCARLET forced (first pass top label) | 4 | 97.6 (0.0) |
| SCARLET first pass best | 4 | 99.6 (7.1) |

| output | experiments | median F1 (% unlabelled) |
| :--- | :--- | :--- |
| SCARLET first pass unique | 4 | 99.9 (18.2) |

| SCARLET output | Seurat predicted.id | SingleR labels | CellTypist best match | SVM | Harmony kNN (transductive) |
| :--- | :--- | :--- | :--- | :--- | :--- |
| SCARLET forced (first pass top label) | 1 / 0 / 3 (0.0 vs 0.0) | 2 / 0 / 2 (0.0 vs 0.0) | 0 / 0 / 4 (0.0 vs 0.0) | 0 / 0 / 4 (0.0 vs 0.0) | 0 / 0 / 4 (0.0 vs 0.0) |
| SCARLET first pass best | 4 / 0 / 0 (7.1 vs 0.0) | 4 / 0 / 0 (7.1 vs 0.0) | 4 / 0 / 0 (7.1 vs 0.0) | 4 / 0 / 0 (7.1 vs 0.0) | 4 / 0 / 0 (7.1 vs 0.0) |

| SCARLET output | Seurat GTP | SingleR pruned | CellTypist prob match | SVMrejection | scmap-cluster | scmap-cell | scType |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| SCARLET first pass unique | 4 / 0 / 0 (18.2 vs 4.1) | 3 / 1 / 0 (18.2 vs 1.8) | 4 / 0 / 0 (18.2 vs 0.3) | 3 / 0 / 1 (18.2 vs 1.2) | 3 / 0 / 1 (18.2 vs 37.6) | 4 / 0 / 0 (18.2 vs 11.8) | 4 / 0 / 0 (18.2 vs 0.9) |

### Fig. 5a: quick mode

SCARLET with `--draws 1`, 10 separate runs (seeds 1 to 10). For each seed, the median over the 4 experiments; shown: the median of the 10 seeds' values, then the lowest and highest seed. Accurate (the default) beside it.

| SCARLET output | accurate: median F1 (% unlabelled) | quick: median F1 (% unlabelled), median of 10 seeds | quick: median F1, lowest–highest seed |
| :--- | :--- | :--- | :--- |
| SCARLET forced (first pass top label) | 97.6 (0.0) | 97.1 (0.0) | 96.6–97.3 |
| SCARLET first pass best | 99.6 (7.1) | 99.0 (5.3) | 98.7–99.4 |
| SCARLET refined best | 99.6 (11.9) | 99.1 (5.9) | 99.0–99.4 |
| SCARLET first pass unique | 99.9 (18.2) | 99.8 (14.3) | 99.6–99.9 |
| SCARLET refined unique | 99.8 (18.8) | 99.6 (8.9) | 99.4–99.7 |

## Supp. Fig. S8: CellBench

### Supp. Fig. S8: CellBench: forced and first pass

The rows of SCARLET's other outputs, from the tables of the figure page (same columns).

| output | experiments | median F1 (% unlabelled) |
| :--- | :--- | :--- |
| SCARLET forced (first pass top label) | 2 | 99.6 (0.0) |
| SCARLET first pass best | 2 | 100.0 (5.2) |

| output | experiments | median F1 (% unlabelled) |
| :--- | :--- | :--- |
| SCARLET first pass unique | 2 | 100.0 (7.5) |

| SCARLET output | Seurat predicted.id | SingleR labels | CellTypist best match | SVM | Harmony kNN (transductive) |
| :--- | :--- | :--- | :--- | :--- | :--- |
| SCARLET forced (first pass top label) | 0 / 1 / 1 (0.0 vs 0.0) | 1 / 0 / 1 (0.0 vs 0.0) | 0 / 1 / 1 (0.0 vs 0.0) | 0 / 1 / 1 (0.0 vs 0.0) | 0 / 1 / 1 (0.0 vs 0.0) |
| SCARLET first pass best | 0 / 1 / 1 (5.2 vs 0.0) | 1 / 0 / 1 (5.2 vs 0.0) | 0 / 1 / 1 (5.2 vs 0.0) | 0 / 1 / 1 (5.2 vs 0.0) | 0 / 1 / 1 (5.2 vs 0.0) |

| SCARLET output | Seurat GTP | SingleR pruned | CellTypist prob match | SVMrejection | scmap-cluster | scmap-cell | scType |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| SCARLET first pass unique | 0 / 2 / 0 (7.5 vs 0.0) | 0 / 2 / 0 (7.5 vs 2.0) | 2 / 0 / 0 (7.5 vs 0.1) | 0 / 2 / 0 (7.5 vs 0.4) | 0 / 2 / 0 (7.5 vs 2.6) | 0 / 2 / 0 (7.5 vs 0.0) | 0 / 2 / 0 (7.5 vs 0.0) |

### Supp. Fig. S8: quick mode

SCARLET with `--draws 1`, 10 separate runs (seeds 1 to 10). For each seed, the median over the 2 experiments; shown: the median of the 10 seeds' values, then the lowest and highest seed. Accurate (the default) beside it.

| SCARLET output | accurate: median F1 (% unlabelled) | quick: median F1 (% unlabelled), median of 10 seeds | quick: median F1, lowest–highest seed |
| :--- | :--- | :--- | :--- |
| SCARLET forced (first pass top label) | 99.6 (0.0) | 99.6 (0.0) | 99.5–99.7 |
| SCARLET first pass best | 100.0 (5.2) | 100.0 (4.4) | 99.9–100.0 |
| SCARLET refined best | 100.0 (1.6) | 100.0 (1.0) | 100.0–100.0 |
| SCARLET first pass unique | 100.0 (7.5) | 100.0 (6.4) | 100.0–100.0 |
| SCARLET refined unique | 100.0 (2.8) | 100.0 (1.7) | 100.0–100.0 |

## Fig. 6: rejection, and our extension

### Fig. 6: rejection: forced and first pass

The rows of SCARLET's other outputs, from the tables of the figure page (same columns).

| output | a | b | c | d |
| :--- | :--- | :--- | :--- | :--- |
| SCARLET first pass unique, accurate | 98.5 | 94.5 | 98.0 | 100.0 |
| SCARLET first pass best, accurate | 93.1 | 89.3 | 90.4 | 99.8 |

**Fig. 6, Part 1: quick mode**

Quick mode (`--draws 1`, seeds 1 to 10), % refused: the median of the 10 seeds (lowest–highest seed).

| SCARLET output | a | b | c | d |
| :--- | :--- | :--- | :--- | :--- |
| SCARLET first pass unique | 91.0 (88.8–93.7) | 92.9 (90.8–94.2) | 95.6 (93.6–98.1) | 95.2 (91.2–96.7) |
| SCARLET refined unique | 88.6 (85.6–93.3) | 92.3 (90.6–92.9) | 93.7 (88.5–96.9) | 94.8 (91.1–96.7) |
| SCARLET first pass best | 34.3 (30.2–40.7) | 82.5 (77.7–85.3) | 79.0 (74.9–86.0) | 80.6 (73.5–84.2) |
| SCARLET refined best | 66.4 (49.7–71.3) | 67.8 (55.3–71.6) | 71.3 (52.7–78.3) | 42.5 (35.0–47.7) |
