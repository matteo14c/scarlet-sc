# Fig. 6: rejection, and our extension

[Back to the benchmarks](index.md) · [how to read the tables](methods.md#how-to-read-the-tables) · [SCARLET's other outputs](further.md)

## Fig. 6: rejection

The benchmark's Fig. 6 (Parts 1 and 2), and a two-sided extension of ours (Part 3, not in Fig. 6). Code: SCARLET 0.6.3 (see [Code versions](methods.md#code-versions));
every other tool run as in [How each tool was run](methods.md#how-each-tool-was-run).

- **Each tool's documented refusal output.** A cell counts as refused only when that output says it is not
  assigned. CellTypist multi-labels are calls.
- **Negative controls are lower bounds.** The reference and the query share immune, endothelial and stromal cells.
  So some query cells in the "negative" datasets have a true match among the reference's cell types, and accepting
  them is not an error. A perfect tool would therefore refuse fewer than 100 % of these cells.
- **SCARLET.** Part 1 is a transfer, so SCARLET's first pass and refined outputs are both shown, accurate and quick.
  Part 2 and our Part 3 are `learn` + `assign` within one dataset (Zheng 68K, downsampled as in the benchmark), where SCARLET
  has no first pass, refinement or draws.
- In transfer, the refinement can assign cells that the first pass refused (compare the first pass unique and
  refined unique rows in Part 1).

### Fig. 6, Part 1: negative controls, % refused

a Zheng 68K PBMC → Baron Human; b Baron Human → Zheng 68K; c AMB16 → Baron Mouse; d Baron Mouse → AMB16. The
outputs that never refuse are listed below the refusal outputs for completeness.

| output | a | b | c | d |
| :--- | :--- | :--- | :--- | :--- |
| SCARLET unique | 99.4 | 93.6 | 99.2 | 100.0 |
| Seurat GTP | 91.9 | 94.6 | 66.0 | 0.1 |
| CellTypist prob match | 0.0 | 68.4 | 0.0 | 0.0 |
| SingleR pruned | 0.4 | 0.1 | 0.9 | 0.2 |
| scmap-cluster | 99.0 | 100.0 | 100.0 | 100.0 |
| scmap-cell | 93.3 | 95.7 | 99.9 | 26.7 |
| SVMrejection | 99.9 | 98.5 | 100.0 | 23.4 |
| scType | 27.5 | 80.9 | 45.1 | 42.6 |
| SCARLET best | 97.8 | 80.3 | 96.3 | 99.8 |
| Harmony kNN (transductive) | 0.0 | 0.0 | 0.0 | 0.0 |
| Seurat predicted.id | 0.0 | 0.0 | 0.0 | 0.0 |
| CellTypist best match | 0.0 | 0.0 | 0.0 | 0.0 |
| SingleR labels | 0.0 | 0.0 | 0.0 | 0.0 |
| SVM | 0.0 | 0.0 | 0.0 | 0.0 |

### Fig. 6, Part 2: unseen populations, % refused, and the labels that absorb the accepted unseen cells

Held out of the training set: i all T cells, ii CD4+ T cells, iii CD4+/CD45RO+ memory T cells; test = the held-out
cells. The absorbing labels are the labels (up to three) most often given to the unseen cells that were not
refused, each with its share of the unseen cells.

| output | i refused | i absorbing labels (% of unseen cells) | ii refused | ii absorbing labels (% of unseen cells) | iii refused | iii absorbing labels (% of unseen cells) |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| SCARLET unique | 54.1 | CD56+ NK 30.2; CD19+ B 15.6; CD34+ 0.1 | 58.2 | CD8+/CD45RA+ Naive Cytotoxic 22.7; CD8+ Cytotoxic T 18.9; CD19+ B 0.2 | 88.7 | CD4+/CD25 T Reg 9.5; CD8+ Cytotoxic T 1.5; CD4+/CD45RA+/CD25- Naive T 0.3 |
| Seurat GTP | 6.0 | CD19+ B 74.0; CD56+ NK 19.9; CD34+ 0.1 | 84.5 | CD8+ Cytotoxic T 11.2; CD8+/CD45RA+ Naive Cytotoxic 4.2; CD19+ B 0.1 | 84.3 | CD4+/CD25 T Reg 8.6; CD8+ Cytotoxic T 7.1 |
| CellTypist prob match | 7.1 | CD19+ B 73.3; CD56+ NK 18.7; CD19+ B\ | CD56+ NK 0.7 | 18.0 | CD8+ Cytotoxic T 44.9; CD8+/CD45RA+ Naive Cytotoxic 30.9; CD8+ Cytotoxic T\ | CD8+/CD45RA+ Naive Cytotoxic 5.8 | 36.5 | CD8+ Cytotoxic T 33.5; CD4+/CD25 T Reg 16.3; CD4+/CD25 T Reg\ | CD8+ Cytotoxic T 5.9 |
| SingleR pruned | 2.1 | CD19+ B 87.9; CD56+ NK 9.5; CD34+ 0.5 | 0.8 | CD8+/CD45RA+ Naive Cytotoxic 91.4; CD8+ Cytotoxic T 7.2; CD34+ 0.2 | 0.9 | CD4+/CD25 T Reg 70.9; CD4+/CD45RA+/CD25- Naive T 11.6; CD4+ T Helper2 7.4 |
| scmap-cluster | 62.3 | CD19+ B 27.9; CD56+ NK 9.7; Dendritic 0.2 | 9.9 | CD8+ Cytotoxic T 45.4; CD8+/CD45RA+ Naive Cytotoxic 44.5; Dendritic 0.2 | 9.2 | CD4+/CD25 T Reg 68.2; CD8+ Cytotoxic T 11.0; CD4+/CD45RA+/CD25- Naive T 6.8 |
| scmap-cell | 44.6 | CD19+ B 40.4; CD56+ NK 14.9; Dendritic 0.1 | 71.7 | CD8+ Cytotoxic T 18.8; CD8+/CD45RA+ Naive Cytotoxic 8.8; CD19+ B 0.4 | 86.1 | CD8+ Cytotoxic T 7.4; CD4+/CD25 T Reg 4.5; CD8+/CD45RA+ Naive Cytotoxic 2.1 |
| SVMrejection | 8.2 | CD19+ B 75.9; CD56+ NK 15.7; Dendritic 0.2 | 75.6 | CD8+ Cytotoxic T 14.9; CD8+/CD45RA+ Naive Cytotoxic 9.2; Dendritic 0.2 | 89.6 | CD8+ Cytotoxic T 9.2; CD4+/CD25 T Reg 1.2 |
| scType | 70.3 | CD56+ NK 23.7; Dendritic 6.0 | 24.5 | CD8+/CD45RA+ Naive Cytotoxic 40.5; CD8+ Cytotoxic T 22.3; CD19+ B 12.7 | 0.0 | CD4+/CD25 T Reg 35.0; CD8+ Cytotoxic T 33.8; CD8+/CD45RA+ Naive Cytotoxic 31.2 |
| SCARLET best | 2.7 | CD56+ NK 53.8; CD19+ B 43.4; CD34+ 0.1 | 4.5 | CD8+/CD45RA+ Naive Cytotoxic 52.7; CD8+ Cytotoxic T 41.9; CD19+ B 0.8 | 14.2 | CD4+/CD25 T Reg 49.9; CD8+ Cytotoxic T 18.4; CD8+/CD45RA+ Naive Cytotoxic 9.5 |
| Harmony kNN (transductive) | 0.0 | CD19+ B 75.7; CD56+ NK 23.5; Dendritic 0.4 | 0.0 | CD8+ Cytotoxic T 52.8; CD8+/CD45RA+ Naive Cytotoxic 45.1; CD19+ B 1.3 | 0.0 | CD4+/CD25 T Reg 51.6; CD8+ Cytotoxic T 32.3; CD8+/CD45RA+ Naive Cytotoxic 14.8 |
| Seurat predicted.id | 0.0 | CD19+ B 77.5; CD56+ NK 22.2; Dendritic 0.1 | 0.0 | CD8+ Cytotoxic T 50.5; CD8+/CD45RA+ Naive Cytotoxic 48.7; Dendritic 0.3 | 0.0 | CD8+ Cytotoxic T 41.8; CD4+/CD25 T Reg 38.6; CD8+/CD45RA+ Naive Cytotoxic 19.6 |
| CellTypist best match | 0.0 | CD19+ B 77.9; CD56+ NK 21.9; Dendritic 0.1 | 0.0 | CD8+ Cytotoxic T 58.0; CD8+/CD45RA+ Naive Cytotoxic 40.6; CD19+ B 0.8 | 0.0 | CD8+ Cytotoxic T 50.1; CD4+/CD25 T Reg 32.3; CD8+/CD45RA+ Naive Cytotoxic 14.5 |
| SingleR labels | 0.0 | CD19+ B 90.0; CD56+ NK 9.5; CD34+ 0.5 | 0.0 | CD8+/CD45RA+ Naive Cytotoxic 92.2; CD8+ Cytotoxic T 7.2; CD34+ 0.2 | 0.0 | CD4+/CD25 T Reg 71.8; CD4+/CD45RA+/CD25- Naive T 11.6; CD4+ T Helper2 7.4 |
| SVM | 0.0 | CD19+ B 80.2; CD56+ NK 19.4; Dendritic 0.4 | 0.0 | CD8+ Cytotoxic T 54.7; CD8+/CD45RA+ Naive Cytotoxic 41.9; CD19+ B 2.3 | 0.0 | CD8+ Cytotoxic T 46.6; CD4+/CD25 T Reg 39.5; CD8+/CD45RA+ Naive Cytotoxic 10.4 |

### Our extension, not in Fig. 6 — Part 3, two-sided: unseen refused / seen refused / median F1 on assigned seen cells (% unlabelled)

The same hold-outs; test = the unseen cells plus a stratified sample of seen cells; refusals are not counted as
wrong.

| output | i | ii | iii |
| :--- | :--- | :--- | :--- |
| SCARLET unique | 51.8 / 19.6 / 98.3 (19.6) | 50.3 / 25.2 / 87.6 (25.2) | 81.3 / 63.2 / 89.1 (63.2) |
| Seurat GTP | 6.2 / 11.2 / 98.4 (11.2) | 82.1 / 49.0 / 86.4 (49.0) | 84.3 / 54.3 / 91.5 (54.3) |
| CellTypist prob match | 5.3 / 3.2 / 81.9 (3.2) | 20.2 / 11.8 / 76.7 (11.8) | 40.7 / 24.5 / 71.7 (24.5) |
| SingleR pruned | 0.9 / 1.0 / 84.0 (1.0) | 1.3 / 1.6 / 75.0 (1.6) | 0.0 / 0.4 / 59.3 (0.4) |
| scmap-cluster | 82.2 / 28.0 / 85.7 (28.0) | 10.0 / 8.7 / 67.4 (8.7) | 11.3 / 10.7 / 49.9 (10.7) |
| scmap-cell | 34.7 / 15.9 / 98.7 (15.9) | 69.2 / 58.7 / 88.5 (58.7) | 81.0 / 67.4 / 51.1 (67.4) |
| SVMrejection | 8.6 / 7.2 / 90.6 (7.2) | 74.2 / 51.3 / 87.4 (51.3) | 89.6 / 64.0 / 84.5 (64.0) |
| scType | 72.8 / 12.7 / 61.9 (12.7) | 0.0 / 0.0 / 71.1 (0.0) | 0.0 / 0.0 / 53.1 (0.0) |
| SCARLET best | 2.9 / 5.5 / 93.7 (5.5) | 5.3 / 2.5 / 80.9 (2.5) | 12.5 / 8.4 / 75.2 (8.4) |
| Harmony kNN (transductive) | 0.0 / 0.0 / 81.2 (0.0) | 0.0 / 0.0 / 73.1 (0.0) | 0.0 / 0.0 / 70.3 (0.0) |
| Seurat predicted.id | 0.0 / 0.0 / 98.0 (0.0) | 0.0 / 0.0 / 80.0 (0.0) | 0.0 / 0.0 / 67.8 (0.0) |
| CellTypist best match | 0.0 / 0.0 / 80.8 (0.0) | 0.0 / 0.0 / 79.4 (0.0) | 0.0 / 0.0 / 67.7 (0.0) |
| SingleR labels | 0.0 / 0.0 / 84.0 (0.0) | 0.0 / 0.0 / 75.0 (0.0) | 0.0 / 0.0 / 59.3 (0.0) |
| SVM | 0.0 / 0.0 / 82.0 (0.0) | 0.0 / 0.0 / 76.2 (0.0) | 0.0 / 0.0 / 69.8 (0.0) |

### Agreement with the benchmark's Fig. 6 (% unlabelled, the same experiments)

For the three tools with a refusal output in Fig. 6, our % unlabelled against the published values (read from the
figure). The caveat above applies: the negative controls are lower bounds, for the published runs as for ours. No
SCARLET code is involved in this table.

| experiment | tool | ours | published | difference |
| :--- | :--- | :--- | :--- | :--- |
| negative control: train Zheng 68K -> test Baron Human | scmap-cluster | 99.0 | 99.0 | -0.0 |
| negative control: train Zheng 68K -> test Baron Human | scmap-cell | 93.3 | 92.4 | 0.9 |
| negative control: train Zheng 68K -> test Baron Human | SVMrejection | 99.9 | 99.9 | 0.0 |
| negative control: train Baron Human -> test Zheng 68K | scmap-cluster | 100.0 | 100.0 | 0.0 |
| negative control: train Baron Human -> test Zheng 68K | scmap-cell | 95.7 | 95.9 | -0.2 |
| negative control: train Baron Human -> test Zheng 68K | SVMrejection | 98.5 | 98.5 | -0.0 |
| negative control: train AMB16 -> test Baron Mouse | scmap-cluster | 100.0 | 100.0 | 0.0 |
| negative control: train AMB16 -> test Baron Mouse | scmap-cell | 99.9 | 100.0 | -0.1 |
| negative control: train AMB16 -> test Baron Mouse | SVMrejection | 100.0 | 100.0 | 0.0 |
| negative control: train Baron Mouse -> test AMB16 | scmap-cluster | 100.0 | 100.0 | 0.0 |
| negative control: train Baron Mouse -> test AMB16 | scmap-cell | 26.7 | 18.8 | 7.9 |
| negative control: train Baron Mouse -> test AMB16 | SVMrejection | 23.4 | 26.2 | -2.8 |
| unseen population (Zheng 68K, left out): T cells | scmap-cluster | 62.3 | 62.34 | -0.0 |
| unseen population (Zheng 68K, left out): T cells | scmap-cell | 44.6 | 46.87 | -2.3 |
| unseen population (Zheng 68K, left out): T cells | SVMrejection | 8.2 | 8.17 | 0.0 |
| unseen population (Zheng 68K, left out): CD4+ T cells | scmap-cluster | 9.9 | 9.88 | -0.0 |
| unseen population (Zheng 68K, left out): CD4+ T cells | scmap-cell | 71.7 | 71.09 | 0.7 |
| unseen population (Zheng 68K, left out): CD4+ T cells | SVMrejection | 75.6 | 75.63 | -0.0 |
| unseen population (Zheng 68K, left out): CD4+/CD45RO+ Memory T cells | scmap-cluster | 9.2 | 9.2 | -0.0 |
| unseen population (Zheng 68K, left out): CD4+/CD45RO+ Memory T cells | scmap-cell | 86.1 | 83.68 | 2.4 |
| unseen population (Zheng 68K, left out): CD4+/CD45RO+ Memory T cells | SVMrejection | 89.6 | 89.91 | -0.3 |
