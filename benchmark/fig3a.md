# Fig. 3a: PBMC, between protocols and between samples

[Back to the benchmarks](index.md) · [how to read the tables](methods.md#how-to-read-the-tables) · [SCARLET's other outputs](further.md)

## Fig. 3a, left of the red line: protocols, sample pbmc1

> Classification performance across the PbmcBench datasets. … Results shown to the left of the red line represent the comparison between different protocols using sample pbmc1.
>
> — Abdelaal et al. (2019), the legend of Fig. 3

One row per experiment: the 42 PbmcBench experiments between protocols on sample pbmc1, training protocol → test protocol. Each value: median F1 (% unlabelled). Code: SCARLET 0.6.3, `transfer`, accurate (see [Transfer settings](methods.md#transfer-settings-b-and-abstention-part-1)); the other tools as in [How each tool was run](methods.md#how-each-tool-was-run).

<div class="views" markdown>

<div markdown>

**Without abstention: outputs that label every cell**

| training → test | SCARLET best | Seurat predicted.id | SingleR labels | CellTypist best match | SVM | Harmony kNN (transductive) |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| 10XV2 → SM2 (pbmc1) | 100.0 (51.4) | 97.7 (0.0) | 92.4 (0.0) | 93.7 (0.0) | 91.1 (0.0) | 84.7 (0.0) |
| 10XV2 → 10XV3 (pbmc1) | 95.9 (6.3) | 96.6 (0.0) | 93.4 (0.0) | 96.0 (0.0) | 94.1 (0.0) | 94.2 (0.0) |
| 10XV2 → CL (pbmc1) | 96.3 (51.8) | 84.2 (0.0) | 76.5 (0.0) | 78.0 (0.0) | 81.0 (0.0) | 77.2 (0.0) |
| 10XV2 → DR (pbmc1) | 92.0 (21.3) | 85.7 (0.0) | 84.0 (0.0) | 78.3 (0.0) | 82.1 (0.0) | 83.3 (0.0) |
| 10XV2 → iD (pbmc1) | 91.0 (32.0) | 85.5 (0.0) | 69.0 (0.0) | 72.2 (0.0) | 64.0 (0.0) | 12.1 (0.0) |
| 10XV2 → SW (pbmc1) | 91.4 (31.4) | 86.1 (0.0) | 63.8 (0.0) | 77.3 (0.0) | 69.4 (0.0) | 76.2 (0.0) |
| 10XV3 → SM2 (pbmc1) | 98.5 (21.7) | 98.5 (0.0) | 98.5 (0.0) | 95.2 (0.0) | 93.9 (0.0) | 97.8 (0.0) |
| 10XV3 → 10XV2 (pbmc1) | 92.8 (8.2) | 93.2 (0.0) | 90.1 (0.0) | 83.9 (0.0) | 92.4 (0.0) | 88.4 (0.0) |
| 10XV3 → CL (pbmc1) | 97.6 (73.1) | 81.1 (0.0) | 76.5 (0.0) | 62.9 (0.0) | 85.0 (0.0) | 71.8 (0.0) |
| 10XV3 → DR (pbmc1) | 92.5 (27.7) | 75.8 (0.0) | 85.6 (0.0) | 37.9 (0.0) | 69.1 (0.0) | 68.5 (0.0) |
| 10XV3 → iD (pbmc1) | 97.3 (49.8) | 78.6 (0.0) | 81.4 (0.0) | 9.9 (0.0) | 32.8 (0.0) | 11.9 (0.0) |
| 10XV3 → SW (pbmc1) | 91.7 (36.4) | 68.1 (0.0) | 66.4 (0.0) | 48.4 (0.0) | 64.4 (0.0) | 59.0 (0.0) |
| CEL-Seq2 (CL) → SM2 (pbmc1) | 93.2 (11.5) | 94.1 (0.0) | 89.2 (0.0) | 91.0 (0.0) | 93.4 (0.0) | 94.6 (0.0) |
| CEL-Seq2 (CL) → 10XV2 (pbmc1) | 97.7 (59.6) | 86.3 (0.0) | 87.9 (0.0) | 63.3 (0.0) | 79.3 (0.0) | 80.0 (0.0) |
| CEL-Seq2 (CL) → 10XV3 (pbmc1) | 86.5 (8.6) | 85.3 (0.0) | 93.3 (0.0) | 84.6 (0.0) | 93.8 (0.0) | 85.8 (0.0) |
| CEL-Seq2 (CL) → DR (pbmc1) | 91.5 (73.4) | 39.7 (0.0) | 83.9 (0.0) | 10.7 (0.0) | 75.3 (0.0) | 79.9 (0.0) |
| CEL-Seq2 (CL) → iD (pbmc1) | 77.3 (65.3) | 44.3 (0.0) | 73.6 (0.0) | 4.7 (0.0) | 26.8 (0.0) | 8.0 (0.0) |
| CEL-Seq2 (CL) → SW (pbmc1) | 81.2 (69.2) | 77.3 (0.0) | 65.6 (0.0) | 33.8 (0.0) | 65.2 (0.0) | 56.4 (0.0) |
| Drop-seq (DR) → SM2 (pbmc1) | 93.2 (78.7) | 86.7 (0.0) | 83.9 (0.0) | 91.4 (0.0) | 85.1 (0.0) | 85.1 (0.0) |
| Drop-seq (DR) → 10XV2 (pbmc1) | 86.8 (5.3) | 88.3 (0.0) | 66.9 (0.0) | 86.1 (0.0) | 87.1 (0.0) | 82.9 (0.0) |
| Drop-seq (DR) → 10XV3 (pbmc1) | 90.7 (10.4) | 84.0 (0.0) | 75.7 (0.0) | 67.2 (0.0) | 86.3 (0.0) | 78.1 (0.0) |
| Drop-seq (DR) → CL (pbmc1) | 96.9 (71.5) | 80.4 (0.0) | 75.0 (0.0) | 55.6 (0.0) | 83.9 (0.0) | 78.6 (0.0) |
| Drop-seq (DR) → iD (pbmc1) | 88.4 (38.2) | 76.8 (0.0) | 65.1 (0.0) | 76.8 (0.0) | 74.2 (0.0) | 21.9 (0.0) |
| Drop-seq (DR) → SW (pbmc1) | 74.5 (8.9) | 76.6 (0.0) | 42.5 (0.0) | 78.9 (0.0) | 69.2 (0.0) | 59.6 (0.0) |
| Seq-Well (SW) → SM2 (pbmc1) | 96.3 (26.1) | 97.0 (0.0) | 94.8 (0.0) | 98.5 (0.0) | 91.3 (0.0) | 97.4 (0.0) |
| Seq-Well (SW) → 10XV2 (pbmc1) | 96.2 (8.0) | 97.4 (0.0) | 89.9 (0.0) | 90.2 (0.0) | 84.9 (0.0) | 93.2 (0.0) |
| Seq-Well (SW) → 10XV3 (pbmc1) | 96.8 (16.0) | 93.9 (0.0) | 91.7 (0.0) | 89.3 (0.0) | 90.9 (0.0) | 96.0 (0.0) |
| Seq-Well (SW) → CL (pbmc1) | 92.5 (36.0) | 97.4 (0.0) | 86.3 (0.0) | 87.5 (0.0) | 85.2 (0.0) | 95.0 (0.0) |
| Seq-Well (SW) → DR (pbmc1) | 93.8 (22.5) | 82.1 (0.0) | 88.2 (0.0) | 81.7 (0.0) | 77.0 (0.0) | 88.1 (0.0) |
| Seq-Well (SW) → iD (pbmc1) | 94.2 (21.5) | 91.3 (0.0) | 87.2 (0.0) | 78.9 (0.0) | 86.0 (0.0) | 38.5 (0.0) |
| Smart-seq2 (SM2) → 10XV2 (pbmc1) | 92.7 (9.4) | 92.3 (0.0) | 92.9 (0.0) | 76.7 (0.0) | 82.3 (0.0) | 89.2 (0.0) |
| Smart-seq2 (SM2) → 10XV3 (pbmc1) | 98.0 (3.0) | 94.9 (0.0) | 97.0 (0.0) | 85.0 (0.0) | 90.4 (0.0) | 96.8 (0.0) |
| Smart-seq2 (SM2) → CL (pbmc1) | 95.8 (37.4) | 88.8 (0.0) | 84.9 (0.0) | 44.5 (0.0) | 72.5 (0.0) | 93.4 (0.0) |
| Smart-seq2 (SM2) → DR (pbmc1) | 94.2 (34.8) | 81.4 (0.0) | 85.8 (0.0) | 36.6 (0.0) | 55.7 (0.0) | 78.7 (0.0) |
| Smart-seq2 (SM2) → iD (pbmc1) | 96.1 (59.9) | 68.9 (0.0) | 84.4 (0.0) | 11.5 (0.0) | 46.8 (0.0) | 21.3 (0.0) |
| Smart-seq2 (SM2) → SW (pbmc1) | 98.9 (45.5) | 69.3 (0.0) | 84.1 (0.0) | 60.3 (0.0) | 60.1 (0.0) | 60.9 (0.0) |
| inDrop (iD) → SM2 (pbmc1) | 99.0 (36.0) | 98.6 (0.0) | 99.3 (0.0) | 93.4 (0.0) | 94.0 (0.0) | 86.0 (0.0) |
| inDrop (iD) → 10XV2 (pbmc1) | 99.3 (19.9) | 97.4 (0.0) | 91.5 (0.0) | 90.9 (0.0) | 90.1 (0.0) | 74.4 (0.0) |
| inDrop (iD) → 10XV3 (pbmc1) | 97.1 (20.3) | 97.4 (0.0) | 97.6 (0.0) | 90.4 (0.0) | 97.1 (0.0) | 41.2 (0.0) |
| inDrop (iD) → CL (pbmc1) | 48.6 (56.4) | 68.9 (0.0) | 80.0 (0.0) | 69.2 (0.0) | 84.7 (0.0) | 7.7 (0.0) |
| inDrop (iD) → DR (pbmc1) | 98.2 (39.1) | 88.3 (0.0) | 84.3 (0.0) | 88.6 (0.0) | 89.0 (0.0) | 86.9 (0.0) |
| inDrop (iD) → SW (pbmc1) | 92.4 (39.7) | 74.4 (0.0) | 72.0 (0.0) | 80.1 (0.0) | 77.5 (0.0) | 76.0 (0.0) |

</div>

<div markdown>

**With abstention: outputs that can refuse a cell**

| training → test | SCARLET unique | Seurat GTP | SingleR pruned | CellTypist prob match | SVMrejection | scmap-cluster | scmap-cell | scType |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| 10XV2 → SM2 (pbmc1) | 100.0 (60.9) | 98.4 (9.9) | 92.4 (0.0) | 20.1 (0.0) | 7.7 (84.6) | 75.9 (5.9) | 96.9 (35.6) | 21.4 (0.0) |
| 10XV2 → 10XV3 (pbmc1) | 96.6 (10.2) | 98.3 (9.9) | 93.5 (0.1) | 92.3 (0.4) | 98.5 (19.3) | 86.8 (7.1) | 97.7 (28.5) | 93.8 (0.0) |
| 10XV2 → CL (pbmc1) | 100.0 (90.5) | 88.2 (17.8) | 76.5 (2.8) | 0.0 (0.0) | 91.7 (33.2) | 75.0 (3.2) | 80.0 (45.5) | 13.3 (0.0) |
| 10XV2 → DR (pbmc1) | 93.4 (38.8) | 93.6 (22.1) | 86.4 (1.3) | 85.0 (11.4) | 96.6 (45.1) | 89.8 (8.2) | 92.5 (31.9) | 65.0 (4.6) |
| 10XV2 → iD (pbmc1) | 93.3 (67.9) | 91.7 (17.7) | 69.0 (1.8) | 92.0 (42.1) | 82.4 (72.2) | 83.5 (62.0) | 88.9 (42.4) | 63.3 (42.0) |
| 10XV2 → SW (pbmc1) | 93.3 (46.6) | 90.6 (21.3) | 68.2 (4.6) | 72.6 (18.5) | 96.0 (53.9) | 57.7 (22.3) | 75.0 (35.5) | 67.2 (10.8) |
| 10XV3 → SM2 (pbmc1) | 98.5 (30.0) | 98.4 (7.5) | 98.5 (0.0) | 57.7 (0.0) | 79.0 (36.8) | 88.9 (5.5) | 98.4 (23.3) | 83.0 (0.0) |
| 10XV3 → 10XV2 (pbmc1) | 95.5 (18.3) | 96.9 (13.2) | 90.3 (0.2) | 84.1 (9.4) | 96.5 (28.8) | 90.4 (8.6) | 94.3 (24.2) | 89.7 (2.1) |
| 10XV3 → CL (pbmc1) | 100.0 (86.2) | 88.2 (23.7) | 76.5 (0.0) | 5.7 (0.4) | 94.9 (20.6) | 81.0 (2.8) | 80.0 (34.0) | 33.3 (0.0) |
| 10XV3 → DR (pbmc1) | 95.6 (42.3) | 90.9 (27.4) | 85.8 (0.9) | 32.3 (10.2) | 94.4 (68.4) | 87.2 (9.9) | 90.8 (26.6) | 85.4 (3.8) |
| 10XV3 → iD (pbmc1) | 96.9 (68.8) | 90.2 (49.8) | 82.2 (0.7) | 6.7 (1.8) | 11.8 (80.7) | 90.2 (54.1) | 88.6 (32.0) | 76.0 (42.9) |
| 10XV3 → SW (pbmc1) | 91.7 (49.9) | 90.8 (37.9) | 67.3 (1.9) | 36.5 (8.6) | 88.6 (80.9) | 69.6 (25.0) | 74.6 (29.7) | 87.7 (5.1) |
| CEL-Seq2 (CL) → SM2 (pbmc1) | 97.0 (18.6) | 95.8 (13.8) | 89.2 (0.0) | 78.3 (36.8) | 96.5 (27.3) | 64.0 (15.8) | 91.3 (28.1) | 56.6 (0.0) |
| CEL-Seq2 (CL) → 10XV2 (pbmc1) | 98.3 (66.7) | 90.3 (14.0) | 87.9 (0.0) | 62.5 (79.4) | 25.0 (99.6) | 93.3 (67.5) | 87.3 (45.3) | 94.2 (1.1) |
| CEL-Seq2 (CL) → 10XV3 (pbmc1) | 87.4 (13.5) | 89.8 (13.0) | 93.3 (0.0) | 90.7 (46.8) | 99.2 (59.3) | 54.6 (14.3) | 89.5 (34.2) | 98.0 (0.0) |
| CEL-Seq2 (CL) → DR (pbmc1) | 94.5 (97.9) | 90.3 (26.4) | 83.9 (0.0) | 0.0 (86.4) | 0.0 (99.9) | 75.0 (85.0) | 95.0 (42.8) | 84.8 (3.1) |
| CEL-Seq2 (CL) → iD (pbmc1) | 50.0 (97.4) | 79.5 (52.7) | 73.7 (0.1) | 6.3 (46.4) | 0.0 (100.0) | 10.9 (97.4) | 73.3 (89.5) | 79.1 (41.9) |
| CEL-Seq2 (CL) → SW (pbmc1) | 100.0 (98.4) | 84.1 (33.2) | 65.6 (0.0) | 25.8 (74.2) | 0.0 (100.0) | 74.1 (82.9) | 68.3 (61.0) | 93.2 (4.4) |
| Drop-seq (DR) → SM2 (pbmc1) | 98.0 (83.8) | 96.5 (23.3) | 83.9 (0.0) | 75.2 (0.0) | 37.5 (88.9) | 0.0 (86.2) | 87.7 (48.6) | 18.3 (0.0) |
| Drop-seq (DR) → 10XV2 (pbmc1) | 91.0 (17.1) | 94.5 (18.7) | 69.0 (2.3) | 72.7 (5.4) | 94.2 (27.4) | 61.3 (14.8) | 83.2 (37.6) | 67.7 (0.9) |
| Drop-seq (DR) → 10XV3 (pbmc1) | 92.8 (14.6) | 93.5 (18.9) | 75.7 (0.1) | 69.2 (0.2) | 94.4 (32.9) | 85.6 (9.7) | 69.4 (38.5) | 63.6 (0.0) |
| Drop-seq (DR) → CL (pbmc1) | 0.0 (94.5) | 88.5 (25.3) | 75.8 (2.4) | 0.0 (2.4) | 87.5 (39.5) | 0.0 (81.0) | 25.0 (63.6) | 10.6 (0.0) |
| Drop-seq (DR) → iD (pbmc1) | 94.7 (62.0) | 90.3 (43.6) | 66.9 (4.2) | 93.6 (43.5) | 69.6 (84.5) | 75.0 (93.7) | 84.8 (74.2) | 49.5 (43.8) |
| Drop-seq (DR) → SW (pbmc1) | 75.1 (21.5) | 88.1 (25.2) | 46.8 (8.0) | 71.1 (21.0) | 71.3 (53.5) | 38.5 (77.1) | 60.0 (52.4) | 44.0 (3.5) |
| Seq-Well (SW) → SM2 (pbmc1) | 100.0 (51.4) | 97.0 (14.7) | 94.8 (0.0) | 92.3 (0.0) | 92.3 (86.1) | 0.0 (85.3) | 100.0 (36.3) | 97.0 (0.0) |
| Seq-Well (SW) → 10XV2 (pbmc1) | 96.2 (15.6) | 99.4 (17.7) | 93.9 (0.8) | 86.2 (1.9) | 95.2 (25.4) | 89.1 (11.5) | 92.9 (42.0) | 92.7 (1.1) |
| Seq-Well (SW) → 10XV3 (pbmc1) | 98.1 (22.3) | 97.3 (21.5) | 91.8 (0.1) | 89.6 (0.2) | 95.8 (31.4) | 89.0 (11.1) | 96.0 (41.7) | 96.2 (0.0) |
| Seq-Well (SW) → CL (pbmc1) | 100.0 (77.5) | 97.4 (26.6) | 89.3 (1.4) | 13.3 (0.9) | 90.8 (32.0) | 61.1 (46.8) | 89.5 (51.8) | 72.0 (1.8) |
| Seq-Well (SW) → DR (pbmc1) | 95.5 (35.1) | 89.3 (12.5) | 89.2 (2.0) | 85.7 (9.2) | 85.9 (24.5) | 90.7 (23.6) | 90.3 (43.2) | 90.4 (3.0) |
| Seq-Well (SW) → iD (pbmc1) | 98.2 (37.2) | 94.8 (19.2) | 89.2 (3.1) | 87.1 (22.8) | 92.9 (59.2) | 52.8 (90.5) | 96.2 (48.2) | 80.8 (46.6) |
| Smart-seq2 (SM2) → 10XV2 (pbmc1) | 94.3 (19.4) | 98.9 (21.4) | 93.9 (0.4) | 68.7 (58.7) | 0.0 (100.0) | 48.5 (98.0) | 92.9 (57.8) | 97.7 (1.2) |
| Smart-seq2 (SM2) → 10XV3 (pbmc1) | 98.2 (6.9) | 97.8 (9.7) | 97.0 (0.0) | 85.7 (31.1) | 0.0 (100.0) | 97.4 (81.8) | 94.6 (31.7) | 98.4 (0.1) |
| Smart-seq2 (SM2) → CL (pbmc1) | 95.6 (40.7) | 94.4 (20.2) | 86.3 (0.8) | 26.3 (15.2) | 0.0 (100.0) | 57.6 (76.5) | 91.5 (37.0) | 55.2 (0.4) |
| Smart-seq2 (SM2) → DR (pbmc1) | 96.0 (51.9) | 92.4 (30.8) | 86.6 (1.3) | 38.9 (62.3) | 0.0 (100.0) | 0.0 (99.9) | 78.4 (83.3) | 89.0 (3.5) |
| Smart-seq2 (SM2) → iD (pbmc1) | 99.8 (81.8) | 91.6 (54.5) | 85.2 (1.6) | 46.7 (93.1) | 0.0 (100.0) | 50.0 (98.9) | 49.5 (86.3) | 80.3 (43.5) |
| Smart-seq2 (SM2) → SW (pbmc1) | 99.1 (57.1) | 88.2 (38.4) | 84.1 (1.9) | 67.5 (75.7) | 0.0 (100.0) | 0.0 (98.4) | 89.1 (84.0) | 79.4 (11.4) |
| inDrop (iD) → SM2 (pbmc1) | 50.0 (84.6) | 100.0 (59.7) | 99.3 (0.4) | 76.7 (0.0) | 96.6 (27.7) | 0.0 (100.0) | 0.0 (57.7) | 87.1 (0.0) |
| inDrop (iD) → 10XV2 (pbmc1) | 99.6 (31.0) | 99.5 (30.7) | 92.3 (0.5) | 72.2 (2.1) | 93.3 (14.9) | 0.0 (70.9) | 65.6 (73.7) | 98.7 (0.0) |
| inDrop (iD) → 10XV3 (pbmc1) | 91.5 (26.4) | 99.1 (52.6) | 97.6 (0.1) | 76.2 (0.2) | 98.7 (8.3) | 49.7 (82.4) | 57.4 (72.5) | 98.2 (0.0) |
| inDrop (iD) → CL (pbmc1) | 0.0 (98.8) | 74.8 (75.3) | 80.4 (0.4) | 0.9 (0.4) | 89.1 (9.1) | 0.0 (98.8) | 12.5 (63.4) | 1.9 (0.0) |
| inDrop (iD) → DR (pbmc1) | 98.9 (54.0) | 98.1 (47.1) | 84.3 (1.7) | 74.7 (3.6) | 95.4 (20.5) | 0.0 (97.4) | 52.8 (85.6) | 89.8 (0.0) |
| inDrop (iD) → SW (pbmc1) | 96.9 (52.8) | 90.9 (42.9) | 76.9 (5.3) | 77.1 (12.0) | 79.9 (31.8) | 0.0 (96.1) | 65.8 (71.5) | 85.9 (0.0) |

</div>

</div>

### Fig. 3a, left of the red line: medians over the 42 experiments

<div class="views" markdown>

<div markdown>

**Without abstention**

| output | experiments | median F1 (% unlabelled) |
| :--- | :--- | :--- |
| SCARLET best | 42 | 94.0 (33.4) |
| Seurat predicted.id | 42 | 85.9 (0.0) |
| SingleR labels | 42 | 84.7 (0.0) |
| CellTypist best match | 42 | 78.6 (0.0) |
| SVM | 42 | 84.3 (0.0) |
| Harmony kNN (transductive) | 42 | 79.3 (0.0) |

</div>

<div markdown>

**With abstention**

| output | experiments | median F1 (% unlabelled) |
| :--- | :--- | :--- |
| SCARLET unique | 42 | 96.1 (50.7) |
| Seurat GTP | 42 | 92.9 (22.7) |
| SingleR pruned | 42 | 86.1 (0.6) |
| CellTypist prob match | 42 | 71.6 (9.3) |
| SVMrejection | 42 | 88.9 (53.7) |
| scmap-cluster | 42 | 62.6 (69.2) |
| scmap-cell | 42 | 88.2 (43.0) |
| scType | 42 | 81.9 (1.0) |

</div>

</div>

### Fig. 3a, left of the red line: better / equal / worse

<div class="views" markdown>

<div markdown>

**Without abstention**

| SCARLET output | Seurat predicted.id | SingleR labels | CellTypist best match | SVM | Harmony kNN (transductive) |
| :--- | :--- | :--- | :--- | :--- | :--- |
| SCARLET best | 31 / 1 / 10 (33.4 vs 0.0) | 36 / 1 / 5 (33.4 vs 0.0) | 38 / 0 / 4 (33.4 vs 0.0) | 38 / 0 / 4 (33.4 vs 0.0) | 39 / 0 / 3 (33.4 vs 0.0) |

</div>

<div markdown>

**With abstention**

| SCARLET output | Seurat GTP | SingleR pruned | CellTypist prob match | SVMrejection | scmap-cluster | scmap-cell | scType |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| SCARLET unique | 28 / 0 / 14 (50.7 vs 22.7) | 35 / 1 / 6 (50.7 vs 0.6) | 38 / 1 / 3 (50.7 vs 9.3) | 31 / 0 / 11 (50.7 vs 53.7) | 40 / 2 / 0 (50.7 vs 69.2) | 35 / 1 / 6 (50.7 vs 43.0) | 34 / 0 / 8 (50.7 vs 1.0) |

</div>

</div>

## Fig. 3a, right of the red line: pbmc1 → pbmc2, same protocol

> Classification performance across the PbmcBench datasets. … Results shown to the right of the red line represent the comparison between different samples using the same protocol, with pbmc 1 used for training and pbmc2 used for testing.
>
> — Abdelaal et al. (2019), the legend of Fig. 3

One row per experiment: the 6 PbmcBench experiments from sample pbmc1 to sample pbmc2 with the same protocol (10Xv3 has no pbmc2 sample). Each value: median F1 (% unlabelled). Code: SCARLET 0.6.3, `transfer`, accurate (see [Transfer settings](methods.md#transfer-settings-b-and-abstention-part-1)); the other tools as in [How each tool was run](methods.md#how-each-tool-was-run).

<div class="views" markdown>

<div markdown>

**Without abstention: outputs that label every cell**

| training → test | SCARLET best | Seurat predicted.id | SingleR labels | CellTypist best match | SVM | Harmony kNN (transductive) |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| 10XV2 → 10XV2 (pbmc2) | 94.9 (5.7) | 96.7 (0.0) | 90.7 (0.0) | 94.9 (0.0) | 93.3 (0.0) | 95.5 (0.0) |
| CEL-Seq2 (CL) → CL (pbmc2) | 60.7 (27.8) | 71.0 (0.0) | 53.3 (0.0) | 68.7 (0.0) | 70.7 (0.0) | 70.1 (0.0) |
| Drop-seq (DR) → DR (pbmc2) | 93.6 (8.4) | 90.2 (0.0) | 84.5 (0.0) | 77.8 (0.0) | 85.4 (0.0) | 87.2 (0.0) |
| Seq-Well (SW) → SW (pbmc2) | 95.5 (14.0) | 91.4 (0.0) | 78.1 (0.0) | 85.6 (0.0) | 78.6 (0.0) | 58.9 (0.0) |
| Smart-seq2 (SM2) → SM2 (pbmc2) | 90.0 (5.4) | 77.3 (0.0) | 94.5 (0.0) | 83.9 (0.0) | 84.5 (0.0) | 91.0 (0.0) |
| inDrop (iD) → iD (pbmc2) | 95.8 (29.1) | 86.9 (0.0) | 69.8 (0.0) | 84.7 (0.0) | 80.4 (0.0) | 78.7 (0.0) |

</div>

<div markdown>

**With abstention: outputs that can refuse a cell**

| training → test | SCARLET unique | Seurat GTP | SingleR pruned | CellTypist prob match | SVMrejection | scmap-cluster | scmap-cell | scType |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| 10XV2 → 10XV2 (pbmc2) | 94.9 (12.1) | 98.3 (11.2) | 91.2 (0.1) | 81.0 (1.6) | 96.3 (12.2) | 88.8 (1.6) | 97.1 (17.7) | 82.1 (0.5) |
| CEL-Seq2 (CL) → CL (pbmc2) | 63.2 (38.8) | 73.5 (13.6) | 53.3 (0.0) | 74.0 (54.2) | 75.4 (27.8) | 72.5 (4.4) | 73.2 (29.7) | 42.1 (0.0) |
| Drop-seq (DR) → DR (pbmc2) | 95.8 (15.3) | 94.4 (14.1) | 84.7 (0.4) | 70.8 (3.3) | 90.7 (16.4) | 76.7 (39.8) | 78.0 (29.0) | 73.0 (0.9) |
| Seq-Well (SW) → SW (pbmc2) | 99.0 (26.5) | 95.6 (19.6) | 78.6 (1.6) | 93.8 (22.3) | 94.6 (45.6) | 73.0 (58.1) | 91.1 (45.6) | 84.0 (11.8) |
| Smart-seq2 (SM2) → SM2 (pbmc2) | 90.0 (8.5) | 88.5 (19.3) | 94.5 (0.0) | 91.0 (13.5) | 94.6 (25.1) | 0.0 (82.5) | 79.1 (47.1) | 7.4 (0.0) |
| inDrop (iD) → iD (pbmc2) | 98.2 (34.2) | 97.2 (26.3) | 72.1 (3.7) | 78.6 (4.8) | 86.1 (15.2) | 0.0 (73.3) | 70.0 (60.4) | 72.0 (0.0) |

</div>

</div>

### Fig. 3a, right of the red line: medians over the 6 experiments

<div class="views" markdown>

<div markdown>

**Without abstention**

| output | experiments | median F1 (% unlabelled) |
| :--- | :--- | :--- |
| SCARLET best | 6 | 94.3 (11.2) |
| Seurat predicted.id | 6 | 88.5 (0.0) |
| SingleR labels | 6 | 81.3 (0.0) |
| CellTypist best match | 6 | 84.3 (0.0) |
| SVM | 6 | 82.4 (0.0) |
| Harmony kNN (transductive) | 6 | 83.0 (0.0) |

</div>

<div markdown>

**With abstention**

| output | experiments | median F1 (% unlabelled) |
| :--- | :--- | :--- |
| SCARLET unique | 6 | 95.4 (20.9) |
| Seurat GTP | 6 | 95.0 (16.7) |
| SingleR pruned | 6 | 81.7 (0.2) |
| CellTypist prob match | 6 | 79.8 (9.1) |
| SVMrejection | 6 | 92.6 (20.7) |
| scmap-cluster | 6 | 72.7 (48.9) |
| scmap-cell | 6 | 78.5 (37.6) |
| scType | 6 | 72.5 (0.3) |

</div>

</div>

### Fig. 3a, right of the red line: better / equal / worse

<div class="views" markdown>

<div markdown>

**Without abstention**

| SCARLET output | Seurat predicted.id | SingleR labels | CellTypist best match | SVM | Harmony kNN (transductive) |
| :--- | :--- | :--- | :--- | :--- | :--- |
| SCARLET best | 4 / 0 / 2 (11.2 vs 0.0) | 5 / 0 / 1 (11.2 vs 0.0) | 4 / 0 / 2 (11.2 vs 0.0) | 5 / 0 / 1 (11.2 vs 0.0) | 3 / 0 / 3 (11.2 vs 0.0) |

</div>

<div markdown>

**With abstention**

| SCARLET output | Seurat GTP | SingleR pruned | CellTypist prob match | SVMrejection | scmap-cluster | scmap-cell | scType |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| SCARLET unique | 4 / 0 / 2 (20.9 vs 16.7) | 5 / 0 / 1 (20.9 vs 0.2) | 4 / 0 / 2 (20.9 vs 9.1) | 3 / 0 / 3 (20.9 vs 20.7) | 5 / 0 / 1 (20.9 vs 48.9) | 4 / 0 / 2 (20.9 vs 37.6) | 6 / 0 / 0 (20.9 vs 0.3) |

</div>

</div>

**10x → Smart-seq2 (pbmc1).** SCARLET refuses about half of these cells after refinement. The chain: the 10x reference has natural killer, dendritic and plasmacytoid dendritic cells, which the Smart-seq2 annotation folds into other labels; so in the first pass many T cells are undetermined across draws (CD4+ T 32 of 58, cytotoxic T 57 of 117); so few T cells become stable anchors (7 and 13); so the refinement scores the query against a background pooled over its cells, 69 % of them T cells, and the T cell types lose to it (refused for background). The calls made are correct (median F1 100.0). The run is faithful to the benchmark, whose rule only drops query cells whose label the reference lacks. Planned for 0.6.4: a check on stable anchors before refining a cell type, and a background balanced across the first pass's cell types.
