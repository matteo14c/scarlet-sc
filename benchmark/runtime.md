[Back to the benchmarks](index.md)

# Run time and memory, as a user runs it

**This section is our measurement, not Abdelaal et al.'s.** It reports how long each tool takes, and how much
memory it needs, when a user runs it once on one reference → query pair. It was measured by the SCARLET authors on
nodes of one cluster: SCARLET's own numbers on 2026-11-27, the other tools' on 2026-11-19. The scores in the rest of
this page do not depend on it.

- **How each tool ran:** at its documented defaults, one process per reference → query pair, from the tool's own
  input files to its written output files (reading and writing included). SCARLET: the command line of the last
  release candidate of SCARLET 0.6.3, whose calls the release does not change (the release also writes a table of
  per-gene detection rates into the model files it saves, a step these numbers do not include), every option at its
  default, including `blocklist auto` (the
  scores above use `blocklist none`); no report. The other tools: the steps of the benchmark's runner for one pair.
  On the three transfer tasks SCARLET ran `transfer` in accurate mode (the default, 10 draws) and quick mode
  (`--draws 1`, one core only); on zheng_68k it ran `learn` then `assign`, which have no draws and no modes.
- **Cores:** each tool on 1 and on 5 physical CPU cores, the run confined to them (every process of the run pinned
  to one thread of each of its cores; the cores' other threads left idle). On 1 core every tool at its defaults; on
  5 cores each tool's documented threads option set to 5 (listed in the tables; scmap and scType have none). The
  thread variables (`OMP_NUM_THREADS`, `MKL_NUM_THREADS`, `OPENBLAS_NUM_THREADS`) were not set, as in a user's
  shell, so each library sized its threads to the cores it could see.
- **Time:** wall-clock time of the whole process. **Memory:** the peak, over the run, of the memory of the whole
  process tree (proportional set size summed over every process, sampled every 0.2 s). **Cores used:** CPU time
  (user + system, every process) divided by wall time.
- **Hardware:** three nodes of one type, each with two Intel Xeon Gold 6442Y processors (48 cores, 96 hardware
  threads in all) and 504 GB of memory; one run at a time on each node, no other run of ours on the same node.
  Other users' jobs were not under our control; the node's load average at the start and end of every run is
  recorded with the measurement.
- **zheng_68k is one fold** (train on four folds, predict the fifth); the benchmark's five-fold cross-validation
  takes about five times as long for every tool.

**All tools, all four tasks** (time, peak memory; 1 core / 5 cores; the tables below give the threads option and the cores used)

| tool | PbmcBench 10Xv2 | Brain Labels34, fold 2 | PBMC → Zheng 68K | A/zheng_68k, fold 1 (learn + assign) |
| :--- | :--- | :--- | :--- | :--- |
| SCARLET, accurate (on zheng_68k: `learn` + `assign`) | 82 s, 0.5 GB / 28 s, 2.3 GB | 3.7 min, 3.1 GB / 71 s, 10.4 GB | 42 s, 1.9 GB / 22 s, 3.9 GB | 7 s, 0.5 GB / 4 s, 0.6 GB |
| SCARLET, quick (`--draws 1`, 1 core) | 16 s, 0.5 GB | 41 s, 3.0 GB | 10 s, 1.9 GB | — |
| Seurat | 33 s, 1.8 GB / 20 s, 2.3 GB | 69 s, 13.3 GB / 61 s, 13.7 GB | 81 s, 4.2 GB / 42 s, 4.5 GB | 117 s, 5.3 GB / 62 s, 5.0 GB |
| SingleR | 88 s, 1.7 GB / 28 s, 1.7 GB | 4.5 min, 9.6 GB / 119 s, 9.7 GB | 4.0 min, 2.5 GB / 58 s, 2.5 GB | 17.5 min, 2.0 GB / 4.1 min, 2.1 GB |
| CellTypist | 58 s, 6.0 GB / 28 s, 6.4 GB | 12.4 min, 9.6 GB / 3.7 min, 10.0 GB | 72 s, 17.8 GB / 41 s, 17.7 GB | 6.7 h, 17.3 GB / 84.6 min, 17.3 GB |
| SVM and SVMrejection (one run) | 7 s, 0.7 GB / 4 s, 1.4 GB | 2.1 min, 5.7 GB / 105 s, 5.7 GB | 8 s, 1.2 GB / 5 s, 1.9 GB | 4.1 min, 1.5 GB / 2.4 min, 3.5 GB |
| Harmony kNN (transductive) | 20 s, 11.4 GB / 15 s, 11.6 GB | 36 s, 22.0 GB / 30 s, 21.7 GB | 49 s, 28.6 GB / 30 s, 29.1 GB | 69 s, 30.5 GB / 52 s, 30.4 GB |
| scmap-cluster | 14 s, 3.5 GB / 12 s, 3.5 GB | 38 s, 9.0 GB / 32 s, 8.9 GB | 22 s, 3.9 GB / 20 s, 3.9 GB | 34 s, 19.3 GB / 30 s, 19.4 GB |
| scmap-cell | 2.5 min, 6.8 GB / 2.2 min, 6.8 GB | 116 s, 12.8 GB / 102 s, 12.8 GB | 43.1 min, 22.1 GB / 31.0 min, 22.1 GB | 5.2 min, 18.4 GB / 4.1 min, 18.4 GB |
| scType | 36 s, 2.4 GB / 20 s, 2.2 GB | 80 s, 10.1 GB / 68 s, 10.2 GB | 2.4 min, 7.2 GB / 75 s, 7.2 GB | 32 s, 2.2 GB / 23 s, 2.2 GB |

**PbmcBench 10Xv2 (reference 10x v2 pbmc1, 6,444 cells; query the 7 other PBMC datasets, 16,710 cells; 22,280 genes)**

| tool | 1 core | 5 cores | threads option at 5 cores | cores used (1 / 5) |
| :--- | :--- | :--- | :--- | :--- |
| SCARLET, accurate (10 draws, the default) | 82 s, 0.5 GB | 28 s, 2.3 GB | --ncpu 5 (default 1) | 0.9 / 3.4 |
| SCARLET, quick (`--draws 1`) | 16 s, 0.5 GB | — | 1 core only | 0.9 /  |
| Seurat | 33 s, 1.8 GB | 20 s, 2.3 GB | future::plan("multicore", workers = 5) (default: no plan, sequential) | 1.0 / 1.9 |
| SingleR | 88 s, 1.7 GB | 28 s, 1.7 GB | num.threads = 5 (default 1) | 1.0 / 2.9 |
| CellTypist | 58 s, 6.0 GB | 28 s, 6.4 GB | celltypist.train(n_jobs = 5) (default None = 1) | 1.0 / 2.1 |
| SVM and SVMrejection (one run) | 7 s, 0.7 GB | 4 s, 1.4 GB | CalibratedClassifierCV(n_jobs = 5) (default None); LinearSVC has none | 0.7 / 1.8 |
| Harmony kNN (transductive) | 20 s, 11.4 GB | 15 s, 11.6 GB | run_harmony(ncores = 5) (default 0), KNeighborsClassifier(n_jobs = 5) (default None) | 1.0 / 1.7 |
| scmap-cluster | 14 s, 3.5 GB | 12 s, 3.5 GB | none | 1.0 / 1.1 |
| scmap-cell | 2.5 min, 6.8 GB | 2.2 min, 6.8 GB | none | 1.0 / 1.1 |
| scType | 36 s, 2.4 GB | 20 s, 2.2 GB | none | 1.0 / 1.7 |

**Brain Labels34, fold 2 (reference mouse ALM, 8,128 cells; query mouse V1, 12,552 cells; 41,718 genes; about 11,000 genes detected per cell)**

| tool | 1 core | 5 cores | threads option at 5 cores | cores used (1 / 5) |
| :--- | :--- | :--- | :--- | :--- |
| SCARLET, accurate (10 draws, the default) | 3.7 min, 3.1 GB | 71 s, 10.4 GB | --ncpu 5 (default 1) | 1.0 / 4.0 |
| SCARLET, quick (`--draws 1`) | 41 s, 3.0 GB | — | 1 core only | 0.9 /  |
| Seurat | 69 s, 13.3 GB | 61 s, 13.7 GB | future::plan("multicore", workers = 5) (default: no plan, sequential) | 1.0 / 1.7 |
| SingleR | 4.5 min, 9.6 GB | 119 s, 9.7 GB | num.threads = 5 (default 1) | 1.0 / 2.2 |
| CellTypist | 12.4 min, 9.6 GB | 3.7 min, 10.0 GB | celltypist.train(n_jobs = 5) (default None = 1) | 1.0 / 3.9 |
| SVM and SVMrejection (one run) | 2.1 min, 5.7 GB | 105 s, 5.7 GB | CalibratedClassifierCV(n_jobs = 5) (default None); LinearSVC has none | 1.0 / 1.1 |
| Harmony kNN (transductive) | 36 s, 22.0 GB | 30 s, 21.7 GB | run_harmony(ncores = 5) (default 0), KNeighborsClassifier(n_jobs = 5) (default None) | 0.9 / 1.7 |
| scmap-cluster | 38 s, 9.0 GB | 32 s, 8.9 GB | none | 1.0 / 1.1 |
| scmap-cell | 116 s, 12.8 GB | 102 s, 12.8 GB | none | 1.0 / 1.2 |
| scType | 80 s, 10.1 GB | 68 s, 10.2 GB | none | 1.0 / 1.3 |

**PBMC → Zheng 68K (reference 10x v2 pbmc1, 6,444 cells, 22,280 genes; query all 65,943 Zheng 68K cells, 20,387 genes)**

| tool | 1 core | 5 cores | threads option at 5 cores | cores used (1 / 5) |
| :--- | :--- | :--- | :--- | :--- |
| SCARLET, accurate (10 draws, the default) | 42 s, 1.9 GB | 22 s, 3.9 GB | --ncpu 5 (default 1) | 0.9 / 2.4 |
| SCARLET, quick (`--draws 1`) | 10 s, 1.9 GB | — | 1 core only | 0.9 /  |
| Seurat | 81 s, 4.2 GB | 42 s, 4.5 GB | future::plan("multicore", workers = 5) (default: no plan, sequential) | 1.0 / 2.2 |
| SingleR | 4.0 min, 2.5 GB | 58 s, 2.5 GB | num.threads = 5 (default 1) | 1.0 / 3.8 |
| CellTypist | 72 s, 17.8 GB | 41 s, 17.7 GB | celltypist.train(n_jobs = 5) (default None = 1) | 1.0 / 1.9 |
| SVM and SVMrejection (one run) | 8 s, 1.2 GB | 5 s, 1.9 GB | CalibratedClassifierCV(n_jobs = 5) (default None); LinearSVC has none | 0.7 / 1.7 |
| Harmony kNN (transductive) | 49 s, 28.6 GB | 30 s, 29.1 GB | run_harmony(ncores = 5) (default 0), KNeighborsClassifier(n_jobs = 5) (default None) | 1.0 / 1.9 |
| scmap-cluster | 22 s, 3.9 GB | 20 s, 3.9 GB | none | 1.0 / 1.1 |
| scmap-cell | 43.1 min, 22.1 GB | 31.0 min, 22.1 GB | none | 1.0 / 1.5 |
| scType | 2.4 min, 7.2 GB | 75 s, 7.2 GB | none | 1.0 / 2.1 |

**A/zheng_68k, fold 1 (train 52,780 cells, test 13,163 cells; 20,387 genes); SCARLET: `learn` + `assign` only (no transfer, no draws, no quick mode)**

| tool | 1 core | 5 cores | threads option at 5 cores | cores used (1 / 5) |
| :--- | :--- | :--- | :--- | :--- |
| SCARLET (`learn` + `assign`) | 7 s, 0.5 GB | 4 s, 0.6 GB | --ncpu 5 (default 1) | 0.6 / 1.3 |
| Seurat | 117 s, 5.3 GB | 62 s, 5.0 GB | future::plan("multicore", workers = 5) (default: no plan, sequential) | 1.0 / 2.4 |
| SingleR | 17.5 min, 2.0 GB | 4.1 min, 2.1 GB | num.threads = 5 (default 1) | 1.0 / 4.2 |
| CellTypist | 6.7 h, 17.3 GB | 84.6 min, 17.3 GB | celltypist.train(n_jobs = 5) (default None = 1) | 1.0 / 4.5 |
| SVM and SVMrejection (one run) | 4.1 min, 1.5 GB | 2.4 min, 3.5 GB | CalibratedClassifierCV(n_jobs = 5) (default None); LinearSVC has none | 1.0 / 1.9 |
| Harmony kNN (transductive) | 69 s, 30.5 GB | 52 s, 30.4 GB | run_harmony(ncores = 5) (default 0), KNeighborsClassifier(n_jobs = 5) (default None) | 1.0 / 2.7 |
| scmap-cluster | 34 s, 19.3 GB | 30 s, 19.4 GB | none | 1.0 / 1.2 |
| scmap-cell | 5.2 min, 18.4 GB | 4.1 min, 18.4 GB | none | 1.0 / 1.2 |
| scType | 32 s, 2.2 GB | 23 s, 2.2 GB | none | 1.0 / 1.9 |

**Conversion for the R tools** (Seurat, SingleR, scmap, scType read an R object; the .h5ad files converted once per task, 1 core; not included in the tools' rows above)

| task | conversion |
| :--- | :--- |
| PbmcBench 10Xv2 | 17 s, 1.0 GB |
| Brain Labels34, fold 2 | 3.7 min, 7.7 GB |
| PBMC → Zheng 68K | 27 s, 1.9 GB |
| A/zheng_68k, fold 1 | 28 s, 1.8 GB |

- **No tool used more than its cores:** at most 1.0 core on 1 core, at most 4.5 on 5. The tools that use the 5
  cores most are SingleR (2.2–4.2), CellTypist (1.9–4.5) and SCARLET (1.3–4.0); Seurat, Harmony kNN and scType
  (through their libraries' threads) use about 2; scmap barely more than 1.
- **SCARLET's memory grows with the draws it runs side by side:** on Brain 3.1 GB with `--ncpu 1`, 6.0 GB with 2, 10.4 GB with 5 and 17.7 GB with 10 (3.7 min, 2.2 min, 71 s and 52 s); the transfer guide gives these figures for each task and how to choose `--ncpu` by the free memory.
- **Two predictions we made before these runs did not hold.** Against an earlier measurement in which 5 CPUs meant 5 hardware threads on 2.5 cores, Harmony kNN took longer on 5 cores (Brain 30 s against 22 s; zheng_68k 52 s against 38 s), and SingleR and Seurat were no faster on Brain (119 s against 119 s; 61 s against 63 s), while both were 11–20 % faster on the other tasks. We have not established why.
- **CellTypist on zheng_68k** (more than 50,000 training cells) trains with a single-threaded solver, one logistic
  regression per cell type: 6.7 h on 1 core, 85 min on 5.
- **SVM on Brain:** SVMrejection was not run on this task (as in the benchmark: a cell type has fewer than 3
  reference cells, too few for its 3-fold calibration), so only LinearSVC ran, which has no threads option.
- **Not included:** installing the tools, and the conversion of the .h5ad files for the R tools (table above).


