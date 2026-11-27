# SCARLET

SCARLET annotates single-cell RNA-seq data from a labelled reference. For each cell type it learns the marker genes
its cells consistently detect, and calls each new cell the cell type that explains it best, **or refuses it** when
no cell type explains it better than the background, or two explain it equally well. Every refusal has a reason.

> **Public beta, v0.6.3.** Stable: `learn`, `assign`, `transfer` and `report`. Coming later: tools to merge, split and
> discover cell identities. Not in scope: cross-species annotation.

**Install** (Python 3.11 to 3.13): `pip install "scarlet-sc[report] @ git+https://github.com/matteo14c/scarlet-sc@v0.6.3"`;
more in [Getting started](getting-started.md). **Tutorials:** [three notebooks](tutorials.md).

| command | what it does | guide |
| :--- | :--- | :--- |
| `learn` | builds a model (one file) from a labelled reference | [learn](guide/learn.md) |
| `assign` | annotates new cells with a model | [assign](guide/assign.md) |
| `transfer` | annotates a query from a reference in one step, across samples, depths and technologies | [transfer](guide/transfer.md) |
| `report` | writes the HTML report of a finished run | [reports](reports/index.md) |

Each command is available from the command line (`scarlet learn ...`) and from Python (`scarlet.learn(...)`), with the
same options and the same defaults.

## How SCARLET was benchmarked

We reran the benchmark of Abdelaal et al. (*Genome Biology*, 2019): annotation within a dataset, transfer across
technologies and across samples (PBMC, mouse brain, pancreas, CellBench), and refusal of cells the reference lacks;
every tool at its defaults, every score beside the share of cells it left unlabelled; plus our own test on
populations missing from the reference. In our view SCARLET holds up well. SCARLET's defaults were developed on these
data, so the results are not an independent test.

- **Strengths.** Accuracy in the range of the best tools (within a dataset, median F1 98.9 with 1.8 % unlabelled,
  against 97.7 to 98.9 for the tools that label every cell). Among the most reliable at refusing cells it cannot place
  (it refuses 93.6 to 100 % of the cells of a tissue its reference does not have; the only other tool above 90 % in
  all four such tests is scmap-cluster (99.0 to 100 %), which also refuses most of the cells it could place). Bit-identical within a fixed environment:
  with the same operating system, Python, maths libraries and CPU type, the same inputs, options and seeds give the
  same output, byte for byte, on any number of CPUs, and every code change is checked against this. Results are
  identical on machines with the same CPU vector instructions; in our tests, a laptop without AVX-512 matched a server
  with it switched off byte for byte. Where the instructions differ, `assign` still gives identical calls, and
  `transfer` can differ in a few cells (21 of 2,990 in our PBMC example): almost all at the 8-of-10 agreement
  threshold, and occasionally (1 of the 21) a cell beyond it, because refinement learns from the stable anchors. `learn` + `assign` in seconds
  (6 s and 0.4 GB for 53,000 reference and 13,000 query cells, on one core). Low memory on one or two CPUs.
- **Limitations.** Across technologies it refuses more cells than most tools (a median of 33 % unlabelled for its best
  call across the 42 PBMC protocol pairs): no call rather than a wrong one. A population missing from the reference
  but close to one it has is mostly called as that relative, not refused. In our view this is a limit of
  reference-based annotation itself; the other tools also call many such cells as the relative, though some refuse
  more of them than SCARLET, and with them more of the cells they could place (CD4+ T cells held out of a reference:
  Seurat refuses 82.1 % of them and 49.0 % of the other cells, SVMrejection 74.2 % and 51.3 %, SCARLET 50.3 % and
  25.2 %). It is the reason for SCARLET's planned split, merge and discover tools. A transfer in accurate mode
  takes mid-pack time (3.7 min on one core for 8,000 → 12,500 deep mouse-brain cells; the tools ranged from 36 s to
  12 min), and its memory grows with the CPUs it uses (3.1 GB on one, 10.4 GB on five).

The full benchmarks, figure by figure, with run time and memory: [Benchmarks](https://github.com/matteo14c/scarlet-sc/blob/v0.6.3/benchmark/index.md). What SCARLET
does not do well yet: [Known limitations](limitations.md).

## Citing SCARLET

A paper is in preparation. Until then, please cite the software and the datasets you use:

Traversa D, Kaftalli J, Canevotti A, Rossi F, Pesole G, Chiara M. *SCARLET: reliable and interpretable cell-type
annotation for single-cell RNA-seq.* Software, version 0.6.3, 2026. https://github.com/matteo14c/scarlet-sc

`CITATION.cff` in the repository holds the same.
