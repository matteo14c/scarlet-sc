# Benchmarks

SCARLET was benchmarked by rerunning one published benchmark of cell-annotation methods: Abdelaal et al., *Genome
Biology* 20, 194 (2019), "the benchmark" below, with its datasets as distributed on Zenodo (record 3357167). Every
tool ran at its own defaults, and every score is shown beside the share of cells the tool left unlabelled, because
a tool can raise its score on the cells it labels by refusing more cells.

**These results are not independent of SCARLET's defaults.** SCARLET's default settings were developed on these
tasks. Read the tables as a description of the defaults on the data they were developed on, not as an independent
test.

## The experiments

Each page follows one figure of the benchmark, in the benchmark's order, and shows SCARLET 0.6.3 once in each view:
**SCARLET best**, which labels every cell it can place (refusing cells its background explains better than any cell
type), beside the tools that label every cell, and **SCARLET unique**, which also refuses ambiguous cells, beside the
tools that can refuse. For `transfer`, these are the final calls (after refinement, accurate mode, the defaults).

| page | the benchmark's figure | what is transferred to what | scored |
| :--- | :--- | :--- | :--- |
| [Within a dataset](fig1.md) | Fig. 1 | 12 datasets, 5-fold cross-validation: learn on four fifths, annotate the fifth | median F1 over cell types (% unlabelled) |
| [PBMC, between protocols and samples](fig3a.md) | Fig. 3a | 42 pairs of seven protocols on one sample (pbmc1), and 6 pairs from one sample to another (pbmc1 → pbmc2) | as above, per pair |
| [Mouse brain](fig4.md) | Fig. 4a and 4b | primary visual cortex ↔ anterior lateral motor cortex, at 3 and at 34 cell populations | as above, per direction |
| [Pancreas](fig5a.md) | Fig. 5a | four pancreas datasets, each annotated from the other three | as above, per held-out dataset |
| [CellBench](figS8.md) | Supp. Fig. S8 | two protocols on the same five cell lines, each direction | as above |
| [Rejection](fig6.md) | Fig. 6, and our extension | Part 1: a reference of one tissue, a query of another (should be refused); Part 2: a population held out of the reference; our Part 3: both sides at once | % refused; Part 3 also the F1 of the cells that should be called |

Also: [how the benchmark was run](methods.md) (how to read the tables, every tool's settings and versions), [run time
and memory](runtime.md) (every tool, 1 and 5 cores, four tasks), [SCARLET's other outputs](further.md) (forced,
first pass, quick mode) and [how to reproduce](reproduce.md) (with every input file's SHA-256).

## At a glance

SCARLET's medians over each figure's experiments, median F1 (% unlabelled), from the figure pages (where every tool
is shown):

| figure | experiments | SCARLET best | SCARLET unique |
| :--- | :--- | :--- | :--- |
| Fig. 1, within a dataset | 12 | 98.9 (1.8) | 99.0 (3.2) |
| Fig. 3a, between protocols (pbmc1) | 42 | 94.0 (33.4) | 96.1 (50.7) |
| Fig. 3a, between samples (pbmc1 → pbmc2) | 6 | 94.3 (11.2) | 95.4 (20.9) |
| Fig. 4a, brain, 3 populations | 2 | 100.0 (1.3) | 100.0 (1.8) |
| Fig. 4b, brain, 34 populations | 2 | 91.9 (1.7) | 92.6 (3.2) |
| Fig. 5a, pancreas | 4 | 99.6 (11.9) | 99.8 (18.8) |
| Supp. Fig. S8, CellBench | 2 | 100.0 (1.6) | 100.0 (2.8) |

On the rejection controls of Fig. 6 Part 1, where nearly every query cell should be refused, SCARLET unique refuses
93.6–100 % of the cells. Across technologies SCARLET refuses more cells than most other tools (Fig. 3a: only
scmap-cluster and SVMrejection refuse more): it gives no call rather than a wrong one, and the share it refuses is
always shown beside its score.

## The gene blocklist

The benchmark ran every tool, SCARLET included, **without a gene blocklist**, so that all tools saw the same genes
(none of the other tools removes genes). In normal use keep SCARLET's blocklist on: its default, `blocklist auto`,
keeps groups of genes that describe a cell's state rather than its type (ribosomal, mitochondrial, stress, sex and
cell-cycle genes) out of the markers. What it changes in these tables, measured on SCARLET 0.6.3 (little, within about 1 point everywhere), what the
list holds and how to choose it: [the blocklist](../docs/blocklist.md#what-the-blocklist-changes-in-the-benchmark).

## References

1. Abdelaal T. et al. A comparison of automatic cell identification methods for single-cell RNA sequencing data.
   *Genome Biology* 20, 194 (2019). doi:[10.1186/s13059-019-1795-z](https://doi.org/10.1186/s13059-019-1795-z).
   Data: Zenodo record 3357167 (CC BY 4.0).
