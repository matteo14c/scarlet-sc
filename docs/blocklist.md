# The blocklist

Before it searches for markers, SCARLET keeps some genes out of them: genes that describe what happened to a cell
or what state it is in, rather than what type it is. A ribosomal or a dissociation-stress gene can be detected more
often in one cell type of a reference only because of how that sample was handled; as a marker it would then call
cells by their handling instead of their type. The blocklist is applied by `learn` and `transfer`; it removes genes
from the marker search only (the counts are not changed, and `assign` scores the model's markers as they are).

Where it is set: `--blocklist` in [learn](guide/learn.md#blocklist) and [transfer](guide/transfer.md);
`--blocklist-modules` in [Advanced options](advanced.md#blocklist_modules). Its messages are explained in
[Troubleshooting](troubleshooting.md#the-blocklist).

## What is on it

The bundled lists (human and mouse) are divided into six gene groups (modules):

| group | human genes | applied | what it holds | how the genes were chosen |
| :--- | :--- | :--- | :--- | :--- |
| `core_stress` | 21 | always | immediate-early and heat-shock genes (FOS, FOSB, JUN, JUNB, JUND, EGR1, ATF3, DUSP1, IER2, IER3, ZFP36, PPP1R15A, HSPA1A, HSPA1B, HSPA6, DNAJB1, DNAJA1, HSP90AA1, HSPH1, HSPB1, BAG3) | induced by tissue handling in at least two of three published dissociation signatures (van den Brink et al. 2017; O'Flanagan et al. 2019; Denisenko et al. 2020), or strictly heat-inducible |
| `core_ribosomal` | 265 | always | the 79 ribosome-subunit genes and 186 ribosomal-protein pseudogenes | HGNC gene families and locus types; paralogs with their own functions (e.g. RPL3L, RPS4Y2), the S6 kinases and read-through transcripts are not on it |
| `core_mitochondrial` | 49 | always | the mitochondrially encoded genes and MTRNR2L1–12 | the mitochondrial chromosome (HGNC) |
| `core_sex` | 8 | always | XIST, TSIX and six Y-linked genes (DDX3Y, EIF1AY, KDM5D, USP9Y, UTY, RPS4Y1) | HGNC |
| `cell_cycle` | 97 | unless a label names dividing cells | the S and G2/M phase genes of Tirosh et al. 2016, with current gene symbols | as distributed with Seurat (`cc.genes`) |
| `stress_extended` | 14 | never by default | activation and quiescence genes with weaker evidence of handling stress (EGR2, EGR3, NR4A1, NFKBIA, SOCS3, DUSP2, BTG1, BTG2, ZFP36L1, ZFP36L2, IER5, DNAJA4, DNAJB4, RHOB) | in only one or two of the studies above, or markers of real cell states |

The **mouse** list mirrors the human groups gene for gene through orthology (MGI homology classes, with Ensembl as a
cross-check), plus the mouse ribosomal-protein pseudogenes, which have no human counterparts, and the withdrawn
symbols MGI lists for them so that older gene names still match. Its groups hold 51 stress, 1,195 ribosomal, 100
mitochondrial, 16 sex-linked and 292 cell-cycle rows, and 38 rows of `stress_extended`, counting the aliases.

Genes are matched by Ensembl id or by symbol (case-insensitive), whichever the data have.

## How SCARLET chooses

**The species** (`blocklist`, `--blocklist`; default `auto`): from the gene names. Mostly upper-case symbols (`CD3E`)
or `ENSG` ids mean human; mostly capitalised symbols (`Cd3e`) or `ENSMUSG` ids mean mouse. When the names do not tell,
SCARLET stops with an error and asks you to choose.

| value | what it does |
| :--- | :--- |
| `auto` (default) | the bundled human or mouse list, chosen from the gene names |
| `human`, `mouse` | the bundled list for that species |
| `none` | no gene is kept out |
| a file path | your own tab-separated list, with a `Symbol` and/or an `ENSG_ID` column (with a `Category` column holding the group names above, `blocklist_modules` acts on it too); a path that does not exist is an error |

**The groups** (`blocklist_modules`, `--blocklist-modules`; default `auto`): the four `core_` groups always apply and
`stress_extended` never does. The `cell_cycle` group applies unless one of the reference's labels names a dividing
or proliferative state or stage (for example "cycling", "proliferating", "mitotic", "progenitor", "erythroblast"):
in a reference that labels such cells, the cell-cycle genes are what tells them apart, so they must stay available
as markers. The rule reads the labels you learn on, never a result. The pre-flight, the run's summary and the report
say which list and which groups were used, and, when the cell-cycle genes stayed in, which labels kept them.

## Changing it

- **Another set of groups:** `blocklist_modules="all"` applies every group (`stress_extended` included); a
  comma-separated list applies exactly those groups, e.g. `--blocklist-modules core_stress,core_ribosomal`.
- **No blocklist:** `blocklist="none"` (`--blocklist none`).
- **Your own list:** `blocklist="my_list.tsv"`, a tab-separated file with a `Symbol` and/or `ENSG_ID` column. To add
  genes to the bundled list, copy it from the package (`scarlet/scarlet_blocklist_universal.tsv` for human,
  `scarlet/scarlet_blocklist_mouse.tsv` for mouse), add your rows and pass the file.

Any of these is recorded in the model and the summary, and the report notes that the run is not at the package
defaults.

## In the benchmarks

The benchmarks ran every tool, SCARLET included, **without** a blocklist (`blocklist none`), so that all tools saw the
same genes: none of the other tools removes genes. In normal use keep it on. See [Benchmarks](https://github.com/matteo14c/scarlet-sc/blob/v0.6.3/benchmark/index.md).

### What the blocklist changes in the benchmark

The same benchmark run again with SCARLET at its package default (`blocklist auto`), everything else as benchmarked.
Each cell: median F1 over the experiments, with the % of cells left unlabelled in brackets (as in the benchmark pages).

| figure | experiments | SCARLET best, no blocklist | with blocklist | SCARLET unique, no blocklist | with blocklist |
| :--- | :--- | :--- | :--- | :--- | :--- |
| Fig. 1, within a dataset | 12 | 98.9 (1.8) | 98.9 (1.8) | 99.0 (3.2) | 99.0 (3.2) |
| Fig. 3a, between protocols | 42 | 94.0 (33.4) | 94.3 (33.3) | 96.1 (50.7) | 96.0 (50.4) |
| Fig. 3a, between samples | 6 | 94.3 (11.2) | 94.0 (11.0) | 95.4 (20.9) | 95.5 (21.0) |
| Fig. 4a, brain, 3 populations | 2 | 100.0 (1.3) | 100.0 (1.3) | 100.0 (1.8) | 100.0 (1.7) |
| Fig. 4b, brain, 34 populations | 2 | 91.9 (1.7) | 92.1 (1.7) | 92.6 (3.2) | 92.6 (3.3) |
| Fig. 5a, pancreas | 4 | 99.6 (11.9) | 99.5 (11.9) | 99.8 (18.8) | 99.8 (19.0) |
| Supp. Fig. S8, CellBench | 2 | 100.0 (1.6) | 100.0 (1.4) | 100.0 (2.8) | 100.0 (2.8) |

Fig. 6, % of the query refused. Part 1 (a reference of one tissue, a query of another; controls a, b, c, d): SCARLET
unique 99.4, 93.6, 99.2, 100.0 without the blocklist and 99.4, 93.8, 99.2, 100.0 with it; SCARLET best 97.8, 80.3,
96.3, 99.8 and 97.3, 81.2, 96.4, 99.8. Part 2 (a population held out of the reference; hold-outs i, ii, iii): SCARLET
unique 54.1, 58.2, 88.7 and 54.2, 58.9, 89.0. Part 3 (held-out cells refused / other cells refused, SCARLET unique):
51.8 / 19.6, 50.3 / 25.2, 81.3 / 63.2 and 51.8 / 19.6, 50.4 / 25.4, 81.0 / 63.1.

**Reading.** On these data the blocklist makes little difference: every median F1 moves by 0.3 points or less, every
share unlabelled by 0.3 points or less, and every share refused in Fig. 6 by 0.9 points or less, in both directions.
These differences are of the size that the order of tied marker genes alone produces (up to about 1 point), so the
benchmark neither shows a gain nor a cost from the blocklist. Its purpose is elsewhere: it keeps markers that describe
how a sample was handled (stress, dissociation) or a cell's state out of a reference, which matters most when the
reference and the query were handled differently, which the benchmark was not designed to test.

## References

- van den Brink SC, et al. Single-cell sequencing reveals dissociation-induced gene expression in tissue
  subpopulations. *Nature Methods* 14, 935–936 (2017).
- O'Flanagan CH, et al. Dissociation of solid tumor tissues with cold active protease for single-cell RNA-seq
  minimizes conserved collagenase-associated stress responses. *Genome Biology* 20, 210 (2019).
- Denisenko E, et al. Systematic assessment of tissue dissociation and storage biases in single-cell and
  single-nucleus RNA-seq workflows. *Genome Biology* 21, 130 (2020).
- Tirosh I, et al. Dissecting the multicellular ecosystem of metastatic melanoma by single-cell RNA-seq. *Science*
  352, 189–196 (2016).
