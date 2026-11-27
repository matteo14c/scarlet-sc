# Reports

A report is **one self-contained HTML page** (open it in any browser; no network or server needed) and, beside it,
a **bundle** of tables (`bundle/tables/*.tsv`) holding every value the page shows. `report` reads a finished run's
outputs; it changes no call and no output of the run, and shows only values SCARLET computed (UMAPs and the
comparison with labels you supply are marked as display only). It needs the `report` extra:
`pip install "scarlet-sc[report]"`.

| what you report on | the report | page |
| :--- | :--- | :--- |
| a model file (`.scarlet`), written by `learn` | the **learn report**: the reference, the markers, and how well the model calls its own reference (cross-validation) | [learn report](learn.md) |
| an annotated `.h5ad`, written by `assign` | the **assign report**: the calls and the refusals of a query without labels | [assign report](assign.md) |
| a transfer output folder | the **transfer report**: both sides, the first pass and the final calls, refinement, agreement across draws | [transfer report](transfer.md) |

```bash
scarlet report reference.scarlet                                   # learn report
scarlet report query_annotated.h5ad --model reference.scarlet      # assign report
scarlet report transfer_out/                                       # transfer report
```
```python
import scarlet
scarlet.report("reference.scarlet")
```

`learn(..., report=True)` and `transfer(..., report=True)` (`--report` on the command line) write the report when
the run ends.

## Why a report can take time

- **The learn report re-runs SCARLET.** Its self-assessment splits the reference into five folds and, for each,
  learns on the other four and assigns the fold left out: five learns and five assigns. A transfer report includes
  the learn report of its reference model (its *supplement*), so it does the same. Switch it off with
  `self_assessment=False` (`--no-self-assessment`) when you only need the rest.
- **The UMAPs** are computed for the page, by default twice: on the model's own genes and on all genes. On large data
  each takes tens of seconds to minutes; `umap_all_genes=False` (`--umap-all-genes no`) skips the all-genes one, and
  `umap_cells` caps the cells shown (50,000 by default; a seeded subsample above it).

At the default verbosity a report prints one line per fold and one when the page is written.

## Common to every report

- **The page opens on a summary** of the result and **ends with the provenance**: the run's own summary, the
  SCARLET version, the inputs with the SHA-256 the run recorded (the report checks it against the file it reads and
  says whether they are the same), the technology preset, the blocklist and its gene groups, and whether the run
  used the package defaults. When it did not, every option that differs is listed; when one is an
  [advanced option](../advanced.md), the report adds that the published benchmarks used the defaults of the advanced
  options and do not apply to the run.
- **Colours:** one fixed colour per cell type in every panel; refused cells are grey. The palette is
  colour-blind-safe for up to 24 cell types (a footnote on the page says which palette and when it is not).
- **Tables** have readable headers; hovering a header shows the column's name in the bundle's TSV and, where
  useful, what it means. Percentages have one decimal, counts are whole numbers, FDRs are in scientific notation.
  Wide tables scroll sideways.
- **Refusal reasons:** every refused cell has one reason, in the order SCARLET applies its steps: *background*, then
  *flip test* (when the flip test refuses), then *ambiguity* (the unique call refused while the best call is kept),
  and in a transfer *undetermined across draws*. The report stops rather than leave a refusal unexplained.
- **Labels you supply** (transfer and assign reports): `query_labels` (`--query-labels`) names a column of the
  query's `obs` that holds labels you already know, as in the tutorials. The report then adds a section, marked as
  labels supplied by you, comparing SCARLET's calls with them: the share refused and the median F1 per label on the
  called cells, the same per label, and a map of labels against calls. SCARLET never reads these labels; without
  `query_labels` the section is absent.

## The UMAPs

Display only: never used for a call. Each UMAP has its own switches:

- **Output** (where there are several): which calls colour the cells.
- **View** (where there are several): *label*; *margin* (learn report: the call's colour darker the larger its
  margin) or *agreement* (transfer: darker the more draws agree); *ambiguity* (the boundary cells in full colour);
  *anchors* (transfer: stable anchors outlined); *unlabelled* (labelled cells faint, refused cells grey with one
  marker per reason: circle background, diamond ambiguity, x undetermined across draws, triangle flip test).
- **Highlight:** one cell type in full colour and the others faint; none by default, so every cell type is in full
  colour.
- **Space:** the genes the UMAP is computed on. By default the **model's own genes** (its universe: the marker
  genes SCARLET decides on), read as SCARLET reads them, detected or not, with neighbours by Jaccard distance; the
  second space is **all genes** (scanpy's standard steps: normalise, log, highly variable genes, PCA, neighbours).
  Each title names its space.

**Which space to look at.** On the model's genes the cells are placed by what SCARLET compares, so cells it cannot
tell apart sit together. All genes is the place to look for **populations the model lacks**: such cells sit with
their nearest relative on the model's genes, as SCARLET calls them, and may stand apart on all genes. In one
donor's PBMC (a Smart-seq2 reference without natural killer cells, a 10x v2 query with 429 of them), the natural
killer cells sit inside the cytotoxic T cells on the model's genes (SCARLET calls 416 of them cytotoxic T), but form
their own group at the edge of the T cells on all genes.

Without WebGL in the browser (switched off, or not available), the UMAPs are drawn as SVG on at most 10,000 cells,
sampled within each cell type so that small populations stay whole; the title says so.
