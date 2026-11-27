# transfer

`transfer` annotates a query with the labels of a reference of the same species, in one step: it matches the
depths of the two, learns on the reference, annotates the query (the **first pass**), learns again on the query's
own stably anchored cells and annotates the query with that model (the **refinement**), and checks each cell type.
The whole transfer is repeated in **draws** (10 by default), each with its own seed, and each cell's final call is
the outcome at least 8 of the 10 draws agree on. It is the command for a query of another depth, protocol or
technology than the reference.

```bash
scarlet transfer --reference reference.h5ad --labels cell_type --query query.h5ad --output transfer/
```
```python
import scarlet

if __name__ == "__main__":   # needed in a script when ncpu is above 1; harmless otherwise
    res = scarlet.transfer("reference.h5ad", "cell_type", "query.h5ad", "transfer/", ncpu=4)
```

In a Python **script**, call `transfer` under `if __name__ == "__main__":` when `ncpu` is above 1: its draws run in
new Python processes, which import the script; without the guard SCARLET stops with a message saying so. Notebooks
and the command line need nothing.

## Inputs

- `reference`: the labelled reference, an `.h5ad` file (or an AnnData) of counts. Reference cells without a label
  are not used.
- `labels`: the `obs` column of the reference with the cell types.
- `query`: the cells to annotate, an `.h5ad` file (or an AnnData) of counts, with the same kind of gene names.
- `output`: the output folder.

Depth matching reads expression values, so both inputs must be non-negative counts.

## Origins

A query often pools cells of several batches, samples or technologies. `--origin-col` names a column of the query's
`obs` that says where each cell comes from (never the labels). Depth matching, the refinement and the transfer check
are then done separately for each origin. Without it the whole query is one origin, called `all`.

## The steps

### Accurate and quick mode

Depth matching and the marker search are random steps: a different seed thins different detections and draws
different cells for the markers, and some cells change call from one seed to another. `transfer` therefore repeats
depth matching, the first pass and the refinement in **draws**, each with its own seed from `--seeds` (by default 1,
2, ... up to the number of draws).

- **Accurate mode** (the default, `--draws 10`). Each cell's final call is the **consensus** of the 10 draws,
  computed separately for the best and the unique call: the outcome most draws give, where "unassigned" counts as an
  outcome, is the final call when **at least 8 of the 10 draws** give it. Otherwise the cell is `unassigned` with
  the reason "undetermined across draws". This is a fixed rule: if a cell's two leading outcomes were equally likely
  in every draw, 8 or more of 10 draws would agree on one of them by chance in about 5.5 % of cells (one-sided p =
  0.055). Each output column `scarlet_best_agreement` / `scarlet_unique_agreement` says how many draws gave the
  cell's outcome. With another number of draws, a final call needs 80 % of them, rounded up: 16 of 20, 4 of 5, 2 of
  2.
- **Quick mode** (`--draws 1`). The transfer runs once: the same method, one draw. Its calls are that draw's, so
  there is no consensus, no per-cell agreement (the agreement columns are not written) and no cell is "undetermined
  across draws". What quick mode gives up: its calls depend on the seed (`--seeds 2` gives another draw, and some
  cells change call), and its stable anchors are that draw's unique calls, so the refinement learns from more,
  less reliable anchors. On the published benchmarks, quick mode refused far fewer cells of the negative controls
  (cells of a tissue absent from the reference) than accurate mode after refinement: 42–71 % against 81–100 %
  (refined best), and the median F1 of a test sample moved by more than 2 points from one seed to another in about
  three quarters of the samples. It takes a fifth to a quarter of the time with one CPU (see
  [Run time, memory and CPUs](#run-time-memory-and-cpus)).

### Reproducibility

The same inputs, options and seeds give the same output, byte for byte, whatever `--ncpu`, on the same operating
system, Python, maths libraries and CPU type. Results are identical on machines with the same CPU vector
instructions; in our tests, a laptop without AVX-512 matched a server with it switched off byte for byte. Where the
instructions differ, `assign` still gives identical calls, and `transfer` can differ in a few cells (21 of 2,990 in
our PBMC example): almost all at the 8-of-10 agreement threshold, and occasionally (1 of the 21) a cell beyond it,
because refinement learns from the stable anchors. Another seed list gives
comparable results, not the same: on the published benchmarks' transfer tasks, two seed lists (1–10 and 1001–1010)
gave the same final call for 77 % to 99.8 % of the cells of a task, and the task's median F1 differed by at most 1.9
points (% unlabelled by at most 4 points). Small references vary more: on test samples whose cell types have only a
few hundred reference cells each, the median F1 moved by more than 2 points between the two seed lists in about one
sample in five.

### Depth matching

Reference and query of different depth differ in which genes they detect at all, and SCARLET reads detection. With
`--depth-match yes` (the default), the deeper side is thinned towards the shallower one before the first pass,
separately for each origin. The depth of a side is the median number of genes detected per cell, counted on the
genes both sides detect; with equal depth nothing is thinned. The thinning fits, across the shared genes, an
increasing curve from each gene's mean expression on the deep side to its detection rate on the shallow side, and
keeps each detection on the deep side with the probability that brings its gene to the fitted rate. It does not use
the labels. Each draw thins with its own seed (per origin), so each draw sees a different thinning.

Thinning affects the first pass only: the refinement learns on the query's own cells, unthinned. When the reference
is thinned, one reference model is learned per origin. With depth matching on, each origin is annotated
separately. `--depth-match no` skips this step.

### Technology preset

With `--technology auto` (the default) the technology preset is chosen from the **shallower side**: `dense` only
when the reference and every query origin detect at least 5,000 genes per cell on average, `sparse` otherwise.
After depth matching both sides have about the shallower side's density, so this is the density the models see. The
pre-flight prints both densities, the gene overlap, the choice and a plain recommendation. `--technology sparse` or
`dense` is used as given; what the presets do is described in [learn](learn.md#technology-preset).

### First pass

In each draw, the reference is learned and the query annotated, as [`learn`](learn.md) and [`assign`](assign.md)
would. The flip test (see [assign](assign.md#flags)) is recorded as a flag (`--flip flag`, the default) and the calls
are kept: across datasets the accuracy of flipped calls varies widely; with `--flip refuse`, flipped cells are
refused, cannot be stable anchors and stay refused in that draw's final calls.
The consensus of the draws' first-pass calls is written to `first_pass.tsv`.

### Refinement

With `--refine yes` (the default), the first pass's stably called cells teach SCARLET what each cell type looks like
in the query itself. Per origin:

- the **stable anchors** of a cell type are the cells the first pass called uniquely to that cell type in every draw
  of the run (in quick mode: the draw's unique calls). They are found once and are the same in every draw;
- a cell type with at least 3 stable anchors is **refined**: in each draw SCARLET learns again on the stable anchors
  of the refined cell types, with their anchor labels and the draw's seed, and annotates the origin's cells with that
  model. Its **background**, the profile a cell must fit better than to be called at all, is the detection rate of
  each gene over all the origin's query cells;
- a cell type with fewer than 3 stable anchors (zero included) is **not refined**: in each draw, the cells whose
  first-pass best call is that cell type keep that draw's first-pass calls; `labels.tsv` notes "not refined: not
  anchored stably" in `scarlet_refine_note` for the cells that kept their first-pass calls in more than half the
  draws;
- an origin with fewer than 2 refined cell types is **not refined**: all its cells keep their first-pass calls in
  every draw, and `labels.tsv` notes "origin not refined: fewer than 2 cell types anchored stably";
- every stable anchor keeps its anchor label in the refined calls (flagged `scarlet_stable_anchor` in
  `first_pass.tsv` and `labels.tsv`), also in an origin that is not refined.

One round per draw. The refinement does not apply the flip test. The consensus of the draws' refined calls is
written to `labels.tsv`; the stable anchors of every cell type, and whether it was refined, are in
`transfer_record.json`. `--refine no` makes the first-pass consensus final.

### The transfer check

With `--qc yes` (the default), SCARLET learns markers on the query with the final calls (cell types with at least 3
cells; at least 2 such cell types per origin) and measures, for each cell type, the share of the reference's markers
that the query confirms (`share_conserved` in `qc.tsv`). A cell type below 0.2 is flagged `doubtful_transfer`. The
flag points to cell types whose transfer deserves a closer look; it changes no call. How reliable it is:
[Reading the output](../output.md#qctsv).

## Options

| option | default | what it does |
| :--- | :--- | :--- |
| `--origin-col` | none | the query's `obs` column naming each cell's origin |
| `--refine` | `yes` | refine on the query's own anchors |
| `--depth-match` | `yes` | match the depths of reference and query before the first pass |
| `--qc` | `yes` | the transfer check per cell type |
| `--technology` | `auto` | `auto` (from the shallower side), `sparse` or `dense` |
| `--blocklist` | `auto` | as in [learn](learn.md#blocklist) |
| `--imbalanced` | `auto` | as in [learn](learn.md#unequal-class-sizes) |
| `--gene-col` | `index` | the `var` column with the gene names; `index` uses the `var` names |
| `--draws` | `10` | how many times the transfer is repeated, each time with its own seed; `1` is quick mode |
| `--seeds` | `1,2,...` | the draws' seeds, separated by commas (by default 1 up to the number of draws) |
| `--flip` | `flag` | the flip test of the first pass: `flag` (recorded, the calls kept) or `refuse` |
| `--chunk-size` | `10000` | cells scored at a time; lower it to use less memory; it does not change the result |
| `--ncpu` | `1` | the draws run side by side, up to this many at a time; it does not change the result |
| `--report` | off | also write the [transfer report](../reports/transfer.md) when the run ends |
| `--verbosity` | `normal` | what the command prints: `normal` (the pre-flight, one progress bar, one line per draw, the result), `quiet` (nothing), `full` (every step of every draw); it does not change the result |

`--seeds` is an [advanced option](../advanced.md#seeds), as are `--n-markers`, `--marker-shortlist`,
`--n-iterations`, `--n-cells`, `--blocklist-modules` and `--alpha`; the published benchmarks hold at their defaults
only.

## Run time, memory and CPUs

Measured with SCARLET 0.6.3 on two of the published benchmarks' transfer tasks (Abdelaal et al. 2019), one fold
each, and on a large query, as a user runs it: the command line, every option at its default (including
`blocklist auto`), `--ncpu k` on k physical CPU cores (the run confined to them, one thread per core); the whole
`transfer` run, from reading the input files to writing its outputs, without the report, one run each. Peak memory
is that of all SCARLET's processes together.

**PBMC, 10x v2 → six other protocols**: 6,444 reference / 16,710 query cells, about 900 genes detected per cell,
sparse preset.

| mode | `--ncpu 1` | `--ncpu 2` | `--ncpu 5` | `--ncpu 10` |
| :--- | :--- | :--- | :--- | :--- |
| accurate (10 draws, the default) | 82 s, 0.5 GB | 51 s, 1.3 GB | 28 s, 2.3 GB | 25 s, 3.7 GB |
| quick (`--draws 1`) | 16 s, 0.5 GB | | | |

**A large query: the same PBMC 10x v2 reference → all 65,943 cells of the Zheng 68K dataset**: about 550 genes
detected per query cell, sparse preset (not a benchmark task).

| mode | `--ncpu 1` | `--ncpu 2` | `--ncpu 5` | `--ncpu 10` |
| :--- | :--- | :--- | :--- | :--- |
| accurate (10 draws, the default) | 42 s, 1.9 GB | 30 s, 2.4 GB | 22 s, 3.9 GB | 21 s, 5.9 GB |
| quick (`--draws 1`) | 10 s, 1.9 GB | | | |

**Mouse brain, one region → the other**: 8,128 reference / 12,552 query cells, about 11,000 genes detected per
cell, dense preset.

| mode | `--ncpu 1` | `--ncpu 2` | `--ncpu 5` | `--ncpu 10` |
| :--- | :--- | :--- | :--- | :--- |
| accurate (10 draws, the default) | 3.7 min, 3.1 GB | 2.2 min, 6.0 GB | 71 s, 10.4 GB | 52 s, 17.7 GB |
| quick (`--draws 1`) | 41 s, 3.0 GB | | | |

- **Deep data need more memory.** The brain task has fewer cells than the PBMC task but detects about 12 times as
  many genes per cell, and needed about 6 times the memory: a deep dataset of brain size needs about 3 GB with one
  CPU.
- **Accurate mode.** With `--ncpu k` the draws run side by side, up to k at a time; they share the input data and
  each adds its own working memory: the run gets faster and needs more memory (with 5 CPUs: brain 3.1 times faster
  with 3.3 times the memory, PBMC 3.2 times faster with 4.1 times the memory). The steps before and after the draws
  run once, with one CPU: the less time the draws take next to them, the less extra CPUs help (the large query: 1.9
  times faster with 10 CPUs). The second CPU added at most the 1-CPU memory plus 0.5 GB, each further CPU at most
  half the 1-CPU memory plus 1 GB.
- **Quick mode** takes a fifth to a quarter of accurate mode's time with one CPU. With one draw,
  `--ncpu` only runs steps inside it side by side: little or no faster, and about the same memory as with one CPU.
- **The report is a separate step** (`scarlet report`, or `--report`), with its own time and memory: with one CPU
  49 s and 0.9 GB on the PBMC transfer, 70 s and 1.2 GB on the large query, 1.7 min and 3.7 GB on the brain one,
  more than the brain transfer itself; `--ncpu` changes it little.
- **The result does not depend on `--ncpu`:** the same inputs, options and seeds give the same output, byte for
  byte.

### Choosing `--ncpu`

1. **Choose a number of CPUs that divides the number of draws:** with the default 10 draws, 1, 2, 5 or 10. The draws
   run in rounds of up to `--ncpu` at a time, so the run takes about as long as its number of rounds: 5 CPUs run the
   10 draws in two rounds, and so do 6 to 9 CPUs, so values between 5 and 10 bring no further speed-up and only add
   memory (3 or 4 CPUs, in four or three rounds, fall between 2 and 5).
2. **Memory decides how high to go.** Take the 1-CPU peak memory M of the table closer to your data, by genes
   detected per cell first, then by the number of cells. The second CPU added at most M + 0.5 GB and each further
   CPU at most M ÷ 2 + 1 GB (measured up to 10 CPUs), so k CPUs (k ≥ 2) need at most
   2 × M + 0.5 GB + (k − 2) × (M ÷ 2 + 1 GB); take the largest of 1, 2, 5 and 10 that fits in the free memory.
   For a deep dataset of brain size (M about 3.1 GB): 2 CPUs with 6.7 GB free, 5 with 14.4 GB, 10 with 27 GB
   (measured: 6.0, 10.4 and 17.7 GB). For PBMC size (M about 0.5 GB): 2 CPUs with 1.5 GB, 5 with 5.3 GB, 10 with
   11.5 GB (measured: 1.3, 2.3 and 3.7 GB).
3. Quick mode needs about the 1-CPU memory whatever `--ncpu`; extra CPUs save little or no time.
4. If you also want the report, it needs its own memory (above: up to 4.2 GB on the brain transfer).

If a run stops for lack of memory, lower `--ncpu` first. In a Python script, call `transfer` with `ncpu` above 1
from code under `if __name__ == "__main__":` (the draws run in new Python processes, which import your script
again); in a notebook nothing is needed.

## Hard transfers

Some transfers are hard, and SCARLET then refuses more cells rather than guess: a query of another technology than
the reference (for example a plate-based query with a droplet reference, or the reverse), a reference and a query
of very different sizes (a reference with a few hundred cells per cell type), and closely related populations
(subtypes that differ by few genes). The refused cells are reported with their reason (background, ambiguity,
undetermined across draws), and the report shows which cell types they come from.

## What transfer writes

Everything goes into the output folder:

| file | what it is |
| :--- | :--- |
| `first_pass.tsv` | the consensus of the draws' first-pass calls (the reference model's calls on the query) |
| `labels.tsv` | the final calls: the consensus of the draws' refined calls |
| `qc.tsv` | with `--qc yes`: the transfer check per origin and cell type |
| `transfer_record.json` | the full record of the run: options, the settings resolved, the pre-flight, the draws and their seeds, depth matching per draw and origin, the refinement per origin, each draw's result, the check, the summary, versions |
| `summary.txt` | the run's [summary](../output.md#the-run-summary) |
| `reference.scarlet` (or `reference_<origin>.scarlet`) | the first draw's reference model(s), usable with `assign` |
| `refine_<origin>.scarlet` | the first draw's refinement model of each refined origin |
| `report/` | with `--report`: the transfer report |

Columns: [Reading the output](../output.md#transfer). In Python, `transfer` returns a result with
`reference_models`, `query_models` (the first draw's refinement models, by origin), `first_pass`, `labels`, `qc` and
`record`
(the content of `transfer_record.json`). A refined model can annotate another query:

```python
import scarlet

if __name__ == "__main__":
    res = scarlet.transfer("reference.h5ad", "cell_type", "query.h5ad", "transfer/", origin_col="batch")
    scarlet.assign(res.query_models["batch1"], "query2.h5ad", "query2_annotated.h5ad")
```

## What can go wrong

- A label or origin column that is not in `obs`, negative values, little gene overlap: see
  [Troubleshooting](../troubleshooting.md).
- **Query cell types the reference does not have.** SCARLET calls only the cell types it learned. A cell of another
  type is refused when the background explains it at least as well; otherwise it receives the most likely learned
  cell type. The refinement can also turn cells refused in the first pass into calls; the report shows these cells,
  and the first-pass calls stay in `first_pass.tsv`.

See also: [Reading the output](../output.md), [Advanced options](../advanced.md),
[Troubleshooting](../troubleshooting.md).
