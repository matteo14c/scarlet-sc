# Advanced options

Seven options that most users never need to change. Each has a default chosen for the method, and the benchmarks
were run at those defaults.

**The [benchmarks](https://github.com/matteo14c/scarlet-sc/blob/v0.6.3/benchmark/index.md) used the defaults of every option in this section, and ran SCARLET without a
gene blocklist (`blocklist none`), so that all tools saw the same genes; in normal use keep the
[blocklist](blocklist.md) on (`blocklist auto`, the default).**

With any other value of an option in this section, the published benchmarks do not apply to your run.

Every run's [summary](output.md#the-run-summary) and its [report](reports/index.md) list every option set to a
value other than its default (for example "Options not at their defaults: alpha 0.05 (default: 0.01)."), and the
report adds that the benchmarks used the defaults of these advanced options.

| option | command line | used by | default |
| :--- | :--- | :--- | :--- |
| [`n_markers`](#n_markers) | `--n-markers` | learn, transfer | the technology preset's: 50 |
| [`marker_shortlist`](#marker_shortlist) | `--marker-shortlist` | learn, transfer | the technology preset's: 150 |
| [`n_iterations`](#n_iterations) | `--n-iterations` | learn, transfer | 100 |
| [`n_cells`](#n_cells) | `--n-cells` | learn, transfer | 100 |
| [`blocklist_modules`](#blocklist_modules) | `--blocklist-modules` | learn, transfer | `auto` |
| [`alpha`](#alpha) | `--alpha` | assign, transfer, report | 0.01 |
| [`seeds`](#seeds) | `--seeds` | transfer | 1, 2, ... up to `--draws` |

## n_markers

**Markers per cell type (default 50).** SCARLET describes each cell type by a short list of marker genes, the genes
that are most consistently detected in that type and rarely elsewhere, and annotates a cell by comparing its
detected genes with every type's markers. `n_markers` is the length of each type's list; the default, 50, is the
same for sparse and dense data. Raising it adds genes that are less specific to the type: calls then rest on more
genes, but the extra genes carry less information each, and models and annotation get larger and slower. Lowering
it keeps only the most specific genes: annotation is faster, but each call rests on fewer genes, so a cell whose few
markers went undetected is more easily left unassigned or miscalled. You might raise it when many closely related
cell types differ only by small amounts in many genes, or lower it when the data have few genes.

## marker_shortlist

**The shortlist size used to rank markers (default 150).** Learning repeats a marker search many times, each time
on a fresh random draw of cells (see [`n_iterations`](#n_iterations)). In each round the genes are ranked for each
cell type, and the top `marker_shortlist` form that round's shortlist; a gene's stability is the number of rounds in
which it made the shortlist, and the markers are the [`n_markers`](#n_markers) genes with the highest stability.
Raising `marker_shortlist` makes the shortlist easier to enter, so more genes reach full stability and the final
ranking among them depends more on how strongly they score; lowering it towards `n_markers` keeps only genes that
are near the top in almost every round. Keep it larger than `n_markers`. Few users need to change it.

## n_iterations

**Rounds of the marker search (default 100).** Each round draws a fresh random sample of cells from each cell type
and finds that sample's best genes; repeating it shows which genes are good markers consistently rather than by
chance. More rounds make the marker lists more reproducible from one seed to another; learning time grows roughly in
proportion, and a transfer, which learns several times (on the reference, then on the query for the refinement and
for the transfer check), grows with it. Fewer rounds are faster, but the markers depend more on the seed. You might
lower it for a quick first look at a large reference; use the default for any result you intend to report.

## n_cells

**Cells drawn per cell type in each round (default 100).** In every round each cell type contributes the same
number of cells, drawn at random with replacement, so that a large cell type does not swamp a small one; a type
with fewer cells than `n_cells` contributes as many cells as it has, again drawn with replacement. Raising it makes
each round's view of a type closer to the whole type, which captures more of a varied type's diversity, but rounds
become more alike, so stability says less, and learning takes longer. Lowering it makes rounds noisier and faster,
and genes detected in only part of a type become less stable. You might raise it when your cell types are large and
varied.

## blocklist_modules

**Which gene groups of the blocklist to apply (default `auto`).** The bundled human and mouse blocklists are
divided into six groups: `core_stress`, `core_ribosomal`, `core_mitochondrial`, `core_sex` (stress-response,
ribosomal, mitochondrial and sex-linked genes), `cell_cycle` and `stress_extended` (a wider set of stress genes).
The genes in an applied group cannot become markers. With `auto`, SCARLET applies the four `core_` groups, leaves
`stress_extended` out, and applies `cell_cycle` unless one of your labels names a dividing state (such as cycling,
proliferating, dividing, stem cell, precursor, or some progenitor and blast stages), so that cell-cycle genes can
mark such cells. A few non-dividing lineages whose names contain such words (oligodendrocyte precursor cells, for
example) do not count, unless the label also says cycling, proliferating, dividing or mitotic. `all` applies all
six groups; a comma-separated list applies only the groups named. Applying more groups keeps more technical and
cell-state genes out of the markers, so calls rely less on dissociation stress or the cell cycle, but a cell type
genuinely defined by such genes loses them; applying fewer groups allows them. The option only acts on the bundled
lists: with `blocklist none`, or your own file without these group names, it has no effect. You might choose `all`
when dissociation stress varies between your samples, or name the groups yourself when your labels describe dividing
cells in words the rule does not recognise.

Every run's summary and its report state whether the cell-cycle genes were blocked and, if they were not, which
label or labels kept them usable (for example "Cell-cycle genes: not blocked, because the label 'cycling T cells'
names dividing cells."). In a transfer, each of its learns decides on its own labels; the summary and the record also
name any learn of the transfer that decided otherwise than the reference's.

## alpha

**The threshold for a confident ("unique") call (default 0.01).** For each cell, SCARLET tests whether its best cell
type fits clearly better than the strongest alternative, the second-best type or the background; the p-values are
adjusted for the number of cells tested (Benjamini–Hochberg), giving `scarlet_fdr`, and when `scarlet_fdr` is below
`alpha` the cell gets a unique call (otherwise `scarlet_unique` is `unassigned`). `alpha` also decides which other
types are listed in `scarlet_all` (those not clearly worse than the best). It does not change the most likely type
(`scarlet_top1`) or `scarlet_best`. Raising it (for example to 0.05) gives more unique calls, with more errors among
them; lowering it (for example to 0.001) gives fewer, each more reliable. In a transfer, `alpha` also selects the
query cells that the refinement learns from (the [stable anchors](glossary.md#stable-anchors): unique calls in every
draw), so it changes the final calls as well. You might raise it when you need calls for more cells and can accept more errors, or lower it when every call
must be reliable.

`assign` records the `alpha` it used in the annotated file (`uns['scarlet_options']`), and the report reads it from
there. `report --alpha` only sets the threshold of the learn report's self-assessment.

## seeds

**The seeds of a transfer's draws (default 1, 2, ... up to the number of draws).** A transfer repeats its random
steps, depth matching and the marker search, in [draws](guide/transfer.md#accurate-and-quick-mode) (10 by default; 1
in quick mode), each with its own seed, and gives each cell the outcome at least 8 of the 10 draws agree on.

- `seeds`: the seeds themselves, one per draw, all different (on the command line separated by commas, for example
  `--seeds 11,12,13,14,15,16,17,18,19,20`). Given, their number sets `draws`. The same inputs, options and seeds give
  the same output, byte for byte, whatever `--ncpu`; another seed list gives another sample of the same random steps
  and comparable results ([how comparable](guide/transfer.md#reproducibility)).
- In quick mode (`--draws 1`) the one seed decides the result: run again with `--seeds 2` (or any other) to see how
  much your cells depend on it.

You might set `seeds` to check how much a result depends on the seeds.

## Unequal class sizes (imbalanced)

`imbalanced` (`--imbalanced`, learn and transfer) is a plain choice, not an advanced option: `auto` (the default),
`yes` or `no`. It decides how SCARLET handles a reference where one cell type holds a large share of the cells.

- **What is measured.** Each learn computes the largest label's share of the cells it learns from. In a transfer,
  every learn decides for itself on its own labels: the learn on the reference, the refinement's learns on the
  query's cells and the transfer check's learns.
- **The choice.** `auto` chooses **balanced** when the largest label holds more than a third of the cells, and
  **standard** otherwise; `yes` always chooses balanced, `no` always standard.
- **What balanced changes.** Two things.
  1. *The marker filter.* A gene can score for a cell type only if its detection rate across all cells is at most a
     fixed fraction of its rate in that type: one third under standard. A gene detected only in a type that holds
     more than a third of the cells can therefore never pass (its rate across all cells is that type's share of its
     rate in the type), which is why the trigger is one third. Under balanced, the allowed fraction becomes the
     largest label's share plus 5 % of the remainder, so such genes can pass.
  2. *The classifier's background.* Under balanced, the background weights every cell type equally, instead of
     following the natural mix of cells.

  The model's marker list and its background file keep the natural mix of cells.
- **These are fixed rules, not options.** The one-third trigger and the 5 % allowance are part of the method; no
  option changes them.
- **Not the same as `alpha`.** The marker filter's fixed fraction (one third under standard) is a rule of the marker
  search; [`alpha`](#alpha) is the threshold of a unique call when cells are annotated. The two are unrelated.
- **What is recorded.** The pre-flight prints the choice and the largest share; the model records it as "standard"
  or "balanced for unequal class sizes"; the run's summary says whether a correction was made ("Unequal identity
  sizes: corrected (the largest identity holds 45% of the cells)." or "no correction needed").
