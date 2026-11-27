# Glossary

The terms the documentation uses, in alphabetical order.

### Agreement

In a transfer, how many [draws](#draw) gave a cell's most frequent outcome ("unassigned" counts as an outcome):
`scarlet_best_agreement` and `scarlet_unique_agreement` (for example 9 of 10). Not written in [quick
mode](#quick-mode).

### Anchors

See [stable anchors](#stable-anchors).

### Best call

`scarlet_best`: the most likely cell type, or `unassigned` when the background explains the cell at least as well
(or the flip test refuses it). See [Reading the output](output.md#the-calls).

### Boundary cell

A cell with a [best call](#best-call) but no [unique call](#unique-call): its best cell type is not separated from
the [runner-up](#runner-up). The report's boundary map counts them by best call × runner-up.

### Blocklist

Genes kept out of the markers before the marker search: technical and broadly expressed genes (stress, ribosomal,
mitochondrial, sex-linked and cell-cycle genes). SCARLET ships a human and a mouse list, divided into
[gene groups](#gene-groups); `auto` chooses between them from the gene names. See [learn](guide/learn.md#blocklist).

### Consensus

In a transfer, each cell's final call: the most frequent outcome across the [draws](#draw), kept when at least 8 of
the 10 draws give it (80 % of the draws, rounded up, for another number of draws); otherwise the cell is unassigned,
[undetermined across draws](#undetermined-across-draws). A fixed rule: by chance alone, 8 or more of 10 draws agree
in about 5.5 % of cases.

### Depth matching

In a transfer, the deeper of reference and query is thinned towards the shallower one before the first pass, so
that both detect comparable numbers of genes. See [transfer](guide/transfer.md#depth-matching).

### Draw

One repetition of a transfer (depth matching, the first pass and the refinement) with its own seed. A transfer runs
10 draws by default (`--draws`, `--seeds`); its calls are the [consensus](#consensus) of the draws. One draw is
[quick mode](#quick-mode).

### Doubtful transfer

The [transfer check](#transfer-check)'s flag on a cell type whose reference markers the query confirms poorly (less
than a fifth of them). It changes no call. See [Reading the output](output.md#qctsv).

### First pass

In a transfer, the query annotated by the model learned on the reference (`first_pass.tsv`), before the
[refinement](#refined-call).

### Flip test

Every cell is also scored without the floor that keeps detection rates away from 0 and 1. When the most likely cell
type differs between the two scorings, the call depends on genes the candidate cell type never showed in the
reference: the cell "flips". `assign` refuses such cells by default; a transfer's first pass only flags them
(`scarlet_flip`) by default. See [assign](guide/assign.md#flags).

### Forced call

`scarlet_top1`: the most likely cell type for every cell, whatever the background and the test say. It never
refuses.

### Gene groups

The six groups the bundled blocklists are divided into (`core_stress`, `core_ribosomal`, `core_mitochondrial`,
`core_sex`, `cell_cycle`, `stress_extended`). A fixed rule on the labels decides which apply; the cell-cycle group
is left out when a label names dividing cells. See [Advanced options](advanced.md#blocklist_modules).

### Identity

A cell type of the reference as SCARLET learns it: one per label. The columns and tables SCARLET writes call cell
types identities (for example `identity` in `qc.tsv`).

### Marker

A gene that a cell type's cells consistently detect and that is rarely detected elsewhere. Each cell type has 50
markers by default; assign scores cells on them. See [learn](guide/learn.md#how-learn-works).

### Marker shortlist

The best-ranked genes of a cell type in one round of the marker search (150 by default). A gene's stability is the
number of rounds in which it made the shortlist; the markers are the most stable genes. See
[Advanced options](advanced.md#marker_shortlist).

### Model file

The one file `learn` writes (`.scarlet`): the markers, the detection rates assign scores, and the run's settings,
provenance and summary. See [Reading the output](output.md#the-model-file).

### Origin

Where a query cell comes from (a batch, a sample, a technology), named by a column of the query given with
`--origin-col`. A transfer matches depths, refines and checks per origin. See
[transfer](guide/transfer.md#origins).

### Positive evidence

A flag on a call whose cell detects fewer of the cell type's markers than a fifth (rounded down) of the number a
cell of that type is expected to detect (`scarlet_pe_low`). It changes no call. See
[assign](guide/assign.md#flags).

### Quick mode

A transfer run once (`--draws 1`): the same method, one [draw](#draw). Its calls are that draw's: no
[consensus](#consensus), no [agreement](#agreement), and the calls depend on the seed. About a fifth of the
time of the default (10 draws) with one CPU (measured: 16 s against 82 s on a PBMC transfer, 41 s against 3.7 min on a
deep mouse-brain one). See [transfer](guide/transfer.md#accurate-and-quick-mode).

### Refined call

In a transfer, the final call (`labels.tsv`): the consensus of the [draws](#draw) in which the query is annotated
again by models learned on the query's own [stable anchors](#stable-anchors). A cell type with fewer than 3 stable
anchors is not refined and its cells keep their first-pass calls; so do all the cells of an origin with fewer than 2
refined cell types.

### Refused cell

A cell left `unassigned` in its best call, its unique call or both, with a reason: background, ambiguity, flip test or (transfer) undetermined across draws.
See [Reading the output](output.md#the-calls).

### Runner-up

The cell type ranked next after a cell's best call: in a transfer the most frequent across the [draws](#draw)
(`scarlet_runner_up`, with `scarlet_runner_up_draws`); in the learn report's self-assessment the second-ranked cell
type.

### Self-assessment

In the learn report, the model's own reference annotated by cross-validation: each fold by a model learned on the
other folds. It shows how well the model separates its own cell types. See
[report](reports/learn.md#the-self-assessment).

### Stable anchors

In a [transfer](guide/transfer.md#refinement), the query cells the first pass called uniquely to the same cell type in
every [draw](#draw) of the run (in [quick mode](#quick-mode), the draw's unique calls). The refinement learns each
cell type again from its stable anchors, and each stable anchor keeps its label; a cell type with fewer than 3
stable anchors is not refined.

### Technology preset

The marker search's setting for the data's density: `sparse` (droplet data such as 10x) or `dense` (plate data such
as Smart-seq2). Both read a gene as detected when its count is above 0; they differ in the background penalty of the
marker search (6 sparse, 0 dense). `auto` chooses from the data. See [learn](guide/learn.md#technology-preset).

### Transfer check

In a transfer, the share of each cell type's reference markers that markers learned on the query confirm
(`qc.tsv`); below 0.2 the cell type is flagged [doubtful transfer](#doubtful-transfer). See
[transfer](guide/transfer.md#the-transfer-check).

### Unassigned

The value of a refused call.

### Undetermined across draws

The reason of a transfer's unassigned cell when no outcome won clearly across the [draws](#draw): the most frequent
outcome was given by fewer than 8 of the 10 draws (the [consensus](#consensus) rule).

### Unequal class sizes

How SCARLET handles a reference where one cell type holds a large share of the cells (`imbalanced`): **standard**,
or **balanced for unequal class sizes** when the largest cell type holds more than a third of the cells. See
[Advanced options](advanced.md#unequal-class-sizes-imbalanced).

### Unique call

`scarlet_unique`: the best cell type only when it is clearly better than its competitor (`scarlet_fdr` below
`alpha`, 0.01), else `unassigned`. See [Reading the output](output.md#the-calls).
