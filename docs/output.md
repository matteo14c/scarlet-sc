# Reading the output

What every column, file and summary line SCARLET writes means. In the tables and columns SCARLET writes, a cell type
of the reference is called an **identity** (one per label of the reference).

## The calls

SCARLET gives each cell three calls, from the most to the least permissive:

| call | column | read it as |
| :--- | :--- | :--- |
| forced | `scarlet_top1` | the most likely cell type, whatever the background and the test say; never refuses |
| best | `scarlet_best` | the most likely cell type, or `unassigned` when the background explains the cell at least as well (the cell looks like none of the reference's cell types), or when the flip test refuses it |
| unique | `scarlet_unique` | the best cell type only when it is clearly better than its competitor (`scarlet_fdr` below `alpha`, 0.01), else `unassigned` |

`unique` leaves more cells unlabelled, and the cells it labels are more reliable; `best` labels more cells. A
refused cell has `unassigned` in `scarlet_best`, in `scarlet_unique` or in both. Its reason (in the
[report](reports/index.md)):

| reason | meaning | refuses |
| :--- | :--- | :--- |
| background | the background explains the cell at least as well as the best cell type: it looks like none of them | best and unique |
| ambiguity | a best cell type exists but is not clearly better than its competitor (the second cell type or the background) | unique only |
| flip test | the call depends on genes the candidate cell type never showed in the reference (the [flip test](glossary.md#flip-test), when it refuses) | best and unique |
| undetermined across draws | transfer only: no outcome was given by at least 8 of the 10 draws | best and unique (each on its own) |

## assign

### assign: the calls table

`assign` returns (in Python) and writes with `--predictions` (tab-separated, first column `cell`) one row per cell:

| column | content |
| :--- | :--- |
| `scarlet_best` | the best call, or `unassigned` |
| `scarlet_unique` | the unique call, or `unassigned` |
| `is_scarlet_unique` | whether the call is unique (`True` / `False`) |
| `scarlet_fdr` | the adjusted p-value of the best cell type against its competitor (Benjamini-Hochberg over the cells of the run) |
| `scarlet_flip` | whether the cell flips (the flip test) |
| `scarlet_reason` | the flip test's reason text for a flipped cell, else empty |
| `scarlet_pe_low` | the low-positive-evidence flag |
| `scarlet_pe_note` | the flag's explanation, with the numbers |

### assign: the annotated file

The annotated `.h5ad` holds the input as given (`X` unchanged), plus:

| where | content |
| :--- | :--- |
| `obs['scarlet_best']`, `obs['is_scarlet_unique']`, `obs['scarlet_fdr']` | as in the calls table (the unique call is `scarlet_best` where `is_scarlet_unique` is true) |
| `obs['scarlet_all']` | every cell type not clearly worse than the best (p-value at least `alpha`), separated by `;`; `unassigned` for a refused cell |
| `obs['scarlet_top1']`, `obs['scarlet_top2']` | the two most likely cell types; `scarlet_top1` is the forced call |
| `obs['scarlet_ll_top1']`, `obs['scarlet_ll_top2']`, `obs['scarlet_ll_bg']` | their log-likelihoods, and the background's |
| `obs['scarlet_flip']`, `obs['scarlet_reason']` | the flip test: the flag and the reason |
| `obs['scarlet_top1_j']`, `obs['scarlet_top1_eps']` | the flip test's two forced calls: with the floor (the scoring of every call) and without it |
| `obs['scarlet_pe_detected']` | how many of the forced cell type's markers the cell detects |
| `obs['scarlet_pe_markers']` | how many markers that cell type has (in the data) |
| `obs['scarlet_pe_expected']` | how many a cell of that type is expected to detect (the sum of its detection rates over its markers) |
| `obs['scarlet_pe_required']` | a fifth of the expected number, rounded down: below it the call is flagged |
| `obs['scarlet_pe_typical']` | how many a typical reference cell of that type detects (the median over its reference cells) |
| `obs['scarlet_pe_low']`, `obs['scarlet_pe_note']` | the low-positive-evidence flag and its explanation |
| `obsm['scarlet_probs']`, `uns['scarlet_identity_names']` | each cell's probability for every cell type (the likelihoods normalised over the cell types), and the columns' order (with `--save-probs yes`) |
| `uns['scarlet_assign']` | the flip mode used (`refuse` or `flag`) and the positive-evidence flag |
| `uns['scarlet_floor']` | the floor of each cell type's detection rates and of the background's (the scoring's margin from 0 and 1) |
| `uns['scarlet_options']` | the `alpha` the run used, and its options not at their defaults (a JSON text) |
| `uns['scarlet_summary']` | the run's [summary](#the-run-summary) |

Next to the annotated file, `<output stem>.summary.txt` holds the same summary.

## The model file

`learn` writes one model file; `transfer` writes one per model it learns. It is a zip archive (`.scarlet` by
convention) that `scarlet.load_model(path)` reads and checks (every member's checksum):

| member | content |
| :--- | :--- |
| `manifest.json` | the format, its version and every member's checksum and size |
| `markers.tsv` | the ranked markers per cell type: `idx`, `gene`, `stab` (stability: rounds on the shortlist), `score`, `cell_type` |
| `programmes.json` | per cell type, its marker genes in rank order |
| `universe.tsv` | what assign scores: one row per marker gene, one column per cell type (its detection rate) and `bg` (the background's) |
| `background_universe.tsv` | `gene_symbol`, `expression_probability`: the detection rate of every gene over the reference |
| `in_class.tsv` | the same table as `universe.tsv` for every gene the marker search scored (all genes after the blocklist), not only the markers: per cell type, the share of its cells detecting the gene as learn estimated it (the average over the marker search's bootstrap), and `bg`. `universe.tsv` is its marker rows. The transfer report's drill-down reads it. It is the larger part of the file: a mouse-brain model (41,000 genes, 34 cell types) about 5 MB instead of 0.5 MB, a PBMC model (22,000 genes, 9 cell types) about 0.9 MB instead of 0.2 MB. Models written before 0.6.3 do not have it |
| `marker_finder_params.json` | every setting of the marker search as used; its input (learn: the reference file, or `in-memory AnnData` and its shape; a learn of transfer: the input file's path and sha256 and how its cells were used, e.g. `in memory, thinned, draw 1, seed 1835504127 (origin X)` or `in memory, not thinned, draw 1`); the blocklist record (file, checksum, genes removed), the class-size handling and the scoring statistics (per cell type, the number of distinct cells drawn) |
| `positive_evidence.json` | per cell type: markers, the expected number detected, the flag's threshold and the typical value |
| `settings.json` | the settings the run resolved (technology preset, marker sizes, blocklist and gene groups, class-size handling) and the run's options, with its options not at their defaults |
| `provenance.json` | command, time, host name, command line, versions (SCARLET, Python and the main libraries), inputs (path, cells, genes, label column, cell counts per type) |
| `summary.txt` | the run's [summary](#the-run-summary) |

`provenance.json` records the host name, the command line and the input paths: check it before sharing a model
file. In Python: `model.identities`, `model.markers`, `model.programmes`, `model.universe`, `model.background`,
`model.in_class`, `model.settings`, `model.provenance`. From the shell: `unzip -p reference.scarlet markers.tsv`.

## transfer

### first_pass.tsv and labels.tsv

One row per query cell (first column `cell`). Both files hold the consensus of the
[draws](guide/transfer.md#accurate-and-quick-mode) (in quick mode, the one draw's calls): `first_pass.tsv` of the
reference model's calls, `labels.tsv` of the final calls, after the refinement. The columns marked "not in quick
mode" are not written with `--draws 1`:

| column | content |
| :--- | :--- |
| `origin` | the cell's origin (`all` without `--origin-col`) |
| `scarlet_best`, `scarlet_unique` | the best and unique calls (the consensus of the draws) |
| `scarlet_best_agreement`, `scarlet_unique_agreement` | how many draws gave the cell's most frequent outcome ("unassigned" counts as an outcome), for example 9 (of 10); the number of draws is in `transfer_record.json`; not in quick mode |
| `scarlet_best_reason`, `scarlet_unique_reason` | empty for a call; for `unassigned`: `background` (the background fits at least as well as the best cell type), `ambiguity` (unique only: the best cell type does not fit clearly better than the alternative), `undetermined across draws` (no outcome was given by at least 8 of the 10 draws; never in quick mode), or `flip test` (with `--flip refuse`) |
| `scarlet_flip`, `scarlet_flip_draws` | the first pass's flip test, raised in at least half the draws, and in how many draws (in both files; `scarlet_flip_draws` not in quick mode) |
| `scarlet_pe_low`, `scarlet_pe_low_draws`, `scarlet_pe_note` | the positive-evidence flag of the call in that file, raised in at least half the draws; in how many (not in quick mode); the note of the first draw that raised it |
| `scarlet_top1` | the cell type ranked first most often across the draws (in both files), whatever the call |
| `scarlet_runner_up`, `scarlet_runner_up_draws` | the cell type ranked next after the cell's best call (after `scarlet_top1` where the cell is unassigned), the most frequent across the draws, and in how many draws (in both files; empty and 0 when the model has one cell type; `scarlet_runner_up_draws` not in quick mode); for a boundary cell (a best call, no unique call) it is the cell type SCARLET could not separate it from |
| `scarlet_refine_note` | `labels.tsv` only: "not refined: not anchored stably" where the cell kept its first-pass calls in more than half the draws (its first-pass cell type was not refined), "origin not refined: fewer than 2 cell types anchored stably" for the cells of an origin that was not refined, else empty |
| `scarlet_stable_anchor` | `True` for a stable anchor (called uniquely to the same cell type in every draw; it keeps that label, so its cell type is its call), else `False`; in both files, with the refinement on |

A refused cell's reason is a note: "not refined" is not a reason, because a cell that was not refined keeps its
first-pass outcome, which carries its own reason.

### qc.tsv

The transfer check, one row per origin and cell type:

| column | content |
| :--- | :--- |
| `origin` | the origin |
| `identity` | the cell type |
| `n_cells_transferred` | query cells with that final call |
| `in_query_programme` | whether markers could be learned for it on the query (at least 3 cells, and at least 2 such cell types in the origin) |
| `reference_markers`, `query_markers` | its markers in the reference model and on the query |
| `conserved` | the reference markers found again on the query |
| `share_conserved` | `conserved` / `reference_markers` |
| `doubtful_transfer` | `True` when `share_conserved` is below 0.2; empty when it could not be measured |

**How reliable the flag is.** When the check was designed, it was measured against transfers whose true quality
was known: among the cell types it flagged, 61 % had in fact been transferred poorly (an F1 below 0.5), where 30 %
of all cell types were transferred that poorly; it caught 49 % of the poorly transferred cell types; and the share
of markers confirmed followed the true F1 with a rank correlation of 0.58. A flag is therefore a reason to look
closer, not a verdict; and an unflagged cell type can still be transferred poorly. The flag changes no call. These
values are also recorded in `transfer_record.json` (`qc_performance`).

### transfer_record.json

The full record of the run: its options and the options not at their defaults, the settings resolved, the
pre-flight (densities, gene overlap, the technology preset and why, the recommendation), the draws (their number,
seeds, the mode, the consensus rule, each draw's result and the final counts by reason), depth matching per
draw and origin (which side was thinned and how much), the refinement per origin (stable anchors per cell type,
refined or not, the origins not refined, the background of each origin and draw), whether each learn of the first
draw blocked the cell-cycle genes, the transfer check's summary and its measured performance, the run's summary and
the versions.

### The models

`reference.scarlet` (or one `reference_<origin>.scarlet` per origin when depth matching thinned the reference) and
`refine_<origin>.scarlet` for each refined origin: ordinary [model files](#the-model-file), usable with `assign`.

## The run summary

Every run writes a few lines in plain words: `<model stem>.summary.txt` next to a model (and the model's member
`summary.txt`), `<output stem>.summary.txt` next to an annotated file (and `uns['scarlet_summary']`), `summary.txt`
in a transfer folder (and in `transfer_record.json`). The report shows it under "What this run did". The lines:

| line | what it says |
| :--- | :--- |
| `SCARLET 0.6.3 learn, <date and time>` | the version, the command and when it ran |
| `What it did: ...` | the command in one sentence |
| `Reference: ...`, `Data: ...`, `Query: ...` | each input: path (or "an AnnData object in memory"), cells x genes, genes detected per cell on average, the label column and how many cell types |
| `Technology preset: ...` | the preset, what it does, and whether it was chosen automatically (with the number it was based on) or given |
| `Blocklist (genes left out of the markers): ...` | which list and why (for `auto`: the share of human- or mouse-style gene names), and the gene groups left out |
| `Cell-cycle genes: ...` | `blocked`; or `not blocked, because the label(s) ... name dividing cells`; or `not blocked (no blocklist applied)`; or `not blocked (blocklist_modules ...)` when the groups were given; or `as the given list says (it has no gene groups)` for your own file |
| `Markers: 50 per cell type, ranked on a shortlist of 150 ...` | the marker sizes, and whether they are the technology preset's or given |
| `Unequal identity sizes: ...` | whether the class-size correction was applied, and the largest cell type's share |
| `Settings: ...` | for `learn` and `assign` the seed, for `transfer` the number of draws (their seeds are in `transfer_record.json`), and, where they apply, alpha, refinement, depth matching and the flip test |
| `All options at their defaults.` or `Options not at their defaults: ...` | every option set to another value than its default, with the default; options that change no call (`ncpu`, `chunk_size`, `report`, `predictions`, `save_probs`, `verbosity`) are not listed |
| notes | for assign, the model used and, with a dense model, how the query was reduced; for transfer, the number of query cells with a unique call and any learn that decided the cell-cycle genes otherwise |
| `Output: ...` | each output written, and what it is |
| `Versions: ...` | SCARLET, Python and the main libraries |

See also: [Glossary](glossary.md), [Troubleshooting](troubleshooting.md).
