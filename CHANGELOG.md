# Changelog

User-visible changes of SCARLET. The public package holds the commands `learn`, `assign`, `transfer` and `report`.
The documentation is in `docs/`.

## 0.6.3 (2026-11-27)

0.6.3 makes the public package exactly the four documented commands, with plain options, and makes every run say
what it did.

- **Python 3.11 to 3.13** (Python 3.10 is no longer supported). With Python 3.12 or later SCARLET requires scanpy
  1.12.6 or later; with Python 3.11, umap-learn below 0.6 and pynndescent below 0.7 (scanpy up to 1.11 does not work
  with their newer releases in the report's UMAP).

- **Dependency floors** are the oldest versions the tests ran with: anndata 0.12.14, scanpy 1.11.5, numpy 2.4.4,
  pandas 2.3.3, scipy 1.17.1, scikit-learn 1.8.0.

### From the beta testers' feedback

- **The benchmark is its own section of the repository** (`benchmark/`, linked from the README and the documentation),
  no longer part of the documentation site.
- **Reports that never stop a run.** With `--report yes`, `learn` and `transfer` keep their outputs and end normally
  when the report cannot be written (plotly missing, or an error in the report): a warning says so. The learn report's
  self-assessment reads the genes from the model's gene column (`--gene-col`), as learn did; repeated gene names no
  longer stop a report.
- **A fixed order among tied marker genes.** When genes tie exactly at the marker search's cut-off, the more
  specific gene (the lower score in every other cell type) is taken first, then the gene name; before, the choice
  depended on the CPU (numpy's vector code). Calls change in a few cells; the benchmark numbers are recomputed.

- **Less output.** `verbosity` (`--verbosity`) on every command: `normal` (the default) prints the pre-flight once,
  one progress bar, one line per draw or fold and the result; `quiet` prints nothing; `full` prints every step, as
  before. It changes no result and is not recorded as a setting.
- **Example data from the package.** `scarlet.datasets.fetch` (and `python -m scarlet.datasets`) downloads the
  tutorials' data and checks their SHA-256; while the repository is private it falls back to the GitHub CLI. The data
  are now the inputs of SCARLET's benchmark, with the benchmark's own selection of cells. Each file shows its progress
  while it downloads.
- **Install extras.** `report` (the HTML reports), `tutorials` (the reports, JupyterLab and the widgets the notebooks
  need) and `docs` (the pinned packages to build the documentation locally).
- **A clear error for scripts.** `transfer` with `ncpu` above 1, called from a Python script without
  `if __name__ == "__main__":`, stops with one line saying so instead of a broken process pool, and the script exits
  with an error; also when the script calls `learn` or `assign` first (before, the script could then wait forever).
- **Input checksums.** `learn` records the SHA-256 of its reference file in the model, `transfer` of both its inputs
  in `transfer_record.json`; the reports show them and check them against the files they read.
- **Reports.** Each page opens on a summary of the result and ends with the provenance. The UMAPs have their own
  switches (output, view, highlight) and are computed by default on the model's own genes, with all genes as a second
  space (`umap_all_genes`). `query_labels` (`--query-labels`) compares the calls of a transfer or an assign with
  labels you know, in a section marked as yours. Readable table headers and numbers; per-cell tables with a filter
  (the cells worth a look first) and sortable columns; agreement histograms; the depth before and after matching;
  the drill-down heatmaps in one gene order; the transfer named by its datasets; folds numbered 1 to 5; certainty
  relative to the run; without WebGL, the UMAPs are drawn as SVG. In the transfer report: the drill-down in three
  parts (markers common to both models, only in A's, only in B's) with each model's own estimate for every one of
  them, read from the models' in-class tables, and each marker's rank in both models; *what refinement changed* on a
  page of its own (`refinement.html`); the per-cell data counted first, the table folded and shown 25 rows at a time.
  In the transfer and assign reports, the bar chart of the calls has no numbers above the bars (in the hover instead).
- **The model keeps learn's in-class table** (`in_class.tsv`): for every gene the marker search scored, each cell
  type's share of cells detecting it, the same bootstrap estimate the universe holds for the markers. The transfer
  report reads it. A mouse-brain model grows from about 0.5 to 5 MB, a PBMC model from 0.2 to 0.9 MB; calls do not
  change. Earlier versions of SCARLET read such models (they ignore the new member).
- **Tutorials.** The notebooks' cells carry ids (nbformat 4.5); the tutorials page says what notebook 03's memory
  numbers measure and on how many CPUs, and its time on a laptop. The data download's progress line ends at 100 %
  also through the GitHub CLI.
- **Troubleshooting** lists the messages in tables: message, cause, fix.
- **Documentation.** New pages: tutorials, the blocklist, the reports (one per kind), the benchmarks figure by figure
  (with run time and memory, and every input's SHA-256), known limitations.

### Transfer: 10 draws, the consensus, quick mode

- **Draws.** A transfer repeats depth matching, the first pass and the refinement in 10 draws by default (`--draws`),
  each with its own seed (`--seeds`, by default 1 up to the number of draws) for depth matching and the marker
  search. Each cell's best and unique calls are the outcome most draws give ("unassigned" counts as one) when at
  least 8 of the 10 draws give it, otherwise `unassigned` with the reason "undetermined across draws". This is a
  fixed rule (by chance alone, 8 or more of 10 agree with one-sided p 0.055); for another number of draws a call
  needs 80 % of them, rounded up. The same inputs, options and seeds give the same output, byte for byte, whatever
  `--ncpu`.
- **Quick mode.** `--draws 1` runs the same method once: the draw's own calls, no consensus and no per-cell
  agreement; a fifth to a quarter of the time with one CPU (measured).
- **CPUs.** With `--ncpu k` the draws run side by side, up to k at a time, sharing the input data: faster, with
  more memory (measured on two benchmark tasks with 5 CPUs: 3.1–3.2 times faster, 3.3–4.1 times the memory). The
  default is 1. Choose a
  number that divides the number of draws (with 10 draws: 1, 2, 5 or 10); values between 5 and 10 bring no further
  speed-up. The transfer guide gives the measured run time and peak memory per mode and CPU count, and how to
  choose `--ncpu` by the available memory: a deep dataset of mouse-brain size needs about 3 GB with one CPU. The
  marker search adds its bootstrap results in the same order whatever `--ncpu`, so the result does not depend on
  it, bit for bit.
- **Memory.** The same outputs, byte for byte, with less memory: input counts stored as 64-bit integers are held as
  32-bit integers when every value fits; the per-cell summaries read the data in bounded chunks; detection patterns
  share one read-only array of their (all true) values; learn keeps only the detection of an input it read; nothing
  of a finished step is kept for the next. Measured with one CPU (through the Python API, without the blocklist): a
  deep mouse-brain transfer 4.7 → 3.0 GB, a learn on 44,000 cells 1.7 → 1.0 GB, the brain transfer's report
  6.2 → 3.6 GB (its UMAP works in 32-bit floats, as scanpy does, and its self-assessment on the detection pattern).
  Each marker search frees its data when it ends: with the default blocklist its data were kept until Python's
  garbage collector ran, so a transfer's twenty or so marker searches piled them up (the deep mouse-brain transfer
  as a user runs it, one CPU: 5.8 → 3.1 GB).
- **Stable anchors.** The refinement learns each cell type from its stable anchors: the cells the first pass called
  uniquely to that cell type in every draw of the run. A cell type is refined with at least 3 stable anchors; an
  origin needs at least 2 refined cell types. Every stable anchor keeps its anchor label.
- **The refinement's background** is the detection rate of each gene over all the origin's query cells.
- **Options:** `draws` and `seeds` (advanced) replace transfer's `seed`. `learn` and `assign` keep `seed` and run
  once.
- **Outputs:** `first_pass.tsv` and `labels.tsv` hold the consensus calls, with each call's agreement (how many
  draws gave the most frequent outcome; not in quick mode) and, for a refused cell, its reason (`background`,
  `ambiguity`, `undetermined across draws`, `flip test`); `scarlet_top1` and `scarlet_runner_up` with
  `scarlet_runner_up_draws`; `scarlet_stable_anchor` (a flag, in both files); in `labels.tsv` the refine note "not
  refined: not anchored stably" or "origin not refined: fewer than 2 cell types anchored stably". The flip and
  positive-evidence flags are set when raised in at least half the draws, with the draw counts. The models saved are
  the first draw's. The record (`transfer_record.json`, schema 4) holds the draws, their seeds, the mode, the
  consensus rule, depth matching per draw, the stable anchors and the backgrounds, and each draw's result.

### Detection only, the presets and the flip test

- **No binarisation.** Both technology presets read a gene as detected when its count is above 0, in learn, assign
  and transfer. The dense preset no longer reduces each cell to a target number of detected genes; it differs from
  the sparse preset only in the marker search's background penalty (0 dense, 6 sparse). The model's settings record
  (schema 3) names the detection rule.
- **The flip test** refuses a flipped call in `assign` with every model (`--flip refuse`, the default; `--flip
  flag` keeps it); `auto` is removed. `transfer` flags flipped calls by default (`--flip refuse` is available); its
  refinement does not apply the test.

### The report

- **Transfer report**: the query UMAP has two switches, the output and the view (label; certainty: intensity =
  agreement k of N, full colour in quick mode; stable anchors outlined; unlabelled cells by reason); hover shows the
  call, agreement, runner-up, reason and anchor flag. Boundary maps (best call × runner-up, first pass and refined),
  a per-identity table (calls, stable anchors, refined or not and why, % full agreement, % boundary, % unlabelled by
  reason) and a per-cell table downloadable from the page.
- **Learn report**: the same certainty and boundary views on the reference, from the self-assessment (FDR, margin,
  runner-up), and a downloadable per-cell table.
- **All reports**: a colour-blind-safe palette (Okabe-Ito, Paul Tol's muted and light) for up to 24 cell types.
- **Fixes**: `transfer(report=True)` with AnnData inputs writes its report; `report=True` and `report` check for
  plotly (the `report` extra) before any work.
- **Removed:** the root files `API.md` and `DEFAULTS.md`, left from earlier versions; the documentation in `docs/`
  covers what they described of the four commands.

### The commands and their options

- **Four commands.** The package holds `learn`, `assign`, `transfer` and `report`, plus `scarlet.load_model` and
  `Model`. The experimental commands are removed from the public package; tools to merge, split and discover cell
  identities will be added in later releases. With them go two dependencies used only for clustering,
  python-igraph and leidenalg, which are no longer installed.
- **Fewer options.** Every option that existed only for development is removed, with no alias: the profile option
  and the earlier-behaviour profile; from `learn` and `transfer` the explicit blocklist categories, the external
  universe, the bootstrap order, active-gene normalisation and its targets, the minimum genes per cell, the marker
  filter's settings, the classifier background mode, the universe mode and temper, the marker score's beta, gamma
  and eps, the positive-evidence switch and the raw marker-search arguments; from `assign` `eps`,
  `save_probs_to_obs`, `scoring` and `positive_evidence`; from `transfer` `depth_match_seed`, `eps` and `scoring`;
  from `report` `folds`. The method uses, fixed, the values every default run used.
- **`marker_shortlist` replaces `n_background`** (command line `--marker-shortlist`): the shortlist size used to
  rank markers. Same default (the technology preset's, 150), same behaviour.
- **Technology preset:** `auto` (the default), `sparse` or `dense`. The former names are not accepted.
- **`transfer`:** `refine` and `depth_match` take yes / no (`True` / `False`): `refine yes` is the refinement where
  anchored, `depth_match yes` the depth matching by capture, as in 0.6.2's defaults; the alternative refinement and
  depth-matching methods are removed. `flip` takes `flag` (the default) or `refuse`.
- **`assign`:** `flip` takes `refuse` (the default) or `flag`.
- **Advanced options.** `n_markers`, `marker_shortlist`, `n_iterations`, `n_cells`, `blocklist_modules` and `alpha`
  keep their defaults and are documented as advanced options (the documentation's "Advanced options" page): the published benchmarks hold at
  their defaults only.

### The blocklist

- **`blocklist`: `auto` (the default), `human`, `mouse`, `none` or a file.** `auto` chooses the bundled human or
  mouse list from the gene names, says how in the pre-flight, the summary and the report, and stops with an error
  when it cannot tell the species. In 0.6.2 the default was the human list for all data.
- On human data the calls at the defaults are the same as 0.6.2's (checked cell for cell on the benchmark's
  datasets). On mouse data they can differ, since the mouse list is now applied.
- **A blocklist path that does not exist is an error** (in 0.6.2: a warning, and no blocklist).
- The earlier lists that were shipped next to the current ones are no longer in the package. A list's version is
  recorded only in the run's provenance.

### What a run says and records

- **Run summary.** Every run writes a few lines in plain words: what it did, on which data, the technology preset
  and blocklist and why, the marker sizes, the class-size handling, the options not at their defaults, the outputs
  and the versions. `learn`: `<model stem>.summary.txt` next to the model and the model member `summary.txt`;
  `assign`: `<output stem>.summary.txt` and `uns['scarlet_summary']`; `transfer`: `summary.txt` in the output folder
  and in `transfer_record.json`. The report shows it under "What this run did".
- **Options not at their defaults.** Every run records them (model `settings.json`, the annotated file's
  `uns['scarlet_options']`, `transfer_record.json`); the summary and the report state them, and the report adds,
  for an advanced option, that the published benchmarks used the defaults.
- **Cell-cycle genes.** The summary and the report state whether the cell-cycle genes were blocked and, if not, which
  labels kept them usable; `transfer_record.json` records the decision of each learn of a transfer.
- **What each learn read.** A learn takes its cells in memory (no temporary file). The model's
  `marker_finder_params.json` records its input: the reference file (`learn`), `in-memory AnnData` and its shape
  (`learn` from Python), or, for each learn of a transfer, the input file's path and sha256 and how the cells were
  used (thinned in draw i with depth-matching seed s, or not thinned; the refinement and the transfer check on the
  query's cells).
- **`assign` records its `alpha`** (`uns['scarlet_options']`); the report's assign view reads it, so it no longer
  needs to be given.
- **Renamed records:** `uns['scarlet_profile']` is now `uns['scarlet_assign']` (the flip mode and the positive-evidence
  flag); `uns['scarlet_scoring']` is now `uns['scarlet_floor']`; the model's settings
  record the marker sizes as `n_markers` and `marker_shortlist`, and the class-size handling as "standard" or
  "balanced for unequal class sizes". Records no longer carry a profile.
- **Older models.** `assign` refuses a model without scoring statistics (written before 0.6.0) and asks to learn it
  again; models written before 0.3.0 are read as written, without conversion.
- The pre-flight summary, the messages and the transfer check's message are in plain words.

## 0.6.2

- **New: the `report` command.** `scarlet report RUN` (Python: `scarlet.report(run)`) writes the HTML report of a
  finished run: one self-contained page plus a bundle of tables, figures and a manifest. A model file gives the
  learn report (with a cross-validated self-assessment), an annotated `.h5ad` the assign view, a transfer folder the
  transfer report. `learn --report` and `transfer --report` write it at the end of the run. It needs the `report`
  extra (plotly).
- **`learn` chooses its technology preset from the reference** (`auto`, the default): dense when the reference
  detects at least 5,000 genes per cell on average, else sparse, the rule `transfer` uses. The pre-flight prints the
  choice and the model records it.
- **`assign` with a dense model reduces the query cells** exactly as `learn` reduced the reference (the model's
  targets and its learn seed), as `transfer` already did. The annotated `.h5ad` keeps the input's values.

## Earlier versions

Versions before 0.6.3 were not public. Among the changes of the earlier versions that 0.6.3 users meet: the floor on detection rates, the
flip test and the positive-evidence flag (0.6.0); depth matching and the refinement where anchored in `transfer`
(0.6.0); 50 markers ranked on a shortlist of 150 for dense data, as for sparse data (0.6.1); the model file, the
pre-flight summary, the bundled human and mouse blocklists with their gene groups, the class-size handling and the
label-order-independent bootstrap (0.2.0); the reduced memory use of learning and annotation (0.2.1).
