# Command-line reference

Generated from the code (tests/test_docs.py); every option is the Python parameter of the same name with
the same default ([Python reference](python.md)). yes/no options take `yes` or `no`.

## scarlet learn

Learn a model of the cell types of a labelled reference.

| option | values | default | description |
| :--- | :--- | :--- | :--- |
| `--reference` | text | required | The labelled reference: an .h5ad file of raw counts (cells x genes); in Python also an AnnData object. |
| `--labels` | text | required | The column of the reference's obs that holds the cell-type labels to learn. |
| `--model` | text | required | The model file to write (.scarlet). In Python, None returns the model without writing it. |
| `--technology` | auto, sparse, dense | auto | The technology preset: auto chooses dense when the reference detects at least 5,000 genes per cell on average, else sparse; sparse suits droplet data (e.g. 10x), dense plate data (e.g. Smart-seq2). |
| `--blocklist` | text | auto | Genes kept out of the markers (technical and stress genes): auto picks the bundled human or mouse list from the gene names; human, mouse, none, or the path of your own tab-separated file with a Symbol and/or ENSG_ID column. |
| `--blocklist-modules` | text | auto | Advanced. Which gene groups of the bundled lists to apply: auto (a fixed rule on the labels), all, or a comma-separated list of groups. See Advanced options. |
| `--imbalanced` | auto, yes, no | auto | Unequal class sizes: auto corrects when the largest cell type holds more than a third of the cells, yes always, no never. See Advanced options. |
| `--gene-col` | text | index | The var column with the gene names; index uses the var names. |
| `--n-markers` | integer | none | Advanced. Markers per cell type; by default the technology preset's (50). |
| `--marker-shortlist` | integer | none | Advanced. The shortlist the markers are ranked on: in each round of the marker search a gene counts towards its stability when it is among the cell type's top marker_shortlist genes; by default the technology preset's (150). Keep it larger than n_markers. |
| `--n-iterations` | integer | 100 | Advanced. Rounds of the marker search. |
| `--n-cells` | integer | 100 | Advanced. Cells drawn per cell type in each round of the marker search. |
| `--seed` | integer | 42 | Random seed; the same inputs, options and seed give the same model. |
| `--ncpu` | integer | 1 | Number of CPU cores to use; it does not change the result. |
| `--report` | yes, no | no | Also write the run's HTML report when the model is written (as scarlet report on it). |
| `--verbosity` | quiet, normal, full | normal | What the command prints: normal, its pre-flight once, one progress bar and one line per fold or draw; quiet, nothing; full, everything each step prints. It does not change the result. |

## scarlet assign

Annotate cells with a model.

| option | values | default | description |
| :--- | :--- | :--- | :--- |
| `--model` | text | required | The model: a .scarlet file written by learn; in Python also a Model. |
| `--data` | text | required | The cells to annotate: an .h5ad file of raw counts (cells x genes); in Python also an AnnData object. |
| `--output` | text | required | The annotated .h5ad to write (the input with the calls in obs). In Python, None returns the calls without writing them. |
| `--predictions` | text | none | Also write the calls per cell as a tab-separated table to this file. |
| `--gene-col` | text | index | The var column with the gene names; index uses the var names. |
| `--alpha` | number | 0.01 | Advanced. The threshold of a unique call: a cell's best call is unique when its adjusted p-value (scarlet_fdr) is below alpha. See Advanced options. |
| `--flip` | refuse, flag | refuse | The flip test, for calls that depend on genes the cell type never showed: refuse (the default) leaves such a cell unassigned, flag keeps the call and marks it in scarlet_flip. |
| `--save-probs` | yes, no | yes | Store each cell's probability for every cell type in the annotated file (obsm['scarlet_probs']). |
| `--chunk-size` | integer | 10000 | Cells scored at a time; lower it to use less memory. It does not change the result. |
| `--ncpu` | integer | 1 | Number of CPU cores to use; it does not change the result. |
| `--verbosity` | quiet, normal, full | normal | What the command prints: normal, its pre-flight once, one progress bar and one line per fold or draw; quiet, nothing; full, everything each step prints. It does not change the result. |

## scarlet transfer

Transfer the labels of a reference to a query.

| option | values | default | description |
| :--- | :--- | :--- | :--- |
| `--reference` | text | required | The labelled reference: an .h5ad file of raw counts (cells x genes); in Python also an AnnData object. |
| `--labels` | text | required | The column of the reference's obs that holds the cell-type labels. |
| `--query` | text | required | The cells to annotate: an .h5ad file of raw counts; in Python also an AnnData object. |
| `--output` | text | required | The output folder; it is created if it does not exist. |
| `--origin-col` | text | none | A column of the query's obs naming where each cell comes from (a batch, a technology or a dataset; never the labels). Depth matching, refinement and the transfer check are done per origin. |
| `--refine` | yes, no | yes | Refine the first-pass calls on the query's own stably anchored cells. |
| `--depth-match` | yes, no | yes | Match the sequencing depth of reference and query before the first pass. |
| `--qc` | yes, no | yes | Run the transfer check: per cell type, the share of its reference markers that the query confirms; below 0.2 the cell type is flagged doubtful transfer. |
| `--technology` | auto, sparse, dense | auto | The technology preset: auto chooses dense when the shallower of reference and query (after depth matching) detects at least 5,000 genes per cell on average, else sparse. |
| `--blocklist` | text | auto | Genes kept out of the markers: auto picks the bundled human or mouse list from the gene names; human, mouse, none, or the path of your own tab-separated file. |
| `--blocklist-modules` | text | auto | Advanced. Which gene groups of the bundled lists to apply: auto, all, or a comma-separated list of groups. See Advanced options. |
| `--imbalanced` | auto, yes, no | auto | Unequal class sizes: auto corrects when the largest cell type holds more than a third of the cells, yes always, no never. |
| `--gene-col` | text | index | The var column with the gene names (in both datasets); index uses the var names. |
| `--n-markers` | integer | none | Advanced. Markers per cell type; by default the technology preset's (50). |
| `--marker-shortlist` | integer | none | Advanced. The shortlist the markers are ranked on; by default the technology preset's (150). |
| `--n-iterations` | integer | 100 | Advanced. Rounds of the marker search, in each learn of the transfer. |
| `--n-cells` | integer | 100 | Advanced. Cells drawn per cell type in each round of the marker search. |
| `--draws` | integer | none | How many times the transfer is repeated, each time with its own seed: by default 10 (accurate mode), or the number of seeds; 1 is quick mode (one draw, no consensus). A final call needs its outcome in at least 80 % of the draws, rounded up (8 of 10). See Advanced options. |
| `--seeds` | text | none | Advanced. The seeds of the draws, one per draw (on the command line separated by commas); by default 1, 2, ... up to draws. The same inputs, options and seeds give the same result. |
| `--alpha` | number | 0.01 | Advanced. The threshold of a unique call; in a transfer it also selects the stable anchors the refinement learns from. See Advanced options. |
| `--flip` | flag, refuse | flag | The flip test on the first pass: flag records it per cell and keeps the calls, refuse refuses a flipped cell. |
| `--chunk-size` | integer | 10000 | Cells scored at a time; lower it to use less memory. It does not change the result. |
| `--ncpu` | integer | 1 | Number of CPU cores to use: the draws run side by side, up to ncpu at a time (faster, with more memory: each draw holds its own copy of the data). It does not change the result. |
| `--report` | yes, no | no | Also write the run's HTML report when the transfer ends. |
| `--verbosity` | quiet, normal, full | normal | What the command prints: normal, its pre-flight once, one progress bar and one line per fold or draw; quiet, nothing; full, everything each step prints. It does not change the result. |

## scarlet report

Write the HTML report of a finished run.

| option | values | default | description |
| :--- | :--- | :--- | :--- |
| `RUN` | text | required | What to report on: a model file (learn), a transfer output folder, or an annotated .h5ad (assign). |
| `--reference` | text | none | The reference a learn or transfer run used, when its recorded path is no longer readable. |
| `--query` | text | none | The query a transfer run used, when its recorded path is no longer readable. |
| `--model` | text | none | The model an assign run used (needed for an assign report). |
| `--query-labels` | text | none | Transfer and assign reports: an obs column of the query (transfer) or of the annotated data (assign) holding labels you already know for its cells, e.g. in a tutorial. The report then adds a section comparing SCARLET's calls with them (share refused, F1 on the called cells, per label), marked as labels supplied by you; SCARLET never reads them. |
| `--output` | text | none | The report folder; by default next to the run (<stem>_report/, or report/ in a transfer folder). |
| `--self-assessment` | yes, no | yes | Learn and transfer reports: estimate how well the model calls its own reference, by cross-validation (it learns once per fold, so it takes longer). |
| `--n-folds` | integer | 5 | Folds of the self-assessment. |
| `--seed` | integer | 0 | Seed of the self-assessment's folds and of the UMAP layout. |
| `--umap-cells` | integer | 50000 | At most this many cells in the report's UMAP plots. |
| `--umap-all-genes` | yes, no | yes | Also compute each UMAP on all genes, a second space the reader can switch to (the place to look for populations the model lacks); the default space is the model's own genes. No skips it (faster on very large data); the report then says so. |
| `--alpha` | number | 0.01 | Advanced. The threshold of a unique call in the learn report's self-assessment. An assign report uses the alpha the assign run recorded (a different value given here is an error). |
| `--ncpu` | integer | 1 | Number of CPU cores to use; it does not change the result. |
| `--verbosity` | quiet, normal, full | normal | What the command prints: normal, its pre-flight once, one progress bar and one line per fold or draw; quiet, nothing; full, everything each step prints. It does not change the result. |
