# Tutorials

Three notebooks in `notebooks/` show SCARLET on the data of its benchmark: what each step prints, what it writes,
and how to read the result. Every code cell is explained, and every setting is at the package default.

| notebook | what it teaches | data (benchmark task) | download | time, cluster node |
| :--- | :--- | :--- | :--- | :--- |
| [01_learn_assign_sparse](https://github.com/matteo14c/scarlet-sc/blob/main/notebooks/01_learn_assign_sparse.ipynb) | `learn` and `assign` on droplet data; the pre-flight; refusals and their reasons; the assign and learn reports; scoring the calls against held-out labels | Baron et al. 2016, human pancreas (inDrop), 8,562 cells (`A/pancreas_baron_human`) | 37 MB | about 1 min |
| [02_learn_assign_dense](https://github.com/matteo14c/scarlet-sc/blob/main/notebooks/02_learn_assign_dense.ipynb) | the same on plate data, where SCARLET chooses the dense technology preset | Segerstolpe et al. 2016, human pancreas (Smart-seq2), 2,111 cells (`A/pancreas_segerstolpe`) | 38 MB | under 1 min |
| [03_transfer](https://github.com/matteo14c/scarlet-sc/blob/main/notebooks/03_transfer.ipynb) | `transfer` across samples (mouse brain, VISp ↔ ALM, 34 cell types) and across technologies (PBMC, Smart-seq2 → 10x v3); the transfer report; agreement across draws; what happens to populations the reference lacks | Tasic et al. 2018, mouse cortex (`B/Brain_MouseV1_MouseALM_Labels34`); Ding et al. 2020, PBMC (`B/PbmcBench_Smart-Seq2`, `B/PbmcBench_10Xv3`) | 701 MB | about 8 min (2 CPUs) |

Times were measured on a cluster node with each notebook's own setting (4 CPUs for 01 and 02, 2 for 03). On a
laptop (14 GB of memory) notebook 03 took about 33 minutes with its 2 CPUs. Each notebook
downloads its data once, checking each file's SHA-256 (`scarlet.datasets.fetch`). The files are the inputs of
SCARLET's benchmark with the benchmark's own selection of cells, compressed; their content is proven identical to
the benchmark's ([Data sources and licences](#data-sources-and-licences)). The notebooks run SCARLET at its package
defaults (the benchmark ran without the blocklist), on one split or pair, so their numbers can differ from the
[Benchmarks](https://github.com/matteo14c/scarlet-sc/blob/v0.6.3/benchmark/index.md).

## Data sources and licences

The ten files of the [v0.6.3 release](https://github.com/matteo14c/scarlet-sc/releases/tag/v0.6.3). The seven data
files are benchmark inputs of Abdelaal T, et al. *A comparison of automatic cell identification methods for
single-cell RNA sequencing data.* Genome Biology 20, 194 (2019), doi:10.1186/s13059-019-1795-z (data: Zenodo record
3357167, doi:10.5281/zenodo.3357167, CC BY 4.0), each with its benchmark task's own selection of cells and nothing
else changed.

| file | study | accession | licence |
| :--- | :--- | :--- | :--- |
| `baron_human.h5ad` (notebook 01) | Baron M, et al. Cell Systems 3, 346–360 (2016) | GEO GSE84133 | CC BY 4.0, as distributed by the benchmark (Zenodo record 3357167); the original data remain under their repository's terms |
| `segerstolpe.h5ad` (notebook 02) | Segerstolpe A, et al. Cell Metabolism 24, 593–607 (2016) | ArrayExpress E-MTAB-5061 | CC BY 4.0, as distributed by the benchmark (Zenodo record 3357167); the original data remain under their repository's terms |
| `pbmc1_smartseq2.h5ad` (getting started, notebook 03) | Ding J, et al. Nature Biotechnology 38, 737–746 (2020) | GEO GSE132044; Single Cell Portal SCP424 | CC BY 4.0, as distributed by the benchmark (Zenodo record 3357167); the original data remain under their repository's terms |
| `pbmc1_10xv3_query.h5ad` (getting started, notebook 03) | Ding J, et al. (2020), as above | GEO GSE132044; Single Cell Portal SCP424 | CC BY 4.0, as distributed by the benchmark (Zenodo record 3357167); the original data remain under their repository's terms |
| `pbmc1_10xv3.h5ad` (notebook 03) | Ding J, et al. (2020), as above | GEO GSE132044; Single Cell Portal SCP424 | CC BY 4.0, as distributed by the benchmark (Zenodo record 3357167); the original data remain under their repository's terms |
| `amb34_visp.h5ad` (notebook 03) | Tasic B, et al. Nature 563, 72–78 (2018) | GEO GSE115746 | CC BY 4.0, as distributed by the benchmark (Zenodo record 3357167); the original data remain under their repository's terms |
| `amb34_alm.h5ad` (notebook 03) | Tasic B, et al. (2018), as above | GEO GSE115746 | CC BY 4.0, as distributed by the benchmark (Zenodo record 3357167); the original data remain under their repository's terms |
| `MANIFEST.tsv` | written by the SCARLET authors: per data file its size, SHA-256, the SHA-256 of its content (counts, cell and gene names, labels), the benchmark file and task it comes from with that file's SHA-256, the cell selection applied, and its cells, genes and labels | none | BSD 3-Clause, as SCARLET |
| `SHA256SUMS` | written by the SCARLET authors: the SHA-256 of each data file, as `sha256sum -c SHA256SUMS` reads it | none | BSD 3-Clause, as SCARLET |
| `ATTRIBUTION.md` | written by the SCARLET authors: the sources above, to cite with the data | none | BSD 3-Clause, as SCARLET |

If you use these data, cite the study and the benchmark.

## Run them on your computer

1. **Get the repository at the release** (the notebooks are in it; `--branch` takes the release's own version of
   them):

    ```bash
    git clone --branch v0.6.3 https://github.com/matteo14c/scarlet-sc
    cd scarlet-sc
    ```

2. **Make a fresh environment** with the Python you choose (3.11, 3.12 or 3.13), and activate it. Name the version
   in the command, so the environment gets that one:

    ```bash
    python3.12 -m venv .venv
    source .venv/bin/activate        # Windows: py -3.12 -m venv .venv, then .venv\Scripts\activate
    python --version                 # the version you chose
    ```

    With conda instead: `conda create -n scarlet python=3.12` and `conda activate scarlet`.

    Each new terminal starts without the environment: activate it again there first (`source .venv/bin/activate`
    from the `scarlet-sc` folder, or `conda activate scarlet`).

3. **Install SCARLET with the `tutorials` extra** (the reports, Jupyter, and the widgets of the progress bars):

    ```bash
    pip install -e ".[tutorials]"
    ```

    To install everything in one go (the reports, the tutorials and the tools to build this documentation):
    `pip install -e ".[tutorials,docs]"`.

4. **Start Jupyter**:

    ```bash
    jupyter lab notebooks/
    ```

    Jupyter keeps this terminal busy while it runs, and opens a browser tab (if it does not, open the address it
    prints). Open a notebook and run its cells in order. To stop Jupyter, press Ctrl+C in the terminal.

    Run one notebook at a time: a notebook keeps its memory until it is shut down. Stopping its cells or closing
    its tab does not free it: close the others first with File → Close and Shut Down Notebook.

    Notebook 03 uses 2 CPUs, which fits a laptop. Its memory, measured with 2 CPUs, SCARLET alone (Jupyter and
    the browser come on top): in a system monitor's process list, the notebook's own process (Python, or
    ipykernel) peaks at about 5 GB of resident memory (5.0 GB on a cluster node, 5.3 GB on a laptop); during each
    transfer two more Python processes run the draws beside it. All of SCARLET's processes together, with the
    memory they share counted once, peak at about 7.5 GB: about what the system monitor's total memory in use
    rises by.

The data download themselves on the first run. To download them beforehand, for example on a machine with a better
connection: `python -m scarlet.datasets --dest notebooks/data` (each notebook looks for its files in its own
`data/` folder). While the repository is private, its release assets need a login: the downloader then uses the
GitHub CLI (`gh auth login` first).

## Run them in Google Colab

Planned for the public release: each notebook will open in Colab from a badge on this page, with nothing installed
on your computer (a setup cell installs SCARLET there).
