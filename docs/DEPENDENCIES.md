# Dependencies and their licences

SCARLET itself is distributed under the BSD 3-Clause licence ([LICENSE](https://github.com/matteo14c/scarlet-sc/blob/main/LICENSE)). It does not
include the code of the packages below; pip or conda installs them. They are the packages `pyproject.toml` requires
(and the `report` extra), with their licences as declared by each package and the versions SCARLET 0.6.3's tests ran
with. For the packages SCARLET requires directly (anndata, scanpy, numpy, pandas, scipy, scikit-learn) these are also
the lowest versions it accepts: nothing older is tested.

| package | tested version | licence | used for |
| :--- | :--- | :--- | :--- |
| anndata | 0.12.14 | BSD-3-Clause | reading and writing .h5ad |
| scanpy | 1.11.5 | BSD-3-Clause | the report's UMAPs (normalisation, PCA, neighbours, UMAP) |
| numpy | 2.4.4 | BSD-3-Clause (with bundled 0BSD, MIT, Zlib, CC0-1.0 parts) | arrays |
| pandas | 2.3.3 | BSD-3-Clause | tables |
| scipy | 1.17.1 | BSD-3-Clause | sparse matrices, statistics |
| scikit-learn | 1.8.0 | BSD-3-Clause | the self-assessment's folds |
| statsmodels | 0.14.6 | BSD-3-Clause | multiple-testing correction (Benjamini-Hochberg) |
| networkx | 3.6.1 | BSD-3-Clause | imported by the marker search |
| h5py | 3.16.0 | BSD-3-Clause | HDF5 access |
| numba | 0.65.1 | BSD-2-Clause | compiled loops |
| tqdm | 4.67.3 | MPL-2.0 and MIT | progress bars |
| matplotlib | 3.10.9 | Matplotlib licence (PSF-based) | the report's static figures (PNG, SVG); installed with scanpy |
| umap-learn | 0.5.12 | BSD-3-Clause | the report's UMAPs, through scanpy (below 0.6 with Python 3.11) |
| pynndescent | 0.5.13 | BSD-2-Clause | the neighbours of the report's UMAPs, through scanpy (below 0.7 with Python 3.11) |
| plotly (extra `report`) | 7.1.0 | MIT | the interactive HTML report |

scanpy computes the report's UMAPs with umap-learn and pynndescent, which it installs. scanpy up to 1.12.4 does not
work with umap-learn 0.6 and pynndescent 0.7 (the report's UMAP fails on some data); scanpy 1.12.6 does, and needs
Python 3.12 or later. So SCARLET requires, with Python 3.12 or later, scanpy 1.12.6 or later, and with Python 3.11
umap-learn below 0.6 and pynndescent below 0.7. SCARLET runs on Python 3.11 to 3.13 (tested on each). SCARLET does no clustering and needs no GPL-licensed package: none of the packages above is under the GPL.

The `tutorials` extra adds JupyterLab, ipykernel and ipywidgets (to run the notebooks); the `docs` extra the pinned mkdocs packages (to build the documentation locally).

The example data are not part of the repository; their sources and licences: [Tutorials](tutorials.md#data-sources-and-licences)
(and `ATTRIBUTION.md` in the release assets). The logo's lettering uses Cormorant Garamond and Inter (SIL Open Font License), converted to outlines.
