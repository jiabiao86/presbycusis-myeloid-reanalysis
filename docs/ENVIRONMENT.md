# Analysis environment

## R

| Component | Version | Notes |
| --- | --- | --- |
| R | 4.6.1 | aarch64-apple-darwin (Apple silicon) |
| limma | 3.68.5 | Ordinal and pairwise severity models |
| CellChat | 1.6.1 | Installed into a project-local library (`Rlibs/`), not the system library |
| Seurat | 5.5.1 | Myeloid subclustering |
| Matrix | 1.7.5 | Sparse matrix operations |
| WGCNA | 1.74 | Version recorded when the network analysis was run; reinstall with `BiocManager::install("WGCNA")` if the module is missing locally |

## Python

The analysis environment recorded at the time of the original run used
Python 3.12 with NumPy 2.5.3, pandas 3.0.6, SciPy 1.18.1,
scikit-learn 1.9.1, statsmodels 0.15.0, and Scanpy 1.12.4.

Document assembly, supplementary workbook generation, and archive packaging use
python-docx, openpyxl, and pandas 2.2.3. PDF rendering checks used the bundled
LibreOffice runtime and pypdf.

## CellChat runtime notes

CellChat 1.6.1 calls `future.apply::future_sapply` for every coreceptor
evaluation inside `computeCommunProb`. With a multisession `future` plan this
dominates the runtime, so the exhaustive database analysis runs each chunk with
`future::plan("sequential")` and parallelises across chunks instead. Twelve
chunks (four per age, 505 interactions each) were executed with a concurrency of
four; each chunk reported 220-310 s of compute for the 3- and 12-month groups.

## Reference builds

Gene symbols follow the Mouse Genome Informatics nomenclature used by the
CellChat mouse database. Two database genes (H2-BI and H2-Ea-ps) are not present
in the official symbol table used by `extractGene()`; they were restored after
gene filtering because several database complexes reference them.
