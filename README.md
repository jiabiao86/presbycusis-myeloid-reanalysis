# Peripheral-central convergence of macrophage and microglial programs in presbycusis severity

This public repository hosts the reproducibility materials for Revision 8 of the
audited cochlea/inferior-colliculus secondary reanalysis.

## R8 reproducibility archive

Download:

- `V14_Code_Archive_R8_public.zip`
- SHA-256:
  `57f6d5d27843b891c0dedd80e7f4f7dc74ad12498eb4f39012bd095f89d4937f`

Extract the archive before running any script. Its `README.md` gives the
software environment, input requirements, execution order, validation scope,
and interpretation limits.

The archive contains:

- R and Python source code for the bulk microarray, sensitivity, WGCNA,
  single-cell, CellChat, external-validation, figure, and submission-QA steps.
- Pinned R and Python environment records.
- Selected result tables, including primary analyses, sample scenarios,
  leave-one-animal-out models, no-sex sensitivity, probe aggregation, WGCNA,
  single-cell descriptive outputs, and external Welch/BH validation.
- GEO source manifests and checksums without redistributing the public source
  matrices.
- Rebuild validation records and complete file-level SHA-256 manifests.
- MIT/CC BY licensing, citation metadata, and a code-to-analysis index.

The included code and result layers support audit and reproduction of the
reported analyses. They do not include cover letters, reviewer responses,
editorial correspondence, internal audit history, or Word/PDF submission
copies.

## Data sources

All source data are public mouse transcriptomic records from GEO:

- GSE49543: mouse cochlea microarray
- GSE49522: mouse inferior colliculus microarray
- GSE233798: mouse cochlea validation data
- GSE153882: mouse cochlear hair-cell data
- GSE154833: mouse cochlea validation data
- GSE274279: mouse cochlea single-nucleus RNA-seq

The original GEO datasets remain subject to the terms of their depositors.

## DOI and repository status

This GitHub repository is the R8 code and result repository. No new R8 Zenodo
DOI is claimed here. The earlier Zenodo DOI
`10.5281/zenodo.23737275` applies only to Revision 3; a separate R8 Zenodo
deposit is still required if an archival DOI is requested.

Files from earlier revisions remain in the Git history and as legacy paths in
the repository. The R8 archive linked above is the authoritative reproducibility
package for the current manuscript.

## License and citation

Code is released under the MIT License. Derived result tables and figure sources
are additionally released under CC BY 4.0. Citation metadata are provided in
`CITATION.cff`.

Corresponding author: Jianming Yang, The Second Affiliated Hospital of Anhui
Medical University, Hefei, Anhui, China.
