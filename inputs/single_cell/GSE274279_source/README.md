# S11 Four-Library GSE274279 Descriptive Single-Nucleus Audit

## Source and scope

GSE274279 is the cochlear single-nucleus dataset reported by Xia et al. (2025,
Cell Reports 44:115781, doi:10.1016/j.celrep.2025.115781). The source mice were
C57BL/6J. This dataset supplies cross-strain descriptive context for the
CBA/CaJ bulk reanalysis and does not constitute validation in the bulk strain.

All four deposited 10x libraries are represented:

| GSM | Library | Age | Raw barcodes | Retained nuclei |
|---|---|---:|---:|---:|
| GSM8446833 | 3M | 3 months | 2,936 | 2,916 |
| GSM8446834 | 12M | 12 months | 5,581 | 5,532 |
| GSM8446835 | 24M | 24 months | 2,968 | 2,954 |
| GSM8446836 | 24M_Cochlea_BU | 24 months | 2,757 | 2,750 |

The four libraries retained 14,152 nuclei and 23,540 genes. The atlas contained
245 marker-defined macrophage/microglial nuclei: 68 at 3 months, 81 at 12
months, and 96 at 24 months. Library totals were 68, 81, 50, and 46.
GSM8446836 is a second public 24-month cochlear library whose animal-level
relationship to GSM8446835 is not documented in GEO.

## Interpretation boundary

Age, library, and 10x batch are not separable. All age-related UMAP,
composition, detection, module-score, bootstrap, and subtype summaries describe
the sampled nuclei only. They do not estimate animal-level age effects.

Cd68, Fcgr3, Mpeg1, and Ptprc overlap the marker panel used to define the
myeloid population. Detection within that population therefore does not
provide independent localization or enrichment relative to other cell types.
Clec7a was not detected in the prepared object.

No cell-level Mann-Whitney P values or BH-FDR estimates are included in the
public-facing S11 summaries. Nuclei are not independent animal-level
observations, and each age is confounded with a 10x library. The legacy
prepared-object file with those columns is retained only in the local audit
archive outside the uploadable package.

The 500-iteration bootstrap samples at most 25 nuclei per age per iteration.
Its intervals describe the sampled nuclei and are not animal-level
biological-effect confidence intervals.

## Main files

- `GSE274279_source_manifest.csv`: accession, tissue, strain, library, and
  public source-file provenance.
- `GSE274279_source_files_sha256.csv`: SHA-256 digest, compressed size, and
  uncompressed size for each source file.
- `gse274279_four_matrix/GSE274279_four_matrix_QC_by_library.csv`: per-library
  barcode and nucleus counts.
- `gse274279_four_matrix/GSE274279_four_matrix_cell_type_counts.csv`: atlas cell
  types by library.
- `gse274279_four_matrix/GSE274279_myeloid_age_descriptive_summary.csv`: gene
  detection by age.
- `gse274279_four_matrix/GSE274279_myeloid_library_descriptive_summary.csv`:
  gene detection by the four libraries.
- `gse274279_four_matrix/myeloid_four_matrix_*`: myeloid cell metadata, UMAP,
  module scores, proportions, and nucleus-level bootstrap samples.
- `gse274279_four_matrix/sessionInfo.txt`: R and package environment.

CellChat is documented separately under S12 and uses only GSM8446833,
GSM8446834, and GSM8446835.
