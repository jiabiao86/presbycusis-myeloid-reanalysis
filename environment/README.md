# Environment record

## R

The primary and sensitivity limma analyses were run with R 4.6.1. Direct
package versions are listed in `R_packages.csv`. The complete session records
for the primary and sensitivity runs are under:

```text
results/primary/sessionInfo.txt
results/sensitivity_exclude_GSM1201568/sessionInfo.txt
```

The archive-wide loaded-package record is `R_sessionInfo_v14.txt`.

## Python

`requirements.txt` records the Python 3.12 dependencies used by the GEO and
supporting preparation scripts. The single-nucleus processing and myeloid
subclustering were implemented in R with Seurat; the original single-cell
matrices are not redistributed in this archive.

The GEO preparation script requires LibreOffice only for the deposited
GSE49543 `.xls` conversion. Set `V14_SOFFICE` when `soffice` is not on `PATH`.
