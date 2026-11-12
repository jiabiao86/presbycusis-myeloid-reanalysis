# Supplementary Table S13 Code Index

The minimal code archive is supplied as
`Analysis_Code.zip` relative to the submission package
root. The archive `README.md` describes dependencies, environment variables,
and the recommended execution order.

The analysis code archive is available from Zenodo at
https://doi.org/10.5281/zenodo.24169355 and from GitHub at
https://github.com/jiabiao86/presbycusis-myeloid-reanalysis.

## Code

| Stage | Files in the archive |
|---|---|
| GEO input preparation | `code/python/prepare_geo_inputs.py` |
| Sample audit | `code/python/build_sample_master.py` |
| Primary, sample-scenario, and no-sex limma | `code/R/run_limma_v14.R`, `code/R/run_no_sex_sensitivity_v14.R`, `code/R/common_v14.R` |
| GSE49543 probe sensitivity | `code/python/prepare_gse49543_probe_sensitivity.py`, `code/R/run_gse49543_probe_sensitivity_v14.R` |
| Leave-one-animal-out and scenario summary | `code/R/run_leave_one_animal_out_v14.R`, `code/R/summarize_sensitivity_v14.R` |
| WGCNA | `code/R/run_wgcna_v14.R` |
| Four-library Seurat processing and marker-defined detection | `code/R/single_cell_context_four_matrices_v14.R` |
| Four-library Seurat myeloid subclustering | `code/R/single_cell_context_four_matrices_v14.R` |
| Candidate-expression extraction | `code/R/extract_myeloid_candidate_expression_v14.R` |
| Cell-ID and count audit | `code/R/create_cell_audit_v14.R` |
| Myeloid subclustering figure | `code/R/generate_myeloid_figure.R` |
| CellChat | `code/R/cellchat_full_analysis.R`, `code/R/cellchat_full_by_age.R`, `code/R/cellchat_full_merge.R`, `code/R/cellchat_full_compare.R`, `code/R/generate_cellchat_full_figure.R` |
| Main figures and graphical abstract | `code/R/generate_figures_v14.R`, `code/R/generate_graphical_abstract_v14.R` |
| Document construction | `code/python/build_submission_docs.py` |
| External matrix preparation | `code/python/prepare_external_inputs.py` |
| External Welch and BH analysis | `code/R/run_external_welch_bh.R`, `code/R/validate_external_welch_bh.R` |
| P-value display | `code/python/add_pvalue_display_columns.py` |
| Integrity manifest | `code/python/build_manifest.py` |

## Environment and validation

- R package versions: `environment/R_packages.csv` and
  `environment/R_sessionInfo_v14.txt`.
- Python dependencies: `environment/requirements.txt` and
  `environment/Python_environment_record.json`.
- Rebuild validation: `validation/VALIDATION_REPORT.md`.
- GEO rebuild outputs: `validation/geo_rebuild/`.
- Primary limma rebuild: `validation/primary_rebuild/`.
- Sensitivity limma rebuild: `validation/sensitivity_rebuild/`.

## Result layers

- Primary analysis: `results/primary/`.
- Sample-scenario analyses:
  `results/sensitivity_exclude_GSM1201568/`,
  `results/sensitivity_exclude_GSM1201602/`, and
  `results/sensitivity_exclude_both_unmatched/`.
- Scenario and leave-one-animal-out summaries:
  `results/sensitivity_summary/` and `results/leave_one_animal_out/`.
- No-sex sensitivity: `results/sensitivity_no_sex/`.
- Probe-aggregation sensitivity: `results/sensitivity_probe_aggregation/`.
- WGCNA: `results/primary/wgcna/`.
- Four-library GSE274279 source records and descriptive outputs:
  `results/single_cell/GSE274279_four_matrix/`.
- Prepared CellChat tables: `inputs/prepared/cellchat_full/`.
- External sample and run mappings, exclusions, normalization metadata, and
  input checksums: selected files under `results/external/`.
- External candidate-level Welch/BH results and validation summaries:
  `results/external/`.
- Numerical-underflow provenance:
  `results/external/external_p_value_underflow.csv`.
- Independent Welch/BH validation:
  `results/external/external_welch_validation*.csv` and
  `results/external/external_welch_validation.json`.

## Integrity

`ANALYSIS_MANIFEST.csv` lists every packaged file with its relative path,
size, and SHA-256 digest. `SHA256SUMS.csv` contains the corresponding
checksum list. The manifest and checksum files do not checksum themselves.

Public GEO workbooks, complete prepared matrices, and complete
transcriptome-wide Welch/BH tables are not duplicated in the minimal archive.
The included preparation and analysis scripts rebuild those files from the
cited public GEO records.
