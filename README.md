# Analysis Code Archive

This archive accompanies the manuscript "Intersection of Myeloid-Associated
Genes in Aging CBA/CaJ Mouse Auditory Tissues: A Secondary Reanalysis with
Cross-Strain Single-Nucleus Context". It preserves the unordered four-cohort
model, the 41-sample primary analysis, mirrored exclusion of the two unmatched
cochlear records, their joint exclusion, and 41 leave-one-animal-out models.
It also records all-four-library processing of the C57BL/6J GSE274279
cochlear single-nucleus context and the GSE49543 probe-aggregation sensitivity.

This analysis code archive is available from Zenodo at
https://doi.org/10.5281/zenodo.24169355. The GitHub repository is
https://github.com/jiabiao86/presbycusis-myeloid-reanalysis.

The author order is Jiabiao Ji, Xiaoqing Yu, Chao Zhang, Long Liu, Jun Song,
Lin Yan, Shan Ren, and Jianming Yang. Jiabiao Ji is the first author, Xiaoqing
Yu is the co-first author, and Jianming Yang is the corresponding author.

## Repository layout

- The repository root contains the analysis code, selected frozen result
  tables, environment records, validation outputs, and integrity manifests.
- `submission/` contains the manuscript, cover letter, highlights, figure
  legends, main and supplementary figures, supplementary materials, analysis
  code archive, and checksum manifests prepared for initial journal submission.
- `submission/Analysis_Code.zip` is the single analysis executable archive
  deposited in Zenodo.

## Reproducible scope

- Preparation of the two bulk microarray matrices from the supplied GEO
  series-family files.
- Animal-level sample audit and construction of the 80-record metadata table.
- Primary and sample-scenario limma/eBayes analyses, including exclusion of
  GSM1201568, exclusion of GSM1201602, exclusion of both unmatched records,
  and leave-one-animal-out resampling.
- The no-sex sensitivity analysis under `~ 0 + cohort_group`, with the
  universe, contrasts, limma/eBayes settings, BH correction, and thresholds
  held fixed.
- GSE49543 probe-to-symbol sensitivity using first-probe, mean-probe, and
  highest-mean-probe aggregation while holding the inferior-colliculus
  result and significance rule fixed.
- WGCNA outputs generated under the unordered cohort model.
- Figure scripts and the compact source tables required by those scripts.
- Seurat processing for all four GSE274279 libraries, plus CellChat
  scripts for the three-library communication subset.
- Candidate-expression extraction for the six genes measurable in the
  marker-defined myeloid population.
- Preparation of the GSE233798, GSE153882, and GSE154833 external expression
  matrices, including biological-sample and sequencing-run mappings,
  missing-value exclusions, duplicate-symbol handling, and normalization.
- Transcriptome-wide two-sided Welch tests with Welch-Satterthwaite degrees of
  freedom, Benjamini-Hochberg correction, numerical-underflow traceability, and
  an independent R `t.test`/`p.adjust` validation.

The original GSE274279 10x matrices are not redistributed. Source file names,
sizes, checksums, library metadata, and the complete four-library descriptive
summaries are supplied so the reported figures can be audited. A complete rerun
of the Seurat layers requires the public GSE274279 matrices and the R
environment described under `environment/`.

Each feature annotation is filtered before concatenation. For duplicated gene
symbols, `collapse_duplicate_symbols()` retains the first occurrence with
`!duplicated(symbols)`, matching the rule reported in Methods 2.6.

GSE274279 contains four libraries: GSM8446833 (3M), GSM8446834 (12M),
GSM8446835 (24M), and GSM8446836 (24M_Cochlea_BU). GSM8446836 is a second
public 24-month cochlear library whose animal-level relationship to GSM8446835
is not documented in GEO. Age, library, and 10x batch are inseparable. UMAP,
composition, marker detection, and subclustering are descriptive rather than
biological age-effect estimates. CellChat uses only GSM8446833, GSM8446834,
and GSM8446835.

This minimal archive includes source-to-run mappings, exclusion and
duplicate-symbol logs, candidate-level Welch/BH results, numerical-underflow
records, and independent validation summaries for GSE233798, GSE153882, and
GSE154833. It does not redistribute the public GEO workbooks, complete
prepared matrices, or complete transcriptome-wide Welch/BH result tables; the
included scripts rebuild them from the cited GEO records. The external
analyses are descriptive and are not independent replications of the primary
contrast.

## Software

- R 4.6.1 and the direct package versions listed in
  `environment/R_packages.csv`.
- Python 3.12 with the packages listed in `environment/requirements.txt`.
- LibreOffice for conversion of the deposited GSE49543 `.xls` workbook.

No machine-specific library or input path is used as a code default. Set the
following optional variables when a local library or alternate input directory
is required:

```text
V14_R_LIB
V14_SOFFICE
V14_GEO_DIR
V14_PREPARED_DIR
V14_PRIMARY_RESULTS_DIR
V14_SENSITIVITY_RESULTS_DIR
V14_PROBE_SENSITIVITY_DIR
V14_PROBE_SENSITIVITY_OUTPUT_DIR
V14_DERIVED_DIR
V14_FIGURE_DIR
V14_SINGLE_CELL_DIR
V14_SINGLE_CELL_OUTPUT_DIR
```

## Reproduction sequence

Run all commands from the archive root.

### 1. Prepare the bulk expression inputs

Download the public GSE49543 and GSE49522 source files described by
`code/python/prepare_geo_inputs.py` and place them under `inputs/geo/`, unless
`V14_GEO_DIR` points to an alternate download directory.

```bash
python3 code/python/prepare_geo_inputs.py
```

The script uses the GSE49543 further-processed workbook, GSE49522 series
matrix, and embedded GPL339 annotation. It writes:

```text
inputs/prepared/expression_GSE49543_log2.csv
inputs/prepared/expression_GSE49522_log2.csv
inputs/prepared/geo_preparation_manifest.json
```

If `soffice` is not available on `PATH`, set `V14_SOFFICE` to its executable.

### 2. Build the animal-level sample audit

```bash
python3 code/python/build_sample_master.py
```

This writes `inputs/prepared/sample_master.csv`, the audit JSON, and the build
log. Animal 13 is recorded as the unmatched SP record (GSM1201602), and animal
16 as the unmatched YC record (GSM1201568).

### 3. Run the primary limma analysis

```bash
Rscript code/R/run_limma_v14.R
```

The script writes complete contrast tables, design matrices and ranks,
significant-gene and intersection summaries, input checksums, and session
information under `results/primary/`.

Run the no-sex sensitivity while holding the 12,437-gene universe, contrasts,
limma/eBayes settings, BH procedure, and thresholds fixed:

```bash
Rscript code/R/run_no_sex_sensitivity_v14.R
```

This writes design matrices and complete result tables under
`results/sensitivity_no_sex/`. All seven primary SP-versus-YC intersection
genes were retained without sex adjustment with concordant directions.

Run the GSE49543 probe-aggregation sensitivity:

```bash
python3 code/python/prepare_gse49543_probe_sensitivity.py
Rscript code/R/run_gse49543_probe_sensitivity_v14.R
```

The three rules select the first valid probe, average all probes for a symbol
in linear signal space, or select the probe with the highest across-sample mean
linear signal. The fixed inferior-colliculus result and threshold rule are used
for each intersection. First-probe, mean-probe, and highest-mean-probe
intersections contain 7, 8, and 9 genes; the seven primary genes are retained
under every rule.

### 4. Run the sample-scenario and leave-one-animal-out analyses

```bash
V14_PRIMARY_RESULTS_DIR=results/primary \
V14_RUN_LABEL=sensitivity_exclude_GSM1201568 \
V14_EXCLUDE_COCHLEA_GSM=GSM1201568 \
V14_FREEZE_GENE_UNIVERSE=true \
V14_RESULTS_DIR=results/sensitivity_exclude_GSM1201568 \
Rscript code/R/run_limma_v14.R
```

This excludes only GSM1201568, retains GSM1201602, and holds the 12,437-gene
primary universe fixed. Repeat the same command with the following changes:

```text
V14_RUN_LABEL=sensitivity_exclude_GSM1201602
V14_EXCLUDE_COCHLEA_GSM=GSM1201602
V14_RESULTS_DIR=results/sensitivity_exclude_GSM1201602
```

```text
V14_RUN_LABEL=sensitivity_exclude_both_unmatched
V14_EXCLUDE_COCHLEA_GSM=GSM1201568,GSM1201602
V14_RESULTS_DIR=results/sensitivity_exclude_both_unmatched
```

Rebuild the 41 leave-one-animal-out models and the scenario summary:

```bash
Rscript code/R/run_leave_one_animal_out_v14.R
Rscript code/R/summarize_sensitivity_v14.R
```

For each of the 41 cochlear animal IDs, the leave-one-animal-out analysis
removes that animal from the cochlear model and, when a paired
inferior-colliculus sample exists, from the inferior-colliculus model. Both
tissue models are then refitted. The frozen universe and significance rule are
held fixed, and the intersection size and candidate-gene membership are
recorded.

### 5. Rebuild WGCNA from the frozen analysis universe

```bash
Rscript code/R/run_wgcna_v14.R
```

WGCNA writes to `results/primary/wgcna/`. Archived validation rebuilt both
networks from the frozen 12,437-gene universe and confirmed the selected powers,
module counts, and module-level results recorded in Table 3.

### 6. Build the four-library single-cell context and figures

```bash
V14_SINGLE_CELL_DIR=/path/to/GSE274279_10x \
Rscript code/R/single_cell_context_four_matrices_v14.R
Rscript code/R/extract_myeloid_candidate_expression_v14.R
Rscript code/R/create_cell_audit_v14.R
Rscript code/R/generate_figures_v14.R
Rscript code/R/generate_myeloid_figure.R
Rscript code/R/generate_graphical_abstract_v14.R
```

The full analysis project writes figure outputs under `figures/source/`. This
minimal code archive does not include generated figure files; use the approved
PNG files from the submission package when rebuilding documents.

To rebuild only Figure 1 after changes to the workflow layout:

```bash
V14_RENDER_FIGURE1_ONLY=true Rscript code/R/generate_figures_v14.R
```

### 7. Rebuild the external GEO analyses

Download the three public processed workbooks named in
`code/python/prepare_external_inputs.py` and place them under
`inputs/external/geo/`.

```bash
python3 code/python/prepare_external_inputs.py
Rscript code/R/run_external_welch_bh.R
Rscript code/R/validate_external_welch_bh.R
python3 code/python/prepare_external_inputs.py
```

When the preparation script is rerun, it records checksums for the completed
Welch/BH, underflow, and independent validation outputs in
`inputs/external/prepared/external_input_manifest.json`. The archival copies
of that manifest and the selected external outputs are under
`results/external/`.

### 8. Rebuild submission documents

The compact code archive intentionally does not include the manuscript Markdown
and supporting-document sources or the generated figure PNGs. To rebuild the
four editable documents, provide a local project root containing the Markdown
files and a directory
containing the approved figure PNGs:

```bash
python3 code/python/build_submission_docs.py \
  --project-root /path/to/project-root \
  --package-root /path/to/document-output \
  --figure-dir /path/to/approved-figures
```

After the package is complete, rebuild the analysis manifest and SHA-256 list:

```bash
python3 code/python/build_manifest.py \
  --package /path/to/analysis-package
```

The code archive has its own file-level integrity records under `manifests/`.
The submission package manifest records the complete submission package,
including the compressed code archive.

## Single-cell scripts

To rerun the four-library descriptive context, provide the public GSE274279 10x
files under
`inputs/single_cell/GSE274279_10x/`:

```text
GSM8446833_3M_Cochlea_matrix.mtx.gz
GSM8446833_3M_Cochlea_barcodes.tsv.gz
GSM8446833_3M_Cochlea_features.tsv.gz
GSM8446834_12M_Cochlea_matrix.mtx.gz
GSM8446834_12M_Cochlea_barcodes.tsv.gz
GSM8446834_12M_Cochlea_features.tsv.gz
GSM8446835_24M_Cochlea_matrix.mtx.gz
GSM8446835_24M_Cochlea_barcodes.tsv.gz
GSM8446835_24M_Cochlea_features.tsv.gz
GSM8446836_24M_Cochlea_BU_matrix.mtx.gz
GSM8446836_24M_Cochlea_BU_barcodes.tsv.gz
GSM8446836_24M_Cochlea_BU_features.tsv.gz
```

Then run:

```bash
Rscript code/R/single_cell_context_four_matrices_v14.R
Rscript code/R/extract_myeloid_candidate_expression_v14.R
Rscript code/R/create_cell_audit_v14.R
Rscript code/R/generate_myeloid_figure.R
```

The prepared `inputs/prepared/cellchat_full/` tables support the archived
CellChat summary and figure. Rebuilding CellChat itself requires the public
single-cell matrices and the CellChat scripts in `code/R/`.

## Validation

`validation/VALIDATION_REPORT.md` records the clean rebuild checks for GEO
preparation and limma, exact sample counts, design ranks, four sample
scenarios, the 41 leave-one-animal-out models, `Tyrobp` values, probe
aggregation, all four GSE274279 libraries, and software records. The minimal
archive retains selected summary outputs supporting those checks; complete
regenerated matrices and large R objects are not duplicated.

## Integrity

`manifests/ANALYSIS_MANIFEST.csv` lists each file in the code archive
with its relative path, size, and SHA-256 digest.
`manifests/SHA256SUMS.csv` is the corresponding checksum list. These
two records do not checksum themselves.

## P-value display

Exact zero P values in public-facing text and display columns are reported as
`<2.2e-16`. Numeric `pval`, `p_value`, and `fdr_bh` fields remain unchanged so
that supplied calculations can be checked. CellChat tables use
`pval_display`; external derived tables use `p_value_display` and
`fdr_bh_display`.

To add or refresh the display columns without changing the numeric fields:

```bash
python3 code/python/add_pvalue_display_columns.py \
  inputs/external/*.csv \
  inputs/prepared/cellchat_full/analysis/CellChat_full_edges_all_ages.csv \
  inputs/prepared/cellchat_full/analysis/CellChat_full_focus_interactions.csv
```

## Interpretation limits

The cohort contrast combines age and source-label hearing status. Individual
animal ages, ABR thresholds, and DPOAE values were not deposited. CellChat
values are model-inferred database interactions and do not measure signal flux
or demonstrate a functional pathway.
