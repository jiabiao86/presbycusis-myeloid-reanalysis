# Analysis Code Archive Validation Report

Validation date: 2026-10-30

All paths in this report are relative to the archive root. The validation used
the archived inputs and scripts, with R 4.6.1 and the package versions recorded
under `environment/`.

The validation covers the cross-strain provenance and interpretation of
GSE274279, complete representation of all four deposited 10x libraries,
public single-nucleus summaries without cell-level inferential statistics,
nucleus-level bootstrap boundaries, a three-probe aggregation sensitivity
analysis for GSE49543, and the first-occurrence duplicate-symbol rule. The four
libraries contributed 14,152 retained nuclei and 245 marker-defined myeloid
nuclei. Age, library, and 10x batch remain inseparable. The archive is
available in Zenodo at https://doi.org/10.5281/zenodo.24169355.

## Bulk expression reconstruction

The GEO preparation script was rerun from the archived GSE49543 and GSE49522
source files. The regenerated matrices matched the frozen prepared matrices
byte for byte.

| Dataset | Rows | Samples | Output SHA-256 |
|---|---:|---:|---|
| GSE49543 | 13,649 | 41 | `6ac19ab7bfddd28750a67a89f7e5db27bb837e24f7746b4cfb94493fe541c622` |
| GSE49522 | 13,049 | 39 | `dd32df2eb8c50d8492926deae662e29b4acb44130d34d27d76f59a7e3d577df5` |

The source checksums and transformation rules are recorded in
`geo_rebuild/geo_preparation_manifest.json`.

## Sample audit

- The archive contains 80 GSM-level records: 41 cochlear and 39
  inferior-colliculus records.
- Thirty-nine animals have concordant records in both tissues.
- Animal 13 (GSM1201602) is the unmatched SP cochlear record.
- Animal 16 (GSM1201568) is the unmatched YC cochlear record.
- Both unmatched records are retained in the 41-sample primary analysis.
- The prespecified sensitivity analysis excludes GSM1201568 and retains
  GSM1201602. The mirrored analysis excludes GSM1201602 and retains
  GSM1201568; the joint analysis excludes both records.
- All sample-scenario and leave-one-animal-out analyses hold the primary
  12,437-gene universe and significance rule fixed.

## Model and design checks

The prespecified model was `~ 0 + cohort_group + sex`. All fitted designs had
five columns and rank five.

| Analysis | Tissue | Samples | Design columns | Rank | Full rank |
|---|---|---:|---:|---:|---|
| Primary | Cochlea | 41 | 5 | 5 | TRUE |
| Primary | Inferior colliculus | 39 | 5 | 5 | TRUE |
| Exclude GSM1201568 | Cochlea | 40 | 5 | 5 | TRUE |
| Exclude GSM1201602 | Cochlea | 40 | 5 | 5 | TRUE |
| Exclude both unmatched | Cochlea | 39 | 5 | 5 | TRUE |
| Inferior colliculus in all scenarios | Inferior colliculus | 39 | 5 | 5 | TRUE |

The primary and sensitivity limma runs used the same thresholds:
Benjamini-Hochberg adjusted P < 0.05 and absolute log2 fold change >= 1.

## Differential-expression rebuild

The rebuilt primary and sample-scenario tables were identical to the canonical
tables under `results/primary/limma/`, the three
`results/sensitivity_exclude_*/limma/` directories, and the
leave-one-animal-out output directory.

| Analysis | Tissue | MA vs YC | SP vs MP | SP vs YC | SP vs YC direction |
|---|---|---:|---:|---:|---|
| Primary | Cochlea | 0 | 0 | 18 | 16 up, 2 down |
| Primary | Inferior colliculus | 0 | 0 | 46 | 46 up, 0 down |
| Exclude GSM1201568 | Cochlea | 0 | 0 | 18 | 16 up, 2 down |
| Exclude GSM1201602 | Cochlea | 0 | 0 | 18 | 16 up, 2 down |
| Exclude both unmatched | Cochlea | 0 | 0 | 18 | 16 up, 2 down |
| Inferior colliculus in all scenarios | Inferior colliculus | 0 | 0 | 46 | 46 up, 0 down |

No gene met both criteria in both tissues for MA versus YC or SP versus MP.

## No-sex sensitivity

The no-sex analysis used `~ 0 + cohort_group` for both tissues and retained
the 12,437-gene universe, three contrasts, limma/eBayes settings,
Benjamini-Hochberg procedure, and thresholds. The sex-adjusted and no-sex
SP-versus-YC intersections both contained seven genes; all seven were retained
with no gain or loss. The other two contrasts remained empty under the
composite criterion.

| Tissue | Genes compared | Spearman rho | Sign concordance | Maximum absolute log2FC difference |
|---|---:|---:|---:|---:|
| Cochlea | 12,437 | 0.9893 | 96.04% | 1.0699 |
| Inferior colliculus | 12,437 | 0.9855 | 95.73% | 2.2389 |

Sex counts in the primary data were:

| Tissue | YC | MA | MP | SP |
|---|---|---|---|---|
| Cochlea | 4F/5M | 10F/7M | 4F/5M | 4F/2M |
| Inferior colliculus | 4F/4M | 10F/7M | 4F/5M | 4F/1M |

The 4F/1M inferior-colliculus SP group limits the identifiability of sex
adjustment. The no-sex analysis supports directional stability but does not
remove this limitation.

## Intersection and Tyrobp

The primary SP-versus-YC analysis yielded seven genes that met the frozen rule
in both separately fitted tissues, all higher in both tissues:

`Cd68`, `Clec7a`, `Fcgr3`, `Lyz1`, `Lyz2`, `Mpeg1`, and `Ptprc`.

Excluding GSM1201568 retained all seven and added `Tyrobp`, producing an
eight-gene intersection. Excluding GSM1201602 produced five genes, with
`Lyz2` and `Ptprc` no longer meeting the rule. Excluding both unmatched records
produced seven genes, with `Ptprc` absent and `Tyrobp` added. Across 41
leave-one-animal-out models, the intersection contained 5-10 genes (median 7;
distribution 5:1, 6:5, 7:27, 8:7, 10:1). `Tyrobp` was present in 6 of 41
models. The exact tables are in `shared_genes.csv` under each scenario
directory and in `results/leave_one_animal_out/`.

## GSE49543 probe-aggregation sensitivity

The primary GSE49543 preparation retains the first probe for each gene symbol.
The probe-aggregation sensitivity analysis rebuilt the cochlear expression
matrix under three
prespecified aggregation rules while holding the GSE49522 matrix, the
12,437-gene analysis universe, the limma model, and the significance criteria
fixed.

| Probe rule | Cochlea significant up/down | Intersection | Primary seven retained | Additional intersection genes |
|---|---:|---:|---:|---|
| First probe | 16/2 | 7 | 7 | None |
| Mean probe | 14/2 | 8 | 7 | `C1qb` |
| Maximum-mean probe | 17/4 | 9 | 7 | `Laptm5`, `Ms4a6d` |

All seven primary-intersection genes remained significant under every probe
rule. The complete intersection was not invariant in size or membership:
mean-probe aggregation added `C1qb`, and maximum-mean aggregation added
`Laptm5` and `Ms4a6d`. The public conclusion is therefore restricted to the
stable seven-gene core, not to a fixed complete eight- or nine-gene set.
Pairwise Spearman correlations of gene-level log2 fold changes ranged from
0.8360 to 0.9680. The input audit, three rebuilt matrices, contrast tables,
membership table, concordance table, and summary are under
`results/sensitivity_probe_aggregation/`.

## WGCNA rebuild

Both signed-hybrid networks were rebuilt from the top 5,000 variable genes in
the frozen common universe.

| Tissue | Samples | Power | Scale-free fit | Biological modules | Unassigned grey genes |
|---|---:|---:|---:|---:|---:|
| Cochlea | 41 | 7 | 0.8241421 | 23 | 109 |
| Inferior colliculus | 39 | 5 | 0.9424200 | 12 | 1,090 |

No biological module met descriptive adjusted P < 0.05 for MA versus YC or SP
versus MP. For SP versus YC, three cochlear modules (`MEgrey60`, `MEpurple`,
and `MEgreenyellow`) and one inferior-colliculus module (`MEgreen`) met the
module-level threshold.

## Single-cell audit

GSE274279 was reported by Xia et al. (2025, Cell Reports 44:115781) and was
generated from C57BL/6J mice. It provides cross-strain descriptive context for
the CBA/CaJ bulk reanalysis and is not validation in the bulk strain. The four
deposited libraries were GSM8446833 (3M), GSM8446834 (12M), GSM8446835 (24M),
and GSM8446836 (24M_Cochlea_BU). They retained 2,916, 5,532, 2,954, and 2,750
nuclei, respectively: 14,152 nuclei and 23,540 genes in total. GSM8446836 is a
second public 24-month cochlear library whose animal-level relationship to
GSM8446835 is not documented in GEO.

The atlas contained 245 marker-defined macrophage/microglial nuclei: 68 at 3
months, 81 at 12 months, and 96 at 24 months. Library totals were 68, 81, 50,
and 46. Six primary-intersection genes were measurable: Cd68, Fcgr3, Lyz1,
Lyz2, Mpeg1, and Ptprc. Clec7a was absent. Cd68, Fcgr3, Mpeg1, and Ptprc also
belonged to the marker panel used to define macrophages/microglia, so the
detection summary is not independent localization. Ptprc detection fractions
were 50.0%, 69.1%, and 47.9% at 3, 12, and 24 months.

Age, library, and 10x batch are inseparable. Age-related UMAP, composition,
detection, module, bootstrap, and subcluster patterns therefore describe the
sampled nuclei only. The public S11 summaries contain no cell-level
Mann-Whitney P values or BH-FDR estimates. The 500-iteration bootstrap samples
at most 25 nuclei per age per iteration; its intervals are not animal-level
biological-effect confidence intervals. CellChat uses only GSM8446833,
GSM8446834, and GSM8446835. The original GSE274279 10x matrices are not
redistributed.

## Environment and reproducibility

- R runtime and direct package versions: `environment/R_packages.csv`.
- Primary and sensitivity package records: `sessionInfo.txt` in each rebuild
  directory.
- Python package ranges and local runtime record: `environment/`.
- Repository-relative checksums: `results/` and
  `validation/primary_rebuild/input_checksums.csv` and
  `validation/sensitivity_rebuild/input_checksums.csv`.
- Document builder: `code/python/build_submission_docs.py`.
- Manifest builder: `code/python/build_manifest.py`.

## External context reproducibility

The archive includes the three deposited processed workbooks and code that
rebuilds the transcript-by-biological-sample matrices for GSE233798,
GSE153882, and GSE154833. The preparation manifest records source hashes,
selected columns, sample and run mappings, excluded rows, duplicate-symbol
handling, the log2(FPKM or RPKM + 1) transformation, and output hashes.

| Dataset | Transcripts | Biological samples | Excluded transcript rows |
|---|---:|---:|---:|
| GSE233798 | 34,395 | 6 | 8,228 |
| GSE153882 | 49,585 | 16 | 0 |
| GSE154833 | 48,709 | 9 | 0 |

The GSE154833 mapping records 36 GSM sequencing runs collapsed to the nine
submitted biological-sample columns; lanes were not treated as independent
replicates. No sample was excluded from any of the three external datasets.

`run_external_welch_bh.R` performs two-sided Welch tests with
Welch-Satterthwaite degrees of freedom and Benjamini-Hochberg correction across
all tested transcripts within each of the four dataset/context comparisons.
`validate_external_welch_bh.R` independently recomputed 28 candidate rows with
R `t.test`, found zero failures, and recomputed BH correction for all four
contexts. The largest absolute BH difference was 2.3315e-15.

`results/external/external_p_value_underflow.csv` retains 118 rows for which R
returned a numeric P value of zero: 32 in GSE233798 and 86 in GSE154833. These
values are labelled `numeric_underflow`, displayed as `<2.2e-16`, and retained
as returned for BH correction. They are not missing values or manually edited
zeros, and the complete transcriptome-wide tables preserve their provenance.
The external analyses remain descriptive tissue context and are not independent
replications of the primary SP-versus-YC contrast.

## Validation scope

The GEO, sample-audit, limma, shared-gene, WGCNA, and GSE49543
probe-aggregation layers were rebuilt from the archived inputs. The Seurat and
CellChat layers were audited against prepared source tables and
their archived scripts. The four public GSE274279 source matrices were checked
against the recorded source sizes and SHA-256 digests; the matrices themselves
and the large prepared RDS objects are not redistributed. Rebuilding the full
single-cell layers requires those public matrices and the R environment listed
under `environment/`.

## P-value display convention

Exact zero P values in public-facing text and display columns are reported as
`<2.2e-16`. Numeric calculation columns retain their original values so that
the supplied tables can be recomputed and checked. The `pval` columns in the
CellChat source tables remain numeric; `pval_display` is the reporting field.
The CellChat interpretation report uses `pval_display` rather than the numeric
field.
