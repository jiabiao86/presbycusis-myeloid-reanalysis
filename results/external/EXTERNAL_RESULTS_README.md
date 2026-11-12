# External Welch and BH results

This directory contains the complete transcriptome-wide Welch tests and
Benjamini-Hochberg adjusted values for the three external GEO datasets.

- GSE233798: old versus young whole cochlea; 3 versus 3 biological samples.
- GSE153882: 26-month versus 9-month IHC and OHC; 4 versus 4 biological samples each.
- GSE154833: 26-month versus 1-month stria vascularis; 3 versus 3 biological samples.

The submitters supplied normalized FPKM or RPKM values. The analysis applies
log2(x + 1), resolves duplicate symbols by the highest mean submitted value,
and performs Welch t tests with Welch-Satterthwaite degrees of freedom.
BH correction is applied across every tested transcript within each stated
dataset/context comparison. No external P value changes the primary analysis.

GSE154833 has 36 GSMs representing four sequencing lanes for each of nine
biological samples. The submitted All values worksheet already supplies nine
biological-sample columns (A/B/C at each age); lanes are not treated as
independent replicates. The 36-run mapping is retained in
inputs/external/prepared/run_mapping_external.csv.

When an extreme Welch statistic exceeds the floating-point resolution of the
student t distribution, R returns p = 0. Such rows are retained without
imputation, labelled p_value_status = numeric_underflow, and displayed as
<2.2e-16. BH correction uses the numeric p values as returned by R. The exact
affected rows are traceable in external_p_value_underflow.csv.

Files:
- *_all_genes.csv: complete results for each dataset/context.
- *_candidate_genes.csv: the seven primary intersection genes; absent genes
  are retained as gene_not_available rather than imputed.
- external_candidate_genes.csv: combined candidate table for Figure S2.
- external_welch_summary.csv: transcript counts and significance totals.
- external_p_value_underflow.csv: traceable rows with numeric p = 0.
- external_welch_validation*.csv and .json: independent R t.test checks for
  all candidate genes and p.adjust checks across every tested transcript.
