# Result index: manuscript items to archived files

## Figures

| Item | Source files |
| --- | --- |
| Figure 1. Study design and quality control | `results/tables/GSE49543_sample_metadata.csv`, `results/tables/limma_summary.csv` |
| Figure 2. Cochlear severity signature | `results/tables/limma_*.csv`, `results/tables/limma_severity_trend.csv`, `results/tables/severity_signature_scores.csv` |
| Figure 3. WGCNA module M12 and preservation | `results/tables/module_assignment.csv`, `results/tables/module_preservation_summary.csv` |
| Figure 4. Peripheral-central convergence | `results/tables/peripheral_central_candidate_comparison.csv`, `results/tables/candidate_panel_prioritization.csv`, `results/tables/candidate_genes_discovery_validation.csv` |
| Figure 5. Myeloid localization | `results/tables/mouse_immune_scores_GSE49543.csv`, `results/tables/mouse_immune_scores_GSE49522.csv`, `results/tables/GSE274279_candidate_expression_by_cell_type.csv`, `results/tables/GSE274279_cell_type_counts.csv` |
| Figure 6. External context and model performance | `results/tables/mouse_immune_deconvolution_external.csv`, `results/tables/nested_cv_hearing_loss_vs_normal_folds.csv`, `results/tables/nested_cv_severe_vs_rest_folds.csv`, `results/tables/candidate_genes_discovery_validation.csv` |
| Figure 7. Severity-focused CellChat and enrichment | `results/cellchat_formal/*.csv`, `results/supplementary/S21-*` to `S24-*` |
| Figure 8. Sensitivity analyses | `results/tables/cochlea_severity_*.csv`, `results/tables/central_severity_*.csv`, `results/tables/candidate_models_sensitivity.csv`, `results/tables/candidate_severity_bootstrap.csv` |
| Figure 9. Myeloid subclustering | `results/tables/myeloid_*.csv` |
| Figure 10. Exhaustive CellChatDB.mouse analysis | `results/cellchat_full/final/*`, `results/cellchat_full/analysis/*` |
| Graphical abstract | `manuscript/Graphical_Abstract/graphical_abstract.{pdf,png,tiff}` |

## Tables

| Item | Source files |
| --- | --- |
| Table 1. Datasets and roles | `results/tables/GSE49543_sample_metadata.csv` |
| Table 2. Severity model output | `results/tables/limma_summary.csv`, `results/tables/limma_severity_trend.csv` |
| Table 3. M12 module content | `results/tables/module_assignment.csv` |
| Table 4. Shared peripheral-central genes | `results/tables/peripheral_central_candidate_comparison.csv`, `results/tables/candidate_panel_prioritization.csv` |
| Table 5. Deconvolution and localization summary | `results/tables/mouse_immune_deconvolution_severity.csv`, `results/tables/mouse_immune_scores_*.csv` |
| Table 6. Strengths and limitations | Not data-derived |

## Supplementary tables

Supplementary Tables S1-S37 are provided in three forms:

- `results/supplementary/Supplementary_Tables_S1-S37.docx` - captions, column
  definitions, and the table index.
- `results/supplementary/Supplementary_Tables_S35-S37_Summary.xlsx` - index and
  previews for the exhaustive CellChat tables.
- `results/supplementary/*.csv` - complete machine-readable tables, including
  the full CellChat network and pathway results for all three ages.

## Exhaustive CellChat outputs

| File | Content |
| --- | --- |
| `results/cellchat_full/final/CellChat_full_GSE274279_{3,12,24}M.rds` | Final CellChat objects, 7 x 7 x 2,019 probability tensors |
| `results/cellchat_full/final/CellChat_full_network_{age}M.csv` | All non-zero ligand-receptor edges |
| `results/cellchat_full/final/CellChat_full_network_significant_{age}M.csv` | Edges with permutation p < 0.05 |
| `results/cellchat_full/final/CellChat_full_pathway_{age}M.csv` | Pathway-level communication |
| `results/cellchat_full/final/CellChat_full_summary_{age}M.csv` | Per-age run summary |
| `results/cellchat_full/analysis/*` | Cross-age comparison, pathway strength, cell-type strength, interpretation report |
| `results/cellchat_full/chunks/*.rds` | Raw per-chunk probability and p-value arrays |
