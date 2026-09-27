# Random seeds

| Analysis | Seed | Where it is set |
| --- | --- | --- |
| Global R seed for all analyses | 2026 | `set.seed(2026)` in `limma_analysis.R`, `wgcna_analysis.R`, `sensitivity_analysis.R`, `myeloid_subclustering.R` |
| CellChat bootstrap permutations (severity-focused and exhaustive runs) | 1 | `seed.use = 1L` in `cellchat_full_chunk.R` and the focused analysis |
| Module preservation (k-means, 100 initialisations) | 2026 | `KMeans(..., random_state=2026)` in `module_preservation.py` |
| Nested cross-validation (outer and inner folds, classifiers) | 2026 | `random_state=2026` in `nested_cross_validation.py` |
| Single-nucleus PCA, k-means (16 clusters, 50 initialisations), UMAP | 2026 | `random_state=2026` in `single_cell_localization.py` |
| Myeloid subclustering bootstrap (500 equal-cell samples) | 2026 | `set.seed(2026)` in `myeloid_subclustering.R` |
| Candidate-gene severity bootstrap (500 samples) | 2026 | `set.seed(2026)` in `sensitivity_analysis.R` |

CellChat draws its permutation matrix from the number of cells and `nboot` only,
so fixing `seed.use = 1` for every chunk guarantees that all chunks use the same
bootstrap permutations and that the merged network equals a single-pass run.
