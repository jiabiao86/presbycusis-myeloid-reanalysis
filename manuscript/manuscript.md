# Peripheral-Central Convergence of Macrophage and Microglial Programs in Presbycusis Severity

**Manuscript type:** Original research

**Short title:** Peripheral-central myeloid programs in presbycusis

**Authors:** Jiabiao Ji^1, Xiaoqing Yu^1, Chao Zhang^1, Long Liu^1, Lin Yan^1, Hongjian Zhang^1, Jianming Yang^1,*

**Affiliations:** ^1 The Second Affiliated Hospital of Anhui Medical University, Hefei, Anhui, China

**Corresponding author:** Jianming Yang, The Second Affiliated Hospital of Anhui Medical University, Hefei, Anhui, China (corresponding email to be confirmed before submission)

**Keywords:** presbycusis; age-related hearing loss; cochlea; inferior colliculus; macrophage; microglia; WGCNA; transcriptomics

## Abstract

**Background:** Presbycusis, or age-related hearing loss, is clinically heterogeneous and involves both peripheral and central auditory dysfunction. Previous transcriptomic studies have mainly analyzed the cochlea, often using single-tissue differential expression and hub-gene strategies. Whether severe hearing loss is associated with shared myeloid programs across the peripheral and central auditory system remains unclear.

**Methods:** We reanalyzed two Affymetrix Mouse Genome 430A microarray datasets from CBA/CaJ mice: GSE49543 (cochlea; 41 samples: 9 young controls, 17 middle-aged, 9 mild presbycusis, and 6 severe presbycusis) and GSE49522 (inferior colliculus; 39 samples: 8 young controls, 17 middle-aged, 9 mild presbycusis, and 5 severe presbycusis). RMA-normalized expression was modeled with limma/eBayes using ordinal severity, pairwise contrasts, and Benjamini-Hochberg correction. Weighted gene co-expression network analysis identified severity-associated modules, followed by module preservation analysis. Peripheral-central convergence was assessed by gene-level correlation and shared significance. Marker-based immune enrichment was used to estimate relative macrophage, microglial, and leukocyte program abundance without interpreting values as absolute cell proportions. Candidate programs were localized using GSE274279 single-nucleus RNA sequencing and evaluated in GSE233798, GSE153882, and GSE154833. Nested cross-validation was used only as a secondary computational analysis.

**Results:** Middle age alone produced no genes meeting FDR < 0.05 and absolute log2 fold change >= 1. In contrast, myeloid and antigen-presentation genes increased during progression from mild to severe presbycusis. The severity model identified 181 genes at FDR < 0.05, including 122 positively and 59 negatively associated genes. WGCNA identified module M12 as the module most strongly correlated with severity (r = 0.860, 66 genes), and M12 was preserved in the inferior colliculus (Zsummary = 8.40, permutation P = 2.28e-21). Among 12,437 genes shared by the cochlea and inferior colliculus, 6,551 changed in the same direction and 48 were significant in both tissues. The shared genes were enriched for immune system processes, macrophages, and microglia. Marker-based enrichment linked macrophage/microglial programs to severity in both the cochlea (r = 0.847, FDR = 2.14e-22) and inferior colliculus (r = 0.854, FDR = 2.14e-22); a microglia-like program was specifically associated with central severity (r = 0.630, FDR = 4.87e-6). In 11,403 single nuclei from GSE274279, the candidate genes were predominantly localized to 147 macrophages/microglia rather than hair cells, supporting cells, fibrocytes, or spiral ganglion neurons. Exhaustive CellChat analysis of all 2,019 CellChatDB.mouse interactions revealed an age-related decline in cochlear communication (455, 360, and 245 significant edges at 3, 12, and 24 months), with 46 ligand-receptor pairs significant at every age; PTN signaling was the strongest pathway throughout, macrophage/microglial PTPRC-MRC1 self-communication was significant at all ages, and APP-CD74 signaling emerged at 12 and 24 months. Transcription-factor enrichment converged on macrophage-associated regulators including IRF8, SPI1, CEBPB, SMRT, and NCOR, while second-pass enrichment supported phagocytosis, cytokine signaling, and antigen presentation. External datasets showed partial directional reproducibility and substantial tissue-specific divergence. Nested cross-validation produced internal areas under the curve of 1.000 for hearing loss versus normal and 0.869 for severe versus non-severe classification, whereas external transfer to GSE233798 yielded areas under the curve of 0.44 to 0.50.

**Conclusion:** Severity-resolved transcriptomic analysis identified a convergent macrophage and microglial program in the aging mouse cochlea and inferior colliculus. The findings support a peripheral-central myeloid axis in presbycusis severity while showing that bulk immune signals are tissue- and model-dependent. The work is computational, and translation to humans will require experimental validation.

## Introduction

Presbycusis is a multifactorial disorder characterized by progressive hearing loss during aging, with contributions from cochlear degeneration, altered central auditory processing, vascular and metabolic dysfunction, and low-grade inflammation. Clinically, patients with similar chronological age can differ substantially in hearing thresholds and speech perception, suggesting that age alone does not capture the biological heterogeneity of presbycusis [1-3]. Recent work has documented chronic inflammatory activation in the aging cochlea and cochlear nucleus, together with oxidative stress, altered cytokine signaling, and transcriptional changes in immune, metabolic, and sensory-cell programs [18-21,24,31]. Transcriptomic studies have therefore sought molecular signatures associated with auditory aging and hearing loss.

Most public transcriptomic analyses of presbycusis have focused on the cochlea. Weighted gene co-expression network analysis and differential expression have repeatedly prioritized immune and myeloid genes, including Mpeg1, Cd68, Fcgr3, Ctss, Clec4d, Ms4a7, and complement genes [4-7]. More recent single-cell studies have shown that immune cells are present in the aging cochlea and that macrophage subsets can acquire inflammatory and antigen-presenting features [2,6,8]. An immune component in presbycusis is therefore well supported, but several questions remain unresolved. First, many studies have treated hearing loss as a binary phenotype, even though mild and severe presbycusis may represent different stages of tissue remodeling. Second, immune changes detected in bulk cochlear tissue have not always been separated from changes in tissue composition, cell abundance, and regional degeneration. Third, analyses of the cochlea alone cannot determine whether similar myeloid programs emerge in central auditory structures.

The inferior colliculus is a major convergence center in the auditory pathway and shows age-related molecular and functional changes in animal models [9]. Structural and physiological studies have identified age-related ultrastructural, volumetric, inhibitory, and neurometric changes in the inferior colliculus, while microglial and interferon-related signaling have also been reported in the aging auditory pathway [32-37]. Peripheral and central auditory aging are not necessarily independent processes. Cochlear damage can alter afferent input and promote central plasticity or neuroinflammation, while central changes may influence auditory processing independently of cochlear threshold sensitivity. However, most presbycusis transcriptomic studies have not directly compared severity-associated programs between the cochlea and the inferior colliculus. A peripheral-central comparison could distinguish broadly shared immune programs from tissue-specific responses and would provide a more systems-level view of presbycusis severity.

We therefore integrated two severity-annotated mouse datasets from the same strain background and study family. GSE49543 contains cochlear microarray profiles, and GSE49522 contains inferior colliculus microarray profiles from CBA/CaJ mice classified as young control, middle-aged, mild presbycusis, or severe presbycusis. We used formal limma/eBayes modeling, ordinal severity regression, WGCNA, module preservation, peripheral-central gene-level comparison, marker-based immune enrichment, single-nucleus RNA sequencing localization, and external dataset evaluation. The analysis was designed to test three hypotheses: first, that severity progression is more informative than chronological age alone; second, that severe presbycusis is associated with shared macrophage and microglial programs across the cochlea and inferior colliculus; and third, that shared programs coexist with tissue-specific transcriptional changes that limit simple biomarker interpretation.

The resulting manuscript presents a computational reanalysis and does not include new animal experiments, qRT-PCR, immunostaining, Western blotting, or other laboratory validation. All inferences are therefore limited to transcriptomic associations and public-data reproducibility.

## Materials and Methods

### Study design and data sources

This study was a secondary analysis of public mouse transcriptomic data. The primary discovery dataset was GSE49543, an Affymetrix Mouse Genome 430A microarray study of the cochlea. The official GEO sample metadata were curated at the sample level. After correcting one severe presbycusis sample that was recorded as "Sever Presbycusis," the dataset contained 41 samples: 9 young controls, 17 middle-aged mice, 9 mice with mild presbycusis, and 6 mice with severe presbycusis. GSE49522, generated on the same microarray platform from the inferior colliculus, contained 39 samples: 8 young controls, 17 middle-aged mice, 9 mice with mild presbycusis, and 5 mice with severe presbycusis. The group labels were encoded for ordinal severity analysis as young control = 0, middle-aged = 1, mild presbycusis = 2, and severe presbycusis = 3.

Independent datasets were used for contextual evaluation rather than for remodeling the discovery model. GSE233798 provided RNA-seq FPKM profiles from aged and young mouse cochlea. GSE153882 provided RNA-seq RPKM profiles from inner hair cells and outer hair cells, and prior analyses of this dataset established that inner and outer hair cells undergo distinct biological aging trajectories [39]. GSE154833 provided RNA-seq RPKM profiles from the stria vascularis, a region whose contribution to age-related hearing loss has been debated [38]. GSE274279 provided 10x single-nucleus RNA sequencing profiles from mouse cochlea at 3, 12, and 24 months of age. Dataset roles and sample sizes are summarized in Table 1.

### Microarray preprocessing and quality control

GSE49543 and GSE49522 were Affymetrix Mouse Genome 430A microarray datasets, not RNA sequencing datasets. Probe-level intensities were normalized with the robust multi-array average method [10]. Expression values were analyzed on the log2 scale. Quality control included sample-level correlation summaries, principal component analysis, and comparison of within-group correlation distributions. No sample was removed solely on the basis of exploratory visualization. One severe sample was reassigned from the misspelled severity label to the severe presbycusis group.

### Differential expression and severity modeling

Differential expression was evaluated using limma with empirical Bayes moderated t-statistics [11]. Six pairwise contrasts were analyzed: middle-aged versus young control, mild presbycusis versus middle-aged, mild presbycusis versus young control, severe presbycusis versus middle-aged, severe presbycusis versus mild presbycusis, and severe presbycusis versus young control. Genes were reported at FDR < 0.05 with absolute log2 fold change >= 1 for the stringent result and absolute log2 fold change >= 0.585 for the moderate-effect result. Benjamini-Hochberg correction was used for multiple testing [12].

To model progression rather than only pairwise differences, ordinal severity was treated as a linear predictor in a limma model. This analysis estimated the direction and strength of monotonic expression change across young control, middle-aged, mild presbycusis, and severe presbycusis. The severity trend was analyzed separately in GSE49543 and GSE49522. Age and severity were not fully separable in the original designs, so severity models were interpreted as progression-associated transcriptomic trends rather than isolated causal effects of hearing loss.

### Weighted gene co-expression network analysis

WGCNA was performed using the 4,000 most variable genes in GSE49543 [13]. A soft-thresholding power of 16 was selected using the scale-free topology criterion. Modules were identified using default dynamic tree cutting followed by merging of highly similar modules. Module eigengenes were correlated with ordinal severity. The module with the strongest positive severity correlation, M12, was selected for downstream characterization.

Module preservation in GSE49522 was evaluated using WGCNA preservation statistics and permutation testing [14]. Preservation was summarized by Zsummary and empirical permutation P values. Modules with strong preservation were interpreted as reproducible co-expression structures, whereas modules without preservation were considered tissue-specific or model-dependent.

### Peripheral-central comparison

Gene-level severity associations from GSE49543 and GSE49522 were matched by gene symbol. The overall relationship between cochlear and central severity correlations was assessed using Pearson correlation. Genes with FDR < 0.05 in both tissues were defined as shared significant severity genes, and direction concordance was calculated among genes that were testable in both datasets. Shared significant genes were used for functional analysis with STRING [15]. Functional enrichment terms were considered supportive when they were coherent with the gene-level evidence and not dependent on a single isolated gene.

### Marker-based immune enrichment

Because human CIBERSORT signatures are not directly appropriate for mouse immune deconvolution, mouse-specific marker panels were used for relative marker-based enrichment scoring. Scores represented the relative abundance of immune transcriptional programs and were not interpreted as absolute cell fractions. The analysis included macrophage/microglia, microglia-like, monocyte, neutrophil, T cell, B cell, natural killer cell, dendritic cell, and mast cell programs. Correlations between immune scores and ordinal severity were calculated separately in GSE49543 and GSE49522. For reference, the previously generated human LM22 CIBERSORT results were compared with the mouse marker-based scores, but the human-signature analysis was not used as a primary result [16].

### Single-nucleus RNA sequencing localization

GSE274279 was processed with Scanpy [17]. Nuclei were filtered using standard quality-control metrics, normalized, log-transformed, and clustered. Cell types were assigned using marker scores for supporting cells, fibrocytes, glia/Schwann cells, inner hair cells, outer hair cells, spiral ganglion neurons, and macrophages/microglia. After quality control, 11,403 nuclei and 23,496 genes were retained. The age distribution was 2,916 nuclei at 3 months, 5,533 nuclei at 12 months, and 2,954 nuclei at 24 months. Macrophages/microglia represented 147 nuclei. Candidate gene localization was summarized by mean log-expression and the percentage of expressing nuclei within each cell type and age group. Age comparisons within the macrophage/microglial compartment were exploratory because of the limited number of nuclei and unequal age composition.

### External validation

Candidate genes and severity-associated programs were evaluated in GSE233798 for cochlear aging, GSE153882 for inner and outer hair cell context, and GSE154833 for stria vascularis context. Directional concordance, P values, and FDR values were calculated where sample structure permitted. Because external datasets differed in strain, model, RNA processing, tissue dissection, and age contrast, they were not treated as independent replications of the same clinical severity phenotype. Genes with concordant effects were classified as reproducible in direction, whereas genes with opposing effects or failed significance were retained as examples of model and tissue heterogeneity.

### Nested cross-validation

As a secondary analysis, nested cross-validation was performed for two classification tasks: hearing loss versus normal and severe presbycusis versus remaining samples. Feature selection and hyperparameter optimization were placed inside the training folds. Logistic regression and support vector machine models were evaluated. The final model was then transferred to GSE233798 as an external test. This analysis was designed to evaluate internal separability and external transfer limits, not to establish a diagnostic assay.

### CellChat and enhancement analyses

Cell-cell communication was evaluated with CellChat 1.6.1 [41]. All 2,019 ligand-receptor pairs of CellChatDB.mouse were evaluated. Database genes that were absent from the GSE274279 matrix were represented as all-zero rows, and two database genes with non-official mouse symbols were restored after gene filtering, so that every database interaction could be tested while complexes were still evaluated from their subunits. Because CellChat dispatches a parallel job for each coreceptor evaluation, the database was analysed in chunks of 505 interactions under a sequential evaluation plan, and the chunks were merged into a single network per age; the permutation seed was fixed at 1 so that every chunk drew the same bootstrap permutations. CellChat was run separately at 3, 12, and 24 months using normalized expression, `triMean` group averages, 100 bootstrap permutations, and a minimum of 10 cells per group. Interactions whose ligand or receptor was undetected return zero probability and a permutation P value of 1, and a severity-focused sensitivity run restricted to 106 curated interactions was retained for comparison. Significant communications were obtained with `subsetCommunication`, and the full interaction probability tensor was exported for age comparison.

An exploratory CellChatDB-informed ligand-receptor score was also calculated to evaluate broader interaction potential beyond the significant CellChat result. This analysis used the same database and complex definitions but substituted a product score based on mean ligand expression, mean receptor expression, and the percentages of expressing sender and receiver cells. These exploratory values represent relative communication potential and are reported separately from the formal CellChat inference.

### Sensitivity and myeloid subclustering analyses

Sex and hearing status were parsed from the GEO series metadata. For each tissue, severity effects were recalculated using limma with either sex adjustment or age-group adjustment. The age-group-adjusted model included severity as a continuous variable and an indicator for old mice with hearing loss. Candidate-gene severity coefficients were also estimated by 500 within-group bootstrap samples.

GSE274279 myeloid nuclei were extracted and reanalyzed with Seurat. Counts were normalized, variable genes were selected, cell-cycle-independent technical effects were controlled by regressing nCount_RNA and mitochondrial percentage, and PCA, nearest-neighbor graph construction, clustering at resolution 0.4, and UMAP were performed. Module scores were calculated for homeostatic, MHC-II, complement, phagocytic, interferon, APP-CD74, and PTPRC-MRC1 programs. To evaluate sensitivity to unequal age counts, 500 bootstrap samples of 25 myeloid nuclei per age were used to estimate age-related module-score changes.

Transcription-factor and miRNA enrichment were performed with Enrichr [42,43] using the ChEA 2022 and TargetScan microRNA 2017 libraries [44,45]. The M12 module and the positive severity gene set were analyzed separately. Drug-signature enrichment used DSigDB [46]. A second functional enrichment analysis used GO Biological Process 2025 and KEGG mouse libraries. Enrichment P values and adjusted P values were reported without additional post hoc filtering. Mouse symbols were converted to uppercase human-style symbols for Enrichr compatibility.

### Statistical analysis and reproducibility

All analyses used public data. P values were adjusted using the Benjamini-Hochberg procedure unless otherwise stated. Correlations were Pearson correlations. Module preservation P values were empirical permutation values. Analyses were performed with R 4.6.1, limma 3.68.5, WGCNA 1.74, Python 3.12, Scanpy 1.12.4, NumPy 2.5.3, pandas 3.0.6, SciPy 1.18.1, scikit-learn 1.9.1, and statsmodels 0.15.0. The random seed was fixed at 2026. Analysis scripts and frozen result tables are listed in the Code Availability section.

## Results

### Cohort structure and quality control

The curated GSE49543 cochlear dataset contained 41 mice. Sample-level correlations were consistently high within groups: mean within-group correlations were 0.972 for young controls, 0.979 for middle-aged mice, 0.979 for mild presbycusis, and 0.976 for severe presbycusis. Principal component analysis did not identify a single sample that justified removal before modeling. The corrected severe group included six samples after merging the misspelled severity label. The central GSE49522 dataset contained 39 mice and was analyzed with the same ordinal severity scale. Complete sample metadata and severity coding are provided in Supplementary Table S1.

The first 10 principal components of GSE49543 explained 19.98%, 14.90%, 7.50%, 6.97%, 5.57%, 5.04%, 4.50%, 3.81%, 2.72%, and 2.51% of variance, respectively. The absence of a dominant first component indicated that group differences were distributed across genes rather than driven by a single global batch effect. These quality-control results supported the use of the full curated dataset for severity modeling.

### Age alone did not produce robust differential expression

At FDR < 0.05 and absolute log2 fold change >= 1, no genes were significant for middle-aged versus young control mice (Table 2). The same threshold produced 9 upregulated and no downregulated genes for mild presbycusis versus middle-aged mice, 11 upregulated and 1 downregulated gene for mild presbycusis versus young controls, 12 upregulated and 3 downregulated genes for severe presbycusis versus middle-aged mice, no genes for severe versus mild presbycusis, and 16 upregulated and 2 downregulated genes for severe presbycusis versus young controls. Count summaries and complete differential expression results are provided in Supplementary Tables S2 and S3.

At the moderate effect threshold of absolute log2 fold change >= 0.585, the same pattern was retained: no genes for middle-aged versus young control, 17 upregulated and no downregulated genes for mild presbycusis versus middle-aged mice, 24 upregulated and 1 downregulated gene for mild presbycusis versus young controls, 29 upregulated and 9 downregulated genes for severe presbycusis versus middle-aged mice, no genes for severe versus mild presbycusis, and 35 upregulated and 12 downregulated genes for severe presbycusis versus young controls. The large number of changes involving the middle-aged-to-mild or middle-aged-to-severe transitions contrasted with the absence of a robust middle-aged versus young effect.

The ordinal severity model identified 181 genes at FDR < 0.05, including 122 genes with positive severity associations and 59 genes with negative severity associations. The strongest positive genes included Ctss, Mpeg1, Ms4a7, Clec4d, Fcgr3, Clec4n, Clec7a, Cd68, Ms4a6d, Il2rg, C1qc, Tyrobp, C3ar1, Fcer1g, Csf1r, H2-Aa, and H2-K1. These genes represent lysosomal proteolysis, phagocytosis, Fc receptor signaling, complement activation, MHC class II antigen presentation, and myeloid immune regulation. The strongest negative genes included Uros, Fgf1, Lynx1, Naalad2, Ttl, Amz2, Otor, Cyth3, Rab3gap2, Abhd10, Hdgfrp3, Aass, Gdpd1, Pfkm, Pfkl, Ankrd40, and Slc25a14. These genes were more heterogeneous and included metabolic, mitochondrial, cytoskeletal, and cochlear maintenance functions rather than a single pathway. The complete cochlear severity trend is provided in Supplementary Table S4.

The positive severity signature therefore represented a coherent myeloid program, while the negative signature represented a broader decline in metabolic and structural programs. This directional imbalance was consistent with previous reports of inflammation-related and mitochondrial gene modules in age-related hearing loss [28,29], and it informed the subsequent emphasis on the myeloid axis while preventing overinterpretation of the negative signature as one unified pathway.

### WGCNA identified a preserved severity-associated myeloid module

WGCNA grouped 4,000 high-variance genes into 16 modules. Module M12 had the strongest correlation with ordinal severity (r = 0.860) (Table 3). M12 contained 66 genes, including Mpeg1, Cd68, Clec7a, Fcgr3, H2-Aa, Ctss, Clec4d, Tyrobp, Csf1r, C1qc, Ms4a7, Lgals3, Cxcl13, Lyz1, Lyz2, Cd14, Fcer1g, C3ar1, Cd84, and Lcp1. The module therefore combined macrophage identity, phagocytic function, lysosomal activity, antigen presentation, complement signaling, and inflammatory chemokines. Module assignments, trait correlations, and preservation statistics are provided in Supplementary Tables S6 and S7.

Preservation analysis in the inferior colliculus supported M12 as a reproducible co-expression structure. M12 had Zsummary = 8.40 and permutation P = 2.28e-21. This result was not simply a consequence of one highly connected hub gene, because M12 contained multiple related myeloid markers and was identified independently of the marker-based immune scoring analysis.

Other modules had different relationships to severity. Module M9 showed a moderate positive correlation (r = 0.437), whereas M5, M11, and M3 showed negative correlations (r = -0.429, -0.471, and -0.162, respectively). The dominance of M12 at severe stages, together with the preservation of its co-expression structure in the central dataset, provided the first evidence for a shared myeloid program across peripheral and central auditory tissues.

### The cochlea and inferior colliculus shared a myeloid-enriched severity signature

GSE49543 and GSE49522 contained 12,437 genes that could be matched by symbol. The correlation between cochlear and central severity effect estimates was positive but modest (r = 0.146, P = 3.03e-61). This small genome-wide correlation indicated that the inferior colliculus was not a simple mirror of the cochlea; most genes either changed in only one tissue or changed in opposing directions. Nevertheless, 6,551 genes changed in the same direction in both tissues, and 48 genes reached FDR < 0.05 in both.

The 48 shared significant genes were strongly enriched for myeloid immune functions. Representative genes included Mpeg1, C1qa, C1qb, C1qc, Clec7a, Ptprc, Fcgr3, Tyrobp, Cd68, Ctss, Cd74, H2-Aa, H2-K1, Cd14, C3ar1, Csf1r, Lgals3, and Icam1 (Table 4). STRING analysis identified immune system process as the leading functional category (FDR = 2.73e-20), followed by macrophage tissue enrichment (FDR = 6.95e-20) and microglia tissue enrichment (FDR = 4.72e-19). Additional enriched functions included antigen presentation, phagocytosis, Fc-gamma receptor signaling, complement activation, and microglial pathogen phagocytosis. The complete gene-level peripheral-central comparison is provided in Supplementary Table S11.

The shared gene list contained both positive and negative severity associations. Uros and Lynx1 were shared significant genes with negative associations, indicating that the peripheral-central commonality was not exclusively a myeloid increase. However, the immune terms were driven primarily by the larger positively associated group, and the strongest shared positive genes were Ctss, Mpeg1, Fcgr3, Clec7a, Cd68, Ms4a6d, Cd84, Il2rg, C3ar1, C1qc, Tyrobp, and H2-K1. The comparison defined a peripheral-central convergent module rather than a single tissue-specific inflammatory response.

### Marker-based enrichment linked severity to macrophage and microglial programs

In the cochlea, macrophage/microglia marker enrichment had the strongest correlation with ordinal severity (r = 0.847, FDR = 2.14e-22). Monocyte (r = 0.536, FDR = 3.34e-4), neutrophil (r = 0.499, FDR = 1.15e-3), T cell (r = 0.424, FDR = 0.0103), and dendritic cell (r = 0.362, FDR = 0.0397) programs also showed weaker positive associations. B cell, natural killer cell, and mast cell programs did not reach FDR < 0.05.

In the inferior colliculus, macrophage/microglia enrichment was similarly associated with severity (r = 0.854, FDR = 2.14e-22). A microglia-like marker program was also significant (r = 0.630, FDR = 4.87e-6). Monocytes, neutrophils, T cells, B cells, natural killer cells, dendritic cells, and mast cells did not show significant severity associations in the central dataset. The contrast suggested that central severity was dominated by a resident macrophage/microglial-like program, whereas the cochlea showed contributions from multiple circulating and tissue-resident myeloid compartments. The complete enrichment analysis is provided in Supplementary Table S10.

The marker-based analysis estimated relative program abundance and not absolute cell proportions. The very similar macrophage/microglia correlations in the two tissues should therefore be interpreted as evidence of shared myeloid transcriptional enrichment, not as proof that the same cell population expands by the same amount in both tissues.

### Single-nucleus RNA sequencing localized candidate genes to macrophages and microglia

After quality control, GSE274279 contained 11,403 nuclei across 3-, 12-, and 24-month-old mice. The dominant cell types were supporting cells (n = 8,387), fibrocytes (n = 1,332), glia/Schwann cells (n = 865), outer hair cells (n = 283), inner hair cells (n = 261), macrophages/microglia (n = 147), and spiral ganglion neurons (n = 128). Macrophages/microglia were captured in two clusters, which is consistent with heterogeneity within the myeloid compartment or with differences in cell state. Cell-type counts and candidate expression summaries are provided in Supplementary Tables S8 and S9.

Candidate genes showed strong myeloid localization. Cd74, H2-Aa, H2-Eb1, Ctss, Cd68, Fcgr3, Tyrobp, C1qa, Mpeg1, and Lgals3 were expressed predominantly in macrophages/microglia. For example, Ctss had mean log-expression of 0.742, 1.460, and 0.999 in macrophages/microglia at 3, 12, and 24 months, compared with mean values generally below 0.04 in supporting cells, fibrocytes, hair cells, and spiral ganglion neurons. Tyrobp had mean log-expression of 0.346, 0.512, and 0.091 in macrophages/microglia, with near-zero expression in most other cell types. Cd74 showed a similar myeloid-restricted pattern, with mean log-expression of 0.461, 1.474, and 0.538 in macrophages/microglia across the three age groups.

The single-nucleus data therefore supported the cellular source of the bulk immune signal. The bulk severity signature was unlikely to reflect intrinsic expression of the same myeloid genes by hair cells or supporting cells. Instead, it was more consistent with myeloid cell abundance, state change, or both. Age-direction estimates within macrophages/microglia were variable: Apoe, C1qa, and C1qc were lower at 24 months than at 3 months, whereas several inflammatory genes were highest at 12 months. Because the number of macrophage/microglial nuclei was limited and the age groups were unbalanced, these within-cell-type age comparisons were treated as exploratory and were not used to define a linear aging trajectory.

### External datasets showed partial reproducibility and strong tissue context

GSE233798 provided a different aged mouse cochlear model. Candidate genes H2-Aa and H2-Eb1 were increased in aged cochlea with FDR < 0.05, and Setd1a, Nfatc2, Cd74, and Apoe also showed significant positive effects. Several myeloid genes had concordant positive directions without FDR significance, including Fcgr3, Cd68, C1qa, C1qb, C3ar1, Ctss, and Ms4a7. In contrast, Tyrobp, Mpeg1, Clec4d, Laptm5, Mif, and Apoe showed directionally inconsistent or significant opposite effects across the discovery and validation contexts. This heterogeneity indicated that not every severe-presbycusis-associated myeloid gene is reproducible across aging models. Complete external validation results are provided in Supplementary Tables S12-S14, and prioritized candidate evidence is provided in Supplementary Table S15.

GSE153882 separated inner hair cells from outer hair cells. In inner hair cells, many myeloid-associated genes were increased with age, including Ctss, Fcgr3, Cd68, Mpeg1, Tyrobp, C1qa, C3ar1, Ms4a7, Clec7a, Setd1a, and H2-Aa. In outer hair cells, the same genes were often increased more strongly, but the interpretation was limited by small group sizes and the possibility of ambient immune signal or altered cell composition. GSE153882 therefore supported the presence of immune-associated transcripts in hair-cell preparations but did not establish cell-autonomous expression.

GSE154833 provided a distinct regional context. In the stria vascularis, many candidate immune genes decreased with age, including Cd74 (log2 fold change = -2.43), Lgals3 (-3.61), H2-Aa (-2.26), Mpeg1 (-2.20), Ctss (-1.82), Cd68 (-1.76), Fcgr3 (-1.44), and Tyrobp (-0.96). This direction was opposite to the average effects in whole cochlea and inferior colliculus. The discrepancy is biologically plausible if bulk cochlear immune enrichment reflects changes in myeloid abundance or infiltration, while stria vascularis aging involves atrophy or altered regional composition. It is also consistent with evidence that sensory-epithelial cells, rather than the stria vascularis alone, carry a major component of age-related genetic risk and with single-cell studies showing region-specific stria vascularis alterations [31,38]. It argued against describing all immune genes as uniformly upregulated across every cochlear region.

Across these datasets, the most reproducible program-level finding was not a single universal biomarker but the recurrent involvement of MHC class II, complement, Fc receptor, lysosomal, and phagocytic genes. Candidate-level reproducibility varied by model and tissue.

### Nested cross-validation revealed internal separability but limited external transfer

Nested cross-validation achieved high internal performance in GSE49543. Logistic regression produced an area under the curve of 1.000 for hearing loss versus normal and 0.869 for severe presbycusis versus remaining samples. Support vector machines produced areas under the curve of 0.997 and 0.857, respectively. Stable features included C3ar1, Cxcl13, Csf1r, Ctss, C1qc, Clec4d, Fcgr3, Tyrobp, Ms4a7, Mpeg1, and Cd68, which aligned with the severity and module analyses. Model summaries, fold-level performance, and feature stability are provided in Supplementary Tables S16-S18.

External transfer of the hearing-loss classifier to GSE233798 produced areas under the curve of 0.444 for logistic regression and 0.500 for support vector machines. The marked decrease in performance indicated that internal discrimination depended strongly on the original study design, group structure, platform, and tissue context. The machine-learning results were therefore retained only as a secondary computational observation. They do not support the claim that the current features form a generalizable diagnostic model.

### CellChat and enhancement analyses identified myeloid communication and regulatory programs

To place the severity-focused analysis in context, we first evaluated the complete CellChatDB.mouse database. All 2,019 ligand-receptor pairs were tested in every age group with 100 bootstrap permutations. The number of significant cell-type-resolved edges decreased monotonically with age, from 455 at 3 months to 360 at 12 months and 245 at 24 months, corresponding to 146, 115, and 64 significant ligand-receptor interactions and to 33, 34, and 27 significant pathways; the total significant communication probability declined in parallel (11.49, 9.39, and 7.44). Forty-six interactions were significant at all three ages, no interaction was significant at 3 months alone, and two were significant at 24 months alone, indicating that presbycusis involves a quantitative contraction of the cochlear communication network rather than a switch between unrelated programs (Figure 10 and Supplementary Tables S35-S37).

The preserved core was dominated by extracellular-matrix and adhesion signaling (laminin, collagen, fibronectin, and SPP1 integrin interactions; CADM1, NCAM1, MPZ, and PTPRM), neurotrophic and guidance pathways (NRG3-ERBB4, NRXN1/3-NLGN1, EFNA5-EPHA3/4/5/7, SEMA3A-NRP1-PLXNA1/2/4, and NTF3-NTRK2/3), and growth-factor signaling (BMP5/6, TGFB2, FGF1, and IGF1). PTN signaling was the strongest pathway at every age, with PTN-PTPRZ1, PTN-ALK, and PTN-SDC2 contributing the largest single-edge probabilities (up to 0.32).

Macrophage/microglial communication behaved differently. PTPRC-MRC1 self-communication was significant at all three ages, with probabilities of 0.064, 0.088, and 0.045, reproducing the result of the severity-focused run (Supplementary Table S25). APP-CD74 signaling, which the exploratory score had highlighted, reached significance for seven myeloid-recipient edges at 12 and 24 months but not at 3 months (maximum probability 0.152 at 12 months and 0.038 at 24 months), and the APP pathway was among the strongest pathways at 12 and 24 months (0.969 and 0.220). CD74-dependent myeloid signaling therefore behaves as a feature acquired during presbycusis progression rather than as a constitutive feature of the young cochlea.

In addition to the exhaustive database run, we calculated graded CellChatDB-informed interaction scores, which are not thresholded by permutation testing. This analysis produced 115,401 ligand-receptor scores across the three ages (Supplementary Table S19). The largest age-related changes involving macrophages/microglia included APP-CD74 interactions from inner hair cells, supporting cells, fibrocytes, and glia/Schwann cells to macrophages/microglia. The APP-CD74 score increased from 0.903 to 1.629 for inner hair cells to macrophages/microglia, from 0.944 to 1.523 for supporting cells to macrophages/microglia, and from 0.867 to 1.438 for fibrocytes to macrophages/microglia. Exploratory myeloid interactions involving PTN-PTPRZ1, SEMA3A-NRP1/PLXNA2, and LAMA2-DAG1 also increased with age. Such scores are hypothesis-generating and are not equivalent to significant CellChat interactions. Complete exploratory scores and age-related changes are provided in Supplementary Tables S19 and S20.

Transcription-factor enrichment of the severity-positive gene set showed strong macrophage-associated programs. The leading terms included SMRT/NCOR and IRF8 ChIP-seq signatures from macrophages or bone marrow-derived macrophages, SPI1, CEBPB, and MECOM. The M12 module produced a similar pattern, with SMRT, IRF8, NCOR, Nerf2, RUNX1, CEBPB, and MECOM among the top enriched regulators. The enriched regulators aligned with the myeloid identity of M12 and pointed to a macrophage-specific component rather than a generic inflammatory transcription program. Complete TF results are provided in Supplementary Table S21.

miRNA target enrichment did not identify any term meeting FDR < 0.05. The strongest nominal terms included miR-4800-3p, miR-4460, and miR-380-5p, but their adjusted P values were close to 1. The absence of a significant miRNA signal may reflect the small size of the M12 and severity-positive gene sets, the use of predicted target libraries, or the limited statistical power of the current design. The negative result is reported in Supplementary Table S22 and is not interpreted as evidence of a specific miRNA mechanism.

The second functional enrichment analysis supported the primary STRING results. GO Biological Process 2025 enrichment of the severity-positive gene set identified positive regulation of phagocytosis, cytokine-mediated signaling, leukocyte cell-cell adhesion, granulocyte chemotaxis, phagocytosis, antifungal innate immune response, and regulation of superoxide anion generation. KEGG mouse enrichment identified phagosome, Staphylococcus aureus infection, viral myocarditis, cell adhesion molecules, antigen processing and presentation, and tuberculosis. The M12 module showed the same broad pattern, including cytokine-cytokine receptor interaction, phagosome, TNF signaling, and C-type lectin receptor signaling. The secondary enrichment analysis provides orthogonal functional support for the interpretation that the severity-associated program combines phagocytic, antigen-presenting, and cytokine-signaling functions (Supplementary Table S23).

Drug-signature enrichment identified multiple perturbation signatures associated with the severity-positive gene set, including mebendazole, nickel sulfate, methotrexate, phorbol 12-myristate 13-acetate, 1-nitropyrene, dexamethasone, and beclomethasone. The M12 module also matched signatures including 1-nitropyrene, phorbol 12-myristate 13-acetate, acetovanillone, mebendazole, and progesterone. The severity signature therefore has detectable pharmacological connectivity, but the analysis does not establish efficacy, causality, or clinical benefit. Drug and chemical terms are therefore reported only as computational repositioning hypotheses in Supplementary Table S24.

### Sensitivity analyses showed robust sex adjustment but incomplete separation from age

Sex adjustment had little effect on the cochlear severity signal. The correlation between severity-only and sex-adjusted effects was 0.997, and the number of genes at FDR < 0.05 changed from 181 to 189. In the inferior colliculus, the corresponding effect correlation was 0.855 and the number of significant genes changed from 228 to 134. Candidate genes retained their direction after sex adjustment in both tissues. Bootstrap confidence intervals for the leading cochlear candidate genes were above zero for all 18 evaluated genes, and the inferior colliculus also showed positive bootstrap intervals for H2-Aa, H2-Eb1, Cd74, Ctss, Fcgr3, Cd68, Tyrobp, Mpeg1, C1qc, C3ar1, and Clec7a. Sex composition therefore does not primarily explain the myeloid severity direction.

In contrast, age-group adjustment substantially attenuated the severity signal. After including an old-hearing-loss indicator, no cochlear genes remained significant at FDR < 0.05, and only Mpeg1 remained significant in the inferior colliculus. The correlation between severity-only and age-adjusted effects decreased to 0.541 in the cochlea and 0.507 in the inferior colliculus. This does not invalidate the ordinal severity model, because severity and age are intrinsically related in presbycusis, but it demonstrates that the current design cannot fully separate age-related and hearing-loss-related effects. The sensitivity results are provided in Supplementary Tables S26-S28 and Figure 8.

### Myeloid subclustering identified two states with age-dependent scores

Reclustering of the 147 macrophage/microglial nuclei identified two myeloid clusters. Myeloid_C0 contained 38 cells at 3 months, 54 at 12 months, and 19 at 24 months, corresponding to 71.7%, 88.5%, and 57.6% of myeloid nuclei, respectively. Myeloid_C1 increased from 28.3% at 3 months to 42.4% at 24 months. Myeloid_C0 had higher APP-CD74, PTPRC-MRC1, MHC-II, complement, and phagocytic scores than Myeloid_C1. Within Myeloid_C0, APP-CD74 and PTPRC-MRC1 scores peaked at 12 months and declined at 24 months. Equal-cell bootstrap analysis showed negative 24-month-minus-3-month changes for APP-CD74 (median = -0.437, 95% CI -0.826 to -0.097) and complement (median = -0.353, 95% CI -0.683 to -0.085). The myeloid states were therefore age-dependent, but the small number of nuclei limits definitive subtype assignment. The subclustering results are provided in Supplementary Tables S29-S34 and Figure 9.

## Discussion

This integrative reanalysis identified a severity-associated macrophage and microglial program shared by the aging mouse cochlea and inferior colliculus. Three observations support the main conclusion. First, middle age alone did not produce robust differential expression at stringent thresholds, whereas progression from middle-aged to mild or severe presbycusis was accompanied by a coherent myeloid signature. Second, WGCNA identified module M12 as strongly correlated with severity and reproducible in the inferior colliculus. Third, 48 genes were significant in both tissues and were enriched for macrophage, microglia, antigen presentation, phagocytosis, complement, and Fc receptor functions. Single-nucleus data localized representative genes to macrophages/microglia.

The severity-resolved design helps explain why previous single-tissue hub-gene studies could produce apparently inconsistent results. Many prior analyses compared aged with young animals, used binary hearing-loss labels, or prioritized network hubs without modeling clinical severity [4-7]. Subsequent studies have linked age-related hearing loss to inflammatory genes, mitochondrial dysfunction, altered m6A regulation, chromatin accessibility, and transcriptional reprogramming [25-30]. In our analysis, severe versus mild presbycusis produced no genes at the stringent threshold, whereas comparisons involving middle-aged mice produced the largest number of changes. This pattern suggests that the transition from middle age to clinically apparent hearing loss may involve nonlinear or threshold-like remodeling. It also suggests that severe presbycusis is not simply an amplified version of mild presbycusis, at least at the transcript level.

The peripheral-central convergence observed here extends the immune interpretation of presbycusis beyond the cochlea. The shared genes included canonical macrophage and microglial markers such as Mpeg1, C1qa, C1qb, C1qc, Cd68, Csf1r, Tyrobp, Fcgr3, Clec7a, and Ptprc. MHC class II and antigen-presentation genes, including Cd74, H2-Aa, and H2-K1, were also shared. The presence of complement and Fc receptor genes links the signature to phagocytic and immune-complex recognition pathways, while the lysosomal protease Ctss connects antigen processing to extracellular matrix and inflammatory signaling. These functions are consistent with activated or reprogrammed myeloid states rather than a single classical inflammatory cytokine.

The central signal should not be interpreted as proof that the inferior colliculus is the primary driver of presbycusis. The modest genome-wide correlation between cochlear and inferior colliculus severity effects (r = 0.146) shows that the two tissues are largely distinct in their transcriptomic responses. The 48 shared genes form a small intersecting core within a larger set of tissue-specific responses. The inferior colliculus may respond to altered peripheral input, undergo independent aging-related glial changes, or combine both processes [32-37]. Causal direction cannot be determined from cross-sectional transcriptomic data.

The tissue-specific findings are equally important. Stria vascularis data showed age-related decreases in several candidate immune genes, whereas the whole cochlea showed increases. Single-nucleus data localized the same genes to macrophages/microglia, with low expression in supporting cells, fibrocytes, hair cells, and spiral ganglion neurons. Together, these findings favor a composition- and cell-state-based explanation for bulk immune enrichment. Regional atrophy, vascular degeneration, altered cell recovery, or immune-cell localization could all contribute [31,38,39]. The data do not support a model in which every cochlear cell increases expression of macrophage genes.

The candidate genes form a prioritized but non-diagnostic panel. H2-Aa, H2-Eb1, Cd74, Setd1a, Ctss, Fcgr3, Cd68, Tyrobp, and Mpeg1 were supported by multiple lines of evidence, including severity association, network membership, peripheral-central sharing, and myeloid localization. However, external reproducibility was partial. H2-Aa and H2-Eb1 were reproducible in GSE233798, while Tyrobp, Mpeg1, Clec4d, and Laptm5 showed context-dependent behavior. These genes may mark biological processes that are common even when individual transcript changes vary. Previous work has implicated individual myeloid genes such as Ctss, Cd68, Fcgr3, Mpeg1, Apoe, and inflammatory signaling in cochlear aging, but no single transcript has yet emerged as a validated human target [6,7,22-24,29,40]. Current data are insufficient to recommend any gene as a human diagnostic or therapeutic target.

The enhancement analyses strengthen the myeloid interpretation in three ways. First, the exhaustive CellChat run over all 2,019 CellChatDB.mouse interactions showed that matrix, adhesion, and neurotrophic signaling dominates the cochlear communication network, that PTN signaling is the strongest pathway at every age, and that the amount of significant communication contracts with age. PTPRC-MRC1 macrophage/microglial communication was significant at all ages with a peak at 12 months, and APP-CD74 signaling became significant at 12 and 24 months but not at 3 months, which converts the earlier exploratory APP-CD74 observation into a formally tested interaction. CD74 is the MHC class II invariant chain and, together with APP, can participate in myeloid activation and antigen-presentation programs [41]. This links the shared H2-Aa, H2-Eb1, and Cd74 signature to candidate intercellular signaling rather than only to intracellular expression, while the age-dependent appearance of APP-CD74 and the declining number of significant interactions at 24 months place myeloid signaling within a broader contraction of cochlear communication. Second, TF enrichment converged on IRF8, SPI1, CEBPB, SMRT, and NCOR, which are consistent with macrophage lineage identity and myeloid enhancer regulation. Third, the GO and KEGG results independently recovered phagocytosis, antigen presentation, cytokine signaling, and phagosome biology. These analyses do not prove signaling activity, but they provide corroborating computational evidence for the myeloid module.

The drug-signature results should be interpreted carefully. Dexamethasone and beclomethasone connect the severity signature to glucocorticoid-related transcriptional responses, whereas methotrexate and several chemical perturbation terms reflect broader immune or stress signatures. None of these matches indicates that a listed drug treats presbycusis or should be used clinically. They identify hypotheses that could be tested experimentally or in appropriately designed pharmacological studies. Similarly, the absence of significant miRNA enrichment is an informative negative result and argues against forcing a miRNA mechanism from the current data.

The sensitivity analyses refine the severity interpretation. Sex adjustment preserved the myeloid direction and had minimal effect on the cochlear effect estimates, whereas age adjustment removed nearly all significant severity genes. Because chronological age and hearing phenotype are strongly linked in the original study design, the fully age-adjusted model is probably over-adjusted and should not be treated as the primary model. The comparison nevertheless disciplines the interpretation: the current analysis supports severity-associated myeloid remodeling, but it cannot establish how much of that signal is due to age alone versus hearing loss.

The myeloid subclustering analysis showed two broad states with different age distributions and module scores. Myeloid_C0 carried stronger APP-CD74, PTPRC-MRC1, MHC-II, complement, and phagocytic signatures, whereas Myeloid_C1 was less immunometabolically active. The apparent decline in APP-CD74 and complement scores at 24 months is consistent with the formal CellChat result, which also showed lower PTPRC-MRC1 probability at 24 months than at 12 months. However, this pattern should not be interpreted as a simple loss of myeloid activation. Differences in cell capture, state composition, and the small number of nuclei could contribute.

The machine-learning analysis provides a useful methodological caution. Nested cross-validation produced perfect internal separation for hearing loss versus normal and good internal discrimination for severe versus non-severe classification, but external transfer was near random. This discrepancy is common when feature selection, sample composition, and batch structure differ between datasets. The failure of external transfer reinforces the decision to treat machine learning as secondary and to avoid presenting the current model as a clinical classifier.

This study has several limitations. It is purely computational and contains no independent experimental validation. The discovery datasets used microarrays from a specific mouse strain and study design, and age, cohort, and hearing phenotype cannot be fully disentangled. Age-adjusted sensitivity analysis removed most significant severity genes, confirming the strong age-severity dependence. Marker-based immune scoring estimates relative program enrichment and is not a substitute for absolute cell counting. The single-nucleus dataset contained relatively few macrophage/microglial nuclei, and its age groups were unbalanced. Two broad myeloid states were detectable, but higher-resolution subtype assignment and lineage inference remain underpowered. The exhaustive CellChat run reports permutation P values without additional multiple-testing correction, so the number of significant edges should be read as a descriptive summary of communication density rather than as a family-wise error-controlled count, and all CellChat probabilities remain model-based estimates. These scores summarize transcript abundance and ligand-receptor database relationships; they are not direct measurements of protein complexes, ligand secretion, receptor activation, or signaling flux. TF and drug enrichments are database-based inferences, and the miRNA analysis was negative after multiple-testing correction. The external datasets differed in tissue region, RNA platform, quantification method, animal model, and age contrast. Human relevance remains indirect because the primary evidence comes from mice. Finally, the ordinal severity scale assumes an approximately monotonic ordering from young control to severe presbycusis, whereas the actual biology may include nonlinear or stage-specific transitions.

Despite these limitations, the analysis provides a reproducible computational framework for severity-resolved study of peripheral-central auditory aging. The principal contribution is not a single hub gene but a cross-tissue myeloid program that links macrophage and microglial biology to presbycusis severity. Future experimental work should determine which myeloid states are causal, which are compensatory, and how central glial responses relate temporally to cochlear degeneration.

## Conclusion

GSE49543 and GSE49522 support a peripheral-central convergence of macrophage and microglial programs in presbycusis severity. The strongest evidence consists of a severity-associated 66-gene WGCNA module, 48 genes shared between the cochlea and inferior colliculus, and myeloid-specific localization in single-nucleus data. Two broad myeloid states displayed different age distributions and immunometabolic scores. Sex-adjusted models preserved the severity direction, whereas age-adjusted models showed that age and hearing loss cannot be fully separated. An exhaustive CellChat analysis of all 2,019 CellChatDB.mouse interactions showed declining cochlear communication with age, a preserved core of 46 ligand-receptor pairs, PTPRC-MRC1 macrophage/microglial communication at all ages, and APP-CD74 signaling from 12 months onwards, while TF enrichment and second-pass functional analysis linked the module to macrophage-associated regulators, phagocytosis, and antigen presentation. The same analyses show substantial tissue specificity and incomplete external reproducibility. The analysis defines a biologically plausible, computationally derived myeloid axis in presbycusis while explicitly excluding claims of a validated diagnostic model or human biomarker.

## Data Availability

All datasets analyzed in this study are publicly available from the Gene Expression Omnibus under accessions GSE49543, GSE49522, GSE233798, GSE153882, GSE154833, and GSE274279. Processed summary tables, frozen analysis outputs, and the full CellChat result objects are archived at https://github.com/jiabiao86/presbycusis-myeloid-reanalysis.

## Code Availability

Analysis scripts include limma_analysis.R, wgcna_analysis.R, reanalysis.py, central_auditory_comparison.py, mouse_immune_deconvolution.py, module_preservation.py, single_cell_localization.py, nested_cross_validation.py, generate_main_figures.R, cellchat_formal_analysis.R, cellchat_enhancement_analysis.R, enrichr_enhancement_analysis.py, generate_enhancement_figure.R, sensitivity_analysis.R, generate_sensitivity_figure.R, myeloid_subclustering.R, generate_myeloid_figure.R, cellchat_full_cache.R, cellchat_full_chunk.R, cellchat_full_merge.R, run_cellchat_full.sh, run_cellchat_full_merge.sh, cellchat_full_compare.R, cellchat_full_analysis.R, and generate_cellchat_full_figure.R. The scripts use public GEO inputs, fixed random seeds where applicable, and the version information recorded in the analysis manifest. The archived code, environment record, and result tables are available at https://github.com/jiabiao86/presbycusis-myeloid-reanalysis.

## Supplementary Materials

Supplementary Tables S1-S37 are provided as Supplementary_Tables_S1-S37.docx, Supplementary_Tables_S1-S11.xlsx, Supplementary_Tables_S12-S18.xlsx, Supplementary_Tables_S19-S25_Summary.xlsx, Supplementary_Tables_S26-S34_Summary.xlsx, Supplementary_Tables_S35-S37_Summary.xlsx, and the complete CSV data files for Supplementary Tables S19-S37. The Word file contains the formal table titles, descriptions, column definitions, and index. The Excel workbooks and CSV files contain the complete machine-readable tables, with Supplementary Table S3 provided across six comparison-specific worksheets (S3A-S3F) and Supplementary Table S13 provided across two cell-type worksheets (S13A-S13B).

## Ethics Statement

This study used only previously published, de-identified public mouse transcriptomic data. No new human participants, human tissue, or identifiable patient data were analyzed.

## Author Contributions

Jiabiao Ji: Conceptualization, Methodology, Formal analysis, Visualization, Writing - original draft. Xiaoqing Yu: Data curation, Formal analysis, Validation. Chao Zhang: Software, Data curation. Long Liu: Investigation, Validation. Lin Yan: Visualization, Writing - review and editing. Hongjian Zhang: Resources, Writing - review and editing. Jianming Yang: Conceptualization, Supervision, Writing - review and editing.

## Funding

This research did not receive any specific grant from funding agencies in the public, commercial, or not-for-profit sectors.

## Conflicts of Interest

The authors declare no competing interests. This statement should be confirmed by all authors before submission.

## Acknowledgements

The authors thank the investigators who generated and deposited the public datasets analyzed in this study (GSE49543, GSE49522, GSE233798, GSE153882, GSE154833, and GSE274279) in the NCBI Gene Expression Omnibus, and acknowledge the Second Affiliated Hospital of Anhui Medical University for institutional support and computing resources. The authors also thank the developers and maintainers of the open-source software and resources used here, including R, limma, WGCNA, Seurat, CellChat, Enrichr, STRING, Scanpy, and the Mouse Genome Informatics gene nomenclature resources.

## References

1. D'Souza M, Zhu X, Frisina RD. Novel approach to select genes from RMA normalized microarray data using functional hearing tests in aging mice. J Neurosci Methods. 2008;171(2):279-287. doi:10.1016/j.jneumeth.2008.02.022. PMID: 18455804.
2. Sun G, et al. Single-cell transcriptomic atlas of mouse cochlear aging. Protein Cell. 2023;14(4):280-300. doi:10.1093/procel/pwac058. PMID: 36933008.
3. Zhao X, et al. Presbycusis: Pathology, signal pathways, and therapeutic strategy. Adv Sci (Weinh). 2025;12(31):e2410413. doi:10.1002/advs.202410413. PMID: 40349177.
4. Cheng Y, et al. Genetic analysis of potential biomarkers and therapeutic targets in age-related hearing loss. Hear Res. 2023;439:108894. doi:10.1016/j.heares.2023.108894. PMID: 37844444.
5. Jiang K, et al. Identification of mRNA expression profiles and their characterization in age-related hearing loss. Cell Mol Biol (Noisy-le-grand). 2024;70(4):235-241. doi:10.14715/cmb/2024.70.4.40. PMID: 38678595.
6. Zhang J, et al. Single cell RNA sequencing provides novel cellular transcriptional profiles and underlying pathogenesis of presbycusis. BMC Med Genomics. 2024;17(1):237. doi:10.1186/s12920-024-02001-7. PMID: 39350266.
7. Luo X, et al. Integrative multi-dataset analysis identifies immune-inflammatory hub genes as diagnostic biomarkers and therapeutic targets for age-related hearing loss. Hear Res. 2025;458:109448. doi:10.1016/j.heares.2025.109448. PMID: 41100952.
8. Wu MY, et al. Single-cell RNA sequencing analysis of mouse cochlea identifies a novel proinflammatory CD74(+)CD14(+) macrophage subset in mice with age-related and noise-induced hearing loss. Hear Res. 2025;466:109376. doi:10.1016/j.heares.2025.109376. PMID: 40743638.
9. Christensen N, D'Souza M, Zhu X, Frisina RD. Age-related hearing loss: aquaporin 4 gene expression changes in the mouse cochlea and auditory midbrain. Brain Res. 2009;1253:27-34. doi:10.1016/j.brainres.2008.11.070. PMID: 19070604.
10. Irizarry RA, Hobbs B, Collin F, et al. Exploration, normalization, and summaries of high density oligonucleotide array probe level data. Biostatistics. 2003;4(2):249-264. doi:10.1093/biostatistics/4.2.249. PMID: 12925520.
11. Ritchie ME, Phipson B, Wu D, et al. limma powers differential expression analyses for RNA-sequencing and microarray studies. Nucleic Acids Res. 2015;43(7):e47. doi:10.1093/nar/gkv007. PMID: 25605792.
12. Benjamini Y, Hochberg Y. Controlling the false discovery rate: a practical and powerful approach to multiple testing. J R Stat Soc Series B Stat Methodol. 1995;57(1):289-300. doi:10.1111/j.2517-6161.1995.tb02031.x.
13. Langfelder P, Horvath S. WGCNA: an R package for weighted correlation network analysis. BMC Bioinformatics. 2008;9:559. doi:10.1186/1471-2105-9-559. PMID: 19114008.
14. Langfelder P, Luo R, Oldham MC, Horvath S. Is my network module preserved and reproducible? PLoS Comput Biol. 2011;7(1):e1001057. doi:10.1371/journal.pcbi.1001057. PMID: 21283776.
15. Szklarczyk D, Kirsch R, Koutrouli M, et al. The STRING database in 2023: protein-protein association networks and functional enrichment analyses for any sequenced genome of interest. Nucleic Acids Res. 2023;51(D1):D638-D646. doi:10.1093/nar/gkac1000. PMID: 36370105.
16. Newman AM, Liu CL, Green MR, et al. Robust enumeration of cell subsets from tissue expression profiles. Nat Methods. 2015;12(5):453-457. doi:10.1038/nmeth.3337. PMID: 25822800.
17. Wolf FA, Angerer P, Theis FJ. SCANPY: large-scale single-cell gene expression data analysis. Genome Biol. 2018;19(1):15. doi:10.1186/s13059-017-1382-0. PMID: 29409532.
18. Wang Y, Zhao H, Zhao K, et al. Chronic inflammation and age-related hearing: based on Mendelian randomization. J Inflamm Res. 2024;17:8921-8934. doi:10.2147/JIR.S486301. PMID: 39569023.
19. Seicol BJ, Lin S, Xie R. Age-related hearing loss is accompanied by chronic inflammation in the cochlea and the cochlear nucleus. Front Aging Neurosci. 2022;14:846804. doi:10.3389/fnagi.2022.846804. PMID: 35418849.
20. Fuentes-Santamaria V, Alvarado JC, Mellado S, et al. Age-related inflammation and oxidative stress in the cochlea are exacerbated by long-term, short-duration noise stimulation. Front Aging Neurosci. 2022;14:853320. doi:10.3389/fnagi.2022.853320. PMID: 35450058.
21. Wu T, Zhou J, Qiu J, et al. Tumor necrosis factor-alpha mediated inflammation versus apoptosis in age-related hearing loss. Front Aging Neurosci. 2022;14:956503. doi:10.3389/fnagi.2022.956503. PMID: 36158549.
22. Chen Y, Huang H, Luo Y, et al. Senolytic treatment alleviates cochlear senescence and delays age-related hearing loss in C57BL/6J mice. Phytomedicine. 2025;142:156772. doi:10.1016/j.phymed.2025.156772. PMID: 40253743.
23. Chen J, Chen H, Wei Q, et al. APOE4 impairs macrophage lipophagy and promotes demyelination of spiral ganglion neurons in mouse cochleae. Cell Death Discov. 2025;11(1):190. doi:10.1038/s41420-025-02454-4. PMID: 40258814.
24. Li Y, Zeng T, Song W, et al. A single-cell transcriptome-wide association study reveals susceptibility genes for age-related hearing loss. Genomics Proteomics Bioinformatics. 2026. doi:10.1093/gpbjnl/qzaf137. PMID: 41495531.
25. Kim J, Lee B, Lee S, Kim JT, Kim BC, Cho HH. Identification and characterization of mRNA and lncRNA expression profiles in age-related hearing loss. Clin Exp Otorhinolaryngol. 2023;16(2):115-124. doi:10.21053/ceo.2022.01235. PMID: 36634670.
26. Feng M, Zhou X, Hu Y, et al. Comprehensive transcriptomic profiling of m6A modification in age-related hearing loss. Biomolecules. 2023;13(10):1537. doi:10.3390/biom13101537. PMID: 37892219.
27. Liu W, Ye B, Cai H, Zou Y, Zou Y. Bioinformatics and experiments reveal the hub genes of age-related hearing loss and the mechanism of PLK1 silencing in the protection of aging cochlear hair cells. Gene. 2025;963:149632. doi:10.1016/j.gene.2025.149632. PMID: 40516835.
28. Ma T, Zeng X, Liu M, et al. Analysis and identification of mitochondria-related genes associated with age-related hearing loss. BMC Genomics. 2025;26(1):218. doi:10.1186/s12864-025-11287-5. PMID: 40045222.
29. Gu X, Chen C, Chen Y, et al. Bioinformatics approach reveals the critical role of inflammation-related genes in age-related hearing loss. Sci Rep. 2025;15(1):2687. doi:10.1038/s41598-024-83428-x. PMID: 39837906.
30. Zhang C, Yang T, Luo X, et al. The chromatin accessibility and transcriptomic landscape of the aging mice cochlea and the identification of potential functional super-enhancers in age-related hearing loss. Clin Epigenetics. 2024;16(1):86. doi:10.1186/s13148-024-01702-1. PMID: 38965562.
31. Gu X, Jiang K, Chen R, et al. Identification of common stria vascularis cellular alteration in sensorineural hearing loss based on scRNA-seq. BMC Genomics. 2024;25(1):213. doi:10.1186/s12864-024-10122-7. PMID: 38413848.
32. Mafi AM, Tokar N, Russ MG, Barat O, Mellott JG. Age-related ultrastructural changes in the lateral cortex of the inferior colliculus. Neurobiol Aging. 2022;120:43-59. doi:10.1016/j.neurobiolaging.2022.08.007. PMID: 36116395.
33. Du EY, Ortega BK, Ninoyu Y, et al. Volumetric analysis of the aging auditory pathway using high resolution magnetic resonance histology. Front Aging Neurosci. 2022;14:1034073. doi:10.3389/fnagi.2022.1034073. PMID: 36437998.
34. Bartlett EL, Han EX, Parthasarathy A. Neurometric amplitude modulation detection in the inferior colliculus of young and aged rats. Hear Res. 2024;447:109028. doi:10.1016/j.heares.2024.109028. PMID: 38733711.
35. Butler T, Wang X, Chiang G, et al. Reduction in constitutively activated auditory brainstem microglia in aging and Alzheimer's disease. J Alzheimers Dis. 2024;99(1):307-319. doi:10.3233/JAD-231312. PMID: 38669537.
36. Liu J, Chen H, Lin X, et al. Age-related activation of cyclic GMP-AMP synthase-stimulator of interferon genes signaling in the auditory system is associated with presbycusis in C57BL/6J male mice. Neuroscience. 2022;481:73-84. doi:10.1016/j.neuroscience.2021.11.031. PMID: 34848262.
37. Ono M, Ito T. Hearing loss-related altered neuronal activity in the inferior colliculus. Hear Res. 2024;449:109033. doi:10.1016/j.heares.2024.109033. PMID: 38797036.
38. Eshel M, Milon B, Hertzano R, Elkon R. The cells of the sensory epithelium, and not the stria vascularis, are the main cochlear cells related to the genetic pathogenesis of age-related hearing loss. Am J Hum Genet. 2024;111(3):614-617. doi:10.1016/j.ajhg.2024.01.008. PMID: 38330941.
39. Liu H, Giffen KP, Chen L, et al. Molecular and cytological profiling of biological aging of mouse cochlear inner and outer hair cells. Cell Rep. 2022;39(2):110665. doi:10.1016/j.celrep.2022.110665. PMID: 35417713.
40. Zhang X, Cao R, Li C, et al. Caffeine ameliorates age-related hearing loss by downregulating the inflammatory pathway in mice. Otol Neurotol. 2024;45(3):227-237. doi:10.1097/MAO.0000000000004098. PMID: 38320571.
41. Jin S, Guerrero-Juarez CF, Zhang L, et al. Inference and analysis of cell-cell communication using CellChat. Nat Commun. 2021;12(1):1088. doi:10.1038/s41467-021-21246-9. PMID: 33597522.
42. Chen EY, Tan CM, Kou Y, et al. Enrichr: interactive and collaborative HTML5 gene list enrichment analysis tool. BMC Bioinformatics. 2013;14:128. doi:10.1186/1471-2105-14-128. PMID: 23586463.
43. Kuleshov MV, Jones MR, Rouillard AD, et al. Enrichr: a comprehensive gene set enrichment analysis web server 2016 update. Nucleic Acids Res. 2016;44(W1):W90-W97. doi:10.1093/nar/gkw377. PMID: 27141961.
44. Lachmann A, Xu H, Krishnan J, et al. ChEA: transcription factor regulation inferred from integrating genome-wide ChIP-X experiments. Bioinformatics. 2010;26(19):2438-2444. doi:10.1093/bioinformatics/btq466. PMID: 20709693.
45. Agarwal V, Bell GW, Nam JW, Bartel DP. Predicting effective microRNA target sites in mammalian mRNAs. Elife. 2015;4:e05005. doi:10.7554/eLife.05005. PMID: 26267216.
46. Yoo M, Shin J, Kim J, et al. DSigDB: drug signatures database for gene set analysis. Bioinformatics. 2015;31(18):3069-3071. doi:10.1093/bioinformatics/btv313. PMID: 25990557.

## Figures

**Figure 1. Study design and quality control.** Overview of GSE49543, GSE49522, and external datasets. The figure shows sample composition, ordinal severity coding, microarray quality-control summaries, and the computational workflow from limma and WGCNA to peripheral-central comparison, single-nucleus localization, and external evaluation.

![Figure 1](figures/figure1_study_design_qc.png)

**Figure 2. Severity-associated transcriptomic remodeling in the cochlea.** Pairwise limma results, ordinal severity effect estimates, and positive and negative severity signatures. The positive axis is enriched for myeloid, phagocytic, complement, and antigen-presentation genes, whereas the negative axis contains metabolic and structural genes.

![Figure 2](figures/figure2_cochlear_severity.png)

**Figure 3. Severity-associated WGCNA module and preservation.** Module-trait relationships, the M12 eigengene trajectory across groups, module membership, and preservation in the inferior colliculus. M12 is the strongest severity-associated module and shows reproducible co-expression structure in the central dataset.

![Figure 3](figures/figure3_wgcna_module12.png)

**Figure 4. Peripheral-central convergence and tissue specificity.** Gene-level comparison of cochlear and central severity associations, direction concordance, the 48 shared significant genes, and STRING enrichment for macrophage and microglial programs. The figure also shows the moderate genome-wide correlation between the two tissues.

![Figure 4](figures/figure4_peripheral_central.png)

**Figure 5. Myeloid program localization.** Marker-based immune enrichment in the cochlea and inferior colliculus, single-nucleus UMAP, major cell-type annotation, and candidate gene expression in macrophages/microglia. The marker-based analysis reports relative program enrichment rather than absolute cell proportions.

![Figure 5](figures/figure5_myeloid_localization.png)

**Figure 6. External context and computational model performance.** Candidate gene directionality in GSE233798, GSE153882, and GSE154833; nested cross-validation results; and external model transfer. The figure highlights partial reproducibility and the failure of the model to generalize to GSE233798.

![Figure 6](figures/figure6_external_validation_models.png)

**Figure 7. Severity-focused CellChat and enhancement analyses.** PTPRC-MRC1 communication in the severity-restricted CellChat run across ages, transcription-factor enrichment, miRNA target enrichment, second-pass GO and KEGG enrichment, and drug-signature enrichment. miRNA results are shown as a negative finding because no term reached FDR < 0.05.

![Figure 7](figures/figure7_enhancement_analysis.png)

**Figure 8. Sensitivity analyses of age, sex, and severity.** Numbers of significant genes under severity-only, sex-adjusted, and age-adjusted models; concordance between severity-only and sex-adjusted effects; concordance between severity-only and age-adjusted effects; and bootstrap confidence intervals for candidate genes.

![Figure 8](figures/figure8_sensitivity_analysis.png)

**Figure 9. Myeloid subclustering and module-score sensitivity.** UMAP of 147 myeloid nuclei, cluster proportions by age, module-score patterns by cluster and age, and equal-cell bootstrap estimates of 24-month-minus-3-month changes.

![Figure 9](figures/figure9_myeloid_subclustering.png)

**Figure 10. Full CellChatDB.mouse communication across presbycusis stages.** Significant cell-type-resolved communication after testing all 2,019 mouse ligand-receptor pairs. Panel A shows the number of significant edges per age with the corresponding numbers of significant ligand-receptor interactions and pathways. Panel B shows the total significant communication probability. Panel C shows pathway-level communication strength for the pathways that are strongest in at least one age group. Panel D compares pathway strength between 3 and 24 months against the identity line.

![Figure 10](figures/figure10_cellchat_full_database.png)

## Tables

**Table 1. Datasets and roles in the analysis.**

| Dataset | Tissue or cell context | Samples or cells | Platform or quantification | Role |
|---|---:|---:|---|---|
| GSE49543 | Cochlea | 41 samples | Affymetrix MOE430A microarray, RMA | Discovery of severity-associated peripheral programs |
| GSE49522 | Inferior colliculus | 39 samples | Affymetrix MOE430A microarray, RMA | Central comparison and module preservation |
| GSE233798 | Cochlea | 6 samples | RNA-seq FPKM | Independent cochlear aging context |
| GSE153882 | Inner and outer hair cells | 16 samples | RNA-seq RPKM | Hair-cell regional context |
| GSE154833 | Stria vascularis | 9 samples | RNA-seq RPKM | Regional tissue context |
| GSE274279 | Cochlea | 11,403 nuclei | 10x snRNA-seq | Myeloid localization |

**Table 2. Pairwise limma results in GSE49543.**

| Comparison | Up, FDR < 0.05 and log2FC >= 1 | Down, FDR < 0.05 and log2FC <= -1 | Up, FDR < 0.05 and log2FC >= 0.585 | Down, FDR < 0.05 and log2FC <= -0.585 |
|---|---:|---:|---:|---:|
| Middle-aged vs young control | 0 | 0 | 0 | 0 |
| Mild presbycusis vs middle-aged | 9 | 0 | 17 | 0 |
| Mild presbycusis vs young control | 11 | 1 | 24 | 1 |
| Severe presbycusis vs middle-aged | 12 | 3 | 29 | 9 |
| Severe presbycusis vs mild presbycusis | 0 | 0 | 0 | 0 |
| Severe presbycusis vs young control | 16 | 2 | 35 | 12 |

**Table 3. Severity-associated WGCNA module and preservation.**

| Module | Genes | Cochlear severity correlation | Preservation dataset | Zsummary | Permutation P | Example genes |
|---|---:|---:|---|---:|---:|---|
| M12 | 66 | 0.860 | GSE49522 inferior colliculus | 8.40 | 2.28e-21 | Mpeg1, Cd68, Clec7a, Fcgr3, H2-Aa, Ctss, Clec4d, Tyrobp, Csf1r, C1qc, Ms4a7 |

**Table 4. Shared significant genes between the cochlea and inferior colliculus.**

| Gene group | Representative genes |
|---|---|
| Macrophage and microglial identity | Mpeg1, Cd68, Csf1r, Tyrobp, Fcgr3, Clec7a, Ptprc, Emr1 |
| Complement and phagocytosis | C1qa, C1qb, C1qc, C3ar1, Lgals3, Cd14, Icam1 |
| Antigen presentation | Cd74, H2-Aa, H2-K1, H2-D1, Psmb8 |
| Fc receptor and myeloid signaling | Fcgr3, Fcer1g, Tyrobp, Nckap1l, Itgb2, Gpr65 |
| Other shared genes | Ms4a6d, Cd84, Il2rg, Ch25h, Csf2rb2, Cd48, Lyz1, Lyz2, Serpina3n, A2m, Lynx1, Uros |

**Table 5. Candidate evidence and external reproducibility.**

| Gene | Cochlear severity direction | Central severity direction | GSE233798 direction | Main interpretation |
|---|---|---|---|---|
| H2-Aa | Positive | Positive | Positive, FDR significant | Reproducible antigen-presentation axis |
| H2-Eb1 | Positive | Positive | Positive, FDR significant | Reproducible antigen-presentation axis |
| Cd74 | Positive | Positive | Positive, FDR significant | Myeloid antigen presentation; reproducible in one external model |
| Ctss | Positive | Positive | Positive, not FDR significant | Strong discovery and peripheral-central gene with model-dependent significance |
| Fcgr3 | Positive | Positive | Positive, not FDR significant | Strong shared myeloid gene |
| Cd68 | Positive | Positive | Positive, not FDR significant | Macrophage marker with stage-specific behavior |
| Tyrobp | Positive | Positive | Negative, FDR significant | Context-dependent myeloid signaling gene |
| Mpeg1 | Positive | Positive | Negative, not significant | Discovery-dominant gene requiring caution |
| Setd1a | Positive | Positive | Positive, FDR significant | Less myeloid-specific but externally reproducible |

**Table 6. Key strengths and limitations of the computational design.**

| Component | Strength | Limitation |
|---|---|---|
| Severity modeling | Uses ordinal progression instead of binary disease status | Age, cohort, and severity are not fully separable |
| Peripheral-central comparison | Directly compares cochlea and inferior colliculus | Cross-sectional data cannot establish causal direction |
| WGCNA | Identifies coordinated modules and preservation | Correlation-based and sensitive to input features |
| Immune enrichment | Uses mouse-specific marker programs | Estimates relative enrichment, not absolute cell proportions |
| Single-nucleus localization | Identifies likely cellular source | Limited macrophage/microglial nuclei and unbalanced ages |
| External validation | Tests directional reproducibility | Different models, platforms, tissues, and age contrasts |
| Machine learning | Prevents internal data leakage with nested cross-validation | External transfer was near random |
| CellChat | Tests the complete 2,019-interaction mouse ligand-receptor database | Permutation P values are not family-wise error controlled |
