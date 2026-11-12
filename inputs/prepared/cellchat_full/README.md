# S12 CellChat Analysis

## Library subset

The CellChat analysis used the three age-specific libraries GSM8446833 (3M),
GSM8446834 (12M), and GSM8446835 (24M). It did not include GSM8446836
(24M_Cochlea_BU). S12 therefore does not summarize all four GSE274279 libraries
or all 14,152 nuclei.

The input mice were C57BL/6J. S12 is cross-strain descriptive context for the
CBA/CaJ bulk reanalysis.

## Statistical and biological boundary

CellChat probabilities are model-inferred from normalized expression and
CellChatDB.mouse. They are not measurements of ligand secretion, receptor
activation, signal flux, protein-complex formation, or functional dependence.
Permutation P values were not additionally family-wise error controlled.

Age is inseparable from library and 10x batch. Cell-type composition and
communication summaries are descriptive. The animal-level relationship between
the two 24-month libraries in S11 is not documented in GEO, and CellChat does
not include the second 24-month library.

## Files

- `CellChat_full_age_summary.csv`: edge, interaction, and pathway counts by age
  in the three-library subset.
- `CellChat_full_edges_all_ages.csv`: cell-pair edges with numeric `pval` and
  display `pval_display`.
- `CellChat_full_interaction_presence.csv`: persistence of interactions across
  ages.
- `CellChat_full_focus_interactions.csv`: PTPRC-MRC1, APP-CD74, and other focus
  interactions.
- `CellChat_full_pathway_edges_all_ages.csv`: pathway-level edge probabilities.
- `CellChat_full_pathway_strength.csv`: pathway summaries.
- `CellChat_full_celltype_strength.csv`: incoming and outgoing strength by cell
  type.
- `CellChat_full_interpretation.md`: database-level descriptive summary.
