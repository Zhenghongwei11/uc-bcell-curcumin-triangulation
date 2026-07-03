# M9 Disease-Context Target-to-Cell Bridge Run

## Scope

This run filters HERB `drug_paper_target` edges to UC/IBD/colitis-relevant references, classifies relationship direction, overlays targets with the bulk UC/IBD module, and quantifies target-gene expression in GSE125527 rectal scRNA cell types.

## Key Outputs

- `results/m9_target_cell_bridge/disease_context_herb_target_edges.tsv`
- `results/m9_target_cell_bridge/gse125527_target_gene_pseudobulk.tsv`
- `results/m9_target_cell_bridge/gse125527_target_gene_contrast.tsv`
- `results/m9_target_cell_bridge/candidate_target_cell_bridge_summary.tsv`
- `results/m9_target_cell_bridge/curcumin_bcell_target_evidence.tsv`

## Acquisition And Filtering Summary

- Disease-context HERB target edge rows: 188
- Disease-context unique target genes: 89
- Target genes present in GSE125527 rectal UMI tables: 47

## Top Candidate Target-to-Cell Bridges

| Ingredient | Score | Disease-context targets | Bulk up hits | B-cell expressed | B-cell disease-increased | M/DC expressed |
|---|---:|---:|---:|---:|---:|---:|
| Curcumin | 63.2 | 22 | 3 | 12 | 11 | 11 |
| Archin | 46 | 18 | 6 | 7 | 5 | 6 |
| Bilobalide | 34.8 | 12 | 2 | 7 | 6 | 7 |
| Berberime | 33.5 | 10 | 2 | 6 | 5 | 6 |
| Geniposide | 31.8 | 11 | 2 | 6 | 6 | 5 |
| Isoarnebin 4 | 29.8 | 9 | 2 | 6 | 5 | 5 |
| Cryptotanshinone | 29.5 | 11 | 2 | 6 | 3 | 6 |
| Polydatin | 27 | 5 | 2 | 4 | 4 | 4 |
| Berberine | 22 | 8 | 1 | 4 | 4 | 4 |
| Taxifolin | 22 | 7 | 1 | 4 | 4 | 4 |
| Trans-resveratrol | 22 | 7 | 1 | 4 | 4 | 4 |
| Honokiol | 21.5 | 6 | 2 | 2 | 2 | 2 |

## Curcumin Focus

- Curcumin disease-context targets: 22.
- Curcumin targets overlapping the bulk UC/IBD up module: 3 (CCL2;IL1B;IL33).
- Curcumin targets expressed in rectal B cells: 12 (BCL6;BLNK;IFNG;IL15;IL1B;IL7;JAK1;PIAS1;SH3KBP1;STAT5A;SYK;TNF).
- Curcumin targets increased in diseased rectal B-cell pseudobulk: 11 (BCL6;BLNK;IL15;IL1B;IL7;JAK1;PIAS1;SH3KBP1;STAT5A;SYK;TNF).

### Figure-Ready Curcumin Target Evidence

| Gene | Tier | Bulk role | B-cell delta | B-cell P | M/DC delta | PMIDs |
|---|---|---|---:|---:|---:|---|
| IL1B | A_module_aligned_b_cell_increased | bulk_up_top150 | 0.00113324 | 1 | 0.516526 | 36196887;36353208 |
| CCL2 | B_module_aligned | bulk_up_top150 |  |  |  | 36196887 |
| IL33 | B_module_aligned | bulk_up_top150 |  |  |  | 36196887 |
| BLNK | C_b_cell_mechanism_increased | not_bulk_top150 | 0.281591 | 0.0424242 | 0.0465185 | 36353208 |
| SYK | C_b_cell_mechanism_increased | not_bulk_top150 | 0.2135 | 0.0424242 | -0.156981 | 36353208 |
| SH3KBP1 | C_b_cell_mechanism_increased | not_bulk_top150 | 0.0850511 | 0.315152 | 0.107259 | 36353208 |
| BCL6 | C_b_cell_mechanism_increased | not_bulk_top150 | 0.0721226 | 0.0179007 | -0.312126 | 36353208 |
| TNF | D_b_cell_increased | not_bulk_top150 | 0.247649 | 0.109091 | 0.0843929 | 36196887;36353208 |
| JAK1 | D_b_cell_increased | not_bulk_top150 | 0.242715 | 0.163636 | 0.0716316 | 33597887 |
| PIAS1 | D_b_cell_increased | not_bulk_top150 | 0.0824401 | 0.230303 | -0.366305 | 33597887 |
| IL7 | D_b_cell_increased | not_bulk_top150 | 0.0748855 | 0.0727273 | 0.0289036 | 33597887 |
| STAT5A | D_b_cell_increased | not_bulk_top150 | 0.0486891 | 0.230303 | 0.0126604 | 33597887 |

## Interpretation Boundary

- The disease-context filter is text-based and should be manually audited for the final figure/table.
- Direction labels are inherited from HERB relationship text and supporting sentences; they are evidence annotations, not new experimental validation.
- A target is considered expressed if it is detected in patient-level pseudobulk for the cell class; disease increase is based on patient-level diseased-vs-healthy pseudobulk contrast.
- Celltype aggregation now deduplicates broad and detailed labels when no cluster label is available, preventing empty-cluster cells from being counted twice in the same celltype key.
- This module supports mechanistic prioritization. It does not prove therapeutic efficacy.
