# M7 TCM Candidate to scRNA Cell-Axis Bridge Run

## Scope

This run links HERB ingredient evidence to the GSE125527 rectal B-cell/M-DC disease-axis result. It uses local HERB evidence records and HERB detail API target relations when available; it does not infer compound-target edges from the standalone HERB target dictionary.

## Inputs

- `results/m2_herb2_uc_ibd_evidence/prioritized_ingredient_candidates.tsv`
- `results/m2_herb2_uc_ibd_evidence/uc_ibd_subject_evidence.tsv`
- `results/m3_lincs_coverage/lincs_candidate_coverage.tsv`
- `results/m6_scrna_cell_state/gse125527_disease_contrast_by_celltype.tsv`
- `results/m6_scrna_cell_state/gse125527_disease_contrast_patient_label_null.tsv` when available
- `results/m8_herb_ingredient_targets/herb_ingredient_target_summary.tsv` when available

## Key Outputs

- `results/m7_tcm_scrna_bridge/tcm_scrna_axis_bridge.tsv`
- `results/m7_tcm_scrna_bridge/tcm_scrna_axis_evidence_details.tsv`

## scRNA Anchor

- Rectal B-cell disease axis: delta=0.538393, Mann-Whitney P=0.00606061, patient-label null P=0.00229977.
- Rectal M/DC disease axis: delta=0.253899, Mann-Whitney P=0.383333, patient-label null P=0.20438; supportive only.

## Top Bridge Candidates

| Role | Ingredient | Score | B-cell evidence | Myeloid/DC evidence | Literature targets | UC/IBD up-target hits | LINCS covered |
|---|---|---:|---:|---:|---:|---:|---|
| primary_lead_for_b_cell_axis | Curcumin | 121.0 | 2 | 1 | 307 | 23 | yes |
| secondary_lead_for_myeloid_axis | Berberime | 74.35 | 0 | 2 | 77 | 13 | yes |
| secondary_lead_for_myeloid_axis | Quercetin | 59.0 | 0 | 1 | 189 | 23 | yes |
| secondary_lead_for_myeloid_axis | Wogonin | 48.65 | 0 | 1 | 73 | 16 | no |
| secondary_lead_for_myeloid_axis | Cis-resveratrol | 41.5 | 0 | 1 | 152 | 9 | no |
| secondary_lead_for_myeloid_axis | Taxifolin | 38.5 | 0 | 1 | 50 | 6 | yes |
| secondary_lead_for_myeloid_axis | Honokiol | 38.45 | 0 | 1 | 99 | 7 | yes |
| secondary_lead_for_myeloid_axis | Paeoniflorin | 35.35 | 0 | 1 | 77 | 7 | no |
| secondary_lead_for_myeloid_axis | Rhein | 32.75 | 0 | 1 | 35 | 6 | no |
| secondary_lead_for_myeloid_axis | Isoarnebin 4 | 31.7 | 0 | 1 | 64 | 5 | no |
| secondary_lead_for_myeloid_axis | Bilobalide | 28.35 | 0 | 1 | 27 | 4 | no |
| secondary_lead_for_myeloid_axis | Asperuloside | 23.55 | 0 | 1 | 11 | 0 | no |

## Interpretation

- Curcumin is the strongest current primary lead because it combines clinical/meta-analysis HERB evidence, direct B-cell-related colitis literature, HERB detail API target support, and exact-ID LINCS coverage, even though the L1000CDS2 reversal query did not rank it as a top reversal hit.
- Berberine remains a high-priority supporting lead because it has clinical/meta-analysis evidence, HERB detail API target support, and myeloid/macrophage-related mechanistic evidence, but it is better framed as secondary to the B-cell axis unless direct B-cell evidence is added.
- Very large database-integrated target sets are flagged as target-count-biased and should not be allowed to dominate the story without literature-supported targets and independent disease evidence.
- Compounds flagged as low specificity should not drive the TCM narrative even if HERB evidence counts are high.
- This bridge is literature/cell-axis triangulation, not wet-lab validation and not a compound-target proof.
