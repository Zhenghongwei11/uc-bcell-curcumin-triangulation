# M6 GSE125527 scRNA Cell-State Localization Run

## Scope

This run scores the GEO bulk-derived UC/IBD consensus up/down modules in processed GSE125527 single-cell UMI data. It uses author-provided sample/tissue/disease/cell-type labels and patient-level pseudobulk contrasts.

## Key Outputs

- `results/m6_scrna_cell_state/gse125527_module_gene_coverage.tsv`
- `results/m6_scrna_cell_state/gse125527_celltype_module_summary.tsv`
- `results/m6_scrna_cell_state/gse125527_pseudobulk_module_summary.tsv`
- `results/m6_scrna_cell_state/gse125527_disease_contrast_by_celltype.tsv`
- `results/m6_scrna_cell_state/gse125527_disease_contrast_by_celltype_detail.tsv`
- `results/m6_scrna_cell_state/gse125527_disease_contrast_patient_label_null.tsv`

## Gene Coverage

| Module | Query genes | Present in scRNA | Coverage |
|---|---:|---:|---:|
| bulk_consensus_up | 150 | 52 | 0.346667 |
| bulk_consensus_down | 150 | 20 | 0.133333 |

## Rectal Tissue Disease Axis Contrast

| Cell type | Diseased samples | Healthy samples | Delta axis | MW P value | Label-null empirical P | Diseased cells | Healthy cells |
|---|---:|---:|---:|---:|---:|---:|---:|
| B | 7 | 4 | 0.538393 | 0.00606061 | 0.00229977 | 3350 | 1425 |
| M/DC | 7 | 3 | 0.253899 | 0.383333 | 0.20438 | 201 | 54 |
| NK | 7 | 4 | 0.0307754 | 0.527273 | 0.931807 | 261 | 52 |
| T | 7 | 4 | 0.00586036 | 1 | 0.956204 | 4237 | 1853 |
| unknown | 7 | 8 | -0.0423548 | 0.53582 | 0.539346 | 10593 | 10128 |

## Actionable Interpretation

- Main publishable signal: rectal B cells show the strongest positive disease-axis shift (delta=0.538393, Mann-Whitney P=0.00606061, patient-label null P=0.00229977) using patient/sample-level pseudobulk scores.
- Supportive direction: rectal M/DC cells are also positive (delta=0.253899, Mann-Whitney P=0.383333), but the healthy comparator has fewer samples/cells, so this should be treated as secondary evidence.
- Recommended manuscript framing: public UC/IBD mucosal transcriptional injury localizes most robustly to a rectal B-cell inflammatory axis; TCM candidates should be connected to this axis by target/pathway and genetic evidence before making mechanism claims.

## Interpretation Boundary

- The primary contrast is patient/sample-level pseudobulk, not cell-level p-value inflation.
- The patient-label null analysis permutes disease labels within each tissue/celltype pseudobulk contrast; it is a sensitivity check for the cell-state disease-axis association, not experimental validation.
- A positive delta means the bulk UC/IBD up module is relatively higher than the down module in diseased samples for that cell class.
- This analysis localizes the public bulk disease signature; it does not validate any TCM compound mechanism by itself.
