# Statistical Decision Rules

This public package reports disease-axis localization and candidate prioritization from public derived tables.

## Bulk transcriptomic signature

- GEO disease-control contrasts used processed expression matrices and platform annotation.
- Genes were summarized by log2 fold change and false-discovery-adjusted support.
- Figure 2 highlights genes passing FDR < 0.05 and absolute log2 fold change >= 1.
- The consensus signature was used as a disease anchor for cell-state localization and perturbational queries, not as a therapeutic mechanism by itself.

## Single-cell disease-axis localization

- Single-cell expression was aggregated to patient-level pseudobulk units before testing.
- Rectal cell-type contrasts use diseased-minus-healthy disease-axis deltas.
- Mann-Whitney P values are reported for pseudobulk group comparisons.
- Patient-label null testing used 10,000 permutations for disease-label specificity.
- Bootstrap confidence intervals are descriptive uncertainty summaries for key deltas.

## Candidate prioritization

- Bridge scores are rule-based prioritization scores integrating disease-context target evidence, disease-signature overlap, and cell-state localization.
- Scores are not effect sizes and do not measure clinical efficacy.
- Open Targets supports target-disease plausibility, not compound action.
- L1000CDS2 exact-identifier non-hits are interpreted as context-limited boundary results.
