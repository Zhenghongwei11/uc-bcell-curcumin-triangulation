# Statistical Confidence-Interval Supplement

Date: 2026-07-01

## Method

Patient-level pseudobulk values were resampled within disease groups using a fixed seed (20260701) and 20,000 bootstrap iterations. The estimand was the diseased-minus-healthy mean difference. Percentile 95% intervals are reported as uncertainty descriptors only; they do not replace the manuscript Mann-Whitney tests or patient-label null tests.

## Main Results

- B disease-axis delta: 0.538; bootstrap 95% CI 0.418 to 0.664 (n=7/4).
- M/DC disease-axis delta: 0.254; bootstrap 95% CI -0.036 to 0.565 (n=7/3).
- BCL6 B-cell log1p-CPM delta: 0.072; bootstrap 95% CI 0.042 to 0.096 (n=7/4).
- BLNK B-cell log1p-CPM delta: 0.282; bootstrap 95% CI 0.087 to 0.511 (n=7/4).
- SYK B-cell log1p-CPM delta: 0.213; bootstrap 95% CI 0.046 to 0.376 (n=7/4).

## Interpretation Boundaries

- Small healthy-group sample sizes make bootstrap intervals descriptive and sometimes wide; they should be used for transparency, not as a new validation layer.
- Gene-level intervals for BCL6, BLNK, and SYK support prioritization of a B-cell mechanism hypothesis only.
- Cytokine-context genes are included in the TSV for completeness but should remain boundary-framed in the manuscript.

## Outputs

- `tables/supplementary/statistical_ci_supplement.tsv`
- `figures/source_data/fig_statistical_ci_supplement.tsv`
