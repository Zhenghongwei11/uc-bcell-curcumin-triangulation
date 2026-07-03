# M5 LINCS/L1000CDS2 Reversal Run

## Scope

This run used the GEO-derived UC/IBD consensus signature and queried the public L1000CDS2 `query2` endpoint in reverse mode. Candidate interpretation is restricted to the HERB-derived priority ingredients with exact LINCS perturbagen identifier coverage.

## Summary

| Metric | Value | Notes |
|---|---:|---|
| query_url | https://maayanlab.cloud/L1000CDS2/query2 |  |
| db_versions | cpcd-gse70138-v1.0;cpcd-gse70138-lm-v1.0 |  |
| up_genes_submitted | 100 | results/m4_geo_uc_ibd_signature/lincs_query_up_genes.txt |
| down_genes_submitted | 100 | results/m4_geo_uc_ibd_signature/lincs_query_down_genes.txt |
| run_utc | 2026-07-03T02:08:12+00:00 |  |
| run_epoch_seconds | 1783044492 |  |
| http_status.cpcd-gse70138-v1.0 | 200 |  |
| query_status.cpcd-gse70138-v1.0 | success |  |
| top_results_parsed.cpcd-gse70138-v1.0 | 50 | results/m5_lincs_reversal/l1000cds2_response.cpcd-gse70138-v1.0.json |
| http_status.cpcd-gse70138-lm-v1.0 | 200 |  |
| query_status.cpcd-gse70138-lm-v1.0 | success |  |
| top_results_parsed.cpcd-gse70138-lm-v1.0 | 50 | results/m5_lincs_reversal/l1000cds2_response.cpcd-gse70138-lm-v1.0.json |
| query_status | success | cpcd-gse70138-v1.0:success; cpcd-gse70138-lm-v1.0:success |
| top_results_parsed | 100 | results/m5_lincs_reversal/l1000cds2_top_results.tsv |
| covered_candidate_hit_rows | 0 | results/m5_lincs_reversal/covered_candidate_reversal.tsv |
| covered_candidate_tests | 15 | results/m5_lincs_reversal/candidate_reversal_status.tsv |
| covered_candidate_ingredient_hits | 0 |  |
| covered_candidate_pert_id_hits | 0 |  |

## Candidate Hits

No exact-ID HERB priority ingredient appeared in the returned L1000CDS2 top results. This is a negative public-API screen, not evidence that the compounds lack activity in the full LINCS matrix. See `results/m5_lincs_reversal/candidate_reversal_status.tsv` for the tested candidate set.

## Interpretation Boundary

- A hit means the public L1000CDS2 top-result list contains a LINCS perturbation that maps exactly to a prequalified HERB ingredient.
- A non-hit does not rule out reversal in the full LINCS Level 5 matrix because the public API returns a ranked subset and uses its own database/version filters.
- This module is suitable for prioritization and reviewer-facing triangulation, not for causal or efficacy claims.
