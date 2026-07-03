# M11 ETCM2 Candidate Cross-Validation Run

## Scope

This run attempts to query ETCM2 ingredient detail payloads for the top M9/M10 candidates, extract ETCM2 `ingredient_target` rows when available, extract `basic_network` target nodes for ingredient details, and compare gene-symbol-level targets with HERB disease-context targets, bulk-aligned targets, and Open Targets genetic-support genes where possible.

## Key Outputs

- `results/m11_etcm2_cross_validation/etcm2_candidate_query_audit.tsv`
- `results/m11_etcm2_cross_validation/etcm2_ingredient_target_edges.tsv`
- `results/m11_etcm2_cross_validation/candidate_etcm2_cross_validation_summary.tsv`

## Acquisition Summary

- Candidate ingredients attempted: 12.
- Browse pages requested: 3 with page size 500.
- Successful/cached candidate detail payloads: 9.
- Failed candidate detail attempts: 6.
- ETCM2 target edge rows extracted: 64.

## Top Cross-Validation Summary

| Ingredient | Query statuses | ETCM2 targets | HERB/M9 overlap | Bulk-aligned overlap | Open Targets genetic overlap |
|---|---|---:|---:|---:|---:|
| Taxifolin | cached | 39 | 0 | 0 | 0 |
| Berberine | cached | 9 | 0 | 0 | 0 |
| Curcumin | cached | 6 | 0 | 0 | 0 |
| Cryptotanshinone | cached | 6 | 0 | 0 | 0 |
| Bilobalide | cached | 2 | 0 | 0 | 0 |
| Geniposide | cached | 1 | 0 | 0 | 0 |
| Isoarnebin 4 | cached | 1 | 0 | 0 | 0 |
| Archin | cached;failed | 0 | 0 | 0 | 0 |
| Berberime | cached;failed | 0 | 0 | 0 | 0 |
| Polydatin | cached | 0 | 0 | 0 | 0 |
| Alpinetin | cached;failed | 0 | 0 | 0 | 0 |
| Asperuloside | cached | 0 | 0 | 0 | 0 |

## Interpretation Boundary

- ETCM2 detail `ingredient_target` rows and `basic_network` target nodes are treated as complementary target annotations, not proof of compound efficacy.
- If the ETCM2 backend times out, the audit table records failed attempts and the command can be rerun later with the same inputs.
- The script uses curl with `--noproxy www.tcmip.cn`, because this endpoint previously stalled through local proxy routing.
- Ingredient-detail target nodes are target names, not necessarily gene symbols; do not claim gene-level overlap unless a target-name mapping is added.
- The detail payload can expose a limited target preview; use extracted rows as cross-validation signals, not an exhaustive target universe unless full ETCM2 export access is obtained.
