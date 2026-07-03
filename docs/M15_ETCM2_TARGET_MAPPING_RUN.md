# M15 ETCM2 Target-Name To Gene-Symbol Mapping Run

## Scope

This run maps ETCM2 ingredient-detail target names to human gene symbols using UniProt reviewed protein-name search and MyGene human gene queries. Original ETCM2 target names are preserved; only exact/high-confidence mappings are used for gene-level overlap.

## Key Outputs

- `results/m15_etcm2_target_mapping/etcm2_target_name_gene_mapping.tsv`
- `results/m15_etcm2_target_mapping/etcm2_target_edges_mapped.tsv`
- `results/m15_etcm2_target_mapping/candidate_etcm2_mapped_overlap_summary.tsv`

## Mapping Summary

- Unique ETCM2 target names: 59.
- Exact mappings: 53.
- High-confidence mappings: 4.
- Medium-confidence mappings: 0.
- Low-confidence mappings: 0.
- Ambiguous mappings: 2.
- Unmapped target names: 0.

## Candidate Gene-Level Overlap After Mapping

| Ingredient | Accepted mapped genes | HERB/M9 overlap | Bulk-aligned overlap | Open Targets genetic overlap |
|---|---:|---:|---:|---:|
| Taxifolin | 36 | 1 | 0 | 1 |
| Berberine | 8 | 0 | 0 | 0 |
| Cryptotanshinone | 6 | 0 | 0 | 0 |
| Curcumin | 5 | 0 | 0 | 0 |
| Bilobalide | 2 | 0 | 0 | 0 |
| Geniposide | 1 | 0 | 0 | 0 |
| Isoarnebin 4 | 1 | 0 | 0 | 0 |

## Interpretation Boundary

- Gene-level overlap is computed only from `exact` and `high` mappings.
- `medium`, `low`, `ambiguous`, and `unmapped` target-name mappings should remain in supplementary audit tables, not primary mechanistic claims.
- Mapping target names to genes can introduce ambiguity for family-level names such as `Cytochrome P450 1A` or `Cholinesterase`; these are explicitly downgraded.
