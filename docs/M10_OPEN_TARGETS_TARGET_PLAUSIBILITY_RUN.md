# M10 Open Targets Target Plausibility Run

## Scope

This run adds an external disease-target plausibility layer from the Open Targets Platform GraphQL API for UC, IBD, and Crohn disease. It filters Open Targets associated targets to the M9 HERB disease-context target genes.

## Key Outputs

- `data/raw/open_targets/open_targets_ulcerative_colitis_MONDO_0005101.json`
- `data/raw/open_targets/open_targets_inflammatory_bowel_disease_MONDO_0005265.json`
- `data/raw/open_targets/open_targets_crohn_disease_MONDO_0005011.json`
- `results/m10_open_targets/open_targets_ibd_target_evidence.tsv`
- `results/m10_open_targets/candidate_open_targets_summary.tsv`

## Acquisition Summary

- ulcerative colitis (MONDO_0005101): 7161 associated targets cached; 7161 rows retrieved.
- inflammatory bowel disease (MONDO_0005265): 7629 associated targets cached; 7629 rows retrieved.
- Crohn disease (MONDO_0005011): 6289 associated targets cached; 6289 rows retrieved.
- M9 disease-context target genes queried against Open Targets: 89.
- Filtered Open Targets disease-target evidence rows retained: 242.

## Top Candidate External Disease-Target Support

| Ingredient | Disease-context targets | OT-supported | Genetic-supported | Clinical-supported | Literature-supported | Max OT score | Max genetic score |
|---|---:|---:|---:|---:|---:|---:|---:|
| Curcumin | 22 | 21 | 6 | 7 | 21 | 0.763493 | 0.88307 |
| Archin | 18 | 17 | 4 | 7 | 17 | 0.645758 | 0.808057 |
| Cryptotanshinone | 11 | 11 | 4 | 4 | 11 | 0.645758 | 0.760492 |
| Alpinetin | 6 | 6 | 4 | 1 | 6 | 0.645758 | 0.808057 |
| Bilobalide | 12 | 11 | 3 | 4 | 11 | 0.645758 | 0.808057 |
| Geniposide | 11 | 10 | 3 | 5 | 10 | 0.763493 | 0.88307 |
| Berberine | 8 | 7 | 3 | 3 | 7 | 0.645758 | 0.744156 |
| Asperuloside | 7 | 7 | 3 | 1 | 7 | 0.545518 | 0.846731 |
| Polydatin | 5 | 5 | 3 | 1 | 5 | 0.567615 | 0.808057 |
| Isoarnebin 4 | 9 | 9 | 2 | 5 | 9 | 0.763493 | 0.88307 |
| Taxifolin | 7 | 7 | 2 | 5 | 7 | 0.645758 | 0.808057 |
| Berberime | 10 | 9 | 1 | 4 | 9 | 0.645758 | 0.760492 |

## Curcumin Focus

- Open Targets-supported curcumin targets: 21 (BCL6;BLNK;CCL2;IFNG;IL10;IL12A;IL13;IL15;IL1B;IL22;IL33;IL4;IL6;IL7;INS;JAK1;PIAS1;SH3KBP1;STAT5A;SYK;TNF).
- Genetic-supported curcumin targets: 6 (CCL2;IFNG;IL10;IL33;INS;SYK).
- Clinical-supported curcumin targets: 7 (IFNG;IL12A;IL13;IL1B;IL6;JAK1;TNF).
- Max Open Targets score among curcumin targets: 0.763493; max genetic association score: 0.88307.

## Interpretation Boundary

- Open Targets supports disease-target plausibility, not compound efficacy.
- Genetic or clinical target support strengthens target relevance, but it does not validate that a TCM ingredient modulates that target in patients.
- Scores are used as orthogonal prioritization evidence and should be reported with source/date rather than treated as experimental results.
