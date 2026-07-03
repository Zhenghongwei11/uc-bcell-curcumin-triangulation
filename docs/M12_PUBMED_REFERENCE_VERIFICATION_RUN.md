# M12 PubMed Reference Verification Run

## Scope

This run independently verifies high-value Curcumin target references from M9 against official NCBI PubMed efetch XML records.

## Key Outputs

- `data/raw/pubmed/curcumin_pubmed_33597887_36196887_36353208.xml`
- `results/m12_pubmed_verification/curcumin_pubmed_article_verification.tsv`
- `results/m12_pubmed_verification/curcumin_target_edge_pubmed_verification.tsv`

## Verification Summary

- High-value PMIDs checked: 3 (33597887;36196887;36353208).
- PubMed records retrieved: 3.
- Article-level HERB title matches: 3/3.
- Edge-level title matches: 25/25.
- Edge-level target terms found in PubMed title/abstract: 23/25.

## Article-Level Verification

| PMID | Title match | Journal | Year | Disease/context terms |
|---|---|---|---:|---|
| 33597887 | yes | Frontiers in pharmacology | 2020 | curcumin;colitis;dss;memory_cell |
| 36196887 | yes | Phytotherapy research : PTR | 2023 | curcumin;colitis;ulcerative_colitis;dss;b_cell |
| 36353208 | yes | World journal of gastroenterology | 2022 | curcumin;colitis;dss;b_cell;memory_cell |

## Interpretation Boundary

- PubMed verification confirms bibliographic identity and title/abstract-level context, not full-text experimental validity.
- Target-term detection is conservative and based on title/abstract text plus common cytokine/gene aliases.
- Final manuscript tables should still manually inspect any high-impact claim before submission.
