# UC Rectal B-lineage and Chinese Medicine Compound Mapping

This repository contains the reproducibility package for a public-data computational pharmacology study of ulcerative colitis. The analysis defines a replicated mucosal disease program, localizes the program in rectal patient-level single-cell pseudobulk profiles, and maps Chinese medicine-related compound annotations in that disease-cell context.

## Repository Contents

- `scripts/`: table and figure rebuild scripts.
- `data/derived/`: lightweight analysis-ready tables used by the rebuild scripts.
- `tables/`: regenerated main and supplementary tables.
- `figures/source_data/`: source data written during figure generation.
- `plots/`: publication-ready PDF and PNG figure exports.
- `docs/`: data manifest, provenance files, statistical decision rules, and compute notes.

## Quick Reproduction

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
bash scripts/reproduce_one_click.sh
```

The one-click script rebuilds the Chinese medicine-focused main tables and figures from the derived tables in `data/derived/`.

## Expected Runtime

On a current laptop, the table rebuild usually finishes in under one minute and the figure rebuild in several minutes. No large raw public datasets are downloaded by the default reproduction path.

## Data Availability

All raw sources are public resources. The default reproduction path uses derived tables to keep the repository lightweight. Public accession IDs and resource URLs are listed in `docs/DATA_MANIFEST.tsv`.
