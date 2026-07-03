# Reproducibility Notes

The public package separates two use cases.

## Fast figure/table rebuild

The committed result tables are derived from public databases and are sufficient to regenerate all main figures and Tables 1-3 quickly:

```bash
SKIP_INSTALL=1 bash scripts/reproduce_one_click.sh
```

The script sets `PUBLIC_REPRO_BUILD=1` so only public-facing table and figure artifacts are generated.

## Full public-data pipeline

To reacquire public inputs and regenerate the analysis outputs before rebuilding all figures and tables, run:

```bash
bash scripts/run_full_public_pipeline.sh
```

For repeat runs after public files already exist locally:

```bash
REUSE_EXISTING=1 SKIP_INSTALL=1 bash scripts/run_full_public_pipeline.sh
```

This pipeline downloads or queries public GEO, HERB 2.0, LINCS, L1000CDS2, Open Targets, PubMed/PMC, ETCM2, MyGene, and UniProt resources, then rebuilds `results/`, `figures/source_data/`, `figures/output/`, `tables/manuscript/`, and `tables/supplementary/`.

## Restricted full-text boundary

Restricted publisher full text is not redistributed. Rows derived from local restricted full-text inspection are retained only as sanitized metadata without copied text snippets. The full public pipeline preserves that sanitized evidence table by default so the submitted figures and tables can be regenerated without distributing restricted text. Set `REBUILD_FULLTEXT_AUDIT=1` only if legal local full-text inputs are available.
