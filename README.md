# UC B-cell curcumin public-data triangulation

This repository contains the public figure-and-table reproducibility package for a disease-first public-data study of ulcerative colitis/inflammatory bowel disease (UC/IBD), rectal B-cell disease-axis localization, and curcumin-linked target prioritization.

The package is designed to regenerate all main figures and manuscript tables from public-data-derived result tables, and to provide the supplementary TSV tables used in the submission package. It does not contain journal submission files, manuscript drafts, cover letters, credentials, private full text, raw downloaded public datasets, or scripts whose only purpose is to reacquire and reprocess the raw databases.

## What is included

- `scripts/`: figure/table reproduction entrypoints.
- `figures/scripts/build_main_figures.py`: rebuilds Figures 1-5.
- `figures/source_data/`: source-data tables used by figure panels.
- `figures/output/`: regenerated PNG, PDF, and SVG figure files. TIFF submission files are intentionally omitted because they are large and can be regenerated locally.
- `tables/source_data/` and `tables/manuscript/`: source and display tables for Tables 1-3.
- `tables/supplementary/`: supplementary TSV tables, including sensitivity analyses, immune-composition proxy analysis, target-overlap baseline comparison, and the GSE182270 B-lineage replication analysis.
- `results/`: lightweight derived result tables required to rebuild the figures and tables.
- `docs/DATA_MANIFEST.tsv`: public data-source manifest and acquisition records.
- `docs/FIGURE_PROVENANCE.tsv`: figure/table to script/source/output map.

## Quick reproduction

```bash
bash scripts/reproduce_one_click.sh
```

For a faster run when the required Python packages are already installed:

```bash
SKIP_INSTALL=1 bash scripts/reproduce_one_click.sh
```

Expected outputs:

- `figures/output/fig1_study_design_evidence_gating.{png,pdf,svg}`
- `figures/output/fig2_geo_disease_signature.{png,pdf,svg}`
- `figures/output/fig3_scrna_bcell_localization.{png,pdf,svg}`
- `figures/output/fig4_curcumin_target_bridge.{png,pdf,svg}`
- `figures/output/fig5_evidence_boundaries.{png,pdf,svg}`
- `tables/manuscript/table1_dataset_and_resource_inventory.tsv`
- `tables/manuscript/table2_prioritized_tcm_candidates.tsv`
- `tables/manuscript/table3_curcumin_target_evidence.tsv`
- `tables/supplementary/*.tsv`

## Data policy

The repository includes derived result tables sufficient to regenerate figures and tables. Raw GEO, HERB, ETCM2, LINCS, Open Targets, and PubMed/PMC inputs are not committed; acquisition sources are documented in `docs/DATA_MANIFEST.tsv`. Restricted full text is not redistributed. Full-text verification tables in this repository are sanitized to remove extracted private-text snippets.

## Compute expectations

Rebuilding figures and display tables from included result tables should take less than a minute on a standard laptop after Python dependencies are installed. Full raw-data reanalysis is larger and requires downloading the public resources listed in `docs/DATA_MANIFEST.tsv`.

## Version and archive

GitHub repository: https://github.com/Zhenghongwei11/uc-bcell-curcumin-triangulation

This repository is released through GitHub and archived by Zenodo. Cite the Zenodo version DOI associated with the exact GitHub release used for reproduction.

Zenodo version DOI for v1.0.1: https://doi.org/10.5281/zenodo.21148932

Zenodo concept DOI: https://doi.org/10.5281/zenodo.21121528
