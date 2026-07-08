# UC B-cell curcumin public-data triangulation

This repository contains the public computational reproducibility package for a disease-first public-data study of ulcerative colitis/inflammatory bowel disease (UC/IBD), rectal B-cell disease-axis localization, and curcumin-linked target prioritization.

The package is designed to support two levels of reproduction: a fast rebuild of all manuscript figures and tables from tracked source tables, and a full public-data pipeline that reacquires public inputs and regenerates the analysis outputs used to build the figures and tables. It does not contain journal submission files, manuscript drafts, cover letters, credentials, private full text, or raw downloaded public datasets.

## What is included

- `scripts/`: public-data download, analysis, and figure/table reproduction entrypoints.
- `figures/scripts/build_main_figures.py`: rebuilds Figures 1-5.
- `figures/source_data/`: source-data tables used by figure panels.
- `figures/output/`: regenerated PNG, PDF, and SVG figure files. TIFF submission files are intentionally omitted because they are large and can be regenerated locally.
- `tables/source_data/` and `tables/manuscript/`: source and display tables for Tables 1-4.
- `tables/supplementary/`: supplementary TSV tables, including sensitivity analyses, immune-composition proxy analysis, target-overlap baseline comparison, therapeutic-direction support, candidate-specificity controls, and the GSE182270 B-lineage replication analysis.
- `results/`: derived result tables required to rebuild the figures and tables; these are also useful for fast verification without rerunning all public downloads.
- `docs/DATA_MANIFEST.tsv`: public data-source manifest and acquisition records.
- `docs/FIGURE_PROVENANCE.tsv`: figure/table to script/source/output map.

## Quick figure/table rebuild

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
- `figures/output/supplementary_figure_s1_robustness_replication.{png,pdf,svg}`
- `tables/manuscript/table1_dataset_and_resource_inventory.tsv`
- `tables/manuscript/table2_prioritized_tcm_candidates.tsv`
- `tables/manuscript/table3_curcumin_target_evidence.tsv`
- `tables/manuscript/table4_therapeutic_direction_consistency.tsv`
- `tables/supplementary/*.tsv`

## Full public-data reproduction

To reacquire public inputs and rerun the analysis pipeline before regenerating all figures and tables:

```bash
bash scripts/run_full_public_pipeline.sh
```

For repeated local runs after files have already been downloaded:

```bash
REUSE_EXISTING=1 SKIP_INSTALL=1 bash scripts/run_full_public_pipeline.sh
```

The full pipeline downloads or queries public GEO, HERB 2.0, LINCS, L1000CDS2, Open Targets, PubMed/PMC, ETCM2, MyGene, and UniProt resources, regenerates the main analysis outputs, rebuilds main figures and manuscript tables, and synchronizes supplementary TSV tables.

One legal boundary is explicit: restricted publisher full text is not redistributed. The repository therefore ships a sanitized full-text evidence table sufficient to reproduce the submitted figures and tables. Set `REBUILD_FULLTEXT_AUDIT=1` only if you have legal local access to the restricted full-text extraction expected by `scripts/audit_curcumin_target_fulltext_evidence.py`.

## Data policy

Raw GEO, HERB, ETCM2, LINCS, Open Targets, PubMed/PMC, MyGene, and UniProt inputs are not committed, but public acquisition is scripted or documented. Restricted full text is not redistributed. Full-text verification tables in this repository are sanitized to remove extracted private-text snippets.

## Compute expectations

Rebuilding figures and display tables from included result tables should take less than a minute on a standard laptop after Python dependencies are installed. Full public-data reanalysis downloads several hundred megabytes, mainly GEO single-cell supplementary files, and can take substantially longer depending on network access to GEO, HERB, ETCM2, L1000CDS2, Open Targets, MyGene, UniProt, and NCBI E-utilities.

## Version and archive

GitHub repository: https://github.com/Zhenghongwei11/uc-bcell-curcumin-triangulation

This repository is released through GitHub and archived by Zenodo. Cite the Zenodo version DOI associated with the exact GitHub release used for reproduction.

Zenodo concept DOI: https://doi.org/10.5281/zenodo.21121528
