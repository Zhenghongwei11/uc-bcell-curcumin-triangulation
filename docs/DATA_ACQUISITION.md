# Public Source Data Acquisition

This repository provides two reproducibility layers.

1. The default layer rebuilds all reported tables and figures from analysis-ready public-data-derived tables in `data/derived/`.
2. The optional source-data layer documents where the public source files were obtained and how to place them locally before rebuilding derived tables.

The source-data layer is documented rather than forced into a single download command because public databases differ in access stability, file size, and interface behavior. GEO FTP files are stable and can usually be downloaded automatically. HERB and ETCM2 may require direct browser access or manual export if the website changes or local network rules interfere.

## Expected Directory Layout

Place acquired source files under:

```text
data/raw/
  geo/
    GSE125527/
    GSE182270/
    GSE59071/
    GSE75214/
    GSE87466/
    platforms/
  herb2/
  etcm2/
  lincs/
  open_targets/
  pubmed/
```

Run the checker after placing files:

```bash
python3 scripts/00_prepare_raw_inputs.py
```

Use `--strict` to return a non-zero exit status when required source files are missing.

## GEO

Bulk microarray inputs:

- GSE75214 series matrix: `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE75nnn/GSE75214/matrix/GSE75214_series_matrix.txt.gz`
- GSE59071 series matrix: `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE59nnn/GSE59071/matrix/GSE59071_series_matrix.txt.gz`
- GSE87466 series matrix: `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE87nnn/GSE87466/matrix/GSE87466_series_matrix.txt.gz`
- GPL6244 annotation: `https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL6nnn/GPL6244/annot/GPL6244.annot.gz`
- GPL13158 annotation: `https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL13nnn/GPL13158/annot/GPL13158.annot.gz`

Single-cell inputs:

- GSE125527 GEO record: `https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE125527`
- GSE182270 GEO record: `https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE182270`

For GSE125527, download the supplementary RAW archive and associated metadata files listed in `docs/RAW_INPUT_MANIFEST.tsv`. The analysis uses patient-level pseudobulk summaries, so metadata files assigning patient, tissue, disease status, and cell type are required.

## HERB 2.0

HERB 2.0 source files used by this project are the V2 tab-delimited tables:

- `HERB_herb_info_v2.txt`
- `HERB_ingredient_info_v2.txt`
- `HERB_formula_info_v2.txt`
- `HERB_target_info_v2.txt`
- `HERB_disease_info_v2.txt`
- `HERB_meta_info_v2.txt`
- `HERB_clinical_trials_v2.txt`
- `HERB_reference_info_v2.txt`
- `HERB_experiment_info_v2.txt`
- `probe2gene.R`

The public website is `http://herb.ac.cn/`; during this project the direct host `http://47.92.70.12/` was also reachable. Store the V2 files in `data/raw/herb2/`.

Ingredient-detail pages were used only to verify compound-to-herb source chains for the prioritized records. If repeating that step, store detail JSON files under `data/raw/herb2/ingredient_details/`.

## ETCM2

ETCM2 was used as a cross-resource target-name context rather than as the primary source of compound-target evidence. The public site is:

`https://www.tcmip.cn/ETCM2/front/`

Place downloaded or exported ingredient search and detail JSON files under:

```text
data/raw/etcm2/ingredient_search/
data/raw/etcm2/ingredient_details/
```

ETCM2 target names may be protein names rather than HGNC symbols. Gene-symbol mapping should therefore be treated as an evidence-normalization step and recorded in derived mapping tables.

## Open Targets

Open Targets was used for target-disease plausibility context. Store JSON exports under `data/raw/open_targets/`:

- `open_targets_ulcerative_colitis_MONDO_0005101.json`
- `open_targets_inflammatory_bowel_disease_MONDO_0005265.json`
- `open_targets_crohn_disease_MONDO_0005011.json`

## LINCS and L1000CDS2

LINCS perturbation metadata used for coverage checks:

- `GSE70138_Broad_LINCS_pert_info.txt.gz`
- `GSE92742_Broad_LINCS_pert_info.txt.gz`

L1000CDS2 was queried through the public service at `https://maayanlab.cloud/L1000CDS2/`. Store downloaded responses under `data/raw/lincs/` or keep the processed response summaries in `data/derived/`.

## PubMed and PMC

PubMed and PMC XML files were used for bibliographic verification of selected curcumin-associated records. Store downloaded XML files under `data/raw/pubmed/`.

## From Source Files To Figures

After source files are obtained and inspected, regenerate or replace the derived tables under `data/derived/`, then run:

```bash
bash scripts/reproduce_one_click.sh
```

The one-click command is intentionally independent of live web access. This prevents review-time failures caused by temporary database downtime or website-interface changes.
