# Compute Plan

## Lightweight reproduction

The default public reproduction command rebuilds tables and figures from included derived results. Expected resources:

- CPU: 2 cores or more
- Memory: <4 GB
- Runtime: usually <1 minute after dependency installation
- Disk: <1 GB for the public repository

## Full raw-data reanalysis

Full reanalysis requires reacquiring public raw or processed inputs from GEO, HERB 2.0, ETCM2, LINCS, Open Targets, and PubMed/PMC according to `docs/DATA_MANIFEST.tsv`. Large public raw files and restricted full text are intentionally excluded from this repository.

Recommended resources for full reanalysis:

- CPU: 4-8 cores
- Memory: 16-32 GB for single-cell matrix processing
- Disk: 20-50 GB depending on retained public raw inputs
- Network: stable access to NCBI GEO, HERB, ETCM2, Open Targets, and PubMed/PMC
