# Compute Plan

The default reproduction path uses lightweight derived tables and does not download raw GEO, HERB, ETCM2, LINCS, PubMed, or Open Targets payloads.

Expected local requirements:

- Python 3.11 or later
- Less than 2 GB free disk space for the default rebuild
- Less than 8 GB RAM for the default table and figure rebuild

Large raw public datasets are not redistributed in this repository. Accession IDs and resource URLs are provided in `docs/DATA_MANIFEST.tsv`.
