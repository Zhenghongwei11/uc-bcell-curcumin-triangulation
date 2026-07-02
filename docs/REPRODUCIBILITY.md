# Reproducibility Notes

The public package separates lightweight figure/table reproduction from full raw-data acquisition. The committed result tables are derived from public databases and are sufficient to regenerate all main figures and Tables 1-3.

To rebuild outputs, run:

```bash
SKIP_INSTALL=1 bash scripts/reproduce_one_click.sh
```

The script sets `PUBLIC_REPRO_BUILD=1` so only public-facing table and figure artifacts are generated.

Restricted full text is not redistributed. Rows derived from local restricted full-text inspection are retained only as sanitized metadata without copied text snippets.
