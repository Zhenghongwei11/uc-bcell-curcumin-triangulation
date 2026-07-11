# HERB and ETCM2 Access Notes

These notes record practical access details for the Chinese medicine knowledge bases used in this study.

## HERB 2.0

- Main website: `http://herb.ac.cn/`
- Direct host observed during this project: `http://47.92.70.12/`
- Browse page observed during this project: `http://herb.ac.cn/Browse/`
- The direct host was useful when the main domain did not render reliably in the controlled browser.
- Allowing third-party cookies helped the Browse interface load in the controlled browser environment.
- If a proxy is enabled system-wide, test whether `herb.ac.cn` and `47.92.70.12` should bypass it. Some local environments can open the direct host while the domain route fails.
- Avoid relying on a single interactive page state. Downloaded V2 text files are more stable than scraping visual pages.
- When recording compound-herb source chains, preserve the original HERB record ID and then standardize the displayed compound name separately. This prevents database spelling variants from entering final tables.

## ETCM2

- Public website: `https://www.tcmip.cn/ETCM2/front/`
- ETCM2 target names are not guaranteed to be HGNC gene symbols. Treat target-name-to-gene-symbol conversion as a mapping step and keep the mapping evidence.
- Ingredient search results and ingredient detail responses should be stored separately because search names, standardized names, and target names can differ.
- ETCM2 is best used as cross-resource context for this project, not as the sole source of target evidence.

## Recommended Record Keeping

For each manual download or export, record:

- Resource name and URL.
- Access date.
- Search term or record ID.
- File path in `data/raw/`.
- Whether the file was obtained automatically, manually downloaded, or exported from an interactive page.
- Any name standardization performed before analysis.
