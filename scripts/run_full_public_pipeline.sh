#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT_DIR"

if [[ "${SKIP_INSTALL:-0}" != "1" ]]; then
  python3 -m venv .venv
  # shellcheck disable=SC1091
  source .venv/bin/activate
  python -m pip install --upgrade pip
  python -m pip install -r requirements.txt
else
  PYTHON_BIN="${PYTHON_BIN:-python3}"
fi

PYTHON_BIN="${PYTHON_BIN:-python}"
REUSE_EXISTING="${REUSE_EXISTING:-0}"

echo "[1/22] Downloading public inputs required before analysis"
if [[ "$REUSE_EXISTING" == "1" ]]; then
  "$PYTHON_BIN" scripts/download_public_inputs.py --reuse-existing
else
  "$PYTHON_BIN" scripts/download_public_inputs.py
fi

echo "[2/22] Auditing GEO metadata"
"$PYTHON_BIN" scripts/audit_geo_uc_ibd_metadata.py

echo "[3/22] Downloading HERB 2.0 full public tables"
if [[ "$REUSE_EXISTING" == "1" ]]; then
  "$PYTHON_BIN" scripts/download_herb2_full.py --reuse-existing
else
  "$PYTHON_BIN" scripts/download_herb2_full.py
fi

echo "[4/22] Building HERB UC/IBD evidence layer"
"$PYTHON_BIN" scripts/build_herb2_uc_ibd_evidence.py

echo "[5/22] Prioritizing HERB ingredients"
"$PYTHON_BIN" scripts/prioritize_herb2_uc_ibd_ingredients.py

echo "[6/22] Building bulk GEO UC/IBD disease signatures"
"$PYTHON_BIN" scripts/build_geo_uc_ibd_signature.py

echo "[7/22] Exporting LINCS query signature"
"$PYTHON_BIN" scripts/export_lincs_query_signature.py

echo "[8/22] Checking LINCS perturbagen coverage"
"$PYTHON_BIN" scripts/check_lincs_candidate_coverage.py

echo "[9/22] Querying L1000CDS2 reversal endpoint"
"$PYTHON_BIN" scripts/query_l1000cds2_reversal.py

echo "[10/22] Scoring GSE125527 single-cell disease-axis localization"
"$PYTHON_BIN" scripts/analyze_gse125527_scrna_modules.py

echo "[11/22] Downloading HERB ingredient-detail target edges"
if [[ "$REUSE_EXISTING" == "1" ]]; then
  "$PYTHON_BIN" scripts/download_herb_ingredient_targets.py --reuse-existing
else
  "$PYTHON_BIN" scripts/download_herb_ingredient_targets.py
fi

echo "[12/22] Linking TCM candidates to the scRNA disease axis"
"$PYTHON_BIN" scripts/link_tcm_candidates_to_scrna_axis.py

echo "[13/22] Building disease-context target-cell bridge"
"$PYTHON_BIN" scripts/build_disease_context_target_cell_bridge.py

echo "[14/22] Querying Open Targets disease relevance"
if [[ "$REUSE_EXISTING" == "1" ]]; then
  "$PYTHON_BIN" scripts/query_open_targets_ibd_evidence.py
else
  "$PYTHON_BIN" scripts/query_open_targets_ibd_evidence.py --no-reuse-existing
fi

echo "[15/22] Verifying curcumin references through PubMed"
if [[ "$REUSE_EXISTING" == "1" ]]; then
  "$PYTHON_BIN" scripts/verify_curcumin_pubmed_references.py
else
  "$PYTHON_BIN" scripts/verify_curcumin_pubmed_references.py --no-reuse-existing
fi

echo "[16/22] Handling full-text evidence layer"
if [[ "${REBUILD_FULLTEXT_AUDIT:-0}" == "1" ]]; then
  "$PYTHON_BIN" scripts/audit_curcumin_target_fulltext_evidence.py
else
  test -s results/m16_curcumin_fulltext_evidence/curcumin_high_value_target_fulltext_audit.tsv
  echo "Preserved sanitized full-text evidence table. Set REBUILD_FULLTEXT_AUDIT=1 only when legal full-text inputs are available locally."
fi

echo "[17/22] Querying ETCM2 candidate cross-resource support"
if [[ "$REUSE_EXISTING" == "1" ]]; then
  "$PYTHON_BIN" scripts/query_etcm2_candidate_cross_validation.py
else
  "$PYTHON_BIN" scripts/query_etcm2_candidate_cross_validation.py --no-reuse-existing
fi

echo "[18/22] Mapping ETCM2 target names to accepted gene symbols"
"$PYTHON_BIN" scripts/map_etcm2_targets_to_genes.py

echo "[19/22] Running robustness analyses"
"$PYTHON_BIN" scripts/run_no_download_robustness.py

echo "[20/22] Running reviewer-facing lightweight extensions"
"$PYTHON_BIN" scripts/run_reviewer_lightweight_extensions.py

echo "[21/22] Running GSE182270 B-lineage replication analysis"
"$PYTHON_BIN" scripts/analyze_gse182270_bcell_replication.py

echo "[22/22] Building manuscript tables, supplementary tables, and figures"
PUBLIC_REPRO_BUILD=1 "$PYTHON_BIN" scripts/build_manuscript_tables_and_legends.py
PUBLIC_REPRO_BUILD=1 "$PYTHON_BIN" scripts/build_cjnm_strengthening_tables.py
"$PYTHON_BIN" figures/scripts/build_main_figures.py
"$PYTHON_BIN" scripts/build_statistical_ci_supplement.py
"$PYTHON_BIN" figures/scripts/build_supplementary_figures.py
"$PYTHON_BIN" scripts/sync_supplementary_tables.py

"$PYTHON_BIN" - <<'PY'
from pathlib import Path
expected = [
    "figures/output/fig1_study_design_evidence_gating.png",
    "figures/output/fig2_geo_disease_signature.png",
    "figures/output/fig3_scrna_bcell_localization.png",
    "figures/output/fig4_curcumin_target_bridge.png",
    "figures/output/fig5_evidence_boundaries.png",
    "figures/output/supplementary_figure_s1_robustness_replication.png",
    "tables/manuscript/table1_dataset_and_resource_inventory.tsv",
    "tables/manuscript/table2_prioritized_tcm_candidates.tsv",
    "tables/manuscript/table3_curcumin_target_evidence.tsv",
    "tables/manuscript/table4_therapeutic_direction_consistency.tsv",
    "tables/supplementary/cjnm_candidate_specificity_controls.tsv",
    "tables/supplementary/cjnm_therapeutic_direction_consistency_full.tsv",
    "tables/supplementary/gse182270_bcell_group_comparison.tsv",
    "tables/supplementary/standard_overlap_baseline_candidate_ranking.tsv",
    "tables/supplementary/scrna_signature_sensitivity_bcell_summary.tsv",
]
missing = [item for item in expected if not Path(item).exists()]
if missing:
    raise SystemExit("Missing expected full-pipeline outputs: " + ", ".join(missing))
print("Full public pipeline completed. Outputs are ready under figures/output, tables/manuscript, and tables/supplementary.")
PY
