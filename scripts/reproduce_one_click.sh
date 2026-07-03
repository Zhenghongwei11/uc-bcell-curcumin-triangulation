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
PUBLIC_REPRO_BUILD=1 "$PYTHON_BIN" scripts/build_manuscript_tables_and_legends.py
"$PYTHON_BIN" figures/scripts/build_main_figures.py

"$PYTHON_BIN" - <<'PY'
from pathlib import Path
expected = [
    'figures/output/fig1_study_design_evidence_gating.png',
    'figures/output/fig2_geo_disease_signature.png',
    'figures/output/fig3_scrna_bcell_localization.png',
    'figures/output/fig4_curcumin_target_bridge.png',
    'figures/output/fig5_evidence_boundaries.png',
    'tables/manuscript/table1_dataset_and_resource_inventory.tsv',
    'tables/manuscript/table2_prioritized_tcm_candidates.tsv',
    'tables/manuscript/table3_curcumin_target_evidence.tsv',
    'tables/supplementary/Table_S1_full_curcumin_target_evidence.tsv',
    'tables/supplementary/gse182270_bcell_group_comparison.tsv',
    'tables/supplementary/standard_overlap_baseline_candidate_ranking.tsv',
    'tables/supplementary/scrna_signature_sensitivity_bcell_summary.tsv',
]
missing = [p for p in expected if not Path(p).exists()]
if missing:
    raise SystemExit('Missing expected outputs: ' + ', '.join(missing))
supplementary = sorted(Path('tables/supplementary').glob('*.tsv'))
if len(supplementary) < 30:
    raise SystemExit(f'Expected supplementary TSV tables, found only {len(supplementary)}')
print(f'Reproduction completed. Main figures, main tables, and {len(supplementary)} supplementary TSV tables are available.')
PY
