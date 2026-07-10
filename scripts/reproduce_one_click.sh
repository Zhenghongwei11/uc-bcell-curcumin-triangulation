#!/usr/bin/env bash
set -euo pipefail
ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT_DIR"
python3 scripts/build_chinese_medicine_rebuild_tables.py
python3 figures/scripts/build_chinese_medicine_figures.py
mkdir -p plots
cp figures/output/chinese_medicine/*.pdf plots/
cp figures/output/chinese_medicine/*.png plots/
python3 scripts/write_checksums.py
