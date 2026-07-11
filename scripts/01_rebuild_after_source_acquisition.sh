#!/usr/bin/env bash
set -euo pipefail
ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT_DIR"

python3 scripts/00_prepare_raw_inputs.py --strict

cat <<'MSG'

The public release rebuilds final tables and figures from data/derived/.

To regenerate data/derived/ from newly acquired public source files, follow
docs/DATA_ACQUISITION.md and keep each processing output under the documented
data/derived/ filenames. After updating data/derived/, run:

    bash scripts/reproduce_one_click.sh

This repository intentionally keeps live website acquisition separate from
the table/figure rebuild so that temporary HERB, ETCM2, LINCS, or GEO access
changes do not break the review-time reproduction path.
MSG
