#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT / "docs" / "RAW_INPUT_MANIFEST.tsv"
REPORT = ROOT / "docs" / "RAW_INPUT_STATUS.tsv"


def path_ready(path: Path) -> bool:
    if path.is_dir():
        return any(child.is_file() for child in path.rglob("*"))
    return path.is_file() and path.stat().st_size > 0


def main() -> int:
    parser = argparse.ArgumentParser(description="Check whether public source files have been placed under the expected data/raw layout.")
    parser.add_argument("--strict", action="store_true", help="Return a non-zero status when required source-rebuild files are missing.")
    args = parser.parse_args()

    rows = list(csv.DictReader(MANIFEST.open(encoding="utf-8"), delimiter="\t"))
    out_rows = []
    missing_required = []
    for row in rows:
        expected = ROOT / row["expected_local_path"]
        ready = path_ready(expected)
        status = "present" if ready else "missing"
        if row["required_for"] == "source_rebuild" and not ready:
            missing_required.append(row["resource_id"])
        out_rows.append(
            {
                "resource_id": row["resource_id"],
                "required_for": row["required_for"],
                "expected_local_path": row["expected_local_path"],
                "status": status,
                "acquisition_mode": row["acquisition_mode"],
                "public_url_or_accession": row["public_url_or_accession"],
            }
        )

    REPORT.parent.mkdir(parents=True, exist_ok=True)
    fields = ["resource_id", "required_for", "expected_local_path", "status", "acquisition_mode", "public_url_or_accession"]
    with REPORT.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(out_rows)

    present = sum(1 for row in out_rows if row["status"] == "present")
    print(f"Checked {len(out_rows)} public source resources.")
    print(f"Present: {present}; missing: {len(out_rows) - present}.")
    print(f"Report: {REPORT.relative_to(ROOT)}")
    if missing_required:
        print("Missing required source-rebuild resources:")
        for item in missing_required:
            print(f"  - {item}")
    return 1 if args.strict and missing_required else 0


if __name__ == "__main__":
    raise SystemExit(main())
