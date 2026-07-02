#!/usr/bin/env python3
"""Download HERB 2.0 V2 analysis tables and register them in data/manifest.tsv."""

from __future__ import annotations

import argparse
import csv
import hashlib
import shutil
import sys
import time
import urllib.parse
import urllib.request
from pathlib import Path
from typing import Dict, List


ACCESS_DATE = "2026-06-27"
SOURCE = "HERB 2.0"
VERSION = "V2"
DOWNLOAD_BASE = "http://47.92.70.12/download/file/"
REMOTE_ROOT = "/www/wwwroot/47.92.70.12/HERB_web/static/download_data/V2"

FILES = [
    "HERB_herb_info_v2.txt",
    "HERB_ingredient_info_v2.txt",
    "HERB_formula_info_v2.txt",
    "HERB_target_info_v2.txt",
    "HERB_disease_info_v2.txt",
    "HERB_meta_info_v2.txt",
    "HERB_clinical_trials_v2.txt",
    "HERB_reference_info_v2.txt",
    "HERB_experiment_info_v2.txt",
    "probe2gene.R",
]

MANIFEST_FIELDS = [
    "dataset_id",
    "source",
    "version_or_date",
    "file_path",
    "source_url_or_api",
    "access_date",
    "file_size_bytes",
    "md5",
    "status",
    "notes",
]


def direct_opener() -> urllib.request.OpenerDirector:
    return urllib.request.build_opener(urllib.request.ProxyHandler({}))


def source_url(filename: str) -> str:
    query = urllib.parse.urlencode({"file_path": f"{REMOTE_ROOT}/{filename}"})
    return f"{DOWNLOAD_BASE}?{query}"


def md5sum(path: Path) -> str:
    digest = hashlib.md5()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def text_audit(path: Path) -> Dict[str, str]:
    if path.suffix != ".txt":
        return {"line_count": "", "field_count": "", "header": ""}
    with path.open(encoding="utf-8", errors="replace", newline="") as handle:
        first = handle.readline().rstrip("\n\r")
        field_count = len(first.split("\t")) if first else 0
        line_count = 1 + sum(1 for _ in handle) if first else 0
    return {"line_count": str(line_count), "field_count": str(field_count), "header": " | ".join(first.split("\t"))}


def audit_existing_file(filename: str, out_dir: Path) -> Dict[str, str]:
    destination = out_dir / filename
    if not destination.exists() or destination.stat().st_size == 0:
        return {
            "filename": filename,
            "file_path": str(destination),
            "source_url_or_api": source_url(filename),
            "file_size_bytes": "",
            "md5": "",
            "status": "failed: missing existing file",
            "line_count": "",
            "field_count": "",
            "header": "",
        }
    audit = text_audit(destination)
    return {
        "filename": filename,
        "file_path": str(destination),
        "source_url_or_api": source_url(filename),
        "file_size_bytes": str(destination.stat().st_size),
        "md5": md5sum(destination),
        "status": "analysis_ready",
        **audit,
    }


def download_file(opener: urllib.request.OpenerDirector, filename: str, out_dir: Path, retries: int, timeout: int) -> Dict[str, str]:
    out_dir.mkdir(parents=True, exist_ok=True)
    destination = out_dir / filename
    tmp = out_dir / f".{filename}.tmp"
    url = source_url(filename)
    last_error = ""

    for attempt in range(retries + 1):
        try:
            request = urllib.request.Request(url, headers={"User-Agent": "zyy-herb2-download/1.0"})
            with opener.open(request, timeout=timeout) as response, tmp.open("wb") as handle:
                shutil.copyfileobj(response, handle, length=1024 * 1024)
            if tmp.stat().st_size == 0:
                raise RuntimeError("empty download")
            tmp.replace(destination)
            audit = text_audit(destination)
            return {
                "filename": filename,
                "file_path": str(destination),
                "source_url_or_api": url,
                "file_size_bytes": str(destination.stat().st_size),
                "md5": md5sum(destination),
                "status": "analysis_ready",
                **audit,
            }
        except Exception as exc:  # noqa: BLE001 - command-line acquisition should report exact failure.
            last_error = f"{type(exc).__name__}: {exc}"
            if tmp.exists():
                tmp.unlink()
            if attempt < retries:
                time.sleep(2 * (attempt + 1))

    return {
        "filename": filename,
        "file_path": str(destination),
        "source_url_or_api": url,
        "file_size_bytes": "",
        "md5": "",
        "status": f"failed: {last_error}",
        "line_count": "",
        "field_count": "",
        "header": "",
    }


def dataset_id(filename: str) -> str:
    stem = filename.replace(".txt", "").replace(".R", "")
    return f"herb2_full_{stem.lower()}"


def update_manifest(manifest_path: Path, audit_rows: List[Dict[str, str]]) -> None:
    existing: List[Dict[str, str]] = []
    if manifest_path.exists():
        with manifest_path.open(encoding="utf-8", newline="") as handle:
            existing = list(csv.DictReader(handle, delimiter="\t"))

    by_id = {row["dataset_id"]: row for row in existing if row.get("dataset_id")}
    for row in audit_rows:
        if row["status"] != "analysis_ready":
            continue
        by_id[dataset_id(row["filename"])] = {
            "dataset_id": dataset_id(row["filename"]),
            "source": SOURCE,
            "version_or_date": VERSION,
            "file_path": row["file_path"],
            "source_url_or_api": row["source_url_or_api"],
            "access_date": ACCESS_DATE,
            "file_size_bytes": row["file_size_bytes"],
            "md5": row["md5"],
            "status": "analysis_ready",
            "notes": f"{row['line_count']} lines; {row['field_count']} tab-delimited fields; HERB 2.0 full V2 acquisition",
        }

    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    with manifest_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=MANIFEST_FIELDS, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for key in sorted(by_id):
            writer.writerow(by_id[key])


def write_audit(path: Path, rows: List[Dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = ["filename", "file_path", "source_url_or_api", "file_size_bytes", "md5", "status", "line_count", "field_count", "header"]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out-dir", default="data/raw/herb2", help="Destination directory for HERB 2.0 full files.")
    parser.add_argument("--audit", default="results/herb2_full_download/file_audit.tsv", help="Output audit TSV.")
    parser.add_argument("--manifest", default="data/manifest.tsv", help="Manifest TSV to update.")
    parser.add_argument("--timeout", type=int, default=360)
    parser.add_argument("--retries", type=int, default=2)
    parser.add_argument("--reuse-existing", action="store_true", help="Audit existing files without downloading again.")
    parser.add_argument("--only", nargs="*", default=FILES, help="Optional subset of filenames to download.")
    args = parser.parse_args()

    opener = direct_opener()
    rows: List[Dict[str, str]] = []
    for filename in args.only:
        if args.reuse_existing:
            print(f"[herb2] auditing existing {filename}", file=sys.stderr)
            rows.append(audit_existing_file(filename, Path(args.out_dir)))
        else:
            print(f"[herb2] downloading {filename}", file=sys.stderr)
            rows.append(download_file(opener, filename, Path(args.out_dir), args.retries, args.timeout))

    write_audit(Path(args.audit), rows)
    update_manifest(Path(args.manifest), rows)

    failures = [row for row in rows if row["status"] != "analysis_ready"]
    if failures:
        for row in failures:
            print(f"[herb2] failed {row['filename']}: {row['status']}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
