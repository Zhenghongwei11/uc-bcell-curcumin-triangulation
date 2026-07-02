#!/usr/bin/env python3
"""Build UC/IBD bulk GEO disease signatures from selected public datasets."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import math
import subprocess
import sys
from collections import defaultdict
from pathlib import Path
from typing import Dict, Iterable, List, Tuple

import numpy as np
from scipy import stats


DATASETS = {
    "GSE75214": {
        "platform": "GPL6244",
        "role": "discovery_or_replication",
        "contrast": "active_uc_colon_vs_normal_colon",
        "evidence_group": "GPL6244_active_uc_colon",
    },
    "GSE59071": {
        "platform": "GPL6244",
        "role": "non_independent_sensitivity",
        "contrast": "active_uc_colon_vs_normal_colon",
        "evidence_group": "GPL6244_active_uc_colon",
    },
    "GSE87466": {
        "platform": "GPL13158",
        "role": "uc_active_validation",
        "contrast": "uc_colon_vs_normal_colon",
        "evidence_group": "GPL13158_uc_colon",
    },
}

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


def gse_prefix(accession: str) -> str:
    return f"GSE{accession[3:-3]}nnn"


def gpl_prefix(platform: str) -> str:
    return f"GPL{platform[3:-3]}nnn"


def series_matrix_url(accession: str) -> str:
    return f"https://ftp.ncbi.nlm.nih.gov/geo/series/{gse_prefix(accession)}/{accession}/matrix/{accession}_series_matrix.txt.gz"


def platform_annot_url(platform: str) -> str:
    return f"https://ftp.ncbi.nlm.nih.gov/geo/platforms/{gpl_prefix(platform)}/{platform}/annot/{platform}.annot.gz"


def run_curl(url: str, destination: Path, timeout: int) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    tmp = destination.with_suffix(destination.suffix + ".tmp")
    cmd = [
        "curl",
        "-L",
        "-sS",
        "--fail",
        "--retry",
        "3",
        "--max-time",
        str(timeout),
        url,
        "-o",
        str(tmp),
    ]
    subprocess.run(cmd, check=True)
    subprocess.run(["gzip", "-t", str(tmp)], check=True)
    tmp.replace(destination)


def md5sum(path: Path) -> str:
    digest = hashlib.md5()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def update_manifest(manifest_path: Path, rows: List[Dict[str, str]]) -> None:
    existing: Dict[str, Dict[str, str]] = {}
    if manifest_path.exists():
        with manifest_path.open(encoding="utf-8", newline="") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                existing[row["dataset_id"]] = row
    for row in rows:
        existing[row["dataset_id"]] = row
    with manifest_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=MANIFEST_FIELDS, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for key in sorted(existing):
            writer.writerow(existing[key])


def read_metadata(path: Path) -> List[Dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def has_term(row: Dict[str, str], field: str, term: str) -> bool:
    return term in {part.strip() for part in row.get(field, "").split(";") if part.strip()}


def in_colon(row: Dict[str, str]) -> bool:
    return has_term(row, "tissue_terms", "colon") or has_term(row, "tissue_terms", "colonic")


def is_active(row: Dict[str, str]) -> bool:
    return has_term(row, "inflammation_terms", "active") or has_term(row, "inflammation_terms", "inflamed")


def is_uc(row: Dict[str, str]) -> bool:
    return has_term(row, "disease_terms", "uc") or has_term(row, "disease_terms", "ulcerative colitis")


def assign_group(accession: str, row: Dict[str, str]) -> str:
    if row.get("accession") != accession:
        return "unused"
    if accession in {"GSE75214", "GSE59071"}:
        if is_uc(row) and in_colon(row) and is_active(row) and row.get("is_control_like") != "yes":
            return "case"
        if row.get("is_control_like") == "yes" and in_colon(row):
            return "control"
        return "excluded"
    if accession == "GSE87466":
        if is_uc(row) and in_colon(row) and row.get("is_control_like") != "yes":
            return "case"
        if row.get("is_control_like") == "yes" and in_colon(row):
            return "control"
        return "excluded"
    return "excluded"


def load_platform_annotation(path: Path) -> Dict[str, Dict[str, str]]:
    annotation: Dict[str, Dict[str, str]] = {}
    with gzip.open(path, "rt", encoding="utf-8", errors="replace", newline="") as handle:
        in_table = False
        reader = None
        for line in handle:
            if line.startswith("!platform_table_begin"):
                in_table = True
                header = next(handle).rstrip("\n")
                reader = csv.DictReader(handle, fieldnames=header.split("\t"), delimiter="\t")
                continue
            if not in_table or reader is None:
                continue
            if line.startswith("!platform_table_end"):
                break
            row = next(csv.DictReader([line], fieldnames=reader.fieldnames, delimiter="\t"))
            probe_id = row.get("ID", "")
            if not probe_id:
                continue
            gene_symbol = first_token(row.get("Gene symbol", ""))
            gene_title = first_token(row.get("Gene title", ""))
            gene_id = first_token(row.get("Gene ID", ""))
            annotation[probe_id] = {
                "gene_symbol": gene_symbol,
                "gene_title": gene_title,
                "gene_id": gene_id,
            }
    return annotation


def first_token(value: str) -> str:
    value = value or ""
    if value in {"NA", "---"}:
        return ""
    return value.split("///")[0].strip()


def iter_matrix_rows(path: Path):
    with gzip.open(path, "rt", encoding="utf-8", errors="replace", newline="") as handle:
        in_table = False
        header: List[str] | None = None
        for line in handle:
            if line.startswith("!series_matrix_table_begin"):
                in_table = True
                continue
            if not in_table:
                continue
            if line.startswith("!series_matrix_table_end"):
                break
            parts = next(csv.reader([line.rstrip("\n")], delimiter="\t"))
            if header is None:
                header = [part.strip('"') for part in parts]
                yield "HEADER", header
                continue
            yield parts[0].strip('"'), [part.strip('"') for part in parts[1:]]


def bh_fdr(p_values: List[float]) -> List[float]:
    n = len(p_values)
    order = sorted(range(n), key=lambda idx: p_values[idx])
    fdr = [1.0] * n
    prev = 1.0
    for rank_from_end, idx in enumerate(reversed(order), start=1):
        rank = n - rank_from_end + 1
        value = min(prev, p_values[idx] * n / rank)
        fdr[idx] = min(value, 1.0)
        prev = value
    return fdr


def safe_float(value: str) -> float:
    try:
        return float(value)
    except ValueError:
        return math.nan


def build_signature(
    accession: str,
    matrix_path: Path,
    metadata: List[Dict[str, str]],
    annotation: Dict[str, Dict[str, str]],
    out_dir: Path,
) -> Tuple[List[Dict[str, str]], Dict[str, str]]:
    case_gsms = {row["gsm"] for row in metadata if assign_group(accession, row) == "case"}
    control_gsms = {row["gsm"] for row in metadata if assign_group(accession, row) == "control"}
    excluded = sum(1 for row in metadata if row.get("accession") == accession and assign_group(accession, row) == "excluded")

    header = None
    rows: List[Dict[str, str]] = []
    case_idx: List[int] = []
    control_idx: List[int] = []

    for probe_id, values in iter_matrix_rows(matrix_path):
        if probe_id == "HEADER":
            header = values
            sample_ids = header[1:]
            case_idx = [idx for idx, gsm in enumerate(sample_ids) if gsm in case_gsms]
            control_idx = [idx for idx, gsm in enumerate(sample_ids) if gsm in control_gsms]
            continue
        if header is None:
            raise RuntimeError(f"{accession}: matrix header not found")
        if len(case_idx) < 3 or len(control_idx) < 3:
            raise RuntimeError(f"{accession}: insufficient groups: case={len(case_idx)}, control={len(control_idx)}")
        numeric = np.array([safe_float(value) for value in values], dtype=float)
        case_values = numeric[case_idx]
        control_values = numeric[control_idx]
        case_values = case_values[~np.isnan(case_values)]
        control_values = control_values[~np.isnan(control_values)]
        if len(case_values) < 3 or len(control_values) < 3:
            continue
        mean_case = float(np.mean(case_values))
        mean_control = float(np.mean(control_values))
        stat = stats.ttest_ind(case_values, control_values, equal_var=False, nan_policy="omit")
        p_value = float(stat.pvalue) if not math.isnan(float(stat.pvalue)) else 1.0
        ann = annotation.get(probe_id, {})
        rows.append(
            {
                "accession": accession,
                "probe_id": probe_id,
                "gene_symbol": ann.get("gene_symbol", ""),
                "gene_title": ann.get("gene_title", ""),
                "gene_id": ann.get("gene_id", ""),
                "n_case": str(len(case_values)),
                "n_control": str(len(control_values)),
                "mean_case": f"{mean_case:.6g}",
                "mean_control": f"{mean_control:.6g}",
                "log2fc": f"{mean_case - mean_control:.6g}",
                "t_stat": f"{float(stat.statistic):.6g}",
                "p_value": f"{p_value:.6g}",
            }
        )

    fdr_values = bh_fdr([float(row["p_value"]) for row in rows])
    for row, fdr in zip(rows, fdr_values):
        row["fdr"] = f"{fdr:.6g}"
    rows.sort(key=lambda row: (float(row["fdr"]), -abs(float(row["log2fc"]))))

    probe_fields = [
        "accession",
        "probe_id",
        "gene_symbol",
        "gene_title",
        "gene_id",
        "n_case",
        "n_control",
        "mean_case",
        "mean_control",
        "log2fc",
        "t_stat",
        "p_value",
        "fdr",
    ]
    write_tsv(out_dir / f"{accession}_probe_signature.tsv", rows, probe_fields)

    gene_rows = collapse_to_gene(rows)
    write_tsv(
        out_dir / f"{accession}_gene_signature.tsv",
        gene_rows,
        ["accession", "gene_symbol", "gene_title", "gene_id", "best_probe_id", "log2fc", "t_stat", "p_value", "fdr", "probe_count"],
    )
    audit = {
        "accession": accession,
        "case_samples": str(len(case_idx)),
        "control_samples": str(len(control_idx)),
        "excluded_samples": str(excluded),
        "probe_rows": str(len(rows)),
        "gene_rows": str(len(gene_rows)),
        "case_group_rule": DATASETS[accession]["contrast"],
        "evidence_group": DATASETS[accession]["evidence_group"],
        "role": DATASETS[accession]["role"],
    }
    return gene_rows, audit


def collapse_to_gene(probe_rows: List[Dict[str, str]]) -> List[Dict[str, str]]:
    grouped: Dict[str, List[Dict[str, str]]] = defaultdict(list)
    for row in probe_rows:
        gene = row.get("gene_symbol", "")
        if not gene:
            continue
        grouped[gene].append(row)
    gene_rows: List[Dict[str, str]] = []
    for gene, rows in grouped.items():
        best = min(rows, key=lambda row: (float(row["fdr"]), -abs(float(row["log2fc"]))))
        gene_rows.append(
            {
                "accession": best["accession"],
                "gene_symbol": gene,
                "gene_title": best.get("gene_title", ""),
                "gene_id": best.get("gene_id", ""),
                "best_probe_id": best["probe_id"],
                "log2fc": best["log2fc"],
                "t_stat": best["t_stat"],
                "p_value": best["p_value"],
                "fdr": best["fdr"],
                "probe_count": str(len(rows)),
            }
        )
    gene_rows.sort(key=lambda row: (float(row["fdr"]), -abs(float(row["log2fc"]))))
    return gene_rows


def write_tsv(path: Path, rows: List[Dict[str, str]], fields: List[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def build_consensus(all_gene_rows: Dict[str, List[Dict[str, str]]], out_dir: Path) -> None:
    by_gene: Dict[str, List[Dict[str, str]]] = defaultdict(list)
    for rows in all_gene_rows.values():
        for row in rows:
            by_gene[row["gene_symbol"]].append(row)
    consensus: List[Dict[str, str]] = []
    for gene, rows in by_gene.items():
        rows_by_group: Dict[str, List[Dict[str, str]]] = defaultdict(list)
        for row in rows:
            rows_by_group[DATASETS[row["accession"]]["evidence_group"]].append(row)
        representative_rows = [min(group_rows, key=lambda row: (float(row["fdr"]), -abs(float(row["log2fc"])))) for group_rows in rows_by_group.values()]
        if len(representative_rows) < 2:
            continue
        lfc = [float(row["log2fc"]) for row in representative_rows]
        up = sum(1 for value in lfc if value > 0)
        down = sum(1 for value in lfc if value < 0)
        direction = "up" if up >= down else "down"
        consistency = max(up, down) / len(representative_rows)
        if consistency < 1.0:
            continue
        mean_lfc = float(np.mean(lfc))
        signed_score = abs(mean_lfc) * consistency * len(representative_rows)
        best_fdr = min(float(row["fdr"]) for row in representative_rows)
        consensus.append(
            {
                "gene_symbol": gene,
                "gene_title": next((row["gene_title"] for row in representative_rows if row.get("gene_title")), ""),
                "datasets": ";".join(row["accession"] for row in rows),
                "evidence_groups": ";".join(sorted(rows_by_group)),
                "evidence_group_count": str(len(representative_rows)),
                "direction": direction,
                "direction_consistency": f"{consistency:.3f}",
                "mean_log2fc": f"{mean_lfc:.6g}",
                "min_log2fc": f"{min(lfc):.6g}",
                "max_log2fc": f"{max(lfc):.6g}",
                "best_fdr": f"{best_fdr:.6g}",
                "consensus_score": f"{signed_score:.6g}",
            }
        )
    consensus.sort(key=lambda row: (-float(row["consensus_score"]), float(row["best_fdr"])))
    write_tsv(
        out_dir / "consensus_gene_signature.tsv",
        consensus,
        [
            "gene_symbol",
            "gene_title",
            "datasets",
            "evidence_groups",
            "evidence_group_count",
            "direction",
            "direction_consistency",
            "mean_log2fc",
            "min_log2fc",
            "max_log2fc",
            "best_fdr",
            "consensus_score",
        ],
    )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--metadata", default="results/m1_uc_ibd_metadata_audit/geo_sample_metadata.tsv")
    parser.add_argument("--raw-dir", default="data/raw/geo")
    parser.add_argument("--out-dir", default="results/m4_geo_uc_ibd_signature")
    parser.add_argument("--manifest", default="data/manifest.tsv")
    parser.add_argument("--timeout", type=int, default=1200)
    parser.add_argument("--accessions", nargs="*", default=list(DATASETS))
    args = parser.parse_args()

    raw_dir = Path(args.raw_dir)
    out_dir = Path(args.out_dir)
    metadata = read_metadata(Path(args.metadata))
    manifest_rows: List[Dict[str, str]] = []
    all_gene_rows: Dict[str, List[Dict[str, str]]] = {}
    audit_rows: List[Dict[str, str]] = []

    for accession in args.accessions:
        if accession not in DATASETS:
            raise ValueError(f"Unsupported accession: {accession}")
        platform = DATASETS[accession]["platform"]
        matrix_path = raw_dir / accession / f"{accession}_series_matrix.txt.gz"
        platform_path = raw_dir / "platforms" / f"{platform}.annot.gz"

        if not matrix_path.exists():
            print(f"[geo-signature] downloading {accession} series matrix", file=sys.stderr)
            run_curl(series_matrix_url(accession), matrix_path, args.timeout)
        if not platform_path.exists():
            print(f"[geo-signature] downloading {platform} annotation", file=sys.stderr)
            run_curl(platform_annot_url(platform), platform_path, args.timeout)

        manifest_rows.append(manifest_row(f"geo_{accession.lower()}_series_matrix", "NCBI GEO", accession, matrix_path, series_matrix_url(accession), "analysis_ready", "GEO series matrix used for UC/IBD bulk expression signature"))
        manifest_rows.append(manifest_row(f"geo_{platform.lower()}_annotation", "NCBI GEO GPL", platform, platform_path, platform_annot_url(platform), "analysis_ready", "GEO platform annotation used for probe-to-gene mapping"))

        print(f"[geo-signature] parsing annotation {platform}", file=sys.stderr)
        annotation = load_platform_annotation(platform_path)
        print(f"[geo-signature] building {accession}", file=sys.stderr)
        gene_rows, audit = build_signature(accession, matrix_path, metadata, annotation, out_dir)
        all_gene_rows[accession] = gene_rows
        audit_rows.append(audit)

    write_tsv(out_dir / "group_audit.tsv", audit_rows, ["accession", "role", "evidence_group", "case_samples", "control_samples", "excluded_samples", "probe_rows", "gene_rows", "case_group_rule"])
    build_consensus(all_gene_rows, out_dir)
    update_manifest(Path(args.manifest), manifest_rows)
    return 0


def manifest_row(dataset_id: str, source: str, version: str, path: Path, url: str, status: str, notes: str) -> Dict[str, str]:
    return {
        "dataset_id": dataset_id,
        "source": source,
        "version_or_date": version,
        "file_path": str(path),
        "source_url_or_api": url,
        "access_date": "2026-06-29",
        "file_size_bytes": str(path.stat().st_size),
        "md5": md5sum(path),
        "status": status,
        "notes": notes,
    }


if __name__ == "__main__":
    raise SystemExit(main())
