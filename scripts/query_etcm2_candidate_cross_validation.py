#!/usr/bin/env python3
"""Query ETCM2 candidate ingredient details and compare targets with HERB/M9/M10."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
import subprocess
import time
import urllib.error
import urllib.parse
import urllib.request
from collections import defaultdict
from pathlib import Path
from typing import Any, Dict, Iterable, List, Sequence, Set, Tuple


BROWSE_API = "http://www.tcmip.cn:18124/home/browse/"
DETAIL_API = "http://www.tcmip.cn:18124/home/detail/"
ACCESS_DATE = "2026-06-29"
RAW_DIR = Path("data/raw/etcm2")
OUT_DIR = Path("results/m11_etcm2_cross_validation")
DOC_PATH = Path("docs/M11_ETCM2_CROSS_VALIDATION_RUN.md")
MANIFEST_PATH = Path("data/manifest.tsv")
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


def read_tsv(path: Path) -> List[Dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: Sequence[Dict[str, Any]], fields: Sequence[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(fields), delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: "" if row.get(field) is None else row.get(field, "") for field in fields})


def md5_file(path: Path) -> str:
    digest = hashlib.md5()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def opener() -> urllib.request.OpenerDirector:
    return urllib.request.build_opener(urllib.request.ProxyHandler({}))


def normalize_name(value: str) -> str:
    return re.sub(r"[^a-z0-9]+", "", value.lower())


def clean_id(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", value).strip("_") or "unknown"


def as_list(value: Any) -> List[Any]:
    if value is None:
        return []
    if isinstance(value, list):
        return value
    return [value]


def first_list_value(value: Any) -> str:
    values = as_list(value)
    return str(values[0]) if values else ""


def post_json(url: str, payload: Dict[str, Any], timeout: int) -> Dict[str, Any]:
    result = subprocess.run(
        [
            "curl",
            "-L",
            "-sS",
            "--fail",
            "--connect-timeout",
            "20",
            "--max-time",
            str(timeout),
            "--retry",
            "1",
            "--noproxy",
            "www.tcmip.cn",
            url,
            "-H",
            "Content-Type: application/json",
            "--data",
            json.dumps(payload, ensure_ascii=False),
        ],
        check=True,
        capture_output=True,
        text=True,
    )
    return json.loads(result.stdout)


def get_json(url: str, params: Dict[str, str], timeout: int) -> Dict[str, Any]:
    command = [
        "curl",
        "-L",
        "-sS",
        "--fail",
        "--connect-timeout",
        "20",
        "--max-time",
        str(timeout),
        "--retry",
        "1",
        "--noproxy",
        "www.tcmip.cn",
        "--get",
        url,
    ]
    for key, value in params.items():
        command.extend(["--data-urlencode", f"{key}={value}"])
    result = subprocess.run(command, check=True, capture_output=True, text=True)
    if not result.stdout.strip():
        raise RuntimeError("empty response body")
    return json.loads(result.stdout)


def retry_call(func, retries: int, sleep_seconds: float) -> Tuple[str, Any]:
    last_error = ""
    for attempt in range(retries + 1):
        try:
            return "ok", func()
        except (urllib.error.URLError, TimeoutError, json.JSONDecodeError, RuntimeError, OSError, subprocess.CalledProcessError) as exc:
            last_error = f"{type(exc).__name__}: {exc}"
            if attempt < retries and sleep_seconds:
                time.sleep(sleep_seconds)
    return last_error, None


def candidate_rows(m9_path: Path, m10_path: Path, limit: int) -> List[Dict[str, str]]:
    merged: Dict[str, Dict[str, str]] = {}
    for rank_source, path in [("m9_bridge", m9_path), ("m10_open_targets", m10_path)]:
        if not path.exists():
            continue
        for rank, row in enumerate(read_tsv(path), start=1):
            ingredient_id = row.get("ingredient_id", "")
            if not ingredient_id:
                continue
            entry = merged.setdefault(
                ingredient_id,
                {
                    "ingredient_id": ingredient_id,
                    "ingredient_name": row.get("ingredient_name", ""),
                    "candidate_priority_class": row.get("candidate_priority_class", ""),
                    "candidate_priority_score": row.get("candidate_priority_score", ""),
                    "m9_rank": "",
                    "m10_rank": "",
                },
            )
            if rank_source == "m9_bridge":
                entry["m9_rank"] = str(rank)
                entry["m9_bridge_score"] = row.get("bridge_score", "")
            else:
                entry["m10_rank"] = str(rank)
                entry["m10_genetic_supported_genes"] = row.get("n_genetic_supported_genes", "")
                entry["m10_open_targets_supported_genes"] = row.get("n_open_targets_supported_genes", "")

    def sort_key(row: Dict[str, str]) -> Tuple[int, int, str]:
        m9_rank = int(row["m9_rank"]) if row.get("m9_rank") else 9999
        m10_rank = int(row["m10_rank"]) if row.get("m10_rank") else 9999
        return (min(m9_rank, m10_rank), m9_rank + m10_rank, row["ingredient_name"].lower())

    return sorted(merged.values(), key=sort_key)[:limit]


def browse_ingredient_pages(page_size: int, pages: int, timeout: int, retries: int, sleep_seconds: float) -> Tuple[List[Dict[str, Any]], List[Dict[str, Any]]]:
    rows: List[Dict[str, Any]] = []
    audit: List[Dict[str, Any]] = []
    browse_dir = RAW_DIR / "browse"
    browse_dir.mkdir(parents=True, exist_ok=True)
    for page in range(1, pages + 1):
        raw_path = browse_dir / f"ingredient_page{page}_size{page_size}.json"
        if raw_path.exists() and raw_path.stat().st_size > 100:
            payload = json.loads(raw_path.read_text(encoding="utf-8"))
            status = "cached"
            message = ""
        else:
            status, payload = retry_call(
                lambda page=page: post_json(
                    BROWSE_API,
                    {"type": "ingredient", "pageNo": page, "pageSize": page_size, "language": "en", "search_key": ""},
                    timeout,
                ),
                retries,
                sleep_seconds,
            )
            message = "" if status == "ok" else status
            if status == "ok":
                raw_path.write_text(json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8")
                status = "downloaded"
        root = (payload.get("data") or [{}])[0] if isinstance(payload, dict) else {}
        page_rows = root.get("data") or []
        rows.extend(page_rows)
        audit.append(
            {
                "query_type": "browse_ingredient",
                "query_id": f"page={page};pageSize={page_size}",
                "status": status,
                "raw_path": str(raw_path) if raw_path.exists() else "",
                "n_rows": len(page_rows),
                "message": message,
            }
        )
        if payload is None or not page_rows:
            break
    return rows, audit


def search_candidate_ingredient(
    candidate: Dict[str, str],
    page_size: int,
    timeout: int,
    retries: int,
    sleep_seconds: float,
    reuse_existing: bool,
) -> Tuple[List[Dict[str, Any]], Dict[str, Any]]:
    search_dir = RAW_DIR / "ingredient_search"
    search_dir.mkdir(parents=True, exist_ok=True)
    query = candidate["ingredient_name"]
    raw_path = search_dir / f"{clean_id(query)}.json"
    if reuse_existing and raw_path.exists() and raw_path.stat().st_size > 100:
        payload = json.loads(raw_path.read_text(encoding="utf-8"))
        status = "cached"
        message = ""
    else:
        status, payload = retry_call(
            lambda: post_json(
                BROWSE_API,
                {"type": "ingredient", "pageNo": 1, "pageSize": page_size, "language": "en", "search_key": query},
                timeout,
            ),
            retries,
            sleep_seconds,
        )
        message = "" if status == "ok" else status
        if status == "ok":
            raw_path.write_text(json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8")
            status = "downloaded"
    root = (payload.get("data") or [{}])[0] if isinstance(payload, dict) else {}
    rows = root.get("data") or []
    return rows, {
        **candidate,
        "query_type": "ingredient_search",
        "query_id": query,
        "status": status,
        "raw_path": str(raw_path) if raw_path.exists() else "",
        "n_rows": len(rows),
        "message": message,
    }


def detail_id_candidates(candidate_name: str, browse_rows: Sequence[Dict[str, Any]]) -> List[str]:
    wanted = normalize_name(candidate_name)
    ids: List[str] = []
    for row in browse_rows:
        ingredient_name = first_list_value(row.get("Ingredient Name"))
        if normalize_name(ingredient_name) == wanted:
            ids.append(ingredient_name)
    ids.extend(
        [
            candidate_name,
            candidate_name.replace(" ", ""),
            candidate_name.replace(" ", "-"),
            candidate_name.lower(),
        ]
    )
    deduped: List[str] = []
    for value in ids:
        if value and value not in deduped:
            deduped.append(value)
    return deduped


def download_detail(detail_id: str, timeout: int, retries: int, sleep_seconds: float, reuse_existing: bool) -> Tuple[str, Dict[str, Any] | None, Path, str]:
    detail_dir = RAW_DIR / "ingredient_details"
    detail_dir.mkdir(parents=True, exist_ok=True)
    raw_path = detail_dir / f"{clean_id(detail_id)}.json"
    if reuse_existing and raw_path.exists() and raw_path.stat().st_size > 100:
        return "cached", json.loads(raw_path.read_text(encoding="utf-8")), raw_path, ""

    status, payload = retry_call(
        lambda: get_json(DETAIL_API, {"id": detail_id, "type": "ingredient", "language": "en"}, timeout),
        retries,
        sleep_seconds,
    )
    if status != "ok":
        return "failed", None, raw_path, status
    if not isinstance(payload, dict) or payload.get("code") != 1:
        return "failed", payload if isinstance(payload, dict) else None, raw_path, f"unexpected payload code: {payload}"
    raw_path.write_text(json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8")
    return "downloaded", payload, raw_path, ""


def extract_section(payload: Dict[str, Any], section_id: str) -> Any:
    for section in payload.get("data") or []:
        if section.get("id") == section_id:
            return section.get("value")
    return None


def extract_etcm_edges(candidate: Dict[str, str], payload: Dict[str, Any], detail_id: str) -> List[Dict[str, Any]]:
    rows: List[Dict[str, Any]] = []
    value = extract_section(payload, "ingredient_target")
    if isinstance(value, dict):
        for row in value.get("data") or []:
            gene = first_list_value(row.get("Gene Symbol")).upper()
            if not gene:
                continue
            linkformat = row.get("linkformat") or {}
            reference_links = ((linkformat.get("References") or {}).get("param") or []) if isinstance(linkformat, dict) else []
            rows.append(
                {
                    "ingredient_id": candidate["ingredient_id"],
                    "ingredient_name": candidate["ingredient_name"],
                    "etcm_detail_id": detail_id,
                    "etcm_ingredient_name": first_list_value(row.get("Ingredient Name")),
                    "edge_source": "ingredient_target_table",
                    "etcm_target_name": gene,
                    "gene_symbol": gene,
                    "similar_score": row.get("Similar Score", ""),
                    "activity": row.get("Activity", ""),
                    "reference_labels": ";".join(str(value) for value in as_list(row.get("References"))),
                    "reference_links": ";".join(str(value) for value in as_list(reference_links)),
                }
            )

    network = extract_section(payload, "basic_network")
    if isinstance(network, dict):
        target_node_ids = {
            node.get("id")
            for node in network.get("nodes", [])
            if ((node.get("color") or {}).get("background") == "#9f7733")
        }
        for edge in network.get("edges", []):
            source = edge.get("from")
            target = edge.get("to")
            if target in target_node_ids:
                target_name = str(target)
            elif source in target_node_ids:
                target_name = str(source)
            else:
                continue
            rows.append(
                {
                    "ingredient_id": candidate["ingredient_id"],
                    "ingredient_name": candidate["ingredient_name"],
                    "etcm_detail_id": detail_id,
                    "etcm_ingredient_name": detail_id,
                    "edge_source": "basic_network_target_node",
                    "etcm_target_name": target_name,
                    "gene_symbol": "",
                    "similar_score": "",
                    "activity": "",
                    "reference_labels": "",
                    "reference_links": "",
                }
            )
    return rows


def load_gene_sets(m9_edge_path: Path, m10_path: Path) -> Tuple[Dict[str, Set[str]], Dict[str, Set[str]], Dict[str, Set[str]]]:
    herb_m9: Dict[str, Set[str]] = defaultdict(set)
    bulk_aligned: Dict[str, Set[str]] = defaultdict(set)
    for row in read_tsv(m9_edge_path):
        ingredient_id = row["ingredient_id"]
        gene = row["gene_symbol"].upper()
        herb_m9[ingredient_id].add(gene)
        if row.get("therapeutic_alignment", "").startswith("aligned"):
            bulk_aligned[ingredient_id].add(gene)
    open_targets_genetic: Dict[str, Set[str]] = defaultdict(set)
    if m10_path.exists():
        for row in read_tsv(m10_path):
            ingredient_id = row["ingredient_id"]
            open_targets_genetic[ingredient_id].update(gene for gene in row.get("genetic_supported_genes", "").split(";") if gene)
    return herb_m9, bulk_aligned, open_targets_genetic


def summarize(
    candidates: Sequence[Dict[str, str]],
    edges: Sequence[Dict[str, Any]],
    audit_rows: Sequence[Dict[str, Any]],
    herb_m9: Dict[str, Set[str]],
    bulk_aligned: Dict[str, Set[str]],
    open_targets_genetic: Dict[str, Set[str]],
) -> List[Dict[str, Any]]:
    by_candidate_edges: Dict[str, List[Dict[str, Any]]] = defaultdict(list)
    for edge in edges:
        by_candidate_edges[edge["ingredient_id"]].append(edge)
    by_candidate_audit: Dict[str, List[Dict[str, Any]]] = defaultdict(list)
    for row in audit_rows:
        by_candidate_audit[row.get("ingredient_id", "")].append(row)

    rows: List[Dict[str, Any]] = []
    for candidate in candidates:
        ingredient_id = candidate["ingredient_id"]
        etcm_targets = {edge.get("etcm_target_name") or edge.get("gene_symbol") for edge in by_candidate_edges.get(ingredient_id, [])}
        etcm_genes = {edge["gene_symbol"] for edge in by_candidate_edges.get(ingredient_id, []) if edge.get("gene_symbol")}
        rows.append(
            {
                **candidate,
                "query_statuses": ";".join(sorted({row.get("status", "") for row in by_candidate_audit.get(ingredient_id, []) if row.get("status", "")})),
                "n_etcm2_target_edges_extracted": len(by_candidate_edges.get(ingredient_id, [])),
                "n_etcm2_unique_targets_extracted": len(etcm_targets),
                "etcm2_unique_targets_extracted": ";".join(sorted(etcm_targets)),
                "n_overlap_with_herb_m9_disease_context_targets": len(etcm_genes & herb_m9.get(ingredient_id, set())),
                "overlap_with_herb_m9_disease_context_targets": ";".join(sorted(etcm_genes & herb_m9.get(ingredient_id, set()))),
                "n_overlap_with_bulk_aligned_targets": len(etcm_genes & bulk_aligned.get(ingredient_id, set())),
                "overlap_with_bulk_aligned_targets": ";".join(sorted(etcm_genes & bulk_aligned.get(ingredient_id, set()))),
                "n_overlap_with_open_targets_genetic_supported": len(etcm_genes & open_targets_genetic.get(ingredient_id, set())),
                "overlap_with_open_targets_genetic_supported": ";".join(sorted(etcm_genes & open_targets_genetic.get(ingredient_id, set()))),
            }
        )
    rows.sort(
        key=lambda row: (
            -int(row["n_overlap_with_herb_m9_disease_context_targets"]),
            -int(row["n_etcm2_unique_targets_extracted"]),
            int(row["m9_rank"]) if row.get("m9_rank") else 9999,
        )
    )
    return rows


def update_manifest(success_paths: Sequence[Tuple[str, str, Path]]) -> None:
    existing: List[Dict[str, str]] = []
    if MANIFEST_PATH.exists():
        with MANIFEST_PATH.open(encoding="utf-8", newline="") as handle:
            existing = list(csv.DictReader(handle, delimiter="\t"))
    by_id = {row["dataset_id"]: row for row in existing if row.get("dataset_id")}
    for entity_type, entity_id, path in success_paths:
        by_id[f"etcm2_{entity_type}_{clean_id(entity_id).lower()}"] = {
            "dataset_id": f"etcm2_{entity_type}_{clean_id(entity_id).lower()}",
            "source": "ETCM2 public API",
            "version_or_date": ACCESS_DATE,
            "file_path": str(path),
            "source_url_or_api": DETAIL_API if entity_type == "ingredient_detail" else BROWSE_API,
            "access_date": ACCESS_DATE,
            "file_size_bytes": str(path.stat().st_size),
            "md5": md5_file(path),
            "status": "analysis_ready",
            "notes": f"ETCM2 {entity_type} payload for {entity_id}; acquired for candidate cross-validation",
        }
    with MANIFEST_PATH.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=MANIFEST_FIELDS, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for key in sorted(by_id):
            writer.writerow(by_id[key])


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--m9-summary", default="results/m9_target_cell_bridge/candidate_target_cell_bridge_summary.tsv")
    parser.add_argument("--m9-edges", default="results/m9_target_cell_bridge/disease_context_herb_target_edges.tsv")
    parser.add_argument("--m10-summary", default="results/m10_open_targets/candidate_open_targets_summary.tsv")
    parser.add_argument("--max-candidates", type=int, default=12)
    parser.add_argument("--browse-pages", type=int, default=3)
    parser.add_argument("--page-size", type=int, default=500)
    parser.add_argument("--search-page-size", type=int, default=25)
    parser.add_argument("--timeout", type=int, default=20)
    parser.add_argument("--retries", type=int, default=1)
    parser.add_argument("--sleep-seconds", type=float, default=1.0)
    parser.add_argument("--max-detail-ids-per-candidate", type=int, default=4)
    parser.add_argument("--no-reuse-existing", action="store_true")
    args = parser.parse_args()

    candidates = candidate_rows(Path(args.m9_summary), Path(args.m10_summary), args.max_candidates)
    herb_m9, bulk_aligned, open_targets_genetic = load_gene_sets(Path(args.m9_edges), Path(args.m10_summary))
    browse_rows, browse_audit = browse_ingredient_pages(args.page_size, args.browse_pages, args.timeout, args.retries, args.sleep_seconds)

    audit_rows: List[Dict[str, Any]] = list(browse_audit)
    success_paths: List[Tuple[str, str, Path]] = []
    for row in browse_audit:
        if row.get("status") in {"cached", "downloaded"} and row.get("raw_path"):
            success_paths.append(("ingredient_browse", clean_id(row["query_id"]), Path(row["raw_path"])))

    edge_rows: List[Dict[str, Any]] = []
    for candidate in candidates:
        search_rows, search_audit = search_candidate_ingredient(
            candidate,
            page_size=args.search_page_size,
            timeout=args.timeout,
            retries=args.retries,
            sleep_seconds=args.sleep_seconds,
            reuse_existing=not args.no_reuse_existing,
        )
        audit_rows.append(search_audit)
        if search_audit.get("status") in {"cached", "downloaded"} and search_audit.get("raw_path"):
            success_paths.append(("ingredient_search", candidate["ingredient_name"], Path(search_audit["raw_path"])))
        detail_ids = detail_id_candidates(candidate["ingredient_name"], search_rows or browse_rows)[: args.max_detail_ids_per_candidate]
        candidate_success = False
        candidate_messages: List[str] = []
        for detail_id in detail_ids:
            status, payload, raw_path, message = download_detail(
                detail_id,
                timeout=args.timeout,
                retries=args.retries,
                sleep_seconds=args.sleep_seconds,
                reuse_existing=not args.no_reuse_existing,
            )
            audit_rows.append(
                {
                    **candidate,
                    "query_type": "ingredient_detail",
                    "query_id": detail_id,
                    "status": status,
                    "raw_path": str(raw_path) if raw_path.exists() else "",
                    "n_rows": "",
                    "message": message,
                }
            )
            if status in {"cached", "downloaded"} and payload:
                candidate_edges = extract_etcm_edges(candidate, payload, detail_id)
                edge_rows.extend(candidate_edges)
                success_paths.append(("ingredient_detail", detail_id, raw_path))
                candidate_success = True
                break
            candidate_messages.append(f"{detail_id}: {message}")
        if not candidate_success and not detail_ids:
            audit_rows.append(
                {
                    **candidate,
                    "query_type": "ingredient_detail",
                    "query_id": "",
                    "status": "not_matched",
                    "raw_path": "",
                    "n_rows": "",
                    "message": "No ETCM2 browse-name match or direct candidate ID available",
                }
            )

    update_manifest(success_paths)
    summary_rows = summarize(candidates, edge_rows, audit_rows, herb_m9, bulk_aligned, open_targets_genetic)

    audit_fields = [
        "ingredient_id",
        "ingredient_name",
        "m9_rank",
        "m10_rank",
        "query_type",
        "query_id",
        "status",
        "raw_path",
        "n_rows",
        "message",
    ]
    edge_fields = [
        "ingredient_id",
        "ingredient_name",
        "etcm_detail_id",
        "etcm_ingredient_name",
        "edge_source",
        "etcm_target_name",
        "gene_symbol",
        "similar_score",
        "activity",
        "reference_labels",
        "reference_links",
    ]
    summary_fields = [
        "ingredient_id",
        "ingredient_name",
        "candidate_priority_class",
        "candidate_priority_score",
        "m9_rank",
        "m9_bridge_score",
        "m10_rank",
        "m10_genetic_supported_genes",
        "m10_open_targets_supported_genes",
        "query_statuses",
        "n_etcm2_target_edges_extracted",
        "n_etcm2_unique_targets_extracted",
        "etcm2_unique_targets_extracted",
        "n_overlap_with_herb_m9_disease_context_targets",
        "overlap_with_herb_m9_disease_context_targets",
        "n_overlap_with_bulk_aligned_targets",
        "overlap_with_bulk_aligned_targets",
        "n_overlap_with_open_targets_genetic_supported",
        "overlap_with_open_targets_genetic_supported",
    ]

    write_tsv(OUT_DIR / "etcm2_candidate_query_audit.tsv", audit_rows, audit_fields)
    write_tsv(OUT_DIR / "etcm2_ingredient_target_edges.tsv", edge_rows, edge_fields)
    write_tsv(OUT_DIR / "candidate_etcm2_cross_validation_summary.tsv", summary_rows, summary_fields)

    downloaded_details = [row for row in audit_rows if row.get("query_type") == "ingredient_detail" and row.get("status") in {"cached", "downloaded"}]
    failed_details = [row for row in audit_rows if row.get("query_type") == "ingredient_detail" and row.get("status") == "failed"]
    report_lines = [
        "# M11 ETCM2 Candidate Cross-Validation Run",
        "",
        "## Scope",
        "",
        "This run attempts to query ETCM2 ingredient detail payloads for the top M9/M10 candidates, extract ETCM2 `ingredient_target` rows when available, extract `basic_network` target nodes for ingredient details, and compare gene-symbol-level targets with HERB disease-context targets, bulk-aligned targets, and Open Targets genetic-support genes where possible.",
        "",
        "## Key Outputs",
        "",
        "- `results/m11_etcm2_cross_validation/etcm2_candidate_query_audit.tsv`",
        "- `results/m11_etcm2_cross_validation/etcm2_ingredient_target_edges.tsv`",
        "- `results/m11_etcm2_cross_validation/candidate_etcm2_cross_validation_summary.tsv`",
        "",
        "## Acquisition Summary",
        "",
        f"- Candidate ingredients attempted: {len(candidates)}.",
        f"- Browse pages requested: {args.browse_pages} with page size {args.page_size}.",
        f"- Successful/cached candidate detail payloads: {len(downloaded_details)}.",
        f"- Failed candidate detail attempts: {len(failed_details)}.",
        f"- ETCM2 target edge rows extracted: {len(edge_rows)}.",
        "",
        "## Top Cross-Validation Summary",
        "",
        "| Ingredient | Query statuses | ETCM2 targets | HERB/M9 overlap | Bulk-aligned overlap | Open Targets genetic overlap |",
        "|---|---|---:|---:|---:|---:|",
    ]
    for row in summary_rows[:12]:
        report_lines.append(
            f"| {row['ingredient_name']} | {row['query_statuses']} | {row['n_etcm2_unique_targets_extracted']} | {row['n_overlap_with_herb_m9_disease_context_targets']} | {row['n_overlap_with_bulk_aligned_targets']} | {row['n_overlap_with_open_targets_genetic_supported']} |"
        )
    report_lines.extend(
        [
            "",
            "## Interpretation Boundary",
            "",
            "- ETCM2 detail `ingredient_target` rows and `basic_network` target nodes are treated as complementary target annotations, not proof of compound efficacy.",
            "- If the ETCM2 backend times out, the audit table records failed attempts and the command can be rerun later with the same inputs.",
            "- The script uses curl with `--noproxy www.tcmip.cn`, because this endpoint previously stalled through local proxy routing.",
            "- Ingredient-detail target nodes are target names, not necessarily gene symbols; do not claim gene-level overlap unless a target-name mapping is added.",
            "- The detail payload can expose a limited target preview; use extracted rows as cross-validation signals, not an exhaustive target universe unless full ETCM2 export access is obtained.",
            "",
        ]
    )
    DOC_PATH.write_text("\n".join(report_lines), encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
