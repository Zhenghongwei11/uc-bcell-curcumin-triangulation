#!/usr/bin/env python3
"""Query L1000CDS2 for UC/IBD reversal hits among LINCS-covered TCM candidates."""

from __future__ import annotations

import argparse
import csv
import json
import re
import sys
import time
import urllib.error
import urllib.request
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, Iterable, List, Sequence, Tuple


DEFAULT_URL = "https://maayanlab.cloud/L1000CDS2/query2"
DEFAULT_DB_VERSION = "cpcd-gse70138-v1.0,cpcd-gse70138-lm-v1.0"


def read_genes(path: Path, limit: int) -> List[str]:
    genes: List[str] = []
    seen = set()
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            gene = line.strip().upper()
            if not gene or gene in seen:
                continue
            genes.append(gene)
            seen.add(gene)
            if len(genes) >= limit:
                break
    return genes


def read_tsv(path: Path) -> Iterable[Dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        yield from csv.DictReader(handle, delimiter="\t")


def write_tsv(path: Path, rows: Sequence[Dict[str, Any]], fields: Sequence[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(fields), delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: "" if row.get(field) is None else row.get(field, "") for field in fields})


def normalize_name(value: str) -> str:
    return re.sub(r"[^a-z0-9]+", " ", (value or "").lower()).strip()


def build_payload(up_genes: Sequence[str], down_genes: Sequence[str], db_version: str, tag: str) -> Dict[str, Any]:
    return {
        "data": {
            "upGenes": list(up_genes),
            "dnGenes": list(down_genes),
        },
        "config": {
            "aggravate": False,
            "share": False,
            "combination": False,
            "db-version": db_version,
            "includeLessSignificant": False,
            "searchMethod": "geneSet",
        },
        "meta": [
            {
                "key": "Tag",
                "value": tag,
            }
        ],
    }


def post_json(url: str, payload: Dict[str, Any], timeout: int) -> Tuple[int, str, str]:
    data = json.dumps(payload, separators=(",", ":")).encode("utf-8")
    request = urllib.request.Request(
        url,
        data=data,
        method="POST",
        headers={
            "Content-Type": "application/json",
            "Accept": "application/json,text/plain,*/*",
            "User-Agent": "zyy-public-data-research/1.0",
        },
    )
    try:
        with urllib.request.urlopen(request, timeout=timeout) as response:
            return response.status, response.read().decode("utf-8", errors="replace"), ""
    except urllib.error.HTTPError as exc:
        body = exc.read().decode("utf-8", errors="replace")
        return exc.code, body, str(exc)
    except urllib.error.URLError as exc:
        return 0, "", str(exc)


def find_top_meta(response: Any) -> List[Dict[str, Any]]:
    if isinstance(response, dict):
        direct = response.get("topMeta")
        if isinstance(direct, list):
            return [item for item in direct if isinstance(item, dict)]
        nested = response.get("result")
        if isinstance(nested, dict):
            nested_meta = nested.get("topMeta")
            if isinstance(nested_meta, list):
                return [item for item in nested_meta if isinstance(item, dict)]
        for key in ("entries", "results", "data"):
            value = response.get(key)
            if isinstance(value, list) and all(isinstance(item, dict) for item in value):
                if any("pert_id" in item or "pert_desc" in item or "sig_id" in item for item in value):
                    return list(value)
    if isinstance(response, list) and all(isinstance(item, dict) for item in response):
        if any("pert_id" in item or "pert_desc" in item or "sig_id" in item for item in response):
            return list(response)
    return []


def safe_float(value: Any) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return 0.0


def load_covered_candidates(path: Path) -> Tuple[Dict[str, List[Dict[str, str]]], Dict[str, List[Dict[str, str]]], List[Dict[str, str]]]:
    by_pert: Dict[str, List[Dict[str, str]]] = {}
    by_name: Dict[str, List[Dict[str, str]]] = {}
    by_candidate_key: Dict[Tuple[str, str], Dict[str, str]] = {}
    for row in read_tsv(path):
        if row.get("priority_class") != "priority_for_lincs_mapping":
            continue
        if row.get("match_status") != "covered":
            continue
        pert_id = row.get("pert_id", "").strip()
        if not pert_id:
            continue
        candidate = {
            "priority_score": row.get("priority_score", ""),
            "ingredient_id": row.get("ingredient_id", ""),
            "ingredient_name": row.get("ingredient_name", ""),
            "pubchem_id": row.get("pubchem_id", ""),
            "inchikey": row.get("inchikey", ""),
            "lincs_source": row.get("lincs_source", ""),
            "pert_id": pert_id,
            "pert_iname": row.get("pert_iname", ""),
        }
        by_pert.setdefault(pert_id.upper(), []).append(candidate)
        for name in (candidate["ingredient_name"], candidate["pert_iname"]):
            normalized = normalize_name(name)
            if normalized:
                by_name.setdefault(normalized, []).append(candidate)
        key = (candidate["ingredient_id"], candidate["pert_id"])
        if key not in by_candidate_key:
            by_candidate_key[key] = {**candidate, "lincs_sources": candidate["lincs_source"]}
        else:
            sources = set(filter(None, by_candidate_key[key]["lincs_sources"].split(";")))
            sources.add(candidate["lincs_source"])
            by_candidate_key[key]["lincs_sources"] = ";".join(sorted(sources))
    candidates = sorted(
        by_candidate_key.values(),
        key=lambda row: (-safe_float(row["priority_score"]), row["ingredient_name"].lower(), row["pert_id"]),
    )
    return by_pert, by_name, candidates


def safe_filename(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", value)


def result_records(top_meta: Sequence[Dict[str, Any]], db_version: str) -> List[Dict[str, Any]]:
    records: List[Dict[str, Any]] = []
    for idx, item in enumerate(top_meta, start=1):
        record = {
            "db_version": db_version,
            "rank": item.get("rank", idx),
            "score": item.get("score", ""),
            "pert_id": item.get("pert_id", ""),
            "pert_desc": item.get("pert_desc", ""),
            "sig_id": item.get("sig_id", ""),
            "cell_id": item.get("cell_id", ""),
            "pert_dose": item.get("pert_dose", ""),
            "pert_dose_unit": item.get("pert_dose_unit", ""),
            "pert_time": item.get("pert_time", ""),
            "pert_time_unit": item.get("pert_time_unit", ""),
            "pubchem_id": item.get("pubchem_id", ""),
            "drugbank_id": item.get("drugbank_id", ""),
        }
        records.append(record)
    records.sort(key=lambda row: int(row["rank"]) if str(row["rank"]).isdigit() else 10**9)
    return records


def match_results(
    records: Sequence[Dict[str, Any]],
    by_pert: Dict[str, List[Dict[str, str]]],
    by_name: Dict[str, List[Dict[str, str]]],
) -> List[Dict[str, Any]]:
    matches: List[Dict[str, Any]] = []
    seen = set()
    for record in records:
        candidates: List[Dict[str, str]] = []
        method = ""
        pert_id = str(record.get("pert_id", "")).strip().upper()
        if pert_id and pert_id in by_pert:
            candidates = by_pert[pert_id]
            method = "pert_id_exact"
        else:
            pert_desc = normalize_name(str(record.get("pert_desc", "")))
            if pert_desc and pert_desc in by_name:
                candidates = by_name[pert_desc]
                method = "pert_desc_exact_normalized"
        for candidate in candidates:
            key = (record.get("rank"), record.get("sig_id"), candidate["ingredient_id"], candidate["pert_id"], method)
            if key in seen:
                continue
            seen.add(key)
            matches.append(
                {
                    **record,
                    "match_method": method,
                    "priority_score": candidate["priority_score"],
                    "ingredient_id": candidate["ingredient_id"],
                    "ingredient_name": candidate["ingredient_name"],
                    "ingredient_pubchem_id": candidate["pubchem_id"],
                    "ingredient_inchikey": candidate["inchikey"],
                    "candidate_lincs_source": candidate["lincs_source"],
                    "candidate_pert_id": candidate["pert_id"],
                    "candidate_pert_iname": candidate["pert_iname"],
                }
            )
    matches.sort(key=lambda row: (int(row["rank"]) if str(row["rank"]).isdigit() else 10**9, -safe_float(row["score"])))
    return matches


def candidate_status_rows(candidates: Sequence[Dict[str, str]], matches: Sequence[Dict[str, Any]]) -> List[Dict[str, Any]]:
    best_by_key: Dict[Tuple[str, str], Dict[str, Any]] = {}
    for match in matches:
        key = (str(match.get("ingredient_id", "")), str(match.get("candidate_pert_id", "")))
        current = best_by_key.get(key)
        rank = int(match["rank"]) if str(match.get("rank", "")).isdigit() else 10**9
        current_rank = int(current["rank"]) if current and str(current.get("rank", "")).isdigit() else 10**9
        if current is None or rank < current_rank:
            best_by_key[key] = match

    rows: List[Dict[str, Any]] = []
    for candidate in candidates:
        key = (candidate["ingredient_id"], candidate["pert_id"])
        hit = best_by_key.get(key)
        rows.append(
            {
                "priority_score": candidate["priority_score"],
                "ingredient_id": candidate["ingredient_id"],
                "ingredient_name": candidate["ingredient_name"],
                "ingredient_pubchem_id": candidate["pubchem_id"],
                "ingredient_inchikey": candidate["inchikey"],
                "candidate_pert_id": candidate["pert_id"],
                "candidate_pert_iname": candidate["pert_iname"],
                "candidate_lincs_sources": candidate.get("lincs_sources", candidate.get("lincs_source", "")),
                "hit_status": "top_result_hit" if hit else "no_top_result_hit",
                "best_db_version": hit.get("db_version", "") if hit else "",
                "best_rank": hit.get("rank", "") if hit else "",
                "best_score": hit.get("score", "") if hit else "",
                "best_sig_id": hit.get("sig_id", "") if hit else "",
                "best_cell_id": hit.get("cell_id", "") if hit else "",
                "best_pert_desc": hit.get("pert_desc", "") if hit else "",
            }
        )
    return rows


def write_report(path: Path, summary: Sequence[Dict[str, Any]], matches: Sequence[Dict[str, Any]]) -> None:
    status = next((row["value"] for row in summary if row["metric"] == "query_status"), "")
    lines = [
        "# M5 LINCS/L1000CDS2 Reversal Run",
        "",
        "## Scope",
        "",
        "This run used the GEO-derived UC/IBD consensus signature and queried the public L1000CDS2 `query2` endpoint in reverse mode. Candidate interpretation is restricted to the HERB-derived priority ingredients with exact LINCS perturbagen identifier coverage.",
        "",
        "## Summary",
        "",
        "| Metric | Value | Notes |",
        "|---|---:|---|",
    ]
    for row in summary:
        value = row.get("value", "")
        notes = row.get("notes", "")
        lines.append(f"| {row.get('metric', '')} | {value} | {notes} |")
    lines.extend(["", "## Candidate Hits", ""])
    if status != "success":
        lines.append("The L1000CDS2 request did not return a parseable success response, so candidate reversal ranking could not be interpreted from this run.")
    elif matches:
        lines.append("| DB version | Rank | Score | Ingredient | LINCS pert_id | Perturbagen | Cell | Dose | Time | Signature |")
        lines.append("|---|---:|---:|---|---|---|---|---|---|---|")
        for row in matches[:50]:
            dose = f"{row.get('pert_dose', '')}{row.get('pert_dose_unit', '')}".strip()
            pert_time = f"{row.get('pert_time', '')}{row.get('pert_time_unit', '')}".strip()
            lines.append(
                "| {db_version} | {rank} | {score} | {ingredient} | {pert_id} | {pert_desc} | {cell} | {dose} | {time} | {sig} |".format(
                    db_version=row.get("db_version", ""),
                    rank=row.get("rank", ""),
                    score=row.get("score", ""),
                    ingredient=row.get("ingredient_name", ""),
                    pert_id=row.get("candidate_pert_id", ""),
                    pert_desc=row.get("pert_desc", ""),
                    cell=row.get("cell_id", ""),
                    dose=dose,
                    time=pert_time,
                    sig=row.get("sig_id", ""),
                )
            )
    else:
        lines.append("No exact-ID HERB priority ingredient appeared in the returned L1000CDS2 top results. This is a negative public-API screen, not evidence that the compounds lack activity in the full LINCS matrix. See `results/m5_lincs_reversal/candidate_reversal_status.tsv` for the tested candidate set.")
    lines.extend(
        [
            "",
            "## Interpretation Boundary",
            "",
            "- A hit means the public L1000CDS2 top-result list contains a LINCS perturbation that maps exactly to a prequalified HERB ingredient.",
            "- A non-hit does not rule out reversal in the full LINCS Level 5 matrix because the public API returns a ranked subset and uses its own database/version filters.",
            "- This module is suitable for prioritization and reviewer-facing triangulation, not for causal or efficacy claims.",
            "",
        ]
    )
    path.write_text("\n".join(lines), encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--up-genes", default="results/m4_geo_uc_ibd_signature/lincs_query_up_genes.txt")
    parser.add_argument("--down-genes", default="results/m4_geo_uc_ibd_signature/lincs_query_down_genes.txt")
    parser.add_argument("--coverage", default="results/m3_lincs_coverage/lincs_candidate_coverage.tsv")
    parser.add_argument("--out-dir", default="results/m5_lincs_reversal")
    parser.add_argument("--report", default="docs/M5_LINCS_REVERSAL_RUN.md")
    parser.add_argument("--url", default=DEFAULT_URL)
    parser.add_argument("--db-version", default=DEFAULT_DB_VERSION, help="Comma-separated L1000CDS2 database version(s).")
    parser.add_argument("--top-n-genes", type=int, default=100)
    parser.add_argument("--tag", default="UC_IBD_consensus_signature")
    parser.add_argument("--timeout", type=int, default=180)
    parser.add_argument("--fail-on-api-error", action="store_true")
    args = parser.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    up_genes = read_genes(Path(args.up_genes), args.top_n_genes)
    down_genes = read_genes(Path(args.down_genes), args.top_n_genes)
    db_versions = [item.strip() for item in args.db_version.split(",") if item.strip()]
    summary: List[Dict[str, Any]] = [
        {"metric": "query_url", "value": args.url, "notes": ""},
        {"metric": "db_versions", "value": ";".join(db_versions), "notes": ""},
        {"metric": "up_genes_submitted", "value": len(up_genes), "notes": str(Path(args.up_genes))},
        {"metric": "down_genes_submitted", "value": len(down_genes), "notes": str(Path(args.down_genes))},
        {"metric": "run_utc", "value": datetime.now(timezone.utc).isoformat(timespec="seconds"), "notes": ""},
        {"metric": "run_epoch_seconds", "value": int(time.time()), "notes": ""},
    ]

    records: List[Dict[str, Any]] = []
    matches: List[Dict[str, Any]] = []
    status_notes: List[str] = []

    for db_version in db_versions:
        payload = build_payload(up_genes, down_genes, db_version, args.tag)
        version_slug = safe_filename(db_version)
        (out_dir / f"l1000cds2_query_payload.{version_slug}.json").write_text(json.dumps(payload, indent=2, ensure_ascii=False), encoding="utf-8")

        status_code, body, error = post_json(args.url, payload, args.timeout)
        response_path = out_dir / f"l1000cds2_response.{version_slug}.json"
        if body:
            response_path.write_text(body, encoding="utf-8")
        else:
            response_path.write_text(json.dumps({"error": error}, indent=2), encoding="utf-8")

        version_status = "failed"
        parse_note = error
        parsed_count = 0
        if status_code and 200 <= status_code < 300 and body:
            try:
                response_obj = json.loads(body)
                top_meta = find_top_meta(response_obj)
                version_records = result_records(top_meta, db_version)
                records.extend(version_records)
                parsed_count = len(version_records)
                version_status = "success" if version_records else "no_parseable_top_results"
                parse_note = ""
            except json.JSONDecodeError as exc:
                parse_note = f"json_decode_error: {exc}"
                version_status = "json_decode_error"
        summary.append({"metric": f"http_status.{db_version}", "value": status_code, "notes": error})
        summary.append({"metric": f"query_status.{db_version}", "value": version_status, "notes": parse_note})
        summary.append({"metric": f"top_results_parsed.{db_version}", "value": parsed_count, "notes": str(response_path)})
        status_notes.append(f"{db_version}:{version_status}")

    by_pert, by_name, covered_candidates = load_covered_candidates(Path(args.coverage))
    matches = match_results(records, by_pert, by_name)
    candidate_status = candidate_status_rows(covered_candidates, matches)
    query_status = "success" if records else "failed"
    parse_note = "; ".join(status_notes)

    result_fields = [
        "db_version",
        "rank",
        "score",
        "pert_id",
        "pert_desc",
        "sig_id",
        "cell_id",
        "pert_dose",
        "pert_dose_unit",
        "pert_time",
        "pert_time_unit",
        "pubchem_id",
        "drugbank_id",
    ]
    match_fields = [
        *result_fields,
        "match_method",
        "priority_score",
        "ingredient_id",
        "ingredient_name",
        "ingredient_pubchem_id",
        "ingredient_inchikey",
        "candidate_lincs_source",
        "candidate_pert_id",
        "candidate_pert_iname",
    ]
    write_tsv(out_dir / "l1000cds2_top_results.tsv", records, result_fields)
    write_tsv(out_dir / "covered_candidate_reversal.tsv", matches, match_fields)
    candidate_status_fields = [
        "priority_score",
        "ingredient_id",
        "ingredient_name",
        "ingredient_pubchem_id",
        "ingredient_inchikey",
        "candidate_pert_id",
        "candidate_pert_iname",
        "candidate_lincs_sources",
        "hit_status",
        "best_db_version",
        "best_rank",
        "best_score",
        "best_sig_id",
        "best_cell_id",
        "best_pert_desc",
    ]
    write_tsv(out_dir / "candidate_reversal_status.tsv", candidate_status, candidate_status_fields)

    unique_candidate_ids = {row.get("ingredient_id") for row in matches if row.get("ingredient_id")}
    unique_candidate_pert_ids = {row.get("candidate_pert_id") for row in matches if row.get("candidate_pert_id")}
    summary.extend(
        [
            {"metric": "query_status", "value": query_status, "notes": parse_note},
            {"metric": "top_results_parsed", "value": len(records), "notes": str(out_dir / "l1000cds2_top_results.tsv")},
            {"metric": "covered_candidate_hit_rows", "value": len(matches), "notes": str(out_dir / "covered_candidate_reversal.tsv")},
            {"metric": "covered_candidate_tests", "value": len(candidate_status), "notes": str(out_dir / "candidate_reversal_status.tsv")},
            {"metric": "covered_candidate_ingredient_hits", "value": len(unique_candidate_ids), "notes": ""},
            {"metric": "covered_candidate_pert_id_hits", "value": len(unique_candidate_pert_ids), "notes": ""},
        ]
    )
    write_tsv(out_dir / "summary.tsv", summary, ["metric", "value", "notes"])
    write_report(Path(args.report), summary, matches)

    if query_status != "success" and args.fail_on_api_error:
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
