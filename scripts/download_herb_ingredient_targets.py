#!/usr/bin/env python3
"""Download HERB ingredient detail payloads and extract compound-target edges."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import sys
import time
import urllib.error
import urllib.request
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Dict, Iterable, List, Sequence, Set


API_URL = "http://47.92.70.12/chedi/api/?"
RAW_DIR = Path("data/raw/herb2/ingredient_details")
OUT_DIR = Path("results/m8_herb_ingredient_targets")
DOC_PATH = Path("docs/M8_HERB_INGREDIENT_TARGETS_RUN.md")
ACCESS_DATE = "2026-06-29"
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


def write_tsv(path: Path, rows: Iterable[Dict[str, object]], fields: Sequence[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(fields), delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def update_manifest(path: Path, audit_rows: List[Dict[str, object]]) -> None:
    existing: List[Dict[str, str]] = []
    if path.exists():
        with path.open(encoding="utf-8", newline="") as handle:
            existing = list(csv.DictReader(handle, delimiter="\t"))

    by_id = {row["dataset_id"]: row for row in existing if row.get("dataset_id")}
    for row in audit_rows:
        if row.get("status") == "failed":
            continue
        ingredient_id = str(row["ingredient_id"])
        by_id[f"herb2_detail_ingredient_{ingredient_id.lower()}"] = {
            "dataset_id": f"herb2_detail_ingredient_{ingredient_id.lower()}",
            "source": "HERB 2.0 detail API",
            "version_or_date": "V2",
            "file_path": str(row["raw_path"]),
            "source_url_or_api": API_URL,
            "access_date": ACCESS_DATE,
            "file_size_bytes": str(row["file_size_bytes"]),
            "md5": str(row["md5"]),
            "status": "analysis_ready",
            "notes": f"Ingredient detail API payload for {ingredient_id}; includes HERB detail relation fields such as drug_paper_target and ingredient_target",
        }

    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=MANIFEST_FIELDS, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for key in sorted(by_id):
            writer.writerow(by_id[key])


def md5_text(text: str) -> str:
    return hashlib.md5(text.encode("utf-8")).hexdigest()


def opener() -> urllib.request.OpenerDirector:
    return urllib.request.build_opener(urllib.request.ProxyHandler({}))


def fetch_detail(ingredient_id: str, timeout: int, reuse_existing: bool) -> Dict[str, object]:
    RAW_DIR.mkdir(parents=True, exist_ok=True)
    raw_path = RAW_DIR / f"{ingredient_id}.json"
    if reuse_existing and raw_path.exists() and raw_path.stat().st_size > 4:
        text = raw_path.read_text(encoding="utf-8")
        return {
            "ingredient_id": ingredient_id,
            "status": "cached",
            "raw_path": str(raw_path),
            "file_size_bytes": raw_path.stat().st_size,
            "md5": md5_text(text),
            "message": "",
        }

    payload = json.dumps(
        {
            "func_name": "detail_api",
            "label": "Ingredient",
            "v": ingredient_id,
            "key_id": ingredient_id,
        }
    ).encode("utf-8")
    request = urllib.request.Request(
        API_URL,
        data=payload,
        headers={
            "Content-Type": "application/json",
            "Accept": "application/json",
            "User-Agent": "zyy-herb-target-acquisition/1.0",
        },
        method="POST",
    )
    try:
        with opener().open(request, timeout=timeout) as response:
            text = response.read().decode("utf-8", errors="replace")
        parsed = json.loads(text)
        if not isinstance(parsed, dict):
            raise RuntimeError(f"unexpected payload type: {type(parsed).__name__}")
        raw_path.write_text(json.dumps(parsed, ensure_ascii=False, indent=2), encoding="utf-8")
        return {
            "ingredient_id": ingredient_id,
            "status": "downloaded",
            "raw_path": str(raw_path),
            "file_size_bytes": raw_path.stat().st_size,
            "md5": md5_text(raw_path.read_text(encoding="utf-8")),
            "message": "",
        }
    except (urllib.error.URLError, TimeoutError, json.JSONDecodeError, RuntimeError) as exc:
        return {
            "ingredient_id": ingredient_id,
            "status": "failed",
            "raw_path": str(raw_path),
            "file_size_bytes": "",
            "md5": "",
            "message": f"{type(exc).__name__}: {exc}",
        }


def link_title(value: object) -> str:
    if isinstance(value, dict):
        return str(value.get("title", ""))
    return str(value or "")


def link_href(value: object) -> str:
    if isinstance(value, dict):
        return str(value.get("link", ""))
    return ""


def load_gene_set(path: Path) -> Set[str]:
    if not path.exists():
        return set()
    return {line.strip().upper() for line in path.read_text(encoding="utf-8").splitlines() if line.strip()}


def extract_edges(ingredient: Dict[str, str], detail: Dict[str, object], up_genes: Set[str], down_genes: Set[str]) -> List[Dict[str, object]]:
    rows: List[Dict[str, object]] = []
    ingredient_id = ingredient["ingredient_id"]
    ingredient_name = ingredient["ingredient_name"]

    for row in detail.get("drug_paper_target", [])[1:]:
        if not isinstance(row, list) or len(row) < 4:
            continue
        target_id = link_title(row[0])
        gene_symbol = str(row[1] or "").upper()
        protein_name = str(row[2] or "")
        references = row[3] if isinstance(row[3], list) else []
        ref_rows = references[1:] if references else [[]]
        for ref in ref_rows:
            ref = ref if isinstance(ref, list) else []
            rows.append(
                {
                    "ingredient_id": ingredient_id,
                    "ingredient_name": ingredient_name,
                    "target_id": target_id,
                    "gene_symbol": gene_symbol,
                    "protein_name": protein_name,
                    "edge_source": "HERB_detail.drug_paper_target",
                    "evidence_level": "literature_mined_pubmed_sentence",
                    "target_source": "",
                    "reference_id": link_title(ref[0]) if len(ref) > 0 else "",
                    "reference_url": link_href(ref[0]) if len(ref) > 0 else "",
                    "pubmed_id": link_title(ref[1]) if len(ref) > 1 else "",
                    "pubmed_url": link_href(ref[1]) if len(ref) > 1 else "",
                    "reference_title": ref[2] if len(ref) > 2 else "",
                    "relationship": ref[3] if len(ref) > 3 else "",
                    "grade": ref[4] if len(ref) > 4 else "",
                    "supporting_sentence": ref[5] if len(ref) > 5 else "",
                    "in_bulk_uc_ibd_up150": "yes" if gene_symbol in up_genes else "no",
                    "in_bulk_uc_ibd_down150": "yes" if gene_symbol in down_genes else "no",
                    "candidate_priority_class": ingredient.get("priority_class", ""),
                    "candidate_priority_score": ingredient.get("priority_score", ""),
                }
            )

    for row in detail.get("ingredient_target", [])[1:]:
        if not isinstance(row, list) or len(row) < 6:
            continue
        gene_symbol = str(row[1] or "").upper()
        rows.append(
            {
                "ingredient_id": ingredient_id,
                "ingredient_name": ingredient_name,
                "target_id": link_title(row[0]),
                "gene_symbol": gene_symbol,
                "protein_name": row[3],
                "edge_source": "HERB_detail.ingredient_target",
                "evidence_level": "database_integrated_target",
                "target_source": row[5],
                "reference_id": "",
                "reference_url": "",
                "pubmed_id": "",
                "pubmed_url": "",
                "reference_title": "",
                "relationship": "",
                "grade": "",
                "supporting_sentence": "",
                "in_bulk_uc_ibd_up150": "yes" if gene_symbol in up_genes else "no",
                "in_bulk_uc_ibd_down150": "yes" if gene_symbol in down_genes else "no",
                "candidate_priority_class": ingredient.get("priority_class", ""),
                "candidate_priority_score": ingredient.get("priority_score", ""),
            }
        )
    return rows


def summarize_edges(ingredients: List[Dict[str, str]], edges: List[Dict[str, object]]) -> List[Dict[str, object]]:
    grouped: Dict[str, List[Dict[str, object]]] = {}
    for row in edges:
        grouped.setdefault(str(row["ingredient_id"]), []).append(row)

    summary_rows: List[Dict[str, object]] = []
    for ingredient in ingredients:
        rows = grouped.get(ingredient["ingredient_id"], [])
        literature = [row for row in rows if row["edge_source"] == "HERB_detail.drug_paper_target"]
        database = [row for row in rows if row["edge_source"] == "HERB_detail.ingredient_target"]
        literature_genes = sorted({str(row["gene_symbol"]) for row in literature if row["gene_symbol"]})
        database_genes = sorted({str(row["gene_symbol"]) for row in database if row["gene_symbol"]})
        all_genes = sorted(set(literature_genes) | set(database_genes))
        up_hits = sorted({str(row["gene_symbol"]) for row in rows if row["in_bulk_uc_ibd_up150"] == "yes"})
        down_hits = sorted({str(row["gene_symbol"]) for row in rows if row["in_bulk_uc_ibd_down150"] == "yes"})
        summary_rows.append(
            {
                "ingredient_id": ingredient["ingredient_id"],
                "ingredient_name": ingredient["ingredient_name"],
                "priority_class": ingredient.get("priority_class", ""),
                "priority_score": ingredient.get("priority_score", ""),
                "literature_edge_rows": len(literature),
                "literature_unique_targets": len(literature_genes),
                "database_edge_rows": len(database),
                "database_unique_targets": len(database_genes),
                "all_unique_targets": len(all_genes),
                "bulk_uc_ibd_up150_target_hits": len(up_hits),
                "bulk_uc_ibd_down150_target_hits": len(down_hits),
                "up150_hit_genes": ";".join(up_hits),
                "down150_hit_genes": ";".join(down_hits),
                "literature_target_genes": ";".join(literature_genes),
                "database_target_genes": ";".join(database_genes),
            }
        )
    summary_rows.sort(
        key=lambda row: (
            -int(row["bulk_uc_ibd_up150_target_hits"]),
            -int(row["literature_unique_targets"]),
            -float(row["priority_score"] or 0),
        )
    )
    return summary_rows


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--candidate-table", default="results/m2_herb2_uc_ibd_evidence/prioritized_ingredient_candidates.tsv")
    parser.add_argument("--priority-class", default="priority_for_lincs_mapping")
    parser.add_argument("--limit", type=int, default=30)
    parser.add_argument("--workers", type=int, default=2)
    parser.add_argument("--timeout", type=int, default=120)
    parser.add_argument("--manifest", default="data/manifest.tsv")
    parser.add_argument("--reuse-existing", action="store_true")
    args = parser.parse_args()

    candidates = read_tsv(Path(args.candidate_table))
    ingredients = [row for row in candidates if not args.priority_class or row.get("priority_class") == args.priority_class]
    if args.limit:
        ingredients = ingredients[: args.limit]

    up_genes = load_gene_set(Path("results/m4_geo_uc_ibd_signature/lincs_query_up_genes.txt"))
    down_genes = load_gene_set(Path("results/m4_geo_uc_ibd_signature/lincs_query_down_genes.txt"))

    audit_rows: List[Dict[str, object]] = []
    print(f"[herb-target] requesting {len(ingredients)} ingredient detail payloads", file=sys.stderr)
    with ThreadPoolExecutor(max_workers=max(1, args.workers)) as pool:
        futures = {
            pool.submit(fetch_detail, row["ingredient_id"], args.timeout, args.reuse_existing): row["ingredient_id"]
            for row in ingredients
        }
        for future in as_completed(futures):
            audit = future.result()
            audit_rows.append(audit)
            print(f"[herb-target] {audit['ingredient_id']} {audit['status']} {audit.get('message', '')}", file=sys.stderr)
            time.sleep(0.2)

    by_id = {row["ingredient_id"]: row for row in ingredients}
    edge_rows: List[Dict[str, object]] = []
    for audit in audit_rows:
        if audit["status"] == "failed":
            continue
        ingredient_id = str(audit["ingredient_id"])
        raw_path = Path(str(audit["raw_path"]))
        detail = json.loads(raw_path.read_text(encoding="utf-8"))
        edge_rows.extend(extract_edges(by_id[ingredient_id], detail, up_genes, down_genes))

    edge_fields = [
        "ingredient_id",
        "ingredient_name",
        "target_id",
        "gene_symbol",
        "protein_name",
        "edge_source",
        "evidence_level",
        "target_source",
        "reference_id",
        "reference_url",
        "pubmed_id",
        "pubmed_url",
        "reference_title",
        "relationship",
        "grade",
        "supporting_sentence",
        "in_bulk_uc_ibd_up150",
        "in_bulk_uc_ibd_down150",
        "candidate_priority_class",
        "candidate_priority_score",
    ]
    summary_fields = [
        "ingredient_id",
        "ingredient_name",
        "priority_class",
        "priority_score",
        "literature_edge_rows",
        "literature_unique_targets",
        "database_edge_rows",
        "database_unique_targets",
        "all_unique_targets",
        "bulk_uc_ibd_up150_target_hits",
        "bulk_uc_ibd_down150_target_hits",
        "up150_hit_genes",
        "down150_hit_genes",
        "literature_target_genes",
        "database_target_genes",
    ]
    audit_fields = ["ingredient_id", "status", "raw_path", "file_size_bytes", "md5", "message"]
    summary_rows = summarize_edges(ingredients, edge_rows)

    write_tsv(OUT_DIR / "herb_ingredient_target_edges.tsv", edge_rows, edge_fields)
    write_tsv(OUT_DIR / "herb_ingredient_target_summary.tsv", summary_rows, summary_fields)
    write_tsv(OUT_DIR / "herb_detail_download_audit.tsv", audit_rows, audit_fields)
    update_manifest(Path(args.manifest), audit_rows)

    successful = [row for row in audit_rows if row["status"] != "failed"]
    failed = [row for row in audit_rows if row["status"] == "failed"]
    literature_edges = [row for row in edge_rows if row["edge_source"] == "HERB_detail.drug_paper_target"]
    database_edges = [row for row in edge_rows if row["edge_source"] == "HERB_detail.ingredient_target"]
    top_rows = summary_rows[:12]
    report_lines = [
        "# M8 HERB Ingredient-Target Acquisition Run",
        "",
        "## Scope",
        "",
        "This run downloads HERB ingredient detail payloads through the public JSON detail API and extracts compound-target relationships. It separates literature-mined PubMed sentence-supported targets from database-integrated ingredient targets.",
        "",
        "## API",
        "",
        f"- Endpoint: `{API_URL}`",
        "- Method: POST JSON body with `func_name=detail_api`, `label=Ingredient`, `v=<Ingredient_id>`, and `key_id=<Ingredient_id>`.",
        "",
        "## Outputs",
        "",
        "- `data/raw/herb2/ingredient_details/*.json`",
        "- `results/m8_herb_ingredient_targets/herb_ingredient_target_edges.tsv`",
        "- `results/m8_herb_ingredient_targets/herb_ingredient_target_summary.tsv`",
        "- `results/m8_herb_ingredient_targets/herb_detail_download_audit.tsv`",
        "",
        "## Acquisition Summary",
        "",
        f"- Requested ingredients: {len(ingredients)}",
        f"- Successful detail payloads: {len(successful)}",
        f"- Failed detail payloads: {len(failed)}",
        f"- Literature/PubMed target edge rows: {len(literature_edges)}",
        f"- Database-integrated target edge rows: {len(database_edges)}",
        "",
        "## Top Candidates By Disease-Module Target Overlap",
        "",
        "| Ingredient | Literature targets | Database targets | UC/IBD up150 target hits | UC/IBD down150 target hits | Up-hit genes |",
        "|---|---:|---:|---:|---:|---|",
    ]
    for row in top_rows:
        report_lines.append(
            f"| {row['ingredient_name']} | {row['literature_unique_targets']} | {row['database_unique_targets']} | {row['bulk_uc_ibd_up150_target_hits']} | {row['bulk_uc_ibd_down150_target_hits']} | {row['up150_hit_genes']} |"
        )
    report_lines.extend(
        [
            "",
            "## Interpretation Boundary",
            "",
            "- `drug_paper_target` is the preferred HERB target layer because it carries PubMed IDs, relationship labels, grades, and supporting sentences.",
            "- `ingredient_target` is retained as a secondary annotation layer because it is database-integrated and may include predicted or inherited targets.",
            "- Target overlap with the bulk UC/IBD module is hypothesis-generating. It should be combined with scRNA localization, GWAS/MR or independent disease-gene evidence before mechanism claims.",
            "",
        ]
    )
    DOC_PATH.write_text("\n".join(report_lines), encoding="utf-8")
    return 1 if failed and not successful else 0


if __name__ == "__main__":
    raise SystemExit(main())
