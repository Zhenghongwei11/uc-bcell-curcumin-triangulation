#!/usr/bin/env python3
"""Query Open Targets disease-target evidence for M9 HERB target genes."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import time
from collections import defaultdict
from pathlib import Path
from typing import Any, Dict, Iterable, List, Sequence, Set, Tuple

import requests


ENDPOINT = "https://api.platform.opentargets.org/api/v4/graphql"
ACCESS_DATE = "2026-06-29"
RAW_DIR = Path("data/raw/open_targets")
OUT_DIR = Path("results/m10_open_targets")
DOC_PATH = Path("docs/M10_OPEN_TARGETS_TARGET_PLAUSIBILITY_RUN.md")
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

DISEASES = [
    ("MONDO_0005101", "ulcerative_colitis"),
    ("MONDO_0005265", "inflammatory_bowel_disease"),
    ("MONDO_0005011", "crohn_disease"),
]

DATATYPE_COLUMNS = [
    "genetic_association",
    "clinical",
    "literature",
    "rna_expression",
    "animal_model",
    "genetic_literature",
    "somatic_mutation",
    "known_drug",
    "affected_pathway",
]

QUERY = """
query DiseaseTargets($efoId:String!, $index:Int!, $size:Int!) {
  disease(efoId:$efoId) {
    id
    name
    associatedTargets(page:{index:$index,size:$size}) {
      count
      rows {
        target {
          id
          approvedSymbol
          approvedName
        }
        score
        datatypeScores {
          id
          score
        }
        datasourceScores {
          id
          score
        }
      }
    }
  }
}
"""


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


def compact_scores(scores: Iterable[Dict[str, Any]]) -> Dict[str, float]:
    out: Dict[str, float] = {}
    for item in scores or []:
        if item.get("id"):
            out[str(item["id"])] = float(item.get("score") or 0.0)
    return out


def format_scores(scores: Dict[str, float]) -> str:
    return ";".join(f"{key}:{value:.6g}" for key, value in sorted(scores.items()))


def fetch_disease_targets(
    disease_id: str,
    disease_slug: str,
    page_size: int,
    timeout: int,
    sleep_seconds: float,
    reuse_existing: bool,
) -> Tuple[Dict[str, Any], Path]:
    RAW_DIR.mkdir(parents=True, exist_ok=True)
    raw_path = RAW_DIR / f"open_targets_{disease_slug}_{disease_id}.json"
    if reuse_existing and raw_path.exists() and raw_path.stat().st_size > 100:
        return json.loads(raw_path.read_text(encoding="utf-8")), raw_path

    session = requests.Session()
    all_rows: List[Dict[str, Any]] = []
    disease_name = disease_slug
    total_count = None
    page_index = 0
    while True:
        payload = {
            "query": QUERY,
            "variables": {"efoId": disease_id, "index": page_index, "size": page_size},
        }
        response = session.post(
            ENDPOINT,
            json=payload,
            timeout=timeout,
            headers={"User-Agent": "zyy-open-targets-acquisition/1.0"},
        )
        response.raise_for_status()
        parsed = response.json()
        if parsed.get("errors"):
            raise RuntimeError(json.dumps(parsed["errors"], ensure_ascii=False))
        disease = (parsed.get("data") or {}).get("disease")
        if not disease:
            raise RuntimeError(f"Open Targets returned no disease for {disease_id}")
        associated = disease.get("associatedTargets") or {}
        total_count = int(associated.get("count") or 0)
        disease_name = disease.get("name") or disease_slug
        rows = associated.get("rows") or []
        all_rows.extend(rows)
        if len(all_rows) >= total_count or not rows:
            break
        page_index += 1
        if sleep_seconds:
            time.sleep(sleep_seconds)

    payload = {
        "source": "Open Targets Platform GraphQL API",
        "endpoint": ENDPOINT,
        "access_date": ACCESS_DATE,
        "disease_id": disease_id,
        "disease_name": disease_name,
        "associated_target_count": total_count,
        "rows": all_rows,
    }
    raw_path.write_text(json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8")
    return payload, raw_path


def update_manifest(raw_paths: Sequence[Tuple[str, str, Path]]) -> None:
    existing: List[Dict[str, str]] = []
    if MANIFEST_PATH.exists():
        with MANIFEST_PATH.open(encoding="utf-8", newline="") as handle:
            existing = list(csv.DictReader(handle, delimiter="\t"))

    by_id = {row["dataset_id"]: row for row in existing if row.get("dataset_id")}
    for disease_id, disease_slug, path in raw_paths:
        by_id[f"open_targets_{disease_slug}_{disease_id.lower()}"] = {
            "dataset_id": f"open_targets_{disease_slug}_{disease_id.lower()}",
            "source": "Open Targets Platform GraphQL API",
            "version_or_date": ACCESS_DATE,
            "file_path": str(path),
            "source_url_or_api": ENDPOINT,
            "access_date": ACCESS_DATE,
            "file_size_bytes": str(path.stat().st_size),
            "md5": md5_file(path),
            "status": "analysis_ready",
            "notes": f"Associated target evidence for {disease_slug}; cached from disease(efoId).associatedTargets",
        }

    with MANIFEST_PATH.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=MANIFEST_FIELDS, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for key in sorted(by_id):
            writer.writerow(by_id[key])


def load_m9_context(edge_path: Path, curcumin_path: Path) -> Tuple[Set[str], Dict[str, Set[str]], Dict[str, Dict[str, str]], Dict[str, str]]:
    m9_edges = read_tsv(edge_path)
    target_genes = {row["gene_symbol"].upper() for row in m9_edges if row.get("gene_symbol")}
    gene_to_candidates: Dict[str, Set[str]] = defaultdict(set)
    candidates: Dict[str, Dict[str, str]] = {}
    for row in m9_edges:
        gene = row["gene_symbol"].upper()
        ingredient_id = row["ingredient_id"]
        ingredient_name = row["ingredient_name"]
        gene_to_candidates[gene].add(ingredient_name)
        candidates[ingredient_id] = {
            "ingredient_id": ingredient_id,
            "ingredient_name": ingredient_name,
            "candidate_priority_class": row.get("candidate_priority_class", ""),
            "candidate_priority_score": row.get("candidate_priority_score", ""),
        }

    curcumin_tiers: Dict[str, str] = {}
    if curcumin_path.exists():
        for row in read_tsv(curcumin_path):
            curcumin_tiers[row["gene_symbol"].upper()] = row.get("evidence_tier", "")
    return target_genes, gene_to_candidates, candidates, curcumin_tiers


def flatten_filtered_rows(
    payloads: Sequence[Dict[str, Any]],
    target_genes: Set[str],
    gene_to_candidates: Dict[str, Set[str]],
    curcumin_tiers: Dict[str, str],
) -> List[Dict[str, Any]]:
    rows: List[Dict[str, Any]] = []
    for payload in payloads:
        for item in payload.get("rows", []):
            target = item.get("target") or {}
            gene = str(target.get("approvedSymbol") or "").upper()
            if gene not in target_genes:
                continue
            datatype_scores = compact_scores(item.get("datatypeScores") or [])
            datasource_scores = compact_scores(item.get("datasourceScores") or [])
            row: Dict[str, Any] = {
                "disease_id": payload["disease_id"],
                "disease_name": payload["disease_name"],
                "target_ensembl_id": target.get("id", ""),
                "gene_symbol": gene,
                "target_approved_name": target.get("approvedName", ""),
                "overall_score": f"{float(item.get('score') or 0.0):.6g}",
                "candidate_ingredients": ";".join(sorted(gene_to_candidates.get(gene, set()))),
                "is_curcumin_target": "yes" if gene in curcumin_tiers else "no",
                "curcumin_evidence_tier": curcumin_tiers.get(gene, ""),
                "all_datatype_scores": format_scores(datatype_scores),
                "all_datasource_scores": format_scores(datasource_scores),
            }
            for column in DATATYPE_COLUMNS:
                row[f"{column}_score"] = f"{datatype_scores.get(column, 0.0):.6g}"
            rows.append(row)
    rows.sort(key=lambda row: (row["gene_symbol"], row["disease_name"], -float(row["overall_score"])))
    return rows


def summarize_candidates(
    edge_path: Path,
    evidence_rows: Sequence[Dict[str, Any]],
) -> List[Dict[str, Any]]:
    by_gene: Dict[str, List[Dict[str, Any]]] = defaultdict(list)
    for row in evidence_rows:
        by_gene[row["gene_symbol"]].append(row)

    by_candidate: Dict[str, Dict[str, Any]] = {}
    for edge in read_tsv(edge_path):
        ingredient_id = edge["ingredient_id"]
        gene = edge["gene_symbol"].upper()
        entry = by_candidate.setdefault(
            ingredient_id,
            {
                "ingredient_id": ingredient_id,
                "ingredient_name": edge["ingredient_name"],
                "candidate_priority_class": edge.get("candidate_priority_class", ""),
                "candidate_priority_score": edge.get("candidate_priority_score", ""),
                "disease_context_targets": set(),
                "open_targets_supported_genes": set(),
                "genetic_supported_genes": set(),
                "clinical_supported_genes": set(),
                "literature_supported_genes": set(),
                "max_overall_score": 0.0,
                "max_genetic_association_score": 0.0,
            },
        )
        entry["disease_context_targets"].add(gene)
        for evidence in by_gene.get(gene, []):
            overall = float(evidence["overall_score"])
            genetic = float(evidence["genetic_association_score"])
            clinical = float(evidence["clinical_score"])
            literature = float(evidence["literature_score"])
            if overall > 0:
                entry["open_targets_supported_genes"].add(gene)
            if genetic > 0:
                entry["genetic_supported_genes"].add(gene)
            if clinical > 0:
                entry["clinical_supported_genes"].add(gene)
            if literature > 0:
                entry["literature_supported_genes"].add(gene)
            entry["max_overall_score"] = max(entry["max_overall_score"], overall)
            entry["max_genetic_association_score"] = max(entry["max_genetic_association_score"], genetic)

    rows: List[Dict[str, Any]] = []
    for entry in by_candidate.values():
        disease_targets = sorted(entry["disease_context_targets"])
        supported = sorted(entry["open_targets_supported_genes"])
        genetic = sorted(entry["genetic_supported_genes"])
        clinical = sorted(entry["clinical_supported_genes"])
        literature = sorted(entry["literature_supported_genes"])
        rows.append(
            {
                "ingredient_id": entry["ingredient_id"],
                "ingredient_name": entry["ingredient_name"],
                "candidate_priority_class": entry["candidate_priority_class"],
                "candidate_priority_score": entry["candidate_priority_score"],
                "n_disease_context_targets": len(disease_targets),
                "n_open_targets_supported_genes": len(supported),
                "open_targets_supported_genes": ";".join(supported),
                "n_genetic_supported_genes": len(genetic),
                "genetic_supported_genes": ";".join(genetic),
                "n_clinical_supported_genes": len(clinical),
                "clinical_supported_genes": ";".join(clinical),
                "n_literature_supported_genes": len(literature),
                "literature_supported_genes": ";".join(literature),
                "max_overall_score": f"{entry['max_overall_score']:.6g}",
                "max_genetic_association_score": f"{entry['max_genetic_association_score']:.6g}",
            }
        )
    rows.sort(
        key=lambda row: (
            -int(row["n_genetic_supported_genes"]),
            -int(row["n_open_targets_supported_genes"]),
            -float(row["max_overall_score"]),
            row["ingredient_name"].lower(),
        )
    )
    return rows


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--m9-edge-table", default="results/m9_target_cell_bridge/disease_context_herb_target_edges.tsv")
    parser.add_argument("--curcumin-table", default="results/m9_target_cell_bridge/curcumin_bcell_target_evidence.tsv")
    parser.add_argument("--page-size", type=int, default=1000)
    parser.add_argument("--timeout", type=int, default=60)
    parser.add_argument("--sleep-seconds", type=float, default=0.2)
    parser.add_argument("--no-reuse-existing", action="store_true")
    args = parser.parse_args()

    target_genes, gene_to_candidates, _, curcumin_tiers = load_m9_context(Path(args.m9_edge_table), Path(args.curcumin_table))
    payloads: List[Dict[str, Any]] = []
    raw_paths: List[Tuple[str, str, Path]] = []
    for disease_id, disease_slug in DISEASES:
        payload, raw_path = fetch_disease_targets(
            disease_id=disease_id,
            disease_slug=disease_slug,
            page_size=args.page_size,
            timeout=args.timeout,
            sleep_seconds=args.sleep_seconds,
            reuse_existing=not args.no_reuse_existing,
        )
        payloads.append(payload)
        raw_paths.append((disease_id, disease_slug, raw_path))
    update_manifest(raw_paths)

    evidence_rows = flatten_filtered_rows(payloads, target_genes, gene_to_candidates, curcumin_tiers)
    summary_rows = summarize_candidates(Path(args.m9_edge_table), evidence_rows)

    evidence_fields = [
        "disease_id",
        "disease_name",
        "target_ensembl_id",
        "gene_symbol",
        "target_approved_name",
        "overall_score",
        *[f"{column}_score" for column in DATATYPE_COLUMNS],
        "candidate_ingredients",
        "is_curcumin_target",
        "curcumin_evidence_tier",
        "all_datatype_scores",
        "all_datasource_scores",
    ]
    summary_fields = [
        "ingredient_id",
        "ingredient_name",
        "candidate_priority_class",
        "candidate_priority_score",
        "n_disease_context_targets",
        "n_open_targets_supported_genes",
        "open_targets_supported_genes",
        "n_genetic_supported_genes",
        "genetic_supported_genes",
        "n_clinical_supported_genes",
        "clinical_supported_genes",
        "n_literature_supported_genes",
        "literature_supported_genes",
        "max_overall_score",
        "max_genetic_association_score",
    ]

    write_tsv(OUT_DIR / "open_targets_ibd_target_evidence.tsv", evidence_rows, evidence_fields)
    write_tsv(OUT_DIR / "candidate_open_targets_summary.tsv", summary_rows, summary_fields)

    curcumin_summary = next((row for row in summary_rows if row["ingredient_name"].lower() == "curcumin"), None)
    top_rows = summary_rows[:12]
    report_lines = [
        "# M10 Open Targets Target Plausibility Run",
        "",
        "## Scope",
        "",
        "This run adds an external disease-target plausibility layer from the Open Targets Platform GraphQL API for UC, IBD, and Crohn disease. It filters Open Targets associated targets to the M9 HERB disease-context target genes.",
        "",
        "## Key Outputs",
        "",
        "- `data/raw/open_targets/open_targets_ulcerative_colitis_MONDO_0005101.json`",
        "- `data/raw/open_targets/open_targets_inflammatory_bowel_disease_MONDO_0005265.json`",
        "- `data/raw/open_targets/open_targets_crohn_disease_MONDO_0005011.json`",
        "- `results/m10_open_targets/open_targets_ibd_target_evidence.tsv`",
        "- `results/m10_open_targets/candidate_open_targets_summary.tsv`",
        "",
        "## Acquisition Summary",
        "",
    ]
    for payload in payloads:
        report_lines.append(
            f"- {payload['disease_name']} ({payload['disease_id']}): {payload['associated_target_count']} associated targets cached; {len(payload.get('rows', []))} rows retrieved."
        )
    report_lines.extend(
        [
            f"- M9 disease-context target genes queried against Open Targets: {len(target_genes)}.",
            f"- Filtered Open Targets disease-target evidence rows retained: {len(evidence_rows)}.",
            "",
            "## Top Candidate External Disease-Target Support",
            "",
            "| Ingredient | Disease-context targets | OT-supported | Genetic-supported | Clinical-supported | Literature-supported | Max OT score | Max genetic score |",
            "|---|---:|---:|---:|---:|---:|---:|---:|",
        ]
    )
    for row in top_rows:
        report_lines.append(
            f"| {row['ingredient_name']} | {row['n_disease_context_targets']} | {row['n_open_targets_supported_genes']} | {row['n_genetic_supported_genes']} | {row['n_clinical_supported_genes']} | {row['n_literature_supported_genes']} | {row['max_overall_score']} | {row['max_genetic_association_score']} |"
        )
    report_lines.extend(["", "## Curcumin Focus", ""])
    if curcumin_summary:
        report_lines.extend(
            [
                f"- Open Targets-supported curcumin targets: {curcumin_summary['n_open_targets_supported_genes']} ({curcumin_summary['open_targets_supported_genes']}).",
                f"- Genetic-supported curcumin targets: {curcumin_summary['n_genetic_supported_genes']} ({curcumin_summary['genetic_supported_genes']}).",
                f"- Clinical-supported curcumin targets: {curcumin_summary['n_clinical_supported_genes']} ({curcumin_summary['clinical_supported_genes']}).",
                f"- Max Open Targets score among curcumin targets: {curcumin_summary['max_overall_score']}; max genetic association score: {curcumin_summary['max_genetic_association_score']}.",
            ]
        )
    report_lines.extend(
        [
            "",
            "## Interpretation Boundary",
            "",
            "- Open Targets supports disease-target plausibility, not compound efficacy.",
            "- Genetic or clinical target support strengthens target relevance, but it does not validate that a TCM ingredient modulates that target in patients.",
            "- Scores are used as orthogonal prioritization evidence and should be reported with source/date rather than treated as experimental results.",
            "",
        ]
    )
    DOC_PATH.write_text("\n".join(report_lines), encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
