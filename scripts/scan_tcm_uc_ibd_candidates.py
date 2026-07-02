#!/usr/bin/env python3
"""Scan local HERB 2.0 smoke-test tables for UC/IBD-relevant TCM entities."""

from __future__ import annotations

import argparse
import csv
import re
from collections import Counter
from pathlib import Path
from typing import Dict, Iterable, List, Tuple


DIRECT_TERMS = [
    "ulcerative colitis",
    "inflammatory bowel disease",
    "crohn",
    "colitis",
    "enteritis",
    "proctitis",
    "溃疡性结肠炎",
    "炎症性肠病",
    "克罗恩",
    "结肠炎",
    "肠炎",
    "直肠炎",
]

PROXY_TERMS = [
    "diarrhea",
    "dysentery",
    "bloody stool",
    "hematochezia",
    "abdominal pain",
    "leukorrheal",
    "泄泻",
    "腹泻",
    "痢疾",
    "便血",
    "下痢",
    "腹痛",
    "肠",
]

FORMULA_SCAN_FIELDS = [
    "Formula_id",
    "Formula_pinyin_name",
    "Formula_cn_name",
    "Formula_en_name",
    "Herbs_in_Chinese",
    "Herbs_in_pinyin",
    "Syndromes_in_Chinese",
    "Syndromes_in_English",
    "Indications_in_Chinese",
    "Indications_in_English",
    "Source",
    "ETCM_id",
]

HERB_SCAN_FIELDS = [
    "Herb_id",
    "Herb_pinyin_name",
    "Herb_cn_name",
    "Herb_en_name",
    "Herb_latin_name",
    "Function",
    "Indication",
    "Therapeutic_en_class",
    "Therapeutic_cn_class",
    "SymMap_id",
    "TCMID_id",
    "TCMSP_id",
    "TCM_ID_id",
]


def normalize(text: str) -> str:
    return " ".join((text or "").lower().split())


def match_terms(text: str) -> Tuple[List[str], List[str], float]:
    haystack = normalize(text)
    direct = [term for term in DIRECT_TERMS if normalize(term) in haystack]
    proxy = [term for term in PROXY_TERMS if normalize(term) in haystack]
    score = len(set(direct)) * 3.0 + len(set(proxy)) * 1.0
    return sorted(set(direct)), sorted(set(proxy)), score


def read_tsv(path: Path) -> Iterable[Dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        yield from csv.DictReader(handle, delimiter="\t")


def scan_table(path: Path, entity_type: str, fields: List[str], id_field: str, name_fields: List[str]) -> List[Dict[str, str]]:
    rows: List[Dict[str, str]] = []
    for row in read_tsv(path):
        scan_text = " | ".join(row.get(field, "") for field in fields)
        direct, proxy, score = match_terms(scan_text)
        if not direct and not proxy:
            continue
        evidence_class = "direct_disease_term" if direct else "symptom_proxy_only"
        entity_name = " / ".join(row.get(field, "") for field in name_fields if row.get(field, ""))
        rows.append(
            {
                "source": "HERB 2.0",
                "entity_type": entity_type,
                "entity_id": row.get(id_field, ""),
                "entity_name": entity_name,
                "evidence_class": evidence_class,
                "direct_terms": ";".join(direct),
                "proxy_terms": ";".join(proxy),
                "match_score": f"{score:.1f}",
                "indications_or_function": " | ".join(
                    row.get(field, "")
                    for field in [
                        "Indications_in_English",
                        "Indications_in_Chinese",
                        "Indication",
                        "Function",
                        "Syndromes_in_English",
                        "Syndromes_in_Chinese",
                    ]
                    if row.get(field, "")
                ),
                "composition_or_taxonomy": " | ".join(
                    row.get(field, "")
                    for field in ["Herbs_in_pinyin", "Herbs_in_Chinese", "Herb_latin_name"]
                    if row.get(field, "")
                ),
                "cross_database_ids": " | ".join(
                    f"{field}={row.get(field, '')}"
                    for field in ["ETCM_id", "SymMap_id", "TCMID_id", "TCMSP_id", "TCM_ID_id"]
                    if row.get(field, "") and row.get(field, "") != "NA"
                ),
                "next_step": "manual_indication_review_then_ingredient_mapping" if direct else "keep_as_negative_or_secondary_unless_supported_by_direct_evidence",
            }
        )
    rows.sort(key=lambda item: (-float(item["match_score"]), item["evidence_class"], item["entity_id"]))
    return rows


def write_tsv(path: Path, rows: List[Dict[str, str]], fieldnames: List[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def summary_rows(formula_rows: List[Dict[str, str]], herb_rows: List[Dict[str, str]]) -> List[Dict[str, str]]:
    rows: List[Dict[str, str]] = []
    for entity_type, candidates in [("formula", formula_rows), ("herb", herb_rows)]:
        by_class = Counter(row["evidence_class"] for row in candidates)
        rows.append(
            {
                "source": "HERB 2.0",
                "entity_type": entity_type,
                "total_matched": str(len(candidates)),
                "direct_disease_term": str(by_class.get("direct_disease_term", 0)),
                "symptom_proxy_only": str(by_class.get("symptom_proxy_only", 0)),
                "top_direct_candidates": "; ".join(
                    f"{row['entity_id']}:{row['entity_name']}"
                    for row in candidates
                    if row["evidence_class"] == "direct_disease_term"
                )[:1000],
                "notes": "Direct terms are candidates for manual review; symptom-only matches are too broad for primary claims.",
            }
        )
    return rows


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--herb-formula", default="data/raw/test_downloads/herb2/HERB_formula_info_v2.txt")
    parser.add_argument("--herb-herb", default="data/raw/test_downloads/herb2/HERB_herb_info_v2.txt")
    parser.add_argument("--out-dir", default="results/m1_tcm_candidate_scan")
    args = parser.parse_args()

    out_dir = Path(args.out_dir)
    formula_rows = scan_table(
        Path(args.herb_formula),
        "formula",
        FORMULA_SCAN_FIELDS,
        "Formula_id",
        ["Formula_pinyin_name", "Formula_cn_name", "Formula_en_name"],
    )
    herb_rows = scan_table(
        Path(args.herb_herb),
        "herb",
        HERB_SCAN_FIELDS,
        "Herb_id",
        ["Herb_pinyin_name", "Herb_cn_name", "Herb_en_name", "Herb_latin_name"],
    )

    fieldnames = [
        "source",
        "entity_type",
        "entity_id",
        "entity_name",
        "evidence_class",
        "direct_terms",
        "proxy_terms",
        "match_score",
        "indications_or_function",
        "composition_or_taxonomy",
        "cross_database_ids",
        "next_step",
    ]
    write_tsv(out_dir / "herb2_formula_uc_ibd_candidates.tsv", formula_rows, fieldnames)
    write_tsv(out_dir / "herb2_herb_uc_ibd_candidates.tsv", herb_rows, fieldnames)
    write_tsv(
        out_dir / "candidate_scan_summary.tsv",
        summary_rows(formula_rows, herb_rows),
        ["source", "entity_type", "total_matched", "direct_disease_term", "symptom_proxy_only", "top_direct_candidates", "notes"],
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
