#!/usr/bin/env python3
"""Prioritize HERB UC/IBD ingredient evidence for downstream LINCS mapping."""

from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path
from typing import Dict, Iterable, List


LOW_TCM_SPECIFICITY = re.compile(
    r"\b(vitamin|progesterone|hydrocortisone|caffeine|calcium|starch|taxol|paclitaxel|indomethacin|caproate|"
    r"cannabidiol|cortisol|cyanocobalamin|tocopherol)\b",
    re.I,
)

HIGH_VALUE_TERMS = ["ulcerative colitis", "inflammatory bowel disease", "colitis", "crohn", "enteritis"]


def read_tsv(path: Path) -> Iterable[Dict[str, str]]:
    with path.open(encoding="utf-8", errors="replace", newline="") as handle:
        yield from csv.DictReader(handle, delimiter="\t")


def priority_score(row: Dict[str, str]) -> float:
    reference_count = int(row.get("reference_count") or 0)
    clinical_count = int(row.get("clinical_trial_count") or 0)
    meta_count = int(row.get("meta_analysis_count") or 0)
    terms = row.get("matched_terms", "")
    score = reference_count * 4.0 + clinical_count * 5.0 + min(meta_count, 3) * 1.0
    if "ulcerative colitis" in terms:
        score += 4.0
    if "inflammatory bowel disease" in terms:
        score += 3.0
    if "crohn" in terms:
        score += 2.0
    if row.get("identifier_status") == "pubchem_and_inchikey":
        score += 2.0
    if LOW_TCM_SPECIFICITY.search(" ".join([row.get("ingredient_name", ""), row.get("alias", "")])):
        score -= 10.0
    return score


def classify(row: Dict[str, str], score: float) -> str:
    name_blob = " ".join([row.get("ingredient_name", ""), row.get("alias", "")])
    if LOW_TCM_SPECIFICITY.search(name_blob):
        return "deprioritize_low_tcm_specificity_or_conventional"
    if row.get("identifier_status") != "pubchem_and_inchikey":
        return "hold_missing_identifier"
    if int(row.get("reference_count") or 0) >= 2 and score >= 12:
        return "priority_for_lincs_mapping"
    if int(row.get("clinical_trial_count") or 0) >= 1 and int(row.get("reference_count") or 0) >= 1:
        return "priority_for_lincs_mapping"
    if int(row.get("reference_count") or 0) >= 1:
        return "secondary_literature_supported"
    return "background_or_control"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", default="results/m2_herb2_uc_ibd_evidence/uc_ibd_ingredient_evidence.tsv")
    parser.add_argument("--output", default="results/m2_herb2_uc_ibd_evidence/prioritized_ingredient_candidates.tsv")
    args = parser.parse_args()

    rows: List[Dict[str, str]] = []
    for row in read_tsv(Path(args.input)):
        score = priority_score(row)
        row = dict(row)
        row["priority_score"] = f"{score:.1f}"
        row["priority_class"] = classify(row, score)
        row["lincs_mapping_key"] = row.get("inchikey") or row.get("pubchem_id") or row.get("ingredient_name", "")
        row["review_note"] = (
            "manual safety and TCM-source review required before claims"
            if row["priority_class"] in {"priority_for_lincs_mapping", "secondary_literature_supported"}
            else "do not use as primary TCM candidate without special rationale"
        )
        rows.append(row)

    rows.sort(
        key=lambda item: (
            item["priority_class"] != "priority_for_lincs_mapping",
            -float(item["priority_score"]),
            item["ingredient_name"],
        )
    )

    fields = [
        "priority_class",
        "priority_score",
        "ingredient_id",
        "ingredient_name",
        "evidence_count",
        "reference_count",
        "clinical_trial_count",
        "meta_analysis_count",
        "matched_terms",
        "pubchem_id",
        "inchikey",
        "lincs_mapping_key",
        "cas_id",
        "drugbank_id",
        "identifier_status",
        "review_note",
        "top_evidence_titles",
    ]
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    with Path(args.output).open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
