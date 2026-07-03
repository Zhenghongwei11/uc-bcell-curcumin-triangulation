#!/usr/bin/env python3
"""Build UC/IBD evidence tables from downloaded HERB 2.0 V2 files."""

from __future__ import annotations

import argparse
import csv
import re
from collections import Counter, defaultdict
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

EVIDENCE_FILES = {
    "reference": "HERB_reference_info_v2.txt",
    "clinical_trial": "HERB_clinical_trials_v2.txt",
    "meta_analysis": "HERB_meta_info_v2.txt",
    "experiment": "HERB_experiment_info_v2.txt",
    "disease": "HERB_disease_info_v2.txt",
}

REFERENCE_TEXT_FIELDS = [
    "Paper_title",
    "Paper_abstract",
    "Experiment_subject",
    "Experiment_type",
    "Phenotype_related",
]

CLINICAL_TEXT_FIELDS = [
    "NCT_title",
    "Study_condition",
    "Study_type",
    "Study_design",
    "Intervention",
    "Outcome_measure",
]

META_TEXT_FIELDS = [
    "CRD_title",
    "Review_question",
    "Condition_being_studied",
    "Participant",
    "Human_disease_modelled",
    "Intervention",
    "Comparator_control",
    "Main_outcome",
    "Outcome_measure",
    "Additional_outcome",
    "Keyword",
]

EXPERIMENT_TEXT_FIELDS = [
    "Subject_disease_name",
    "Experiment_type",
    "Sequence_type ",
    "Experiment_subject",
    "Experiment_detail",
    "Control_condition",
    "Experiment_subject_detail",
    "Data_type",
    "Tissue",
    "Cell_type",
    "Cell_line",
]

DISEASE_TEXT_FIELDS = [
    "Disease_name",
    "Disease_alias_name",
    "MeSH_disease_class",
    "HPO_disease_class",
    "DO_disease_class",
    "ICD10_id",
]


def normalize(text: str) -> str:
    return " ".join((text or "").lower().split())


def match_terms(text: str) -> Tuple[List[str], float]:
    haystack = normalize(text)
    terms = sorted({term for term in DIRECT_TERMS if normalize(term) in haystack})
    return terms, float(len(terms) * 3)


def read_tsv(path: Path) -> Iterable[Dict[str, str]]:
    with path.open(encoding="utf-8", errors="replace", newline="") as handle:
        yield from csv.DictReader(handle, delimiter="\t")


def text_from(row: Dict[str, str], fields: List[str]) -> str:
    return " | ".join(row.get(field, "") for field in fields if row.get(field, ""))


def load_entity_maps(data_dir: Path) -> Dict[str, Dict[str, Dict[str, str]]]:
    maps: Dict[str, Dict[str, Dict[str, str]]] = {"Ingredient": {}, "Herb": {}, "Formula": {}, "Disease": {}}
    for row in read_tsv(data_dir / "HERB_ingredient_info_v2.txt"):
        maps["Ingredient"][row["Ingredient_id"]] = row
    for row in read_tsv(data_dir / "HERB_herb_info_v2.txt"):
        maps["Herb"][row["Herb_id"]] = row
    for row in read_tsv(data_dir / "HERB_formula_info_v2.txt"):
        maps["Formula"][row["Formula_id"]] = row
    for row in read_tsv(data_dir / "HERB_disease_info_v2.txt"):
        maps["Disease"][row["Disease_id"]] = row
    return maps


def entity_name(subject_type: str, subject_id: str, subject_name: str, maps: Dict[str, Dict[str, Dict[str, str]]]) -> str:
    row = maps.get(subject_type, {}).get(subject_id, {})
    if subject_type == "Ingredient" and row:
        return row.get("Ingredient_name", subject_name)
    if subject_type == "Herb" and row:
        return " / ".join(v for v in [row.get("Herb_pinyin_name", ""), row.get("Herb_cn_name", ""), row.get("Herb_en_name", "")] if v and v != "NA")
    if subject_type == "Formula" and row:
        return " / ".join(v for v in [row.get("Formula_pinyin_name", ""), row.get("Formula_cn_name", ""), row.get("Formula_en_name", "")] if v and v != "NA")
    if subject_type == "Disease" and row:
        return row.get("Disease_name", subject_name)
    return subject_name


def scan_reference(data_dir: Path, maps: Dict[str, Dict[str, Dict[str, str]]]) -> List[Dict[str, str]]:
    rows: List[Dict[str, str]] = []
    for row in read_tsv(data_dir / EVIDENCE_FILES["reference"]):
        terms, score = match_terms(text_from(row, REFERENCE_TEXT_FIELDS))
        if not terms:
            continue
        subject_type = row.get("Subject_type", "")
        subject_id = row.get("Subject_id", "")
        rows.append(
            {
                "evidence_source": "HERB_reference",
                "evidence_id": row.get("Reference_id", ""),
                "subject_id": subject_id,
                "subject_type": subject_type,
                "subject_name": entity_name(subject_type, subject_id, row.get("Subject_name", ""), maps),
                "matched_terms": ";".join(terms),
                "match_score": f"{score:.1f}",
                "external_id": row.get("PubMed_id", ""),
                "title_or_condition": row.get("Paper_title", ""),
                "evidence_detail": row.get("Experiment_type", ""),
                "journal_or_registry": row.get("Journal", ""),
                "year_or_date": row.get("Publish_date", ""),
                "doi_or_url": row.get("DOI", ""),
            }
        )
    return rows


def scan_clinical(data_dir: Path, maps: Dict[str, Dict[str, Dict[str, str]]]) -> List[Dict[str, str]]:
    rows: List[Dict[str, str]] = []
    for row in read_tsv(data_dir / EVIDENCE_FILES["clinical_trial"]):
        terms, score = match_terms(text_from(row, CLINICAL_TEXT_FIELDS))
        if not terms:
            continue
        subject_type = row.get("Subject_type", "")
        subject_id = row.get("Subject_id", "")
        rows.append(
            {
                "evidence_source": "HERB_clinical_trial",
                "evidence_id": row.get("Clinical_trial_id", ""),
                "subject_id": subject_id,
                "subject_type": subject_type,
                "subject_name": entity_name(subject_type, subject_id, row.get("Subject_name", ""), maps),
                "matched_terms": ";".join(terms),
                "match_score": f"{score:.1f}",
                "external_id": row.get("NCT_id", ""),
                "title_or_condition": row.get("NCT_title", "") or row.get("Study_condition", ""),
                "evidence_detail": " | ".join(v for v in [row.get("Status", ""), row.get("Phase ", ""), row.get("Study_result", ""), row.get("Study_type", "")] if v and v != "NA"),
                "journal_or_registry": "ClinicalTrials.gov",
                "year_or_date": row.get("First_posted", ""),
                "doi_or_url": row.get("URL", ""),
            }
        )
    return rows


def scan_meta(data_dir: Path, maps: Dict[str, Dict[str, Dict[str, str]]]) -> List[Dict[str, str]]:
    rows: List[Dict[str, str]] = []
    for row in read_tsv(data_dir / EVIDENCE_FILES["meta_analysis"]):
        terms, score = match_terms(text_from(row, META_TEXT_FIELDS))
        if not terms:
            continue
        subject_type = row.get("Subject_type", "")
        subject_id = row.get("Subject_id", "")
        rows.append(
            {
                "evidence_source": "HERB_meta_analysis",
                "evidence_id": row.get("Meta-analysis_id", ""),
                "subject_id": subject_id,
                "subject_type": subject_type,
                "subject_name": entity_name(subject_type, subject_id, row.get("Subject_name", ""), maps),
                "matched_terms": ";".join(terms),
                "match_score": f"{score:.1f}",
                "external_id": row.get("CRD_id", ""),
                "title_or_condition": row.get("CRD_title", "") or row.get("Condition_being_studied", ""),
                "evidence_detail": row.get("Review_type", "") or row.get("Study_type_included", ""),
                "journal_or_registry": "PROSPERO/CRD",
                "year_or_date": row.get("Registration_date", "") or row.get("First_submission_date", ""),
                "doi_or_url": row.get("Final_publication", ""),
            }
        )
    return rows


def scan_experiment(data_dir: Path, maps: Dict[str, Dict[str, Dict[str, str]]]) -> List[Dict[str, str]]:
    rows: List[Dict[str, str]] = []
    for row in read_tsv(data_dir / EVIDENCE_FILES["experiment"]):
        terms, score = match_terms(text_from(row, EXPERIMENT_TEXT_FIELDS))
        if not terms:
            continue
        subject_type = row.get("Subject_type", "")
        subject_id = row.get("Subject_disease_id", "")
        rows.append(
            {
                "evidence_source": "HERB_high_throughput_experiment",
                "evidence_id": row.get("EXP_id", ""),
                "subject_id": subject_id,
                "subject_type": subject_type,
                "subject_name": entity_name(subject_type, subject_id, row.get("Subject_disease_name", ""), maps),
                "matched_terms": ";".join(terms),
                "match_score": f"{score:.1f}",
                "external_id": row.get("GSE_id", ""),
                "title_or_condition": row.get("Subject_disease_name", ""),
                "evidence_detail": " | ".join(v for v in [row.get("Organism", ""), row.get("Experiment_type", ""), row.get("Data_type", ""), row.get("Tissue", "")] if v and v != "NA"),
                "journal_or_registry": "GEO",
                "year_or_date": "",
                "doi_or_url": row.get("Original_id ", ""),
            }
        )
    return rows


def scan_disease(data_dir: Path) -> List[Dict[str, str]]:
    rows: List[Dict[str, str]] = []
    for row in read_tsv(data_dir / EVIDENCE_FILES["disease"]):
        terms, score = match_terms(text_from(row, DISEASE_TEXT_FIELDS))
        if not terms:
            continue
        rows.append(
            {
                "disease_id": row.get("Disease_id", ""),
                "disease_name": row.get("Disease_name", ""),
                "matched_terms": ";".join(terms),
                "match_score": f"{score:.1f}",
                "mesh_id": row.get("MeSH_id", ""),
                "do_id": row.get("DO_id", ""),
                "icd10_id": row.get("ICD10_id", ""),
                "alias": row.get("Disease_alias_name", ""),
            }
        )
    return rows


def ingredient_summary(evidence_rows: List[Dict[str, str]], maps: Dict[str, Dict[str, Dict[str, str]]]) -> List[Dict[str, str]]:
    grouped: Dict[str, List[Dict[str, str]]] = defaultdict(list)
    for row in evidence_rows:
        if row["subject_type"] == "Ingredient" and row["subject_id"]:
            grouped[row["subject_id"]].append(row)

    out: List[Dict[str, str]] = []
    for ingredient_id, rows in grouped.items():
        info = maps["Ingredient"].get(ingredient_id, {})
        evidence_counts = Counter(row["evidence_source"] for row in rows)
        out.append(
            {
                "ingredient_id": ingredient_id,
                "ingredient_name": info.get("Ingredient_name", rows[0].get("subject_name", "")),
                "alias": info.get("Ingredient_alias_name", ""),
                "evidence_count": str(len(rows)),
                "reference_count": str(evidence_counts.get("HERB_reference", 0)),
                "clinical_trial_count": str(evidence_counts.get("HERB_clinical_trial", 0)),
                "meta_analysis_count": str(evidence_counts.get("HERB_meta_analysis", 0)),
                "matched_terms": ";".join(sorted({term for row in rows for term in row["matched_terms"].split(";") if term})),
                "pubchem_id": info.get("PubChem_id", ""),
                "inchikey": info.get("InChIKey", ""),
                "canonical_smiles": info.get("Canonical_smiles", ""),
                "cas_id": info.get("CAS_id", ""),
                "drugbank_id": info.get("DrugBank_id", ""),
                "identifier_status": identifier_status(info),
                "top_evidence_titles": " || ".join(row["title_or_condition"] for row in rows[:5]),
            }
        )

    out.sort(key=lambda item: (-int(item["evidence_count"]), item["ingredient_name"]))
    return out


def identifier_status(info: Dict[str, str]) -> str:
    pubchem = info.get("PubChem_id", "")
    inchikey = info.get("InChIKey", "")
    smiles = info.get("Canonical_smiles", "")
    if pubchem and pubchem != "NA" and inchikey and inchikey != "NA":
        return "pubchem_and_inchikey"
    if pubchem and pubchem != "NA":
        return "pubchem_only"
    if inchikey and inchikey != "NA":
        return "inchikey_only"
    if smiles and smiles != "NA":
        return "smiles_only"
    return "missing_core_chemical_id"


def write_tsv(path: Path, rows: List[Dict[str, str]], fields: List[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def write_summary(path: Path, evidence_rows: List[Dict[str, str]], ingredient_rows: List[Dict[str, str]], disease_rows: List[Dict[str, str]]) -> None:
    counts = Counter(row["evidence_source"] for row in evidence_rows)
    subject_counts = Counter(row["subject_type"] for row in evidence_rows)
    identifier_counts = Counter(row["identifier_status"] for row in ingredient_rows)
    rows = [
        {"metric": "total_evidence_rows", "value": str(len(evidence_rows)), "notes": "reference + clinical trial + meta-analysis + high-throughput rows matching direct UC/IBD terms"},
        {"metric": "disease_records", "value": str(len(disease_rows)), "notes": "HERB disease records matching direct UC/IBD terms"},
        {"metric": "ingredient_candidates", "value": str(len(ingredient_rows)), "notes": "Ingredient subjects with direct UC/IBD evidence"},
    ]
    rows += [{"metric": f"evidence_source:{key}", "value": str(value), "notes": ""} for key, value in sorted(counts.items())]
    rows += [{"metric": f"subject_type:{key}", "value": str(value), "notes": ""} for key, value in sorted(subject_counts.items())]
    rows += [{"metric": f"ingredient_identifier:{key}", "value": str(value), "notes": ""} for key, value in sorted(identifier_counts.items())]
    write_tsv(path, rows, ["metric", "value", "notes"])


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", default="data/raw/herb2")
    parser.add_argument("--out-dir", default="results/m2_herb2_uc_ibd_evidence")
    args = parser.parse_args()

    data_dir = Path(args.data_dir)
    out_dir = Path(args.out_dir)
    maps = load_entity_maps(data_dir)

    evidence_rows = []
    evidence_rows.extend(scan_reference(data_dir, maps))
    evidence_rows.extend(scan_clinical(data_dir, maps))
    evidence_rows.extend(scan_meta(data_dir, maps))
    evidence_rows.extend(scan_experiment(data_dir, maps))
    evidence_rows.sort(key=lambda item: (item["subject_type"], item["subject_id"], item["evidence_source"], item["evidence_id"]))

    disease_rows = scan_disease(data_dir)
    ingredient_rows = ingredient_summary(evidence_rows, maps)

    evidence_fields = [
        "evidence_source",
        "evidence_id",
        "subject_id",
        "subject_type",
        "subject_name",
        "matched_terms",
        "match_score",
        "external_id",
        "title_or_condition",
        "evidence_detail",
        "journal_or_registry",
        "year_or_date",
        "doi_or_url",
    ]
    ingredient_fields = [
        "ingredient_id",
        "ingredient_name",
        "alias",
        "evidence_count",
        "reference_count",
        "clinical_trial_count",
        "meta_analysis_count",
        "matched_terms",
        "pubchem_id",
        "inchikey",
        "canonical_smiles",
        "cas_id",
        "drugbank_id",
        "identifier_status",
        "top_evidence_titles",
    ]
    disease_fields = ["disease_id", "disease_name", "matched_terms", "match_score", "mesh_id", "do_id", "icd10_id", "alias"]

    write_tsv(out_dir / "uc_ibd_subject_evidence.tsv", evidence_rows, evidence_fields)
    write_tsv(out_dir / "uc_ibd_ingredient_evidence.tsv", ingredient_rows, ingredient_fields)
    write_tsv(out_dir / "uc_ibd_disease_records.tsv", disease_rows, disease_fields)
    write_summary(out_dir / "summary.tsv", evidence_rows, ingredient_rows, disease_rows)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
