#!/usr/bin/env python3
"""Check LINCS perturbagen coverage for prioritized HERB UC/IBD ingredients."""

from __future__ import annotations

import argparse
import csv
import gzip
from collections import Counter, defaultdict
from pathlib import Path
from typing import Dict, Iterable, List, Tuple


def open_text(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8", errors="replace", newline="")
    return path.open(encoding="utf-8", errors="replace", newline="")


def read_tsv(path: Path) -> Iterable[Dict[str, str]]:
    with open_text(path) as handle:
        yield from csv.DictReader(handle, delimiter="\t")


def norm(value: str) -> str:
    return (value or "").strip().lower()


def load_lincs_pert(paths: List[Tuple[str, Path]]) -> Tuple[List[Dict[str, str]], Dict[str, List[Dict[str, str]]], Dict[str, List[Dict[str, str]]]]:
    rows: List[Dict[str, str]] = []
    by_inchikey: Dict[str, List[Dict[str, str]]] = defaultdict(list)
    by_pubchem: Dict[str, List[Dict[str, str]]] = defaultdict(list)
    for source, path in paths:
        for row in read_tsv(path):
            item = {
                "lincs_source": source,
                "pert_id": row.get("pert_id", ""),
                "pert_iname": row.get("pert_iname", ""),
                "pert_type": row.get("pert_type", ""),
                "inchi_key": row.get("inchi_key", ""),
                "pubchem_cid": row.get("pubchem_cid", ""),
                "canonical_smiles": row.get("canonical_smiles", ""),
            }
            rows.append(item)
            if item["inchi_key"] and item["inchi_key"] not in {"-666", "NA"}:
                by_inchikey[item["inchi_key"]].append(item)
            if item["pubchem_cid"] and item["pubchem_cid"] not in {"-666", "NA"}:
                by_pubchem[item["pubchem_cid"]].append(item)
    return rows, by_inchikey, by_pubchem


def load_sig_counts(path: Path) -> Tuple[Counter, Dict[str, set]]:
    counts: Counter = Counter()
    cells: Dict[str, set] = defaultdict(set)
    if not path.exists():
        return counts, cells
    try:
        for row in read_tsv(path):
            pert_id = row.get("pert_id", "")
            if not pert_id:
                continue
            counts[pert_id] += 1
            if row.get("cell_id"):
                cells[pert_id].add(row["cell_id"])
    except (EOFError, OSError, gzip.BadGzipFile):
        return Counter(), defaultdict(set)
    return counts, cells


def match_candidate(row: Dict[str, str], by_inchikey: Dict[str, List[Dict[str, str]]], by_pubchem: Dict[str, List[Dict[str, str]]]) -> Tuple[str, List[Dict[str, str]]]:
    inchikey = row.get("inchikey", "")
    pubchem = row.get("pubchem_id", "")
    if inchikey and inchikey in by_inchikey:
        return "inchikey_exact", by_inchikey[inchikey]
    if pubchem and pubchem in by_pubchem:
        return "pubchem_exact", by_pubchem[pubchem]
    return "no_exact_identifier_match", []


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--candidates", default="results/m2_herb2_uc_ibd_evidence/prioritized_ingredient_candidates.tsv")
    parser.add_argument("--gse92742-pert", default="data/raw/lincs/GSE92742_Broad_LINCS_pert_info.txt.gz")
    parser.add_argument("--gse70138-pert", default="data/raw/lincs/GSE70138_Broad_LINCS_pert_info.txt.gz")
    parser.add_argument("--gse92742-sig", default="data/raw/lincs/GSE92742_Broad_LINCS_sig_info.txt.gz")
    parser.add_argument("--out-dir", default="results/m3_lincs_coverage")
    args = parser.parse_args()

    _, by_inchikey, by_pubchem = load_lincs_pert(
        [
            ("GSE92742", Path(args.gse92742_pert)),
            ("GSE70138", Path(args.gse70138_pert)),
        ]
    )
    sig_counts, sig_cells = load_sig_counts(Path(args.gse92742_sig))

    rows: List[Dict[str, str]] = []
    for candidate in read_tsv(Path(args.candidates)):
        match_method, matches = match_candidate(candidate, by_inchikey, by_pubchem)
        if matches:
            for match in matches:
                pert_id = match["pert_id"]
                rows.append(
                    {
                        "priority_class": candidate.get("priority_class", ""),
                        "priority_score": candidate.get("priority_score", ""),
                        "ingredient_id": candidate.get("ingredient_id", ""),
                        "ingredient_name": candidate.get("ingredient_name", ""),
                        "pubchem_id": candidate.get("pubchem_id", ""),
                        "inchikey": candidate.get("inchikey", ""),
                        "match_status": "covered",
                        "match_method": match_method,
                        "lincs_source": match["lincs_source"],
                        "pert_id": pert_id,
                        "pert_iname": match["pert_iname"],
                        "pert_type": match["pert_type"],
                        "gse92742_signature_count": str(sig_counts.get(pert_id, 0)),
                        "gse92742_cell_count": str(len(sig_cells.get(pert_id, set()))),
                        "gse92742_cells": ";".join(sorted(sig_cells.get(pert_id, set()))[:50]),
                    }
                )
        else:
            rows.append(
                {
                    "priority_class": candidate.get("priority_class", ""),
                    "priority_score": candidate.get("priority_score", ""),
                    "ingredient_id": candidate.get("ingredient_id", ""),
                    "ingredient_name": candidate.get("ingredient_name", ""),
                    "pubchem_id": candidate.get("pubchem_id", ""),
                    "inchikey": candidate.get("inchikey", ""),
                    "match_status": "not_covered_by_exact_id",
                    "match_method": match_method,
                    "lincs_source": "",
                    "pert_id": "",
                    "pert_iname": "",
                    "pert_type": "",
                    "gse92742_signature_count": "0",
                    "gse92742_cell_count": "0",
                    "gse92742_cells": "",
                }
            )

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    fields = [
        "priority_class",
        "priority_score",
        "ingredient_id",
        "ingredient_name",
        "pubchem_id",
        "inchikey",
        "match_status",
        "match_method",
        "lincs_source",
        "pert_id",
        "pert_iname",
        "pert_type",
        "gse92742_signature_count",
        "gse92742_cell_count",
        "gse92742_cells",
    ]
    with (out_dir / "lincs_candidate_coverage.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)

    unique_candidates = {(row["ingredient_id"], row["match_status"]) for row in rows}
    covered_priority = {
        row["ingredient_id"]
        for row in rows
        if row["priority_class"] == "priority_for_lincs_mapping" and row["match_status"] == "covered"
    }
    all_priority = {
        row["ingredient_id"]
        for row in rows
        if row["priority_class"] == "priority_for_lincs_mapping"
    }
    summary = [
        {"metric": "coverage_rows", "value": str(len(rows)), "notes": "Multiple rows per ingredient possible when present in both GSE92742 and GSE70138"},
        {"metric": "unique_candidate_status_pairs", "value": str(len(unique_candidates)), "notes": ""},
        {"metric": "priority_candidates", "value": str(len(all_priority)), "notes": ""},
        {"metric": "priority_candidates_covered", "value": str(len(covered_priority)), "notes": "Exact InChIKey or PubChem match"},
    ]
    with (out_dir / "summary.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["metric", "value", "notes"], delimiter="\t")
        writer.writeheader()
        writer.writerows(summary)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
