#!/usr/bin/env python3
"""Bridge HERB ingredient candidates to the GSE125527 scRNA cell-state axis.

The bridge combines HERB disease evidence, HERB detail API compound-target
relations when available, LINCS identifier coverage, and the GSE125527
patient-level scRNA cell-state localization result.
"""

from __future__ import annotations

import csv
import re
from collections import defaultdict
from pathlib import Path
from typing import Dict, Iterable, List


OUT_DIR = Path("results/m7_tcm_scrna_bridge")
DOC_PATH = Path("docs/M7_TCM_SCRNA_BRIDGE_RUN.md")

B_CELL_PATTERNS = [
    r"\bb[- ]?cell",
    r"\bb[- ]?cells",
    r"regulatory b",
    r"memory b",
    r"\bbcl[- ]?6\b",
    r"\bsyk\b",
    r"\bblnk\b",
    r"\bcd79",
    r"\bms4a1\b",
]

MYELOID_PATTERNS = [
    r"macrophage",
    r"monocyte",
    r"dendritic",
    r"myeloid",
    r"\bm1\b",
    r"\bm2\b",
    r"nlrp3",
    r"inflammasome",
    r"nf[- ]?κb",
    r"nf[- ]?kb",
    r"tlr4",
    r"myd88",
]

GENERIC_LOW_SPECIFICITY = {
    "vitamin d",
    "progesterone",
    "caffeine",
    "starch",
    "calcium",
    "hydrocortisone",
    "caproate",
}


def read_tsv(path: Path) -> List[Dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: Iterable[Dict[str, object]], fields: List[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def pattern_hits(text: str, patterns: List[str]) -> List[str]:
    hits: List[str] = []
    for pattern in patterns:
        if re.search(pattern, text, flags=re.IGNORECASE):
            hits.append(pattern)
    return hits


def evidence_type_counts(rows: List[Dict[str, str]]) -> Dict[str, int]:
    counts: Dict[str, int] = defaultdict(int)
    for row in rows:
        source = row.get("evidence_source", "")
        if source:
            counts[source] += 1
    return dict(counts)


def first_nonempty(*values: str) -> str:
    for value in values:
        if value and value != "NA":
            return value
    return ""


def optional_target_summary() -> Dict[str, Dict[str, str]]:
    path = Path("results/m8_herb_ingredient_targets/herb_ingredient_target_summary.tsv")
    if not path.exists():
        return {}
    return {row["ingredient_id"]: row for row in read_tsv(path)}


def optional_scrna_null_summary() -> Dict[str, Dict[str, str]]:
    path = Path("results/m6_scrna_cell_state/gse125527_disease_contrast_patient_label_null.tsv")
    if not path.exists():
        return {}
    return {row["celltype"]: row for row in read_tsv(path) if row.get("tissue_assignment") == "R"}


def main() -> int:
    candidates = read_tsv(Path("results/m2_herb2_uc_ibd_evidence/prioritized_ingredient_candidates.tsv"))
    evidence_rows = read_tsv(Path("results/m2_herb2_uc_ibd_evidence/uc_ibd_subject_evidence.tsv"))
    lincs_rows = read_tsv(Path("results/m3_lincs_coverage/lincs_candidate_coverage.tsv"))
    sc_contrast = read_tsv(Path("results/m6_scrna_cell_state/gse125527_disease_contrast_by_celltype.tsv"))
    target_summary = optional_target_summary()
    sc_null = optional_scrna_null_summary()

    evidence_by_subject: Dict[str, List[Dict[str, str]]] = defaultdict(list)
    for row in evidence_rows:
        if row.get("subject_type") == "Ingredient":
            evidence_by_subject[row["subject_id"]].append(row)

    lincs_by_subject: Dict[str, List[Dict[str, str]]] = defaultdict(list)
    for row in lincs_rows:
        lincs_by_subject[row["ingredient_id"]].append(row)

    sc_by_cell = {row["celltype"]: row for row in sc_contrast}
    b_axis_delta = sc_by_cell.get("B", {}).get("delta_axis_diseased_minus_healthy", "")
    b_axis_p = sc_by_cell.get("B", {}).get("p_axis_mannwhitney", "")
    b_axis_null_p = sc_null.get("B", {}).get("empirical_p_two_sided", "")
    mdc_axis_delta = sc_by_cell.get("M/DC", {}).get("delta_axis_diseased_minus_healthy", "")
    mdc_axis_p = sc_by_cell.get("M/DC", {}).get("p_axis_mannwhitney", "")
    mdc_axis_null_p = sc_null.get("M/DC", {}).get("empirical_p_two_sided", "")

    bridge_rows: List[Dict[str, object]] = []
    evidence_detail_rows: List[Dict[str, object]] = []

    for candidate in candidates:
        ingredient_id = candidate["ingredient_id"]
        name = candidate["ingredient_name"]
        candidate_evidence = evidence_by_subject.get(ingredient_id, [])

        b_records: List[Dict[str, str]] = []
        myeloid_records: List[Dict[str, str]] = []
        for evidence in candidate_evidence:
            text = " ".join(
                [
                    evidence.get("title_or_condition", ""),
                    evidence.get("evidence_detail", ""),
                    evidence.get("journal_or_registry", ""),
                ]
            )
            b_hits = pattern_hits(text, B_CELL_PATTERNS)
            myeloid_hits = pattern_hits(text, MYELOID_PATTERNS)
            if b_hits:
                b_records.append(evidence)
                evidence_detail_rows.append(
                    {
                        "ingredient_id": ingredient_id,
                        "ingredient_name": name,
                        "axis": "B_cell",
                        "evidence_id": evidence.get("evidence_id", ""),
                        "external_id": evidence.get("external_id", ""),
                        "title_or_condition": evidence.get("title_or_condition", ""),
                        "evidence_detail": evidence.get("evidence_detail", ""),
                        "year_or_date": evidence.get("year_or_date", ""),
                        "doi_or_url": evidence.get("doi_or_url", ""),
                        "matched_patterns": ";".join(b_hits),
                    }
                )
            if myeloid_hits:
                myeloid_records.append(evidence)
                evidence_detail_rows.append(
                    {
                        "ingredient_id": ingredient_id,
                        "ingredient_name": name,
                        "axis": "myeloid_dc",
                        "evidence_id": evidence.get("evidence_id", ""),
                        "external_id": evidence.get("external_id", ""),
                        "title_or_condition": evidence.get("title_or_condition", ""),
                        "evidence_detail": evidence.get("evidence_detail", ""),
                        "year_or_date": evidence.get("year_or_date", ""),
                        "doi_or_url": evidence.get("doi_or_url", ""),
                        "matched_patterns": ";".join(myeloid_hits),
                    }
                )

        counts = evidence_type_counts(candidate_evidence)
        lincs_records = lincs_by_subject.get(ingredient_id, [])
        lincs_covered = any(row.get("match_status") == "covered" for row in lincs_records)
        target_stats = target_summary.get(ingredient_id, {})
        literature_targets = int(target_stats.get("literature_unique_targets") or 0)
        database_targets = int(target_stats.get("database_unique_targets") or 0)
        up150_target_hits = int(target_stats.get("bulk_uc_ibd_up150_target_hits") or 0)
        down150_target_hits = int(target_stats.get("bulk_uc_ibd_down150_target_hits") or 0)
        target_count_bias = "yes" if database_targets > 500 else "no"
        low_specificity = name.strip().lower() in GENERIC_LOW_SPECIFICITY or candidate["priority_class"].startswith("deprioritize")
        priority_score = float(candidate.get("priority_score") or 0)
        bridge_score = (
            priority_score
            + 12 * len(b_records)
            + 4 * len(myeloid_records)
            + min(up150_target_hits, 20) * 1.5
            + min(literature_targets, 120) * 0.05
            + (2 if lincs_covered else 0)
        )
        if low_specificity:
            bridge_score -= 20
        if target_count_bias == "yes":
            bridge_score -= 8

        direct_b_cell_bridge = "yes" if b_records else "no"
        myeloid_dc_support = "yes" if myeloid_records else "no"
        if low_specificity:
            recommended_role = "exclude_or_background"
        elif b_records:
            recommended_role = "primary_lead_for_b_cell_axis"
        elif myeloid_records:
            recommended_role = "secondary_lead_for_myeloid_axis"
        elif candidate["priority_class"] == "priority_for_lincs_mapping":
            recommended_role = "backup_tcm_candidate"
        else:
            recommended_role = "low_priority"

        bridge_rows.append(
            {
                "recommended_role": recommended_role,
                "bridge_score": round(bridge_score, 3),
                "ingredient_id": ingredient_id,
                "ingredient_name": name,
                "priority_class": candidate.get("priority_class", ""),
                "priority_score": candidate.get("priority_score", ""),
                "evidence_count": candidate.get("evidence_count", ""),
                "reference_count": candidate.get("reference_count", ""),
                "clinical_trial_count": candidate.get("clinical_trial_count", ""),
                "meta_analysis_count": candidate.get("meta_analysis_count", ""),
                "b_cell_evidence_count": len(b_records),
                "myeloid_dc_evidence_count": len(myeloid_records),
                "lincs_exact_id_covered": "yes" if lincs_covered else "no",
                "herb_literature_unique_targets": literature_targets,
                "herb_database_unique_targets": database_targets,
                "bulk_uc_ibd_up150_target_hits": up150_target_hits,
                "bulk_uc_ibd_down150_target_hits": down150_target_hits,
                "target_count_bias_flag": target_count_bias,
                "up150_hit_genes": target_stats.get("up150_hit_genes", ""),
                "b_axis_delta": b_axis_delta,
                "b_axis_p": b_axis_p,
                "b_axis_patient_label_null_p": b_axis_null_p,
                "mdc_axis_delta": mdc_axis_delta,
                "mdc_axis_p": mdc_axis_p,
                "mdc_axis_patient_label_null_p": mdc_axis_null_p,
                "low_specificity_flag": "yes" if low_specificity else "no",
                "best_b_cell_evidence": first_nonempty(*(row.get("title_or_condition", "") for row in b_records[:2])),
                "best_myeloid_dc_evidence": first_nonempty(*(row.get("title_or_condition", "") for row in myeloid_records[:2])),
                "herb_evidence_sources": ";".join(f"{key}:{value}" for key, value in sorted(counts.items())),
            }
        )

    bridge_rows.sort(
        key=lambda row: (
            row["recommended_role"] != "primary_lead_for_b_cell_axis",
            row["recommended_role"] != "secondary_lead_for_myeloid_axis",
            -float(row["bridge_score"]),
        )
    )

    bridge_fields = [
        "recommended_role",
        "bridge_score",
        "ingredient_id",
        "ingredient_name",
        "priority_class",
        "priority_score",
        "evidence_count",
        "reference_count",
        "clinical_trial_count",
        "meta_analysis_count",
        "b_cell_evidence_count",
        "myeloid_dc_evidence_count",
        "lincs_exact_id_covered",
        "herb_literature_unique_targets",
        "herb_database_unique_targets",
        "bulk_uc_ibd_up150_target_hits",
        "bulk_uc_ibd_down150_target_hits",
        "target_count_bias_flag",
        "up150_hit_genes",
        "b_axis_delta",
        "b_axis_p",
        "b_axis_patient_label_null_p",
        "mdc_axis_delta",
        "mdc_axis_p",
        "mdc_axis_patient_label_null_p",
        "low_specificity_flag",
        "best_b_cell_evidence",
        "best_myeloid_dc_evidence",
        "herb_evidence_sources",
    ]
    detail_fields = [
        "ingredient_id",
        "ingredient_name",
        "axis",
        "evidence_id",
        "external_id",
        "title_or_condition",
        "evidence_detail",
        "year_or_date",
        "doi_or_url",
        "matched_patterns",
    ]

    write_tsv(OUT_DIR / "tcm_scrna_axis_bridge.tsv", bridge_rows, bridge_fields)
    write_tsv(OUT_DIR / "tcm_scrna_axis_evidence_details.tsv", evidence_detail_rows, detail_fields)

    top_rows = bridge_rows[:12]
    report_lines = [
        "# M7 TCM Candidate to scRNA Cell-Axis Bridge Run",
        "",
        "## Scope",
        "",
        "This run links HERB ingredient evidence to the GSE125527 rectal B-cell/M-DC disease-axis result. It uses local HERB evidence records and HERB detail API target relations when available; it does not infer compound-target edges from the standalone HERB target dictionary.",
        "",
        "## Inputs",
        "",
        "- `results/m2_herb2_uc_ibd_evidence/prioritized_ingredient_candidates.tsv`",
        "- `results/m2_herb2_uc_ibd_evidence/uc_ibd_subject_evidence.tsv`",
        "- `results/m3_lincs_coverage/lincs_candidate_coverage.tsv`",
        "- `results/m6_scrna_cell_state/gse125527_disease_contrast_by_celltype.tsv`",
        "- `results/m6_scrna_cell_state/gse125527_disease_contrast_patient_label_null.tsv` when available",
        "- `results/m8_herb_ingredient_targets/herb_ingredient_target_summary.tsv` when available",
        "",
        "## Key Outputs",
        "",
        "- `results/m7_tcm_scrna_bridge/tcm_scrna_axis_bridge.tsv`",
        "- `results/m7_tcm_scrna_bridge/tcm_scrna_axis_evidence_details.tsv`",
        "",
        "## scRNA Anchor",
        "",
        f"- Rectal B-cell disease axis: delta={b_axis_delta}, Mann-Whitney P={b_axis_p}, patient-label null P={b_axis_null_p}.",
        f"- Rectal M/DC disease axis: delta={mdc_axis_delta}, Mann-Whitney P={mdc_axis_p}, patient-label null P={mdc_axis_null_p}; supportive only.",
        "",
        "## Top Bridge Candidates",
        "",
        "| Role | Ingredient | Score | B-cell evidence | Myeloid/DC evidence | Literature targets | UC/IBD up-target hits | LINCS covered |",
        "|---|---|---:|---:|---:|---:|---:|---|",
    ]
    for row in top_rows:
        report_lines.append(
            f"| {row['recommended_role']} | {row['ingredient_name']} | {row['bridge_score']} | {row['b_cell_evidence_count']} | {row['myeloid_dc_evidence_count']} | {row['herb_literature_unique_targets']} | {row['bulk_uc_ibd_up150_target_hits']} | {row['lincs_exact_id_covered']} |"
        )
    report_lines.extend(
        [
            "",
            "## Interpretation",
            "",
            "- Curcumin is the strongest current primary lead because it combines clinical/meta-analysis HERB evidence, direct B-cell-related colitis literature, HERB detail API target support, and exact-ID LINCS coverage, even though the L1000CDS2 reversal query did not rank it as a top reversal hit.",
            "- Berberine remains a high-priority supporting lead because it has clinical/meta-analysis evidence, HERB detail API target support, and myeloid/macrophage-related mechanistic evidence, but it is better framed as secondary to the B-cell axis unless direct B-cell evidence is added.",
            "- Very large database-integrated target sets are flagged as target-count-biased and should not be allowed to dominate the story without literature-supported targets and independent disease evidence.",
            "- Compounds flagged as low specificity should not drive the TCM narrative even if HERB evidence counts are high.",
            "- This bridge is literature/cell-axis triangulation, not wet-lab validation and not a compound-target proof.",
            "",
        ]
    )
    DOC_PATH.write_text("\n".join(report_lines), encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
