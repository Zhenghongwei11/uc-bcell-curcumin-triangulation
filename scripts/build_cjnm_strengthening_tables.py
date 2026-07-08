#!/usr/bin/env python3
"""Build CJNM-oriented strengthening tables.

The outputs focus on two review risks for a public-data natural-product
pharmacology manuscript: therapeutic directionality and candidate specificity.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]


def read_tsv(relpath: str) -> pd.DataFrame:
    return pd.read_csv(ROOT / relpath, sep="\t", dtype=str).fillna("")


def write_tsv(df: pd.DataFrame, relpath: str) -> None:
    path = ROOT / relpath
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, sep="\t", index=False)


def fmt_delta(value: str) -> str:
    if value == "" or value.lower() == "nan":
        return "NA"
    return f"{float(value):.4g}"


def disease_direction(row: pd.Series) -> tuple[str, str]:
    delta_text = row.get("b_delta_log1p_cpm_diseased_minus_healthy", "")
    if delta_text != "" and delta_text.lower() != "nan":
        delta = float(delta_text)
        if delta > 0:
            return "Increased in diseased rectal B cells", fmt_delta(delta_text)
        if delta < 0:
            return "Decreased in diseased rectal B cells", fmt_delta(delta_text)
        return "No rectal B-cell change", fmt_delta(delta_text)
    if row.get("bulk_module_role") == "bulk_up_top150":
        return "Increased in replicated bulk UC/IBD signature", "Not resolved"
    return "Not resolved in rectal B-cell pseudobulk", "Not resolved"


def consistency_call(gene: str, disease: str, relation: str, role: str) -> tuple[str, str, str]:
    relation = relation.lower()
    if "down" in relation and ("Increased" in disease or "bulk UC/IBD" in disease):
        expected = "Suppression of a disease-increased target or cytokine-context signal"
        call = "Directionally consistent"
        interpretation = (
            "Curcumin-reported downregulation aligns with a disease-increased signal; "
            "this supports a normalization-oriented hypothesis, not proven target engagement."
        )
    elif "up" in relation and "Increased" in disease:
        expected = "Context modulation rather than simple normalization"
        call = "Contextual, not normalizing"
        interpretation = (
            "The reported curcumin-associated increase occurs in a B-cell signaling context "
            "that is also higher in diseased B cells; treat as pathway-state evidence rather "
            "than a claim of direct reversal."
        )
    elif "up" in relation:
        expected = "Potential restoration or immune-regulatory context"
        call = "Context-dependent"
        interpretation = (
            "The reported upregulation is compatible with immune-regulatory context but is "
            "not used as a direct reversal claim."
        )
    else:
        expected = "No directionality inference"
        call = "Insufficient directionality"
        interpretation = "Directionality could not be resolved from the current evidence layer."

    if gene in {"CCL2", "IL33", "IL1B", "TNF"}:
        interpretation += " This target is retained as cytokine-context support."
    elif gene in {"BCL6", "BLNK", "SYK"}:
        interpretation += " This target is retained in the B-cell mechanism-candidate layer."
    elif role:
        interpretation += f" Evidence role: {role}."
    return expected, call, interpretation


def build_directionality_table() -> pd.DataFrame:
    table3 = read_tsv("tables/source_data/table3_curcumin_target_evidence.tsv")
    focus = ["BCL6", "BLNK", "SYK", "CCL2", "IL33", "IL1B", "TNF"]
    table3 = table3[table3["gene_symbol"].isin(focus)].copy()

    role_map = {
        "BCL6": "B-cell mechanism candidate",
        "BLNK": "B-cell mechanism candidate",
        "SYK": "B-cell mechanism candidate",
        "CCL2": "Cytokine-context target",
        "IL33": "Cytokine-context target",
        "IL1B": "Cytokine-context target",
        "TNF": "Cytokine-context target",
    }
    order = {gene: idx for idx, gene in enumerate(focus)}
    rows: list[dict[str, str]] = []
    for _, row in table3.sort_values("gene_symbol", key=lambda s: s.map(order)).iterrows():
        gene = row["gene_symbol"]
        disease, delta = disease_direction(row)
        relation = row["relationship_directions"] or "not resolved"
        expected, call, interpretation = consistency_call(gene, disease, relation, role_map[gene])
        rows.append(
            {
                "gene": gene,
                "evidence_role": role_map[gene],
                "disease_direction": disease,
                "rectal_B_delta_log1p_CPM": delta,
                "curcumin_reported_regulation": relation,
                "expected_therapeutic_direction": expected,
                "directional_consistency": call,
                "evidence_source": row.get("pubmed_ids", ""),
                "interpretation_for_review": interpretation,
            }
        )
    return pd.DataFrame(rows)


def build_specificity_control_table() -> pd.DataFrame:
    table2 = read_tsv("tables/source_data/table2_prioritized_tcm_candidates.tsv")
    m7 = read_tsv("results/m7_tcm_scrna_bridge/tcm_scrna_axis_bridge.tsv")
    specificity = read_tsv("results/m17_no_download_robustness/herb_candidate_specificity_ranking.tsv")
    baseline = read_tsv("results/m18_reviewer_lightweight_extensions/standard_overlap_baseline_candidate_ranking.tsv")
    role_map = {
        "primary_lead_for_b_cell_axis": "Lead B-cell-axis candidate",
        "secondary_lead_for_myeloid_axis": "Secondary myeloid-axis candidate",
        "backup_tcm_candidate": "Lower-consensus candidate",
    }
    m7_meta = m7[
        [
            "ingredient_id",
            "recommended_role",
            "evidence_count",
            "reference_count",
            "clinical_trial_count",
            "meta_analysis_count",
            "lincs_exact_id_covered",
        ]
    ].copy()
    m7_meta["m7_display_role"] = m7_meta["recommended_role"].map(role_map).fillna(m7_meta["recommended_role"])

    merged = (
        specificity.merge(
            baseline[
                [
                    "ingredient_id",
                    "standard_overlap_rank",
                    "standard_overlap_hits",
                    "standard_overlap_up_hits",
                    "standard_overlap_genes",
                    "disease_cell_aware_rank",
                    "rank_shift_overlap_minus_disease_cell",
                ]
            ],
            on="ingredient_id",
            how="left",
        )
        .merge(
            table2[
                [
                    "ingredient_id",
                    "ingredient_name",
                    "display_role",
                    "evidence_count",
                    "reference_count",
                    "clinical_trial_count",
                    "meta_analysis_count",
                    "lincs_exact_id_covered",
                ]
            ],
            on="ingredient_id",
            how="left",
            suffixes=("", "_table2"),
        )
        .merge(m7_meta, on="ingredient_id", how="left", suffixes=("", "_m7"))
    )
    table2_name = merged["ingredient_name_table2"].fillna("")
    merged["ingredient_name"] = table2_name.where(table2_name.ne(""), merged["ingredient_name"])
    role = merged["display_role"].fillna("")
    role = role.where(role.ne(""), merged["m7_display_role"].fillna(""))
    merged["display_role"] = role.where(role.ne(""), "Specificity-control candidate")
    for col in [
        "evidence_count",
        "reference_count",
        "clinical_trial_count",
        "meta_analysis_count",
        "lincs_exact_id_covered",
    ]:
        fallback = merged[f"{col}_m7"].fillna("") if f"{col}_m7" in merged.columns else ""
        merged[col] = merged[col].fillna("").where(merged[col].fillna("").ne(""), fallback)
    selected_names = {
        "Curcumin",
        "Archin",
        "Berberine",
        "Bilobalide",
        "Geniposide",
        "Quercetin",
        "Trans-resveratrol",
        "Honokiol",
        "Polydatin",
        "Taxifolin",
    }
    selected = merged[
        merged["ingredient_name"].isin(selected_names)
        | merged["target_set_specificity_rank"].astype(int).le(8)
        | merged["standard_overlap_rank"].astype(float).le(8)
    ].copy()
    selected = selected.sort_values(
        ["target_set_specificity_rank", "standard_overlap_rank"],
        key=lambda s: pd.to_numeric(s, errors="coerce"),
    )
    cols = [
        "ingredient_name",
        "display_role",
        "evidence_count",
        "reference_count",
        "clinical_trial_count",
        "meta_analysis_count",
        "disease_context_targets",
        "target_set_specificity_rank",
        "target_set_specificity_score_without_direction_alignment",
        "standard_overlap_rank",
        "standard_overlap_hits",
        "standard_overlap_up_hits",
        "disease_cell_aware_rank",
        "rank_shift_overlap_minus_disease_cell",
        "b_cell_increased",
        "mdc_expressed",
        "lincs_exact_id_covered",
        "standard_overlap_genes",
    ]
    out = selected[cols].copy()
    out = out.rename(
        columns={
            "ingredient_name": "candidate",
            "display_role": "candidate_role",
            "evidence_count": "HERB_records",
            "reference_count": "HERB_references",
            "clinical_trial_count": "clinical_trial_records",
            "meta_analysis_count": "meta_analysis_records",
            "disease_context_targets": "disease_context_targets",
            "target_set_specificity_rank": "B_cell_axis_specificity_rank",
            "target_set_specificity_score_without_direction_alignment": "B_cell_axis_specificity_score",
            "standard_overlap_rank": "standard_overlap_rank",
            "standard_overlap_hits": "standard_overlap_hits",
            "standard_overlap_up_hits": "standard_overlap_up_hits",
            "disease_cell_aware_rank": "disease_cell_aware_rank",
            "rank_shift_overlap_minus_disease_cell": "rank_shift_overlap_minus_disease_cell",
            "b_cell_increased": "B_cell_increased_targets",
            "mdc_expressed": "MDC_expressed_targets",
            "lincs_exact_id_covered": "LINCS_exact_ID_covered",
            "standard_overlap_genes": "standard_overlap_genes",
        }
    )
    return out


def build_specificity_summary() -> pd.DataFrame:
    null = read_tsv("results/m17_no_download_robustness/herb_decoy_null_summary.tsv")
    row = null.iloc[0]
    return pd.DataFrame(
        [
            {
                "analysis": "Real-candidate ranking",
                "result": "Curcumin ranked first among 25 HERB disease-context candidates by B-cell-axis specificity score.",
                "interpretation": "Supports lead-candidate status within the real candidate landscape.",
            },
            {
                "analysis": "Standard target-overlap baseline",
                "result": "Curcumin ranked second by simple bulk-module overlap, whereas Archin ranked first.",
                "interpretation": "The disease-cell-aware framework changes the ranking away from generic inflammatory overlap.",
            },
            {
                "analysis": "Matched random target-set null",
                "result": (
                    f"Observed score {row['observed_target_set_score_without_direction_alignment']}; "
                    f"null mean {row['null_mean']}; empirical P>={row['empirical_p_ge_observed']}."
                ),
                "interpretation": "Does not eliminate target-count or literature-density bias; this limit should remain explicit.",
            },
            {
                "analysis": "Target-count association",
                "result": (
                    f"Spearman rho {float(row['target_count_score_spearman_rho_among_candidates']):.3f}; "
                    f"P={float(row['target_count_score_spearman_p']):.3g}."
                ),
                "interpretation": "Candidate scores remain strongly coupled to disease-context target count.",
            },
        ]
    )


def manuscript_directionality_table(direction: pd.DataFrame) -> pd.DataFrame:
    cols = [
        "gene",
        "evidence_role",
        "disease_direction",
        "rectal_B_delta_log1p_CPM",
        "curcumin_reported_regulation",
        "expected_therapeutic_direction",
        "directional_consistency",
        "evidence_source",
    ]
    out = direction[cols].copy()
    out = out.rename(
        columns={
            "evidence_role": "evidence role",
            "disease_direction": "disease direction",
            "rectal_B_delta_log1p_CPM": "rectal B delta",
            "curcumin_reported_regulation": "curcumin-reported regulation",
            "expected_therapeutic_direction": "expected therapeutic direction",
            "directional_consistency": "directional consistency",
            "evidence_source": "PMID",
        }
    )
    return out


def main() -> None:
    direction = build_directionality_table()
    direction_main = manuscript_directionality_table(direction)
    controls = build_specificity_control_table()
    summary = build_specificity_summary()

    write_tsv(direction, "tables/source_data/table4_therapeutic_direction_consistency.tsv")
    write_tsv(direction_main, "tables/manuscript/table4_therapeutic_direction_consistency.tsv")
    write_tsv(direction, "tables/supplementary/cjnm_therapeutic_direction_consistency_full.tsv")

    write_tsv(controls, "tables/supplementary/cjnm_candidate_specificity_controls.tsv")
    write_tsv(summary, "tables/supplementary/cjnm_specificity_summary.tsv")
    print("Wrote CJNM strengthening tables.")


if __name__ == "__main__":
    main()
