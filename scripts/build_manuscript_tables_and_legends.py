#!/usr/bin/env python3
"""Build manuscript-facing tables and figure legend drafts from analysis outputs."""

from __future__ import annotations

import os
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
TABLE_DIR = ROOT / "tables" / "source_data"
MANUSCRIPT_TABLE_DIR = ROOT / "tables" / "manuscript"
DOCS_DIR = ROOT / "docs"


def read_tsv(relpath: str) -> pd.DataFrame:
    return pd.read_csv(ROOT / relpath, sep="\t")


def write_tsv(df: pd.DataFrame, relpath: str) -> None:
    path = ROOT / relpath
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, sep="\t", index=False)


def write_curated_tsv(df: pd.DataFrame, relpath: str) -> None:
    path = ROOT / relpath
    if path.exists():
        print(f"Preserved curated {path.relative_to(ROOT)}")
        return
    write_tsv(df, relpath)


def standardize_candidate_names(names: pd.Series) -> list[str]:
    standardized = names.astype(str).replace({"Berberime": "Berberine"})
    seen: dict[str, int] = {}
    out: list[str] = []
    for name in standardized:
        seen[name] = seen.get(name, 0) + 1
        out.append(name if seen[name] == 1 else f"{name} (record {seen[name]})")
    return out


def build_manuscript_table_2(table: pd.DataFrame) -> pd.DataFrame:
    table = table.copy()
    table["ingredient_name"] = standardize_candidate_names(table["ingredient_name"])
    cols = [
        "display_role",
        "ingredient_name",
        "bridge_score",
        "evidence_count",
        "reference_count",
        "clinical_trial_count",
        "meta_analysis_count",
        "b_cell_evidence_count",
        "myeloid_dc_evidence_count",
        "lincs_exact_id_covered",
        "n_disease_context_unique_targets",
        "n_bulk_up_target_hits",
        "n_b_cell_expressed_targets",
        "n_b_cell_disease_increased_targets",
        "n_mdc_disease_increased_targets",
    ]
    out = table[cols].copy()
    out = out.rename(
        columns={
            "display_role": "candidate_role",
            "ingredient_name": "candidate",
            "bridge_score": "candidate_priority_score",
            "evidence_count": "HERB_evidence_records",
            "reference_count": "HERB_references",
            "clinical_trial_count": "HERB_clinical_trial_records",
            "meta_analysis_count": "HERB_meta_analysis_records",
            "b_cell_evidence_count": "B_cell_literature_records",
            "myeloid_dc_evidence_count": "myeloid_DC_literature_records",
            "lincs_exact_id_covered": "LINCS_exact_ID_covered",
            "n_disease_context_unique_targets": "disease_context_targets",
            "n_bulk_up_target_hits": "bulk_up_module_targets",
            "n_b_cell_expressed_targets": "rectal_B_cell_expressed_targets",
            "n_b_cell_disease_increased_targets": "rectal_B_cell_increased_targets",
            "n_mdc_disease_increased_targets": "rectal_MDC_increased_targets",
        }
    )
    return out


def build_manuscript_table_3(table: pd.DataFrame) -> pd.DataFrame:
    figure_use_map = {
        "main_figure_candidate": "Main B-cell mechanism",
        "cytokine_context_candidate": "Cytokine-context evidence",
        "main_or_supplementary": "Main or supplementary inflammatory context",
        "supplementary_candidate": "Supplementary context",
    }
    fulltext_level_map = {
        "open_fulltext_body_match": "Open full-text body match",
        "user_provided_fulltext_body_match": "Restricted full-text body match",
    }

    def map_figure_use(value: str) -> str:
        if pd.isna(value) or str(value).strip() == "":
            return "Not used in main mechanism figure"
        parts = str(value).split(";")
        return "; ".join(figure_use_map.get(part, part) for part in parts)

    def map_fulltext_level(value: str) -> str:
        if pd.isna(value) or str(value).strip() == "":
            return "No full-text body match in high-value verification"
        parts = str(value).split(";")
        return "; ".join(fulltext_level_map.get(part, part) for part in parts)

    out = table[
        [
            "gene_symbol",
            "display_evidence_layer",
            "display_bulk_role",
            "relationship_directions",
            "pubmed_ids",
            "b_delta_log1p_cpm_diseased_minus_healthy",
            "b_p_log1p_cpm_mannwhitney",
            "mdc_delta_log1p_cpm_diseased_minus_healthy",
            "mdc_p_log1p_cpm_mannwhitney",
            "open_targets_diseases",
            "max_open_targets_overall",
            "max_open_targets_genetic",
            "max_open_targets_clinical",
            "pubmed_title_matches",
            "pubmed_target_terms_found",
            "fulltext_evidence_levels",
            "recommended_figure_use",
            "interpretation_boundary",
        ]
    ].copy()
    out["recommended_figure_use"] = out["recommended_figure_use"].map(map_figure_use)
    out["fulltext_evidence_levels"] = out["fulltext_evidence_levels"].map(map_fulltext_level)
    out = out.rename(
        columns={
            "gene_symbol": "gene",
            "display_evidence_layer": "evidence_layer",
            "display_bulk_role": "bulk_signature_role",
            "relationship_directions": "HERB_relation_direction",
            "pubmed_ids": "supporting_PMIDs",
            "b_delta_log1p_cpm_diseased_minus_healthy": "rectal_B_cell_delta_log1p_CPM",
            "b_p_log1p_cpm_mannwhitney": "rectal_B_cell_P",
            "mdc_delta_log1p_cpm_diseased_minus_healthy": "rectal_MDC_delta_log1p_CPM",
            "mdc_p_log1p_cpm_mannwhitney": "rectal_MDC_P",
            "open_targets_diseases": "Open_Targets_disease_contexts",
            "max_open_targets_overall": "max_Open_Targets_overall_score",
            "max_open_targets_genetic": "max_Open_Targets_genetic_score",
            "max_open_targets_clinical": "max_Open_Targets_clinical_score",
            "pubmed_title_matches": "PubMed_title_matches",
            "pubmed_target_terms_found": "PubMed_target_terms_found",
            "fulltext_evidence_levels": "full_text_evidence_level",
            "recommended_figure_use": "figure_evidence_role",
            "interpretation_boundary": "interpretation_boundary",
        }
    )
    return out


def fmt_float(value, digits: int = 3) -> str:
    if pd.isna(value):
        return "NA"
    return f"{float(value):.{digits}g}"


def yes_no_count(series: pd.Series, yes_value: str = "yes") -> int:
    return int((series.astype(str) == yes_value).sum())


def build_table_1() -> pd.DataFrame:
    dataset = read_tsv("docs/DATASET_LANDSCAPE.tsv")
    manifest = read_tsv("data/manifest.tsv")

    def dataset_row(dataset_id: str, role: str, limitation: str) -> dict:
        row = dataset[dataset["dataset_id_or_name"] == dataset_id].iloc[0]
        return {
            "resource": dataset_id,
            "type": row["modality"],
            "sample_or_record_count": row["sample_size"],
            "tissue_or_scope": row["tissue_or_site"],
            "manuscript_role": role,
            "primary_limitation": limitation,
            "source_url_or_record": row["access_link"],
        }

    herb_files = manifest[manifest["dataset_id"].astype(str).str.contains("herb", case=False, na=False)]
    etcm_files = manifest[manifest["dataset_id"].astype(str).str.contains("etcm2", case=False, na=False)]
    rows = [
        dataset_row(
            "GSE75214",
            "Bulk mucosal disease-signature discovery evidence group",
            "GSE75214 and GSE59071 are not treated as independent validation because they share a platform/evidence group.",
        ),
        dataset_row(
            "GSE59071",
            "Non-independent bulk sensitivity evidence",
            "Used as sensitivity within the same GPL6244 evidence group rather than as an independent cohort.",
        ),
        dataset_row(
            "GSE87466",
            "Independent active UC mucosal validation evidence group",
            "Processed microarray analysis; detailed treatment covariates remain limited.",
        ),
        dataset_row(
            "GSE125527",
            "Primary scRNA pseudobulk cell-state localization resource",
            "Patient count is modest after tissue/cell-type stratification; claims remain associative.",
        ),
        {
            "resource": "HERB 2.0",
            "type": "TCM compound-target and disease-context knowledgebase",
            "sample_or_record_count": f"{len(herb_files)} local manifest files",
            "tissue_or_scope": "Ingredients, herbs, formulas, targets, references",
            "manuscript_role": "Primary TCM disease-context target source",
            "primary_limitation": "Database relations are evidence annotations, not experimental validation in this study.",
            "source_url_or_record": "http://herb.ac.cn/",
        },
        {
            "resource": "ETCM2",
            "type": "TCM compound and target-name resource",
            "sample_or_record_count": f"{len(etcm_files)} local manifest files",
            "tissue_or_scope": "Ingredient details and target-name annotations",
            "manuscript_role": "Supplementary cross-resource check",
            "primary_limitation": "Target names require gene-symbol mapping and did not support the Curcumin-B-cell target bridge.",
            "source_url_or_record": "http://www.tcmip.cn/ETCM2/front/",
        },
        {
            "resource": "Open Targets",
            "type": "Disease-target association resource",
            "sample_or_record_count": "242 disease-target rows in local output",
            "tissue_or_scope": "UC, Crohn disease, and IBD target associations",
            "manuscript_role": "Orthogonal target disease-relevance check",
            "primary_limitation": "Supports target-disease plausibility, not compound action.",
            "source_url_or_record": "https://platform.opentargets.org/",
        },
        {
            "resource": "PubMed/PMC and full-text verification",
            "type": "Bibliographic and full-text evidence verification",
            "sample_or_record_count": "3 high-value PMIDs; 14 target-level full-text rows",
            "tissue_or_scope": "Curcumin colitis target evidence",
            "manuscript_role": "Reference identity and target-context verification",
            "primary_limitation": "Restricted full text was used only for local evidence inspection and is not redistributed.",
            "source_url_or_record": "https://pubmed.ncbi.nlm.nih.gov/",
        },
        {
            "resource": "L1000CDS2",
            "type": "Perturbational signature reversal query",
            "sample_or_record_count": "15 exact-ID candidate checks",
            "tissue_or_scope": "LINCS/L1000 compound signatures",
            "manuscript_role": "Negative boundary-setting perturbation layer",
            "primary_limitation": "No exact-ID top-result hit does not rule out activity in untested contexts or raw LINCS matrices.",
            "source_url_or_record": "https://maayanlab.cloud/L1000CDS2/",
        },
    ]
    return pd.DataFrame(rows)


def build_table_2() -> pd.DataFrame:
    m7 = read_tsv("results/m7_tcm_scrna_bridge/tcm_scrna_axis_bridge.tsv")
    m9 = read_tsv("results/m9_target_cell_bridge/candidate_target_cell_bridge_summary.tsv")
    merged = m7.merge(
        m9[
            [
                "ingredient_id",
                "n_disease_context_unique_targets",
                "n_bulk_up_target_hits",
                "n_b_cell_expressed_targets",
                "n_b_cell_disease_increased_targets",
                "n_mdc_disease_increased_targets",
            ]
        ],
        on="ingredient_id",
        how="left",
    )
    role_map = {
        "primary_lead_for_b_cell_axis": "Lead B-cell-axis candidate",
        "secondary_lead_for_myeloid_axis": "Secondary myeloid-axis candidate",
        "backup_tcm_candidate": "Lower-consensus candidate",
    }
    merged["display_role"] = merged["recommended_role"].map(role_map).fillna(merged["recommended_role"])
    merged["ingredient_name"] = standardize_candidate_names(merged["ingredient_name"])
    cols = [
        "display_role",
        "recommended_role",
        "bridge_score",
        "ingredient_id",
        "ingredient_name",
        "evidence_count",
        "reference_count",
        "clinical_trial_count",
        "meta_analysis_count",
        "b_cell_evidence_count",
        "myeloid_dc_evidence_count",
        "lincs_exact_id_covered",
        "n_disease_context_unique_targets",
        "n_bulk_up_target_hits",
        "n_b_cell_expressed_targets",
        "n_b_cell_disease_increased_targets",
        "n_mdc_disease_increased_targets",
        "target_count_bias_flag",
        "low_specificity_flag",
    ]
    table = merged[cols].sort_values("bridge_score", ascending=False).head(15).copy()
    table["bridge_score"] = table["bridge_score"].round(3)
    return table


def build_table_3() -> pd.DataFrame:
    cur = read_tsv("results/m9_target_cell_bridge/curcumin_bcell_target_evidence.tsv")
    ot = read_tsv("results/m10_open_targets/open_targets_ibd_target_evidence.tsv")
    m12 = read_tsv("results/m12_pubmed_verification/curcumin_target_edge_pubmed_verification.tsv")
    m16 = read_tsv("results/m16_curcumin_fulltext_evidence/curcumin_high_value_target_fulltext_audit.tsv")

    ot_cur = ot[ot["is_curcumin_target"] == "yes"].copy()
    ot_gene = (
        ot_cur.groupby("gene_symbol", as_index=False)
        .agg(
            open_targets_diseases=("disease_name", lambda x: ";".join(sorted(set(map(str, x))))),
            max_open_targets_overall=("overall_score", "max"),
            max_open_targets_genetic=("genetic_association_score", "max"),
            max_open_targets_clinical=("clinical_score", "max"),
            max_open_targets_literature=("literature_score", "max"),
        )
    )

    m12_gene = (
        m12.groupby("gene_symbol", as_index=False)
        .agg(
            pubmed_edges_checked=("pmid", "size"),
            pubmed_pmids=("pmid", lambda x: ";".join(map(str, sorted(set(x))))),
            pubmed_title_matches=("title_match", yes_no_count),
            pubmed_target_terms_found=("target_term_found_in_pubmed_title_or_abstract", yes_no_count),
        )
    )

    m16_gene = (
        m16.groupby("gene_symbol", as_index=False)
        .agg(
            fulltext_pmids=("pmid", lambda x: ";".join(map(str, sorted(set(x))))),
            fulltext_evidence_levels=("evidence_level", lambda x: ";".join(sorted(set(map(str, x))))),
            recommended_figure_use=("recommended_figure_use", lambda x: ";".join(sorted(set(map(str, x))))),
            interpretation_boundary=("interpretation_boundary", lambda x: " | ".join(dict.fromkeys(map(str, x)))),
        )
    )

    table = cur.merge(ot_gene, on="gene_symbol", how="left").merge(m12_gene, on="gene_symbol", how="left").merge(m16_gene, on="gene_symbol", how="left")
    tier_map = {
        "A_module_aligned_b_cell_increased": "Bulk up-module and B-cell increased",
        "B_module_aligned": "Bulk up-module cytokine-context target",
        "C_b_cell_mechanism_increased": "B-cell mechanism candidate",
        "D_b_cell_increased": "B-cell increased supplementary target",
        "F_context_target_only": "Context-only literature target",
    }
    bulk_map = {
        "bulk_up_top150": "UC/IBD up-module",
        "not_bulk_top150": "Not in top bulk module",
    }
    table["display_evidence_layer"] = table["evidence_tier"].map(tier_map).fillna(table["evidence_tier"])
    table["display_bulk_role"] = table["bulk_module_role"].map(bulk_map).fillna(table["bulk_module_role"])
    table["b_delta_log1p_cpm_diseased_minus_healthy"] = table["b_delta_log1p_cpm_diseased_minus_healthy"].map(lambda x: fmt_float(x, 4))
    table["b_p_log1p_cpm_mannwhitney"] = table["b_p_log1p_cpm_mannwhitney"].map(lambda x: fmt_float(x, 3))
    table["mdc_delta_log1p_cpm_diseased_minus_healthy"] = table["mdc_delta_log1p_cpm_diseased_minus_healthy"].map(lambda x: fmt_float(x, 4))
    table["mdc_p_log1p_cpm_mannwhitney"] = table["mdc_p_log1p_cpm_mannwhitney"].map(lambda x: fmt_float(x, 3))
    for col in ["max_open_targets_overall", "max_open_targets_genetic", "max_open_targets_clinical", "max_open_targets_literature"]:
        table[col] = table[col].map(lambda x: fmt_float(x, 3))
    for col in ["pubmed_edges_checked", "pubmed_title_matches", "pubmed_target_terms_found"]:
        table[col] = table[col].fillna(0).astype(int)
    keep = [
        "gene_symbol",
        "display_evidence_layer",
        "display_bulk_role",
        "evidence_tier",
        "bulk_module_role",
        "relationship_directions",
        "pubmed_ids",
        "b_delta_log1p_cpm_diseased_minus_healthy",
        "b_p_log1p_cpm_mannwhitney",
        "mdc_delta_log1p_cpm_diseased_minus_healthy",
        "mdc_p_log1p_cpm_mannwhitney",
        "open_targets_diseases",
        "max_open_targets_overall",
        "max_open_targets_genetic",
        "max_open_targets_clinical",
        "pubmed_edges_checked",
        "pubmed_title_matches",
        "pubmed_target_terms_found",
        "fulltext_evidence_levels",
        "recommended_figure_use",
        "interpretation_boundary",
    ]
    return table[keep]


def write_legend_doc() -> None:
    out_path = DOCS_DIR / "MANUSCRIPT_FIGURE_LEGENDS_AND_STAT_NOTES.md"
    if out_path.exists():
        print(f"Preserved curated {out_path.relative_to(ROOT)}")
        return
    text = """# Draft Figure Legends And Statistical Notes

Last updated: 2026-07-01

These legends are manuscript-facing drafts. They avoid claiming therapeutic efficacy or direct human compound-target validation.

## Figure 1. Disease-first public-data triangulation and evidence boundaries

	(a) Disease-first synthesis used to move from public UC/IBD mucosal disease signatures to rectal single-cell localization, TCM candidate prioritization, orthogonal target support, and explicit interpretive limits. (b) Candidate-role convergence showing that Curcumin is the single primary B-cell-axis candidate rather than one of many equivalent target-overlap hits. (c) Claim map separating supported wording from claims not supported by the current public-data evidence.

Statistical notes: Figure 1 summarizes results from the downstream quantitative analyses and evidence verification. The figure defines the manuscript logic and claim hierarchy; it does not add an independent efficacy test.

Source data: `figures/source_data/fig1a_evidence_gate_summary.tsv`, `figures/source_data/fig1b_candidate_role_summary.tsv`, and `figures/source_data/fig1c_claim_boundary_summary.tsv`.

## Figure 2. Replicated UC/IBD mucosal disease signature

(a) GEO cohort inventory showing the discovery/sensitivity evidence group and the independent active UC validation group. GSE75214 and GSE59071 are retained within the same GPL6244 evidence group rather than treated as independent validation. (b) Gene-level disease-versus-control volcano plots for the discovery and independent validation groups. Points indicate genes; red and blue mark genes passing FDR < 0.05 with absolute log2 fold change at least 1. (c) Consensus mucosal signature genes ranked by mean log2 fold change across the two independent evidence groups.

Statistical notes: bulk transcriptomic contrasts were computed from processed GEO matrices. Gene-level P values and FDR values are reported in the corresponding source-data files. The downstream consensus signature is used as a disease-first anchor and should not be interpreted as a treatment mechanism.

Source data: `figures/source_data/fig2a_geo_group_audit.tsv`, `figures/source_data/fig2b_volcano_source.tsv`, and `figures/source_data/fig2c_consensus_signature.tsv`.

## Figure 3. Rectal B cells show the strongest scRNA disease-axis localization

(a) Disease-axis score differences between diseased and healthy rectal pseudobulk samples in GSE125527. B cells show the largest positive shift. (b) Patient-label null intervals for each cell type, with the observed disease-axis delta shown as a point. The B-cell shift exceeds the patient-label null expectation, whereas M/DC remains secondary. (c) Cell and patient support for each rectal pseudobulk cell-type contrast.

Statistical notes: cell-level counts were aggregated to patient/sample pseudobulk before disease contrasts. Panel a uses Mann-Whitney tests across pseudobulk samples. Panel b uses 10,000 patient-label permutations with a fixed seed. The B-cell empirical null P value is 0.00229977; the M/DC empirical null P value is 0.20438.

Source data: `figures/source_data/fig3a_celltype_disease_axis.tsv`, `figures/source_data/fig3b_patient_label_null.tsv`, and `figures/source_data/fig3c_pseudobulk_counts.tsv`.

## Figure 4. Curcumin target evidence bridges the disease signature to the rectal B-cell axis

(a) Candidate cell-state linkage ranking from disease-context HERB evidence and rectal scRNA localization. Curcumin is the leading B-cell-axis candidate. (b) Curcumin target evidence layers separating B-cell mechanism candidates, cytokine-context targets, and supplementary full-text-supported targets. (c) Rectal B-cell and M/DC pseudobulk expression deltas for selected Curcumin-linked targets.

Statistical notes: cell-state linkage scores are prioritization scores and do not measure therapeutic efficacy. BCL6, BLNK, and SYK are presented as full-text-supported B-cell mechanism candidates from PMID 36353208. CCL2, IL33, IL1B, and TNF are presented as cytokine-context targets supported by PMID 36196887 and should not be described as direct B-cell-intrinsic target modulation.

Source data: `figures/source_data/fig4a_candidate_bridge_ranking.tsv`, `figures/source_data/fig4b_curcumin_target_evidence_tiers.tsv`, and `figures/source_data/fig4c_curcumin_target_cell_contrast.tsv`.

## Figure 5. Orthogonal evidence supports target plausibility while defining boundary conditions

(a) Open Targets support for Curcumin disease-context targets, including genetic and clinical evidence subsets. (b) PubMed and full-text verification matrix for high-value Curcumin references. (c) ETCM2 mapped target overlap check, showing that Curcumin has mapped ETCM2 genes but no accepted overlap with the HERB disease-context target set. (d) L1000CDS2 exact-ID candidate screen, showing no top-result hit among covered candidate perturbagens.

Statistical notes: Open Targets scores support target-disease plausibility and do not establish compound action. PubMed/PMC and full-text verification supports bibliographic and target-context traceability. ETCM2 and LINCS define resource-specific limits; the ETCM2 non-overlap and L1000CDS2 non-hit should not be interpreted as proof that Curcumin lacks biological activity.

Source data: `figures/source_data/fig5a_open_targets_curcumin_support.tsv`, `figures/source_data/fig5b_pubmed_fulltext_verification.tsv`, `figures/source_data/fig5c_etcm2_overlap_boundary.tsv`, `figures/source_data/fig5d_lincs_candidate_reversal_status.tsv`, and `figures/source_data/fig5d_lincs_query_summary.tsv`.
"""
    out_path.write_text(text, encoding="utf-8")


def write_table_doc(table1: pd.DataFrame, table2: pd.DataFrame, table3: pd.DataFrame) -> None:
    out_path = DOCS_DIR / "MANUSCRIPT_TABLE_PACKAGE.md"
    if out_path.exists():
        print(f"Preserved curated {out_path.relative_to(ROOT)}")
        return
    text = f"""# Manuscript Table Package

Last updated: 2026-07-01

## Table 1. Public datasets and evidence resources

Source files: `tables/source_data/table1_dataset_and_resource_inventory.tsv`; manuscript-ready version `tables/manuscript/table1_dataset_and_resource_inventory.tsv`

Rows: {len(table1)}

Purpose: document the public data and evidence resources used in the manuscript, their role in the evidence chain, and the limitation attached to each resource.

## Table 2. Prioritized TCM-related candidates and evidence boundaries

Source files: `tables/source_data/table2_prioritized_tcm_candidates.tsv`; manuscript-ready version `tables/manuscript/table2_prioritized_tcm_candidates.tsv`

Rows: {len(table2)}

Purpose: summarize the candidate-ranking layer while making target-count bias, LINCS coverage, and B-cell/myeloid evidence boundaries visible.

## Table 3. Curcumin disease-context target evidence

Source files: `tables/source_data/table3_curcumin_target_evidence.tsv`; manuscript-ready version `tables/manuscript/table3_curcumin_target_evidence.tsv`

Rows: {len(table3)}

Purpose: provide gene-by-gene support for the Curcumin target bridge, integrating HERB evidence tier, rectal pseudobulk expression, Open Targets disease plausibility, PubMed verification, and full-text evidence tier.

## Submission Boundary

These tables support prioritization and evidence traceability. They do not claim validated therapeutic efficacy, direct compound-target engagement in human UC tissue, or prospective clinical utility.
"""
    out_path.write_text(text, encoding="utf-8")


def main() -> None:
    TABLE_DIR.mkdir(parents=True, exist_ok=True)
    MANUSCRIPT_TABLE_DIR.mkdir(parents=True, exist_ok=True)
    table1 = build_table_1()
    table2 = build_table_2()
    table3 = build_table_3()
    write_tsv(table1, "tables/source_data/table1_dataset_and_resource_inventory.tsv")
    write_tsv(table2, "tables/source_data/table2_prioritized_tcm_candidates.tsv")
    write_tsv(table3, "tables/source_data/table3_curcumin_target_evidence.tsv")
    write_curated_tsv(table1, "tables/manuscript/table1_dataset_and_resource_inventory.tsv")
    write_curated_tsv(build_manuscript_table_2(table2), "tables/manuscript/table2_prioritized_tcm_candidates.tsv")
    write_curated_tsv(build_manuscript_table_3(table3), "tables/manuscript/table3_curcumin_target_evidence.tsv")
    if os.environ.get("PUBLIC_REPRO_BUILD") != "1":
        write_legend_doc()
        write_table_doc(table1, table2, table3)
    print("Wrote manuscript table source data.")


if __name__ == "__main__":
    main()
