#!/usr/bin/env python3
"""Sync manuscript supplementary TSV tables from regenerated analysis outputs."""

from __future__ import annotations

import shutil
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SUPP = ROOT / "tables" / "supplementary"

COPIES = {
    "tables/manuscript/table3_curcumin_target_evidence.tsv": "Table_S1_full_curcumin_target_evidence.tsv",
    "figures/source_data/fig1a_evidence_gate_summary.tsv": "fig1a_evidence_gate_summary.tsv",
    "figures/source_data/fig1b_candidate_role_summary.tsv": "fig1b_candidate_role_summary.tsv",
    "figures/source_data/fig1c_claim_boundary_summary.tsv": "fig1c_claim_boundary_summary.tsv",
    "figures/source_data/fig2a_geo_group_audit.tsv": "fig2a_geo_group_audit.tsv",
    "figures/source_data/fig2b_volcano_source.tsv": "fig2b_volcano_source.tsv",
    "figures/source_data/fig2c_consensus_signature.tsv": "fig2c_consensus_signature.tsv",
    "figures/source_data/fig3_unclassified_cell_qc.tsv": "fig3_unclassified_cell_qc.tsv",
    "figures/source_data/fig3a_celltype_disease_axis.tsv": "fig3a_celltype_disease_axis.tsv",
    "figures/source_data/fig3b_patient_label_null.tsv": "fig3b_patient_label_null.tsv",
    "figures/source_data/fig3c_pseudobulk_counts.tsv": "fig3c_pseudobulk_counts.tsv",
    "figures/source_data/fig4a_candidate_bridge_ranking.tsv": "fig4a_candidate_bridge_ranking.tsv",
    "figures/source_data/fig4b_curcumin_target_evidence_tiers.tsv": "fig4b_curcumin_target_evidence_tiers.tsv",
    "figures/source_data/fig4c_curcumin_target_cell_contrast.tsv": "fig4c_curcumin_target_cell_contrast.tsv",
    "figures/source_data/fig5a_open_targets_curcumin_support.tsv": "fig5a_open_targets_curcumin_support.tsv",
    "figures/source_data/fig5b_pubmed_fulltext_verification.tsv": "fig5b_pubmed_fulltext_verification.tsv",
    "figures/source_data/fig5c_etcm2_overlap_boundary.tsv": "fig5c_etcm2_overlap_boundary.tsv",
    "figures/source_data/fig5d_lincs_candidate_reversal_status.tsv": "fig5d_lincs_candidate_reversal_status.tsv",
    "figures/source_data/fig5d_lincs_query_summary.tsv": "fig5d_lincs_query_summary.tsv",
    "figures/source_data/fig_statistical_ci_supplement.tsv": "fig_statistical_ci_supplement.tsv",
    "tables/supplementary/statistical_ci_supplement.tsv": "statistical_ci_supplement.tsv",
    "results/m17_no_download_robustness/bulk_b_plasma_marker_gene_directions.tsv": "bulk_b_plasma_marker_gene_directions.tsv",
    "results/m17_no_download_robustness/bulk_b_plasma_marker_summary.tsv": "bulk_b_plasma_marker_summary.tsv",
    "results/m17_no_download_robustness/herb_candidate_specificity_ranking.tsv": "herb_candidate_specificity_ranking.tsv",
    "results/m17_no_download_robustness/herb_decoy_null_summary.tsv": "herb_decoy_null_summary.tsv",
    "results/m17_no_download_robustness/scrna_signature_sensitivity_bcell_summary.tsv": "scrna_signature_sensitivity_bcell_summary.tsv",
    "results/m17_no_download_robustness/scrna_signature_sensitivity_by_celltype.tsv": "scrna_signature_sensitivity_by_celltype.tsv",
    "results/m17_no_download_robustness/scrna_signature_variant_gene_coverage.tsv": "scrna_signature_variant_gene_coverage.tsv",
    "results/m18_reviewer_lightweight_extensions/bulk_scrna_reference_deconvolution_proxy_by_sample.tsv": "bulk_scrna_reference_deconvolution_proxy_by_sample.tsv",
    "results/m18_reviewer_lightweight_extensions/bulk_scrna_reference_deconvolution_proxy_summary.tsv": "bulk_scrna_reference_deconvolution_proxy_summary.tsv",
    "results/m18_reviewer_lightweight_extensions/scrna_reference_marker_genes.tsv": "scrna_reference_marker_genes.tsv",
    "results/m18_reviewer_lightweight_extensions/standard_overlap_baseline_candidate_ranking.tsv": "standard_overlap_baseline_candidate_ranking.tsv",
    "results/m18_reviewer_lightweight_extensions/standard_overlap_baseline_summary.tsv": "standard_overlap_baseline_summary.tsv",
    "results/m19_gse182270_bcell_replication/gse182270_bcell_cell_module_scores_compact.tsv": "gse182270_bcell_cell_module_scores_compact.tsv",
    "results/m19_gse182270_bcell_replication/gse182270_bcell_group_comparison.tsv": "gse182270_bcell_group_comparison.tsv",
    "results/m19_gse182270_bcell_replication/gse182270_bcell_sample_module_scores.tsv": "gse182270_bcell_sample_module_scores.tsv",
    "results/m19_gse182270_bcell_replication/gse182270_sample_metadata_interpretation.tsv": "gse182270_sample_metadata_interpretation.tsv",
}


def main() -> int:
    SUPP.mkdir(parents=True, exist_ok=True)
    missing: list[str] = []
    for source, dest_name in COPIES.items():
        source_path = ROOT / source
        if not source_path.exists():
            missing.append(source)
            continue
        dest_path = SUPP / dest_name
        if source_path.resolve() != dest_path.resolve():
            shutil.copy2(source_path, dest_path)
    if missing:
        raise SystemExit("Missing supplementary source tables:\n" + "\n".join(missing))
    print(f"Synchronized {len(COPIES)} supplementary TSV tables.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
