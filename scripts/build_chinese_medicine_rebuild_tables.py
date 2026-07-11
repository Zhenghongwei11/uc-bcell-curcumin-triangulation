#!/usr/bin/env python3
"""Build Chinese Medicine main tables from existing reproducible outputs."""

from __future__ import annotations

from pathlib import Path
import math
import re

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]

PUBCHEM_IDENTITY_OVERRIDES = {
    "969516": {
        "cas": "458-37-7",
        "note": "CAS corrected from HERB CAS field using PubChem CID 969516 synonyms.",
    },
    "3220": {
        "cas": "518-82-1",
        "note": "CAS added from PubChem.",
    },
    "479503": {
        "display_name": "Shikonin / Isoarnebin 4",
        "cas": "517-89-5",
        "note": "HERB name 'Isoarnebin 4' maps to PubChem CID 479503, whose primary synonym is shikonin; alkannin is the stereochemical counterpart with a different PubChem CID/InChIKey.",
    },
    "1548910": {
        "cas": "61434-67-1",
        "note": "CAS corrected for cis-resveratrol using PubChem CID 1548910; 501-36-0 corresponds to trans-resveratrol/resveratrol.",
    },
    "445154": {
        "cas": "501-36-0",
        "note": "CAS added/verified for trans-resveratrol using PubChem CID 445154.",
    },
}


def read_tsv(path: str) -> pd.DataFrame:
    return pd.read_csv(ROOT / path, sep="\t")


def write_tsv(df: pd.DataFrame, path: str) -> None:
    out = ROOT / path
    out.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out, sep="\t", index=False)


def clean_join(values) -> str:
    vals = []
    for value in values:
        if pd.isna(value):
            continue
        text = str(value).strip()
        if not text or text.lower() == "nan":
            continue
        vals.append(text)
    return "; ".join(sorted(set(vals)))


def split_genes(value) -> set[str]:
    if pd.isna(value):
        return set()
    text = str(value).strip()
    if not text or text.lower() == "nan":
        return set()
    return {x.strip() for x in text.split(";") if x.strip()}


def fmt_int(value) -> int:
    if pd.isna(value):
        return 0
    return int(round(float(value)))


def display_compound_name(source_name: str, pubchem_id) -> str:
    name = str(source_name).strip()
    pubchem_key = ""
    if not pd.isna(pubchem_id):
        pubchem_key = str(int(float(pubchem_id)))
    if pubchem_key in PUBCHEM_IDENTITY_OVERRIDES and "display_name" in PUBCHEM_IDENTITY_OVERRIDES[pubchem_key]:
        return PUBCHEM_IDENTITY_OVERRIDES[pubchem_key]["display_name"]
    if name == "Berberime":
        return "Berberine"
    if name == "Archin" and pubchem_key == "3220":
        return "Emodin"
    return name


def bh_fdr(pvalues: list[float]) -> list[float]:
    n = len(pvalues)
    order = sorted(range(n), key=lambda i: pvalues[i])
    adjusted = [math.nan] * n
    prev = 1.0
    for rank_from_end, i in enumerate(reversed(order), start=1):
        rank = n - rank_from_end + 1
        value = min(prev, pvalues[i] * n / rank)
        prev = value
        adjusted[i] = min(value, 1.0)
    return adjusted


def build_compound_tables() -> None:
    metadata = read_tsv("data/derived/herb_uc_ibd_ingredient_candidates.tsv")
    bridge = read_tsv("data/derived/compound_cell_context_summary.tsv")
    specificity = read_tsv("data/derived/compound_specificity_ranking.tsv")
    overlap = read_tsv("data/derived/standard_overlap_compound_ranking.tsv")
    null_summary = read_tsv("data/derived/curcumin_matched_target_null.tsv")

    base = bridge.merge(
        metadata[
            [
                "ingredient_id",
                "ingredient_name",
                "evidence_count",
                "reference_count",
                "clinical_trial_count",
                "meta_analysis_count",
                "pubchem_id",
                "inchikey",
                "cas_id",
                "identifier_status",
            ]
        ],
        on=["ingredient_id", "ingredient_name"],
        how="left",
    ).merge(
        specificity[
            [
                "ingredient_id",
                "target_set_specificity_rank",
                "target_set_specificity_score_without_direction_alignment",
            ]
        ],
        on="ingredient_id",
        how="left",
    ).merge(
        overlap[
            [
                "ingredient_id",
                "standard_overlap_rank",
                "standard_overlap_hits",
                "standard_overlap_up_hits",
                "standard_overlap_genes",
            ]
        ],
        on="ingredient_id",
        how="left",
    )

    base["display_name"] = [
        display_compound_name(n, p) for n, p in zip(base["ingredient_name"], base["pubchem_id"])
    ]
    base["canonical_key"] = base.apply(
        lambda row: str(row["inchikey"]).strip()
        if not pd.isna(row["inchikey"]) and str(row["inchikey"]).strip()
        else (
            f"PUBCHEM:{int(float(row['pubchem_id']))}"
            if not pd.isna(row["pubchem_id"])
            else re.sub(r"[^A-Za-z0-9]+", "_", str(row["display_name"]).lower()).strip("_")
        ),
        axis=1,
    )

    rows = []
    full_rows = []
    for _, group in base.groupby("canonical_key", sort=False):
        first = group.iloc[0]
        disease_targets = set()
        bulk_up = set()
        b_expr = set()
        b_inc = set()
        mdc_inc = set()
        overlap_genes = set()
        for _, row in group.iterrows():
            disease_targets |= split_genes(row.get("disease_context_unique_targets"))
            bulk_up |= split_genes(row.get("bulk_up_target_hits"))
            b_expr |= split_genes(row.get("b_cell_expressed_targets"))
            b_inc |= split_genes(row.get("b_cell_disease_increased_targets"))
            mdc_inc |= split_genes(row.get("mdc_disease_increased_targets"))
            overlap_genes |= split_genes(row.get("standard_overlap_genes"))

        duplicate_status = (
            f"merged_{len(group)}_source_records"
            if len(group) > 1
            else "single_source_record"
        )
        display_name = clean_join(group["display_name"]).split("; ")[0]
        if "Curcumin" in set(group["display_name"]):
            display_name = "Curcumin"
        elif "Berberine" in set(group["display_name"]):
            display_name = "Berberine"
        elif "Emodin" in set(group["display_name"]):
            display_name = "Emodin"

        source_names = clean_join(group["ingredient_name"])
        pubchem_ids = [
            str(int(float(x))) for x in group["pubchem_id"] if not pd.isna(x)
        ]
        primary_pubchem = sorted(set(pubchem_ids))[0] if pubchem_ids else ""
        source_note = ""
        if "Berberime" in set(group["ingredient_name"]):
            source_note = "HERB source spelling 'Berberime' standardized to Berberine."
        if "Archin" in set(group["ingredient_name"]):
            source_note = "HERB source name 'Archin' standardized to Emodin using PubChem CID 3220/InChIKey"
        cas_value = clean_join(group["cas_id"])
        if primary_pubchem in PUBCHEM_IDENTITY_OVERRIDES:
            override = PUBCHEM_IDENTITY_OVERRIDES[primary_pubchem]
            cas_value = override.get("cas", cas_value)
            note = override.get("note", "")
            source_note = "; ".join(x for x in [source_note, note] if x)

        row = {
            "Compound": display_name,
            "HERB source names": source_names,
            "HERB IDs": clean_join(group["ingredient_id"]),
            "PubChem CID": clean_join(pubchem_ids),
            "InChIKey": clean_join(group["inchikey"]),
            "CAS": cas_value,
            "Disease-context targets": len(disease_targets),
            "Bulk up-module targets": len(bulk_up),
            "Rectal B-lineage expressed targets": len(b_expr),
            "Rectal B-lineage increased targets": len(b_inc),
            "M/DC increased targets": len(mdc_inc),
            "HERB evidence records": fmt_int(group["evidence_count"].max()),
            "HERB references": fmt_int(group["reference_count"].max()),
            "Clinical trial records": fmt_int(group["clinical_trial_count"].max()),
            "Meta-analysis records": fmt_int(group["meta_analysis_count"].max()),
            "Target-set score": round(float(group["target_set_specificity_score_without_direction_alignment"].max()), 2)
            if not group["target_set_specificity_score_without_direction_alignment"].isna().all()
            else math.nan,
            "Target-set rank": fmt_int(group["target_set_specificity_rank"].min()),
            "Overlap baseline rank": fmt_int(group["standard_overlap_rank"].min()),
            "Overlap baseline genes": ";".join(sorted(overlap_genes)),
            "Duplicate status": duplicate_status,
            "Source note": source_note,
        }
        full_rows.append(row)

    full = pd.DataFrame(full_rows)
    full = full.sort_values(
        ["Rectal B-lineage increased targets", "Disease-context targets", "HERB evidence records"],
        ascending=[False, False, False],
    ).reset_index(drop=True)
    full.insert(0, "Table order", range(1, len(full) + 1))
    herb_source_path = ROOT / "tables/supplementary/chinese_medicine_compound_to_herb_sources.tsv"
    if herb_source_path.exists():
        herb_source = pd.read_csv(herb_source_path, sep="\t")
        herb_source["Representative HERB source"] = herb_source.apply(
            lambda row: f"{row['Representative herb Chinese']} ({row['Representative herb pinyin']}; {row['Representative herb Latin']})"
            if not pd.isna(row.get("Representative herb Chinese")) and str(row.get("Representative herb Chinese")).strip()
            else "",
            axis=1,
        )
        full = full.merge(
            herb_source[["Compound", "Representative HERB source"]],
            on="Compound",
            how="left",
        )
    if "Representative HERB source" not in full.columns:
        full["Representative HERB source"] = ""
    full = full.rename(columns={"Representative HERB source": "HERB-linked herbal record"})

    main_cols = [
        "Table order",
        "Compound",
        "HERB-linked herbal record",
        "PubChem CID",
        "Disease-context targets",
        "Rectal B-lineage expressed targets",
        "Rectal B-lineage increased targets",
        "M/DC increased targets",
        "HERB evidence records",
    ]
    main = full.loc[:, main_cols].head(12)

    ledger = full.loc[
        :,
        [
            "Table order",
            "Compound",
            "HERB source names",
            "HERB IDs",
            "PubChem CID",
            "InChIKey",
            "CAS",
            "Duplicate status",
            "Source note",
        ],
    ]

    null_summary = null_summary.rename(
        columns={
            "analysis": "Analysis",
            "result": "Result",
            "interpretation": "Interpretation",
        }
    )

    write_tsv(main, "tables/main/table2_chinese_medicine_compound_mapping.tsv")
    write_tsv(full, "tables/supplementary/chinese_medicine_compound_mapping_full.tsv")
    write_tsv(ledger, "tables/supplementary/chinese_medicine_compound_ledger.tsv")
    write_tsv(null_summary, "tables/supplementary/chinese_medicine_candidate_specificity_context.tsv")


def build_dataset_table() -> None:
    rows = [
        {
            "Accession/resource": "GSE75214/GSE59071",
            "Tissue/data type": "Colonic mucosal bulk transcriptomics",
            "Platform": "GPL6244 microarray",
            "Age/tissue context": "Adult or mixed-age mucosal biopsies; colon",
            "Clinical groups or records": "Active UC colon and normal control colon mucosa",
            "Source total": "194 samples in GSE75214; 116-sample GSE59071 colon subset",
            "Included in primary contrast": "74 active disease vs 11 controls",
            "Analytical role": "Shared bulk evidence group",
            "Query/access date": "2026-06-29",
            "Source note": "GSE59071 selected contrast samples are contained in GSE75214 and are counted once.",
        },
        {
            "Accession/resource": "GSE87466",
            "Tissue/data type": "Colonic mucosal bulk transcriptomics",
            "Platform": "GPL13158 microarray",
            "Age/tissue context": "Colonic mucosal biopsies; colon",
            "Clinical groups or records": "Active UC and normal control mucosa",
            "Source total": "108 samples",
            "Included in primary contrast": "87 active UC vs 21 controls",
            "Analytical role": "Independent bulk evidence group",
            "Query/access date": "2026-06-29",
            "Source note": "Independent platform bulk evidence group.",
        },
        {
            "Accession/resource": "GSE125527",
            "Tissue/data type": "Rectal single-cell RNA-seq",
            "Platform": "10x Genomics scRNA-seq",
            "Age/tissue context": "Pediatric rectal biopsies",
            "Clinical groups or records": "Pediatric UC and healthy controls",
            "Source total": "103 samples/captures in source record",
            "Included in primary contrast": "Rectal B-lineage pseudobulk: 7 diseased vs 4 healthy patients",
            "Analytical role": "Primary cell-state context analysis",
            "Query/access date": "2026-06-29",
            "Source note": "Used for patient-level rectal cell-state context analysis.",
        },
        {
            "Accession/resource": "GSE182270",
            "Tissue/data type": "B-lineage-focused single-cell RNA-seq",
            "Platform": "scRNA-seq",
            "Age/tissue context": "Colonic or mucosal single-cell samples",
            "Clinical groups or records": "UC inflamed mucosa and healthy/noninflamed controls",
            "Source total": "9 samples",
            "Included in primary contrast": "5 UC inflamed vs 4 controls",
            "Analytical role": "Exploratory external directional comparison",
            "Query/access date": "2026-07-02",
            "Source note": "Directionally contextual only; not treated as validation.",
        },
        {
            "Accession/resource": "HERB 2.0",
            "Tissue/data type": "Chinese medicine compound-target knowledge base",
            "Platform": "Curated and literature-mined database",
            "Age/tissue context": "Not applicable",
            "Clinical groups or records": "TCM-related ingredient, target, disease, and reference records",
            "Source total": "Queried for UC/IBD-related disease-context compounds",
            "Included in primary contrast": "Deduplicated disease-context compound and target records",
            "Analytical role": "Compound and target annotation",
            "Query/access date": "2026-06-27 to 2026-07-02",
            "Source note": "Knowledge-base annotations; not experimental validation.",
        },
        {
            "Accession/resource": "ETCM2",
            "Tissue/data type": "Chinese medicine compound-target knowledge base",
            "Platform": "Web database",
            "Age/tissue context": "Not applicable",
            "Clinical groups or records": "Ingredient-target annotations",
            "Source total": "Queried for selected candidate compounds",
            "Included in primary contrast": "Mapped human gene symbols after target-name standardization",
            "Analytical role": "External database context check",
            "Query/access date": "2026-06-27 to 2026-07-02",
            "Source note": "Target names require mapping; unmapped or ambiguous records are not over-interpreted.",
        },
        {
            "Accession/resource": "Open Targets Platform",
            "Tissue/data type": "Target-disease evidence platform",
            "Platform": "Public target-disease database",
            "Age/tissue context": "Not applicable",
            "Clinical groups or records": "UC, Crohn's disease, and IBD target-disease evidence",
            "Source total": "Queried for disease-context target genes",
            "Included in primary contrast": "Curcumin-associated disease-context targets",
            "Analytical role": "Target-disease plausibility check",
            "Query/access date": "2026-07-01",
            "Source note": "Supports target-disease plausibility, not compound action.",
        },
        {
            "Accession/resource": "L1000CDS2/LINCS",
            "Tissue/data type": "Perturbational transcriptomic signatures",
            "Platform": "L1000CDS2 query against LINCS signatures",
            "Age/tissue context": "Not applicable",
            "Clinical groups or records": "Compound perturbation signatures",
            "Source total": "Exact-identifier query for covered candidates",
            "Included in primary contrast": "Top-result reversal query status",
            "Analytical role": "Perturbational context check",
            "Query/access date": "2026-06-30",
            "Source note": "Exact-identifier top-result screen; non-hit does not rule out broader LINCS activity.",
        },
    ]
    df = pd.DataFrame(rows)
    write_tsv(df, "tables/main/table1_chinese_medicine_resources.tsv")


def build_curcumin_gene_tables() -> None:
    cur = read_tsv("data/derived/curcumin_blineage_target_evidence.tsv")
    cur["p_for_bh_across_22"] = cur["b_p_log1p_cpm_mannwhitney"].fillna(1.0).astype(float)
    cur["b_fdr_bh_across_22"] = bh_fdr(cur["p_for_bh_across_22"].tolist())
    cur["tested_in_rectal_b_lineage"] = cur["b_p_log1p_cpm_mannwhitney"].notna().map(
        {True: "yes", False: "no"}
    )

    evidence_label = {
        "A_module_aligned_b_cell_increased": "Bulk-aligned and B-lineage increased",
        "B_module_aligned": "Bulk-aligned cytokine-context target",
        "C_b_cell_mechanism_increased": "Selected B-cell-context gene",
        "D_b_cell_increased": "B-lineage increased context gene",
        "F_context_target_only": "Literature context target",
    }
    bulk_label = {
        "bulk_up_top150": "Consensus up-module",
        "bulk_down_top150": "Consensus down-module",
        "not_bulk_top150": "Not in consensus top module",
    }
    relation_label = {
        "suppress/inhibitor": "suppress/inhibit",
        "upregulate/increase": "upregulate/increase",
        "downregulate/decrease": "downregulate/decrease",
        "downregulate/decrease;suppress/inhibitor": "downregulate/suppress",
    }
    out = pd.DataFrame(
        {
            "Gene": cur["gene_symbol"],
            "Evidence context": cur["evidence_tier"].map(evidence_label).fillna(cur["evidence_tier"].str.replace("_", " ", regex=False)),
            "Bulk module role": cur["bulk_module_role"].map(bulk_label).fillna(cur["bulk_module_role"].str.replace("_", " ", regex=False)),
            "Reported curcumin relation": cur["relationships"].map(relation_label).fillna(cur["relationships"].fillna("not resolved")),
            "Reported direction": cur["relationship_directions"].fillna("not resolved"),
            "Rectal B-lineage delta": cur["b_delta_log1p_cpm_diseased_minus_healthy"],
            "Rectal B-lineage P value": cur["b_p_log1p_cpm_mannwhitney"],
            "Rectal B-lineage FDR across 22 genes": cur["b_fdr_bh_across_22"],
            "M/DC delta": cur["mdc_delta_log1p_cpm_diseased_minus_healthy"],
            "M/DC P value": cur["mdc_p_log1p_cpm_mannwhitney"],
            "Tested in rectal B-lineage": cur["tested_in_rectal_b_lineage"],
            "PubMed IDs": cur["pubmed_ids"],
            "Reference context": cur["reference_titles"],
        }
    )
    out.loc[out["Gene"] == "SH3KBP1", "Evidence context"] = "Additional literature-context gene"

    main = out[out["Gene"].isin(["BCL6", "BLNK", "SYK"])].copy()
    main = main.sort_values(
        ["Rectal B-lineage delta", "Gene"], ascending=[False, True]
    ).reset_index(drop=True)

    write_tsv(main, "tables/main/table3_curcumin_gene_followup.tsv")
    write_tsv(out, "tables/supplementary/curcumin_gene_followup_full.tsv")


def build_cross_context_direction_table() -> None:
    src = ROOT / "tables/main/table4_therapeutic_direction_consistency.tsv"
    if not src.exists():
        return
    df = pd.read_csv(src, sep="\t")
    rename = {
        "gene": "Gene",
        "evidence role": "Evidence context",
        "disease direction": "Disease direction in analyzed UC data",
        "rectal B delta": "Rectal B-lineage delta",
        "curcumin-reported regulation": "Reported curcumin regulation",
        "curcumin_reported_regulation_direction": "Reported curcumin regulation",
        "disease_direction_in_rectal_B_cells": "Disease direction in analyzed UC data",
        "directional_consistency": "Cross-context observation",
        "directional consistency": "Cross-context observation",
        "evidence_level": "Evidence level",
    }
    df = df.rename(columns={k: v for k, v in rename.items() if k in df.columns})
    for col in ["Previously stated expected direction", "expected therapeutic direction"]:
        if col in df.columns:
            df = df.drop(columns=[col])
    if "Cross-context observation" in df.columns:
        df["Cross-context observation"] = df["Cross-context observation"].replace(
            {
                "Directionally consistent": "Directionally compatible across heterogeneous contexts",
                "Contextual, not normalizing": "Reported regulation does not indicate simple transcriptional normalization",
            }
        )
    df.insert(
        0,
        "Interpretation note",
        "Cross-context literature comparison only; not evidence of therapeutic normalization.",
    )
    write_tsv(df, "tables/supplementary/cross_context_reported_direction_comparison.tsv")


def main() -> None:
    build_dataset_table()
    build_compound_tables()
    build_curcumin_gene_tables()
    build_cross_context_direction_table()


if __name__ == "__main__":
    main()
