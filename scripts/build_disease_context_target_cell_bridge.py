#!/usr/bin/env python3
"""Build a disease-context HERB target to scRNA cell-axis bridge."""

from __future__ import annotations

import argparse
import csv
import gzip
import math
import re
from collections import defaultdict
from pathlib import Path
from statistics import mean
from typing import Any, Dict, Iterable, List, Sequence, Set, Tuple

from scipy.stats import mannwhitneyu

from analyze_gse125527_scrna_modules import (
    load_cluster_maps,
    load_metadata_map,
    load_patient_id_map,
    parse_sample_from_filename,
)


OUT_DIR = Path("results/m9_target_cell_bridge")
DOC_PATH = Path("docs/M9_DISEASE_CONTEXT_TARGET_CELL_BRIDGE_RUN.md")

DISEASE_PATTERNS = {
    "colitis": r"\bcolitis\b",
    "ulcerative_colitis": r"ulcerative colitis|\buc\b",
    "inflammatory_bowel_disease": r"inflammatory bowel disease|\bibd\b",
    "crohn": r"crohn",
    "dss": r"dextran sulfate sodium|\bdss\b",
    "tnbs": r"\btnbs\b|trinitrobenzene",
    "intestinal_inflammation": r"intestinal inflammation|colonic inflammation|colon inflammation",
    "mucosal_barrier": r"mucosal barrier|intestinal barrier|colonic barrier",
    "colitis_associated": r"colitis-associated|colitis associated",
}

CORE_DISEASE_TERMS = {
    "colitis",
    "ulcerative_colitis",
    "inflammatory_bowel_disease",
    "crohn",
    "dss",
    "tnbs",
    "colitis_associated",
}

UP_RELATION_PATTERNS = [
    r"upregulate",
    r"increase",
    r"agonist",
    r"activator",
    r"activate",
    r"promote",
]

DOWN_RELATION_PATTERNS = [
    r"downregulate",
    r"decrease",
    r"suppress",
    r"inhibitor",
    r"inhibit",
    r"antagonist",
    r"block",
]


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


def read_gene_list(path: Path) -> Set[str]:
    return {line.strip().upper() for line in path.read_text(encoding="utf-8").splitlines() if line.strip()}


def disease_terms(text: str) -> List[str]:
    terms: List[str] = []
    for name, pattern in DISEASE_PATTERNS.items():
        if re.search(pattern, text, flags=re.IGNORECASE):
            terms.append(name)
    return terms


def relation_direction(relationship: str, supporting_sentence: str) -> str:
    rel = relationship.lower()
    if rel and rel != "na":
        rel_up = any(re.search(pattern, rel, flags=re.IGNORECASE) for pattern in UP_RELATION_PATTERNS)
        rel_down = any(re.search(pattern, rel, flags=re.IGNORECASE) for pattern in DOWN_RELATION_PATTERNS)
        if rel_up and not rel_down:
            return "up"
        if rel_down and not rel_up:
            return "down"
        if rel_up and rel_down:
            return "mixed"
        if "bind" in rel or "target" in rel:
            return "binding_or_target"

    text = supporting_sentence.lower()
    up = any(re.search(pattern, text, flags=re.IGNORECASE) for pattern in UP_RELATION_PATTERNS)
    down = any(re.search(pattern, text, flags=re.IGNORECASE) for pattern in DOWN_RELATION_PATTERNS)
    if up and not down:
        return "up"
    if down and not up:
        return "down"
    if up and down:
        return "mixed"
    if "bind" in text or "target" in text:
        return "binding_or_target"
    return "unknown"


def therapeutic_alignment(direction: str, in_up: bool, in_down: bool) -> str:
    if direction == "down" and in_up:
        return "aligned_suppresses_disease_up_target"
    if direction == "up" and in_down:
        return "aligned_restores_disease_down_target"
    if direction == "up" and in_up:
        return "opposes_upregulates_disease_up_target"
    if direction == "down" and in_down:
        return "opposes_suppresses_disease_down_target"
    if in_up or in_down:
        return "module_hit_direction_unclear"
    return "not_in_bulk_module"


def load_disease_context_edges(edge_path: Path, up_genes: Set[str], down_genes: Set[str]) -> List[Dict[str, Any]]:
    rows: List[Dict[str, Any]] = []
    for row in read_tsv(edge_path):
        if row["edge_source"] != "HERB_detail.drug_paper_target":
            continue
        text = " ".join([row.get("reference_title", ""), row.get("supporting_sentence", "")])
        terms = disease_terms(text)
        if not set(terms) & CORE_DISEASE_TERMS:
            continue
        gene = row["gene_symbol"].upper()
        direction = relation_direction(row.get("relationship", ""), row.get("supporting_sentence", ""))
        rows.append(
            {
                **row,
                "gene_symbol": gene,
                "disease_context_terms": ";".join(terms),
                "relationship_direction": direction,
                "in_bulk_uc_ibd_up150": "yes" if gene in up_genes else "no",
                "in_bulk_uc_ibd_down150": "yes" if gene in down_genes else "no",
                "therapeutic_alignment": therapeutic_alignment(direction, gene in up_genes, gene in down_genes),
            }
        )
    rows.sort(
        key=lambda row: (
            row["ingredient_name"].lower(),
            row["gene_symbol"],
            row.get("pubmed_id", ""),
            row.get("reference_id", ""),
        )
    )
    return rows


def safe_pvalue(a: Sequence[float], b: Sequence[float]) -> str:
    if len(a) < 2 or len(b) < 2:
        return ""
    try:
        return f"{mannwhitneyu(a, b, alternative='two-sided').pvalue:.6g}"
    except ValueError:
        return ""


def celltype_group_labels(celltype: str, cluster: str) -> List[str]:
    if cluster:
        return [celltype, f"{celltype}:{cluster}"]
    return [celltype]


def score_target_expression(
    sample_dir: Path,
    raw_dir: Path,
    target_genes: Set[str],
) -> Tuple[List[Dict[str, Any]], Set[str]]:
    metadata_map = load_metadata_map(raw_dir / "GSE125527_cell_metadata.csv.gz")
    cluster_map = load_cluster_maps(raw_dir)
    patient_id_map = load_patient_id_map(raw_dir)

    aggregate: Dict[Tuple[str, str, str, str, str], Dict[str, float]] = defaultdict(lambda: {"sum_counts": 0.0, "detected_cells": 0.0, "total_umi": 0.0, "n_cells": 0.0})
    present_genes: Set[str] = set()

    for path in sorted(sample_dir.glob("*_cell-gene_UMI_table.tsv.gz")):
        old_patient, tissue = parse_sample_from_filename(path)
        patient = patient_id_map.get(old_patient, old_patient)
        with gzip.open(path, "rt", encoding="utf-8", errors="replace", newline="") as handle:
            reader = csv.reader(handle, delimiter="\t")
            header = next(reader)
            genes = [gene.upper() for gene in header[1:]]
            target_indices = [(idx, gene) for idx, gene in enumerate(genes) if gene in target_genes]
            present_genes.update(gene for _, gene in target_indices)
            if not target_indices:
                continue
            for raw in reader:
                if not raw:
                    continue
                barcode = raw[0]
                counts = raw[1:]
                total = sum(float(value) for value in counts if value)
                if total <= 0:
                    continue
                meta = metadata_map.get((patient, tissue, barcode), {})
                disease = meta.get("disease_assignment") or ("diseased" if patient.startswith("U") else "healthy")
                celltype = meta.get("celltype", "unknown")
                cluster = cluster_map.get((patient, tissue, barcode), "")
                group_celltypes = celltype_group_labels(celltype, cluster)
                for idx, gene in target_indices:
                    value = float(counts[idx]) if idx < len(counts) and counts[idx] else 0.0
                    for group_celltype in group_celltypes:
                        key = (patient, tissue, disease, group_celltype, gene)
                        rec = aggregate[key]
                        rec["sum_counts"] += value
                        rec["detected_cells"] += 1.0 if value > 0 else 0.0
                        rec["total_umi"] += total
                        rec["n_cells"] += 1.0

    rows: List[Dict[str, Any]] = []
    for (patient, tissue, disease, celltype, gene), rec in sorted(aggregate.items()):
        n_cells = int(rec["n_cells"])
        total_umi = rec["total_umi"]
        normalized = math.log1p((rec["sum_counts"] / total_umi) * 10000.0) if total_umi > 0 else 0.0
        rows.append(
            {
                "patient_assignment": patient,
                "tissue_assignment": tissue,
                "disease_assignment": disease,
                "celltype": celltype,
                "gene_symbol": gene,
                "n_cells": n_cells,
                "sum_counts": f"{rec['sum_counts']:.6g}",
                "detected_cells": int(rec["detected_cells"]),
                "detection_fraction": f"{rec['detected_cells'] / n_cells:.6g}" if n_cells else "",
                "total_umi": f"{total_umi:.6g}",
                "log1p_cpm": f"{normalized:.6g}",
            }
        )
    return rows, present_genes


def contrast_gene_expression(rows: Sequence[Dict[str, Any]]) -> List[Dict[str, Any]]:
    by_group: Dict[Tuple[str, str, str], List[Dict[str, Any]]] = defaultdict(list)
    for row in rows:
        by_group[(row["tissue_assignment"], row["celltype"], row["gene_symbol"])].append(row)
    contrast_rows: List[Dict[str, Any]] = []
    for (tissue, celltype, gene), group_rows in sorted(by_group.items()):
        diseased = [row for row in group_rows if row["disease_assignment"] == "diseased"]
        healthy = [row for row in group_rows if row["disease_assignment"] == "healthy"]
        if not diseased or not healthy:
            continue
        d_expr = [float(row["log1p_cpm"]) for row in diseased]
        h_expr = [float(row["log1p_cpm"]) for row in healthy]
        d_det = [float(row["detection_fraction"]) for row in diseased]
        h_det = [float(row["detection_fraction"]) for row in healthy]
        contrast_rows.append(
            {
                "tissue_assignment": tissue,
                "celltype": celltype,
                "gene_symbol": gene,
                "n_diseased_samples": len(diseased),
                "n_healthy_samples": len(healthy),
                "diseased_cells": sum(int(row["n_cells"]) for row in diseased),
                "healthy_cells": sum(int(row["n_cells"]) for row in healthy),
                "mean_log1p_cpm_diseased": f"{mean(d_expr):.6g}",
                "mean_log1p_cpm_healthy": f"{mean(h_expr):.6g}",
                "delta_log1p_cpm_diseased_minus_healthy": f"{mean(d_expr) - mean(h_expr):.6g}",
                "p_log1p_cpm_mannwhitney": safe_pvalue(d_expr, h_expr),
                "mean_detection_diseased": f"{mean(d_det):.6g}",
                "mean_detection_healthy": f"{mean(h_det):.6g}",
                "delta_detection_diseased_minus_healthy": f"{mean(d_det) - mean(h_det):.6g}",
                "p_detection_mannwhitney": safe_pvalue(d_det, h_det),
            }
        )
    contrast_rows.sort(key=lambda row: (row["tissue_assignment"], row["celltype"], row["gene_symbol"]))
    return contrast_rows


def summarize_candidate_bridge(edges: Sequence[Dict[str, Any]], contrast_rows: Sequence[Dict[str, Any]]) -> List[Dict[str, Any]]:
    contrast = {(row["celltype"], row["gene_symbol"]): row for row in contrast_rows if row["tissue_assignment"] == "R"}
    by_candidate: Dict[str, Dict[str, Any]] = {}
    for edge in edges:
        key = edge["ingredient_id"]
        entry = by_candidate.setdefault(
            key,
            {
                "ingredient_id": edge["ingredient_id"],
                "ingredient_name": edge["ingredient_name"],
                "candidate_priority_class": edge["candidate_priority_class"],
                "candidate_priority_score": edge["candidate_priority_score"],
                "disease_context_edge_rows": 0,
                "disease_context_unique_targets": set(),
                "bulk_up_target_hits": set(),
                "bulk_down_target_hits": set(),
                "aligned_target_genes": set(),
                "b_cell_expressed_targets": set(),
                "b_cell_disease_increased_targets": set(),
                "mdc_expressed_targets": set(),
                "mdc_disease_increased_targets": set(),
                "top_b_cell_target_evidence": [],
            },
        )
        gene = edge["gene_symbol"]
        entry["disease_context_edge_rows"] += 1
        entry["disease_context_unique_targets"].add(gene)
        if edge["in_bulk_uc_ibd_up150"] == "yes":
            entry["bulk_up_target_hits"].add(gene)
        if edge["in_bulk_uc_ibd_down150"] == "yes":
            entry["bulk_down_target_hits"].add(gene)
        if edge["therapeutic_alignment"].startswith("aligned"):
            entry["aligned_target_genes"].add(gene)
        for celltype, expressed_key, increased_key in [
            ("B", "b_cell_expressed_targets", "b_cell_disease_increased_targets"),
            ("M/DC", "mdc_expressed_targets", "mdc_disease_increased_targets"),
        ]:
            row = contrast.get((celltype, gene))
            if not row:
                continue
            if float(row["mean_detection_diseased"]) > 0 or float(row["mean_log1p_cpm_diseased"]) > 0:
                entry[expressed_key].add(gene)
            if float(row["delta_log1p_cpm_diseased_minus_healthy"]) > 0:
                entry[increased_key].add(gene)
        if edge["ingredient_name"].lower() == "curcumin" and gene in {"IL10", "IL1B", "PTGS2", "MMP9", "CXCL8", "CCL2", "SYK", "BLNK", "BCL6"}:
            entry["top_b_cell_target_evidence"].append(
                f"{gene}|PMID:{edge.get('pubmed_id','')}|{edge.get('relationship','')}|{edge.get('reference_title','')}"
            )

    summary_rows: List[Dict[str, Any]] = []
    for entry in by_candidate.values():
        row = {
            key: value
            for key, value in entry.items()
            if not isinstance(value, set) and key != "top_b_cell_target_evidence"
        }
        for key in [
            "disease_context_unique_targets",
            "bulk_up_target_hits",
            "bulk_down_target_hits",
            "aligned_target_genes",
            "b_cell_expressed_targets",
            "b_cell_disease_increased_targets",
            "mdc_expressed_targets",
            "mdc_disease_increased_targets",
        ]:
            genes = sorted(entry[key])
            row[f"n_{key}"] = len(genes)
            row[key] = ";".join(genes)
        row["top_b_cell_target_evidence"] = " || ".join(entry["top_b_cell_target_evidence"][:8])
        row["bridge_score"] = (
            len(entry["bulk_up_target_hits"]) * 2
            + len(entry["aligned_target_genes"]) * 3
            + len(entry["b_cell_expressed_targets"]) * 1.5
            + len(entry["b_cell_disease_increased_targets"]) * 2
            + len(entry["mdc_expressed_targets"]) * 0.75
        )
        row["bridge_score"] = f"{row['bridge_score']:.3g}"
        summary_rows.append(row)
    summary_rows.sort(key=lambda row: (-float(row["bridge_score"]), row["ingredient_name"].lower()))
    return summary_rows


def safe_float(value: Any) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return 0.0


def joined_unique(values: Iterable[str]) -> str:
    clean = sorted({value for value in values if value})
    return ";".join(clean)


def contrast_lookup(
    contrast_rows: Sequence[Dict[str, Any]],
    celltype: str,
    gene: str,
) -> Dict[str, Any]:
    for row in contrast_rows:
        if row["tissue_assignment"] == "R" and row["celltype"] == celltype and row["gene_symbol"] == gene:
            return row
    return {}


def build_curcumin_target_evidence(
    edges: Sequence[Dict[str, Any]],
    contrast_rows: Sequence[Dict[str, Any]],
) -> List[Dict[str, Any]]:
    by_gene: Dict[str, List[Dict[str, Any]]] = defaultdict(list)
    for edge in edges:
        if edge["ingredient_name"].lower() == "curcumin":
            by_gene[edge["gene_symbol"]].append(edge)

    b_mechanism_genes = {"BCL6", "BLNK", "SYK", "SH3KBP1"}
    rows: List[Dict[str, Any]] = []
    for gene, gene_edges in by_gene.items():
        b_row = contrast_lookup(contrast_rows, "B", gene)
        mdc_row = contrast_lookup(contrast_rows, "M/DC", gene)
        alignments = {edge["therapeutic_alignment"] for edge in gene_edges}
        directions = {edge["relationship_direction"] for edge in gene_edges}
        in_up = any(edge["in_bulk_uc_ibd_up150"] == "yes" for edge in gene_edges)
        in_down = any(edge["in_bulk_uc_ibd_down150"] == "yes" for edge in gene_edges)
        aligned = any(str(alignment).startswith("aligned") for alignment in alignments)
        b_delta = safe_float(b_row.get("delta_log1p_cpm_diseased_minus_healthy"))
        mdc_delta = safe_float(mdc_row.get("delta_log1p_cpm_diseased_minus_healthy"))
        b_detected = safe_float(b_row.get("mean_detection_diseased"))
        mdc_detected = safe_float(mdc_row.get("mean_detection_diseased"))
        b_status = "no_rectal_B_contrast" if not b_row else ("detected_in_diseased_B" if b_detected > 0 else "not_detected_in_diseased_B")
        mdc_status = "no_rectal_MDC_contrast" if not mdc_row else ("detected_in_diseased_MDC" if mdc_detected > 0 else "not_detected_in_diseased_MDC")

        if aligned and b_delta > 0:
            tier = "A_module_aligned_b_cell_increased"
        elif aligned:
            tier = "B_module_aligned"
        elif gene in b_mechanism_genes and b_delta > 0:
            tier = "C_b_cell_mechanism_increased"
        elif b_delta > 0 and b_detected > 0:
            tier = "D_b_cell_increased"
        elif mdc_delta > 0 and mdc_detected > 0:
            tier = "E_mdc_increased"
        else:
            tier = "F_context_target_only"

        if in_up and in_down:
            bulk_role = "bulk_up_and_down_top150"
        elif in_up:
            bulk_role = "bulk_up_top150"
        elif in_down:
            bulk_role = "bulk_down_top150"
        else:
            bulk_role = "not_bulk_top150"

        rows.append(
            {
                "gene_symbol": gene,
                "evidence_tier": tier,
                "bulk_module_role": bulk_role,
                "therapeutic_alignments": joined_unique(alignments),
                "relationship_directions": joined_unique(directions),
                "relationships": joined_unique(edge.get("relationship", "") for edge in gene_edges),
                "pubmed_ids": joined_unique(edge.get("pubmed_id", "") for edge in gene_edges),
                "reference_titles": " || ".join(sorted({edge.get("reference_title", "") for edge in gene_edges if edge.get("reference_title", "")})),
                "b_mean_log1p_cpm_diseased": b_row.get("mean_log1p_cpm_diseased", ""),
                "b_mean_log1p_cpm_healthy": b_row.get("mean_log1p_cpm_healthy", ""),
                "b_delta_log1p_cpm_diseased_minus_healthy": b_row.get("delta_log1p_cpm_diseased_minus_healthy", ""),
                "b_p_log1p_cpm_mannwhitney": b_row.get("p_log1p_cpm_mannwhitney", ""),
                "b_mean_detection_diseased": b_row.get("mean_detection_diseased", ""),
                "b_mean_detection_healthy": b_row.get("mean_detection_healthy", ""),
                "b_expression_status": b_status,
                "mdc_mean_log1p_cpm_diseased": mdc_row.get("mean_log1p_cpm_diseased", ""),
                "mdc_mean_log1p_cpm_healthy": mdc_row.get("mean_log1p_cpm_healthy", ""),
                "mdc_delta_log1p_cpm_diseased_minus_healthy": mdc_row.get("delta_log1p_cpm_diseased_minus_healthy", ""),
                "mdc_p_log1p_cpm_mannwhitney": mdc_row.get("p_log1p_cpm_mannwhitney", ""),
                "mdc_mean_detection_diseased": mdc_row.get("mean_detection_diseased", ""),
                "mdc_mean_detection_healthy": mdc_row.get("mean_detection_healthy", ""),
                "mdc_expression_status": mdc_status,
                "first_supporting_sentence": next((edge.get("supporting_sentence", "") for edge in gene_edges if edge.get("supporting_sentence", "")), ""),
            }
        )

    rows.sort(
        key=lambda row: (
            row["evidence_tier"],
            -safe_float(row["b_delta_log1p_cpm_diseased_minus_healthy"]),
            -safe_float(row["mdc_delta_log1p_cpm_diseased_minus_healthy"]),
            row["gene_symbol"],
        )
    )
    return rows


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--edge-table", default="results/m8_herb_ingredient_targets/herb_ingredient_target_edges.tsv")
    parser.add_argument("--up-genes", default="results/m4_geo_uc_ibd_signature/lincs_query_up_genes.txt")
    parser.add_argument("--down-genes", default="results/m4_geo_uc_ibd_signature/lincs_query_down_genes.txt")
    parser.add_argument("--raw-dir", default="data/raw/geo/GSE125527")
    parser.add_argument("--sample-dir", default="data/raw/geo/GSE125527/raw_rectal")
    parser.add_argument("--out-dir", default=str(OUT_DIR))
    args = parser.parse_args()

    out_dir = Path(args.out_dir)
    up_genes = read_gene_list(Path(args.up_genes))
    down_genes = read_gene_list(Path(args.down_genes))
    disease_edges = load_disease_context_edges(Path(args.edge_table), up_genes, down_genes)
    target_genes = {row["gene_symbol"] for row in disease_edges}
    expression_rows, present_genes = score_target_expression(Path(args.sample_dir), Path(args.raw_dir), target_genes)
    contrast_rows = contrast_gene_expression(expression_rows)
    summary_rows = summarize_candidate_bridge(disease_edges, contrast_rows)
    curcumin_target_rows = build_curcumin_target_evidence(disease_edges, contrast_rows)

    edge_fields = [
        "ingredient_id",
        "ingredient_name",
        "target_id",
        "gene_symbol",
        "protein_name",
        "reference_id",
        "pubmed_id",
        "pubmed_url",
        "reference_title",
        "relationship",
        "relationship_direction",
        "grade",
        "supporting_sentence",
        "disease_context_terms",
        "in_bulk_uc_ibd_up150",
        "in_bulk_uc_ibd_down150",
        "therapeutic_alignment",
        "candidate_priority_class",
        "candidate_priority_score",
    ]
    expression_fields = [
        "patient_assignment",
        "tissue_assignment",
        "disease_assignment",
        "celltype",
        "gene_symbol",
        "n_cells",
        "sum_counts",
        "detected_cells",
        "detection_fraction",
        "total_umi",
        "log1p_cpm",
    ]
    contrast_fields = [
        "tissue_assignment",
        "celltype",
        "gene_symbol",
        "n_diseased_samples",
        "n_healthy_samples",
        "diseased_cells",
        "healthy_cells",
        "mean_log1p_cpm_diseased",
        "mean_log1p_cpm_healthy",
        "delta_log1p_cpm_diseased_minus_healthy",
        "p_log1p_cpm_mannwhitney",
        "mean_detection_diseased",
        "mean_detection_healthy",
        "delta_detection_diseased_minus_healthy",
        "p_detection_mannwhitney",
    ]
    summary_fields = [
        "ingredient_id",
        "ingredient_name",
        "candidate_priority_class",
        "candidate_priority_score",
        "bridge_score",
        "disease_context_edge_rows",
        "n_disease_context_unique_targets",
        "disease_context_unique_targets",
        "n_bulk_up_target_hits",
        "bulk_up_target_hits",
        "n_bulk_down_target_hits",
        "bulk_down_target_hits",
        "n_aligned_target_genes",
        "aligned_target_genes",
        "n_b_cell_expressed_targets",
        "b_cell_expressed_targets",
        "n_b_cell_disease_increased_targets",
        "b_cell_disease_increased_targets",
        "n_mdc_expressed_targets",
        "mdc_expressed_targets",
        "n_mdc_disease_increased_targets",
        "mdc_disease_increased_targets",
        "top_b_cell_target_evidence",
    ]
    curcumin_target_fields = [
        "gene_symbol",
        "evidence_tier",
        "bulk_module_role",
        "therapeutic_alignments",
        "relationship_directions",
        "relationships",
        "pubmed_ids",
        "reference_titles",
        "b_mean_log1p_cpm_diseased",
        "b_mean_log1p_cpm_healthy",
        "b_delta_log1p_cpm_diseased_minus_healthy",
        "b_p_log1p_cpm_mannwhitney",
        "b_mean_detection_diseased",
        "b_mean_detection_healthy",
        "b_expression_status",
        "mdc_mean_log1p_cpm_diseased",
        "mdc_mean_log1p_cpm_healthy",
        "mdc_delta_log1p_cpm_diseased_minus_healthy",
        "mdc_p_log1p_cpm_mannwhitney",
        "mdc_mean_detection_diseased",
        "mdc_mean_detection_healthy",
        "mdc_expression_status",
        "first_supporting_sentence",
    ]

    write_tsv(out_dir / "disease_context_herb_target_edges.tsv", disease_edges, edge_fields)
    write_tsv(out_dir / "gse125527_target_gene_pseudobulk.tsv", expression_rows, expression_fields)
    write_tsv(out_dir / "gse125527_target_gene_contrast.tsv", contrast_rows, contrast_fields)
    write_tsv(out_dir / "candidate_target_cell_bridge_summary.tsv", summary_rows, summary_fields)
    write_tsv(out_dir / "curcumin_bcell_target_evidence.tsv", curcumin_target_rows, curcumin_target_fields)

    top_rows = summary_rows[:12]
    curcumin = next((row for row in summary_rows if row["ingredient_name"].lower() == "curcumin"), None)
    report_lines = [
        "# M9 Disease-Context Target-to-Cell Bridge Run",
        "",
        "## Scope",
        "",
        "This run filters HERB `drug_paper_target` edges to UC/IBD/colitis-relevant references, classifies relationship direction, overlays targets with the bulk UC/IBD module, and quantifies target-gene expression in GSE125527 rectal scRNA cell types.",
        "",
        "## Key Outputs",
        "",
        "- `results/m9_target_cell_bridge/disease_context_herb_target_edges.tsv`",
        "- `results/m9_target_cell_bridge/gse125527_target_gene_pseudobulk.tsv`",
        "- `results/m9_target_cell_bridge/gse125527_target_gene_contrast.tsv`",
        "- `results/m9_target_cell_bridge/candidate_target_cell_bridge_summary.tsv`",
        "- `results/m9_target_cell_bridge/curcumin_bcell_target_evidence.tsv`",
        "",
        "## Acquisition And Filtering Summary",
        "",
        f"- Disease-context HERB target edge rows: {len(disease_edges)}",
        f"- Disease-context unique target genes: {len(target_genes)}",
        f"- Target genes present in GSE125527 rectal UMI tables: {len(present_genes)}",
        "",
        "## Top Candidate Target-to-Cell Bridges",
        "",
        "| Ingredient | Score | Disease-context targets | Bulk up hits | B-cell expressed | B-cell disease-increased | M/DC expressed |",
        "|---|---:|---:|---:|---:|---:|---:|",
    ]
    for row in top_rows:
        report_lines.append(
            f"| {row['ingredient_name']} | {row['bridge_score']} | {row['n_disease_context_unique_targets']} | {row['n_bulk_up_target_hits']} | {row['n_b_cell_expressed_targets']} | {row['n_b_cell_disease_increased_targets']} | {row['n_mdc_expressed_targets']} |"
        )
    report_lines.extend(["", "## Curcumin Focus", ""])
    if curcumin:
        report_lines.extend(
            [
                f"- Curcumin disease-context targets: {curcumin['n_disease_context_unique_targets']}.",
                f"- Curcumin targets overlapping the bulk UC/IBD up module: {curcumin['n_bulk_up_target_hits']} ({curcumin['bulk_up_target_hits']}).",
                f"- Curcumin targets expressed in rectal B cells: {curcumin['n_b_cell_expressed_targets']} ({curcumin['b_cell_expressed_targets']}).",
                f"- Curcumin targets increased in diseased rectal B-cell pseudobulk: {curcumin['n_b_cell_disease_increased_targets']} ({curcumin['b_cell_disease_increased_targets']}).",
            ]
        )
    report_lines.extend(
        [
            "",
            "### Figure-Ready Curcumin Target Evidence",
            "",
            "| Gene | Tier | Bulk role | B-cell delta | B-cell P | M/DC delta | PMIDs |",
            "|---|---|---|---:|---:|---:|---|",
        ]
    )
    for row in curcumin_target_rows[:12]:
        report_lines.append(
            f"| {row['gene_symbol']} | {row['evidence_tier']} | {row['bulk_module_role']} | {row['b_delta_log1p_cpm_diseased_minus_healthy']} | {row['b_p_log1p_cpm_mannwhitney']} | {row['mdc_delta_log1p_cpm_diseased_minus_healthy']} | {row['pubmed_ids']} |"
        )
    report_lines.extend(
        [
            "",
            "## Interpretation Boundary",
            "",
            "- The disease-context filter is text-based and should be manually audited for the final figure/table.",
            "- Direction labels are inherited from HERB relationship text and supporting sentences; they are evidence annotations, not new experimental validation.",
            "- A target is considered expressed if it is detected in patient-level pseudobulk for the cell class; disease increase is based on patient-level diseased-vs-healthy pseudobulk contrast.",
            "- Celltype aggregation now deduplicates broad and detailed labels when no cluster label is available, preventing empty-cluster cells from being counted twice in the same celltype key.",
            "- This module supports mechanistic prioritization. It does not prove therapeutic efficacy.",
            "",
        ]
    )
    DOC_PATH.write_text("\n".join(report_lines), encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
