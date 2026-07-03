#!/usr/bin/env python3
"""Run no-new-download robustness analyses for the UC B-cell/curcumin study."""

from __future__ import annotations

import csv
import gzip
import math
from collections import defaultdict
from pathlib import Path
from statistics import mean
from typing import Any, Dict, Iterable, List, Sequence, Tuple

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, spearmanr


ROOT = Path(__file__).resolve().parents[1]
OUT_DIR = ROOT / "results" / "m17_no_download_robustness"
RAW_DIR = ROOT / "data" / "raw" / "geo" / "GSE125527"
SAMPLE_DIR = RAW_DIR / "raw_rectal"


B_MARKERS = [
    "MS4A1",
    "CD19",
    "CD79A",
    "CD79B",
    "CD74",
    "BANK1",
    "BLK",
    "PAX5",
    "BCL6",
    "BLNK",
    "SYK",
    "TNFRSF13C",
]

PLASMA_MARKERS = [
    "MZB1",
    "XBP1",
    "JCHAIN",
    "SDC1",
    "PRDM1",
    "TNFRSF17",
    "IGHG1",
    "IGHG3",
    "IGHA1",
    "IGKC",
    "DERL3",
]


def read_gzip_lines(path: Path) -> Iterable[str]:
    with gzip.open(path, "rt", encoding="utf-8", errors="replace", newline="") as handle:
        yield from handle


def write_tsv(path: Path, rows: Sequence[Dict[str, Any]], fields: Sequence[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(fields), delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: "" if row.get(field) is None else row.get(field, "") for field in fields})


def split_genes(value: Any) -> set[str]:
    if pd.isna(value) or str(value).strip() == "":
        return set()
    return {part.strip().upper() for part in str(value).split(";") if part.strip()}


def fmt(value: float, digits: int = 6) -> str:
    if value is None or not math.isfinite(float(value)):
        return ""
    return f"{float(value):.{digits}g}"


def bh_fdr(p_values: Sequence[float]) -> list[float]:
    p = np.array(p_values, dtype=float)
    adjusted = np.full(len(p), np.nan)
    valid = np.isfinite(p)
    indices = np.where(valid)[0]
    if not len(indices):
        return adjusted.tolist()
    ordered = indices[np.argsort(p[indices])]
    running = 1.0
    m = len(ordered)
    for rank, idx in reversed(list(enumerate(ordered, start=1))):
        running = min(running, p[idx] * m / rank)
        adjusted[idx] = min(running, 1.0)
    return adjusted.tolist()


def safe_mannwhitney(a: Sequence[float], b: Sequence[float]) -> float:
    if len(a) < 2 or len(b) < 2:
        return math.nan
    try:
        return float(mannwhitneyu(a, b, alternative="two-sided").pvalue)
    except ValueError:
        return math.nan


def load_metadata_map() -> dict[tuple[str, str, str], dict[str, str]]:
    mapping: dict[tuple[str, str, str], dict[str, str]] = {}
    reader = csv.reader(read_gzip_lines(RAW_DIR / "GSE125527_cell_metadata.csv.gz"))
    next(reader)
    for raw in reader:
        if len(raw) < 6:
            continue
        row = {
            "row_id": raw[0],
            "cell_barcode": raw[1],
            "patient_assignment": raw[2],
            "tissue_assignment": raw[3],
            "disease_assignment": raw[4],
            "celltype": raw[5],
        }
        mapping[(row["patient_assignment"], row["tissue_assignment"], row["cell_barcode"])] = row
    return mapping


def load_patient_id_map() -> dict[str, str]:
    path = RAW_DIR / "GSE125527_oldPatientId-newPatientId.csv.gz"
    if not path.exists():
        return {}
    reader = csv.DictReader(read_gzip_lines(path))
    return {row["old_id"]: row["new_id"] for row in reader}


def parse_sample_from_filename(path: Path) -> tuple[str, str]:
    parts = path.name.split("_")
    if len(parts) < 4:
        return "", ""
    return parts[1], parts[2]


def signature_from_consensus(n: int, mode: str) -> dict[str, Any]:
    consensus = pd.read_csv(ROOT / "results/m4_geo_uc_ibd_signature/consensus_gene_signature.tsv", sep="\t")
    up = consensus[consensus["direction"].eq("up")].sort_values("consensus_score", ascending=False).head(n)
    down = consensus[consensus["direction"].eq("down")].sort_values("consensus_score", ascending=False).head(n)
    return {
        "variant": f"consensus_top{n}_{mode}",
        "source": "two_evidence_group_consensus",
        "mode": mode,
        "up_weights": dict(zip(up["gene_symbol"].str.upper(), up["mean_log2fc"].abs().astype(float))),
        "down_weights": dict(zip(down["gene_symbol"].str.upper(), down["mean_log2fc"].abs().astype(float))),
    }


def signature_from_dataset(accession: str, n: int = 150) -> dict[str, Any]:
    sig = pd.read_csv(ROOT / f"results/m4_geo_uc_ibd_signature/{accession}_gene_signature.tsv", sep="\t")
    up = sig[sig["log2fc"] > 0].assign(abs_fc=lambda x: x["log2fc"].abs()).sort_values(["fdr", "abs_fc"], ascending=[True, False]).head(n)
    down = sig[sig["log2fc"] < 0].assign(abs_fc=lambda x: x["log2fc"].abs()).sort_values(["fdr", "abs_fc"], ascending=[True, False]).head(n)
    return {
        "variant": f"{accession}_top{n}_avg_logcpm",
        "source": accession,
        "mode": "avg_logcpm",
        "up_weights": dict(zip(up["gene_symbol"].str.upper(), up["log2fc"].abs().astype(float))),
        "down_weights": dict(zip(down["gene_symbol"].str.upper(), down["log2fc"].abs().astype(float))),
    }


def load_signature_variants() -> list[dict[str, Any]]:
    variants: list[dict[str, Any]] = []
    for n in [50, 100, 150, 300]:
        variants.append(signature_from_consensus(n, "logsum_cpm"))
        variants.append(signature_from_consensus(n, "avg_logcpm"))
    variants.append(signature_from_consensus(150, "weighted_avg_logcpm"))
    variants.append(signature_from_consensus(150, "weighted_logsum_cpm"))
    for accession in ["GSE75214", "GSE87466"]:
        dataset_variant = signature_from_dataset(accession, 150)
        variants.append({**dataset_variant, "variant": f"{accession}_top150_logsum_cpm", "mode": "logsum_cpm"})
        variants.append(signature_from_dataset(accession, 150))
    return variants


def score_gene_module(counts: list[float], idx_weights: list[tuple[int, float]], scale: float, mode: str) -> float:
    if not idx_weights:
        return math.nan
    if mode == "logsum_cpm":
        return math.log1p(sum(counts[idx] for idx, _ in idx_weights) * scale)
    if mode == "weighted_logsum_cpm":
        denom = sum(weight for _, weight in idx_weights)
        if denom <= 0:
            return math.nan
        weighted_counts = sum(counts[idx] * weight for idx, weight in idx_weights) / denom * len(idx_weights)
        return math.log1p(weighted_counts * scale)
    if mode == "weighted_avg_logcpm":
        denom = sum(weight for _, weight in idx_weights)
        if denom <= 0:
            return math.nan
        return sum(math.log1p(counts[idx] * scale) * weight for idx, weight in idx_weights) / denom
    return sum(math.log1p(counts[idx] * scale) for idx, _ in idx_weights) / len(idx_weights)


def summarize_variant_pseudobulk() -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    variants = load_signature_variants()
    metadata = load_metadata_map()
    patient_map = load_patient_id_map()
    grouped: dict[tuple[str, str, str, str, str], dict[str, float]] = defaultdict(lambda: {"score_sum": 0.0, "n_cells": 0.0})
    coverage_rows: list[dict[str, Any]] = []

    first_sample = next(iter(sorted(SAMPLE_DIR.glob("*_cell-gene_UMI_table.tsv.gz"))))
    with gzip.open(first_sample, "rt", encoding="utf-8", errors="replace", newline="") as handle:
        header = next(csv.reader(handle, delimiter="\t"))
    gene_to_idx = {gene.upper(): idx for idx, gene in enumerate(header[1:])}

    prepared: list[dict[str, Any]] = []
    for variant in variants:
        up_weights = variant["up_weights"]
        down_weights = variant["down_weights"]
        up_idx = [(gene_to_idx[gene], weight) for gene, weight in up_weights.items() if gene in gene_to_idx]
        down_idx = [(gene_to_idx[gene], weight) for gene, weight in down_weights.items() if gene in gene_to_idx]
        prepared.append({**variant, "up_idx": up_idx, "down_idx": down_idx})
        coverage_rows.append(
            {
                "variant": variant["variant"],
                "source": variant["source"],
                "mode": variant["mode"],
                "up_query_genes": len(up_weights),
                "down_query_genes": len(down_weights),
                "up_present_genes": len(up_idx),
                "down_present_genes": len(down_idx),
                "up_coverage_fraction": fmt(len(up_idx) / len(up_weights) if up_weights else math.nan),
                "down_coverage_fraction": fmt(len(down_idx) / len(down_weights) if down_weights else math.nan),
            }
        )

    for path in sorted(SAMPLE_DIR.glob("*_cell-gene_UMI_table.tsv.gz")):
        old_patient, tissue = parse_sample_from_filename(path)
        patient = patient_map.get(old_patient, old_patient)
        with gzip.open(path, "rt", encoding="utf-8", errors="replace", newline="") as handle:
            reader = csv.reader(handle, delimiter="\t")
            next(reader)
            for raw in reader:
                if not raw:
                    continue
                barcode = raw[0]
                counts = [float(value) if value else 0.0 for value in raw[1:]]
                total = sum(counts)
                if total <= 0:
                    continue
                meta = metadata.get((patient, tissue, barcode), {})
                disease = meta.get("disease_assignment") or ("diseased" if patient.startswith("U") else "healthy")
                celltype = meta.get("celltype", "unknown")
                scale = 10000.0 / total
                for variant in prepared:
                    up_score = score_gene_module(counts, variant["up_idx"], scale, variant["mode"])
                    down_score = score_gene_module(counts, variant["down_idx"], scale, variant["mode"])
                    if not math.isfinite(up_score) or not math.isfinite(down_score):
                        continue
                    key = (variant["variant"], patient, tissue, disease, celltype)
                    grouped[key]["score_sum"] += up_score - down_score
                    grouped[key]["n_cells"] += 1.0

    pseudobulk_rows: list[dict[str, Any]] = []
    for (variant, patient, tissue, disease, celltype), values in sorted(grouped.items()):
        if values["n_cells"] <= 0:
            continue
        pseudobulk_rows.append(
            {
                "variant": variant,
                "patient_assignment": patient,
                "tissue_assignment": tissue,
                "disease_assignment": disease,
                "celltype": celltype,
                "n_cells": int(values["n_cells"]),
                "mean_disease_axis_score": values["score_sum"] / values["n_cells"],
            }
        )
    return pseudobulk_rows, coverage_rows


def contrast_and_null(
    pseudobulk_rows: Sequence[dict[str, Any]],
    n_permutations: int = 10000,
    seed: int = 20260702,
) -> list[dict[str, Any]]:
    rng = np.random.default_rng(seed)
    rows: list[dict[str, Any]] = []
    by_variant_cell: dict[tuple[str, str, str], list[dict[str, Any]]] = defaultdict(list)
    for row in pseudobulk_rows:
        if row["tissue_assignment"] != "R" or row["celltype"] == "unknown":
            continue
        by_variant_cell[(row["variant"], row["tissue_assignment"], row["celltype"])].append(row)

    for (variant, tissue, celltype), group in sorted(by_variant_cell.items()):
        diseased = [float(row["mean_disease_axis_score"]) for row in group if row["disease_assignment"] == "diseased"]
        healthy = [float(row["mean_disease_axis_score"]) for row in group if row["disease_assignment"] == "healthy"]
        if not diseased or not healthy:
            continue
        values = np.array([float(row["mean_disease_axis_score"]) for row in group])
        labels = np.array([row["disease_assignment"] for row in group])
        observed = float(np.mean(values[labels == "diseased"]) - np.mean(values[labels == "healthy"]))
        diseased_count = int(np.sum(labels == "diseased"))
        null_values = np.empty(n_permutations, dtype=float)
        all_idx = np.arange(len(values))
        for i in range(n_permutations):
            perm_d = rng.choice(all_idx, size=diseased_count, replace=False)
            mask = np.zeros(len(values), dtype=bool)
            mask[perm_d] = True
            null_values[i] = float(np.mean(values[mask]) - np.mean(values[~mask]))
        empirical_p = (int(np.sum(np.abs(null_values) >= abs(observed))) + 1.0) / (n_permutations + 1.0)
        rows.append(
            {
                "variant": variant,
                "tissue_assignment": tissue,
                "celltype": celltype,
                "n_diseased_samples": len(diseased),
                "n_healthy_samples": len(healthy),
                "delta_axis_diseased_minus_healthy": observed,
                "p_mannwhitney": safe_mannwhitney(diseased, healthy),
                "empirical_p_two_sided": empirical_p,
                "null_q025_delta": float(np.quantile(null_values, 0.025)),
                "null_q975_delta": float(np.quantile(null_values, 0.975)),
                "n_permutations": n_permutations,
                "seed": seed,
            }
        )

    out: list[dict[str, Any]] = []
    for variant in sorted({row["variant"] for row in rows}):
        idxs = [idx for idx, row in enumerate(rows) if row["variant"] == variant]
        q_mw = bh_fdr([rows[idx]["p_mannwhitney"] for idx in idxs])
        q_emp = bh_fdr([rows[idx]["empirical_p_two_sided"] for idx in idxs])
        for local_idx, idx in enumerate(idxs):
            row = rows[idx].copy()
            row["q_mannwhitney_bh"] = q_mw[local_idx]
            row["q_empirical_bh"] = q_emp[local_idx]
            out.append(row)
    out.sort(key=lambda row: (row["variant"], -float(row["delta_axis_diseased_minus_healthy"])))
    return out


def run_signature_sensitivity() -> None:
    pseudobulk_rows, coverage_rows = summarize_variant_pseudobulk()
    contrast_rows = contrast_and_null(pseudobulk_rows)
    fields = [
        "variant",
        "tissue_assignment",
        "celltype",
        "n_diseased_samples",
        "n_healthy_samples",
        "delta_axis_diseased_minus_healthy",
        "p_mannwhitney",
        "q_mannwhitney_bh",
        "empirical_p_two_sided",
        "q_empirical_bh",
        "null_q025_delta",
        "null_q975_delta",
        "n_permutations",
        "seed",
    ]
    formatted = [{k: fmt(v) if isinstance(v, float) else v for k, v in row.items()} for row in contrast_rows]
    write_tsv(OUT_DIR / "scrna_signature_sensitivity_by_celltype.tsv", formatted, fields)
    write_tsv(
        OUT_DIR / "scrna_signature_variant_gene_coverage.tsv",
        coverage_rows,
        [
            "variant",
            "source",
            "mode",
            "up_query_genes",
            "down_query_genes",
            "up_present_genes",
            "down_present_genes",
            "up_coverage_fraction",
            "down_coverage_fraction",
        ],
    )
    b_rows = [row for row in formatted if row["celltype"] == "B"]
    write_tsv(OUT_DIR / "scrna_signature_sensitivity_bcell_summary.tsv", b_rows, fields)


def run_herb_decoy_null() -> None:
    m9 = pd.read_csv(ROOT / "results/m9_target_cell_bridge/candidate_target_cell_bridge_summary.tsv", sep="\t")
    cur = m9[m9["ingredient_name"].str.lower().eq("curcumin")].iloc[0]
    disease_edges = pd.read_csv(ROOT / "results/m9_target_cell_bridge/disease_context_herb_target_edges.tsv", sep="\t")
    contrast = pd.read_csv(ROOT / "results/m9_target_cell_bridge/gse125527_target_gene_contrast.tsv", sep="\t")
    rectal = contrast[(contrast["tissue_assignment"].eq("R")) & (contrast["celltype"].isin(["B", "M/DC"]))].copy()
    bulk_up = set(disease_edges.loc[disease_edges["in_bulk_uc_ibd_up150"].eq("yes"), "gene_symbol"].dropna().str.upper())
    b_expr = set(rectal[(rectal["celltype"].eq("B")) & ((rectal["mean_detection_diseased"] > 0) | (rectal["mean_log1p_cpm_diseased"] > 0))]["gene_symbol"].str.upper())
    b_inc = set(rectal[(rectal["celltype"].eq("B")) & (rectal["delta_log1p_cpm_diseased_minus_healthy"] > 0)]["gene_symbol"].str.upper())
    mdc_expr = set(rectal[(rectal["celltype"].eq("M/DC")) & ((rectal["mean_detection_diseased"] > 0) | (rectal["mean_log1p_cpm_diseased"] > 0))]["gene_symbol"].str.upper())
    universe = sorted(set(disease_edges["gene_symbol"].dropna().str.upper()) & (set(rectal["gene_symbol"].dropna().str.upper()) | bulk_up))

    target_count = int(cur["n_disease_context_unique_targets"])
    cur_targets = split_genes(cur["disease_context_unique_targets"])

    def score(targets: set[str]) -> tuple[float, int, int, int, int]:
        n_bulk = len(targets & bulk_up)
        n_b_expr = len(targets & b_expr)
        n_b_inc = len(targets & b_inc)
        n_mdc_expr = len(targets & mdc_expr)
        return n_bulk * 2.0 + n_b_expr * 1.5 + n_b_inc * 2.0 + n_mdc_expr * 0.75, n_bulk, n_b_expr, n_b_inc, n_mdc_expr

    observed_score, observed_bulk, observed_b_expr, observed_b_inc, observed_mdc_expr = score(cur_targets)
    rng = np.random.default_rng(20260702)
    n_perm = 20000
    scores = np.empty(n_perm, dtype=float)
    for i in range(n_perm):
        sampled = set(rng.choice(universe, size=target_count, replace=False))
        scores[i] = score(sampled)[0]
    empirical_p = (int(np.sum(scores >= observed_score)) + 1.0) / (n_perm + 1.0)

    rank_rows: list[dict[str, Any]] = []
    for _, row in m9.iterrows():
        targets = split_genes(row["disease_context_unique_targets"])
        target_set_score, n_bulk, n_b_expr, n_b_inc, n_mdc_expr = score(targets)
        rank_rows.append(
            {
                "ingredient_id": row["ingredient_id"],
                "ingredient_name": "Berberine" if row["ingredient_name"] == "Berberime" else row["ingredient_name"],
                "disease_context_targets": int(row["n_disease_context_unique_targets"]),
                "manuscript_bridge_score": row["bridge_score"],
                "target_set_specificity_score_without_direction_alignment": fmt(target_set_score),
                "bulk_up_hits": n_bulk,
                "b_cell_expressed": n_b_expr,
                "b_cell_increased": n_b_inc,
                "mdc_expressed": n_mdc_expr,
            }
        )
    rank_rows.sort(key=lambda r: (-float(r["target_set_specificity_score_without_direction_alignment"]), r["ingredient_name"]))
    for idx, row in enumerate(rank_rows, start=1):
        row["target_set_specificity_rank"] = idx

    rho, rho_p = spearmanr(
        [row["disease_context_targets"] for row in rank_rows],
        [float(row["target_set_specificity_score_without_direction_alignment"]) for row in rank_rows],
    )
    summary_rows = [
        {
            "analysis": "matched_target_count_random_target_set_null",
            "observed_candidate": "Curcumin",
            "target_count": target_count,
            "observed_target_set_score_without_direction_alignment": fmt(observed_score),
            "observed_bulk_up_hits": observed_bulk,
            "observed_b_cell_expressed": observed_b_expr,
            "observed_b_cell_increased": observed_b_inc,
            "observed_mdc_expressed": observed_mdc_expr,
            "null_iterations": n_perm,
            "null_mean": fmt(float(np.mean(scores))),
            "null_sd": fmt(float(np.std(scores, ddof=1))),
            "null_q95": fmt(float(np.quantile(scores, 0.95))),
            "null_q99": fmt(float(np.quantile(scores, 0.99))),
            "empirical_p_ge_observed": fmt(empirical_p),
            "target_count_score_spearman_rho_among_candidates": fmt(float(rho)),
            "target_count_score_spearman_p": fmt(float(rho_p)),
        }
    ]
    write_tsv(
        OUT_DIR / "herb_decoy_null_summary.tsv",
        summary_rows,
        [
            "analysis",
            "observed_candidate",
            "target_count",
            "observed_target_set_score_without_direction_alignment",
            "observed_bulk_up_hits",
            "observed_b_cell_expressed",
            "observed_b_cell_increased",
            "observed_mdc_expressed",
            "null_iterations",
            "null_mean",
            "null_sd",
            "null_q95",
            "null_q99",
            "empirical_p_ge_observed",
            "target_count_score_spearman_rho_among_candidates",
            "target_count_score_spearman_p",
        ],
    )
    write_tsv(
        OUT_DIR / "herb_candidate_specificity_ranking.tsv",
        rank_rows,
        [
            "target_set_specificity_rank",
            "ingredient_id",
            "ingredient_name",
            "disease_context_targets",
            "manuscript_bridge_score",
            "target_set_specificity_score_without_direction_alignment",
            "bulk_up_hits",
            "b_cell_expressed",
            "b_cell_increased",
            "mdc_expressed",
        ],
    )


def run_bulk_marker_overlap() -> None:
    rows: list[dict[str, Any]] = []
    marker_group = {gene: "B-cell marker" for gene in B_MARKERS}
    marker_group.update({gene: "Plasma-cell marker" for gene in PLASMA_MARKERS})
    datasets = ["GSE75214", "GSE59071", "GSE87466"]
    for accession in datasets:
        sig = pd.read_csv(ROOT / f"results/m4_geo_uc_ibd_signature/{accession}_gene_signature.tsv", sep="\t")
        sig["gene_symbol"] = sig["gene_symbol"].str.upper()
        for gene, group in marker_group.items():
            hit = sig[sig["gene_symbol"].eq(gene)]
            if hit.empty:
                rows.append({"accession": accession, "marker_group": group, "gene_symbol": gene, "present": "no"})
                continue
            row = hit.iloc[0]
            rows.append(
                {
                    "accession": accession,
                    "marker_group": group,
                    "gene_symbol": gene,
                    "present": "yes",
                    "log2fc": fmt(float(row["log2fc"])),
                    "fdr": fmt(float(row["fdr"])),
                    "direction": "up" if float(row["log2fc"]) > 0 else "down",
                    "fdr_lt_0_05": "yes" if float(row["fdr"]) < 0.05 else "no",
                }
            )
    summary: list[dict[str, Any]] = []
    df = pd.DataFrame(rows)
    for (accession, group), sub in df[df["present"].eq("yes")].groupby(["accession", "marker_group"]):
        summary.append(
            {
                "accession": accession,
                "marker_group": group,
                "markers_present": len(sub),
                "markers_up": int(((sub["direction"] == "up")).sum()),
                "markers_up_fdr_lt_0_05": int(((sub["direction"] == "up") & (sub["fdr_lt_0_05"] == "yes")).sum()),
                "median_log2fc": fmt(float(pd.to_numeric(sub["log2fc"]).median())),
            }
        )
    write_tsv(
        OUT_DIR / "bulk_b_plasma_marker_gene_directions.tsv",
        rows,
        ["accession", "marker_group", "gene_symbol", "present", "log2fc", "fdr", "direction", "fdr_lt_0_05"],
    )
    write_tsv(
        OUT_DIR / "bulk_b_plasma_marker_summary.tsv",
        summary,
        ["accession", "marker_group", "markers_present", "markers_up", "markers_up_fdr_lt_0_05", "median_log2fc"],
    )


def main() -> int:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    run_signature_sensitivity()
    run_herb_decoy_null()
    run_bulk_marker_overlap()
    print(f"Wrote robustness outputs to {OUT_DIR.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
