#!/usr/bin/env python3
"""Reviewer-facing lightweight extensions that do not require new downloads.

Outputs:
- scRNA-reference marker-constrained bulk deconvolution proxy
- standard target-overlap baseline comparison for TCM candidate ranking
"""

from __future__ import annotations

import csv
import gzip
import math
from collections import defaultdict
from pathlib import Path
from typing import Any, Iterable

import numpy as np
import pandas as pd
from scipy.optimize import nnls
from scipy.stats import mannwhitneyu, spearmanr

import build_geo_uc_ibd_signature as geo


ROOT = Path(__file__).resolve().parents[1]
OUT_DIR = ROOT / "results" / "m18_reviewer_lightweight_extensions"
RAW_SCRNA = ROOT / "data" / "raw" / "geo" / "GSE125527"
RAW_RECTAL = RAW_SCRNA / "raw_rectal"
GEO_RAW = ROOT / "data" / "raw" / "geo"
METADATA = ROOT / "results" / "m1_uc_ibd_metadata_audit" / "geo_sample_metadata.tsv"

CURATED_IMMUNE_MARKERS = {
    "T": [
        "CD3D",
        "CD3E",
        "CD3G",
        "TRAC",
        "CD2",
        "CD247",
        "IL7R",
        "CD4",
        "CD8A",
        "CD8B",
        "LTB",
        "CCR7",
        "FOXP3",
    ],
    "B": [
        "MS4A1",
        "CD19",
        "CD79A",
        "CD79B",
        "CD74",
        "BANK1",
        "BLK",
        "PAX5",
        "TNFRSF13C",
        "IGKC",
        "IGHM",
        "IGHD",
    ],
    "M/DC": [
        "LYZ",
        "LST1",
        "FCER1G",
        "TYROBP",
        "CD14",
        "FCGR3A",
        "ITGAM",
        "ITGAX",
        "HLA-DRA",
        "HLA-DPA1",
        "C1QA",
        "C1QB",
        "C1QC",
        "S100A8",
        "S100A9",
    ],
    "NK": [
        "NKG7",
        "GNLY",
        "KLRD1",
        "KLRF1",
        "GZMB",
        "PRF1",
        "FGFBP2",
        "XCL1",
        "XCL2",
        "TRDC",
    ],
}


def read_gzip_lines(path: Path) -> Iterable[str]:
    with gzip.open(path, "rt", encoding="utf-8", errors="replace", newline="") as handle:
        yield from handle


def write_tsv(path: Path, rows: list[dict[str, Any]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: "" if row.get(field) is None else row.get(field, "") for field in fields})


def fmt(value: Any, digits: int = 6) -> str:
    try:
        v = float(value)
    except (TypeError, ValueError):
        return "" if value is None else str(value)
    if not math.isfinite(v):
        return ""
    return f"{v:.{digits}g}"


def split_genes(value: Any) -> set[str]:
    if pd.isna(value) or not str(value).strip():
        return set()
    return {part.strip().upper() for part in str(value).split(";") if part.strip()}


def load_metadata_map() -> dict[tuple[str, str, str], dict[str, str]]:
    mapping: dict[tuple[str, str, str], dict[str, str]] = {}
    reader = csv.reader(read_gzip_lines(RAW_SCRNA / "GSE125527_cell_metadata.csv.gz"))
    next(reader)
    for raw in reader:
        if len(raw) < 6:
            continue
        row = {
            "cell_id": raw[0],
            "cell_barcode": raw[1],
            "patient_assignment": raw[2],
            "tissue_assignment": raw[3],
            "disease_assignment": raw[4],
            "celltype": raw[5],
        }
        mapping[(row["patient_assignment"], row["tissue_assignment"], row["cell_barcode"])] = row
    return mapping


def load_patient_id_map() -> dict[str, str]:
    reader = csv.DictReader(read_gzip_lines(RAW_SCRNA / "GSE125527_oldPatientId-newPatientId.csv.gz"))
    return {row["old_id"]: row["new_id"] for row in reader}


def parse_sample_from_filename(path: Path) -> tuple[str, str]:
    parts = path.name.split("_")
    return (parts[1], parts[2]) if len(parts) >= 4 else ("", "")


def build_rectal_scrna_reference() -> tuple[pd.DataFrame, pd.DataFrame]:
    """Return log1p-CPM pseudobulk reference profiles and marker genes."""
    metadata = load_metadata_map()
    patient_map = load_patient_id_map()
    first = next(iter(sorted(RAW_RECTAL.glob("*_cell-gene_UMI_table.tsv.gz"))))
    with gzip.open(first, "rt", encoding="utf-8", errors="replace", newline="") as handle:
        genes = next(csv.reader(handle, delimiter="\t"))[1:]
    genes_upper = [gene.upper() for gene in genes]
    n_genes = len(genes_upper)

    sum_counts: dict[str, np.ndarray] = defaultdict(lambda: np.zeros(n_genes, dtype=np.float64))
    total_counts: dict[str, float] = defaultdict(float)
    n_cells: dict[str, int] = defaultdict(int)

    for path in sorted(RAW_RECTAL.glob("*_cell-gene_UMI_table.tsv.gz")):
        old_patient, tissue = parse_sample_from_filename(path)
        patient = patient_map.get(old_patient, old_patient)
        with gzip.open(path, "rt", encoding="utf-8", errors="replace", newline="") as handle:
            reader = csv.reader(handle, delimiter="\t")
            next(reader)
            for raw in reader:
                if not raw:
                    continue
                barcode = raw[0]
                meta = metadata.get((patient, tissue, barcode))
                if not meta or meta["tissue_assignment"] != "R":
                    continue
                celltype = meta["celltype"]
                if celltype == "unknown":
                    continue
                values = np.array([float(v) if v else 0.0 for v in raw[1:]], dtype=np.float64)
                total = float(values.sum())
                if total <= 0:
                    continue
                sum_counts[celltype] += values
                total_counts[celltype] += total
                n_cells[celltype] += 1

    profiles = {}
    for celltype, counts in sum_counts.items():
        cpm = counts / max(total_counts[celltype], 1.0) * 1_000_000.0
        profiles[celltype] = np.log1p(cpm)
    ref = pd.DataFrame(profiles, index=genes_upper)
    ref = ref.loc[~ref.index.duplicated(keep="first")]

    marker_rows: list[dict[str, Any]] = []
    for celltype in ref.columns:
        own = ref[celltype]
        others = ref.drop(columns=[celltype]).mean(axis=1)
        specificity = own - others
        curated = [gene for gene in CURATED_IMMUNE_MARKERS.get(celltype, []) if gene in ref.index]
        for rank, gene in enumerate(curated, start=1):
            marker_rows.append(
                {
                    "celltype": celltype,
                    "marker_rank": rank,
                    "gene_symbol": gene,
                    "marker_source": "curated_immune_marker",
                    "reference_log1p_cpm": fmt(own.loc[gene]),
                    "specificity_delta_vs_other_celltypes": fmt(specificity.loc[gene]),
                    "n_reference_cells": n_cells[celltype],
                }
            )
    markers = pd.DataFrame(marker_rows)
    return ref, markers


def load_bulk_gene_expression(accession: str) -> tuple[pd.DataFrame, dict[str, str]]:
    meta = geo.read_metadata(METADATA)
    sample_groups = {
        row["gsm"]: geo.assign_group(accession, row)
        for row in meta
        if row.get("accession") == accession and geo.assign_group(accession, row) in {"case", "control"}
    }
    matrix_path = GEO_RAW / accession / f"{accession}_series_matrix.txt.gz"
    platform = geo.DATASETS[accession]["platform"]
    annot = geo.load_platform_annotation(GEO_RAW / "platforms" / f"{platform}.annot.gz")

    sample_ids: list[str] = []
    selected_idx: list[int] = []
    per_gene: dict[str, tuple[float, np.ndarray]] = {}

    for probe_id, values in geo.iter_matrix_rows(matrix_path):
        if probe_id == "HEADER":
            sample_ids = values[1:]
            selected_idx = [idx for idx, gsm in enumerate(sample_ids) if gsm in sample_groups]
            selected_samples = [sample_ids[idx] for idx in selected_idx]
            continue
        if not selected_idx:
            continue
        ann = annot.get(probe_id, {})
        gene = ann.get("gene_symbol", "").upper()
        if not gene:
            continue
        numeric = np.array([geo.safe_float(value) for value in values], dtype=np.float64)[selected_idx]
        if np.isnan(numeric).all():
            continue
        numeric = np.where(np.isnan(numeric), np.nanmedian(numeric), numeric)
        variance = float(np.var(numeric))
        if gene not in per_gene or variance > per_gene[gene][0]:
            per_gene[gene] = (variance, numeric)

    expr = pd.DataFrame({gene: vals for gene, (_, vals) in per_gene.items()}, index=selected_samples).T
    groups = {gsm: sample_groups[gsm] for gsm in selected_samples}
    return expr, groups


def minmax_rows(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy().astype(float)
    mins = out.min(axis=1)
    maxs = out.max(axis=1)
    denom = (maxs - mins).replace(0, np.nan)
    return out.sub(mins, axis=0).div(denom, axis=0).fillna(0.0)


def run_deconvolution_proxy() -> None:
    ref, markers = build_rectal_scrna_reference()
    marker_genes = sorted(set(markers["gene_symbol"]))
    write_tsv(
        OUT_DIR / "scrna_reference_marker_genes.tsv",
        markers.to_dict("records"),
        [
            "celltype",
            "marker_rank",
            "gene_symbol",
            "marker_source",
            "reference_log1p_cpm",
            "specificity_delta_vs_other_celltypes",
            "n_reference_cells",
        ],
    )

    sample_rows: list[dict[str, Any]] = []
    summary_rows: list[dict[str, Any]] = []
    for accession in ["GSE75214", "GSE59071", "GSE87466"]:
        bulk, groups = load_bulk_gene_expression(accession)
        common = sorted(set(marker_genes) & set(ref.index) & set(bulk.index))
        if len(common) < 20:
            continue
        x = minmax_rows(ref.loc[common])
        y = minmax_rows(bulk.loc[common])
        celltypes = list(x.columns)
        x_mat = x.to_numpy(dtype=float)
        if np.allclose(x_mat, 0):
            continue
        for sample in y.columns:
            coef, residual = nnls(x_mat, y[sample].to_numpy(dtype=float))
            total = float(coef.sum())
            fractions = coef / total if total > 0 else np.zeros_like(coef)
            for celltype, fraction in zip(celltypes, fractions):
                sample_rows.append(
                    {
                        "accession": accession,
                        "gsm": sample,
                        "group": groups.get(sample, ""),
                        "celltype": celltype,
                        "fraction_proxy": fmt(float(fraction)),
                        "nnls_residual": fmt(float(residual)),
                        "n_common_marker_genes": len(common),
                    }
                )
    sample_df = pd.DataFrame(sample_rows)
    if not sample_df.empty:
        for (accession, celltype), sub in sample_df.groupby(["accession", "celltype"]):
            case = pd.to_numeric(sub.loc[sub["group"].eq("case"), "fraction_proxy"])
            control = pd.to_numeric(sub.loc[sub["group"].eq("control"), "fraction_proxy"])
            if len(case) < 3 or len(control) < 3:
                continue
            p_value = float(mannwhitneyu(case, control, alternative="two-sided").pvalue)
            summary_rows.append(
                {
                    "accession": accession,
                    "celltype": celltype,
                    "n_case": len(case),
                    "n_control": len(control),
                    "mean_case_fraction_proxy": fmt(case.mean()),
                    "mean_control_fraction_proxy": fmt(control.mean()),
                    "delta_case_minus_control": fmt(case.mean() - control.mean()),
                    "median_case_fraction_proxy": fmt(case.median()),
                    "median_control_fraction_proxy": fmt(control.median()),
                    "mannwhitney_p": fmt(p_value),
                    "n_common_marker_genes": int(sub["n_common_marker_genes"].iloc[0]),
                    "interpretation": "B_lineage_focus" if celltype in {"B", "Plasma"} else "context",
                }
            )
    write_tsv(
        OUT_DIR / "bulk_scrna_reference_deconvolution_proxy_by_sample.tsv",
        sample_rows,
        ["accession", "gsm", "group", "celltype", "fraction_proxy", "nnls_residual", "n_common_marker_genes"],
    )
    write_tsv(
        OUT_DIR / "bulk_scrna_reference_deconvolution_proxy_summary.tsv",
        summary_rows,
        [
            "accession",
            "celltype",
            "n_case",
            "n_control",
            "mean_case_fraction_proxy",
            "mean_control_fraction_proxy",
            "delta_case_minus_control",
            "median_case_fraction_proxy",
            "median_control_fraction_proxy",
            "mannwhitney_p",
            "n_common_marker_genes",
            "interpretation",
        ],
    )


def run_standard_overlap_baseline() -> None:
    m9 = pd.read_csv(ROOT / "results/m9_target_cell_bridge/candidate_target_cell_bridge_summary.tsv", sep="\t")
    disease_edges = pd.read_csv(ROOT / "results/m9_target_cell_bridge/disease_context_herb_target_edges.tsv", sep="\t")
    bulk_up = set(disease_edges.loc[disease_edges["in_bulk_uc_ibd_up150"].eq("yes"), "gene_symbol"].dropna().str.upper())
    bulk_down = set(disease_edges.loc[disease_edges["in_bulk_uc_ibd_down150"].eq("yes"), "gene_symbol"].dropna().str.upper())
    disease_module = bulk_up | bulk_down

    rows: list[dict[str, Any]] = []
    for _, row in m9.iterrows():
        targets = split_genes(row["disease_context_unique_targets"])
        overlap = targets & disease_module
        up_overlap = targets & bulk_up
        baseline_score = len(overlap)
        disease_cell_score = float(row["bridge_score"])
        rows.append(
            {
                "ingredient_id": row["ingredient_id"],
                "ingredient_name": "Berberine" if row["ingredient_name"] == "Berberime" else row["ingredient_name"],
                "n_disease_context_targets": int(row["n_disease_context_unique_targets"]),
                "standard_overlap_hits": baseline_score,
                "standard_overlap_up_hits": len(up_overlap),
                "standard_overlap_genes": ";".join(sorted(overlap)),
                "disease_cell_aware_score": fmt(disease_cell_score),
                "b_cell_increased_targets": int(row["n_b_cell_disease_increased_targets"]),
                "mdc_increased_targets": int(row["n_mdc_disease_increased_targets"]),
            }
        )
    rows.sort(key=lambda r: (-int(r["standard_overlap_hits"]), -int(r["n_disease_context_targets"]), r["ingredient_name"]))
    for idx, row in enumerate(rows, start=1):
        row["standard_overlap_rank"] = idx
    disease_sorted = sorted(rows, key=lambda r: (-float(r["disease_cell_aware_score"]), r["ingredient_name"]))
    disease_rank = {row["ingredient_id"]: idx for idx, row in enumerate(disease_sorted, start=1)}
    for row in rows:
        row["disease_cell_aware_rank"] = disease_rank[row["ingredient_id"]]
        row["rank_shift_overlap_minus_disease_cell"] = int(row["standard_overlap_rank"]) - int(row["disease_cell_aware_rank"])

    rho, rho_p = spearmanr(
        [int(row["standard_overlap_rank"]) for row in rows],
        [int(row["disease_cell_aware_rank"]) for row in rows],
    )
    cur = next(row for row in rows if row["ingredient_name"] == "Curcumin")
    summary = [
        {
            "analysis": "standard_target_overlap_baseline_vs_disease_cell_aware_ranking",
            "n_candidates": len(rows),
            "curcumin_standard_overlap_rank": cur["standard_overlap_rank"],
            "curcumin_disease_cell_aware_rank": cur["disease_cell_aware_rank"],
            "curcumin_standard_overlap_hits": cur["standard_overlap_hits"],
            "curcumin_b_cell_increased_targets": cur["b_cell_increased_targets"],
            "rank_spearman_rho": fmt(float(rho)),
            "rank_spearman_p": fmt(float(rho_p)),
            "interpretation": "Disease-cell-aware ranking preserves Curcumin as lead while exposing how a standard overlap baseline is driven by generic inflammatory overlap.",
        }
    ]
    fields = [
        "standard_overlap_rank",
        "disease_cell_aware_rank",
        "rank_shift_overlap_minus_disease_cell",
        "ingredient_id",
        "ingredient_name",
        "n_disease_context_targets",
        "standard_overlap_hits",
        "standard_overlap_up_hits",
        "standard_overlap_genes",
        "disease_cell_aware_score",
        "b_cell_increased_targets",
        "mdc_increased_targets",
    ]
    write_tsv(OUT_DIR / "standard_overlap_baseline_candidate_ranking.tsv", rows, fields)
    write_tsv(
        OUT_DIR / "standard_overlap_baseline_summary.tsv",
        summary,
        [
            "analysis",
            "n_candidates",
            "curcumin_standard_overlap_rank",
            "curcumin_disease_cell_aware_rank",
            "curcumin_standard_overlap_hits",
            "curcumin_b_cell_increased_targets",
            "rank_spearman_rho",
            "rank_spearman_p",
            "interpretation",
        ],
    )


def main() -> int:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    run_deconvolution_proxy()
    run_standard_overlap_baseline()
    print(f"Wrote reviewer lightweight extension outputs to {OUT_DIR.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
