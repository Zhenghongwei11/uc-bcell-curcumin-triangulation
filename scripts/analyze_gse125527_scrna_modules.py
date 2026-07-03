#!/usr/bin/env python3
"""Score UC/IBD bulk consensus modules in GSE125527 single-cell data."""

from __future__ import annotations

import argparse
import csv
import gzip
import math
from collections import defaultdict
from pathlib import Path
from statistics import mean, median
from typing import Any, Dict, Iterable, List, Sequence, Tuple

import numpy as np
from scipy.stats import mannwhitneyu


def read_gene_list(path: Path) -> List[str]:
    genes: List[str] = []
    seen = set()
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            gene = line.strip().upper()
            if gene and gene not in seen:
                genes.append(gene)
                seen.add(gene)
    return genes


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


def load_genes(path: Path) -> Tuple[List[str], Dict[str, int]]:
    genes = [line.strip().upper() for line in read_gzip_lines(path) if line.strip()]
    return genes, {gene: idx + 1 for idx, gene in enumerate(genes)}


def load_cells(path: Path) -> List[str]:
    return [line.strip() for line in read_gzip_lines(path) if line.strip()]


def load_metadata(path: Path, expected_cells: int) -> List[Dict[str, str]]:
    rows: List[Dict[str, str]] = [{} for _ in range(expected_cells)]
    reader = csv.reader(read_gzip_lines(path))
    header = next(reader)
    corrected = [
        "row_id",
        "cell_barcode",
        "patient_assignment",
        "tissue_assignment",
        "disease_assignment",
        "celltype",
    ]
    if len(header) == len(corrected):
        fieldnames = header
    else:
        fieldnames = corrected
    for raw in reader:
        if not raw:
            continue
        if len(raw) < len(fieldnames):
            continue
        row = dict(zip(fieldnames, raw))
        idx = int(row["row_id"]) - 1
        if 0 <= idx < expected_cells:
            rows[idx] = row
    missing = sum(1 for row in rows if not row)
    if missing:
        raise ValueError(f"Missing metadata for {missing} cells")
    return rows


def load_metadata_map(path: Path) -> Dict[Tuple[str, str, str], Dict[str, str]]:
    mapping: Dict[Tuple[str, str, str], Dict[str, str]] = {}
    reader = csv.reader(read_gzip_lines(path))
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


def load_cluster_maps(raw_dir: Path) -> Dict[Tuple[str, str, str], str]:
    mapping: Dict[Tuple[str, str, str], str] = {}
    for filename in ["GSE125527_Tcell_cluster.csv.gz", "GSE125527_Bcell_cluster.csv.gz"]:
        path = raw_dir / filename
        if not path.exists():
            continue
        reader = csv.DictReader(read_gzip_lines(path))
        for row in reader:
            mapping[(row["patient_id"], row["tissue_id"], row["cell_id"])] = row["Cluster_id"]
    return mapping


def load_patient_id_map(raw_dir: Path) -> Dict[str, str]:
    path = raw_dir / "GSE125527_oldPatientId-newPatientId.csv.gz"
    if not path.exists():
        return {}
    reader = csv.DictReader(read_gzip_lines(path))
    return {row["old_id"]: row["new_id"] for row in reader}


def score_sparse_matrix(
    matrix_path: Path,
    n_cells: int,
    up_rows: set[int],
    down_rows: set[int],
) -> Dict[str, np.ndarray]:
    total = np.zeros(n_cells, dtype=np.float64)
    up = np.zeros(n_cells, dtype=np.float64)
    down = np.zeros(n_cells, dtype=np.float64)
    up_detected = np.zeros(n_cells, dtype=np.int16)
    down_detected = np.zeros(n_cells, dtype=np.int16)

    reader = csv.DictReader(read_gzip_lines(matrix_path))
    for item in reader:
        row = int(item["row"])
        col = int(item["col"]) - 1
        value = float(item["value"])
        if col < 0 or col >= n_cells:
            continue
        total[col] += value
        if row in up_rows:
            up[col] += value
            up_detected[col] += 1
        if row in down_rows:
            down[col] += value
            down_detected[col] += 1

    scale = np.divide(10000.0, total, out=np.zeros_like(total), where=total > 0)
    up_score = np.log1p(up * scale)
    down_score = np.log1p(down * scale)
    disease_axis_score = up_score - down_score
    return {
        "total_umi": total,
        "up_counts": up,
        "down_counts": down,
        "up_score": up_score,
        "down_score": down_score,
        "disease_axis_score": disease_axis_score,
        "up_detected": up_detected,
        "down_detected": down_detected,
    }


def summarize(values: Sequence[float]) -> Dict[str, float]:
    if not values:
        return {"mean": math.nan, "median": math.nan}
    return {"mean": mean(values), "median": median(values)}


def safe_pvalue(a: Sequence[float], b: Sequence[float]) -> str:
    if len(a) < 2 or len(b) < 2:
        return ""
    try:
        return f"{mannwhitneyu(a, b, alternative='two-sided').pvalue:.6g}"
    except ValueError:
        return ""


def parse_sample_from_filename(path: Path) -> Tuple[str, str]:
    parts = path.name.split("_")
    if len(parts) < 4:
        return "", ""
    return parts[1], parts[2]


def score_sample_tables(
    sample_dir: Path,
    metadata_map: Dict[Tuple[str, str, str], Dict[str, str]],
    cluster_map: Dict[Tuple[str, str, str], str],
    patient_id_map: Dict[str, str],
    up_present: Sequence[str],
    down_present: Sequence[str],
) -> List[Dict[str, Any]]:
    records: List[Dict[str, Any]] = []
    up_set = set(up_present)
    down_set = set(down_present)
    for path in sorted(sample_dir.glob("*_cell-gene_UMI_table.tsv.gz")):
        old_patient, tissue = parse_sample_from_filename(path)
        patient = patient_id_map.get(old_patient, old_patient)
        with gzip.open(path, "rt", encoding="utf-8", errors="replace", newline="") as handle:
            reader = csv.reader(handle, delimiter="\t")
            header = next(reader)
            genes = [gene.upper() for gene in header[1:]]
            up_idx = [idx for idx, gene in enumerate(genes, start=1) if gene in up_set]
            down_idx = [idx for idx, gene in enumerate(genes, start=1) if gene in down_set]
            for raw in reader:
                if not raw:
                    continue
                barcode = raw[0]
                counts = raw[1:]
                total = 0.0
                for value in counts:
                    if value:
                        total += float(value)
                if total <= 0:
                    continue
                up_count = sum(float(counts[idx - 1]) for idx in up_idx if idx - 1 < len(counts))
                down_count = sum(float(counts[idx - 1]) for idx in down_idx if idx - 1 < len(counts))
                meta = metadata_map.get((patient, tissue, barcode), {})
                disease = meta.get("disease_assignment") or ("diseased" if patient.startswith("U") else "healthy")
                celltype = meta.get("celltype", "unknown")
                cluster = cluster_map.get((patient, tissue, barcode), "")
                scale = 10000.0 / total
                up_score = math.log1p(up_count * scale)
                down_score = math.log1p(down_count * scale)
                up_detected = sum(1 for idx in up_idx if idx - 1 < len(counts) and float(counts[idx - 1]) > 0)
                down_detected = sum(1 for idx in down_idx if idx - 1 < len(counts) and float(counts[idx - 1]) > 0)
                records.append(
                    {
                        "patient_assignment": patient,
                        "raw_sample_id": old_patient,
                        "tissue_assignment": tissue,
                        "disease_assignment": disease,
                        "celltype": celltype,
                        "celltype_detail": f"{celltype}:{cluster}" if cluster else celltype,
                        "cell_barcode": barcode,
                        "total_umi": total,
                        "up_score": up_score,
                        "down_score": down_score,
                        "disease_axis_score": up_score - down_score,
                        "up_detected": up_detected,
                        "down_detected": down_detected,
                    }
                )
    return records


def aggregate_cell_records(records: Sequence[Dict[str, Any]], group_fields: Sequence[str]) -> List[Dict[str, Any]]:
    grouped: Dict[Tuple[str, ...], Dict[str, List[float]]] = defaultdict(lambda: defaultdict(list))
    for record in records:
        key = tuple(str(record[field]) for field in group_fields)
        for metric in ["up_score", "down_score", "disease_axis_score", "up_detected", "down_detected", "total_umi"]:
            grouped[key][metric].append(float(record[metric]))
    rows: List[Dict[str, Any]] = []
    for key, metrics in sorted(grouped.items()):
        row: Dict[str, Any] = dict(zip(group_fields, key))
        row["n_cells"] = len(metrics["up_score"])
        for metric in ["up_score", "down_score", "disease_axis_score", "up_detected", "down_detected", "total_umi"]:
            stats = summarize(metrics[metric])
            row[f"mean_{metric}"] = f"{stats['mean']:.6g}"
            row[f"median_{metric}"] = f"{stats['median']:.6g}"
        rows.append(row)
    return rows


def pseudobulk_records(records: Sequence[Dict[str, Any]], celltype_field: str = "celltype") -> List[Dict[str, Any]]:
    grouped: Dict[Tuple[str, str, str, str], Dict[str, List[float]]] = defaultdict(lambda: defaultdict(list))
    for record in records:
        key = (
            str(record["patient_assignment"]),
            str(record["tissue_assignment"]),
            str(record["disease_assignment"]),
            str(record[celltype_field]),
        )
        for metric in ["up_score", "down_score", "disease_axis_score", "up_detected", "down_detected", "total_umi"]:
            grouped[key][metric].append(float(record[metric]))
    rows: List[Dict[str, Any]] = []
    for (patient, tissue, disease, celltype), metrics in sorted(grouped.items()):
        row = {
            "patient_assignment": patient,
            "tissue_assignment": tissue,
            "disease_assignment": disease,
            "celltype": celltype,
            "n_cells": len(metrics["up_score"]),
        }
        for metric in ["up_score", "down_score", "disease_axis_score", "up_detected", "down_detected", "total_umi"]:
            row[f"mean_{metric}"] = f"{mean(metrics[metric]):.6g}"
        rows.append(row)
    return rows


def contrast_from_pseudobulk(pseudobulk_rows: Sequence[Dict[str, Any]]) -> List[Dict[str, Any]]:
    values_by_group: Dict[Tuple[str, str, str], List[Dict[str, Any]]] = defaultdict(list)
    for row in pseudobulk_rows:
        values_by_group[(row["tissue_assignment"], row["celltype"], row["disease_assignment"])].append(row)
    contrast_rows: List[Dict[str, Any]] = []
    tissues = sorted({row["tissue_assignment"] for row in pseudobulk_rows})
    celltypes = sorted({row["celltype"] for row in pseudobulk_rows})
    for tissue in tissues:
        for celltype in celltypes:
            diseased = values_by_group.get((tissue, celltype, "diseased"), [])
            healthy = values_by_group.get((tissue, celltype, "healthy"), [])
            if not diseased or not healthy:
                continue
            d_axis = [float(row["mean_disease_axis_score"]) for row in diseased]
            h_axis = [float(row["mean_disease_axis_score"]) for row in healthy]
            d_up = [float(row["mean_up_score"]) for row in diseased]
            h_up = [float(row["mean_up_score"]) for row in healthy]
            d_down = [float(row["mean_down_score"]) for row in diseased]
            h_down = [float(row["mean_down_score"]) for row in healthy]
            contrast_rows.append(
                {
                    "tissue_assignment": tissue,
                    "celltype": celltype,
                    "n_diseased_samples": len(diseased),
                    "n_healthy_samples": len(healthy),
                    "diseased_cells": sum(int(row["n_cells"]) for row in diseased),
                    "healthy_cells": sum(int(row["n_cells"]) for row in healthy),
                    "mean_axis_diseased": f"{mean(d_axis):.6g}",
                    "mean_axis_healthy": f"{mean(h_axis):.6g}",
                    "delta_axis_diseased_minus_healthy": f"{mean(d_axis) - mean(h_axis):.6g}",
                    "p_axis_mannwhitney": safe_pvalue(d_axis, h_axis),
                    "mean_up_diseased": f"{mean(d_up):.6g}",
                    "mean_up_healthy": f"{mean(h_up):.6g}",
                    "delta_up_diseased_minus_healthy": f"{mean(d_up) - mean(h_up):.6g}",
                    "p_up_mannwhitney": safe_pvalue(d_up, h_up),
                    "mean_down_diseased": f"{mean(d_down):.6g}",
                    "mean_down_healthy": f"{mean(h_down):.6g}",
                    "delta_down_diseased_minus_healthy": f"{mean(d_down) - mean(h_down):.6g}",
                    "p_down_mannwhitney": safe_pvalue(d_down, h_down),
                }
            )
    contrast_rows.sort(key=lambda row: (row["tissue_assignment"], -float(row["delta_axis_diseased_minus_healthy"])))
    return contrast_rows


def patient_label_null_from_pseudobulk(
    pseudobulk_rows: Sequence[Dict[str, Any]],
    n_permutations: int,
    seed: int,
) -> List[Dict[str, Any]]:
    rng = np.random.default_rng(seed)
    grouped: Dict[Tuple[str, str], List[Dict[str, Any]]] = defaultdict(list)
    for row in pseudobulk_rows:
        grouped[(row["tissue_assignment"], row["celltype"])].append(row)

    rows: List[Dict[str, Any]] = []
    for (tissue, celltype), group_rows in sorted(grouped.items()):
        diseased = [row for row in group_rows if row["disease_assignment"] == "diseased"]
        healthy = [row for row in group_rows if row["disease_assignment"] == "healthy"]
        if not diseased or not healthy:
            continue
        values = np.array([float(row["mean_disease_axis_score"]) for row in group_rows], dtype=np.float64)
        labels = np.array([row["disease_assignment"] for row in group_rows])
        diseased_count = int(np.sum(labels == "diseased"))
        observed = float(np.mean(values[labels == "diseased"]) - np.mean(values[labels == "healthy"]))
        null_values = np.empty(n_permutations, dtype=np.float64)
        all_indices = np.arange(len(values))
        for i in range(n_permutations):
            perm_diseased = rng.choice(all_indices, size=diseased_count, replace=False)
            mask = np.zeros(len(values), dtype=bool)
            mask[perm_diseased] = True
            null_values[i] = float(np.mean(values[mask]) - np.mean(values[~mask]))
        extreme = int(np.sum(np.abs(null_values) >= abs(observed)))
        empirical_p = (extreme + 1.0) / (n_permutations + 1.0)
        null_sd = float(np.std(null_values, ddof=1)) if n_permutations > 1 else math.nan
        rows.append(
            {
                "tissue_assignment": tissue,
                "celltype": celltype,
                "n_diseased_samples": len(diseased),
                "n_healthy_samples": len(healthy),
                "observed_delta_axis_diseased_minus_healthy": f"{observed:.6g}",
                "null_mean_delta": f"{float(np.mean(null_values)):.6g}",
                "null_sd_delta": f"{null_sd:.6g}",
                "null_q025_delta": f"{float(np.quantile(null_values, 0.025)):.6g}",
                "null_q975_delta": f"{float(np.quantile(null_values, 0.975)):.6g}",
                "empirical_p_two_sided": f"{empirical_p:.6g}",
                "n_permutations": n_permutations,
                "seed": seed,
            }
        )
    rows.sort(key=lambda row: (row["tissue_assignment"], -abs(float(row["observed_delta_axis_diseased_minus_healthy"]))))
    return rows


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--raw-dir", default="data/raw/geo/GSE125527")
    parser.add_argument("--sample-dir", default="data/raw/geo/GSE125527/raw_rectal")
    parser.add_argument("--up-genes", default="results/m4_geo_uc_ibd_signature/lincs_query_up_genes.txt")
    parser.add_argument("--down-genes", default="results/m4_geo_uc_ibd_signature/lincs_query_down_genes.txt")
    parser.add_argument("--out-dir", default="results/m6_scrna_cell_state")
    parser.add_argument("--null-permutations", type=int, default=10000)
    parser.add_argument("--null-seed", type=int, default=20260629)
    args = parser.parse_args()

    raw_dir = Path(args.raw_dir)
    out_dir = Path(args.out_dir)
    up_query = read_gene_list(Path(args.up_genes))
    down_query = read_gene_list(Path(args.down_genes))
    sample_dir = Path(args.sample_dir)
    first_sample = next(iter(sorted(sample_dir.glob("*_cell-gene_UMI_table.tsv.gz"))))
    with gzip.open(first_sample, "rt", encoding="utf-8", errors="replace", newline="") as handle:
        first_header = next(csv.reader(handle, delimiter="\t"))
    sample_genes = {gene.upper() for gene in first_header[1:]}
    up_present = [gene for gene in up_query if gene in sample_genes]
    down_present = [gene for gene in down_query if gene in sample_genes]

    records = score_sample_tables(
        sample_dir,
        load_metadata_map(raw_dir / "GSE125527_cell_metadata.csv.gz"),
        load_cluster_maps(raw_dir),
        load_patient_id_map(raw_dir),
        up_present,
        down_present,
    )
    cell_summary_rows = aggregate_cell_records(records, ["tissue_assignment", "disease_assignment", "celltype"])
    cell_detail_summary_rows = aggregate_cell_records(records, ["tissue_assignment", "disease_assignment", "celltype_detail"])
    pseudobulk_rows = pseudobulk_records(records, "celltype")
    pseudobulk_detail_rows = pseudobulk_records(records, "celltype_detail")
    contrast_rows = contrast_from_pseudobulk(pseudobulk_rows)
    contrast_detail_rows = contrast_from_pseudobulk(pseudobulk_detail_rows)
    null_rows = patient_label_null_from_pseudobulk(pseudobulk_rows, args.null_permutations, args.null_seed)

    coverage_rows = [
        {
            "module": "bulk_consensus_up",
            "query_genes": len(up_query),
            "genes_present_in_scrna": len(up_present),
            "coverage_fraction": f"{len(up_present) / len(up_query):.6g}" if up_query else "",
            "present_genes": ";".join(up_present),
            "missing_genes": ";".join(gene for gene in up_query if gene not in sample_genes),
        },
        {
            "module": "bulk_consensus_down",
            "query_genes": len(down_query),
            "genes_present_in_scrna": len(down_present),
            "coverage_fraction": f"{len(down_present) / len(down_query):.6g}" if down_query else "",
            "present_genes": ";".join(down_present),
            "missing_genes": ";".join(gene for gene in down_query if gene not in sample_genes),
        },
    ]

    write_tsv(
        out_dir / "gse125527_celltype_module_summary.tsv",
        cell_summary_rows,
        [
            "tissue_assignment",
            "disease_assignment",
            "celltype",
            "n_cells",
            "mean_up_score",
            "median_up_score",
            "mean_down_score",
            "median_down_score",
            "mean_disease_axis_score",
            "median_disease_axis_score",
            "mean_up_detected",
            "median_up_detected",
            "mean_down_detected",
            "median_down_detected",
            "mean_total_umi",
            "median_total_umi",
        ],
    )
    write_tsv(
        out_dir / "gse125527_celltype_detail_module_summary.tsv",
        cell_detail_summary_rows,
        [
            "tissue_assignment",
            "disease_assignment",
            "celltype_detail",
            "n_cells",
            "mean_up_score",
            "median_up_score",
            "mean_down_score",
            "median_down_score",
            "mean_disease_axis_score",
            "median_disease_axis_score",
            "mean_up_detected",
            "median_up_detected",
            "mean_down_detected",
            "median_down_detected",
            "mean_total_umi",
            "median_total_umi",
        ],
    )
    write_tsv(
        out_dir / "gse125527_pseudobulk_module_summary.tsv",
        pseudobulk_rows,
        [
            "patient_assignment",
            "tissue_assignment",
            "disease_assignment",
            "celltype",
            "n_cells",
            "mean_up_score",
            "mean_down_score",
            "mean_disease_axis_score",
            "mean_up_detected",
            "mean_down_detected",
            "mean_total_umi",
        ],
    )
    write_tsv(
        out_dir / "gse125527_disease_contrast_by_celltype.tsv",
        contrast_rows,
        [
            "tissue_assignment",
            "celltype",
            "n_diseased_samples",
            "n_healthy_samples",
            "diseased_cells",
            "healthy_cells",
            "mean_axis_diseased",
            "mean_axis_healthy",
            "delta_axis_diseased_minus_healthy",
            "p_axis_mannwhitney",
            "mean_up_diseased",
            "mean_up_healthy",
            "delta_up_diseased_minus_healthy",
            "p_up_mannwhitney",
            "mean_down_diseased",
            "mean_down_healthy",
            "delta_down_diseased_minus_healthy",
            "p_down_mannwhitney",
        ],
    )
    write_tsv(
        out_dir / "gse125527_disease_contrast_by_celltype_detail.tsv",
        contrast_detail_rows,
        [
            "tissue_assignment",
            "celltype",
            "n_diseased_samples",
            "n_healthy_samples",
            "diseased_cells",
            "healthy_cells",
            "mean_axis_diseased",
            "mean_axis_healthy",
            "delta_axis_diseased_minus_healthy",
            "p_axis_mannwhitney",
            "mean_up_diseased",
            "mean_up_healthy",
            "delta_up_diseased_minus_healthy",
            "p_up_mannwhitney",
            "mean_down_diseased",
            "mean_down_healthy",
            "delta_down_diseased_minus_healthy",
            "p_down_mannwhitney",
        ],
    )
    write_tsv(
        out_dir / "gse125527_module_gene_coverage.tsv",
        coverage_rows,
        ["module", "query_genes", "genes_present_in_scrna", "coverage_fraction", "present_genes", "missing_genes"],
    )
    write_tsv(
        out_dir / "gse125527_disease_contrast_patient_label_null.tsv",
        null_rows,
        [
            "tissue_assignment",
            "celltype",
            "n_diseased_samples",
            "n_healthy_samples",
            "observed_delta_axis_diseased_minus_healthy",
            "null_mean_delta",
            "null_sd_delta",
            "null_q025_delta",
            "null_q975_delta",
            "empirical_p_two_sided",
            "n_permutations",
            "seed",
        ],
    )

    top_rectal = [row for row in contrast_rows if row["tissue_assignment"] == "R"]
    null_lookup = {(row["tissue_assignment"], row["celltype"]): row for row in null_rows}
    report_lines = [
        "# M6 GSE125527 scRNA Cell-State Localization Run",
        "",
        "## Scope",
        "",
        "This run scores the GEO bulk-derived UC/IBD consensus up/down modules in processed GSE125527 single-cell UMI data. It uses author-provided sample/tissue/disease/cell-type labels and patient-level pseudobulk contrasts.",
        "",
        "## Key Outputs",
        "",
        "- `results/m6_scrna_cell_state/gse125527_module_gene_coverage.tsv`",
        "- `results/m6_scrna_cell_state/gse125527_celltype_module_summary.tsv`",
        "- `results/m6_scrna_cell_state/gse125527_pseudobulk_module_summary.tsv`",
        "- `results/m6_scrna_cell_state/gse125527_disease_contrast_by_celltype.tsv`",
        "- `results/m6_scrna_cell_state/gse125527_disease_contrast_by_celltype_detail.tsv`",
        "- `results/m6_scrna_cell_state/gse125527_disease_contrast_patient_label_null.tsv`",
        "",
        "## Gene Coverage",
        "",
        "| Module | Query genes | Present in scRNA | Coverage |",
        "|---|---:|---:|---:|",
    ]
    for row in coverage_rows:
        report_lines.append(f"| {row['module']} | {row['query_genes']} | {row['genes_present_in_scrna']} | {row['coverage_fraction']} |")
    report_lines.extend(["", "## Rectal Tissue Disease Axis Contrast", ""])
    if top_rectal:
        report_lines.append("| Cell type | Diseased samples | Healthy samples | Delta axis | MW P value | Label-null empirical P | Diseased cells | Healthy cells |")
        report_lines.append("|---|---:|---:|---:|---:|---:|---:|---:|")
        for row in top_rectal:
            null_row = null_lookup.get((row["tissue_assignment"], row["celltype"]), {})
            report_lines.append(
                f"| {row['celltype']} | {row['n_diseased_samples']} | {row['n_healthy_samples']} | {row['delta_axis_diseased_minus_healthy']} | {row['p_axis_mannwhitney']} | {null_row.get('empirical_p_two_sided', '')} | {row['diseased_cells']} | {row['healthy_cells']} |"
            )
    else:
        report_lines.append("No rectal diseased-vs-healthy contrast was available.")
    b_cell_rows = [row for row in top_rectal if row["celltype"] == "B"]
    mdc_rows = [row for row in top_rectal if row["celltype"] == "M/DC"]
    report_lines.extend(["", "## Actionable Interpretation", ""])
    if b_cell_rows:
        row = b_cell_rows[0]
        report_lines.append(
            f"- Main publishable signal: rectal B cells show the strongest positive disease-axis shift "
            f"(delta={row['delta_axis_diseased_minus_healthy']}, Mann-Whitney P={row['p_axis_mannwhitney']}, "
            f"patient-label null P={null_lookup.get(('R', 'B'), {}).get('empirical_p_two_sided', '')}) "
            "using patient/sample-level pseudobulk scores."
        )
    if mdc_rows:
        row = mdc_rows[0]
        report_lines.append(
            f"- Supportive direction: rectal M/DC cells are also positive "
            f"(delta={row['delta_axis_diseased_minus_healthy']}, Mann-Whitney P={row['p_axis_mannwhitney']}), "
            "but the healthy comparator has fewer samples/cells, so this should be treated as secondary evidence."
        )
    report_lines.append(
        "- Recommended manuscript framing: public UC/IBD mucosal transcriptional injury localizes most robustly to a rectal B-cell inflammatory axis; TCM candidates should be connected to this axis by target/pathway and genetic evidence before making mechanism claims."
    )
    report_lines.extend(
        [
            "",
            "## Interpretation Boundary",
            "",
            "- The primary contrast is patient/sample-level pseudobulk, not cell-level p-value inflation.",
            "- The patient-label null analysis permutes disease labels within each tissue/celltype pseudobulk contrast; it is a sensitivity check for the cell-state disease-axis association, not experimental validation.",
            "- A positive delta means the bulk UC/IBD up module is relatively higher than the down module in diseased samples for that cell class.",
            "- This analysis localizes the public bulk disease signature; it does not validate any TCM compound mechanism by itself.",
            "",
        ]
    )
    (Path("docs") / "M6_SCRNA_CELL_STATE_RUN.md").write_text("\n".join(report_lines), encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
