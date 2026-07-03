#!/usr/bin/env python3
"""Analyze GSE182270 as a small independent B-lineage UC replication dataset.

The script streams filtered 10x matrices from the downloaded GEO tar without
fully extracting all sample folders. It summarizes sample-level B/plasma and
disease-axis module scores for UC inflamed mucosa versus healthy colon.
"""

from __future__ import annotations

import csv
import gzip
import io
import math
import tarfile
from collections import defaultdict
from pathlib import Path
from typing import Any, Iterable

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu


ROOT = Path(__file__).resolve().parents[1]
RAW_TAR = ROOT / "data" / "raw" / "geo" / "GSE182270" / "GSE182270_RAW.tar"
OUT_DIR = ROOT / "results" / "m19_gse182270_bcell_replication"


MODULES = {
    "naive_b_marker": ["MS4A1", "CD19", "CD79A", "CD79B", "BANK1", "TCL1A", "IGHD", "IGHM", "CD74"],
    "plasma_cell_marker": ["MZB1", "XBP1", "JCHAIN", "SDC1", "PRDM1", "TNFRSF17", "IGHG1", "IGHA1", "IGKC", "DERL3"],
    "igg_plasma_marker": ["IGHG1", "IGHG3", "JCHAIN", "MZB1", "XBP1", "PRDM1"],
    "bcl6_blnk_syk_axis": ["BCL6", "BLNK", "SYK"],
    "curcumin_bcell_context_targets": ["BCL6", "BLNK", "SYK", "IL1B", "TNF", "IL7", "JAK1", "PIAS1", "STAT5A"],
}

SAMPLE_STATUS = {
    "GSM5525955": {
        "condition": "UC_inflamed",
        "geo_status": "Inflamed",
        "geo_disease": "Ulcerative colitis",
        "interpretation_note": "UC inflamed by GEO status, disease, description, and file name.",
    },
    "GSM5525956": {
        "condition": "UC_inflamed",
        "geo_status": "Inflamed",
        "geo_disease": "Ulcerative colitis",
        "interpretation_note": "UC inflamed by GEO status, disease, description, and file name.",
    },
    "GSM5525957": {
        "condition": "UC_inflamed",
        "geo_status": "Inflamed",
        "geo_disease": "Ulcerative colitis",
        "interpretation_note": "UC inflamed by GEO status, disease, description, and file name.",
    },
    "GSM5525958": {
        "condition": "UC_inflamed",
        "geo_status": "Healthy control",
        "geo_disease": "Ulcerative colitis",
        "interpretation_note": "Assigned to UC inflamed because GEO disease, description, overall design, and file name support UC; GEO status field is internally discordant.",
    },
    "GSM5525959": {
        "condition": "healthy_control",
        "geo_status": "non inflamed",
        "geo_disease": "Healthy control",
        "interpretation_note": "Assigned to healthy/noninflamed by GEO disease, status, and description despite UC-like file name.",
    },
    "GSM5525960": {
        "condition": "healthy_control",
        "geo_status": "non inflamed",
        "geo_disease": "Healthy control",
        "interpretation_note": "Healthy/noninflamed by GEO disease, status, description, and file name.",
    },
    "GSM5525961": {
        "condition": "healthy_control",
        "geo_status": "non inflamed",
        "geo_disease": "Healthy control",
        "interpretation_note": "Healthy/noninflamed by GEO disease, status, description, and file name.",
    },
    "GSM5525962": {
        "condition": "healthy_control",
        "geo_status": "non inflamed",
        "geo_disease": "Healthy control",
        "interpretation_note": "Healthy/noninflamed by GEO disease, status, description, and file name.",
    },
    "GSM5525963": {
        "condition": "UC_inflamed",
        "geo_status": "Inflamed",
        "geo_disease": "Ulcerative colitis",
        "interpretation_note": "UC inflamed by GEO status, disease, description, and file name.",
    },
}


def fmt(value: Any, digits: int = 6) -> str:
    try:
        v = float(value)
    except (TypeError, ValueError):
        return "" if value is None else str(value)
    if not math.isfinite(v):
        return ""
    return f"{v:.{digits}g}"


def write_tsv(path: Path, rows: list[dict[str, Any]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: "" if row.get(field) is None else row.get(field, "") for field in fields})


def parse_sample(name: str) -> dict[str, str]:
    stem = name.replace(".tar.gz", "")
    gsm, sample = stem.split("_", 1)
    status = SAMPLE_STATUS.get(gsm, {})
    condition = status.get("condition") or ("UC_inflamed" if sample.upper().startswith("UC") else "healthy_control")
    return {
        "gsm": gsm,
        "sample_id": sample,
        "condition": condition,
        "geo_status": status.get("geo_status", ""),
        "geo_disease": status.get("geo_disease", ""),
        "sample_interpretation_note": status.get("interpretation_note", ""),
    }


def find_member(tar: tarfile.TarFile, suffixes: list[str]) -> tarfile.TarInfo:
    lowered = [suffix.lower() for suffix in suffixes]
    for member in tar.getmembers():
        name = member.name.lower()
        if any(name.endswith(suffix) for suffix in lowered):
            return member
    raise FileNotFoundError(f"none of the requested suffixes found: {suffixes}")


def read_text_member(tar: tarfile.TarFile, suffixes: list[str]) -> list[str]:
    member = find_member(tar, suffixes)
    handle = tar.extractfile(member)
    if handle is None:
        return []
    data = handle.read()
    if member.name.lower().endswith(".gz"):
        data = gzip.decompress(data)
    return data.decode("utf-8", errors="replace").splitlines()


def iter_matrix_entries(tar: tarfile.TarFile, suffixes: list[str]) -> tuple[tuple[int, int, int], Iterable[tuple[int, int, float]]]:
    member = find_member(tar, suffixes)
    handle = tar.extractfile(member)
    if handle is None:
        raise RuntimeError(f"missing matrix member {suffixes}")
    if member.name.lower().endswith(".gz"):
        text = io.TextIOWrapper(gzip.GzipFile(fileobj=handle), encoding="utf-8", errors="replace")
    else:
        text = io.TextIOWrapper(handle, encoding="utf-8", errors="replace")

    shape: tuple[int, int, int] | None = None

    def generator() -> Iterable[tuple[int, int, float]]:
        nonlocal shape
        for line in text:
            if not line.strip() or line.startswith("%"):
                continue
            if shape is None:
                a, b, c = line.strip().split()[:3]
                shape = (int(a), int(b), int(c))
                continue
            i, j, value = line.strip().split()[:3]
            yield int(i) - 1, int(j) - 1, float(value)

    entries = generator()
    # Prime the generator so shape is available to callers.
    first_entries: list[tuple[int, int, float]] = []
    for entry in entries:
        first_entries.append(entry)
        break
    if shape is None:
        raise RuntimeError("matrix shape not found")

    def chained() -> Iterable[tuple[int, int, float]]:
        yield from first_entries
        yield from entries

    return shape, chained()


def module_score(sum_counts: dict[str, np.ndarray], total_counts: np.ndarray, genes: list[str]) -> np.ndarray:
    present = [gene for gene in genes if gene in sum_counts]
    if not present:
        return np.full(total_counts.shape, np.nan)
    scale = np.divide(10000.0, total_counts, out=np.zeros_like(total_counts, dtype=float), where=total_counts > 0)
    scores = []
    for gene in present:
        scores.append(np.log1p(sum_counts[gene] * scale))
    return np.vstack(scores).mean(axis=0)


def load_consensus_axis() -> tuple[list[str], list[str]]:
    sig = pd.read_csv(ROOT / "results/m4_geo_uc_ibd_signature/consensus_gene_signature.tsv", sep="\t")
    up = sig[sig["direction"].eq("up")].sort_values("consensus_score", ascending=False).head(150)["gene_symbol"].str.upper().tolist()
    down = sig[sig["direction"].eq("down")].sort_values("consensus_score", ascending=False).head(150)["gene_symbol"].str.upper().tolist()
    return up, down


def analyze_sample(outer: tarfile.TarFile, member: tarfile.TarInfo, axis_up: list[str], axis_down: list[str]) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    meta = parse_sample(member.name)
    raw = outer.extractfile(member)
    if raw is None:
        raise RuntimeError(f"cannot extract {member.name}")
    payload = gzip.decompress(raw.read())
    with tarfile.open(fileobj=io.BytesIO(payload), mode="r:") as sample_tar:
        gene_lines = read_text_member(
            sample_tar,
            ["/filtered/genes.tsv", "/filtered/features.tsv", "/filtered/genes.tsv.gz", "/filtered/features.tsv.gz"],
        )
        genes: list[str] = []
        for line in gene_lines:
            parts = line.split("\t")
            genes.append((parts[1] if len(parts) > 1 else parts[0]).upper())
        shape, entries = iter_matrix_entries(
            sample_tar,
            ["/filtered/matrix.mtx", "/filtered/matrix.mtx.txt", "/filtered/matrix.mtx.gz"],
        )
        n_genes, n_cells, _ = shape
        if n_genes != len(genes):
            raise RuntimeError(f"{member.name}: genes length mismatch")
        wanted = sorted(set(g for values in MODULES.values() for g in values) | set(axis_up) | set(axis_down))
        gene_to_idx = {gene: idx for idx, gene in enumerate(genes)}
        wanted_idx = {gene_to_idx[gene]: gene for gene in wanted if gene in gene_to_idx}
        total_counts = np.zeros(n_cells, dtype=float)
        sum_counts = {gene: np.zeros(n_cells, dtype=float) for gene in wanted_idx.values()}
        for gene_idx, cell_idx, value in entries:
            total_counts[cell_idx] += value
            gene = wanted_idx.get(gene_idx)
            if gene is not None:
                sum_counts[gene][cell_idx] += value

    cell_rows: list[dict[str, Any]] = []
    module_values: dict[str, np.ndarray] = {}
    for module, genes_for_module in MODULES.items():
        module_values[module] = module_score(sum_counts, total_counts, [g.upper() for g in genes_for_module])
    up_score = module_score(sum_counts, total_counts, axis_up)
    down_score = module_score(sum_counts, total_counts, axis_down)
    disease_axis = up_score - down_score
    module_values["consensus_uc_ibd_axis"] = disease_axis
    module_values["axis_up_score"] = up_score
    module_values["axis_down_score"] = down_score

    plasma = module_values["plasma_cell_marker"]
    naive = module_values["naive_b_marker"]
    plasma_skew = plasma - naive
    module_values["plasma_minus_naive_b_score"] = plasma_skew
    plasma_high = plasma_skew > 0

    summary: dict[str, Any] = {
        **meta,
        "n_filtered_cells": int(n_cells),
        "median_total_umi": fmt(float(np.median(total_counts))),
        "mean_total_umi": fmt(float(np.mean(total_counts))),
        "fraction_plasma_skewed_cells": fmt(float(np.mean(plasma_high))),
    }
    for module, values in module_values.items():
        finite = values[np.isfinite(values)]
        summary[f"mean_{module}"] = fmt(float(np.mean(finite))) if len(finite) else ""
        summary[f"median_{module}"] = fmt(float(np.median(finite))) if len(finite) else ""
    for module, genes_for_module in MODULES.items():
        present = [gene for gene in genes_for_module if gene.upper() in gene_to_idx]
        summary[f"{module}_present_genes"] = ";".join(present)
        summary[f"{module}_n_present_genes"] = len(present)
    summary["axis_up_present_genes"] = len([gene for gene in axis_up if gene in gene_to_idx])
    summary["axis_down_present_genes"] = len([gene for gene in axis_down if gene in gene_to_idx])

    # Store compact cell-level rows for the most relevant modules only.
    for idx in range(n_cells):
        cell_rows.append(
            {
                **meta,
                "cell_index": idx + 1,
                "total_umi": fmt(total_counts[idx]),
                "plasma_cell_marker": fmt(plasma[idx]),
                "naive_b_marker": fmt(naive[idx]),
                "plasma_minus_naive_b_score": fmt(plasma_skew[idx]),
                "bcl6_blnk_syk_axis": fmt(module_values["bcl6_blnk_syk_axis"][idx]),
                "consensus_uc_ibd_axis": fmt(disease_axis[idx]),
            }
        )
    return summary, cell_rows


def compare_groups(sample_rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    df = pd.DataFrame(sample_rows)
    metrics = [
        "fraction_plasma_skewed_cells",
        "mean_plasma_cell_marker",
        "mean_igg_plasma_marker",
        "mean_naive_b_marker",
        "mean_plasma_minus_naive_b_score",
        "mean_bcl6_blnk_syk_axis",
        "mean_curcumin_bcell_context_targets",
        "mean_consensus_uc_ibd_axis",
    ]
    out: list[dict[str, Any]] = []
    for metric in metrics:
        uc = pd.to_numeric(df.loc[df["condition"].eq("UC_inflamed"), metric], errors="coerce").dropna()
        hc = pd.to_numeric(df.loc[df["condition"].eq("healthy_control"), metric], errors="coerce").dropna()
        if len(uc) < 2 or len(hc) < 2:
            continue
        p_value = float(mannwhitneyu(uc, hc, alternative="two-sided").pvalue)
        out.append(
            {
                "metric": metric,
                "n_uc": len(uc),
                "n_healthy": len(hc),
                "mean_uc": fmt(uc.mean()),
                "mean_healthy": fmt(hc.mean()),
                "delta_uc_minus_healthy": fmt(uc.mean() - hc.mean()),
                "median_uc": fmt(uc.median()),
                "median_healthy": fmt(hc.median()),
                "mannwhitney_p": fmt(p_value),
                "direction": "higher_in_uc" if uc.mean() > hc.mean() else "lower_in_uc",
            }
        )
    return out


def main() -> int:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    axis_up, axis_down = load_consensus_axis()
    sample_rows: list[dict[str, Any]] = []
    cell_rows: list[dict[str, Any]] = []
    with tarfile.open(RAW_TAR, mode="r:") as outer:
        for member in sorted(outer.getmembers(), key=lambda m: m.name):
            if not member.name.endswith(".tar.gz"):
                continue
            summary, cells = analyze_sample(outer, member, axis_up, axis_down)
            sample_rows.append(summary)
            cell_rows.extend(cells)
            print(f"processed {member.name}: {summary['n_filtered_cells']} cells")

    sample_fields = sorted({key for row in sample_rows for key in row})
    ordered_sample_fields = [
        "gsm",
        "sample_id",
        "condition",
        "n_filtered_cells",
        "median_total_umi",
        "mean_total_umi",
        "fraction_plasma_skewed_cells",
        "mean_plasma_cell_marker",
        "mean_naive_b_marker",
        "mean_plasma_minus_naive_b_score",
        "mean_igg_plasma_marker",
        "mean_bcl6_blnk_syk_axis",
        "mean_curcumin_bcell_context_targets",
        "mean_consensus_uc_ibd_axis",
        "axis_up_present_genes",
        "axis_down_present_genes",
    ]
    ordered_sample_fields += [field for field in sample_fields if field not in ordered_sample_fields]
    write_tsv(OUT_DIR / "gse182270_bcell_sample_module_scores.tsv", sample_rows, ordered_sample_fields)
    write_tsv(
        OUT_DIR / "gse182270_bcell_group_comparison.tsv",
        compare_groups(sample_rows),
        [
            "metric",
            "n_uc",
            "n_healthy",
            "mean_uc",
            "mean_healthy",
            "delta_uc_minus_healthy",
            "median_uc",
            "median_healthy",
            "mannwhitney_p",
            "direction",
        ],
    )
    metadata_rows = []
    for row in sample_rows:
        metadata_rows.append(
            {
                "gsm": row["gsm"],
                "sample_id": row["sample_id"],
                "assigned_condition": row["condition"],
                "geo_status": row.get("geo_status", ""),
                "geo_disease": row.get("geo_disease", ""),
                "interpretation_note": row.get("sample_interpretation_note", ""),
            }
        )
    write_tsv(
        OUT_DIR / "gse182270_sample_metadata_interpretation.tsv",
        metadata_rows,
        ["gsm", "sample_id", "assigned_condition", "geo_status", "geo_disease", "interpretation_note"],
    )
    write_tsv(
        OUT_DIR / "gse182270_bcell_cell_module_scores_compact.tsv",
        cell_rows,
        [
            "gsm",
            "sample_id",
            "condition",
            "cell_index",
            "total_umi",
            "plasma_cell_marker",
            "naive_b_marker",
            "plasma_minus_naive_b_score",
            "bcl6_blnk_syk_axis",
            "consensus_uc_ibd_axis",
        ],
    )
    print(f"Wrote GSE182270 replication outputs to {OUT_DIR.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
