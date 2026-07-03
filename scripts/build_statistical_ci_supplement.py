#!/usr/bin/env python3
"""Build bootstrap confidence-interval supplement from patient-level tables."""

from __future__ import annotations

import csv
from pathlib import Path
from typing import Any, Iterable

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
OUT_TABLE = ROOT / "tables" / "supplementary" / "statistical_ci_supplement.tsv"
OUT_FIG = ROOT / "figures" / "source_data" / "fig_statistical_ci_supplement.tsv"
SEED = 20260701
N_BOOT = 20000


def fmt(value: float) -> str:
    return f"{float(value):.6g}"


def bootstrap_delta(diseased: np.ndarray, healthy: np.ndarray, rng: np.random.Generator) -> tuple[float, float, float]:
    observed = float(np.mean(diseased) - np.mean(healthy))
    draws = np.empty(N_BOOT, dtype=float)
    for i in range(N_BOOT):
        d = rng.choice(diseased, size=len(diseased), replace=True)
        h = rng.choice(healthy, size=len(healthy), replace=True)
        draws[i] = np.mean(d) - np.mean(h)
    low, high = np.percentile(draws, [2.5, 97.5])
    return observed, float(low), float(high)


def row(
    analysis_family: str,
    feature: str,
    measure: str,
    diseased: Iterable[float],
    healthy: Iterable[float],
    rng: np.random.Generator,
    source_table: str,
    interpretation: str,
) -> dict[str, Any]:
    d = np.array(list(diseased), dtype=float)
    h = np.array(list(healthy), dtype=float)
    observed, low, high = bootstrap_delta(d, h, rng)
    return {
        "analysis_family": analysis_family,
        "feature": feature,
        "measure": measure,
        "n_diseased": len(d),
        "n_healthy": len(h),
        "observed_delta": fmt(observed),
        "bootstrap_ci_low": fmt(low),
        "bootstrap_ci_high": fmt(high),
        "bootstrap_iterations": N_BOOT,
        "seed": SEED,
        "source_table": source_table,
        "interpretation": interpretation,
    }


def write_tsv(path: Path, rows: list[dict[str, Any]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    rng = np.random.default_rng(SEED)
    rows: list[dict[str, Any]] = []

    m6_rel = "results/m6_scrna_cell_state/gse125527_pseudobulk_module_summary.tsv"
    m6 = pd.read_csv(ROOT / m6_rel, sep="\t")
    for celltype in ["B", "M/DC", "T", "NK", "unknown"]:
        sub = m6[(m6["tissue_assignment"] == "R") & (m6["celltype"] == celltype)]
        rows.append(
            row(
                "celltype_disease_axis",
                celltype,
                "mean_disease_axis_score",
                sub.loc[sub["disease_assignment"] == "diseased", "mean_disease_axis_score"],
                sub.loc[sub["disease_assignment"] == "healthy", "mean_disease_axis_score"],
                rng,
                m6_rel,
                "patient-level pseudobulk bootstrap CI for diseased-minus-healthy mean difference",
            )
        )

    m9_rel = "results/m9_target_cell_bridge/gse125527_target_gene_pseudobulk.tsv"
    m9 = pd.read_csv(ROOT / m9_rel, sep="\t")
    for gene in ["BCL6", "BLNK", "SYK", "IL1B", "TNF"]:
        sub = m9[(m9["tissue_assignment"] == "R") & (m9["celltype"] == "B") & (m9["gene_symbol"] == gene)]
        rows.append(
            row(
                "curcumin_target_bcell_expression",
                gene,
                "B-cell log1p CPM",
                sub.loc[sub["disease_assignment"] == "diseased", "log1p_cpm"],
                sub.loc[sub["disease_assignment"] == "healthy", "log1p_cpm"],
                rng,
                m9_rel,
                "patient-level pseudobulk bootstrap CI for diseased-minus-healthy mean difference; prioritization support only",
            )
        )

    fields = [
        "analysis_family",
        "feature",
        "measure",
        "n_diseased",
        "n_healthy",
        "observed_delta",
        "bootstrap_ci_low",
        "bootstrap_ci_high",
        "bootstrap_iterations",
        "seed",
        "source_table",
        "interpretation",
    ]
    write_tsv(OUT_TABLE, rows, fields)
    write_tsv(OUT_FIG, [{k: v for k, v in r.items() if k not in {"source_table", "interpretation"}} for r in rows], fields[:-2])
    print(f"Wrote {OUT_TABLE.relative_to(ROOT)} and {OUT_FIG.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
