#!/usr/bin/env python3
"""Build supplementary robustness and replication figures."""

from __future__ import annotations

from pathlib import Path
from textwrap import wrap

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "figures" / "output"

PALETTE = {
    "text": "#2B2B2B",
    "grid": "#D9D9D9",
    "blue": "#3F6FA6",
    "gold": "#C7892B",
    "red": "#B5473E",
    "olive": "#8A6F3D",
    "grey": "#7A869A",
    "light": "#F4F5F6",
    "very_light": "#FAFAFA",
}


def setup_style() -> None:
    mpl.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
            "svg.fonttype": "none",
            "pdf.fonttype": 42,
            "font.size": 7,
            "axes.spines.right": False,
            "axes.spines.top": False,
            "axes.linewidth": 0.8,
            "legend.frameon": False,
            "figure.facecolor": "white",
            "axes.facecolor": "white",
        }
    )


def read_tsv(path: str) -> pd.DataFrame:
    return pd.read_csv(ROOT / path, sep="\t")


def save(fig: mpl.figure.Figure, stem: str) -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    for ext in ["svg", "pdf"]:
        fig.savefig(OUT / f"{stem}.{ext}", bbox_inches="tight")
    fig.savefig(OUT / f"{stem}.png", dpi=240, bbox_inches="tight")
    fig.savefig(OUT / f"{stem}.tiff", dpi=600, bbox_inches="tight")
    plt.close(fig)


def panel_label(ax: mpl.axes.Axes, label: str) -> None:
    ax.text(-0.12, 1.06, label, transform=ax.transAxes, fontsize=9, fontweight="bold", va="top")


def clean_xgrid(ax: mpl.axes.Axes) -> None:
    ax.grid(axis="x", color=PALETTE["grid"], lw=0.5, alpha=0.75)
    ax.tick_params(labelsize=6, length=2)


def build_s1() -> None:
    sens = read_tsv("results/m17_no_download_robustness/scrna_signature_sensitivity_bcell_summary.tsv")
    proxy = read_tsv("results/m18_reviewer_lightweight_extensions/bulk_scrna_reference_deconvolution_proxy_summary.tsv")
    baseline = read_tsv("results/m18_reviewer_lightweight_extensions/standard_overlap_baseline_candidate_ranking.tsv")
    rep = read_tsv("results/m19_gse182270_bcell_replication/gse182270_bcell_group_comparison.tsv")

    fig = plt.figure(figsize=(7.2, 6.9), constrained_layout=True)
    axes = fig.subplot_mosaic(
        [["A", "B"], ["A", "C"], ["D", "D"]],
        width_ratios=[0.9, 1.55],
        height_ratios=[1.0, 1.0, 1.05],
    )

    ax = axes["A"]
    panel_label(ax, "a")
    ax.set_axis_off()
    ax.set_title("Supportive checks", loc="left", fontsize=8, pad=8)
    items = [
        ("Signature sensitivity", "B-cell localization is strongest for the replicated top-150 consensus axis; smaller or single-dataset variants are less stable.", PALETTE["blue"]),
        ("Composition proxy", "Bulk marker proxy does not show a simple replicated total B-cell fraction increase; M/DC proxy is more consistent.", PALETTE["olive"]),
        ("Baseline ranking", "A standard target-overlap baseline keeps curcumin near the top but is driven by generic inflammatory overlap.", PALETTE["gold"]),
        ("B-lineage replication", "GSE182270 shows directionally higher B-lineage module scores in UC, without sample-level significance.", PALETTE["red"]),
    ]
    y = 0.86
    for title, body, color in items:
        ax.add_patch(mpl.patches.Rectangle((0.02, y - 0.045), 0.025, 0.025, color=color, transform=ax.transAxes))
        ax.text(0.07, y, title, transform=ax.transAxes, ha="left", va="center", fontsize=6.4, fontweight="bold")
        ax.text(
            0.07,
            y - 0.085,
            "\n".join(wrap(body, width=34, break_long_words=False)),
            transform=ax.transAxes,
            ha="left",
            va="top",
            fontsize=5.6,
            color=PALETTE["grey"],
            linespacing=1.15,
        )
        y -= 0.23

    ax = axes["B"]
    panel_label(ax, "b")
    keep_variants = [
        "GSE75214_top150_logsum_cpm",
        "GSE87466_top150_logsum_cpm",
        "consensus_top50_logsum_cpm",
        "consensus_top100_logsum_cpm",
        "consensus_top150_logsum_cpm",
        "consensus_top150_weighted_logsum_cpm",
        "consensus_top300_logsum_cpm",
        "consensus_top150_avg_logcpm",
        "consensus_top150_weighted_avg_logcpm",
    ]
    label_map = {
        "GSE75214_top150_logsum_cpm": "GSE75214\nlog-sum",
        "GSE87466_top150_logsum_cpm": "GSE87466\nlog-sum",
        "consensus_top50_logsum_cpm": "Consensus 50\nlog-sum",
        "consensus_top100_logsum_cpm": "Consensus 100\nlog-sum",
        "consensus_top150_logsum_cpm": "Consensus 150\nlog-sum",
        "consensus_top150_weighted_logsum_cpm": "Weighted 150\nlog-sum",
        "consensus_top300_logsum_cpm": "Consensus 300\nlog-sum",
        "consensus_top150_avg_logcpm": "Consensus 150\naverage",
        "consensus_top150_weighted_avg_logcpm": "Weighted 150\naverage",
    }
    sens = sens[sens["variant"].isin(keep_variants)].copy()
    sens["variant_label"] = sens["variant"].map(label_map)
    sens_plot = sens.sort_values("delta_axis_diseased_minus_healthy")
    colors = np.where(sens_plot["q_empirical_bh"] < 0.05, PALETTE["blue"], "#B8BEC5")
    y = np.arange(len(sens_plot))
    ax.hlines(y, sens_plot["null_q025_delta"], sens_plot["null_q975_delta"], color="#CFD4D8", lw=2.2)
    ax.scatter(sens_plot["delta_axis_diseased_minus_healthy"], y, c=colors, s=26, zorder=3)
    ax.axvline(0, color=PALETTE["text"], lw=0.7)
    ax.set_yticks(y)
    ax.set_yticklabels(sens_plot["variant_label"], fontsize=5.4)
    ax.set_xlabel("B-cell disease-axis delta")
    ax.set_title("Signature sensitivity", loc="left", fontsize=8, pad=6)
    clean_xgrid(ax)

    ax = axes["C"]
    panel_label(ax, "c")
    focus = proxy[proxy["celltype"].isin(["B", "M/DC"])].copy()
    focus["label"] = focus["accession"] + " " + focus["celltype"]
    focus = focus.sort_values(["accession", "celltype"], ascending=[True, True])
    y = np.arange(len(focus))
    colors = [PALETTE["blue"] if ct == "B" else PALETTE["olive"] for ct in focus["celltype"]]
    ax.barh(y, focus["delta_case_minus_control"], color=colors, height=0.62)
    ax.axvline(0, color=PALETTE["text"], lw=0.7)
    for i, row in focus.reset_index(drop=True).iterrows():
        delta = float(row["delta_case_minus_control"])
        label_x = delta + 0.010 if delta >= 0 else 0.006
        ax.text(
            label_x,
            i,
            f"P={row['mannwhitney_p']:.3g}",
            va="center",
            ha="left",
            fontsize=5.5,
            color=PALETTE["grey"],
        )
    ax.set_yticks(y)
    ax.set_yticklabels(focus["label"], fontsize=5.8)
    ax.set_xlabel("Case-control proxy fraction delta")
    ax.set_title("Bulk immune-composition proxy", loc="left", fontsize=8, pad=6)
    ax.set_xlim(-0.025, 0.165)
    clean_xgrid(ax)

    ax = axes["D"]
    panel_label(ax, "d")
    rep_metrics = [
        "mean_plasma_cell_marker",
        "mean_igg_plasma_marker",
        "mean_bcl6_blnk_syk_axis",
        "mean_curcumin_bcell_context_targets",
        "mean_consensus_uc_ibd_axis",
    ]
    rep_plot = rep[rep["metric"].isin(rep_metrics)].copy()
    label_map = {
        "mean_plasma_cell_marker": "Plasma-cell\nmarker",
        "mean_igg_plasma_marker": "IgG plasma\nmarker",
        "mean_bcl6_blnk_syk_axis": "BCL6/BLNK/SYK\naxis",
        "mean_curcumin_bcell_context_targets": "Curcumin B-cell\ncontext targets",
        "mean_consensus_uc_ibd_axis": "Consensus UC/IBD\naxis",
    }
    rep_plot["label"] = rep_plot["metric"].map(label_map)
    x = np.arange(len(rep_plot))
    ax.bar(x, rep_plot["delta_uc_minus_healthy"], color=[PALETTE["red"], PALETTE["red"], PALETTE["gold"], PALETTE["gold"], PALETTE["blue"]], width=0.62)
    ax.axhline(0, color=PALETTE["text"], lw=0.7)
    for i, row in rep_plot.reset_index(drop=True).iterrows():
        y_text = max(float(row["delta_uc_minus_healthy"]) + 0.018, 0.035)
        ax.text(i, y_text, f"P={row['mannwhitney_p']:.3g}", ha="center", va="bottom", fontsize=5.5, color=PALETTE["grey"])
    ax.set_xticks(x)
    ax.set_xticklabels(rep_plot["label"], fontsize=5.8)
    ax.set_ylabel("UC minus control score delta")
    ax.set_title("GSE182270 focused B-lineage replication", loc="left", fontsize=8, pad=6)
    ax.text(0.99, 0.94, "5 UC inflamed vs 4 healthy/noninflamed samples", transform=ax.transAxes, ha="right", va="top", fontsize=5.5, color=PALETTE["grey"])
    ax.grid(axis="y", color=PALETTE["grid"], lw=0.5, alpha=0.75)
    ax.tick_params(labelsize=6, length=2)

    save(fig, "supplementary_figure_s1_robustness_replication")


def main() -> None:
    setup_style()
    build_s1()


if __name__ == "__main__":
    main()
