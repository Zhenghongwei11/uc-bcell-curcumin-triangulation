#!/usr/bin/env python3
"""Build Chinese Medicine-specific manuscript figures."""

from __future__ import annotations

from pathlib import Path
from textwrap import wrap

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "figures" / "output" / "chinese_medicine"
SOURCE = ROOT / "figures" / "source_data" / "chinese_medicine"

PALETTE = {
    "text": "#263238",
    "muted": "#69747C",
    "grid": "#D7DDE2",
    "healthy": "#718096",
    "disease": "#B54A42",
    "bcell": "#2F6B9A",
    "myeloid": "#7A6A3A",
    "gold": "#C78F2F",
    "green": "#4F7F62",
    "light": "#F6F8FA",
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
            "axes.edgecolor": PALETTE["text"],
            "axes.labelcolor": PALETTE["text"],
            "xtick.color": PALETTE["text"],
            "ytick.color": PALETTE["text"],
            "text.color": PALETTE["text"],
            "legend.frameon": False,
            "figure.facecolor": "white",
            "axes.facecolor": "white",
        }
    )


def read_tsv(path: str) -> pd.DataFrame:
    return pd.read_csv(ROOT / path, sep="\t")


def save_source(df: pd.DataFrame, name: str) -> None:
    SOURCE.mkdir(parents=True, exist_ok=True)
    df.to_csv(SOURCE / name, sep="\t", index=False)


def save_figure(fig: mpl.figure.Figure, stem: str) -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    for ext in ("svg", "pdf"):
        fig.savefig(OUT / f"{stem}.{ext}", bbox_inches="tight")
    fig.savefig(OUT / f"{stem}.png", dpi=300, bbox_inches="tight")
    fig.savefig(OUT / f"{stem}.tiff", dpi=600, bbox_inches="tight")
    plt.close(fig)


def panel_label(ax: mpl.axes.Axes, label: str) -> None:
    ax.text(-0.10, 1.06, label, transform=ax.transAxes, fontsize=10, fontweight="bold", va="top")


def p_text(p: float) -> str:
    if pd.isna(p):
        return "P = NA"
    if p < 0.001:
        return "P < 0.001"
    return f"P = {p:.3g}"


def q_text(q: float) -> str:
    if pd.isna(q):
        return "q = NA"
    if q < 0.001:
        return "q < 0.001"
    return f"q = {q:.3g}"


def bh_fdr(values: pd.Series) -> pd.Series:
    p = pd.to_numeric(values, errors="coerce")
    valid = p.notna()
    out = pd.Series(np.nan, index=values.index)
    if not valid.any():
        return out
    ordered = p[valid].sort_values()
    m = len(ordered)
    running = 1.0
    for rank, idx in reversed(list(enumerate(ordered.index, start=1))):
        running = min(running, float(p.loc[idx]) * m / rank)
        out.loc[idx] = min(running, 1.0)
    return out


def capped_neglog10(series: pd.Series, cap: float = 45.0) -> pd.Series:
    return (-np.log10(series.clip(lower=np.nextafter(0, 1)))).clip(upper=cap)


def bootstrap_delta_interval(values_case: np.ndarray, values_control: np.ndarray, seed: int = 20260710) -> tuple[float, float]:
    if len(values_case) == 0 or len(values_control) == 0:
        return (np.nan, np.nan)
    rng = np.random.default_rng(seed)
    deltas = []
    for _ in range(2000):
        c = rng.choice(values_case, size=len(values_case), replace=True)
        h = rng.choice(values_control, size=len(values_control), replace=True)
        deltas.append(float(np.mean(c) - np.mean(h)))
    return tuple(np.quantile(deltas, [0.025, 0.975]))


def display_variant_name(value: str) -> str:
    mapping = {
        "consensus_top50_logsum_cpm": "Consensus top 50, log-sum",
        "consensus_top50_avg_logcpm": "Consensus top 50, average",
        "consensus_top100_logsum_cpm": "Consensus top 100, log-sum",
        "consensus_top150_logsum_cpm": "Consensus top 150, log-sum",
        "consensus_top150_weighted_logsum_cpm": "Weighted top 150, log-sum",
        "consensus_top150_avg_logcpm": "Consensus top 150, average",
        "consensus_top150_weighted_avg_logcpm": "Weighted top 150, average",
        "consensus_top300_logsum_cpm": "Consensus top 300, log-sum",
        "consensus_top300_avg_logcpm": "Consensus top 300, average",
        "GSE75214_top150_logsum_cpm": "GSE75214 top 150, log-sum",
        "GSE87466_top150_logsum_cpm": "GSE87466 top 150, log-sum",
        "GSE75214_top150_avg_logcpm": "GSE75214 top 150, average",
        "GSE87466_top150_avg_logcpm": "GSE87466 top 150, average",
        "consensus_top100_avg_logcpm": "Consensus top 100, average",
    }
    return mapping.get(value, value.replace("_", " "))


def build_figure_1() -> None:
    group = read_tsv("data/derived/bulk_evidence_groups.tsv")
    scrna = read_tsv("data/derived/rectal_celltype_disease_contrasts.tsv")
    null = read_tsv("data/derived/rectal_patient_label_null.tsv")
    compound = read_tsv("tables/manuscript/table2_chinese_medicine_compound_mapping.tsv")
    gene = read_tsv("tables/manuscript/table3_curcumin_gene_followup.tsv")

    b = scrna[scrna["celltype"] == "B"].iloc[0]
    b_null = null[null["celltype"] == "B"].iloc[0]
    cur = compound[compound["Compound"] == "Curcumin"].iloc[0]

    flow = pd.DataFrame(
        [
            {
                "step": "Mucosal program",
                "result": f"{group['evidence_group'].nunique()} GEO evidence groups; active UC/IBD vs control",
            },
            {
                "step": "Cell-state localization",
                "result": f"B-lineage delta {b['delta_axis_diseased_minus_healthy']:.3f}; {p_text(b['p_axis_mannwhitney'])}; null {p_text(b_null['empirical_p_two_sided'])}",
            },
            {
                "step": "Compound mapping",
                "result": f"Deduplicated compound table; curcumin: {int(cur['Disease-context targets'])} disease-context targets",
            },
            {
                "step": "Gene follow-up",
                "result": f"{len(gene)} curcumin-associated genes with rectal B-lineage contrast; BH-FDR across 22 screened genes",
            },
        ]
    )
    save_source(flow, "figure1_study_design.tsv")

    fig, ax = plt.subplots(figsize=(7.4, 2.7))
    ax.set_axis_off()
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.text(0.01, 0.94, "Study design", fontsize=9, fontweight="bold", ha="left", va="top")
    ax.text(
        0.01,
        0.83,
        "Disease biology was defined before Chinese medicine-related compound mapping.",
        fontsize=6.5,
        color=PALETTE["muted"],
        ha="left",
        va="top",
    )

    colors = [PALETTE["disease"], PALETTE["bcell"], PALETTE["green"], PALETTE["gold"]]
    x_positions = [0.04, 0.29, 0.54, 0.79]
    width = 0.18
    for i, (_, row) in enumerate(flow.iterrows()):
        x = x_positions[i]
        ax.add_patch(
            mpl.patches.FancyBboxPatch(
                (x, 0.31),
                width,
                0.35,
                boxstyle="round,pad=0.012,rounding_size=0.012",
                linewidth=0.9,
                edgecolor=colors[i],
                facecolor="white",
            )
        )
        ax.add_patch(
            mpl.patches.Rectangle((x, 0.58), width, 0.08, linewidth=0, facecolor=colors[i])
        )
        ax.text(x + 0.012, 0.62, row["step"], fontsize=5.9, color="white", va="center", ha="left")
        ax.text(
            x + 0.012,
            0.53,
            "\n".join(wrap(row["result"], width=28, break_long_words=False)),
            fontsize=5.7,
            color=PALETTE["text"],
            va="top",
            ha="left",
            linespacing=1.15,
        )
        if i < len(x_positions) - 1:
            ax.annotate(
                "",
                xy=(x_positions[i + 1] - 0.025, 0.485),
                xytext=(x + width + 0.018, 0.485),
                arrowprops={"arrowstyle": "->", "lw": 0.9, "color": PALETTE["muted"]},
            )

    ax.text(
        0.04,
        0.16,
        "Primary interpretation: replicated UC/IBD mucosal program with rectal B-lineage localization.",
        fontsize=6.3,
        color=PALETTE["text"],
        ha="left",
    )
    ax.text(
        0.04,
        0.07,
        "Compound interpretation: annotation-based context for prioritizing future pharmacological experiments.",
        fontsize=6.0,
        color=PALETTE["muted"],
        ha="left",
    )
    save_figure(fig, "cm_figure1_study_design")


def build_figure_2() -> None:
    group = read_tsv("data/derived/bulk_evidence_groups.tsv")
    gse75214 = read_tsv("data/derived/gse75214_gene_signature.tsv")
    gse87466 = read_tsv("data/derived/gse87466_gene_signature.tsv")
    consensus = read_tsv("data/derived/consensus_gene_signature.tsv")

    volcano = pd.concat(
        [
            gse75214.assign(evidence_group_label="GPL6244 evidence group"),
            gse87466.assign(evidence_group_label="GPL13158 evidence group"),
        ],
        ignore_index=True,
    )
    volcano["neglog10_fdr"] = capped_neglog10(volcano["fdr"])
    volcano["direction_class"] = "Not selected"
    volcano.loc[(volcano["fdr"] < 0.05) & (volcano["log2fc"] >= 1), "direction_class"] = "Increased"
    volcano.loc[(volcano["fdr"] < 0.05) & (volcano["log2fc"] <= -1), "direction_class"] = "Decreased"
    save_source(volcano[["evidence_group_label", "gene_symbol", "log2fc", "fdr", "neglog10_fdr", "direction_class"]], "figure2_volcano_source.tsv")

    top = pd.concat(
        [
            consensus[consensus["direction"] == "up"].nlargest(10, "consensus_score"),
            consensus[consensus["direction"] == "down"].nlargest(10, "consensus_score"),
        ]
    ).copy()
    top = top.sort_values("mean_log2fc")
    save_source(top[["gene_symbol", "direction", "mean_log2fc", "best_fdr", "consensus_score"]], "figure2_consensus_genes.tsv")
    save_source(group, "figure2_geo_evidence_groups.tsv")

    fig = plt.figure(figsize=(7.1, 5.4), constrained_layout=True)
    axes = fig.subplot_mosaic([["A", "B", "C"], ["D", "D", "D"]], height_ratios=[1.0, 1.0], width_ratios=[0.95, 1.2, 1.2])

    ax = axes["A"]
    panel_label(ax, "A")
    group_plot = group.copy()
    y = np.arange(len(group_plot))
    left = np.zeros(len(group_plot))
    for col, label, color in [
        ("control_samples", "Control", PALETTE["healthy"]),
        ("case_samples", "Active disease", PALETTE["disease"]),
        ("excluded_samples", "Other/excluded", "#C9D0D6"),
    ]:
        ax.barh(y, group_plot[col], left=left, height=0.62, color=color, label=label)
        left += group_plot[col].to_numpy()
    ax.set_yticks(y)
    ax.set_yticklabels(group_plot["accession"])
    ax.invert_yaxis()
    ax.set_xlabel("Samples")
    ax.set_title("Bulk evidence groups", loc="left", fontsize=8)
    ax.legend(fontsize=5.8, loc="lower right")
    ax.grid(axis="x", color=PALETTE["grid"], lw=0.5)

    for key, label, title in [("B", "GPL6244 evidence group", "GPL6244 evidence group"), ("C", "GPL13158 evidence group", "GPL13158 evidence group")]:
        ax = axes[key]
        panel_label(ax, key)
        data = volcano[volcano["evidence_group_label"] == label]
        colors = data["direction_class"].map({"Increased": PALETTE["disease"], "Decreased": PALETTE["bcell"], "Not selected": "#C8D0D7"})
        ax.scatter(data["log2fc"], data["neglog10_fdr"], s=4, c=colors, alpha=0.42, linewidths=0)
        ax.axvline(0, color=PALETTE["text"], lw=0.6)
        ax.axvline(1, color=PALETTE["grid"], lw=0.6, ls="--")
        ax.axvline(-1, color=PALETTE["grid"], lw=0.6, ls="--")
        ax.axhline(-np.log10(0.05), color=PALETTE["grid"], lw=0.6, ls="--")
        ax.set_title(title, loc="left", fontsize=8)
        ax.set_xlabel("log2 fold change")
        ax.set_ylabel("-log10 FDR" if key == "B" else "")
        ax.grid(color=PALETTE["grid"], lw=0.4, alpha=0.5)

    ax = axes["D"]
    panel_label(ax, "D")
    colors = np.where(top["mean_log2fc"] >= 0, PALETTE["disease"], PALETTE["bcell"])
    ax.barh(np.arange(len(top)), top["mean_log2fc"], color=colors, height=0.68)
    ax.axvline(0, color=PALETTE["text"], lw=0.7)
    ax.set_yticks(np.arange(len(top)))
    ax.set_yticklabels(top["gene_symbol"], fontsize=6)
    ax.set_xlabel("Consensus mean log2 fold change")
    ax.set_title("Consensus mucosal disease-program genes", loc="left", fontsize=8)
    ax.grid(axis="x", color=PALETTE["grid"], lw=0.5)

    save_figure(fig, "cm_figure2_bulk_disease_program")


def build_figure_3() -> None:
    pseudo = read_tsv("data/derived/rectal_pseudobulk_module_scores.tsv")
    contrast = read_tsv("data/derived/rectal_celltype_disease_contrasts.tsv")
    null = read_tsv("data/derived/rectal_patient_label_null.tsv")

    keep = ["B", "M/DC", "T", "NK"]
    pseudo = pseudo[(pseudo["tissue_assignment"] == "R") & (pseudo["celltype"].isin(keep))].copy()
    contrast_all = contrast.copy()
    contrast_all["q_axis_bh_fdr"] = bh_fdr(contrast_all["p_axis_mannwhitney"])
    contrast = contrast_all[contrast_all["celltype"].isin(keep)].copy()
    null_all = null.copy()
    null_all["q_empirical_bh_fdr"] = bh_fdr(null_all["empirical_p_two_sided"])
    null = null_all[null_all["celltype"].isin(keep)].copy()

    order = ["B", "M/DC", "T", "NK"]
    save_source(pseudo, "figure3_patient_level_scores.tsv")
    save_source(contrast, "figure3_celltype_effects.tsv")
    save_source(null, "figure3_patient_label_null.tsv")

    fig = plt.figure(figsize=(7.0, 4.8), constrained_layout=True)
    axes = fig.subplot_mosaic([["A", "B"]], width_ratios=[1.4, 1.0])

    ax = axes["A"]
    panel_label(ax, "A")
    x_base = np.arange(len(order))
    rng = np.random.default_rng(20260710)
    for i, celltype in enumerate(order):
        sub = pseudo[pseudo["celltype"] == celltype]
        for status, color, offset in [
            ("healthy", PALETTE["healthy"], -0.16),
            ("diseased", PALETTE["disease"], 0.16),
        ]:
            vals = sub[sub["disease_assignment"] == status]["mean_disease_axis_score"].to_numpy()
            jitter = rng.normal(0, 0.025, size=len(vals))
            ax.scatter(
                np.full(len(vals), i + offset) + jitter,
                vals,
                s=24,
                color=color,
                alpha=0.88,
                edgecolor="white",
                linewidth=0.4,
                label=status.capitalize() if i == 0 else None,
                zorder=3,
            )
            if len(vals):
                ax.plot([i + offset - 0.08, i + offset + 0.08], [np.median(vals), np.median(vals)], color=color, lw=1.4)
        row = contrast[contrast["celltype"] == celltype].iloc[0]
        y = max(sub["mean_disease_axis_score"]) + 0.10
        ax.text(
            i,
            y,
            f"delta {row['delta_axis_diseased_minus_healthy']:.3f}\n{p_text(row['p_axis_mannwhitney'])}; {q_text(row['q_axis_bh_fdr'])}",
            ha="center",
            va="bottom",
            fontsize=5.6,
            color=PALETTE["muted"],
        )
    ax.set_xticks(x_base)
    ax.set_xticklabels(order)
    ax.set_ylabel("Disease-axis score")
    ax.set_title("Patient-level rectal pseudobulk scores", loc="left", fontsize=8, pad=6)
    ax.legend(loc="upper left", fontsize=6)
    ax.set_ylim(pseudo["mean_disease_axis_score"].min() - 0.15, pseudo["mean_disease_axis_score"].max() + 0.38)
    ax.grid(axis="y", color=PALETTE["grid"], lw=0.5, alpha=0.75)
    ax.tick_params(labelsize=6, length=2)

    ax = axes["B"]
    panel_label(ax, "B")
    null_plot = null.set_index("celltype").reindex(order).reset_index()
    y = np.arange(len(null_plot))
    colors = [PALETTE["bcell"], PALETTE["myeloid"], "#9AA5AE", "#9AA5AE"]
    ax.hlines(y, null_plot["null_q025_delta"], null_plot["null_q975_delta"], color="#C8D0D7", lw=2.4)
    ax.scatter(null_plot["observed_delta_axis_diseased_minus_healthy"], y, s=34, c=colors, zorder=3, edgecolor="white", linewidth=0.5)
    ax.axvline(0, color=PALETTE["text"], lw=0.7)
    ax.set_yticks(y)
    ax.set_yticklabels(null_plot["celltype"])
    ax.invert_yaxis()
    for i, row in null_plot.iterrows():
        ax.text(
            0.83,
            i,
            f"emp. {p_text(row['empirical_p_two_sided'])}\n{q_text(row['q_empirical_bh_fdr'])}",
            ha="right",
            va="center",
            fontsize=5.6,
            color=PALETTE["muted"],
        )
    ax.set_xlim(-0.55, 0.88)
    ax.set_xlabel("Observed delta and patient-label null interval")
    ax.set_title("Patient-label permutation", loc="left", fontsize=8, pad=6)
    ax.grid(axis="x", color=PALETTE["grid"], lw=0.5, alpha=0.75)
    ax.tick_params(labelsize=6, length=2)

    save_figure(fig, "cm_figure3_rectal_blineage_localization")


def build_figure_4() -> None:
    compounds = read_tsv("tables/manuscript/table2_chinese_medicine_compound_mapping.tsv")
    null_summary = read_tsv("data/derived/curcumin_matched_target_null.tsv").iloc[0]
    full = read_tsv("tables/supplementary/chinese_medicine_compound_mapping_full.tsv")

    plot = compounds.head(10).copy().sort_values("Rectal B-lineage increased targets")
    save_source(compounds, "figure4_compound_mapping.tsv")
    save_source(pd.DataFrame([null_summary]), "figure4_candidate_specificity_context.tsv")

    fig = plt.figure(figsize=(7.0, 5.3), constrained_layout=True)
    axes = fig.subplot_mosaic([["A", "B"], ["A", "C"]], width_ratios=[1.25, 1.0], height_ratios=[1.0, 1.0])

    ax = axes["A"]
    panel_label(ax, "A")
    y = np.arange(len(plot))
    sizes = 28 + plot["Disease-context targets"] * 5
    colors = [PALETTE["gold"] if x == "Curcumin" else PALETTE["bcell"] for x in plot["Compound"]]
    ax.scatter(plot["Rectal B-lineage increased targets"], y, s=sizes, c=colors, alpha=0.9, edgecolor="white", linewidth=0.5)
    for _, row in plot.iterrows():
        yv = list(plot["Compound"]).index(row["Compound"])
        ax.plot([row["M/DC increased targets"], row["Rectal B-lineage increased targets"]], [yv, yv], color="#B8C1C9", lw=1.0, zorder=0)
        ax.scatter(row["M/DC increased targets"], yv, s=22, color=PALETTE["myeloid"], alpha=0.85, edgecolor="white", linewidth=0.4)
    ax.set_yticks(y)
    ax.set_yticklabels(plot["Compound"], fontsize=6.2)
    ax.set_xlabel("Disease-increased targets")
    ax.set_title("Deduplicated compound mapping", loc="left", fontsize=8)
    ax.text(0.02, 0.02, "Gold/blue: B-lineage targets; brown: M/DC targets; point size: disease-context targets", transform=ax.transAxes, fontsize=5.6, color=PALETTE["muted"], ha="left", va="bottom")
    ax.grid(axis="x", color=PALETTE["grid"], lw=0.5)

    ax = axes["B"]
    panel_label(ax, "B")
    summary_cols = ["Disease-context targets", "Rectal B-lineage expressed targets", "Rectal B-lineage increased targets", "M/DC increased targets"]
    cur = compounds[compounds["Compound"] == "Curcumin"].iloc[0]
    vals = [cur[c] for c in summary_cols]
    labels = ["Disease\ncontext", "B-lineage\nexpressed", "B-lineage\nincreased", "M/DC\nincreased"]
    ax.bar(np.arange(len(vals)), vals, color=[PALETTE["green"], PALETTE["bcell"], PALETTE["gold"], PALETTE["myeloid"]], width=0.62)
    for i, v in enumerate(vals):
        ax.text(i, v + 0.5, str(int(v)), ha="center", fontsize=6)
    ax.set_xticks(np.arange(len(vals)))
    ax.set_xticklabels(labels, fontsize=6)
    ax.set_ylabel("Targets")
    ax.set_title("Curcumin descriptive target context", loc="left", fontsize=8)
    ax.grid(axis="y", color=PALETTE["grid"], lw=0.5)

    ax = axes["C"]
    panel_label(ax, "C")
    observed = float(null_summary["observed_target_set_score_without_direction_alignment"])
    null_mean = float(null_summary["null_mean"])
    null_q95 = float(null_summary["null_q95"])
    null_q99 = float(null_summary["null_q99"])
    ax.scatter([observed], [0], s=45, color=PALETTE["gold"], label="Curcumin observed", zorder=3)
    ax.hlines(0, null_mean, null_q99, color="#C8D0D7", lw=5, label="Null mean to 99th percentile")
    ax.scatter([null_mean, null_q95, null_q99], [0, 0, 0], s=[30, 30, 30], color=[PALETTE["healthy"], PALETTE["muted"], PALETTE["text"]], zorder=3)
    ax.set_yticks([])
    ax.set_xlabel("Target-set score")
    ax.set_title("Matched target-count null", loc="left", fontsize=8)
    ax.text(
        0.02,
        0.80,
        f"Empirical P = {float(null_summary['empirical_p_ge_observed']):.1f}\nscore-target count rho = {float(null_summary['target_count_score_spearman_rho_among_candidates']):.3f}",
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=6,
        color=PALETTE["muted"],
    )
    ax.grid(axis="x", color=PALETTE["grid"], lw=0.5)

    save_figure(fig, "cm_figure4_compound_mapping")


def build_figure_5() -> None:
    main = read_tsv("tables/manuscript/table3_curcumin_gene_followup.tsv")
    pseudo = read_tsv("data/derived/target_gene_pseudobulk_scores.tsv")
    genes = main["Gene"].tolist()
    rows = []
    for gene in genes:
        for celltype in ["B", "M/DC"]:
            sub = pseudo[(pseudo["tissue_assignment"] == "R") & (pseudo["celltype"] == celltype) & (pseudo["gene_symbol"] == gene)]
            diseased = sub[sub["disease_assignment"] == "diseased"]["log1p_cpm"].to_numpy()
            healthy = sub[sub["disease_assignment"] == "healthy"]["log1p_cpm"].to_numpy()
            if len(diseased) == 0 or len(healthy) == 0:
                continue
            delta = float(np.mean(diseased) - np.mean(healthy))
            lo, hi = bootstrap_delta_interval(diseased, healthy, seed=20260710 + len(rows))
            rows.append({"Gene": gene, "Cell type": celltype, "Delta": delta, "CI lower": lo, "CI upper": hi, "n diseased": len(diseased), "n healthy": len(healthy)})
    eff = pd.DataFrame(rows)
    eff = eff.merge(main[["Gene", "Rectal B-lineage FDR across 22 genes", "Reported direction"]], on="Gene", how="left")
    save_source(eff, "figure5_curcumin_gene_effects.tsv")

    order = main.sort_values("Rectal B-lineage delta", ascending=True)["Gene"].tolist()
    fig = plt.figure(figsize=(7.0, 5.2), constrained_layout=True)
    axes = fig.subplot_mosaic([["A", "B"]], width_ratios=[1.15, 1.0])

    ax = axes["A"]
    panel_label(ax, "A")
    b_eff = eff[eff["Cell type"] == "B"].set_index("Gene").reindex(order).reset_index()
    y = np.arange(len(b_eff))
    colors = [PALETTE["gold"] if g in {"BCL6", "BLNK", "SYK"} else PALETTE["bcell"] for g in b_eff["Gene"]]
    ax.hlines(y, b_eff["CI lower"], b_eff["CI upper"], color="#C8D0D7", lw=2.0)
    ax.scatter(b_eff["Delta"], y, c=colors, s=30, zorder=3, edgecolor="white", linewidth=0.5)
    ax.axvline(0, color=PALETTE["text"], lw=0.7)
    ax.set_yticks(y)
    ax.set_yticklabels(b_eff["Gene"], fontsize=6)
    ax.set_xlabel("Diseased minus healthy log1p CPM")
    ax.set_title("Rectal B-lineage expression effects", loc="left", fontsize=8)
    for i, row in b_eff.iterrows():
        ax.text(0.36, i, f"q = {row['Rectal B-lineage FDR across 22 genes']:.3g}", va="center", fontsize=5.5, color=PALETTE["muted"])
    ax.set_xlim(min(b_eff["CI lower"].min(), -0.08), max(b_eff["CI upper"].max(), 0.42))
    ax.grid(axis="x", color=PALETTE["grid"], lw=0.5)

    ax = axes["B"]
    panel_label(ax, "B")
    pivot = eff.pivot(index="Gene", columns="Cell type", values="Delta").reindex(order)
    ax.scatter(pivot["B"], pivot.index, color=PALETTE["bcell"], s=30, label="B-lineage", zorder=3)
    ax.scatter(pivot["M/DC"], pivot.index, color=PALETTE["myeloid"], s=30, label="M/DC", zorder=3)
    for gene in pivot.index:
        ax.plot([pivot.loc[gene, "M/DC"], pivot.loc[gene, "B"]], [gene, gene], color="#C8D0D7", lw=1.0, zorder=1)
    ax.axvline(0, color=PALETTE["text"], lw=0.7)
    ax.set_xlabel("Expression delta")
    ax.set_title("B-lineage versus M/DC context", loc="left", fontsize=8)
    ax.legend(fontsize=6, loc="lower right")
    ax.grid(axis="x", color=PALETTE["grid"], lw=0.5)

    save_figure(fig, "cm_figure5_curcumin_gene_followup")


def build_supplementary_figures() -> None:
    sens = read_tsv("data/derived/signature_sensitivity_blineage.tsv")
    marker = read_tsv("data/derived/bulk_blineage_marker_summary.tsv")
    proxy = read_tsv("data/derived/bulk_composition_proxy_summary.tsv")
    gse182 = read_tsv("data/derived/external_blineage_sample_scores.tsv")
    gse182_cmp = read_tsv("data/derived/external_blineage_group_comparison.tsv")
    ot = read_tsv("data/derived/open_targets_candidate_summary.tsv")
    etcm = read_tsv("data/derived/etcm2_mapped_overlap_summary.tsv")
    lincs = read_tsv("data/derived/lincs_candidate_reversal_status.tsv")
    pubmed = read_tsv("data/derived/pubmed_target_edge_verification.tsv")

    # S1
    s1 = sens.copy()
    s1["Variant"] = s1["variant"].map(display_variant_name)
    s1 = s1.sort_values("delta_axis_diseased_minus_healthy")
    save_source(s1, "supplementary_figure_s1_signature_sensitivity.tsv")
    fig, ax = plt.subplots(figsize=(6.2, 4.4), constrained_layout=True)
    y = np.arange(len(s1))
    colors = np.where(s1["q_empirical_bh"] < 0.05, PALETTE["bcell"], "#AEB8C2")
    colors = np.where(s1["delta_axis_diseased_minus_healthy"] < 0, PALETTE["disease"], colors)
    ax.hlines(y, s1["null_q025_delta"], s1["null_q975_delta"], color="#CBD3DA", lw=2.0)
    ax.scatter(s1["delta_axis_diseased_minus_healthy"], y, c=colors, s=28, zorder=3, edgecolor="white", linewidth=0.4)
    ax.axvline(0, color=PALETTE["text"], lw=0.7)
    ax.set_yticks(y)
    ax.set_yticklabels(s1["Variant"], fontsize=6)
    ax.set_xlabel("B-lineage disease-axis delta")
    ax.set_title("Figure S1. Disease-signature sensitivity", loc="left", fontsize=8)
    ax.grid(axis="x", color=PALETTE["grid"], lw=0.5)
    save_figure(fig, "cm_supplementary_figure_s1_signature_sensitivity")

    # S2
    save_source(marker, "supplementary_figure_s2_marker_summary.tsv")
    save_source(proxy, "supplementary_figure_s2_composition_proxy.tsv")
    fig = plt.figure(figsize=(7.0, 4.4), constrained_layout=True)
    axes = fig.subplot_mosaic([["A", "B"]], width_ratios=[1.0, 1.2])
    ax = axes["A"]
    panel_label(ax, "A")
    m = marker.copy()
    m["label"] = m["accession"] + "\n" + m["marker_group"].str.replace(" marker", "", regex=False)
    ax.bar(np.arange(len(m)), m["markers_up_fdr_lt_0_05"], color=np.where(m["marker_group"].str.contains("Plasma"), PALETTE["gold"], PALETTE["bcell"]))
    ax.set_xticks(np.arange(len(m)))
    ax.set_xticklabels(m["label"], rotation=45, ha="right", fontsize=5.7)
    ax.set_ylabel("Markers increased at FDR < 0.05")
    ax.set_title("B-lineage marker direction", loc="left", fontsize=8)
    ax.grid(axis="y", color=PALETTE["grid"], lw=0.5)
    ax = axes["B"]
    panel_label(ax, "B")
    p = proxy[proxy["celltype"].isin(["B", "M/DC"])].copy()
    p["label"] = p["accession"] + " " + p["celltype"]
    y = np.arange(len(p))
    ax.barh(y, p["delta_case_minus_control"], color=np.where(p["celltype"] == "B", PALETTE["bcell"], PALETTE["myeloid"]))
    ax.axvline(0, color=PALETTE["text"], lw=0.7)
    ax.set_yticks(y)
    ax.set_yticklabels(p["label"], fontsize=6)
    ax.set_xlabel("Case-control proxy fraction delta")
    ax.set_title("Marker-constrained composition proxy", loc="left", fontsize=8)
    ax.grid(axis="x", color=PALETTE["grid"], lw=0.5)
    save_figure(fig, "cm_supplementary_figure_s2_blineage_markers")

    # S3
    metrics = ["mean_plasma_cell_marker", "mean_igg_plasma_marker", "mean_bcl6_blnk_syk_axis", "mean_curcumin_bcell_context_targets", "mean_consensus_uc_ibd_axis"]
    labels = {
        "mean_plasma_cell_marker": "Plasma marker",
        "mean_igg_plasma_marker": "IgG plasma marker",
        "mean_bcl6_blnk_syk_axis": "BCL6/BLNK/SYK context",
        "mean_curcumin_bcell_context_targets": "Curcumin B-lineage targets",
        "mean_consensus_uc_ibd_axis": "Consensus disease axis",
    }
    save_source(gse182, "supplementary_figure_s3_gse182270_sample_scores.tsv")
    save_source(gse182_cmp, "supplementary_figure_s3_gse182270_group_comparison.tsv")
    fig, ax = plt.subplots(figsize=(6.5, 4.2), constrained_layout=True)
    rng = np.random.default_rng(20260710)
    for i, metric in enumerate(metrics):
        for condition, color, offset in [("healthy_control", PALETTE["healthy"], -0.13), ("UC_inflamed", PALETTE["disease"], 0.13)]:
            vals = gse182[gse182["condition"] == condition][metric].to_numpy()
            ax.scatter(np.full(len(vals), i + offset) + rng.normal(0, 0.02, len(vals)), vals, s=28, color=color, edgecolor="white", linewidth=0.4, label=condition.replace("_", " ") if i == 0 else None)
            ax.plot([i + offset - 0.06, i + offset + 0.06], [np.median(vals), np.median(vals)], color=color, lw=1.2)
        cmp_row = gse182_cmp[gse182_cmp["metric"] == metric].iloc[0]
        ax.text(i, max(gse182[metric]) + 0.06, f"P = {cmp_row['mannwhitney_p']:.3g}", ha="center", fontsize=5.6, color=PALETTE["muted"])
    ax.set_xticks(np.arange(len(metrics)))
    ax.set_xticklabels([labels[m] for m in metrics], rotation=25, ha="right", fontsize=6)
    ax.set_ylabel("Module score")
    ax.set_title("Figure S3. GSE182270 exploratory sample-level comparison", loc="left", fontsize=8)
    ax.legend(fontsize=6, loc="upper right")
    ax.set_ylim(-0.10, max(gse182[metrics].max()) + 0.28)
    ax.grid(axis="y", color=PALETTE["grid"], lw=0.5)
    save_figure(fig, "cm_supplementary_figure_s3_gse182270_exploratory")

    # S4
    cur_ot = ot[ot["ingredient_name"] == "Curcumin"].iloc[0]
    cur_etcm = etcm[etcm["ingredient_name"] == "Curcumin"].iloc[0]
    db = pd.DataFrame(
        [
            {"Resource": "Open Targets", "Metric": "Supported disease-context genes", "Count": cur_ot["n_open_targets_supported_genes"], "Denominator": cur_ot["n_disease_context_targets"]},
            {"Resource": "Open Targets", "Metric": "Genetic evidence genes", "Count": cur_ot["n_genetic_supported_genes"], "Denominator": cur_ot["n_disease_context_targets"]},
            {"Resource": "PubMed", "Metric": "Verified curcumin-target edges", "Count": (pubmed["title_match"] == "yes").sum(), "Denominator": len(pubmed)},
            {"Resource": "ETCM2", "Metric": "Overlap with HERB disease-context genes", "Count": cur_etcm["n_overlap_with_herb_m9_disease_context_targets"], "Denominator": max(cur_etcm["n_accepted_mapped_gene_symbols_exact_high"], 1)},
            {"Resource": "L1000CDS2", "Metric": "Top-result exact-ID hits", "Count": (lincs["hit_status"] != "no_top_result_hit").sum(), "Denominator": len(lincs)},
        ]
    )
    save_source(db, "supplementary_figure_s4_database_context_checks.tsv")
    fig, ax = plt.subplots(figsize=(6.4, 3.6), constrained_layout=True)
    db["Fraction"] = db["Count"] / db["Denominator"]
    y = np.arange(len(db))
    ax.barh(y, db["Fraction"], color=[PALETTE["green"], PALETTE["green"], PALETTE["gold"], PALETTE["myeloid"], PALETTE["healthy"]], height=0.62)
    for i, row in db.iterrows():
        ax.text(min(row["Fraction"] + 0.03, 0.98), i, f"{int(row['Count'])}/{int(row['Denominator'])}", va="center", fontsize=6)
    ax.set_yticks(y)
    ax.set_yticklabels(db["Resource"] + "\n" + db["Metric"], fontsize=6)
    ax.set_xlim(0, 1.08)
    ax.set_xlabel("Fraction")
    ax.set_title("Figure S4. External database context checks", loc="left", fontsize=8)
    ax.grid(axis="x", color=PALETTE["grid"], lw=0.5)
    save_figure(fig, "cm_supplementary_figure_s4_database_context_checks")


def build_graphical_abstract() -> None:
    fig, ax = plt.subplots(figsize=(7.38, 2.95))
    ax.set_axis_off()
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    steps = [
        ("UC/IBD mucosa", "replicated\ntranscriptomic program", PALETTE["disease"]),
        ("Rectal B-lineage", "patient-level\npseudobulk localization", PALETTE["bcell"]),
        ("TCM compounds", "deduplicated\nannotation mapping", PALETTE["green"]),
        ("Curcumin case", "exploratory\nB-cell-context genes", PALETTE["gold"]),
    ]
    x = [0.06, 0.31, 0.56, 0.81]
    for i, (title, body, color) in enumerate(steps):
        ax.add_patch(
            mpl.patches.FancyBboxPatch(
                (x[i], 0.34),
                0.15,
                0.34,
                boxstyle="round,pad=0.014,rounding_size=0.018",
                facecolor="white",
                edgecolor=color,
                linewidth=1.2,
            )
        )
        ax.add_patch(mpl.patches.Rectangle((x[i], 0.60), 0.15, 0.08, facecolor=color, edgecolor=color, linewidth=0))
        ax.text(x[i] + 0.075, 0.64, title, color="white", fontsize=6.0, ha="center", va="center")
        ax.text(x[i] + 0.075, 0.50, body, color=PALETTE["text"], fontsize=6.1, ha="center", va="center", linespacing=1.15)
        if i < 3:
            ax.annotate(
                "",
                xy=(x[i + 1] - 0.025, 0.51),
                xytext=(x[i] + 0.175, 0.51),
                arrowprops={"arrowstyle": "->", "lw": 1.0, "color": PALETTE["muted"]},
            )
    ax.text(
        0.5,
        0.20,
        "Disease-first public transcriptomics places Chinese medicine-related compound mapping into a rectal B-lineage context.",
        ha="center",
        va="center",
        fontsize=6.4,
        color=PALETTE["text"],
    )
    ax.text(
        0.5,
        0.10,
        "Curcumin-associated genes define experimentally testable B-cell-context hypotheses.",
        ha="center",
        va="center",
        fontsize=5.8,
        color=PALETTE["muted"],
    )
    save_figure(fig, "cm_graphical_abstract")


def main() -> None:
    setup_style()
    build_figure_1()
    build_figure_2()
    build_figure_3()
    build_figure_4()
    build_figure_5()
    build_supplementary_figures()
    build_graphical_abstract()


if __name__ == "__main__":
    main()
