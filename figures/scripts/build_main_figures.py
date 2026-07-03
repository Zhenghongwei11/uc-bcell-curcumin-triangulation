#!/usr/bin/env python3
"""Build manuscript-style main figure drafts from tracked source results."""

from __future__ import annotations

import math
from pathlib import Path
from textwrap import wrap

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
SOURCE_DIR = ROOT / "figures" / "source_data"
OUTPUT_DIR = ROOT / "figures" / "output"

PALETTE = {
    "text": "#2B2B2B",
    "grid": "#D9D9D9",
    "healthy": "#7A869A",
    "disease": "#B5473E",
    "bcell": "#3F6FA6",
    "myeloid": "#8A6F3D",
    "curcumin": "#C7892B",
    "boundary": "#5F6B6D",
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


def read_tsv(relpath: str) -> pd.DataFrame:
    return pd.read_csv(ROOT / relpath, sep="\t")


def display_ingredient_name(name: str) -> str:
    return "Berberine" if name == "Berberime" else name


def make_unique_labels(labels: pd.Series) -> list[str]:
    seen: dict[str, int] = {}
    out: list[str] = []
    for label in labels.astype(str):
        seen[label] = seen.get(label, 0) + 1
        out.append(label if seen[label] == 1 else f"{label} (record {seen[label]})")
    return out


def add_bh_fdr(df: pd.DataFrame, p_col: str, out_col: str) -> pd.DataFrame:
    df = df.copy()
    pvals = pd.to_numeric(df[p_col], errors="coerce")
    valid = pvals.notna()
    adjusted = pd.Series(np.nan, index=df.index, dtype=float)
    if valid.any():
        ordered = pvals[valid].sort_values()
        m = len(ordered)
        running = 1.0
        for rank, idx in reversed(list(enumerate(ordered.index, start=1))):
            running = min(running, float(pvals.loc[idx]) * m / rank)
            adjusted.loc[idx] = min(running, 1.0)
    df[out_col] = adjusted
    return df


def save_source(df: pd.DataFrame, name: str) -> None:
    SOURCE_DIR.mkdir(parents=True, exist_ok=True)
    df.to_csv(SOURCE_DIR / name, sep="\t", index=False)


def save_figure(fig: mpl.figure.Figure, stem: str) -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    for ext in ("svg", "pdf"):
        fig.savefig(OUTPUT_DIR / f"{stem}.{ext}", bbox_inches="tight")
    fig.savefig(OUTPUT_DIR / f"{stem}.tiff", dpi=600, bbox_inches="tight")
    fig.savefig(OUTPUT_DIR / f"{stem}.png", dpi=240, bbox_inches="tight")
    plt.close(fig)


def panel_label(ax: mpl.axes.Axes, label: str) -> None:
    ax.text(
        -0.14,
        1.08,
        label,
        transform=ax.transAxes,
        fontsize=9,
        fontweight="bold",
        va="top",
        ha="left",
    )


def clean_axis(ax: mpl.axes.Axes) -> None:
    ax.grid(axis="x", color=PALETTE["grid"], linewidth=0.5, alpha=0.75)
    ax.tick_params(labelsize=6, length=2)


def p_label(p_value: float | int | None) -> str:
    if p_value is None or pd.isna(p_value):
        return "P=NA"
    if p_value < 0.001:
        return "P<0.001"
    return f"P={p_value:.3g}"


def q_label(q_value: float | int | None) -> str:
    if q_value is None or pd.isna(q_value):
        return "q=NA"
    if q_value < 0.001:
        return "q<0.001"
    return f"q={q_value:.3g}"


def capped_neglog10(series: pd.Series, cap: float = 50.0) -> pd.Series:
    values = -np.log10(series.clip(lower=np.nextafter(0, 1)))
    return values.clip(upper=cap)


def short_title(text: str, width: int = 28) -> str:
    return "\n".join(wrap(text, width=width, break_long_words=False))


def annotate_volcano_labels(ax: mpl.axes.Axes, rows: pd.DataFrame, offsets: list[tuple[float, float]]) -> None:
    for (_, row), (dx, dy) in zip(rows.iterrows(), offsets):
        ax.annotate(
            row["gene_symbol"],
            xy=(row["log2fc"], row["neglog10_fdr"]),
            xytext=(dx, dy),
            textcoords="offset points",
            fontsize=5.5,
            ha="center",
            va="bottom",
            arrowprops={"arrowstyle": "-", "linewidth": 0.35, "color": PALETTE["grid"]},
        )


def draw_gate_box(
    ax: mpl.axes.Axes,
    x: float,
    y: float,
    width: float,
    height: float,
    title: str,
    body: str,
    color: str,
    footer: str,
) -> None:
    rect = mpl.patches.Rectangle(
        (x, y),
        width,
        height,
        facecolor="white",
        edgecolor=color,
        linewidth=1.1,
    )
    ax.add_patch(rect)
    ax.add_patch(
        mpl.patches.Rectangle(
            (x, y + height - 0.18),
            width,
            0.18,
            facecolor=color,
            edgecolor=color,
            linewidth=0,
        )
    )
    ax.text(x + 0.04, y + height - 0.09, title, va="center", ha="left", fontsize=6.2, color="white")
    ax.text(x + 0.04, y + height - 0.25, body, va="top", ha="left", fontsize=5.25, color=PALETTE["text"], linespacing=1.15)
    ax.text(x + 0.04, y + 0.06, footer, va="bottom", ha="left", fontsize=5.2, color=PALETTE["boundary"])


def build_figure_1() -> None:
    group = read_tsv("results/m4_geo_uc_ibd_signature/group_audit.tsv")
    scrna = read_tsv("results/m6_scrna_cell_state/gse125527_disease_contrast_by_celltype.tsv")
    scrna_null = read_tsv("results/m6_scrna_cell_state/gse125527_disease_contrast_patient_label_null.tsv")
    m7 = read_tsv("results/m7_tcm_scrna_bridge/tcm_scrna_axis_bridge.tsv")
    m9 = read_tsv("results/m9_target_cell_bridge/candidate_target_cell_bridge_summary.tsv")
    ot_summary = read_tsv("results/m10_open_targets/candidate_open_targets_summary.tsv")
    m12 = read_tsv("results/m12_pubmed_verification/curcumin_target_edge_pubmed_verification.tsv")
    m16 = read_tsv("results/m16_curcumin_fulltext_evidence/curcumin_high_value_target_fulltext_audit.tsv")
    etcm2 = read_tsv("results/m15_etcm2_target_mapping/candidate_etcm2_mapped_overlap_summary.tsv")
    lincs = read_tsv("results/m5_lincs_reversal/candidate_reversal_status.tsv")

    b_row = scrna[scrna["celltype"] == "B"].iloc[0]
    b_null = scrna_null[scrna_null["celltype"] == "B"].iloc[0]
    cur_m9 = m9[m9["ingredient_name"] == "Curcumin"].iloc[0]
    cur_ot = ot_summary[ot_summary["ingredient_name"] == "Curcumin"].iloc[0]
    cur_etcm = etcm2[etcm2["ingredient_name"] == "Curcumin"].iloc[0]
    lincs_hits = int((lincs["hit_status"] != "no_top_result_hit").sum())

    evidence_gates = pd.DataFrame(
        [
            {
                "layer": "Disease signature",
                "source_layer": "Bulk GEO",
                "key_result": f"{group['evidence_group'].nunique()} independent evidence groups",
                "allowed_claim": "UC/IBD mucosal disease signature",
                "boundary": "Not a treatment mechanism",
            },
            {
                "layer": "Cell localization",
                "source_layer": "GSE125527 scRNA pseudobulk",
                "key_result": f"B-cell delta={b_row['delta_axis_diseased_minus_healthy']:.3f}; null P={b_null['empirical_p_two_sided']:.4g}",
                "allowed_claim": "Rectal B-cell disease-axis localization",
                "boundary": "Association, not compound action",
            },
            {
                "layer": "Candidate linkage",
                "source_layer": "HERB disease-context targets",
                "key_result": f"Curcumin linkage={cur_m9['bridge_score']:.1f}; targets={int(cur_m9['n_disease_context_unique_targets'])}",
                "allowed_claim": "Curcumin prioritized as B-cell-axis candidate",
                "boundary": "Prioritization, not efficacy",
            },
            {
                "layer": "Orthogonal support",
                "source_layer": "Open Targets and PubMed/full text",
                "key_result": f"Open Targets {int(cur_ot['n_open_targets_supported_genes'])}/{int(cur_ot['n_disease_context_targets'])}; PubMed title matches {int((m12['title_match'] == 'yes').sum())}/{len(m12)}",
                "allowed_claim": "Target disease plausibility and evidence traceability",
                "boundary": "Does not prove compound-target engagement",
            },
            {
                "layer": "Resource limits",
                "source_layer": "ETCM2 and L1000CDS2",
                "key_result": f"ETCM2 overlap={int(cur_etcm['n_overlap_with_herb_m9_disease_context_targets'])}; LINCS hits={lincs_hits}/{len(lincs)}",
                "allowed_claim": "Weak or negative layers reported",
                "boundary": "Non-hits do not prove inactivity",
            },
        ]
    )
    save_source(evidence_gates, "fig1a_evidence_gate_summary.tsv")

    role_map = {
        "primary_lead_for_b_cell_axis": "Lead B-cell-axis",
        "secondary_lead_for_myeloid_axis": "Secondary myeloid-axis",
        "backup_tcm_candidate": "Lower consensus",
        "low_priority": "Low priority",
        "exclude_or_background": "Background/excluded",
    }
    role_summary = (
        m7["recommended_role"]
        .map(role_map)
        .value_counts()
        .rename_axis("candidate_role")
        .reset_index(name="n_candidates")
    )
    role_summary["display_order"] = role_summary["candidate_role"].map(
        {
            "Lead B-cell-axis": 0,
            "Secondary myeloid-axis": 1,
            "Lower consensus": 2,
            "Low priority": 3,
            "Background/excluded": 4,
        }
    )
    role_summary = role_summary.sort_values("display_order")
    save_source(role_summary, "fig1b_candidate_role_summary.tsv")

    claim_boundaries = pd.DataFrame(
        [
            {
                "evidence_class": "Disease association",
                "can_claim": "UC/IBD mucosal disease signature",
                "cannot_claim": "Treatment response or causal disease driver",
            },
            {
                "evidence_class": "Cell-state localization",
                "can_claim": "Rectal B-cell enrichment of disease-axis signal",
                "cannot_claim": "Curcumin acts through B cells in patients",
            },
            {
                "evidence_class": "Target bridge",
                "can_claim": "Curcumin prioritized by disease-context targets",
                "cannot_claim": "Validated direct compound-target engagement",
            },
            {
                "evidence_class": "Full-text support",
                "can_claim": "BCL6/BLNK/SYK as main B-cell mechanism candidates",
                "cannot_claim": "CCL2/IL33/IL1B/TNF as B-cell-intrinsic targets",
            },
            {
                "evidence_class": "Boundary resources",
                "can_claim": "ETCM2/LINCS constrain interpretation",
                "cannot_claim": "ETCM2 or LINCS validates the mechanism",
            },
        ]
    )
    save_source(claim_boundaries, "fig1c_claim_boundary_summary.tsv")

    fig = plt.figure(figsize=(7.4, 6.8), constrained_layout=True)
    axes = fig.subplot_mosaic([["A", "A"], ["B", "C"]], height_ratios=[1.22, 1.2], width_ratios=[0.9, 1.55])

    ax = axes["A"]
    panel_label(ax, "a")
    ax.set_axis_off()
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_title("Disease-first evidence synthesis", loc="left", fontsize=8, pad=8)
    gate_colors = [PALETTE["disease"], PALETTE["bcell"], PALETTE["curcumin"], PALETTE["myeloid"], PALETTE["boundary"]]
    gate_bodies = [
        "2 GEO groups\nUC/IBD signature",
        f"B-cell delta {b_row['delta_axis_diseased_minus_healthy']:.3f}\nnull P={b_null['empirical_p_two_sided']:.4g}",
        f"Curcumin linkage {cur_m9['bridge_score']:.1f}\n22 disease targets",
        f"Open Targets {int(cur_ot['n_open_targets_supported_genes'])}/{int(cur_ot['n_disease_context_targets'])}\nPubMed {int((m12['title_match'] == 'yes').sum())}/{len(m12)}",
        f"ETCM2 overlap {int(cur_etcm['n_overlap_with_herb_m9_disease_context_targets'])}\nLINCS {lincs_hits}/{len(lincs)}",
    ]
    footers = ["association", "localization", "prioritization", "plausibility", "limits"]
    titles = ["Disease", "Cell state", "TCM candidate", "Support", "Limits"]
    x0 = 0.02
    w = 0.17
    gap = 0.025
    for i, (title, body, color, footer) in enumerate(zip(titles, gate_bodies, gate_colors, footers)):
        x = x0 + i * (w + gap)
        draw_gate_box(ax, x, 0.36, w, 0.43, title, body, color, footer)
        if i < 4:
            ax.annotate(
                "",
                xy=(x + w + gap * 0.72, 0.58),
                xytext=(x + w + gap * 0.15, 0.58),
                arrowprops={"arrowstyle": "->", "lw": 0.8, "color": PALETTE["boundary"]},
            )
    ax.text(
        0.02,
        0.19,
        "Central hypothesis: public data triangulation nominates a Curcumin-linked rectal B-cell inflammatory axis in UC/IBD.",
        fontsize=6.5,
        ha="left",
        va="center",
        color=PALETTE["text"],
    )
    ax.text(
        0.02,
        0.08,
        "Interpretation is tiered: association -> localization -> prioritization -> plausibility -> evidentiary limits.",
        fontsize=5.6,
        ha="left",
        va="center",
        color=PALETTE["boundary"],
    )

    ax = axes["B"]
    panel_label(ax, "b")
    ax.set_title("Candidate convergence", loc="left", fontsize=8, pad=16)
    plot_roles = role_summary.sort_values("display_order", ascending=False)
    role_colors = {
        "Lead B-cell-axis": PALETTE["curcumin"],
        "Secondary myeloid-axis": PALETTE["myeloid"],
        "Lower consensus": "#AEB6BE",
        "Low priority": "#D7DBDF",
        "Background/excluded": "#ECEFF1",
    }
    ax.barh(
        plot_roles["candidate_role"],
        plot_roles["n_candidates"],
        color=[role_colors[r] for r in plot_roles["candidate_role"]],
        height=0.62,
    )
    for _, row in plot_roles.iterrows():
        ax.text(row["n_candidates"] + 1.0, row["candidate_role"], str(int(row["n_candidates"])), va="center", fontsize=5.8, color=PALETTE["boundary"])
    ax.set_xlabel("Candidates")
    ax.set_xlim(0, max(plot_roles["n_candidates"]) * 1.22)
    clean_axis(ax)
    ax.text(
        0.02,
        -0.18,
        f"Curcumin: 13 HERB evidence records; {int(cur_m9['n_b_cell_expressed_targets'])} B-cell-expressed targets; {int(cur_m9['n_b_cell_disease_increased_targets'])} increased in diseased B-cell pseudobulk",
        transform=ax.transAxes,
        fontsize=5.4,
        ha="left",
        va="bottom",
        color=PALETTE["boundary"],
    )

    ax = axes["C"]
    panel_label(ax, "c")
    ax.set_title("Claim map", loc="left", fontsize=8, pad=16)
    ax.set_axis_off()
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    compact_boundaries = pd.DataFrame(
        [
            ("Disease", "UC/IBD mucosal\nsignature", "treatment\nresponse"),
            ("Cell state", "rectal B-cell\nlocalization", "patient compound\naction"),
            ("Target linkage", "Curcumin\nprioritization", "validated target\nengagement"),
            ("Full text", "BCL6/BLNK/SYK\nmechanism", "cytokines as\nB-cell-intrinsic"),
            ("Resource limits", "ETCM2/LINCS\nconstraints", "validation by\nnon-hit layers"),
        ],
        columns=["evidence_class", "can_claim", "cannot_claim"],
    )
    row_y = np.linspace(0.80, 0.12, len(compact_boundaries))
    for y, (_, row) in zip(row_y, compact_boundaries.iterrows()):
        ax.text(0.02, y, row["evidence_class"], fontsize=6.0, fontweight="bold", ha="left", va="center", color=PALETTE["text"])
        ax.text(0.33, y, row["can_claim"], fontsize=5.35, ha="left", va="center", color=PALETTE["bcell"], linespacing=1.12)
        ax.text(0.75, y, row["cannot_claim"], fontsize=5.35, ha="left", va="center", color=PALETTE["boundary"], linespacing=1.12)
        ax.plot([0.02, 0.98], [y - 0.083, y - 0.083], color=PALETTE["grid"], lw=0.5)
    ax.text(0.33, 0.93, "Supported wording", fontsize=5.8, color=PALETTE["bcell"], ha="left", va="center")
    ax.text(0.75, 0.93, "Not supported", fontsize=5.8, color=PALETTE["boundary"], ha="left", va="center")

    save_figure(fig, "fig1_study_design_evidence_gating")


def build_figure_2() -> None:
    group = read_tsv("results/m4_geo_uc_ibd_signature/group_audit.tsv")
    gse75214 = read_tsv("results/m4_geo_uc_ibd_signature/GSE75214_gene_signature.tsv")
    gse87466 = read_tsv("results/m4_geo_uc_ibd_signature/GSE87466_gene_signature.tsv")
    consensus = read_tsv("results/m4_geo_uc_ibd_signature/consensus_gene_signature.tsv")

    group_source = group[
        [
            "accession",
            "role",
            "evidence_group",
            "case_samples",
            "control_samples",
            "excluded_samples",
            "gene_rows",
        ]
    ].copy()
    save_source(group_source, "fig2a_geo_group_audit.tsv")

    volcano = pd.concat(
        [
            gse75214.assign(display_group="GSE75214 discovery"),
            gse87466.assign(display_group="GSE87466 validation"),
        ],
        ignore_index=True,
    )
    volcano = volcano[
        ["display_group", "accession", "gene_symbol", "log2fc", "fdr", "p_value", "t_stat"]
    ].copy()
    volcano["neglog10_fdr"] = capped_neglog10(volcano["fdr"])
    volcano["plot_class"] = "not_significant"
    volcano.loc[(volcano["fdr"] < 0.05) & (volcano["log2fc"] >= 1), "plot_class"] = "up"
    volcano.loc[(volcano["fdr"] < 0.05) & (volcano["log2fc"] <= -1), "plot_class"] = "down"
    save_source(volcano, "fig2b_volcano_source.tsv")

    consensus_source = consensus[
        [
            "gene_symbol",
            "direction",
            "evidence_group_count",
            "direction_consistency",
            "mean_log2fc",
            "best_fdr",
            "consensus_score",
        ]
    ].copy()
    save_source(consensus_source, "fig2c_consensus_signature.tsv")

    fig = plt.figure(figsize=(7.2, 6.2), constrained_layout=True)
    mosaic = [["A", "B1", "B2"], ["C", "C", "C"]]
    axes = fig.subplot_mosaic(mosaic, width_ratios=[1.05, 1.25, 1.25], height_ratios=[1.05, 1.0])

    ax = axes["A"]
    panel_label(ax, "a")
    y = np.arange(len(group_source))
    left = np.zeros(len(group_source))
    bars = [
        ("control_samples", "Control", PALETTE["healthy"]),
        ("case_samples", "Active UC", PALETTE["disease"]),
        ("excluded_samples", "Excluded", "#C9CDD2"),
    ]
    for col, label, color in bars:
        ax.barh(y, group_source[col], left=left, color=color, height=0.62, label=label)
        left += group_source[col].to_numpy()
    ax.set_yticks(y)
    ax.set_yticklabels(group_source["accession"])
    ax.invert_yaxis()
    ax.set_xlabel("Samples")
    ax.set_title("GEO evidence groups", loc="left", fontsize=8, pad=6)
    for i, row in group_source.iterrows():
        label = "non-independent\nsensitivity" if "non_independent" in row["role"] else "independent\nvalidation" if "validation" in row["role"] else "discovery"
        ax.text(left[i] + 5, i, label, va="center", fontsize=5.8, color=PALETTE["boundary"])
    ax.legend(loc="lower right", fontsize=5.8, handlelength=1.0)
    ax.set_xlim(0, max(left) * 1.6)
    clean_axis(ax)

    for axis_key, display_group, ax_title in [
        ("B1", "GSE75214 discovery", "Discovery"),
        ("B2", "GSE87466 validation", "Independent validation"),
    ]:
        ax = axes[axis_key]
        if axis_key == "B1":
            panel_label(ax, "b")
        data = volcano[volcano["display_group"] == display_group]
        colors = data["plot_class"].map(
            {"up": PALETTE["disease"], "down": PALETTE["bcell"], "not_significant": "#C8CED4"}
        )
        ax.scatter(data["log2fc"], data["neglog10_fdr"], s=3, c=colors, alpha=0.45, linewidths=0)
        ax.axvline(0, color=PALETTE["text"], linewidth=0.5)
        ax.axvline(1, color=PALETTE["grid"], linewidth=0.5, linestyle="--")
        ax.axvline(-1, color=PALETTE["grid"], linewidth=0.5, linestyle="--")
        ax.axhline(-math.log10(0.05), color=PALETTE["grid"], linewidth=0.5, linestyle="--")
        ax.set_title(ax_title, loc="left", fontsize=8, pad=6)
        ax.set_xlabel("log2 fold change")
        ax.set_ylabel("-log10 FDR" if axis_key == "B1" else "")
        if axis_key == "B1":
            label_genes = ["SLC26A2", "ABCA12", "PADI2", "TCN1"]
            offsets = [(0, 7), (0, 7), (-10, 7), (10, 7)]
        else:
            label_genes = ["RTEL1", "SLC26A2", "DHRS11", "MTMR11"]
            offsets = [(0, 7), (-18, 9), (16, 10), (8, 18)]
        labels = data[data["gene_symbol"].isin(label_genes)].drop_duplicates("gene_symbol")
        labels = labels.set_index("gene_symbol").reindex(label_genes).dropna(subset=["log2fc"]).reset_index()
        annotate_volcano_labels(ax, labels, offsets)
        clean_axis(ax)

    ax = axes["C"]
    panel_label(ax, "c")
    top_up = consensus_source[consensus_source["direction"] == "up"].nlargest(10, "consensus_score")
    top_down = consensus_source[consensus_source["direction"] == "down"].nlargest(10, "consensus_score")
    plot_df = pd.concat([top_down.sort_values("mean_log2fc"), top_up.sort_values("mean_log2fc")])
    colors = np.where(plot_df["mean_log2fc"] >= 0, PALETTE["disease"], PALETTE["bcell"])
    ax.barh(np.arange(len(plot_df)), plot_df["mean_log2fc"], color=colors, height=0.72)
    ax.axvline(0, color=PALETTE["text"], linewidth=0.7)
    ax.set_yticks(np.arange(len(plot_df)))
    ax.set_yticklabels(plot_df["gene_symbol"], fontsize=6)
    ax.set_xlabel("Consensus mean log2 fold change")
    ax.set_title("Consensus mucosal signature genes", loc="left", fontsize=8, pad=6)
    ax.text(
        0.99,
        0.04,
        "Two evidence groups; GSE59071 retained only as non-independent sensitivity",
        transform=ax.transAxes,
        ha="right",
        va="bottom",
        fontsize=6,
        color=PALETTE["boundary"],
    )
    clean_axis(ax)

    save_figure(fig, "fig2_geo_disease_signature")


def build_figure_3() -> None:
    contrast = read_tsv("results/m6_scrna_cell_state/gse125527_disease_contrast_by_celltype.tsv")
    null = read_tsv("results/m6_scrna_cell_state/gse125527_disease_contrast_patient_label_null.tsv")
    pseudo = read_tsv("results/m6_scrna_cell_state/gse125527_pseudobulk_module_summary.tsv")

    contrast = add_bh_fdr(contrast, "p_axis_mannwhitney", "q_axis_bh_fdr")
    null = add_bh_fdr(null, "empirical_p_two_sided", "q_empirical_bh_fdr")
    unclassified_qc = contrast[contrast["celltype"].eq("unknown")].copy()
    contrast = contrast[~contrast["celltype"].eq("unknown")].sort_values("delta_axis_diseased_minus_healthy", ascending=True)
    null = null[~null["celltype"].eq("unknown")].copy()
    save_source(contrast, "fig3a_celltype_disease_axis.tsv")
    save_source(null, "fig3b_patient_label_null.tsv")
    save_source(unclassified_qc, "fig3_unclassified_cell_qc.tsv")

    counts = (
        pseudo.groupby(["celltype", "disease_assignment"], as_index=False)
        .agg(n_patients=("patient_assignment", "nunique"), n_pseudobulk_samples=("patient_assignment", "size"), n_cells=("n_cells", "sum"))
    )
    counts_plot = counts[~counts["celltype"].eq("unknown")].copy()
    save_source(counts, "fig3c_pseudobulk_counts.tsv")

    fig = plt.figure(figsize=(7.1, 5.7), constrained_layout=True)
    axes = fig.subplot_mosaic([["A", "B"], ["C", "C"]], height_ratios=[1.2, 0.9])

    ax = axes["A"]
    panel_label(ax, "a")
    colors = [
        PALETTE["bcell"] if c == "B" else PALETTE["myeloid"] if c == "M/DC" else "#AEB6BE"
        for c in contrast["celltype"]
    ]
    ax.barh(
        contrast["celltype"],
        contrast["delta_axis_diseased_minus_healthy"],
        color=colors,
        height=0.65,
    )
    ax.axvline(0, color=PALETTE["text"], linewidth=0.7)
    for _, row in contrast.iterrows():
        ax.text(
            row["delta_axis_diseased_minus_healthy"] + 0.018,
            row["celltype"],
            f"{p_label(row['p_axis_mannwhitney'])}\n{q_label(row['q_axis_bh_fdr'])}",
            va="center",
            fontsize=5.8,
            color=PALETTE["boundary"],
        )
    ax.set_xlabel("Disease-axis delta, diseased minus healthy")
    ax.set_title("Rectal cell-type localization", loc="left", fontsize=8, pad=6)
    ax.set_xlim(-0.12, 0.68)
    clean_axis(ax)

    ax = axes["B"]
    panel_label(ax, "b")
    null_plot = null.sort_values("observed_delta_axis_diseased_minus_healthy", ascending=True)
    y = np.arange(len(null_plot))
    colors = [
        PALETTE["bcell"] if c == "B" else PALETTE["myeloid"] if c == "M/DC" else "#AEB6BE"
        for c in null_plot["celltype"]
    ]
    ax.hlines(y, null_plot["null_q025_delta"], null_plot["null_q975_delta"], color="#AEB6BE", linewidth=2.2)
    ax.scatter(null_plot["observed_delta_axis_diseased_minus_healthy"], y, c=colors, s=25, zorder=3)
    ax.axvline(0, color=PALETTE["grid"], linewidth=0.7)
    ax.set_yticks(y)
    ax.set_yticklabels(null_plot["celltype"])
    for i, row in null_plot.reset_index(drop=True).iterrows():
        ax.text(
            0.73,
            i,
            f"emp. {q_label(row['q_empirical_bh_fdr'])}",
            va="center",
            ha="right",
            fontsize=5.8,
            color=PALETTE["boundary"],
        )
    ax.set_xlabel("Observed delta vs patient-label null 95% interval")
    ax.set_title("Patient-label null support", loc="left", fontsize=8, pad=6)
    ax.set_xlim(-0.55, 0.82)
    clean_axis(ax)

    ax = axes["C"]
    panel_label(ax, "c")
    pivot_cells = counts_plot.pivot(index="celltype", columns="disease_assignment", values="n_cells").fillna(0)
    pivot_patients = counts_plot.pivot(index="celltype", columns="disease_assignment", values="n_patients").fillna(0)
    order = contrast.sort_values("delta_axis_diseased_minus_healthy", ascending=False)["celltype"].tolist()
    x = np.arange(len(order))
    width = 0.36
    healthy = pivot_cells.reindex(order).get("healthy", pd.Series(0, index=order))
    diseased = pivot_cells.reindex(order).get("diseased", pd.Series(0, index=order))
    ax.bar(x - width / 2, healthy, width=width, color=PALETTE["healthy"], label="Healthy cells")
    ax.bar(x + width / 2, diseased, width=width, color=PALETTE["disease"], label="Diseased cells")
    for i, celltype in enumerate(order):
        h_pat = int(pivot_patients.reindex(order).loc[celltype].get("healthy", 0))
        d_pat = int(pivot_patients.reindex(order).loc[celltype].get("diseased", 0))
        ax.text(i, max(healthy.loc[celltype], diseased.loc[celltype]) * 1.05, f"{d_pat}/{h_pat} pts", ha="center", fontsize=5.8)
    ax.set_xticks(x)
    ax.set_xticklabels(order)
    ax.set_ylabel("Cells contributing to pseudobulk")
    ax.set_title("Pseudobulk support by cell type", loc="left", fontsize=8, pad=6)
    ax.legend(loc="upper left", fontsize=6, ncol=2)
    ax.set_ylim(0, max(diseased.max(), healthy.max()) * 1.34)
    clean_axis(ax)

    save_figure(fig, "fig3_scrna_bcell_localization")


def evidence_category(row: pd.Series, fulltext_use: str | float | None) -> str:
    if isinstance(fulltext_use, str) and fulltext_use == "main_figure_candidate":
        return "B-cell mechanism"
    if isinstance(fulltext_use, str) and fulltext_use == "cytokine_context_candidate":
        return "Cytokine context"
    if isinstance(row["evidence_tier"], str) and row["evidence_tier"].startswith("A_"):
        return "Bulk + B-cell"
    if isinstance(row["evidence_tier"], str) and row["evidence_tier"].startswith("D_"):
        return "B-cell increased"
    return "Supplementary context"


def build_figure_4() -> None:
    bridge = read_tsv("results/m9_target_cell_bridge/candidate_target_cell_bridge_summary.tsv")
    cur = read_tsv("results/m9_target_cell_bridge/curcumin_bcell_target_evidence.tsv")
    contrast = read_tsv("results/m9_target_cell_bridge/gse125527_target_gene_contrast.tsv")
    m16 = read_tsv("results/m16_curcumin_fulltext_evidence/curcumin_high_value_target_fulltext_audit.tsv")

    bridge_source = bridge[
        [
            "ingredient_id",
            "ingredient_name",
            "bridge_score",
            "n_disease_context_unique_targets",
            "n_b_cell_expressed_targets",
            "n_b_cell_disease_increased_targets",
            "n_mdc_disease_increased_targets",
        ]
    ].copy()
    bridge_source["ingredient_name"] = bridge_source["ingredient_name"].map(display_ingredient_name)
    save_source(bridge_source, "fig4a_candidate_bridge_ranking.tsv")

    m16_use = (
        m16.sort_values("recommended_figure_use")
        .drop_duplicates(["gene_symbol"], keep="first")
        [["gene_symbol", "evidence_level", "recommended_figure_use", "pmid"]]
    )
    evidence = cur.merge(m16_use, on="gene_symbol", how="left")
    evidence["display_category"] = evidence.apply(lambda r: evidence_category(r, r.get("recommended_figure_use")), axis=1)
    evidence["has_b_cell_contrast"] = evidence["b_delta_log1p_cpm_diseased_minus_healthy"].notna()
    evidence_source = evidence[
        [
            "gene_symbol",
            "evidence_tier",
            "bulk_module_role",
            "b_delta_log1p_cpm_diseased_minus_healthy",
            "b_p_log1p_cpm_mannwhitney",
            "mdc_delta_log1p_cpm_diseased_minus_healthy",
            "pubmed_ids",
            "evidence_level",
            "recommended_figure_use",
            "display_category",
            "has_b_cell_contrast",
        ]
    ].copy()
    save_source(evidence_source, "fig4b_curcumin_target_evidence_tiers.tsv")

    selected_genes = [
        "BCL6",
        "BLNK",
        "SYK",
        "CCL2",
        "IL33",
        "IL1B",
        "TNF",
        "JAK1",
        "PIAS1",
        "STAT5A",
        "IL7",
        "IL15",
    ]
    heat = contrast[
        (contrast["celltype"].isin(["B", "M/DC"])) & (contrast["gene_symbol"].isin(selected_genes))
    ][
        [
            "celltype",
            "gene_symbol",
            "delta_log1p_cpm_diseased_minus_healthy",
            "p_log1p_cpm_mannwhitney",
            "mean_log1p_cpm_diseased",
            "mean_log1p_cpm_healthy",
        ]
    ].copy()
    save_source(heat, "fig4c_curcumin_target_cell_contrast.tsv")

    fig = plt.figure(figsize=(7.2, 6.4), constrained_layout=True)
    axes = fig.subplot_mosaic([["A", "B"], ["A", "C"]], width_ratios=[1.0, 1.45], height_ratios=[1.0, 1.1])

    ax = axes["A"]
    panel_label(ax, "a")
    top = bridge_source.nlargest(10, "bridge_score").sort_values("bridge_score", ascending=True)
    y = np.arange(len(top))
    y_labels = make_unique_labels(top["ingredient_name"])
    colors = [PALETTE["curcumin"] if name == "Curcumin" else "#B8BEC5" for name in top["ingredient_name"]]
    ax.barh(y, top["bridge_score"], color=colors, height=0.65)
    ax.set_yticks(y)
    ax.set_yticklabels(y_labels)
    for i, (_, row) in enumerate(top.iterrows()):
        if row["ingredient_name"] == "Curcumin":
            ax.text(
                row["bridge_score"] + 1.2,
                i,
                "22 disease-context targets\n12 expressed in rectal B cells",
                va="center",
                fontsize=5.8,
                color=PALETTE["boundary"],
            )
    ax.set_xlabel("Cell-state linkage score")
    ax.set_title("Candidate bridge ranking", loc="left", fontsize=8, pad=6)
    ax.set_xlim(0, top["bridge_score"].max() * 1.55)
    clean_axis(ax)

    ax = axes["B"]
    panel_label(ax, "b")
    category_order = ["B-cell mechanism", "Cytokine context", "Bulk + B-cell", "B-cell increased", "Supplementary context"]
    display_genes = ["BCL6", "BLNK", "SYK", "CCL2", "IL33", "IL1B", "TNF", "JAK1", "PIAS1", "STAT5A", "IL7", "IL15"]
    genes_by_priority = evidence_source.set_index("gene_symbol").reindex(display_genes).dropna(how="all").reset_index()
    columns = [
        ("bulk_up", "Bulk\nup"),
        ("b_contrast", "B-cell\ncontrast"),
        ("main_fulltext", "Main\nfull text"),
        ("cytokine_context", "Cytokine\ncontext"),
        ("supp_fulltext", "Supp.\nfull text"),
    ]
    genes_by_priority["bulk_up"] = genes_by_priority["bulk_module_role"].eq("bulk_up_top150")
    genes_by_priority["b_contrast"] = genes_by_priority["has_b_cell_contrast"].fillna(False)
    genes_by_priority["main_fulltext"] = genes_by_priority["recommended_figure_use"].eq("main_figure_candidate")
    genes_by_priority["cytokine_context"] = genes_by_priority["recommended_figure_use"].eq("cytokine_context_candidate")
    genes_by_priority["supp_fulltext"] = genes_by_priority["recommended_figure_use"].isin(["supplementary_candidate", "main_or_supplementary"])
    y = np.arange(len(genes_by_priority))
    ax.set_xlim(-0.5, len(columns) + 1.05)
    ax.set_ylim(-0.5, len(genes_by_priority) - 0.5)
    for j, (col, _) in enumerate(columns):
        for i, (_, row) in enumerate(genes_by_priority.iterrows()):
            active = bool(row[col]) if not pd.isna(row[col]) else False
            ax.scatter(
                j,
                i,
                s=54 if active else 18,
                color=PALETTE["light"] if not active else {
                    "bulk_up": PALETTE["disease"],
                    "b_contrast": PALETTE["bcell"],
                    "main_fulltext": PALETTE["curcumin"],
                    "cytokine_context": PALETTE["myeloid"],
                    "supp_fulltext": PALETTE["boundary"],
                }[col],
                edgecolor="white" if active else PALETTE["grid"],
                linewidth=0.4,
                zorder=3,
            )
    layer_colors = {
        "B-cell mechanism": PALETTE["curcumin"],
        "Cytokine context": PALETTE["myeloid"],
        "Bulk + B-cell": PALETTE["disease"],
        "B-cell increased": PALETTE["bcell"],
        "Supplementary context": PALETTE["boundary"],
    }
    for i, (_, row) in enumerate(genes_by_priority.iterrows()):
        category = row["display_category"]
        ax.text(
            len(columns) - 0.05,
            i,
            category,
            ha="left",
            va="center",
            fontsize=5.7,
            color=layer_colors.get(category, PALETTE["boundary"]),
        )
    ax.set_yticks(y)
    ax.set_yticklabels(genes_by_priority["gene_symbol"], fontsize=6)
    ax.set_xticks(np.arange(len(columns)))
    ax.set_xticklabels([label for _, label in columns], fontsize=6)
    ax.invert_yaxis()
    ax.set_title("Curcumin target evidence layers", loc="left", fontsize=8, pad=6)
    ax.tick_params(length=0)
    for spine in ax.spines.values():
        spine.set_visible(False)

    ax = axes["C"]
    panel_label(ax, "c")
    heat_genes = [g for g in selected_genes if g in set(heat["gene_symbol"])]
    heat_matrix = (
        heat.pivot(index="gene_symbol", columns="celltype", values="delta_log1p_cpm_diseased_minus_healthy")
        .reindex(heat_genes)
        .reindex(columns=["B", "M/DC"])
    )
    masked = np.ma.masked_invalid(heat_matrix.to_numpy(dtype=float))
    cmap = mpl.colors.LinearSegmentedColormap.from_list("delta", [PALETTE["bcell"], "white", PALETTE["disease"]])
    im = ax.imshow(masked, aspect="auto", cmap=cmap, vmin=-0.55, vmax=0.55)
    ax.set_xticks(np.arange(heat_matrix.shape[1]))
    ax.set_xticklabels(heat_matrix.columns)
    ax.set_yticks(np.arange(heat_matrix.shape[0]))
    ax.set_yticklabels(heat_matrix.index, fontsize=6)
    for i in range(heat_matrix.shape[0]):
        for j in range(heat_matrix.shape[1]):
            val = heat_matrix.iloc[i, j]
            if pd.isna(val):
                ax.text(j, i, "NA", ha="center", va="center", fontsize=5.5, color=PALETTE["boundary"])
            else:
                ax.text(j, i, f"{val:.2f}", ha="center", va="center", fontsize=5.5, color=PALETTE["text"])
    ax.set_title("Rectal target expression delta", loc="left", fontsize=8, pad=6)
    cbar = fig.colorbar(im, ax=ax, fraction=0.038, pad=0.03)
    cbar.set_label("Diseased minus healthy log1p CPM", fontsize=6)
    cbar.ax.tick_params(labelsize=5.5, length=2)
    ax.text(
        0,
        -0.18,
        "CCL2 and IL33 are cytokine-context targets without rectal B-cell contrast in this matrix.",
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=5.8,
        color=PALETTE["boundary"],
    )
    ax.tick_params(length=0)
    for spine in ax.spines.values():
        spine.set_visible(False)

    save_figure(fig, "fig4_curcumin_target_bridge")


def build_figure_5() -> None:
    ot_summary = read_tsv("results/m10_open_targets/candidate_open_targets_summary.tsv")
    m12_articles = read_tsv("results/m12_pubmed_verification/curcumin_pubmed_article_verification.tsv")
    m12_edges = read_tsv("results/m12_pubmed_verification/curcumin_target_edge_pubmed_verification.tsv")
    m16 = read_tsv("results/m16_curcumin_fulltext_evidence/curcumin_high_value_target_fulltext_audit.tsv")
    etcm2 = read_tsv("results/m15_etcm2_target_mapping/candidate_etcm2_mapped_overlap_summary.tsv")
    lincs = read_tsv("results/m5_lincs_reversal/candidate_reversal_status.tsv")
    lincs_summary = read_tsv("results/m5_lincs_reversal/summary.tsv")

    cur_ot = ot_summary[ot_summary["ingredient_name"] == "Curcumin"].iloc[0]
    ot_source = pd.DataFrame(
        [
            {
                "metric": "Disease-context targets",
                "count": int(cur_ot["n_disease_context_targets"]),
                "denominator": int(cur_ot["n_disease_context_targets"]),
                "fraction": 1.0,
                "display_layer": "HERB/M9 input",
            },
            {
                "metric": "Open Targets-supported",
                "count": int(cur_ot["n_open_targets_supported_genes"]),
                "denominator": int(cur_ot["n_disease_context_targets"]),
                "fraction": cur_ot["n_open_targets_supported_genes"] / cur_ot["n_disease_context_targets"],
                "display_layer": "disease relevance",
            },
            {
                "metric": "Genetic-supported",
                "count": int(cur_ot["n_genetic_supported_genes"]),
                "denominator": int(cur_ot["n_disease_context_targets"]),
                "fraction": cur_ot["n_genetic_supported_genes"] / cur_ot["n_disease_context_targets"],
                "display_layer": "genetic evidence",
            },
            {
                "metric": "Clinical-supported",
                "count": int(cur_ot["n_clinical_supported_genes"]),
                "denominator": int(cur_ot["n_disease_context_targets"]),
                "fraction": cur_ot["n_clinical_supported_genes"] / cur_ot["n_disease_context_targets"],
                "display_layer": "clinical evidence",
            },
        ]
    )
    save_source(ot_source, "fig5a_open_targets_curcumin_support.tsv")

    article_short = {
        33597887: "Memory T-cell colitis",
        36196887: "Breg/TLR-MyD88 colitis",
        36353208: "Bcl-6/Syk/BLNK colitis",
    }
    edge_counts = (
        m12_edges.groupby("pmid", as_index=False)
        .agg(
            herb_edges=("gene_symbol", "size"),
            title_matches=("title_match", lambda x: int((x == "yes").sum())),
            target_terms_in_record=("target_term_found_in_pubmed_title_or_abstract", lambda x: int((x == "yes").sum())),
        )
    )
    fulltext_counts = (
        m16.groupby("pmid", as_index=False)
        .agg(
            fulltext_body_matches=("gene_symbol", "size"),
            main_figure_candidates=("recommended_figure_use", lambda x: int((x == "main_figure_candidate").sum())),
            cytokine_context_candidates=("recommended_figure_use", lambda x: int((x == "cytokine_context_candidate").sum())),
            supplementary_candidates=("recommended_figure_use", lambda x: int(x.isin(["supplementary_candidate", "main_or_supplementary"]).sum())),
        )
    )
    pubmed_source = (
        m12_articles[["pmid", "publication_year", "journal", "all_herb_titles_match_pubmed", "abstract_available"]]
        .merge(edge_counts, on="pmid", how="left")
        .merge(fulltext_counts, on="pmid", how="left")
        .fillna(0)
    )
    pubmed_source["article_short"] = pubmed_source["pmid"].map(article_short)
    save_source(pubmed_source, "fig5b_pubmed_fulltext_verification.tsv")

    etcm_source = etcm2[
        [
            "ingredient_name",
            "n_accepted_mapped_gene_symbols_exact_high",
            "n_overlap_with_herb_m9_disease_context_targets",
            "n_overlap_with_bulk_aligned_targets",
            "n_overlap_with_open_targets_genetic_supported",
            "accepted_mapped_gene_symbols_exact_high",
            "overlap_with_herb_m9_disease_context_targets",
        ]
    ].copy()
    save_source(etcm_source, "fig5c_etcm2_overlap_boundary.tsv")

    lincs_source = lincs[
        [
            "ingredient_name",
            "priority_score",
            "candidate_pert_iname",
            "candidate_lincs_sources",
            "hit_status",
            "best_rank",
            "best_score",
        ]
    ].copy()
    lincs_source["ingredient_name"] = lincs_source["ingredient_name"].map(display_ingredient_name)
    lincs_source = lincs_source.rename(columns={"priority_score": "initial_lincs_screen_score"})
    save_source(lincs_source, "fig5d_lincs_candidate_reversal_status.tsv")
    save_source(lincs_summary, "fig5d_lincs_query_summary.tsv")

    fig = plt.figure(figsize=(7.2, 6.2), constrained_layout=True)
    axes = fig.subplot_mosaic([["A", "B"], ["C", "D"]], width_ratios=[1.05, 1.25], height_ratios=[1.05, 1.0])

    ax = axes["A"]
    panel_label(ax, "a")
    ot_plot = ot_source.iloc[::-1].reset_index(drop=True)
    y = np.arange(len(ot_plot))
    colors = [PALETTE["boundary"], PALETTE["myeloid"], PALETTE["bcell"], PALETTE["curcumin"]][::-1]
    ax.barh(y, ot_plot["count"], color=colors, height=0.62)
    for i, row in ot_plot.iterrows():
        ax.text(
            row["count"] + 0.35,
            i,
            f"{int(row['count'])}/{int(row['denominator'])}",
            va="center",
            fontsize=6,
            color=PALETTE["boundary"],
        )
    ax.set_yticks(y)
    ax.set_yticklabels(ot_plot["metric"], fontsize=6)
    ax.set_xlabel("Curcumin target genes")
    ax.set_xlim(0, int(cur_ot["n_disease_context_targets"]) + 3)
    ax.set_title("Open Targets disease relevance", loc="left", fontsize=8, pad=6)
    ax.text(
        0.02,
        0.04,
        f"Max overall score={cur_ot['max_overall_score']:.3f}",
        transform=ax.transAxes,
        fontsize=5.8,
        color=PALETTE["boundary"],
    )
    clean_axis(ax)

    ax = axes["B"]
    panel_label(ax, "b")
    pubmed_plot = pubmed_source.sort_values("publication_year").reset_index(drop=True)
    rows = pubmed_plot["article_short"].tolist()
    columns = [
        ("title_matches", "Title\nmatch", PALETTE["healthy"]),
        ("target_terms_in_record", "Target in\nrecord", PALETTE["bcell"]),
        ("fulltext_body_matches", "Full-text\nbody", PALETTE["curcumin"]),
        ("main_figure_candidates", "Main\nmechanism", PALETTE["curcumin"]),
        ("cytokine_context_candidates", "Cytokine\ncontext", PALETTE["myeloid"]),
        ("supplementary_candidates", "Supp.\ncontext", PALETTE["boundary"]),
    ]
    ax.set_xlim(-0.5, len(columns) - 0.5)
    ax.set_ylim(-0.5, len(rows) - 0.5)
    for j, (col, _, color) in enumerate(columns):
        for i, row in pubmed_plot.iterrows():
            value = int(row[col])
            ax.scatter(
                j,
                i,
                s=22 + value * 22,
                color=color if value else PALETTE["light"],
                edgecolor="white" if value else PALETTE["grid"],
                linewidth=0.4,
            )
            if value:
                ax.text(j, i, str(value), ha="center", va="center", fontsize=5.5, color="white")
    ax.set_yticks(np.arange(len(rows)))
    ax.set_yticklabels([f"{int(pmid)}\n{title}" for pmid, title in zip(pubmed_plot["pmid"], rows)], fontsize=5.8)
    ax.set_xticks(np.arange(len(columns)))
    ax.set_xticklabels([label for _, label, _ in columns], fontsize=5.8)
    ax.invert_yaxis()
    ax.set_title("PubMed and full-text evidence tiers", loc="left", fontsize=8, pad=6)
    ax.tick_params(length=0)
    for spine in ax.spines.values():
        spine.set_visible(False)

    ax = axes["C"]
    panel_label(ax, "c")
    etcm_plot = etcm_source.sort_values("n_accepted_mapped_gene_symbols_exact_high", ascending=True)
    y = np.arange(len(etcm_plot))
    ax.barh(
        y,
        etcm_plot["n_accepted_mapped_gene_symbols_exact_high"],
        color="#C8CED4",
        height=0.62,
        label="Accepted mapped genes",
    )
    ax.barh(
        y,
        etcm_plot["n_overlap_with_herb_m9_disease_context_targets"],
        color=PALETTE["disease"],
        height=0.30,
        label="HERB/M9 overlap",
    )
    for i, row in etcm_plot.reset_index(drop=True).iterrows():
        if row["ingredient_name"] == "Curcumin":
            ax.text(
                row["n_accepted_mapped_gene_symbols_exact_high"] + 0.8,
                i,
                "Curcumin: 5 mapped,\n0 disease-context overlap",
                va="center",
                fontsize=5.8,
                color=PALETTE["boundary"],
            )
    ax.set_yticks(y)
    ax.set_yticklabels(etcm_plot["ingredient_name"], fontsize=6)
    ax.set_xlabel("ETCM2 mapped target genes")
    ax.set_title("ETCM2 cross-resource concordance", loc="left", fontsize=8, pad=6)
    ax.legend(loc="lower right", fontsize=5.8)
    ax.set_xlim(0, etcm_plot["n_accepted_mapped_gene_symbols_exact_high"].max() * 1.45)
    clean_axis(ax)

    ax = axes["D"]
    panel_label(ax, "d")
    lincs_plot = lincs_source.sort_values("initial_lincs_screen_score", ascending=True).tail(12)
    y = np.arange(len(lincs_plot))
    y_labels = make_unique_labels(lincs_plot["ingredient_name"])
    colors = [PALETTE["curcumin"] if name == "Curcumin" else "#B8BEC5" for name in lincs_plot["ingredient_name"]]
    ax.barh(y, lincs_plot["initial_lincs_screen_score"], color=colors, height=0.62)
    ax.set_yticks(y)
    ax.set_yticklabels(y_labels)
    n_tests = len(lincs_source)
    n_hits = int((lincs_source["hit_status"] != "no_top_result_hit").sum())
    ax.text(
        0.98,
        0.08,
        f"{n_hits}/{n_tests} exact-ID candidate top-result hits",
        transform=ax.transAxes,
        ha="right",
        va="bottom",
        fontsize=6.2,
        color=PALETTE["boundary"],
    )
    ax.set_xlabel("Initial LINCS-screen score")
    ax.set_title("L1000CDS2 exact-ID coverage", loc="left", fontsize=8, pad=6)
    ax.set_xlim(0, lincs_plot["initial_lincs_screen_score"].max() * 1.18)
    clean_axis(ax)

    save_figure(fig, "fig5_evidence_boundaries")


def main() -> None:
    setup_style()
    build_figure_1()
    build_figure_2()
    build_figure_3()
    build_figure_4()
    build_figure_5()
    print(f"Wrote source data to {SOURCE_DIR.relative_to(ROOT)}")
    print(f"Wrote figure outputs to {OUTPUT_DIR.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
