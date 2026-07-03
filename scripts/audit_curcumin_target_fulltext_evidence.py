#!/usr/bin/env python3
"""Audit high-value Curcumin target evidence against PubMed/PMC text."""

from __future__ import annotations

import argparse
import csv
import re
import textwrap
import xml.etree.ElementTree as ET
from pathlib import Path
from typing import Any, Dict, Iterable, List, Sequence, Tuple


OUT_DIR = Path("results/m16_curcumin_fulltext_evidence")
DOC_PATH = Path("docs/M16_CURCUMIN_FULLTEXT_EVIDENCE_AUDIT.md")

PRIMARY_TARGETS = ["IL1B", "CCL2", "IL33", "BCL6", "BLNK", "SYK", "TNF"]
SECONDARY_TARGETS = ["JAK1", "PIAS1", "IL7", "IL15", "STAT5A"]

TARGET_PATTERNS: Dict[str, List[str]] = {
    "IL1B": [r"\bIL[- ]?1(?:β|beta|B)\b", r"\binterleukin[- ]?1(?:β|beta)\b"],
    "CCL2": [r"\bCCL[- ]?2\b", r"\bMCP[- ]?1\b"],
    "IL33": [r"\bIL[- ]?33\b", r"\binterleukin[- ]?33\b"],
    "BCL6": [r"\bBcl[- ]?6\b", r"\bBCL6\b"],
    "BLNK": [r"\bBLNK\b", r"\bB cell linker\b"],
    "SYK": [r"\bSyk\b", r"\bSYK\b", r"\bp[- ]?Syk\b"],
    "TNF": [r"\bTNF[- ]?α\b", r"\bTNF[- ]?alpha\b", r"\bTNF\b"],
    "JAK1": [r"\bJAK1\b", r"\bJAK[- ]?1\b"],
    "PIAS1": [r"\bPIAS1\b"],
    "IL7": [r"\bIL[- ]?7\b", r"\binterleukin[- ]?7\b"],
    "IL15": [r"\bIL[- ]?15\b", r"\binterleukin[- ]?15\b"],
    "STAT5A": [r"\bSTAT5A\b", r"\bSTAT5\b", r"\bp[- ]?STAT5\b"],
}


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


def clean_text(value: str) -> str:
    return re.sub(r"\s+", " ", value).strip()


def elem_text(elem: ET.Element | None) -> str:
    if elem is None:
        return ""
    return clean_text(" ".join(elem.itertext()))


def parse_pubmed(path: Path) -> Dict[str, Dict[str, str]]:
    root = ET.parse(path).getroot()
    records: Dict[str, Dict[str, str]] = {}
    for article in root.findall(".//PubmedArticle"):
        pmid = article.findtext("./PubmedData/ArticleIdList/ArticleId[@IdType='pubmed']", "")
        title = elem_text(article.find(".//ArticleTitle"))
        abstract = clean_text(" ".join(elem_text(item) for item in article.findall(".//AbstractText")))
        doi = article.findtext("./PubmedData/ArticleIdList/ArticleId[@IdType='doi']", "")
        pmc = article.findtext("./PubmedData/ArticleIdList/ArticleId[@IdType='pmc']", "")
        journal = elem_text(article.find(".//Journal/Title"))
        year = article.findtext(".//PubDate/Year", "")
        if pmid:
            records[pmid] = {
                "pmid": pmid,
                "title": title,
                "abstract": abstract,
                "doi": doi,
                "pmc_id": pmc,
                "journal": journal,
                "year": year,
            }
    return records


def section_title(sec: ET.Element) -> str:
    title = sec.find("title")
    return elem_text(title)


def parse_pmc(path: Path) -> Dict[str, Any]:
    root = ET.parse(path).getroot()
    article = root.find(".//article")
    if article is None:
        return {}
    pmid = elem_text(article.find(".//article-id[@pub-id-type='pmid']"))
    pmc_id = elem_text(article.find(".//article-id[@pub-id-type='pmc']"))
    title = elem_text(article.find(".//article-title"))
    abstract = clean_text(" ".join(elem_text(item) for item in article.findall(".//abstract//p")))
    sections = []
    body = article.find("body")
    if body is not None:
        for sec in body.findall(".//sec"):
            title_text = section_title(sec)
            paragraphs = [elem_text(p) for p in sec.findall("p")]
            text = clean_text(" ".join(paragraphs))
            if text:
                sections.append({"title": title_text, "text": text})
    return {
        "pmid": pmid,
        "pmc_id": f"PMC{pmc_id}" if pmc_id and not pmc_id.startswith("PMC") else pmc_id,
        "title": title,
        "abstract": abstract,
        "sections": sections,
    }


def load_pmc_records(paths: Iterable[Path]) -> Dict[str, Dict[str, Any]]:
    records: Dict[str, Dict[str, Any]] = {}
    for path in paths:
        if not path.exists():
            continue
        record = parse_pmc(path)
        if record.get("pmid"):
            records[str(record["pmid"])] = record
    return records


def parse_private_text_arg(value: str) -> Tuple[str, Path]:
    if "=" not in value:
        raise ValueError(f"Private full-text input must be formatted as PMID=path, got: {value}")
    pmid, path = value.split("=", 1)
    return pmid.strip(), Path(path.strip())


def load_private_fulltext_records(values: Iterable[str]) -> Dict[str, Dict[str, Any]]:
    records: Dict[str, Dict[str, Any]] = {}
    for value in values:
        pmid, path = parse_private_text_arg(value)
        if not path.exists():
            continue
        text = clean_text(path.read_text(encoding="utf-8", errors="replace"))
        records[pmid] = {
            "pmid": pmid,
            "source_type": "user_provided_private_fulltext",
            "sections": [
                {
                    "title": "private_fulltext_extracted_text",
                    "text": text,
                }
            ],
        }
    return records


def find_matches(text: str, patterns: Sequence[str]) -> List[str]:
    hits = []
    for pattern in patterns:
        if re.search(pattern, text, flags=re.IGNORECASE):
            hits.append(pattern)
    return hits


def first_snippet(text: str, patterns: Sequence[str], width: int = 180) -> str:
    for pattern in patterns:
        match = re.search(pattern, text, flags=re.IGNORECASE)
        if match:
            start = max(match.start() - width // 2, 0)
            end = min(match.end() + width // 2, len(text))
            return clean_text(text[start:end])
    return ""


def evidence_level(
    pmid: str,
    gene: str,
    pubmed: Dict[str, Dict[str, str]],
    pmc: Dict[str, Dict[str, Any]],
    private_fulltext: Dict[str, Dict[str, Any]],
) -> Tuple[str, str, str, str, str, str]:
    patterns = TARGET_PATTERNS.get(gene, [rf"\b{re.escape(gene)}\b"])
    pmc_record = pmc.get(pmid)
    private_record = private_fulltext.get(pmid)
    pub_record = pubmed.get(pmid, {})
    abstract_text = pmc_record.get("abstract", "") if pmc_record else pub_record.get("abstract", "")
    abstract_hits = find_matches(abstract_text, patterns)
    section_hits: List[str] = []
    section_snippet = ""
    fulltext_record = pmc_record or private_record
    if fulltext_record:
        for sec in fulltext_record.get("sections", []):
            hits = find_matches(sec["text"], patterns)
            if hits:
                section_hits.append(sec["title"] or "untitled_section")
                if not section_snippet:
                    section_snippet = first_snippet(sec["text"], patterns)
        if section_hits:
            level = "open_fulltext_body_match" if pmc_record else "user_provided_fulltext_body_match"
        elif abstract_hits:
            level = "open_fulltext_abstract_match_only" if pmc_record else "user_provided_fulltext_abstract_match_only"
        else:
            level = "open_fulltext_no_target_term_match" if pmc_record else "user_provided_fulltext_no_target_term_match"
        status = "open_pmc_fulltext" if pmc_record else "user_provided_private_fulltext"
    elif abstract_hits:
        level = "pubmed_abstract_match_no_public_fulltext"
        status = "pubmed_only_or_publisher_restricted"
    else:
        level = "pubmed_record_no_target_term_match"
        status = "pubmed_only_or_publisher_restricted"
    abstract_snippet = first_snippet(abstract_text, patterns)
    return (
        status,
        level,
        ";".join(sorted(set(abstract_hits))),
        ";".join(sorted(set(section_hits))),
        abstract_snippet,
        section_snippet,
    )


def interpretation_for(gene: str, pmid: str, level: str) -> Tuple[str, str]:
    if level in {"open_fulltext_body_match", "user_provided_fulltext_body_match"}:
        if gene in {"BCL6", "BLNK", "SYK"} and pmid == "36353208":
            return "main_figure_candidate", "Directly supports the memory B-cell/Bcl-6-Syk-BLNK axis claim."
        if gene in {"CCL2", "IL33", "IL1B", "TNF"} and pmid == "36196887":
            return "cytokine_context_candidate", "Full text supports cytokine-context evidence in the chronic colitis/Breg study; use as inflammatory-context support unless the claim is explicitly cell-specific."
        if gene in {"IL1B", "TNF"}:
            return "main_or_supplementary", "Supports inflammatory target context; use with disease/bulk/scRNA evidence rather than as standalone mechanism."
        if gene in SECONDARY_TARGETS:
            return "supplementary_candidate", "Supports secondary signaling/cytokine context; do not center the main B-cell mechanism on this target."
        return "supplementary_candidate", "Open full text contains target-context evidence."
    if level == "pubmed_abstract_match_no_public_fulltext":
        if gene in {"CCL2", "IL33"}:
            return "supplementary_until_fulltext_checked", "PubMed abstract supports cytokine modulation, but publisher full text was not publicly accessible in this run."
        return "supplementary_until_fulltext_checked", "Use only as abstract-level support until full text is inspected."
    return "do_not_use_for_primary_claim", "Target term was not located at the audited text level."


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--curcumin-table", default="results/m9_target_cell_bridge/curcumin_bcell_target_evidence.tsv")
    parser.add_argument("--pubmed-xml", default="data/raw/pubmed/curcumin_pubmed_33597887_36196887_36353208.xml")
    parser.add_argument("--pmc-xml", action="append", default=[
        "data/raw/pubmed/curcumin_pmc_PMC7882737.xml",
        "data/raw/pubmed/curcumin_pmc_PMC9639655.xml",
    ])
    parser.add_argument(
        "--private-fulltext-text",
        action="append",
        default=[
            "36196887=data/private/fulltext/PMID36196887_PTR_2023_curcumin_regulatory_B_cells.txt",
        ],
        help="User-provided full-text extraction formatted as PMID=path. These files should remain outside git.",
    )
    parser.add_argument("--out-dir", default=str(OUT_DIR))
    args = parser.parse_args()

    out_dir = Path(args.out_dir)
    curcumin_rows = read_tsv(Path(args.curcumin_table))
    pubmed = parse_pubmed(Path(args.pubmed_xml))
    pmc = load_pmc_records(Path(path) for path in args.pmc_xml)
    private_fulltext = load_private_fulltext_records(args.private_fulltext_text)
    target_set = set(PRIMARY_TARGETS + SECONDARY_TARGETS)
    rows: List[Dict[str, Any]] = []
    for cur_row in curcumin_rows:
        gene = cur_row["gene_symbol"]
        if gene not in target_set:
            continue
        for pmid in [value for value in cur_row["pubmed_ids"].split(";") if value]:
            pub_record = pubmed.get(pmid, {})
            status, level, abstract_hits, section_hits, abstract_snippet, body_snippet = evidence_level(pmid, gene, pubmed, pmc, private_fulltext)
            figure_use, interpretation = interpretation_for(gene, pmid, level)
            rows.append(
                {
                    "gene_symbol": gene,
                    "target_priority": "primary" if gene in PRIMARY_TARGETS else "secondary",
                    "pmid": pmid,
                    "doi": pub_record.get("doi", ""),
                    "pmc_id": pub_record.get("pmc_id", ""),
                    "article_title": pub_record.get("title", ""),
                    "journal": pub_record.get("journal", ""),
                    "year": pub_record.get("year", ""),
                    "fulltext_status": status,
                    "evidence_level": level,
                    "abstract_matched_patterns": abstract_hits,
                    "body_matched_sections": section_hits,
                    "abstract_snippet": abstract_snippet,
                    "body_snippet": body_snippet,
                    "m9_evidence_tier": cur_row["evidence_tier"],
                    "m9_bulk_module_role": cur_row["bulk_module_role"],
                    "m9_b_cell_delta": cur_row["b_delta_log1p_cpm_diseased_minus_healthy"],
                    "m9_b_cell_p": cur_row["b_p_log1p_cpm_mannwhitney"],
                    "recommended_figure_use": figure_use,
                    "interpretation_boundary": interpretation,
                }
            )
    rows.sort(key=lambda row: (row["target_priority"] != "primary", row["gene_symbol"], row["pmid"]))
    fields = [
        "gene_symbol",
        "target_priority",
        "pmid",
        "doi",
        "pmc_id",
        "article_title",
        "journal",
        "year",
        "fulltext_status",
        "evidence_level",
        "abstract_matched_patterns",
        "body_matched_sections",
        "abstract_snippet",
        "body_snippet",
        "m9_evidence_tier",
        "m9_bulk_module_role",
        "m9_b_cell_delta",
        "m9_b_cell_p",
        "recommended_figure_use",
        "interpretation_boundary",
    ]
    write_tsv(out_dir / "curcumin_high_value_target_fulltext_audit.tsv", rows, fields)

    level_counts: Dict[str, int] = {}
    use_counts: Dict[str, int] = {}
    for row in rows:
        level_counts[row["evidence_level"]] = level_counts.get(row["evidence_level"], 0) + 1
        use_counts[row["recommended_figure_use"]] = use_counts.get(row["recommended_figure_use"], 0) + 1

    primary_rows = [row for row in rows if row["target_priority"] == "primary"]
    main_ready = [row for row in primary_rows if row["recommended_figure_use"] == "main_figure_candidate"]
    abstract_only = [row for row in rows if row["evidence_level"] == "pubmed_abstract_match_no_public_fulltext"]

    report = [
        "# M16 Curcumin Full-Text Evidence Audit",
        "",
        "## Scope",
        "",
        "This run audits high-value Curcumin target evidence against PubMed records and publicly accessible PMC full-text XML where available.",
        "",
        "## Key Output",
        "",
        "- `results/m16_curcumin_fulltext_evidence/curcumin_high_value_target_fulltext_audit.tsv`",
        "",
        "## Evidence Availability",
        "",
        "- PubMed records audited: 3.",
        f"- Public PMC full texts available in this run: {len(pmc)} (PMIDs: {';'.join(sorted(pmc))}).",
        f"- User-provided private full texts available in this run: {len(private_fulltext)} (PMIDs: {';'.join(sorted(private_fulltext))}).",
        "- PMID 36196887 has no PMC full text in the NCBI record, but a user-provided private PDF/text extraction is now available for local evidence auditing.",
        "",
        "## Evidence-Level Counts",
        "",
        "| Evidence level | Rows |",
        "|---|---:|",
    ]
    for key, value in sorted(level_counts.items()):
        report.append(f"| {key} | {value} |")
    report.extend(["", "## Figure-Use Counts", "", "| Recommended use | Rows |", "|---|---:|"])
    for key, value in sorted(use_counts.items()):
        report.append(f"| {key} | {value} |")
    report.extend(
        [
            "",
            "## Main-Figure Candidates After Full-Text Audit",
            "",
            "| Gene | PMID | Evidence level | Figure use | Boundary |",
            "|---|---:|---|---|---|",
        ]
    )
    for row in main_ready:
        report.append(
            f"| {row['gene_symbol']} | {row['pmid']} | {row['evidence_level']} | {row['recommended_figure_use']} | {row['interpretation_boundary']} |"
        )
    report.extend(
        [
            "",
            "## Abstract-Only Rows Requiring Caution",
            "",
            "| Gene | PMID | Reason |",
            "|---|---:|---|",
        ]
    )
    for row in abstract_only:
        report.append(f"| {row['gene_symbol']} | {row['pmid']} | {row['interpretation_boundary']} |")
    report.extend(
        [
            "",
            "## Interpretation Boundary",
            "",
            "- Open PMC full-text matches can support figure-level evidence annotations but still require careful wording because most data are mouse DSS-colitis experiments, not human therapeutic validation.",
            "- User-provided private full-text rows can support local evidence upgrading, but the copyrighted PDF/text itself should not be committed or redistributed.",
            "- The strongest B-cell mechanism evidence after this audit is the PMID 36353208 Bcl-6/Syk/BLNK memory B-cell axis.",
            "- Cytokine rows such as IL1B, TNF, CCL2, and IL33 from PMID 36196887 now have private-full-text support for the chronic colitis/Breg study context; draw them as cytokine-context support unless a final manual review confirms cell-specific Breg measurement.",
            "",
        ]
    )
    DOC_PATH.write_text("\n".join(report), encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
