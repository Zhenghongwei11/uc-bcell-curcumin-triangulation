#!/usr/bin/env python3
"""Verify high-value Curcumin target references against PubMed efetch XML."""

from __future__ import annotations

import argparse
import csv
import hashlib
import re
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from collections import defaultdict
from pathlib import Path
from typing import Any, Dict, Iterable, List, Sequence, Set


EFETCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
ACCESS_DATE = "2026-06-29"
RAW_DIR = Path("data/raw/pubmed")
OUT_DIR = Path("results/m12_pubmed_verification")
DOC_PATH = Path("docs/M12_PUBMED_REFERENCE_VERIFICATION_RUN.md")
MANIFEST_PATH = Path("data/manifest.tsv")
MANIFEST_FIELDS = [
    "dataset_id",
    "source",
    "version_or_date",
    "file_path",
    "source_url_or_api",
    "access_date",
    "file_size_bytes",
    "md5",
    "status",
    "notes",
]

TARGET_PATTERNS = {
    "BCL6": [r"\bBCL6\b", r"\bBcl-?6\b"],
    "BLNK": [r"\bBLNK\b", r"\bp-?BLNK\b"],
    "CCL2": [r"\bCCL-?2\b", r"\bMCP-?1\b"],
    "IFNG": [r"\bIFN-?gamma\b", r"\bIFN-?γ\b", r"\bIFNG\b"],
    "IL10": [r"\bIL-?10\b"],
    "IL12A": [r"\bIL-?12A\b", r"\bIL-?35\b"],
    "IL13": [r"\bIL-?13\b"],
    "IL15": [r"\bIL-?15\b"],
    "IL1B": [r"\bIL-?1beta\b", r"\bIL-?1β\b", r"\bIL1B\b"],
    "IL22": [r"\bIL-?22\b"],
    "IL33": [r"\bIL-?33\b"],
    "IL4": [r"\bIL-?4\b"],
    "IL6": [r"\bIL-?6\b"],
    "IL7": [r"\bIL-?7A?\b"],
    "INS": [r"\binsulin\b", r"\bINS\b"],
    "JAK1": [r"\bJAK1\b"],
    "PIAS1": [r"\bPIAS1\b"],
    "SH3KBP1": [r"\bSH3KBP1\b", r"\bCIN85\b"],
    "STAT5A": [r"\bSTAT5A?\b", r"\bp-?STAT5\b"],
    "SYK": [r"\bSYK\b", r"\bp-?SYK\b"],
    "TNF": [r"\bTNF-?alpha\b", r"\bTNF-?α\b", r"\bTNF\b"],
}

DISEASE_PATTERNS = {
    "curcumin": r"\bcurcumin\b|\bCur\b",
    "colitis": r"\bcolitis\b",
    "ulcerative_colitis": r"ulcerative colitis|\bUC\b",
    "dss": r"dextran sulfate sodium|\bDSS\b",
    "b_cell": r"\bB cells?\b|\bregulatory B cells?\b|\bBreg\b",
    "memory_cell": r"\bmemory [BT] cells?\b",
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


def md5_file(path: Path) -> str:
    digest = hashlib.md5()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def normalize_title(value: str) -> str:
    return re.sub(r"[^a-z0-9]+", "", value.lower())


def text_content(node: ET.Element | None) -> str:
    if node is None:
        return ""
    return " ".join("".join(node.itertext()).split())


def join_unique(values: Iterable[str]) -> str:
    return ";".join(sorted({value for value in values if value}))


def fetch_pubmed_xml(pmids: Sequence[str], timeout: int, reuse_existing: bool) -> Path:
    RAW_DIR.mkdir(parents=True, exist_ok=True)
    key = "_".join(pmids)
    raw_path = RAW_DIR / f"curcumin_pubmed_{key}.xml"
    if reuse_existing and raw_path.exists() and raw_path.stat().st_size > 100:
        return raw_path
    params = urllib.parse.urlencode({"db": "pubmed", "id": ",".join(pmids), "retmode": "xml"})
    request = urllib.request.Request(
        f"{EFETCH_URL}?{params}",
        headers={"User-Agent": "zyy-pubmed-verification/1.0"},
        method="GET",
    )
    with urllib.request.urlopen(request, timeout=timeout) as response:
        raw_path.write_bytes(response.read())
    return raw_path


def update_manifest(path: Path) -> None:
    existing: List[Dict[str, str]] = []
    if MANIFEST_PATH.exists():
        with MANIFEST_PATH.open(encoding="utf-8", newline="") as handle:
            existing = list(csv.DictReader(handle, delimiter="\t"))
    by_id = {row["dataset_id"]: row for row in existing if row.get("dataset_id")}
    by_id[f"pubmed_curcumin_reference_verification_{path.stem.replace('curcumin_pubmed_', '')}"] = {
        "dataset_id": f"pubmed_curcumin_reference_verification_{path.stem.replace('curcumin_pubmed_', '')}",
        "source": "NCBI PubMed efetch",
        "version_or_date": ACCESS_DATE,
        "file_path": str(path),
        "source_url_or_api": EFETCH_URL,
        "access_date": ACCESS_DATE,
        "file_size_bytes": str(path.stat().st_size),
        "md5": md5_file(path),
        "status": "analysis_ready",
        "notes": "PubMed XML used to verify high-value Curcumin HERB target references",
    }
    with MANIFEST_PATH.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=MANIFEST_FIELDS, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for key in sorted(by_id):
            writer.writerow(by_id[key])


def parse_pubmed_xml(path: Path) -> Dict[str, Dict[str, Any]]:
    root = ET.fromstring(path.read_text(encoding="utf-8", errors="replace"))
    articles: Dict[str, Dict[str, Any]] = {}
    for article in root.findall(".//PubmedArticle"):
        pmid = text_content(article.find(".//MedlineCitation/PMID"))
        art = article.find(".//Article")
        title = text_content(art.find("./ArticleTitle") if art is not None else None)
        journal = text_content(article.find(".//Journal/Title"))
        year = text_content(article.find(".//PubDate/Year")) or text_content(article.find(".//ArticleDate/Year"))
        abstract = " ".join(text_content(node) for node in article.findall(".//Abstract/AbstractText"))
        doi = ""
        pubmed_ids = article.find("./PubmedData/ArticleIdList")
        if pubmed_ids is not None:
            for article_id in pubmed_ids.findall("./ArticleId"):
                if article_id.attrib.get("IdType") == "doi":
                    doi = text_content(article_id)
        combined = f"{title} {abstract}"
        disease_terms = [name for name, pattern in DISEASE_PATTERNS.items() if re.search(pattern, combined, re.IGNORECASE)]
        articles[pmid] = {
            "pmid": pmid,
            "pubmed_url": f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/",
            "pubmed_title": title,
            "journal": journal,
            "publication_year": year,
            "doi": doi,
            "abstract_available": "yes" if abstract else "no",
            "disease_context_terms_in_pubmed_record": ";".join(disease_terms),
            "pubmed_combined_text": combined,
        }
    return articles


def select_high_value_edges(edge_path: Path, pmids: Set[str]) -> List[Dict[str, str]]:
    rows = []
    for row in read_tsv(edge_path):
        if row.get("ingredient_name", "").lower() != "curcumin":
            continue
        if row.get("pubmed_id", "") in pmids:
            rows.append(row)
    rows.sort(key=lambda row: (row["pubmed_id"], row["gene_symbol"]))
    return rows


def target_terms_in_text(gene: str, text: str) -> List[str]:
    patterns = TARGET_PATTERNS.get(gene.upper(), [rf"\b{re.escape(gene)}\b"])
    return [pattern for pattern in patterns if re.search(pattern, text, re.IGNORECASE)]


def build_article_rows(articles: Dict[str, Dict[str, Any]], herb_titles_by_pmid: Dict[str, Set[str]]) -> List[Dict[str, Any]]:
    rows = []
    for pmid, article in sorted(articles.items()):
        herb_titles = herb_titles_by_pmid.get(pmid, set())
        title_matches = [
            "yes" if normalize_title(title) == normalize_title(article["pubmed_title"]) else "no"
            for title in herb_titles
        ]
        rows.append(
            {
                "pmid": pmid,
                "pubmed_url": article["pubmed_url"],
                "pubmed_title": article["pubmed_title"],
                "herb_reference_titles": " || ".join(sorted(herb_titles)),
                "all_herb_titles_match_pubmed": "yes" if title_matches and all(value == "yes" for value in title_matches) else "no",
                "journal": article["journal"],
                "publication_year": article["publication_year"],
                "doi": article["doi"],
                "abstract_available": article["abstract_available"],
                "disease_context_terms_in_pubmed_record": article["disease_context_terms_in_pubmed_record"],
            }
        )
    return rows


def build_edge_rows(edges: Sequence[Dict[str, str]], articles: Dict[str, Dict[str, Any]]) -> List[Dict[str, Any]]:
    rows = []
    for edge in edges:
        pmid = edge["pubmed_id"]
        article = articles.get(pmid, {})
        combined = article.get("pubmed_combined_text", "")
        gene = edge["gene_symbol"].upper()
        matched_terms = target_terms_in_text(gene, combined)
        rows.append(
            {
                "gene_symbol": gene,
                "pmid": pmid,
                "pubmed_found": "yes" if article else "no",
                "herb_reference_title": edge.get("reference_title", ""),
                "pubmed_title": article.get("pubmed_title", ""),
                "title_match": "yes" if article and normalize_title(edge.get("reference_title", "")) == normalize_title(article.get("pubmed_title", "")) else "no",
                "relationship": edge.get("relationship", ""),
                "relationship_direction": edge.get("relationship_direction", ""),
                "therapeutic_alignment": edge.get("therapeutic_alignment", ""),
                "target_term_found_in_pubmed_title_or_abstract": "yes" if matched_terms else "no",
                "matched_target_patterns": ";".join(matched_terms),
                "disease_context_terms_in_pubmed_record": article.get("disease_context_terms_in_pubmed_record", ""),
                "supporting_sentence": edge.get("supporting_sentence", ""),
            }
        )
    return rows


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--curcumin-table", default="results/m9_target_cell_bridge/curcumin_bcell_target_evidence.tsv")
    parser.add_argument("--m9-edge-table", default="results/m9_target_cell_bridge/disease_context_herb_target_edges.tsv")
    parser.add_argument("--timeout", type=int, default=40)
    parser.add_argument("--no-reuse-existing", action="store_true")
    args = parser.parse_args()

    high_value_pmids: Set[str] = set()
    for row in read_tsv(Path(args.curcumin_table)):
        if row.get("evidence_tier", "")[:1] in {"A", "B", "C", "D"}:
            high_value_pmids.update(pmid for pmid in row.get("pubmed_ids", "").split(";") if pmid)
    pmids = sorted(high_value_pmids)
    raw_path = fetch_pubmed_xml(pmids, args.timeout, not args.no_reuse_existing)
    update_manifest(raw_path)

    articles = parse_pubmed_xml(raw_path)
    edges = select_high_value_edges(Path(args.m9_edge_table), high_value_pmids)
    herb_titles_by_pmid: Dict[str, Set[str]] = defaultdict(set)
    for edge in edges:
        herb_titles_by_pmid[edge["pubmed_id"]].add(edge.get("reference_title", ""))

    article_rows = build_article_rows(articles, herb_titles_by_pmid)
    edge_rows = build_edge_rows(edges, articles)

    article_fields = [
        "pmid",
        "pubmed_url",
        "pubmed_title",
        "herb_reference_titles",
        "all_herb_titles_match_pubmed",
        "journal",
        "publication_year",
        "doi",
        "abstract_available",
        "disease_context_terms_in_pubmed_record",
    ]
    edge_fields = [
        "gene_symbol",
        "pmid",
        "pubmed_found",
        "herb_reference_title",
        "pubmed_title",
        "title_match",
        "relationship",
        "relationship_direction",
        "therapeutic_alignment",
        "target_term_found_in_pubmed_title_or_abstract",
        "matched_target_patterns",
        "disease_context_terms_in_pubmed_record",
        "supporting_sentence",
    ]
    write_tsv(OUT_DIR / "curcumin_pubmed_article_verification.tsv", article_rows, article_fields)
    write_tsv(OUT_DIR / "curcumin_target_edge_pubmed_verification.tsv", edge_rows, edge_fields)

    n_title_matches = sum(1 for row in article_rows if row["all_herb_titles_match_pubmed"] == "yes")
    n_edge_title_matches = sum(1 for row in edge_rows if row["title_match"] == "yes")
    n_target_terms = sum(1 for row in edge_rows if row["target_term_found_in_pubmed_title_or_abstract"] == "yes")
    report_lines = [
        "# M12 PubMed Reference Verification Run",
        "",
        "## Scope",
        "",
        "This run independently verifies high-value Curcumin target references from M9 against official NCBI PubMed efetch XML records.",
        "",
        "## Key Outputs",
        "",
        "- `data/raw/pubmed/curcumin_pubmed_33597887_36196887_36353208.xml`",
        "- `results/m12_pubmed_verification/curcumin_pubmed_article_verification.tsv`",
        "- `results/m12_pubmed_verification/curcumin_target_edge_pubmed_verification.tsv`",
        "",
        "## Verification Summary",
        "",
        f"- High-value PMIDs checked: {len(pmids)} ({';'.join(pmids)}).",
        f"- PubMed records retrieved: {len(articles)}.",
        f"- Article-level HERB title matches: {n_title_matches}/{len(article_rows)}.",
        f"- Edge-level title matches: {n_edge_title_matches}/{len(edge_rows)}.",
        f"- Edge-level target terms found in PubMed title/abstract: {n_target_terms}/{len(edge_rows)}.",
        "",
        "## Article-Level Verification",
        "",
        "| PMID | Title match | Journal | Year | Disease/context terms |",
        "|---|---|---|---:|---|",
    ]
    for row in article_rows:
        report_lines.append(
            f"| {row['pmid']} | {row['all_herb_titles_match_pubmed']} | {row['journal']} | {row['publication_year']} | {row['disease_context_terms_in_pubmed_record']} |"
        )
    report_lines.extend(
        [
            "",
            "## Interpretation Boundary",
            "",
            "- PubMed verification confirms bibliographic identity and title/abstract-level context, not full-text experimental validity.",
            "- Target-term detection is conservative and based on title/abstract text plus common cytokine/gene aliases.",
            "- Final manuscript tables should still manually inspect any high-impact claim before submission.",
            "",
        ]
    )
    DOC_PATH.write_text("\n".join(report_lines), encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
