#!/usr/bin/env python3
"""Lightweight GEO metadata audit for the UC/IBD TCM public-data project.

The script reads only GEO series-matrix metadata headers and supplementary
directory listings. It stops before expression tables, so it is suitable as a
pre-analysis feasibility check.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import io
import re
import sys
import time
import urllib.error
import urllib.request
from pathlib import Path
from typing import Dict, Iterable, List, Tuple


DEFAULT_ACCESSIONS = [
    "GSE75214",
    "GSE11223",
    "GSE59071",
    "GSE87466",
    "GSE92415",
    "GSE16879",
    "GSE73661",
    "GSE125527",
    "GSE134809",
    "GSE282580",
    "GSE270661",
]

SAMPLE_TAGS = [
    "!Sample_geo_accession",
    "!Sample_title",
    "!Sample_source_name_ch1",
    "!Sample_characteristics_ch1",
    "!Sample_treatment_protocol_ch1",
    "!Sample_description",
]

KEYWORD_PATTERNS = {
    "disease": re.compile(r"\b(ulcerative colitis|crohn'?s?|ibd|inflammatory bowel|colitis|control|normal|healthy|uc|cd)\b", re.I),
    "tissue": re.compile(r"\b(colon|colonic|ileum|ileal|rectum|rectal|sigmoid|mucosa|mucosal|biopsy|blood|pbmc|organoid)\b", re.I),
    "inflammation": re.compile(r"\b(inflamed|uninflamed|non[- ]?inflamed|active|inactive|normal|lesion|non[- ]?lesion|inflammatory)\b", re.I),
    "treatment": re.compile(r"\b(infliximab|vedolizumab|golimumab|anti[- ]?tnf|tnf|placebo|baseline|week|treatment|treated|therapy|induction)\b", re.I),
    "response": re.compile(r"\b(responder|non[- ]?responder|response|remission|refractory|resistant|effective|ineffective|mayo|endoscopic)\b", re.I),
}


def gse_prefix(accession: str) -> str:
    match = re.fullmatch(r"(GSE)(\d+)", accession)
    if not match:
        raise ValueError(f"Not a GSE accession: {accession}")
    number = match.group(2)
    return f"GSE{number[:-3]}nnn" if len(number) > 3 else "GSEnnn"


def geo_urls(accession: str) -> Dict[str, str]:
    prefix = gse_prefix(accession)
    base = f"https://ftp.ncbi.nlm.nih.gov/geo/series/{prefix}/{accession}"
    return {
        "geo_page": f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={accession}",
        "matrix": f"{base}/matrix/{accession}_series_matrix.txt.gz",
        "supplementary_dir": f"{base}/suppl/",
    }


def fetch_bytes(url: str, timeout: int, retries: int = 2) -> bytes:
    last_error: Exception | None = None
    headers = {"User-Agent": "zyy-uc-ibd-metadata-audit/1.0"}
    for attempt in range(retries + 1):
        try:
            request = urllib.request.Request(url, headers=headers)
            with urllib.request.urlopen(request, timeout=timeout) as response:
                return response.read()
        except (urllib.error.URLError, TimeoutError) as exc:
            last_error = exc
            if attempt < retries:
                time.sleep(1.5 * (attempt + 1))
    assert last_error is not None
    raise last_error


def open_url(url: str, timeout: int):
    headers = {"User-Agent": "zyy-uc-ibd-metadata-audit/1.0"}
    request = urllib.request.Request(url, headers=headers)
    return urllib.request.urlopen(request, timeout=timeout)


def list_supplementary_files(url: str, timeout: int) -> Tuple[str, List[str]]:
    try:
        html = fetch_bytes(url, timeout=timeout).decode("utf-8", errors="replace")
    except urllib.error.HTTPError as exc:
        return f"HTTP_{exc.code}", []
    except Exception as exc:
        return f"ERROR_{type(exc).__name__}", []

    files = []
    for href in re.findall(r'href="([^"]+)"', html):
        if href in {"../", "/"} or href.endswith("/") or href.startswith("http://") or href.startswith("https://"):
            continue
        files.append(href)
    return "available" if files else "empty", files


def parse_tsv_line(line: str) -> List[str]:
    reader = csv.reader([line.rstrip("\n")], delimiter="\t")
    return next(reader)


def read_series_header(url: str, timeout: int) -> Tuple[str, Dict[str, List[List[str]]], int, str]:
    try:
        response = open_url(url, timeout=timeout)
    except urllib.error.HTTPError as exc:
        return f"HTTP_{exc.code}", {}, 0, ""
    except Exception as exc:
        return f"ERROR_{type(exc).__name__}", {}, 0, ""

    tags: Dict[str, List[List[str]]] = {}
    header_lines = 0
    compressed_bytes = response.headers.get("Content-Length", "")
    try:
        with response:
            with gzip.GzipFile(fileobj=response, mode="rb") as gz_handle:
                text_handle = io.TextIOWrapper(gz_handle, encoding="utf-8", errors="replace")
                for line in text_handle:
                    if line.startswith("!series_matrix_table_begin"):
                        break
                    header_lines += 1
                    if not line.startswith("!"):
                        continue
                    parts = parse_tsv_line(line)
                    tag = parts[0]
                    if tag.startswith("!Series_") or tag in SAMPLE_TAGS:
                        tags.setdefault(tag, []).append(parts[1:])
    except Exception as exc:
        return f"ERROR_{type(exc).__name__}", tags, header_lines, compressed_bytes

    return "available", tags, header_lines, compressed_bytes


def first_or_empty(values: List[str], index: int) -> str:
    if index < len(values):
        return values[index]
    return ""


def sample_rows(accession: str, tags: Dict[str, List[List[str]]]) -> List[Dict[str, str]]:
    geo = tags.get("!Sample_geo_accession", [[]])[0] if tags.get("!Sample_geo_accession") else []
    titles = tags.get("!Sample_title", [[]])[0] if tags.get("!Sample_title") else []
    sources = tags.get("!Sample_source_name_ch1", [[]])[0] if tags.get("!Sample_source_name_ch1") else []
    treatments = tags.get("!Sample_treatment_protocol_ch1", [[]])[0] if tags.get("!Sample_treatment_protocol_ch1") else []
    descriptions = tags.get("!Sample_description", [[]])[0] if tags.get("!Sample_description") else []
    characteristics = tags.get("!Sample_characteristics_ch1", [])

    sample_count = max(len(geo), len(titles), len(sources), len(treatments), len(descriptions), *(len(row) for row in characteristics)) if characteristics or geo or titles else 0
    rows: List[Dict[str, str]] = []
    for index in range(sample_count):
        char_values = [first_or_empty(row, index) for row in characteristics if first_or_empty(row, index)]
        title = first_or_empty(titles, index)
        source = first_or_empty(sources, index)
        treatment = first_or_empty(treatments, index)
        description = first_or_empty(descriptions, index)
        label_specific = " | ".join([title, source, *char_values])
        treatment_context = " | ".join([label_specific, treatment])
        response_context = " | ".join([label_specific, treatment, description])
        rows.append(
            {
                "accession": accession,
                "gsm": first_or_empty(geo, index),
                "sample_title": title,
                "source_name": source,
                "characteristics": " || ".join(char_values),
                "treatment_protocol": treatment,
                "description": description,
                "disease_terms": extract_terms("disease", label_specific),
                "tissue_terms": extract_terms("tissue", label_specific),
                "inflammation_terms": extract_terms("inflammation", label_specific),
                "treatment_terms": extract_terms("treatment", treatment_context),
                "response_terms": extract_terms("response", response_context),
                "is_control_like": "yes" if re.search(r"\b(control|normal|healthy)\b", label_specific, re.I) else "no",
            }
        )
    return rows


def extract_terms(kind: str, text: str) -> str:
    pattern = KEYWORD_PATTERNS[kind]
    terms = sorted({match.group(0).lower() for match in pattern.finditer(text)})
    return ";".join(terms)


def unique_nonempty(values: Iterable[str], limit: int = 30) -> str:
    seen: List[str] = []
    for value in values:
        if not value:
            continue
        for part in value.split(";"):
            if part and part not in seen:
                seen.append(part)
    return ";".join(seen[:limit])


def write_tsv(path: Path, rows: List[Dict[str, str]], fieldnames: List[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def audit(accessions: List[str], out_dir: Path, timeout: int) -> int:
    summary_rows: List[Dict[str, str]] = []
    all_sample_rows: List[Dict[str, str]] = []

    for accession in accessions:
        urls = geo_urls(accession)
        print(f"[geo-audit] {accession}", file=sys.stderr)
        matrix_status, tags, header_lines, compressed_bytes = read_series_header(urls["matrix"], timeout=timeout)
        suppl_status, suppl_files = list_supplementary_files(urls["supplementary_dir"], timeout=timeout)
        rows = sample_rows(accession, tags) if matrix_status == "available" else []
        all_sample_rows.extend(rows)

        series_title = " | ".join(tags.get("!Series_title", [[""]])[0][:3]) if tags.get("!Series_title") else ""
        series_summary = " | ".join(tags.get("!Series_summary", [[""]])[0][:2]) if tags.get("!Series_summary") else ""
        summary_rows.append(
            {
                "accession": accession,
                "geo_page": urls["geo_page"],
                "matrix_url": urls["matrix"],
                "matrix_status": matrix_status,
                "series_title": series_title,
                "series_summary": series_summary[:500],
                "header_lines_read": str(header_lines),
                "compressed_matrix_bytes": compressed_bytes,
                "parsed_sample_count": str(len(rows)),
                "supplementary_status": suppl_status,
                "supplementary_file_count": str(len(suppl_files)),
                "supplementary_files": ";".join(suppl_files[:25]),
                "disease_terms_observed": unique_nonempty(row["disease_terms"] for row in rows),
                "tissue_terms_observed": unique_nonempty(row["tissue_terms"] for row in rows),
                "inflammation_terms_observed": unique_nonempty(row["inflammation_terms"] for row in rows),
                "treatment_terms_observed": unique_nonempty(row["treatment_terms"] for row in rows),
                "response_terms_observed": unique_nonempty(row["response_terms"] for row in rows),
                "control_like_samples": str(sum(1 for row in rows if row["is_control_like"] == "yes")),
            }
        )

    write_tsv(
        out_dir / "geo_series_summary.tsv",
        summary_rows,
        [
            "accession",
            "geo_page",
            "matrix_url",
            "matrix_status",
            "series_title",
            "series_summary",
            "header_lines_read",
            "compressed_matrix_bytes",
            "parsed_sample_count",
            "supplementary_status",
            "supplementary_file_count",
            "supplementary_files",
            "disease_terms_observed",
            "tissue_terms_observed",
            "inflammation_terms_observed",
            "treatment_terms_observed",
            "response_terms_observed",
            "control_like_samples",
        ],
    )
    write_tsv(
        out_dir / "geo_sample_metadata.tsv",
        all_sample_rows,
        [
            "accession",
            "gsm",
            "sample_title",
            "source_name",
            "characteristics",
            "treatment_protocol",
            "description",
            "disease_terms",
            "tissue_terms",
            "inflammation_terms",
            "treatment_terms",
            "response_terms",
            "is_control_like",
        ],
    )
    return 0


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out-dir", default="results/m1_uc_ibd_metadata_audit", help="Output directory for audit TSV files.")
    parser.add_argument("--timeout", type=int, default=90, help="HTTP timeout per GEO request in seconds.")
    parser.add_argument("accessions", nargs="*", default=DEFAULT_ACCESSIONS, help="GEO GSE accessions to audit.")
    args = parser.parse_args()

    return audit(args.accessions, Path(args.out_dir), args.timeout)


if __name__ == "__main__":
    raise SystemExit(main())
