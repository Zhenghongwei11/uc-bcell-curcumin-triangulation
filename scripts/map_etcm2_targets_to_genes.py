#!/usr/bin/env python3
"""Map ETCM2 target names to human gene symbols and recompute overlaps."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
import time
import urllib.parse
import urllib.request
import urllib.error
from collections import defaultdict
from pathlib import Path
from typing import Any, Dict, Iterable, List, Sequence, Set, Tuple


ACCESS_DATE = "2026-06-29"
MYGENE_URL = "https://mygene.info/v3/query"
UNIPROT_URL = "https://rest.uniprot.org/uniprotkb/search"
RAW_DIR = Path("data/raw/target_mapping/etcm2")
OUT_DIR = Path("results/m15_etcm2_target_mapping")
DOC_PATH = Path("docs/M15_ETCM2_TARGET_MAPPING_RUN.md")
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

MANUAL_ALIASES = {
    "Aromatase cytochrome P450": ("CYP19A1", "aromatase", "P11511"),
    "beta-Secretase (BACE-1)": ("BACE1", "beta-secretase 1", "P56817"),
    "Corticotropin-releasing factor-binding protein": ("CRHBP", "corticotropin-releasing hormone-binding protein", "P24387"),
    "Estradiol 17-beta-dehydrogenase 1": ("HSD17B1", "17-beta-hydroxysteroid dehydrogenase type 1", "P14061"),
    "Cytochrome P450 1A": ("", "", ""),
    "Cholinesterase": ("", "", ""),
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


def clean_id(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", value).strip("_") or "unknown"


def normalize(value: str) -> str:
    return re.sub(r"[^a-z0-9]+", "", value.lower())


def parenthetical_alias(value: str) -> str:
    match = re.search(r"\(([A-Za-z0-9 -]+)\)", value)
    if not match:
        return ""
    return re.sub(r"[^A-Za-z0-9]+", "", match.group(1)).upper()


def url_json(url: str, params: Dict[str, str], timeout: int) -> Dict[str, Any]:
    full_url = f"{url}?{urllib.parse.urlencode(params)}"
    request = urllib.request.Request(full_url, headers={"User-Agent": "zyy-target-name-mapping/1.0"})
    with urllib.request.urlopen(request, timeout=timeout) as response:
        return json.loads(response.read().decode("utf-8", errors="replace"))


def safe_url_json(url: str, params: Dict[str, str], timeout: int) -> Dict[str, Any]:
    try:
        return url_json(url, params, timeout)
    except (urllib.error.HTTPError, urllib.error.URLError, TimeoutError, json.JSONDecodeError, OSError) as exc:
        return {"_error": f"{type(exc).__name__}: {exc}", "_params": params}


def load_or_fetch(path: Path, fetcher, reuse_existing: bool, sleep_seconds: float) -> Dict[str, Any]:
    path.parent.mkdir(parents=True, exist_ok=True)
    if reuse_existing and path.exists() and path.stat().st_size > 10:
        return json.loads(path.read_text(encoding="utf-8"))
    payload = fetcher()
    path.write_text(json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8")
    if sleep_seconds:
        time.sleep(sleep_seconds)
    return payload


def sanitize_mygene_query(value: str) -> str:
    cleaned = re.sub(r"[\[\]{}]", " ", value)
    cleaned = re.sub(r"\s+", " ", cleaned).strip()
    return cleaned


def uniprot_gene(result: Dict[str, Any]) -> str:
    genes = result.get("genes") or []
    if not genes:
        return ""
    return ((genes[0].get("geneName") or {}).get("value") or "").upper()


def uniprot_protein_name(result: Dict[str, Any]) -> str:
    desc = result.get("proteinDescription") or {}
    return (((desc.get("recommendedName") or {}).get("fullName") or {}).get("value") or "")


def mygene_ensembl(hit: Dict[str, Any]) -> str:
    ensembl = hit.get("ensembl")
    if isinstance(ensembl, dict):
        return str(ensembl.get("gene", ""))
    if isinstance(ensembl, list) and ensembl:
        return ";".join(sorted({str(item.get("gene", "")) for item in ensembl if item.get("gene")}))
    return ""


def mygene_uniprot(hit: Dict[str, Any]) -> str:
    uniprot = hit.get("uniprot") or {}
    swiss = uniprot.get("Swiss-Prot")
    if isinstance(swiss, list):
        return ";".join(str(value) for value in swiss)
    return str(swiss or "")


def map_one_target(name: str, timeout: int, reuse_existing: bool, sleep_seconds: float) -> Dict[str, Any]:
    target_key = clean_id(name)
    alias_symbol, alias_name, alias_uniprot = MANUAL_ALIASES.get(name, ("", "", ""))

    uniprot_exact = load_or_fetch(
        RAW_DIR / "uniprot_exact" / f"{target_key}.json",
        lambda: url_json(
            UNIPROT_URL,
            {
                "query": f'protein_name:"{name}" AND organism_id:9606 AND reviewed:true',
                "format": "json",
                "fields": "accession,gene_names,protein_name,organism_name",
                "size": "10",
            },
            timeout,
        ),
        reuse_existing,
        sleep_seconds,
    )
    mygene = load_or_fetch(
        RAW_DIR / "mygene" / f"{target_key}.json",
        lambda: safe_url_json(
            MYGENE_URL,
            {
                "q": sanitize_mygene_query(name),
                "species": "human",
                "fields": "symbol,name,entrezgene,ensembl.gene,uniprot.Swiss-Prot",
                "size": "10",
            },
            timeout,
        ),
        reuse_existing,
        sleep_seconds,
    )

    candidates: List[Dict[str, Any]] = []
    for result in uniprot_exact.get("results", []):
        candidates.append(
            {
                "source": "UniProt_exact_protein_name",
                "symbol": uniprot_gene(result),
                "name": uniprot_protein_name(result),
                "uniprot": result.get("primaryAccession", ""),
                "score": "",
            }
        )
    for hit in mygene.get("hits", []):
        candidates.append(
            {
                "source": "MyGene_query",
                "symbol": str(hit.get("symbol", "")).upper(),
                "name": str(hit.get("name", "")),
                "uniprot": mygene_uniprot(hit),
                "entrezgene": hit.get("entrezgene", ""),
                "ensembl_gene": mygene_ensembl(hit),
                "score": hit.get("_score", ""),
            }
        )

    n_name = normalize(name)
    alias = parenthetical_alias(name)
    best: Dict[str, Any] = {}
    confidence = "unmapped"
    rationale = "No confident human gene mapping from UniProt/MyGene"

    if alias_symbol:
        best = {"source": "manual_alias", "symbol": alias_symbol, "name": alias_name, "uniprot": alias_uniprot, "score": ""}
        confidence = "high" if alias_symbol else "ambiguous"
        rationale = "Manual alias override for known target name"
    elif name in MANUAL_ALIASES and not alias_symbol:
        confidence = "ambiguous"
        rationale = "Manual ambiguity override; target name maps to multiple possible genes"
    else:
        exact_uniprot = [
            candidate
            for candidate in candidates
            if candidate["source"] == "UniProt_exact_protein_name"
            and normalize(candidate.get("name", "")) == n_name
            and candidate.get("symbol")
        ]
        if exact_uniprot:
            best = exact_uniprot[0]
            confidence = "exact"
            rationale = "UniProt reviewed human protein-name exact match"
        else:
            exact_mygene = [
                candidate
                for candidate in candidates
                if candidate["source"] == "MyGene_query"
                and normalize(candidate.get("name", "")) == n_name
                and candidate.get("symbol")
            ]
            if exact_mygene:
                best = exact_mygene[0]
                confidence = "exact"
                rationale = "MyGene human gene-name exact match"
            elif alias:
                alias_hits = [candidate for candidate in candidates if normalize(candidate.get("symbol", "")) == normalize(alias)]
                if alias_hits:
                    best = alias_hits[0]
                    confidence = "high"
                    rationale = "Parenthetical target alias matches mapped gene symbol"
            if not best:
                mygene_hits = [candidate for candidate in candidates if candidate["source"] == "MyGene_query" and candidate.get("symbol")]
                if mygene_hits:
                    top = mygene_hits[0]
                    top_name = normalize(top.get("name", ""))
                    if n_name and (n_name in top_name or top_name in n_name):
                        best = top
                        confidence = "medium"
                        rationale = "Top MyGene hit has partial name containment"
                    else:
                        best = top
                        confidence = "low"
                        rationale = "Top MyGene hit selected but name match is weak"

    candidate_preview = " || ".join(
        f"{candidate.get('source')}:{candidate.get('symbol')}:{candidate.get('name')}:{candidate.get('uniprot')}:{candidate.get('score')}"
        for candidate in candidates[:5]
    )
    return {
        "etcm_target_name": name,
        "mapped_gene_symbol": best.get("symbol", ""),
        "mapped_gene_name": best.get("name", ""),
        "mapped_entrezgene": best.get("entrezgene", ""),
        "mapped_ensembl_gene": best.get("ensembl_gene", ""),
        "mapped_uniprot": best.get("uniprot", ""),
        "mapping_confidence": confidence,
        "mapping_source": best.get("source", ""),
        "mapping_rationale": rationale,
        "candidate_preview": candidate_preview,
    }


def update_manifest(paths: Sequence[Path]) -> None:
    existing: List[Dict[str, str]] = []
    if MANIFEST_PATH.exists():
        with MANIFEST_PATH.open(encoding="utf-8", newline="") as handle:
            existing = list(csv.DictReader(handle, delimiter="\t"))
    by_id = {row["dataset_id"]: row for row in existing if row.get("dataset_id")}
    for path in paths:
        by_id[f"etcm2_target_mapping_cache_{path.parent.name}_{path.stem.lower()}"] = {
            "dataset_id": f"etcm2_target_mapping_cache_{path.parent.name}_{path.stem.lower()}",
            "source": "UniProt/MyGene target mapping API",
            "version_or_date": ACCESS_DATE,
            "file_path": str(path),
            "source_url_or_api": f"{UNIPROT_URL}; {MYGENE_URL}",
            "access_date": ACCESS_DATE,
            "file_size_bytes": str(path.stat().st_size),
            "md5": md5_file(path),
            "status": "analysis_ready",
            "notes": "Raw API cache for ETCM2 target-name to human gene-symbol mapping",
        }
    with MANIFEST_PATH.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=MANIFEST_FIELDS, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for key in sorted(by_id):
            writer.writerow(by_id[key])


def build_mapped_edges(edge_rows: Sequence[Dict[str, str]], mapping_rows: Sequence[Dict[str, Any]]) -> List[Dict[str, Any]]:
    by_name = {row["etcm_target_name"]: row for row in mapping_rows}
    rows = []
    for edge in edge_rows:
        mapping = by_name.get(edge["etcm_target_name"], {})
        rows.append({**edge, **mapping})
    return rows


def load_reference_sets(m9_edges: Path, m10_summary: Path) -> Tuple[Dict[str, Set[str]], Dict[str, Set[str]], Dict[str, Set[str]]]:
    herb_m9: Dict[str, Set[str]] = defaultdict(set)
    bulk_aligned: Dict[str, Set[str]] = defaultdict(set)
    for row in read_tsv(m9_edges):
        herb_m9[row["ingredient_id"]].add(row["gene_symbol"].upper())
        if row.get("therapeutic_alignment", "").startswith("aligned"):
            bulk_aligned[row["ingredient_id"]].add(row["gene_symbol"].upper())
    genetic: Dict[str, Set[str]] = defaultdict(set)
    for row in read_tsv(m10_summary):
        genetic[row["ingredient_id"]].update(gene for gene in row.get("genetic_supported_genes", "").split(";") if gene)
    return herb_m9, bulk_aligned, genetic


def summarize(mapped_edges: Sequence[Dict[str, Any]], herb_m9: Dict[str, Set[str]], bulk_aligned: Dict[str, Set[str]], genetic: Dict[str, Set[str]]) -> List[Dict[str, Any]]:
    by_candidate: Dict[Tuple[str, str], List[Dict[str, Any]]] = defaultdict(list)
    for row in mapped_edges:
        by_candidate[(row["ingredient_id"], row["ingredient_name"])].append(row)
    rows = []
    accepted_conf = {"exact", "high"}
    for (ingredient_id, ingredient_name), edges in by_candidate.items():
        accepted_genes = {
            row["mapped_gene_symbol"]
            for row in edges
            if row.get("mapped_gene_symbol") and row.get("mapping_confidence") in accepted_conf
        }
        all_genes = {row["mapped_gene_symbol"] for row in edges if row.get("mapped_gene_symbol")}
        rows.append(
            {
                "ingredient_id": ingredient_id,
                "ingredient_name": ingredient_name,
                "n_etcm2_target_name_edges": len(edges),
                "n_mapped_gene_symbols_all_confidence": len(all_genes),
                "mapped_gene_symbols_all_confidence": ";".join(sorted(all_genes)),
                "n_accepted_mapped_gene_symbols_exact_high": len(accepted_genes),
                "accepted_mapped_gene_symbols_exact_high": ";".join(sorted(accepted_genes)),
                "n_overlap_with_herb_m9_disease_context_targets": len(accepted_genes & herb_m9.get(ingredient_id, set())),
                "overlap_with_herb_m9_disease_context_targets": ";".join(sorted(accepted_genes & herb_m9.get(ingredient_id, set()))),
                "n_overlap_with_bulk_aligned_targets": len(accepted_genes & bulk_aligned.get(ingredient_id, set())),
                "overlap_with_bulk_aligned_targets": ";".join(sorted(accepted_genes & bulk_aligned.get(ingredient_id, set()))),
                "n_overlap_with_open_targets_genetic_supported": len(accepted_genes & genetic.get(ingredient_id, set())),
                "overlap_with_open_targets_genetic_supported": ";".join(sorted(accepted_genes & genetic.get(ingredient_id, set()))),
            }
        )
    rows.sort(
        key=lambda row: (
            -int(row["n_overlap_with_herb_m9_disease_context_targets"]),
            -int(row["n_accepted_mapped_gene_symbols_exact_high"]),
            row["ingredient_name"].lower(),
        )
    )
    return rows


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--etcm-edges", default="results/m11_etcm2_cross_validation/etcm2_ingredient_target_edges.tsv")
    parser.add_argument("--m9-edges", default="results/m9_target_cell_bridge/disease_context_herb_target_edges.tsv")
    parser.add_argument("--m10-summary", default="results/m10_open_targets/candidate_open_targets_summary.tsv")
    parser.add_argument("--timeout", type=int, default=30)
    parser.add_argument("--sleep-seconds", type=float, default=0.2)
    parser.add_argument("--no-reuse-existing", action="store_true")
    args = parser.parse_args()

    edge_rows = read_tsv(Path(args.etcm_edges))
    target_names = sorted({row["etcm_target_name"] for row in edge_rows if row.get("etcm_target_name")})
    mapping_rows = [
        map_one_target(name, args.timeout, not args.no_reuse_existing, args.sleep_seconds)
        for name in target_names
    ]
    mapped_edges = build_mapped_edges(edge_rows, mapping_rows)
    herb_m9, bulk_aligned, genetic = load_reference_sets(Path(args.m9_edges), Path(args.m10_summary))
    summary_rows = summarize(mapped_edges, herb_m9, bulk_aligned, genetic)

    mapping_fields = [
        "etcm_target_name",
        "mapped_gene_symbol",
        "mapped_gene_name",
        "mapped_entrezgene",
        "mapped_ensembl_gene",
        "mapped_uniprot",
        "mapping_confidence",
        "mapping_source",
        "mapping_rationale",
        "candidate_preview",
    ]
    mapped_edge_fields = list(edge_rows[0].keys()) + [
        "mapped_gene_symbol",
        "mapped_gene_name",
        "mapped_entrezgene",
        "mapped_ensembl_gene",
        "mapped_uniprot",
        "mapping_confidence",
        "mapping_source",
        "mapping_rationale",
        "candidate_preview",
    ]
    summary_fields = [
        "ingredient_id",
        "ingredient_name",
        "n_etcm2_target_name_edges",
        "n_mapped_gene_symbols_all_confidence",
        "mapped_gene_symbols_all_confidence",
        "n_accepted_mapped_gene_symbols_exact_high",
        "accepted_mapped_gene_symbols_exact_high",
        "n_overlap_with_herb_m9_disease_context_targets",
        "overlap_with_herb_m9_disease_context_targets",
        "n_overlap_with_bulk_aligned_targets",
        "overlap_with_bulk_aligned_targets",
        "n_overlap_with_open_targets_genetic_supported",
        "overlap_with_open_targets_genetic_supported",
    ]
    write_tsv(OUT_DIR / "etcm2_target_name_gene_mapping.tsv", mapping_rows, mapping_fields)
    write_tsv(OUT_DIR / "etcm2_target_edges_mapped.tsv", mapped_edges, mapped_edge_fields)
    write_tsv(OUT_DIR / "candidate_etcm2_mapped_overlap_summary.tsv", summary_rows, summary_fields)
    update_manifest(sorted(RAW_DIR.glob("*/*.json")))

    counts = defaultdict(int)
    for row in mapping_rows:
        counts[row["mapping_confidence"]] += 1
    report_lines = [
        "# M15 ETCM2 Target-Name To Gene-Symbol Mapping Run",
        "",
        "## Scope",
        "",
        "This run maps ETCM2 ingredient-detail target names to human gene symbols using UniProt reviewed protein-name search and MyGene human gene queries. Original ETCM2 target names are preserved; only exact/high-confidence mappings are used for gene-level overlap.",
        "",
        "## Key Outputs",
        "",
        "- `results/m15_etcm2_target_mapping/etcm2_target_name_gene_mapping.tsv`",
        "- `results/m15_etcm2_target_mapping/etcm2_target_edges_mapped.tsv`",
        "- `results/m15_etcm2_target_mapping/candidate_etcm2_mapped_overlap_summary.tsv`",
        "",
        "## Mapping Summary",
        "",
        f"- Unique ETCM2 target names: {len(target_names)}.",
        f"- Exact mappings: {counts['exact']}.",
        f"- High-confidence mappings: {counts['high']}.",
        f"- Medium-confidence mappings: {counts['medium']}.",
        f"- Low-confidence mappings: {counts['low']}.",
        f"- Ambiguous mappings: {counts['ambiguous']}.",
        f"- Unmapped target names: {counts['unmapped']}.",
        "",
        "## Candidate Gene-Level Overlap After Mapping",
        "",
        "| Ingredient | Accepted mapped genes | HERB/M9 overlap | Bulk-aligned overlap | Open Targets genetic overlap |",
        "|---|---:|---:|---:|---:|",
    ]
    for row in summary_rows:
        report_lines.append(
            f"| {row['ingredient_name']} | {row['n_accepted_mapped_gene_symbols_exact_high']} | {row['n_overlap_with_herb_m9_disease_context_targets']} | {row['n_overlap_with_bulk_aligned_targets']} | {row['n_overlap_with_open_targets_genetic_supported']} |"
        )
    report_lines.extend(
        [
            "",
            "## Interpretation Boundary",
            "",
            "- Gene-level overlap is computed only from `exact` and `high` mappings.",
            "- `medium`, `low`, `ambiguous`, and `unmapped` target-name mappings should remain in supplementary audit tables, not primary mechanistic claims.",
            "- Mapping target names to genes can introduce ambiguity for family-level names such as `Cytochrome P450 1A` or `Cholinesterase`; these are explicitly downgraded.",
            "",
        ]
    )
    DOC_PATH.write_text("\n".join(report_lines), encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
