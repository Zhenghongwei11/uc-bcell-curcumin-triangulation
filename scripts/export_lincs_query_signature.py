#!/usr/bin/env python3
"""Export top UC/IBD consensus genes as LINCS/CMap query gene sets."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path
from typing import Dict, Iterable, List


def read_tsv(path: Path) -> Iterable[Dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        yield from csv.DictReader(handle, delimiter="\t")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--consensus", default="results/m4_geo_uc_ibd_signature/consensus_gene_signature.tsv")
    parser.add_argument("--out-dir", default="results/m4_geo_uc_ibd_signature")
    parser.add_argument("--top-n", type=int, default=150)
    parser.add_argument("--max-fdr", type=float, default=0.05)
    args = parser.parse_args()

    rows = [
        row
        for row in read_tsv(Path(args.consensus))
        if float(row["best_fdr"]) <= args.max_fdr and row["gene_symbol"] and not row["gene_symbol"].startswith("LOC")
    ]
    up = [row for row in rows if row["direction"] == "up"][: args.top_n]
    down = [row for row in rows if row["direction"] == "down"][: args.top_n]

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / "lincs_query_up_genes.txt").write_text("\n".join(row["gene_symbol"] for row in up) + "\n", encoding="utf-8")
    (out_dir / "lincs_query_down_genes.txt").write_text("\n".join(row["gene_symbol"] for row in down) + "\n", encoding="utf-8")

    fields = [
        "query_direction",
        "gene_symbol",
        "gene_title",
        "evidence_groups",
        "mean_log2fc",
        "best_fdr",
        "consensus_score",
    ]
    with (out_dir / "lincs_query_signature.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for query_direction, selected in [("up", up), ("down", down)]:
            for row in selected:
                writer.writerow({"query_direction": query_direction, **row})

    summary = [
        {"metric": "top_n_requested", "value": str(args.top_n), "notes": ""},
        {"metric": "max_fdr", "value": str(args.max_fdr), "notes": ""},
        {"metric": "up_genes_exported", "value": str(len(up)), "notes": str(out_dir / "lincs_query_up_genes.txt")},
        {"metric": "down_genes_exported", "value": str(len(down)), "notes": str(out_dir / "lincs_query_down_genes.txt")},
    ]
    with (out_dir / "lincs_query_signature_summary.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["metric", "value", "notes"], delimiter="\t")
        writer.writeheader()
        writer.writerows(summary)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
