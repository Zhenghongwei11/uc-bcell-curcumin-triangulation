#!/usr/bin/env python3
"""Download public raw inputs that are not acquired by downstream modules.

The downstream analysis scripts download several resources themselves
(HERB 2.0 full tables, bulk GEO matrices, Open Targets, ETCM2, PubMed, and
L1000CDS2). This helper covers public supplementary and metadata files that
must exist before the single-cell and LINCS modules run.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import io
import tarfile
import urllib.parse
import urllib.request
from pathlib import Path
from typing import Iterable


ACCESS_DATE = "2026-07-03"
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

DOWNLOADS = [
    (
        "geo_gse125527_raw_tar",
        "NCBI GEO",
        "GSE125527",
        "data/raw/geo/GSE125527/GSE125527_RAW.tar",
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE125nnn/GSE125527/suppl/GSE125527_RAW.tar",
        "GSE125527 supplementary raw archive; rectal UMI tables are extracted to data/raw/geo/GSE125527/raw_rectal.",
    ),
    (
        "geo_gse125527_gene_id_rownames",
        "NCBI GEO",
        "GSE125527",
        "data/raw/geo/GSE125527/GSE125527_gene_id_rownames.csv.gz",
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE125nnn/GSE125527/suppl/GSE125527_gene_id_rownames.csv.gz",
        "GSE125527 processed feature names.",
    ),
    (
        "geo_gse125527_cell_metadata",
        "NCBI GEO",
        "GSE125527",
        "data/raw/geo/GSE125527/GSE125527_cell_metadata.csv.gz",
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE125nnn/GSE125527/suppl/GSE125527_cell_metadata.csv.gz",
        "GSE125527 processed cell metadata.",
    ),
    (
        "geo_gse125527_cell_id_colnames",
        "NCBI GEO",
        "GSE125527",
        "data/raw/geo/GSE125527/GSE125527_cell_id_colnames.csv.gz",
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE125nnn/GSE125527/suppl/GSE125527_cell_id_colnames.csv.gz",
        "GSE125527 processed cell identifiers.",
    ),
    (
        "geo_gse125527_patient_id_map",
        "NCBI GEO",
        "GSE125527",
        "data/raw/geo/GSE125527/GSE125527_oldPatientId-newPatientId.csv.gz",
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE125nnn/GSE125527/suppl/GSE125527_oldPatientId-newPatientId.csv.gz",
        "GSE125527 old-to-new patient ID map.",
    ),
    (
        "geo_gse125527_bcell_cluster",
        "NCBI GEO",
        "GSE125527",
        "data/raw/geo/GSE125527/GSE125527_Bcell_cluster.csv.gz",
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE125nnn/GSE125527/suppl/GSE125527_Bcell_cluster.csv.gz",
        "GSE125527 B-cell cluster annotations.",
    ),
    (
        "geo_gse125527_tcell_cluster",
        "NCBI GEO",
        "GSE125527",
        "data/raw/geo/GSE125527/GSE125527_Tcell_cluster.csv.gz",
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE125nnn/GSE125527/suppl/GSE125527_Tcell_cluster.csv.gz",
        "GSE125527 T-cell cluster annotations.",
    ),
    (
        "geo_gse182270_raw_tar",
        "NCBI GEO",
        "GSE182270",
        "data/raw/geo/GSE182270/GSE182270_RAW.tar",
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE182nnn/GSE182270/suppl/GSE182270_RAW.tar",
        "GSE182270 supplementary 10x filtered matrices used for B-lineage replication.",
    ),
    (
        "lincs_gse92742_pert_info",
        "LINCS L1000 GSE92742",
        "2017-03-03",
        "data/raw/lincs/GSE92742_Broad_LINCS_pert_info.txt.gz",
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE92nnn/GSE92742/suppl/GSE92742_Broad_LINCS_pert_info.txt.gz",
        "Perturbagen metadata for exact InChIKey/PubChem matching.",
    ),
    (
        "lincs_gse70138_pert_info",
        "LINCS L1000 GSE70138",
        "2017-03-22",
        "data/raw/lincs/GSE70138_Broad_LINCS_pert_info.txt.gz",
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE70nnn/GSE70138/suppl/GSE70138_Broad_LINCS_pert_info.txt.gz",
        "Perturbagen metadata for exact InChIKey matching.",
    ),
]

PMC_DOWNLOADS = [
    (
        "pubmed_pmc_curcumin_33597887",
        "data/raw/pubmed/curcumin_pmc_PMC7882737.xml",
        "7882737",
        "Open PMC full-text XML for PMID 33597887 / PMC7882737.",
    ),
    (
        "pubmed_pmc_curcumin_36353208",
        "data/raw/pubmed/curcumin_pmc_PMC9639655.xml",
        "9639655",
        "Open PMC full-text XML for PMID 36353208 / PMC9639655.",
    ),
]


def md5_file(path: Path) -> str:
    digest = hashlib.md5()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def download(url: str, destination: Path, timeout: int, reuse_existing: bool) -> None:
    if reuse_existing and destination.exists() and destination.stat().st_size > 0:
        return
    destination.parent.mkdir(parents=True, exist_ok=True)
    tmp = destination.with_suffix(destination.suffix + ".tmp")
    request = urllib.request.Request(url, headers={"User-Agent": "uc-bcell-curcumin-repro/1.0"})
    with urllib.request.urlopen(request, timeout=timeout) as response:
        tmp.write_bytes(response.read())
    if destination.suffix == ".gz":
        with gzip.open(tmp, "rb") as handle:
            handle.read(1)
    tmp.replace(destination)


def efetch_pmc(pmc_numeric_id: str, destination: Path, timeout: int, reuse_existing: bool) -> None:
    params = urllib.parse.urlencode({"db": "pmc", "id": pmc_numeric_id, "retmode": "xml"})
    url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?{params}"
    download(url, destination, timeout, reuse_existing)


def extract_gse125527_rectal_tables(raw_tar: Path, out_dir: Path, reuse_existing: bool) -> int:
    existing = sorted(out_dir.glob("*_R_cell-gene_UMI_table.tsv.gz"))
    if reuse_existing and len(existing) >= 15:
        return len(existing)
    out_dir.mkdir(parents=True, exist_ok=True)
    n = 0
    with tarfile.open(raw_tar, "r:") as archive:
        for member in archive.getmembers():
            name = Path(member.name).name
            if not name.endswith("_R_cell-gene_UMI_table.tsv.gz"):
                continue
            extracted = archive.extractfile(member)
            if extracted is None:
                continue
            (out_dir / name).write_bytes(extracted.read())
            n += 1
    if n == 0:
        raise RuntimeError(f"No rectal UMI tables were extracted from {raw_tar}")
    return n


def read_manifest(path: Path) -> dict[str, dict[str, str]]:
    if not path.exists():
        return {}
    with path.open(encoding="utf-8", newline="") as handle:
        return {row["dataset_id"]: row for row in csv.DictReader(handle, delimiter="\t") if row.get("dataset_id")}


def write_manifest(path: Path, rows: Iterable[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    by_id = read_manifest(path)
    for row in rows:
        by_id[row["dataset_id"]] = row
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=MANIFEST_FIELDS, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for dataset_id in sorted(by_id):
            writer.writerow(by_id[dataset_id])


def manifest_row(dataset_id: str, source: str, version: str, path: Path, url: str, notes: str) -> dict[str, str]:
    return {
        "dataset_id": dataset_id,
        "source": source,
        "version_or_date": version,
        "file_path": str(path),
        "source_url_or_api": url,
        "access_date": ACCESS_DATE,
        "file_size_bytes": str(path.stat().st_size),
        "md5": md5_file(path),
        "status": "analysis_ready",
        "notes": notes,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--timeout", type=int, default=1200)
    parser.add_argument("--manifest", default="data/manifest.tsv")
    parser.add_argument("--reuse-existing", action="store_true")
    args = parser.parse_args()

    rows: list[dict[str, str]] = []
    for dataset_id, source, version, relpath, url, notes in DOWNLOADS:
        path = Path(relpath)
        print(f"{dataset_id}: {'reuse' if args.reuse_existing and path.exists() and path.stat().st_size > 0 else 'download'} {path}", flush=True)
        download(url, path, args.timeout, args.reuse_existing)
        rows.append(manifest_row(dataset_id, source, version, path, url, notes))

    for dataset_id, relpath, pmc_id, notes in PMC_DOWNLOADS:
        path = Path(relpath)
        print(f"{dataset_id}: {'reuse' if args.reuse_existing and path.exists() and path.stat().st_size > 0 else 'download'} {path}", flush=True)
        efetch_pmc(pmc_id, path, args.timeout, args.reuse_existing)
        rows.append(
            manifest_row(
                dataset_id,
                "NCBI PMC efetch",
                "2026-07-03",
                path,
                "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi",
                notes,
            )
        )

    n_rectal = extract_gse125527_rectal_tables(
        Path("data/raw/geo/GSE125527/GSE125527_RAW.tar"),
        Path("data/raw/geo/GSE125527/raw_rectal"),
        args.reuse_existing,
    )
    write_manifest(Path(args.manifest), rows)
    print(f"Downloaded public inputs and extracted {n_rectal} GSE125527 rectal UMI tables.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
