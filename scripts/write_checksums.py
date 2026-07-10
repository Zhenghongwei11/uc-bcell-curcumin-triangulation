#!/usr/bin/env python3
from pathlib import Path
import hashlib

root = Path(__file__).resolve().parents[1]
patterns = [
    "data/derived/*.tsv",
    "tables/main/*.tsv",
    "tables/supplementary/*.tsv",
    "figures/source_data/chinese_medicine/*.tsv",
    "plots/*.pdf",
    "plots/*.png",
    "docs/*.tsv",
    "docs/*.md",
]
rows = []
for pattern in patterns:
    for path in sorted(root.glob(pattern)):
        h = hashlib.sha256(path.read_bytes()).hexdigest()
        rows.append(f"{h}  {path.relative_to(root).as_posix()}")
(root / "CHECKSUMS.sha256").write_text("\n".join(rows) + "\n", encoding="utf-8")
print(root / "CHECKSUMS.sha256")
