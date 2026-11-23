#!/usr/bin/env python3
"""Merge GTDB genome metadata with BioSample fields and aggregated SRA runs.

One output row per GTDB genome. Joins are by BioSample accession.
No network requests; supports .tsv and .tsv.gz inputs.

Example:
  python scripts/04_merge_metadata.py \
    --gtdb data/processed/Serratia_gtdb.tsv \
    --biosample data/processed/Serratia_biosample.tsv \
    --sra data/processed/Serratia_sra.tsv \
    --out data/processed/Serratia_merged_metadata.tsv \
    --limit 10
"""

from __future__ import annotations

import argparse
import csv
import gzip
import sys
from collections import defaultdict
from pathlib import Path


MISSING = {"", "NA", "N/A", "NULL", "NONE", "-", "NAN"}


def open_text(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8-sig", newline="")
    return path.open("r", encoding="utf-8-sig", newline="")


def clean(value: str | None) -> str:
    value = (value or "").strip()
    return "" if value.upper() in MISSING else value


def reader_for(path: Path):
    handle = open_text(path)
    reader = csv.DictReader(handle, delimiter="\t")
    if not reader.fieldnames:
        handle.close()
        raise ValueError(f"Empty or headerless TSV: {path}")
    reader.fieldnames = [f.strip() for f in reader.fieldnames]
    return handle, reader


def find_key(columns: list[str], preferred: str, path: Path) -> str:
    by_lower = {column.lower(): column for column in columns}
    key = by_lower.get(preferred.lower())
    if key is None:
        raise ValueError(
            f"{path}: BioSample key '{preferred}' missing. "
            f"Columns: {', '.join(columns)}"
        )
    return key


def read_gtdb(path: Path, key_name: str, limit: int | None):
    handle, reader = reader_for(path)
    try:
        columns = list(reader.fieldnames or [])
        key = find_key(columns, key_name, path)
        rows = []
        for row in reader:
            rows.append(row)
            if limit is not None and len(rows) >= limit:
                break
        return columns, key, rows
    finally:
        handle.close()


def read_biosample(path: Path, key_name: str, wanted: set[str]):
    handle, reader = reader_for(path)
    try:
        columns = list(reader.fieldnames or [])
        key = find_key(columns, key_name, path)
        indexed = {}
        for row in reader:
            accession = clean(row.get(key))
            if accession in wanted and accession not in indexed:
                indexed[accession] = row
        return columns, key, indexed
    finally:
        handle.close()


def read_sra_summary(path: Path, key_name: str, wanted: set[str]):
    handle, reader = reader_for(path)
    try:
        columns = list(reader.fieldnames or [])
        key = find_key(columns, key_name, path)
        col = {name.lower(): name for name in columns}
        summaries = defaultdict(
            lambda: {"runs": 0, "bases": 0, "spots": 0,
                     "platforms": set(), "strategies": set(), "layouts": set()}
        )
        for row in reader:
            accession = clean(row.get(key))
            if accession not in wanted:
                continue
            summary = summaries[accession]
            summary["runs"] += 1
            for source, target in (("bases", "bases"), ("spots", "spots")):
                value = clean(row.get(col.get(source, "")))
                if value:
                    try:
                        summary[target] += int(float(value.replace(",", "")))
                    except ValueError:
                        pass
            for source, target in (("platform", "platforms"),
                                   ("librarystrategy", "strategies"),
                                   ("librarylayout", "layouts")):
                value = clean(row.get(col.get(source, "")))
                if value:
                    summary[target].add(value)
        return summaries
    finally:
        handle.close()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gtdb", type=Path, required=True, help="GTDB genus metadata TSV")
    parser.add_argument("--biosample", type=Path, required=True, help="BioSample TSV")
    parser.add_argument("--sra", type=Path, required=True, help="SRA runs TSV")
    parser.add_argument("--out", type=Path, required=True, help="Merged genome-level TSV")
    parser.add_argument("--gtdb-key", default="ncbi_biosample")
    parser.add_argument("--biosample-key", default="BioSample")
    parser.add_argument("--sra-key", default="BioSample")
    parser.add_argument("--limit", type=int, help="Only first N GTDB rows (quick test)")
    parser.add_argument("--check-only", action="store_true",
                        help="Check TSV headers/keys only; do not merge data")
    args = parser.parse_args()

    if args.limit is not None and args.limit < 1:
        parser.error("--limit must be >= 1")
    for path in (args.gtdb, args.biosample, args.sra):
        if not path.is_file():
            parser.error(f"Input file not found: {path}")

    if args.check_only:
        for path, key_name in ((args.gtdb, args.gtdb_key),
                               (args.biosample, args.biosample_key),
                               (args.sra, args.sra_key)):
            handle, reader = reader_for(path)
            try:
                cols = list(reader.fieldnames or [])
                key = find_key(cols, key_name, path)
                print(f"OK: {path} ({len(cols)} columns; join key: {key})")
            finally:
                handle.close()
        print("Header check passed. No records were processed.")
        return

    gtdb_cols, gtdb_key, gtdb_rows = read_gtdb(args.gtdb, args.gtdb_key, args.limit)
    wanted = {clean(row.get(gtdb_key)) for row in gtdb_rows} - {""}
    biosample_cols, bs_key, biosamples = read_biosample(
        args.biosample, args.biosample_key, wanted
    )
    sra_summaries = read_sra_summary(args.sra, args.sra_key, wanted)

    sra_cols = ["sra_run_count", "sra_total_bases", "sra_total_spots",
                "sra_platforms", "sra_library_strategies", "sra_library_layouts"]
    out_cols = (["BioSample"] + ["gtdb_" + c for c in gtdb_cols]
                + ["biosample_" + c for c in biosample_cols if c != bs_key]
                + sra_cols)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=out_cols, delimiter="\t",
                                lineterminator="\n")
        writer.writeheader()
        for gtdb_row in gtdb_rows:
            accession = clean(gtdb_row.get(gtdb_key))
            bs_row = biosamples.get(accession, {})
            stats = sra_summaries.get(accession)
            out_row = {"BioSample": accession}
            out_row.update({"gtdb_" + c: clean(gtdb_row.get(c)) for c in gtdb_cols})
            out_row.update({"biosample_" + c: clean(bs_row.get(c))
                            for c in biosample_cols if c != bs_key})
            if stats is not None:
                out_row.update({
                    "sra_run_count": stats["runs"],
                    "sra_total_bases": stats["bases"],
                    "sra_total_spots": stats["spots"],
                    "sra_platforms": ";".join(sorted(stats["platforms"])),
                    "sra_library_strategies": ";".join(sorted(stats["strategies"])),
                    "sra_library_layouts": ";".join(sorted(stats["layouts"])),
                })
            else:
                out_row["sra_run_count"] = 0
            writer.writerow(out_row)

    print(f"GTDB genomes: {len(gtdb_rows)}")
    print(f"GTDB BioSample IDs: {len(wanted)}")
    print(f"Matched BioSample records: {len(biosamples)}")
    print(f"BioSample IDs with SRA runs: {len(sra_summaries)}")
    print(f"Output: {args.out}")


if __name__ == "__main__":
    try:
        main()
    except (OSError, ValueError, csv.Error) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        sys.exit(1)
