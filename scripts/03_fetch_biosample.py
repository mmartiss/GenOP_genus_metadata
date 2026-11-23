#!/usr/bin/env python3
"""Fetch and normalize NCBI BioSample metadata for BioSample IDs in an SRA TSV.

Requires Python >= 3.9; standard library only. Raw XML is cached in batches.
"""

import argparse
import csv
import gzip
import hashlib
import os
from collections import Counter, defaultdict
from pathlib import Path
import re
import sys
import time
from urllib.error import HTTPError, URLError
from urllib.parse import urlencode
from urllib.request import Request, urlopen
import xml.etree.ElementTree as ET

EFETCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
VALID_ID = re.compile(r"^SAM[NED]\d+$", re.I)

# Aliases are checked in order; the harmonized NCBI name is preferred.
FIELDS = {
    "sample_name": ("sample_name", "sample_alias"),
    "strain": ("strain",),
    "isolate": ("isolate", "isolate_name_alias"),
    "isolation_source": ("isolation_source", "isolation_site"),
    "geo_loc_name": ("geo_loc_name", "geographic_location", "country"),
    "lat_lon": ("lat_lon", "latitude_and_longitude", "latitude_longitude"),
    "collection_date": ("collection_date", "sampling_date"),
    "host": ("host", "host_scientific_name"),
    "host_disease": ("host_disease",),
    "env_broad_scale": ("env_broad_scale", "environmental_biome"),
    "env_local_scale": ("env_local_scale", "environmental_feature"),
    "env_medium": ("env_medium", "environmental_material"),
    "environment": ("environment",),
    "sample_type": ("sample_type",),
    "tissue": ("tissue",),
    "body_site": ("body_site",),
    "host_tissue_sampled": ("host_tissue_sampled",),
}
ECOLOGY_FIELDS = (
    "isolation_source", "geo_loc_name", "lat_lon", "host",
    "env_broad_scale", "env_local_scale", "env_medium", "environment",
    "tissue", "body_site", "host_tissue_sampled",
)
MISSING_VALUES = {
    "", "missing", "not provided", "not collected", "not available",
    "not applicable", "not determined", "not reported", "not known",
    "not recorded", "not supplied", "unspecified", "unknown", "n/a",
    "na", "null", "no data", "no information", "unavailable",
}
BASE_FIELDS = (
    "BioSample", "organism", "taxid", "title", "description", "package",
    "bioproject", "submission_date", "publication_date", "last_update",
)
OUTPUT_FIELDS = (*BASE_FIELDS, *FIELDS, "has_location", "has_ecology")


def tidy(value):
    return " ".join((value or "").split())


def key(value):
    return re.sub(r"[^a-z0-9]+", "_", tidy(value).lower()).strip("_")


def clean_value(value):
    value = tidy(value)
    return "" if value.lower() in MISSING_VALUES else value


def read_accessions(tsv_path):
    with tsv_path.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames or "BioSample" not in reader.fieldnames:
            raise ValueError(f"{tsv_path}: missing 'BioSample' column")
        ids = set()
        invalid = set()
        for row in reader:
            accession = tidy(row.get("BioSample", "")).upper()
            if not accession or accession in {"NA", "N/A", "-"}:
                continue
            if VALID_ID.fullmatch(accession):
                ids.add(accession)
            else:
                invalid.add(accession)
    if invalid:
        print(f"[WARN] Skipped {len(invalid)} invalid BioSample IDs: "
              + ", ".join(sorted(invalid)[:5]), file=sys.stderr)
    return sorted(ids)


def parse_xml(data):
    root = ET.fromstring(data)
    if root.tag == "BioSample":
        return [root]
    if root.tag != "BioSampleSet":
        raise ValueError(f"Unexpected NCBI response: XML root={root.tag!r}")
    return root.findall("./BioSample")


def download_batch(accessions, email, api_key):
    params = {
        "db": "biosample", "id": ",".join(accessions), "retmode": "xml",
        "tool": "genop_biosample_metadata",
    }
    if email:
        params["email"] = email
    if api_key:
        params["api_key"] = api_key
    request = Request(
        EFETCH, data=urlencode(params).encode("utf-8"),
        headers={"User-Agent": "GenOP-BioSample/1.0", "Accept": "application/xml"},
        method="POST",
    )
    for attempt in range(5):
        try:
            # Less than 3 E-utilities requests/second without an API key.
            time.sleep(0.4 if not api_key else 0.12)
            with urlopen(request, timeout=120) as response:
                data = response.read()
            parse_xml(data)  # Reject HTML/errors/malformed XML before caching.
            return data
        except (HTTPError, URLError, TimeoutError, ET.ParseError, ValueError) as exc:
            if attempt == 4:
                raise RuntimeError(f"NCBI request failed: {exc}") from exc
            pause = 2 ** attempt
            print(f"[WARN] NCBI error ({exc}); retrying in {pause} seconds",
                  file=sys.stderr)
            time.sleep(pause)
    raise AssertionError("Unreachable")


def get_batch(accessions, raw_dir, email, api_key, refresh):
    digest = hashlib.sha256(",".join(accessions).encode()).hexdigest()[:16]
    cache_file = raw_dir / f"biosample_{digest}.xml.gz"
    if cache_file.exists() and not refresh:
        try:
            with gzip.open(cache_file, "rb") as handle:
                data = handle.read()
            return parse_xml(data), True
        except (OSError, ET.ParseError, ValueError) as exc:
            print(f"[WARN] Corrupted cache file {cache_file}: {exc}", file=sys.stderr)

    data = download_batch(accessions, email, api_key)
    records = parse_xml(data)
    tmp_file = cache_file.with_suffix(cache_file.suffix + ".tmp")
    with gzip.open(tmp_file, "wb") as handle:
        handle.write(data)
    os.replace(tmp_file, cache_file)
    return records, False


def decode_sample(sample):
    accession = (sample.get("accession") or "").upper()
    organism = sample.find("./Description/Organism")
    links = []
    for link in sample.findall("./Links/Link"):
        if link.get("target", "").lower() == "bioproject":
            links.append(tidy(link.get("label") or link.text))

    row = {
        "BioSample": accession,
        "organism": (tidy(organism.get("taxonomy_name") or
                           organism.findtext("OrganismName")) if organism is not None else ""),
        "taxid": organism.get("taxonomy_id", "") if organism is not None else "",
        "title": tidy(sample.findtext("./Description/Title")),
        "description": tidy(" ".join(x.strip() for x in
                                      sample.find("./Description/Comment").itertext()))
                       if sample.find("./Description/Comment") is not None else "",
        "package": tidy(sample.findtext("./Package")),
        "bioproject": " | ".join(dict.fromkeys(filter(None, links))),
        "submission_date": sample.get("submission_date", ""),
        "publication_date": sample.get("publication_date", ""),
        "last_update": sample.get("last_update", ""),
    }

    attributes = []
    lookup = defaultdict(list)
    for attribute in sample.findall("./Attributes/Attribute"):
        name = tidy(attribute.get("attribute_name"))
        harmonized = tidy(attribute.get("harmonized_name"))
        value = tidy(" ".join(attribute.itertext()))
        attributes.append((accession, name, harmonized, value))
        usable = clean_value(value)
        if usable:
            for attr_key in {key(name), key(harmonized)} - {""}:
                if usable not in lookup[attr_key]:
                    lookup[attr_key].append(usable)

    for field, aliases in FIELDS.items():
        found = []
        for alias in aliases:
            for value in lookup.get(alias, []):
                if value not in found:
                    found.append(value)
        row[field] = " | ".join(found)

    row["has_location"] = int(bool(row["geo_loc_name"] or row["lat_lon"]))
    row["has_ecology"] = int(any(row[field] for field in ECOLOGY_FIELDS))
    return row, attributes


def write_table(path, fieldnames, rows):
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(fieldnames)
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, type=Path,
                        help="SRA metadata TSV containing a BioSample column")
    parser.add_argument("--genus", help="Output file prefix (e.g., Serratia)")
    parser.add_argument("--raw-dir", type=Path, default=Path("data/raw/biosample"))
    parser.add_argument("--out-dir", type=Path, default=Path("data/processed"))
    parser.add_argument("--batch-size", type=int, default=100)
    parser.add_argument("--limit", type=int, default=None,
                        help="Process at most N unique BioSample IDs (useful for testing)")
    parser.add_argument("--email", default=os.getenv("NCBI_EMAIL", ""))
    parser.add_argument("--api-key", default=os.getenv("NCBI_API_KEY", ""))
    parser.add_argument("--refresh", action="store_true",
                        help="Redownload batches even if the XML cache exists")
    args = parser.parse_args()
    if not 1 <= args.batch_size <= 200:
        parser.error("--batch-size must be between 1 and 200")
    if args.limit is not None and args.limit < 1:
        parser.error("--limit must be a positive integer")
    if not args.input.is_file():
        parser.error(f"Input file not found: {args.input}")
    genus = args.genus or args.input.name.removesuffix("_sra.tsv").removesuffix(".tsv")
    accessions = read_accessions(args.input)
    if not accessions:
        parser.error("No valid BioSample IDs found")
    total_accessions = len(accessions)
    if args.limit is not None:
        accessions = accessions[:args.limit]

    raw_dir = args.raw_dir / genus
    raw_dir.mkdir(parents=True, exist_ok=True)
    args.out_dir.mkdir(parents=True, exist_ok=True)
    print(f"BioSample ID: {len(accessions)} / {total_accessions} unique IDs found "
          f"| batch size: {args.batch_size}"
          + (f" | limit: {args.limit}" if args.limit is not None else ""))
    if not args.email:
        print("[INFO] Consider setting NCBI_EMAIL (recommended).", file=sys.stderr)

    samples = {}
    attrs_by_id = {}
    cache_hits = 0
    batches = [accessions[i:i + args.batch_size] for i in range(0, len(accessions), args.batch_size)]
    for n, batch in enumerate(batches, 1):
        records, cached = get_batch(batch, raw_dir, args.email, args.api_key, args.refresh)
        cache_hits += int(cached)
        for record in records:
            accession = (record.get("accession") or "").upper()
            if accession in batch:
                samples[accession], attrs_by_id[accession] = decode_sample(record)
        print(f"[{n}/{len(batches)}] {'CACHE' if cached else 'NCBI'}: "
              f"{sum(acc in samples for acc in batch)}/{len(batch)}")

    base = args.out_dir / f"{genus}_biosample"
    summary_rows = [samples[a] for a in accessions if a in samples]
    write_table(Path(str(base) + ".tsv"), OUTPUT_FIELDS,
                ([row.get(field, "") for field in OUTPUT_FIELDS] for row in summary_rows))
    write_table(Path(str(base) + "_ecology.tsv"), OUTPUT_FIELDS,
                ([row.get(field, "") for field in OUTPUT_FIELDS] for row in summary_rows
                 if row["has_ecology"]))
    write_table(Path(str(base) + "_attributes.tsv"),
                ("BioSample", "attribute_name", "harmonized_name", "value"),
                (attr for a in accessions for attr in attrs_by_id.get(a, [])))

    counts = Counter()
    for acc in accessions:
        counts.update(set((name, harm) for _, name, harm, _ in attrs_by_id.get(acc, [])))
    write_table(Path(str(base) + "_attribute_counts.tsv"),
                ("attribute_name", "harmonized_name", "n_biosamples", "pct_biosamples"),
                ((name, harm, count, f"{100 * count / len(summary_rows):.1f}")
                 for (name, harm), count in sorted(counts.items(), key=lambda p: (-p[1], p[0]))))
    missing = [a for a in accessions if a not in samples]
    write_table(Path(str(base) + "_missing.tsv"), ("BioSample",), ((a,) for a in missing))
    print(f"Completed. Retrieved: {len(summary_rows)} / {len(accessions)}; "
          f"with ecological metadata: {sum(row['has_ecology'] for row in summary_rows)}; "
          f"missing: {len(missing)}; cached batches: {cache_hits}/{len(batches)}")
    print(f"Output: {base}*.tsv")
    if missing:
        print("[WARN] Some BioSample records were not returned. Check *_missing.tsv; "
              "rerun with --refresh if needed.", file=sys.stderr)


if __name__ == "__main__":
    try:
        main()
    except (OSError, RuntimeError, ValueError) as exc:
        sys.exit(f"ERROR: {exc}")
