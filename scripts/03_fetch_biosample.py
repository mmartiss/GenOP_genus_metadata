#!/usr/bin/env python3
"""Fetch BioSample metadata for accessions in an SRA RunInfo TSV.

Uses NCBI BioSample EFetch in batches and falls back to the EMBL-EBI BioSamples
API for accessions not returned by NCBI. Supports SAMN, SAMEA and SAMD IDs.
Raw XML/JSON is cached to allow reruns without unnecessary downloads.
Python 3.9+; no third-party packages.

Example usage:
    python3 scripts/03_fetch_biosample.py \
    --input data/processed/Serratia_sra.tsv \
    --genus Serratia \
    --limit 10

"""

import argparse
import csv
import gzip
import hashlib
import json
import os
import re
import sys
import time
from collections import Counter, defaultdict
from pathlib import Path
from urllib.error import HTTPError, URLError
from urllib.parse import urlencode
from urllib.request import Request, urlopen
from xml.etree import ElementTree as ET

NCBI_EFETCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
EBI_SAMPLE_API = "https://www.ebi.ac.uk/biosamples/samples/"
# INSDC accessions: SAMN (NCBI), SAMD (DDBJ), SAMEA (EMBL-EBI).
# SAME and SAMEG are also allowed by the EMBL-EBI accession schema.
VALID_ID = re.compile(r"^SAM(?:N|D|E(?:A|G)?)\d+$", re.IGNORECASE)

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
    "metadata_source",
)
OUTPUT_FIELDS = (*BASE_FIELDS, *FIELDS, "has_location", "has_ecology")


def tidy(value):
    return " ".join(str(value or "").split())


def key(value):
    return re.sub(r"[^a-z0-9]+", "_", tidy(value).lower()).strip("_")


def clean_value(value):
    value = tidy(value)
    return "" if value.lower() in MISSING_VALUES else value


def read_accessions(path):
    with path.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames or "BioSample" not in reader.fieldnames:
            raise ValueError(f"{path}: missing BioSample column")
        accessions, invalid = set(), set()
        for row in reader:
            accession = tidy(row.get("BioSample")).upper()
            if not accession or accession in {"NA", "N/A", "-"}:
                continue
            if VALID_ID.fullmatch(accession):
                accessions.add(accession)
            else:
                invalid.add(accession)
    if invalid:
        print(f"[WARN] Ignored {len(invalid)} malformed BioSample IDs: "
              + ", ".join(sorted(invalid)[:5]), file=sys.stderr)
    return sorted(accessions)


def xml_records(data):
    root = ET.fromstring(data)
    if root.tag == "BioSample":
        return [root]
    if root.tag != "BioSampleSet":
        raise ValueError(f"Unexpected NCBI XML root: {root.tag!r}")
    return root.findall("./BioSample")


def request_ncbi(accessions, email, api_key):
    params = {
        "db": "biosample", "id": ",".join(accessions), "retmode": "xml",
        "tool": "GenOP_BioSample",
    }
    if email:
        params["email"] = email
    if api_key:
        params["api_key"] = api_key
    for attempt in range(4):
        request = Request(
            NCBI_EFETCH,
            data=urlencode(params).encode("utf-8"),
            headers={"User-Agent": "GenOP/1.1", "Accept": "application/xml"},
            method="POST",
        )
        try:
            time.sleep(0.4 if not api_key else 0.12)
            with urlopen(request, timeout=90) as response:
                data = response.read()
            xml_records(data)
            return data
        except (HTTPError, URLError, TimeoutError, ET.ParseError, ValueError) as exc:
            if attempt == 3:
                raise RuntimeError(f"NCBI EFetch failed after retries: {exc}") from exc
            time.sleep(2 ** attempt)
    raise AssertionError("Unreachable")


def ncbi_batch(accessions, raw_dir, email, api_key, refresh):
    digest = hashlib.sha256(",".join(accessions).encode()).hexdigest()[:16]
    cache_file = raw_dir / f"biosample_{digest}.xml.gz"
    if cache_file.is_file() and not refresh:
        try:
            with gzip.open(cache_file, "rb") as handle:
                return xml_records(handle.read()), True
        except (OSError, ValueError, ET.ParseError) as exc:
            print(f"[WARN] Invalid NCBI cache ({cache_file}): {exc}", file=sys.stderr)
    data = request_ncbi(accessions, email, api_key)
    records = xml_records(data)
    temp = cache_file.with_suffix(".xml.gz.tmp")
    with gzip.open(temp, "wb") as handle:
        handle.write(data)
    os.replace(temp, cache_file)
    return records, False


def ebi_sample(accession, raw_dir, refresh):
    """Return (sample_dict, cached), (None, False) on HTTP 404; raise on errors."""
    cache_file = raw_dir / f"{accession}.json.gz"
    if cache_file.is_file() and not refresh:
        try:
            with gzip.open(cache_file, "rt", encoding="utf-8") as handle:
                record = json.load(handle)
            if isinstance(record, dict) and tidy(record.get("accession")).upper() == accession:
                return record, True
            print(f"[WARN] Invalid EBI cache for {accession}; redownloading", file=sys.stderr)
        except (OSError, json.JSONDecodeError) as exc:
            print(f"[WARN] Invalid EBI cache for {accession}: {exc}", file=sys.stderr)

    for attempt in range(4):
        request = Request(
            EBI_SAMPLE_API + accession,
            headers={"Accept": "application/json", "User-Agent": "GenOP/1.1"},
        )
        try:
            time.sleep(0.15)
            with urlopen(request, timeout=45) as response:
                data = response.read()
            record = json.loads(data)
            if not isinstance(record, dict) or tidy(record.get("accession")).upper() != accession:
                raise ValueError(f"EBI returned unexpected accession for {accession}")
            temp = cache_file.with_suffix(".json.gz.tmp")
            with gzip.open(temp, "wb") as handle:
                handle.write(data)
            os.replace(temp, cache_file)
            return record, False
        except HTTPError as exc:
            if exc.code == 404:
                return None, False
            error = exc
        except (URLError, TimeoutError, json.JSONDecodeError, ValueError) as exc:
            error = exc
        if attempt == 3:
            raise RuntimeError(f"EBI fetch failed for {accession}: {error}") from error
        time.sleep(2 ** attempt)
    raise AssertionError("Unreachable")


def normalize_attributes(attributes, row):
    lookup = defaultdict(list)
    for _, name, harmonized, value in attributes:
        usable = clean_value(value)
        if usable:
            for attribute_key in {key(name), key(harmonized)} - {""}:
                if usable not in lookup[attribute_key]:
                    lookup[attribute_key].append(usable)
    for field, aliases in FIELDS.items():
        values = []
        for alias in aliases:
            for value in lookup.get(alias, []):
                if value not in values:
                    values.append(value)
        row[field] = " | ".join(values)
    row["has_location"] = int(bool(row["geo_loc_name"] or row["lat_lon"]))
    row["has_ecology"] = int(any(row[field] for field in ECOLOGY_FIELDS))
    return row, attributes


def decode_ncbi(sample):
    accession = tidy(sample.get("accession")).upper()
    organism = sample.find("./Description/Organism")
    comment = sample.find("./Description/Comment")
    links = []
    for link in sample.findall("./Links/Link"):
        if link.get("target", "").lower() == "bioproject":
            links.append(tidy(link.get("label") or link.text))
    row = {
        "BioSample": accession,
        "organism": tidy(organism.get("taxonomy_name") or organism.findtext("OrganismName")) if organism is not None else "",
        "taxid": tidy(organism.get("taxonomy_id")) if organism is not None else "",
        "title": tidy(sample.findtext("./Description/Title")),
        "description": tidy(" ".join(comment.itertext())) if comment is not None else "",
        "package": tidy(sample.findtext("./Package")),
        "bioproject": " | ".join(dict.fromkeys(filter(None, links))),
        "submission_date": tidy(sample.get("submission_date")),
        "publication_date": tidy(sample.get("publication_date")),
        "last_update": tidy(sample.get("last_update")),
        "metadata_source": "NCBI",
    }
    attrs = []
    for attr in sample.findall("./Attributes/Attribute"):
        attrs.append((accession, tidy(attr.get("attribute_name")),
                      tidy(attr.get("harmonized_name")), tidy(" ".join(attr.itertext()))))
    return normalize_attributes(attrs, row)


def decode_ebi(sample):
    accession = tidy(sample["accession"]).upper()
    attrs = []
    for name, entries in (sample.get("characteristics") or {}).items():
        if not isinstance(entries, list):
            entries = [entries]
        for entry in entries:
            value = entry.get("text", "") if isinstance(entry, dict) else entry
            attrs.append((accession, tidy(name), "", tidy(value)))
    attr_index = defaultdict(list)
    for _, name, _, value in attrs:
        if clean_value(value):
            attr_index[key(name)].append(tidy(value))

    def first(*names):
        for name in names:
            if attr_index.get(key(name)):
                return attr_index[key(name)][0]
        return ""

    refs = []
    for ref in sample.get("externalReferences", []) or []:
        url = ref.get("url", "") if isinstance(ref, dict) else ""
        refs.extend(re.findall(r"PRJ(?:NA|EB|DB)\d+", url, re.IGNORECASE))
    row = {
        "BioSample": accession,
        "organism": first("organism", "scientific name") or tidy(sample.get("organism")),
        "taxid": tidy(sample.get("taxId")) or first("NCBI taxon id", "taxon id"),
        "title": first("title", "sample title") or tidy(sample.get("name")),
        "description": first("description", "sample description"),
        "package": first("model", "biosample package"),
        "bioproject": " | ".join(dict.fromkeys(refs)),
        "submission_date": tidy(sample.get("create")),
        "publication_date": tidy(sample.get("release")),
        "last_update": tidy(sample.get("update")),
        "metadata_source": "EBI",
    }
    return normalize_attributes(attrs, row)


def write_tsv(path, columns, rows):
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(columns)
        writer.writerows(rows)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True,
                        help="Input SRA TSV with a BioSample column")
    parser.add_argument("--genus", help="Genus prefix for output files, e.g. Serratia")
    parser.add_argument("--raw-dir", type=Path, default=Path("data/raw/biosample"))
    parser.add_argument("--out-dir", type=Path, default=Path("data/processed"))
    parser.add_argument("--batch-size", type=int, default=100)
    parser.add_argument("--limit", type=int, help="Process only the first N unique BioSample IDs")
    parser.add_argument("--email", default=os.getenv("NCBI_EMAIL", ""))
    parser.add_argument("--api-key", default=os.getenv("NCBI_API_KEY", ""))
    parser.add_argument("--refresh", action="store_true",
                        help="Ignore cached XML/JSON and download again")
    args = parser.parse_args(argv)

    if not 1 <= args.batch_size <= 200:
        parser.error("--batch-size must be between 1 and 200")
    if args.limit is not None and args.limit < 1:
        parser.error("--limit must be positive")
    if not args.input.is_file():
        parser.error(f"Input file not found: {args.input}")

    genus = args.genus or args.input.name.removesuffix("_sra.tsv").removesuffix(".tsv")
    all_accessions = read_accessions(args.input)
    if not all_accessions:
        parser.error("No valid BioSample IDs found")
    accessions = all_accessions[:args.limit] if args.limit else all_accessions
    raw_dir = args.raw_dir / genus
    raw_dir.mkdir(parents=True, exist_ok=True)
    args.out_dir.mkdir(parents=True, exist_ok=True)
    print(f"BioSamples: {len(accessions)} / {len(all_accessions)} unique IDs "
          f"| batch size: {args.batch_size}"
          + (f" | limit: {args.limit}" if args.limit else ""), flush=True)
    if not args.email:
        print("[INFO] NCBI_EMAIL is recommended for NCBI requests.", file=sys.stderr)

    samples, attributes, missing, errors = {}, {}, [], []
    cache_hits = 0
    for start in range(0, len(accessions), args.batch_size):
        batch = accessions[start:start + args.batch_size]
        n = 1 + start // args.batch_size
        n_batches = (len(accessions) + args.batch_size - 1) // args.batch_size
        ncbi_failed = False
        try:
            records, from_cache = ncbi_batch(batch, raw_dir, args.email,
                                             args.api_key, args.refresh)
            cache_hits += int(from_cache)
        except RuntimeError as exc:
            print(f"[WARN] {exc}; attempting EBI fallback for this batch", file=sys.stderr)
            records, from_cache, ncbi_failed = [], False, True

        for record in records:
            accession = tidy(record.get("accession")).upper()
            if accession in batch and accession not in samples:
                samples[accession], attributes[accession] = decode_ncbi(record)
        found_ncbi = sum(acc in samples for acc in batch)
        print(f"[{n}/{n_batches}] NCBI: {found_ncbi}/{len(batch)}"
              + (" (cached)" if from_cache else ""), flush=True)

        for accession in batch:
            if accession in samples:
                continue
            try:
                record, cached = ebi_sample(accession, raw_dir, args.refresh)
                cache_hits += int(cached)
                if record is None:
                    if ncbi_failed:
                        errors.append((accession, "NCBI request failed; EBI returned 404"))
                    else:
                        missing.append(accession)
                else:
                    samples[accession], attributes[accession] = decode_ebi(record)
            except (RuntimeError, OSError, ValueError) as exc:
                errors.append((accession, str(exc)))
                print(f"[ERROR] {accession}: {exc}", file=sys.stderr)
        print(f"[{n}/{n_batches}] Total retrieved: "
              f"{sum(a in samples for a in batch)}/{len(batch)} after EBI fallback",
              flush=True)

    base = args.out_dir / f"{genus}_biosample"
    rows = [samples[acc] for acc in accessions if acc in samples]
    write_tsv(Path(f"{base}.tsv"), OUTPUT_FIELDS,
              ([row.get(f, "") for f in OUTPUT_FIELDS] for row in rows))
    write_tsv(Path(f"{base}_ecology.tsv"), OUTPUT_FIELDS,
              ([row.get(f, "") for f in OUTPUT_FIELDS] for row in rows if row["has_ecology"]))
    write_tsv(Path(f"{base}_attributes.tsv"),
              ("BioSample", "attribute_name", "harmonized_name", "value"),
              (attr for acc in accessions for attr in attributes.get(acc, [])))
    counts = Counter()
    for acc in accessions:
        counts.update({(name, harmonized) for _, name, harmonized, _ in attributes.get(acc, [])})
    write_tsv(Path(f"{base}_attribute_counts.tsv"),
              ("attribute_name", "harmonized_name", "n_biosamples", "pct_biosamples"),
              ((name, harmonized, n, f"{100 * n / len(rows):.1f}")
               for (name, harmonized), n in sorted(counts.items(),
                                                   key=lambda pair: (-pair[1], pair[0]))))
    write_tsv(Path(f"{base}_missing.tsv"), ("BioSample",), ((acc,) for acc in missing))
    write_tsv(Path(f"{base}_errors.tsv"), ("BioSample", "error"), errors)
    from_ncbi = sum(row["metadata_source"] == "NCBI" for row in rows)
    print(f"Completed. Retrieved: {len(rows)}/{len(accessions)} "
          f"(NCBI: {from_ncbi}, EBI: {len(rows) - from_ncbi}); "
          f"ecological metadata: {sum(row['has_ecology'] for row in rows)}; "
          f"missing: {len(missing)}; errors: {len(errors)}; cache hits: {cache_hits}")
    print(f"Output: {base}*.tsv")
    if errors:
        print(f"[ERROR] {len(errors)} records failed; check {base}_errors.tsv.", file=sys.stderr)
        return 1
    if missing:
        print(f"[WARN] {len(missing)} BioSample records absent from NCBI and EBI. "
              "See *_missing.tsv.", file=sys.stderr)
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except (OSError, ValueError, RuntimeError) as error:
        sys.exit(f"ERROR: {error}")
