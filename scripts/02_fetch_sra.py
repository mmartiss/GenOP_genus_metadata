
#!/usr/bin/env python3

# [USAGE]
# python3 scripts/02_fetch_sra.py \
#     --genus Serratia \
#     --input data/filtered/serratia_genera.tsv.gz

# Example
# python3 scripts/02_fetch_sra.py \
#     --genus Serratia \
#     --input data/filtered/serratia_genera.tsv.gz \
#     --limit 10

# Example
# python3 scripts/02_fetch_sra.py \
#     --genus Serratia \
#     --input data/filtered/serratia_genera.tsv.gz \
#     --refresh

import argparse
import csv
import gzip
import io
import json
import os
import time
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path

BASE_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
MISSING = {"", "none", "na", "n/a", "null"}

# Conservative NCBI request interval (without API key)
REQUEST_INTERVAL = 0.45
LAST_REQUEST = 0.0


def ncbi_request(endpoint, params):
    """Request data from NCBI E-utilities with retries."""
    global LAST_REQUEST

    params = {
        **params,
        "tool": "GenOP_genus_metadata",
    }

    url = (
        f"{BASE_URL}/{endpoint}?"
        + urllib.parse.urlencode(params)
    )

    for attempt in range(5):
        try:
            elapsed = time.monotonic() - LAST_REQUEST
            if elapsed < REQUEST_INTERVAL:
                time.sleep(REQUEST_INTERVAL - elapsed)

            LAST_REQUEST = time.monotonic()

            req = urllib.request.Request(
                url,
                headers={"User-Agent": "GenOP/0.1"},
            )

            with urllib.request.urlopen(
                req, timeout=60
            ) as response:
                return response.read().decode("utf-8-sig")

        except (
            urllib.error.URLError,
            TimeoutError,
        ) as exc:
            if attempt == 4:
                raise RuntimeError(
                    f"NCBI request failed: {exc}"
                ) from exc

            time.sleep(2 ** attempt)


def get_biosamples(input_file):
    """Read unique BioSample accessions from filtered GTDB."""
    biosamples = set()
    genome_count = 0
    missing_count = 0

    opener = (
        gzip.open
        if str(input_file).endswith(".gz")
        else open
    )

    with opener(input_file, "rt") as handle:
        reader = csv.DictReader(handle, delimiter="\t")

        if "ncbi_biosample" not in (reader.fieldnames or []):
            raise ValueError(
                "Input file lacks 'ncbi_biosample' column"
            )

        for row in reader:
            genome_count += 1

            sample = row["ncbi_biosample"].strip()

            if sample.lower() in MISSING:
                missing_count += 1
                continue

            biosamples.add(sample)

    print(f"GTDB genome records: {genome_count}")
    print(f"Missing BioSample: {missing_count}")
    print(f"Unique BioSamples: {len(biosamples)}")

    return sorted(biosamples)


def fetch_sra(sample):
    """Retrieve SRA RunInfo metadata for one BioSample."""

    search_text = ncbi_request(
        "esearch.fcgi",
        {
            "db": "sra",
            "term": f'"{sample}"',
            "retmode": "json",
            "retmax": 0,
            "usehistory": "y",
        },
    )

    search = json.loads(search_text)["esearchresult"]

    count = int(search["count"])

    if count == 0:
        return []

    rows = []

    for start in range(0, count, 200):
        runinfo_text = ncbi_request(
            "efetch.fcgi",
            {
                "db": "sra",
                "query_key": search["querykey"],
                "WebEnv": search["webenv"],
                "rettype": "runinfo",
                "retmode": "text",
                "retstart": start,
                "retmax": 200,
            },
        )

        reader = csv.DictReader(io.StringIO(runinfo_text))

        if not reader.fieldnames or "Run" not in reader.fieldnames:
            raise ValueError(
                f"Invalid SRA RunInfo response: {sample}"
            )

        if "BioSample" not in reader.fieldnames:
            raise ValueError(
                f"Missing BioSample field: {sample}"
            )

        for row in reader:
            if row["BioSample"].strip() == sample:
                rows.append(row)

    # Deduplicate by SRA Run accession
    unique_rows = {}

    for row in rows:
        run = row.get("Run", "").strip()

        if run:
            unique_rows[run] = row

    return list(unique_rows.values())


def load_or_fetch(sample, cache_dir, refresh=False):
    """Use cached results or fetch from NCBI."""
    cache_file = cache_dir / f"{sample}.json"

    if cache_file.exists() and not refresh:
        return json.loads(
            cache_file.read_text(encoding="utf-8")
        )

    rows = fetch_sra(sample)

    # Write atomically to avoid partially written cache
    temporary = cache_file.with_suffix(".tmp")

    temporary.write_text(
        json.dumps(rows, ensure_ascii=False),
        encoding="utf-8",
    )

    os.replace(temporary, cache_file)

    return rows


def write_tsv(path, rows, preferred_fields=None):
    """Write a list of dictionaries to TSV."""
    path.parent.mkdir(parents=True, exist_ok=True)

    fields = list(preferred_fields or [])

    for row in rows:
        for key in row:
            if key not in fields:
                fields.append(key)

    with path.open(
        "w", newline="", encoding="utf-8"
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fields,
            delimiter="\t",
            extrasaction="ignore",
        )

        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(
        description="Fetch NCBI SRA metadata for filtered GTDB genomes"
    )

    parser.add_argument(
        "--input",
        required=True,
        help="Filtered GTDB TSV or TSV.GZ file",
    )

    parser.add_argument(
        "--genus",
        required=True,
        help="Genus name used for output filenames",
    )

    parser.add_argument(
        "--limit",
        type=int,
        default=0,
        help="Number of BioSamples to test (0 = all)",
    )

    parser.add_argument(
        "--refresh",
        action="store_true",
        help="Fetch again even if cached",
    )

    args = parser.parse_args()

    input_file = Path(args.input)

    if not input_file.exists():
        parser.error(f"Input file not found: {input_file}")

    if args.limit < 0:
        parser.error("--limit must be >= 0")

    genus = args.genus.strip()

    processed_dir = Path("data/processed")
    cache_dir = Path("data/raw/ncbi/sra_cache") / genus

    processed_dir.mkdir(parents=True, exist_ok=True)
    cache_dir.mkdir(parents=True, exist_ok=True)

    biosamples = get_biosamples(input_file)

    if args.limit:
        biosamples = biosamples[:args.limit]

    all_runs = {}
    status_rows = []

    print(f"\nProcessing {len(biosamples)} BioSamples\n")

    for index, sample in enumerate(biosamples, start=1):
        try:
            rows = load_or_fetch(
                sample,
                cache_dir,
                refresh=args.refresh,
            )

            for row in rows:
                run = row.get("Run", "").strip()

                if run:
                    all_runs[run] = row

            status = "found" if rows else "no_runs"

            status_rows.append({
                "BioSample": sample,
                "status": status,
                "run_count": len(rows),
                "error": "",
            })

            print(
                f"[{index}/{len(biosamples)}] "
                f"{sample}: {status} ({len(rows)} runs)"
            )

        except Exception as exc:
            status_rows.append({
                "BioSample": sample,
                "status": "error",
                "run_count": "",
                "error": str(exc),
            })

            print(
                f"[{index}/{len(biosamples)}] "
                f"{sample}: ERROR - {exc}"
            )

    sra_output = processed_dir / f"{genus}_sra.tsv"
    status_output = (
        processed_dir / f"{genus}_sra_status.tsv"
    )

    write_tsv(
        sra_output,
        list(all_runs.values()),
        preferred_fields=["Run", "BioSample", "BioProject"],
    )

    write_tsv(
        status_output,
        status_rows,
        preferred_fields=[
            "BioSample",
            "status",
            "run_count",
            "error",
        ],
    )

    found = sum(
        row["status"] == "found"
        for row in status_rows
    )

    no_runs = sum(
        row["status"] == "no_runs"
        for row in status_rows
    )

    errors = sum(
        row["status"] == "error"
        for row in status_rows
    )

    print("\n=== Summary ===")
    print(f"BioSamples processed: {len(biosamples)}")
    print(f"BioSamples with SRA: {found}")
    print(f"BioSamples without SRA: {no_runs}")
    print(f"Errors: {errors}")
    print(f"Unique SRA runs: {len(all_runs)}")
    print(f"\nSRA metadata: {sra_output}")
    print(f"Status file: {status_output}")


if __name__ == "__main__":
    main()
