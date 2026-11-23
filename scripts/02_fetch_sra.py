
#!/usr/bin/env python3

import argparse
import csv
import gzip
import io
import json
import time
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path

BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"


def request(endpoint, params):
    url = f"{BASE}/{endpoint}?" + urllib.parse.urlencode(params)

    for attempt in range(5):
        try:
            time.sleep(0.45)
            with urllib.request.urlopen(url, timeout=60) as response:
                return response.read().decode("utf-8-sig")
        except (urllib.error.URLError, TimeoutError) as exc:
            if attempt == 4:
                raise
            time.sleep(2 ** attempt)

    raise RuntimeError("NCBI request failed")


def get_biosamples(gtdb_file, genus, output):
    biosamples = set()
    genomes = []

    with gzip.open(gtdb_file, "rt") as handle:
        reader = csv.DictReader(handle, delimiter="\t")

        for row in reader:
            taxonomy = row["gtdb_taxonomy"].split(";")

            if f"g__{genus}" not in taxonomy:
                continue

            sample = row.get("ncbi_biosample", "").strip()

            genomes.append({
                "accession": row["accession"],
                "gtdb_taxonomy": row["gtdb_taxonomy"],
                "ncbi_biosample": sample,
                "ncbi_bioproject": row["ncbi_bioproject"],
                "ncbi_organism_name": row["ncbi_organism_name"],
            })

            if sample and sample.lower() not in {"none", "na", "n/a"}:
                biosamples.add(sample)

    output.parent.mkdir(parents=True, exist_ok=True)

    with output.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "accession",
                "gtdb_taxonomy",
                "ncbi_biosample",
                "ncbi_bioproject",
                "ncbi_organism_name",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(genomes)

    print(f"GTDB genomes: {len(genomes)}")
    print(f"Unique BioSamples: {len(biosamples)}")

    return sorted(biosamples)


def fetch_sra(sample):
    search_text = request("esearch.fcgi", {
        "db": "sra",
        "term": sample,
        "retmode": "json",
        "retmax": 0,
        "usehistory": "y",
    })

    result = json.loads(search_text)["esearchresult"]
    count = int(result["count"])

    if count == 0:
        return []

    rows = []

    # Read all matching records, not just the first page.
    for start in range(0, count, 200):
        text = request("efetch.fcgi", {
            "db": "sra",
            "query_key": result["querykey"],
            "WebEnv": result["webenv"],
            "rettype": "runinfo",
            "retmode": "text",
            "retstart": start,
            "retmax": 200,
        })

        reader = csv.DictReader(io.StringIO(text))

        if not reader.fieldnames or "Run" not in reader.fieldnames:
            raise ValueError("Unexpected SRA RunInfo response")

        for row in reader:
            # Keep only rows belonging to this BioSample.
            if row.get("BioSample", "").strip() == sample:
                rows.append(row)

    return rows


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--genus", required=True)
    parser.add_argument("--limit", type=int, default=0)
    parser.add_argument(
        "--gtdb",
        default="data/raw/gtdb/bac120_metadata.tsv.gz",
    )
    args = parser.parse_args()

    genus = args.genus.strip()

    processed = Path("data/processed")
    cache = Path("data/raw/ncbi/sra_cache") / genus
    cache.mkdir(parents=True, exist_ok=True)

    gtdb_output = processed / f"{genus}_gtdb.tsv"
    biosamples = get_biosamples(
        args.gtdb, genus, gtdb_output
    )

    if args.limit:
        biosamples = biosamples[:args.limit]

    all_rows = []
    status_rows = []

    for index, sample in enumerate(biosamples, 1):
        cached = cache / f"{sample}.json"

        try:
            if cached.exists():
                rows = json.loads(cached.read_text())
            else:
                rows = fetch_sra(sample)
                cached.write_text(json.dumps(rows))

            all_rows.extend(rows)

            status_rows.append({
                "BioSample": sample,
                "status": "found" if rows else "no_runs",
                "run_count": len(rows),
            })

            print(
                f"[{index}/{len(biosamples)}] "
                f"{sample}: {len(rows)} runs"
            )

        except Exception as exc:
            status_rows.append({
                "BioSample": sample,
                "status": "error",
                "run_count": "",
            })
            print(f"ERROR {sample}: {exc}")

    processed.mkdir(parents=True, exist_ok=True)

    # Preserve all fields returned by NCBI.
    fields = list(dict.fromkeys(
        key for row in all_rows for key in row.keys()
    ))

    output = processed / f"{genus}_sra.tsv"

    with output.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fields or ["Run", "BioSample"],
            delimiter="\t",
            extrasaction="ignore",
        )
        writer.writeheader()
        writer.writerows(all_rows)

    status_output = processed / f"{genus}_sra_status.tsv"

    with status_output.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["BioSample", "status", "run_count"],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(status_rows)

    print(f"\nSaved: {output}")
    print(f"Status: {status_output}")
    print(f"Total SRA runs: {len(all_rows)}")


if __name__ == "__main__":
    main()
