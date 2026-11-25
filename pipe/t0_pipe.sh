#!/usr/bin/env bash

set -euo pipefail

cd "$(dirname "${BASH_SOURCE[0]}")"

GENUS="${1:?Usage: bash run_pipeline.sh <genus>}"

FILTERED="data/filtered/${GENUS,,}_genera.tsv.gz"
PROCESSED="data/processed"

mkdir -p data/filtered "$PROCESSED"

# 00 — Download GTDB metadata if missing
bash scripts/00_download_gtdb.sh

# 01 — Filter GTDB metadata
bash scripts/01_filter_genus_gtdb.sh \
    --input data/raw/gtdb/bac120_metadata.tsv.gz \
    --output "$FILTERED" \
    --genus "$GENUS"

# 02 — Fetch SRA metadata
python3 scripts/02_fetch_sra.py \
    --genus "$GENUS" \
    --input "$FILTERED"

# 03 — Fetch BioSample metadata
python3 scripts/03_fetch_biosample.py \
    --input "$PROCESSED/${GENUS}_sra.tsv" \
    --genus "$GENUS"

# 04 — Merge metadata
python3 scripts/04_merge_metadata.py \
    --gtdb "$FILTERED" \
    --biosample "$PROCESSED/${GENUS}_biosample.tsv" \
    --sra "$PROCESSED/${GENUS}_sra.tsv" \
    --out "$PROCESSED/${GENUS}_merged_metadata.tsv"

# 05 — Clean metadata
Rscript scripts/05_clean_metadata.R \
    "$PROCESSED/${GENUS}_merged_metadata.tsv" \
    "$PROCESSED/${GENUS}_clean_metadata.tsv"

echo "Pipeline completed: $GENUS"
