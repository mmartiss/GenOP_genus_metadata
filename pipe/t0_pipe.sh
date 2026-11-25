#!/usr/bin/env bash

set -euo pipefail

# Set project root directory
PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$PROJECT_DIR"

GENUS="${1:?Usage: bash pipe/t0_pipe.sh <genus>}"

FILTERED="data/filtered/${GENUS,,}_genera.tsv.gz"
PROCESSED="data/processed"

mkdir -p data/filtered "$PROCESSED"

# 00 — Download GTDB metadata
bash scripts/00_download_gtdb.sh

# 01 — Filter GTDB metadata
bash scripts/01_filter_genus_gtdb.sh \
    --input data/raw/gtdb/bac120_metadata.tsv.gz \
    --output "$FILTERED" \
    --genus "$GENUS" \
    --force

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
Rscript scripts/05_clean_metadata.r \
    "$PROCESSED/${GENUS}_merged_metadata.tsv" \
    "$PROCESSED/${GENUS}_clean_metadata.tsv"

# 06 - Standardize metadata
echo "[06] Standardizing metadata..."

Rscript scripts/06_standardize_metadata.r \
    "data/processed/${GENUS}_merged_metadata.tsv" \
    "data/processed/${GENUS}_standardized_metadata.tsv"

# 07 - Species summary
echo "[07] Generating species summary..."

Rscript scripts/07_species_summary.r \
    "data/processed/${GENUS}_standardized_metadata.tsv" \
    "data/processed/${GENUS}_species_summary.tsv"

# 08 - Download GTDB tree
echo "[STEP 08] Downloading GTDB tree..."

bash scripts/08_download_gtdb_tree.sh

echo "[DONE] Pipeline completed for ${GENUS}"
