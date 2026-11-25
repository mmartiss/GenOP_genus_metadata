#!/usr/bin/env bash

# USAGE:
# bash scripts/08_download_gtdb_tree.sh

set -euo pipefail

PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

OUT_DIR="${PROJECT_DIR}/data/raw/gtdb"
OUT_FILE="${OUT_DIR}/bac120.tree.gz"

BASE_URL="${GTDB_BASE_URL:-https://data.gtdb.ecogenomic.org/releases/latest}"
URL="${BASE_URL}/bac120.tree.gz"

mkdir -p "${OUT_DIR}"

echo "[INFO] Checking GTDB bacterial tree..."

if [[ -s "${OUT_FILE}" ]] && gzip -t "${OUT_FILE}" 2>/dev/null; then
    echo "[INFO] GTDB tree already exists."
    echo "[INFO] ${OUT_FILE}"
    exit 0
fi

echo "[INFO] Downloading GTDB bacterial tree..."

curl -fL --retry 3 \
    "${URL}" \
    -o "${OUT_FILE}.tmp"

gzip -t "${OUT_FILE}.tmp"

mv "${OUT_FILE}.tmp" "${OUT_FILE}"

echo "[INFO] Download completed."
echo "[INFO] ${OUT_FILE}"
