
#!/usr/bin/env bash
set -euo pipefail

PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

DATA_DIR="${PROJECT_DIR}/data/raw/gtdb"
FILE_NAME="bac120_metadata.tsv.gz"
URL="https://data.gtdb.ecogenomic.org/releases/latest/${FILE_NAME}"

mkdir -p "$DATA_DIR"

OUTPUT="${DATA_DIR}/${FILE_NAME}"

if [[ -s "$OUTPUT" ]] && gzip -t "$OUTPUT" 2>/dev/null; then
    echo "[INFO] GTDB metadata already exists: $OUTPUT"
    echo "[INFO] Skipping download."
    exit 0
fi

echo "[INFO] Downloading GTDB bacterial metadata..."

TEMP_FILE="${OUTPUT}.tmp"
trap 'rm -f "$TEMP_FILE"' EXIT

curl -fL --retry 3 --connect-timeout 30 \
    -o "$TEMP_FILE" "$URL"

if ! gzip -t "$TEMP_FILE"; then
    echo "[ERROR] Downloaded file is corrupted."
    exit 1
fi

mv "$TEMP_FILE" "$OUTPUT"

echo "[INFO] Download completed: $OUTPUT"
