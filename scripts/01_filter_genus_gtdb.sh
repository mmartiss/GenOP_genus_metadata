#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'USAGE'
Usage:
  filter_gtdb.sh --input FILE --output FILE --genus GENUS [GENUS ...] [--force]

Example:
  bash scripts/01_filter_genus_gtdb.sh \
    --input data/raw/gtdb/bac120_metadata.tsv.gz \
    --output data/filtered/selected_genera.tsv.gz \
    --genus Serratia Bacillus

Input/output: .tsv or .tsv.gz
--force: overwrite existing output file
USAGE
}

input=""
output=""
force=false
genera=()

while (( $# > 0 )); do
  case "$1" in
    --input|-i)
      (( $# >= 2 )) || { echo "Missing value for --input" >&2; exit 2; }
      input=$2; shift 2 ;;
    --output|-o)
      (( $# >= 2 )) || { echo "Missing value for --output" >&2; exit 2; }
      output=$2; shift 2 ;;
    --genus|-g)
      shift
      while (( $# > 0 )) && [[ $1 != -* ]]; do
        genera+=("$1")
        shift
      done ;;
    --force)
      force=true; shift ;;
    --help|-h)
      usage; exit 0 ;;
    *) echo "Unknown argument: $1" >&2; usage >&2; exit 2 ;;
  esac
done

[[ -n $input && -n $output && ${#genera[@]} -gt 0 ]] || { usage >&2; exit 2; }
[[ -f $input && -r $input ]] || { echo "Cannot read input: $input" >&2; exit 1; }
[[ $input == *.tsv || $input == *.tsv.gz ]] || { echo "Input must be .tsv or .tsv.gz" >&2; exit 2; }
[[ $output == *.tsv || $output == *.tsv.gz ]] || { echo "Output must be .tsv or .tsv.gz" >&2; exit 2; }
[[ $input != "$output" ]] || { echo "Input and output must differ" >&2; exit 2; }
if [[ -e $output && $force != true ]]; then
  echo "Output already exists: $output (use --force to overwrite)" >&2
  exit 1
fi

for genus in "${genera[@]}"; do
  [[ $genus =~ ^[A-Za-z0-9_]+$ ]] || {
    echo "Invalid genus '$genus': use a GTDB genus name such as Serratia or Serratia_A" >&2
    exit 2
  }
done

mkdir -p "$(dirname "$output")"
tmp=$(mktemp "$(dirname "$output")/.filter_gtdb.XXXXXXXX")
trap 'rm -f -- "$tmp"' EXIT

genus_list=$(printf '%s\n' "${genera[@]}")

read_input() {
  if [[ $input == *.gz ]]; then
    gzip -dc -- "$input"
  else
    cat -- "$input"
  fi
}

filter_rows() {
  awk -F '\t' -v OFS='\t' -v wanted="$genus_list" '
    BEGIN {
      n = split(wanted, names, "\n")
      for (i = 1; i <= n; i++) if (names[i] != "") selected[names[i]] = 1
    }
    NR == 1 {
      for (i = 1; i <= NF; i++) {
        col = $i
        sub(/\r$/, "", col)
        if (col == "gtdb_taxonomy") tax_col = i
      }
      if (!tax_col) {
        print "ERROR: missing gtdb_taxonomy column" > "/dev/stderr"
        exit 3
      }
      print
      next
    }
    {
      total++
      k = split($tax_col, ranks, ";")
      genus = ""
      for (i = 1; i <= k; i++) {
        if (substr(ranks[i], 1, 3) == "g__") {
          genus = substr(ranks[i], 4)
          sub(/\r$/, "", genus)
          break
        }
      }
      if (genus in selected) {
        print
        matched++
        counts[genus]++
      }
    }
    END {
      if (!tax_col) exit 3
      printf "Processed genomes: %d\nMatched genomes: %d\n", total, matched > "/dev/stderr"
      for (i = 1; i <= n; i++) {
        if (names[i] != "" && !reported[names[i]]++)
          printf "  %s: %d\n", names[i], counts[names[i]] > "/dev/stderr"
      }
    }
  '
}

if [[ $output == *.gz ]]; then
  read_input | filter_rows | gzip -c > "$tmp"
else
  read_input | filter_rows > "$tmp"
fi

mv -f -- "$tmp" "$output"
trap - EXIT
printf 'Saved: %s\n' "$output" >&2
