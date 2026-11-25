# GenOP — stage 06: metadata standardization

Run from repository root:

```bash
Rscript scripts/06_standardization/06_standardize_metadata.R \
  data/processed/Serratia_merged_metadata.tsv \
  data/processed/Serratia_standardized_metadata.tsv
```

Alternatively use a third argument pointing to a directory containing the seven mapping TSV files. Requires base R only. Keep mappings with the script. Original columns and rows are retained. Prefer the cleaned merged TSV if its relevant column names are unchanged. This first version is based on the supplied Serratia file, not a universal ontology. New exact values are preserved and marked unmapped; extend mapping tables for additional genera. `preserved` means no synonym normalization was performed, not validated biological annotation. Counts in mapping files count source-field occurrences, not unique genomes. `std_country` and `std_geography` can contain multiple values separated by ` | `. Do not treat those as single map locations. `std_isolation_source_disagreement` is a textual mismatch, not proof of biological conflict. Coordinates with >0.01 degree component disagreement are withheld and flagged. Dates are checked for shape, not full calendar validity. Clinical categories are conservative, evidence-based rule labels, not diagnostic conclusions; absence of evidence is Unknown. Multiple ecological labels can coexist. Historical USSR and ambiguous INP require review. R was not available in the authoring environment, so execution is not verified. `field_audit.tsv` was generated directly from the uploaded TSV using Python. Runtime audit/review/conflicts outputs are produced by R.
