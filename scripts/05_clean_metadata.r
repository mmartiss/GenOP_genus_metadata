#!/usr/bin/env Rscript

# 05_clean_metadata.R
# 1) Keep a curated set of metadata fields.
# 2) Standardize blank/missing values, without merging sources.
# 3) Report coverage of ALL original fields (kept and discarded).
# Usage: Rscript scripts/05_clean_metadata.R [input.tsv] [output.tsv]

args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 2L) {
  stop("Usage: Rscript scripts/05_clean_metadata.R [input.tsv] [output.tsv]")
}

input <- if (length(args) >= 1L) args[1] else
  "data/processed/Serratia_merged_metadata.tsv"
output <- if (length(args) >= 2L) {
  args[2]
} else if (grepl("_merged_metadata\\.tsv$", input)) {
  sub("_merged_metadata\\.tsv$", "_clean_metadata.tsv", input)
} else {
  sub("\\.tsv$", "_clean_metadata.tsv", input)
}

if (!grepl("\\.tsv$", output)) stop("Output must be a .tsv file.")
if (identical(input, output)) {
  stop("Output must differ from input (expected *_merged_metadata.tsv).")
}
if (!file.exists(input)) stop("Input not found: ", input)

# Grouped whitelist: deliberately preserve parallel GTDB/BioSample fields.
fields <- list(
  identifiers = c(
    "BioSample", "gtdb_accession", "gtdb_ncbi_biosample",
    "gtdb_ncbi_genbank_assembly_accession",
    "gtdb_ncbi_bioproject", "biosample_bioproject"
  ),
  taxonomy = c(
    "gtdb_gtdb_taxonomy", "gtdb_ncbi_taxonomy",
    "gtdb_ncbi_organism_name", "biosample_organism",
    "gtdb_ncbi_species_taxid", "gtdb_ncbi_taxid", "biosample_taxid",
    "gtdb_gtdb_representative"
  ),
  genome_quality = c(
    "gtdb_ncbi_assembly_level", "gtdb_ncbi_assembly_type",
    "gtdb_checkm2_completeness", "gtdb_checkm2_contamination",
    "gtdb_checkm_completeness", "gtdb_checkm_contamination",
    "gtdb_genome_size", "gtdb_gc_percentage", "gtdb_contig_count",
    "gtdb_n50_contigs", "gtdb_scaffold_count", "gtdb_n50_scaffolds"
  ),
  source_and_geography = c(
    "gtdb_ncbi_country", "gtdb_ncbi_date", "gtdb_ncbi_isolate",
    "gtdb_ncbi_isolation_source", "gtdb_ncbi_lat_lon",
    "gtdb_ncbi_strain_identifiers",
    "biosample_sample_name", "biosample_strain", "biosample_isolate",
    "biosample_isolation_source", "biosample_geo_loc_name",
    "biosample_lat_lon", "biosample_collection_date", "biosample_host",
    "biosample_host_disease", "biosample_env_broad_scale",
    "biosample_env_local_scale", "biosample_env_medium",
    "biosample_environment", "biosample_sample_type", "biosample_tissue",
    "biosample_body_site", "biosample_host_tissue_sampled"
  ),
  biosample_context = c(
    "biosample_title", "biosample_description", "biosample_package",
    "biosample_metadata_source", "biosample_has_location",
    "biosample_has_ecology"
  ),
  sequencing = c(
    "sra_run_count", "sra_total_bases", "sra_total_spots",
    "sra_platforms", "sra_library_strategies", "sra_library_layouts"
  )
)
keep <- unlist(fields, use.names = FALSE)
stopifnot(!anyDuplicated(keep))

message("Reading: ", input)
dat <- read.delim(
  input, sep = "\t", header = TRUE, check.names = FALSE,
  colClasses = "character", quote = "\"", comment.char = "",
  fill = FALSE, stringsAsFactors = FALSE
)
if (anyDuplicated(names(dat))) stop("Duplicate column names in input header.")
if (nrow(dat) == 0L) stop("Input contains no rows.")

required <- c("BioSample", "gtdb_accession")
missing_required <- setdiff(required, names(dat))
if (length(missing_required)) {
  stop("Missing required columns: ", paste(missing_required, collapse = ", "))
}

missing_optional <- setdiff(keep, names(dat))
if (length(missing_optional)) {
  warning("Requested columns not found: ", paste(missing_optional, collapse = ", "))
}

# Conservative missing-value normalization. Do not turn "none", "no", or "0" into NA.
missing_tokens <- c(
  "", "na", "n/a", "null", "missing", "unknown",
  "not available", "not provided", "not collected", "not applicable"
)
for (nm in names(dat)) {
  x <- trimws(dat[[nm]])
  x[tolower(x) %in% missing_tokens] <- NA_character_
  dat[[nm]] <- x
}

selected <- keep[keep %in% names(dat)]
clean <- dat[, selected, drop = FALSE]

# Audit every original column so exclusions can be reviewed later.
non_missing <- vapply(dat, function(x) sum(!is.na(x)), integer(1))
unique_values <- vapply(dat, function(x) length(unique(x[!is.na(x)])), integer(1))
report <- data.frame(
  column = names(dat),
  kept = names(dat) %in% selected,
  non_missing = as.integer(non_missing),
  coverage_pct = round(100 * non_missing / nrow(dat), 1),
  unique_non_missing = as.integer(unique_values),
  stringsAsFactors = FALSE
)
report <- report[order(!report$kept, -report$coverage_pct), , drop = FALSE]
report_file <- sub("\\.tsv$", "_column_report.tsv", output)

dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
write.table(clean, output, sep = "\t", quote = FALSE,
            row.names = FALSE, na = "")
write.table(report, report_file, sep = "\t", quote = FALSE,
            row.names = FALSE, na = "")

message("Rows: ", nrow(clean), " | Columns: ", ncol(dat), " -> ", ncol(clean))
message("Repeated GTDB accessions: ", sum(duplicated(clean$gtdb_accession[!is.na(clean$gtdb_accession)])))
if ("gtdb_ncbi_biosample" %in% names(clean)) {
  both <- !is.na(clean$BioSample) & !is.na(clean$gtdb_ncbi_biosample)
  message("BioSample ID disagreements (both present): ",
          sum(clean$BioSample[both] != clean$gtdb_ncbi_biosample[both]),
          " / ", sum(both))
}
message("Cleaned table: ", output)
message("Column report: ", report_file)
