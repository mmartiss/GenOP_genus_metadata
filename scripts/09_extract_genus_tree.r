#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(ape)
    library(data.table)
})

# ==========================================
# Arguments
# ==========================================

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 1) {
    stop(
        "Usage: Rscript scripts/09_extract_genus_tree.R <genus>"
    )
}

genus <- args[1]

if (!grepl("^[A-Za-z0-9_]+$", genus)) {
    stop("Invalid genus name.")
}

# ==========================================
# Project directory
# ==========================================

script_arg <- grep(
    "^--file=",
    commandArgs(trailingOnly = FALSE),
    value = TRUE
)

script_path <- normalizePath(
    sub("^--file=", "", script_arg[1])
)

project_dir <- dirname(dirname(script_path))

# ==========================================
# Input and output files
# ==========================================

tree_file <- file.path(
    project_dir,
    "data/raw/gtdb/bac120.tree.gz"
)

metadata_file <- file.path(
    project_dir,
    "data/filtered",
    paste0(tolower(genus), "_genera.tsv.gz")
)

out_dir <- file.path(
    project_dir,
    "data/trees",
    genus
)

dir.create(
    out_dir,
    recursive = TRUE,
    showWarnings = FALSE
)

if (!file.exists(tree_file)) {
    stop("GTDB tree not found: ", tree_file)
}

if (!file.exists(metadata_file)) {
    stop("Metadata file not found: ", metadata_file)
}

# ==========================================
# Read metadata
# ==========================================

cat("[INFO] Reading GTDB metadata...\n")

required_cols <- c(
    "accession",
    "gtdb_taxonomy"
)

metadata <- fread(
    metadata_file,
    select = required_cols
)

if (!all(required_cols %in% names(metadata))) {
    stop(
        "Required columns missing: ",
        paste(
            setdiff(required_cols, names(metadata)),
            collapse = ", "
        )
    )
}

# Exact GTDB genus matching

genus_pattern <- paste0(
    "(^|;)g__", genus, "(;|$)"
)

metadata <- metadata[
    grepl(genus_pattern, gtdb_taxonomy)
]

if (nrow(metadata) == 0) {
    stop("No GTDB genomes found for genus: ", genus)
}

# ==========================================
# Normalize accession IDs
# ==========================================

normalize_accession <- function(x) {
    sub("^(RS_|GB_)", "", x)
}

metadata[, accession_clean :=
    normalize_accession(accession)]

target_ids <- unique(
    metadata$accession_clean[
        !is.na(metadata$accession_clean) &
        metadata$accession_clean != ""
    ]
)

cat(
    "[INFO] Unique GTDB genome IDs:",
    length(target_ids), "\n"
)

# ==========================================
# Read GTDB tree
# ==========================================

cat("[INFO] Reading GTDB bacterial tree...\n")

con <- gzfile(tree_file, open = "rt")

tree <- tryCatch(
    read.tree(file = con),
    finally = close(con)
)

if (!inherits(tree, "phylo")) {
    stop("Expected a single phylogenetic tree.")
}

cat(
    "[INFO] Total GTDB tree tips:",
    Ntip(tree), "\n"
)

# ==========================================
# Match genome IDs
# ==========================================

tree_ids <- normalize_accession(
    tree$tip.label
)

keep_tips <- tree$tip.label[
    tree_ids %in% target_ids
]

matched_ids <- intersect(
    target_ids,
    tree_ids
)

unmatched_ids <- setdiff(
    target_ids,
    tree_ids
)

cat(
    "[INFO] Matched genome IDs:",
    length(matched_ids), "\n"
)

cat(
    "[INFO] Genome IDs absent from tree:",
    length(unmatched_ids), "\n"
)

if (length(keep_tips) < 2) {
    stop("Fewer than two matching tree tips.")
}

# ==========================================
# Extract genus tree
# ==========================================

cat("[INFO] Extracting genus tree...\n")

genus_tree <- keep.tip(
    tree,
    keep_tips
)

# ==========================================
# Save Newick tree
# ==========================================

out_tree <- file.path(
    out_dir,
    paste0(genus, "_tree.nwk")
)

write.tree(
    genus_tree,
    file = out_tree
)

# ==========================================
# Save tree matching report
# ==========================================

report <- data.table(
    accession = target_ids,
    in_tree = target_ids %in% matched_ids
)

fwrite(
    report,
    file.path(
        out_dir,
        paste0(genus, "_tree_matching.tsv")
    ),
    sep = "\t"
)

# ==========================================
# Save tree summary
# ==========================================

summary <- data.table(
    genus = genus,
    metadata_genomes = length(target_ids),
    gtdb_tree_tips = Ntip(tree),
    matched_genomes = length(matched_ids),
    unmatched_genomes = length(unmatched_ids),
    extracted_tree_tips = Ntip(genus_tree)
)

fwrite(
    summary,
    file.path(
        out_dir,
        paste0(genus, "_tree_summary.tsv")
    ),
    sep = "\t"
)

# ==========================================
# Save PDF visualization
# ==========================================

pdf(
    file.path(
        out_dir,
        paste0(genus, "_tree.pdf")
    ),
    width = 12,
    height = max(6, Ntip(genus_tree) * 0.22)
)

plot(
    genus_tree,
    type = "phylogram",
    cex = 0.7,
    no.margin = TRUE
)

add.scale.bar()

dev.off()

# ==========================================
# Finish
# ==========================================

cat("[INFO] Tree extraction completed.\n")
cat("[INFO] Output directory:", out_dir, "\n")
cat("[INFO] Final tree tips:", Ntip(genus_tree), "\n")
