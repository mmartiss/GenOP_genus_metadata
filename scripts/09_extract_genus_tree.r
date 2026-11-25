#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(ape)
    library(data.table)
    library(phangorn)
})

# ==========================================
# Arguments
# ==========================================

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 1) {
    stop(
        "Usage: Rscript 09_extract_genus_tree.R <genus>"
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

project_dir <- dirname(dirname(normalizePath(
    sub("^--file=", "", script_arg[1])
)))

# ==========================================
# Input and output
# ==========================================

metadata_file <- file.path(
    project_dir,
    "data/filtered",
    paste0(tolower(genus), "_genera.tsv.gz")
)

tree_file <- file.path(
    project_dir,
    "data/raw/gtdb/bac120.tree.gz"
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

for (f in c(metadata_file, tree_file)) {

    if (!file.exists(f)) {
        stop("Input file not found: ", f)
    }
}

# ==========================================
# Read metadata
# ==========================================

cat("[INFO] Reading GTDB metadata...\n")

required_cols <- c(
    "accession",
    "gtdb_taxonomy",
    "gtdb_representative",
    "gtdb_genome_representative"
)

metadata <- fread(
    cmd = paste("gzip -dc", shQuote(metadata_file)),
    select = required_cols
)

missing_cols <- setdiff(
    required_cols,
    names(metadata)
)

if (length(missing_cols) > 0) {

    stop(
        "Required columns missing: ",
        paste(missing_cols, collapse = ", ")
    )
}

# ==========================================
# Filter exact GTDB genus
# ==========================================

genus_pattern <- paste0(
    "(^|;)g__", genus, "(;|$)"
)

metadata <- metadata[
    grepl(genus_pattern, gtdb_taxonomy)
]

if (nrow(metadata) == 0) {
    stop("No genomes found for genus: ", genus)
}

# ==========================================
# Normalize accession IDs
# ==========================================

normalize_accession <- function(x) {

    sub(
        "^(RS_|GB_)",
        "",
        trimws(as.character(x))
    )
}

metadata[, accession_clean :=
    normalize_accession(accession)]

metadata[, representative_clean :=
    normalize_accession(gtdb_genome_representative)]

metadata[, is_representative :=
    tolower(as.character(gtdb_representative)) %in%
    c("t", "true", "1")]

# Representatives must map to themselves.

metadata[
    is_representative == TRUE,
    representative_clean := accession_clean
]

# Remove invalid accession IDs.

metadata <- metadata[
    !is.na(accession_clean) &
    accession_clean != ""
]

# Check conflicting representative assignments.

conflicts <- metadata[
    ,
    .(
        n_representatives = uniqueN(
            representative_clean
        )
    ),
    by = accession_clean
][n_representatives > 1]

if (nrow(conflicts) > 0) {
    stop("Conflicting representative assignments detected.")
}

# One row per genome.

metadata <- unique(
    metadata,
    by = "accession_clean"
)

target_ids <- metadata$accession_clean

representative_ids <- metadata[
    is_representative == TRUE,
    accession_clean
]

cat(
    "[INFO] Total genomes:",
    length(target_ids), "\n"
)

cat(
    "[INFO] Expected representatives:",
    length(representative_ids), "\n"
)

# ==========================================
# Read GTDB reference tree
# ==========================================

cat("[INFO] Reading GTDB reference tree...\n")

con <- gzfile(
    tree_file,
    open = "rt"
)

reference_tree <- tryCatch(
    read.tree(file = con),
    finally = close(con)
)

if (!inherits(reference_tree, "phylo")) {
    stop("Expected a single phylogenetic tree.")
}

reference_tree$tip.label <- normalize_accession(
    reference_tree$tip.label
)

if (anyDuplicated(reference_tree$tip.label)) {
    stop("Duplicate normalized GTDB tree labels.")
}

cat(
    "[INFO] Total GTDB reference tree tips:",
    Ntip(reference_tree), "\n"
)

# ==========================================
# Match representative genomes
# ==========================================

matched_representatives <- intersect(
    representative_ids,
    reference_tree$tip.label
)

missing_representatives <- setdiff(
    representative_ids,
    reference_tree$tip.label
)

cat(
    "[INFO] Matched representatives:",
    length(matched_representatives), "\n"
)

cat(
    "[INFO] Missing representatives:",
    length(missing_representatives), "\n"
)

if (length(missing_representatives) > 0) {

    warning(
        "Missing representatives:\n",
        paste(
            missing_representatives,
            collapse = "\n"
        )
    )
}

if (length(matched_representatives) < 2) {
    stop("Fewer than two matching representatives.")
}

# ==========================================
# Extract representative tree
# ==========================================

cat("[INFO] Extracting representative tree...\n")

representative_tree <- keep.tip(
    reference_tree,
    matched_representatives
)

# Remove original internal GTDB node labels.

representative_tree$node.label <- NULL

# Validate branch lengths.

if (is.null(representative_tree$edge.length)) {
    stop("Reference tree has no branch lengths.")
}

if (any(!is.finite(representative_tree$edge.length)) ||
    any(representative_tree$edge.length < 0)) {

    stop("Invalid reference tree branch lengths.")
}

# ==========================================
# Midpoint rooting
# ==========================================

cat("[INFO] Midpoint rooting...\n")

if (is.rooted(representative_tree)) {

    representative_tree <- unroot(
        representative_tree
    )
}

representative_tree <- phangorn::midpoint(
    representative_tree
)

representative_tree <- ladderize(
    representative_tree
)

# ==========================================
# Save representative tree
# ==========================================

representative_tree_file <- file.path(
    out_dir,
    paste0(
        genus,
        "_representatives_tree.nwk"
    )
)

write.tree(
    representative_tree,
    file = representative_tree_file,
    digits = 12
)

# ==========================================
# Map genomes to representatives
# ==========================================

cat("[INFO] Mapping genomes to representatives...\n")

metadata[, representative_in_tree :=
    representative_clean %in%
    representative_tree$tip.label]

metadata[, matching_status := fifelse(
    is.na(representative_clean) |
        representative_clean == "",
    "missing_representative_assignment",
    fifelse(
        !representative_in_tree,
        "representative_not_in_tree",
        fifelse(
            is_representative,
            "representative",
            "assigned_to_representative"
        )
    )
)]

mapped_metadata <- metadata[
    representative_in_tree == TRUE
]

unmapped_metadata <- metadata[
    representative_in_tree == FALSE
]

cat(
    "[INFO] Successfully mapped genomes:",
    nrow(mapped_metadata), "\n"
)

cat(
    "[INFO] Unmapped genomes:",
    nrow(unmapped_metadata), "\n"
)

if (nrow(unmapped_metadata) > 0) {

    warning(
        nrow(unmapped_metadata),
        " genomes cannot be included in the expanded tree."
    )
}

# ==========================================
# Group genomes by representative
# ==========================================

genome_groups <- split(
    mapped_metadata$accession_clean,
    mapped_metadata$representative_clean
)

genome_groups <- lapply(
    genome_groups,
    unique
)

# ==========================================
# Expand representative tree
# ==========================================
#
# Each representative tip is replaced by:
#
#       original branch
#             |
#             +--- representative (0)
#             +--- genome 1       (0)
#             +--- genome 2       (0)
#
# All additional branches have zero length.
# No within-species topology is inferred.
#
# ==========================================

expand_tree <- function(tree, groups) {

    old_n_tips <- Ntip(tree)

    old_n_nodes <- tree$Nnode

    old_labels <- tree$tip.label

    # Arrange genome groups in original tree order.

    ordered_groups <- lapply(
        old_labels,
        function(rep_id) {

            members <- groups[[rep_id]]

            if (is.null(members)) {
                members <- rep_id
            }

            unique(c(rep_id, members))
        }
    )

    new_labels <- unlist(
        ordered_groups,
        use.names = FALSE
    )

    if (anyDuplicated(new_labels)) {
        stop("Duplicate tips in expanded tree.")
    }

    new_n_tips <- length(new_labels)

    # Identify representatives requiring polytomies.

    expanded_tips <- which(
        lengths(ordered_groups) > 1
    )

    n_added_nodes <- length(expanded_tips)

    # Map old internal nodes to new numbering.

    node_map <- integer(
        old_n_tips + old_n_nodes
    )

    old_internal_nodes <- seq.int(
        old_n_tips + 1L,
        old_n_tips + old_n_nodes
    )

    node_map[old_internal_nodes] <-
        new_n_tips + seq_len(old_n_nodes)

    # Assign new internal nodes to expanded tips.

    added_nodes <- setNames(
        new_n_tips +
            old_n_nodes +
            seq_len(n_added_nodes),
        expanded_tips
    )

    extra_edges <- list()

    extra_index <- 0L

    # Replace original tips with single tips
    # or new polytomy nodes.

    for (i in seq_len(old_n_tips)) {

        members <- ordered_groups[[i]]

        member_indices <- match(
            members,
            new_labels
        )

        if (length(members) == 1L) {

            node_map[i] <- member_indices[1]

        } else {

            new_node <- added_nodes[[as.character(i)]]
            
	    node_map[i] <- new_node

            extra_index <- extra_index + 1L

            extra_edges[[extra_index]] <- cbind(
                rep(
                    new_node,
                    length(member_indices)
                ),
                member_indices
            )
        }
    }

    # Remap original edges without changing
    # their evolutionary branch lengths.

    new_edges <- matrix(
        node_map[tree$edge],
        ncol = 2
    )

    new_lengths <- tree$edge.length

    # Add zero-length genome branches.

    if (length(extra_edges) > 0) {

        extra_edges <- do.call(
            rbind,
            extra_edges
        )

        new_edges <- rbind(
            new_edges,
            extra_edges
        )

        new_lengths <- c(
            new_lengths,
            rep(0, nrow(extra_edges))
        )
    }

    result <- list(
        edge = new_edges,
        edge.length = new_lengths,
        tip.label = new_labels,
        Nnode = as.integer(
            old_n_nodes + n_added_nodes
        )
    )

    storage.mode(result$edge) <- "integer"

    class(result) <- "phylo"

    result <- reorder.phylo(
        result,
        order = "cladewise"
    )

    return(result)
}

# ==========================================
# Construct expanded tree
# ==========================================

cat("[INFO] Constructing expanded tree...\n")

expanded_tree <- expand_tree(
    representative_tree,
    genome_groups
)

expanded_tree <- ladderize(
    expanded_tree
)

# ==========================================
# Validate expanded tree
# ==========================================

cat("[INFO] Validating expanded tree...\n")

if (anyDuplicated(expanded_tree$tip.label)) {
    stop("Expanded tree contains duplicate tips.")
}

if (!setequal(
    expanded_tree$tip.label,
    mapped_metadata$accession_clean
)) {

    stop("Expanded tree tip IDs do not match mapped genomes.")
}

if (Ntip(expanded_tree) != nrow(mapped_metadata)) {
    stop("Expanded tree genome count mismatch.")
}

if (!is.rooted(expanded_tree)) {
    stop("Expanded tree lost its root.")
}

# Check whether original representative distances
# have been preserved.

original_distances <- cophenetic.phylo(
    representative_tree
)

reduced_tree <- keep.tip(
    expanded_tree,
    representative_tree$tip.label
)

expanded_distances <- cophenetic.phylo(
    reduced_tree
)

expanded_distances <- expanded_distances[
    rownames(original_distances),
    colnames(original_distances),
    drop = FALSE
]

distances_preserved <- isTRUE(
    all.equal(
        original_distances,
        expanded_distances,
        tolerance = 1e-10
    )
)

if (!distances_preserved) {
    stop("Original phylogenetic distances were altered.")
}

cat(
    "[INFO] Original phylogenetic distances preserved:",
    distances_preserved, "\n"
)

# ==========================================
# Save expanded tree
# ==========================================

expanded_tree_file <- file.path(
    out_dir,
    paste0(genus, "_expanded_tree.nwk")
)

write.tree(
    expanded_tree,
    file = expanded_tree_file,
    digits = 12
)

# ==========================================
# Save genome matching report
# ==========================================

metadata[, in_representative_tree :=
    accession_clean %in%
    representative_tree$tip.label]

metadata[, in_expanded_tree :=
    accession_clean %in%
    expanded_tree$tip.label]

fwrite(
    metadata,
    file.path(
        out_dir,
        paste0(genus, "_tree_matching.tsv")
    ),
    sep = "\t"
)

# ==========================================
# Save unmatched genomes
# ==========================================

fwrite(
    unmapped_metadata,
    file.path(
        out_dir,
        paste0(genus, "_tree_unmatched.tsv")
    ),
    sep = "\t"
)

# ==========================================
# Representative group summary
# ==========================================

group_summary <- mapped_metadata[
    ,
    .(
        genome_count = .N,
        representative_count = sum(
            is_representative
        ),
        non_representative_count = sum(
            !is_representative
        ),
        gtdb_species = paste(
            sort(unique(
                sub(
                    "^.*;s__",
                    "",
                    gtdb_taxonomy
                )
            )),
            collapse = "; "
        )
    ),
    by = representative_clean
]

setorder(
    group_summary,
    -genome_count
)

fwrite(
    group_summary,
    file.path(
        out_dir,
        paste0(genus, "_tree_groups.tsv")
    ),
    sep = "\t"
)

# ==========================================
# Tree summary
# ==========================================

tree_summary <- data.table(
    genus = genus,
    total_genomes = length(target_ids),
    expected_representatives = length(
        representative_ids
    ),
    matched_representatives = length(
        matched_representatives
    ),
    missing_representatives = length(
        missing_representatives
    ),
    representative_tree_tips = Ntip(
        representative_tree
    ),
    expanded_tree_tips = Ntip(
        expanded_tree
    ),
    mapped_genomes = nrow(mapped_metadata),
    unmapped_genomes = nrow(unmapped_metadata),
    zero_length_branches_added =
        nrow(expanded_tree$edge) -
        nrow(representative_tree$edge),
    representative_rooting = "midpoint",
    expanded_rooting = "inherited_midpoint",
    phylogenetic_distances_preserved =
        distances_preserved,
    expanded_tree_type =
        "representative_based_polytomies"
)

fwrite(
    tree_summary,
    file.path(
        out_dir,
        paste0(genus, "_tree_summary.tsv")
    ),
    sep = "\t"
)

# ==========================================
# Representative tree PDF
# ==========================================

cat("[INFO] Plotting representative tree...\n")

pdf(
    file.path(out_dir, paste0(genus, "_representatives_tree.pdf")),
    width = 12,
    height = max(6, Ntip(representative_tree) * 0.25)
)

par(mar = c(3, 2, 3, 2))

plot(
    representative_tree,
    type = "phylogram",
    direction = "rightwards",
    show.tip.label = TRUE,
    use.edge.length = TRUE,
    cex = 0.7,
    main = paste0(genus, " representative tree (midpoint rooted)")
)

add.scale.bar()

dev.off()

# ==========================================
# Expanded tree PDF
# ==========================================

cat("[INFO] Plotting expanded tree...\n")

pdf(
    file.path(out_dir, paste0(genus, "_expanded_tree.pdf")),
    width = 16,
    height = 30
)

par(mar = c(3, 2, 3, 2))

plot(
    expanded_tree,
    type = "phylogram",
    direction = "rightwards",
    show.tip.label = FALSE,
    use.edge.length = TRUE,
    edge.width = 0.5,
    main = paste0(genus, " expanded tree (",
                  Ntip(expanded_tree), " genomes, midpoint rooted)")
)

add.scale.bar()

dev.off()

# ==========================================
# Final report
# ==========================================

cat("\n")

cat("========================================\n")
cat("TREE EXTRACTION COMPLETED\n")
cat("========================================\n")

cat(
    "Genus:",
    genus, "\n"
)

cat(
    "Total metadata genomes:",
    length(target_ids), "\n"
)

cat(
    "Representative tree tips:",
    Ntip(representative_tree), "\n"
)

cat(
    "Expanded tree tips:",
    Ntip(expanded_tree), "\n"
)

cat(
    "Unmapped genomes:",
    nrow(unmapped_metadata), "\n"
)

cat(
    "Rooting: midpoint\n"
)

cat(
    "Original distances preserved:",
    distances_preserved, "\n"
)

cat(
    "Output directory:",
    out_dir, "\n"
)

cat("========================================\n")
