#!/usr/bin/env Rscript

# USAGE
# Rscript scripts/07_species_summary.r \
#    data/processed/Serratia_standardized_metadata.tsv \
#    data/processed/Serratia_species_summary.tsv

# Summarize unique GTDB/NCBI species pairs; one genome belongs to one group.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop('Usage: Rscript scripts/07_species_summary.R INPUT.tsv OUTPUT.tsv', call. = FALSE)
if (!file.exists(args[1])) stop('Input not found: ', args[1])
if (normalizePath(args[1]) == normalizePath(args[2], mustWork = FALSE)) stop('Input and output must differ')
d <- read.delim(args[1], sep='\t', quote='"', comment.char='', check.names=FALSE, colClasses='character', na.strings=NULL)
required <- c('accession','species','ncbi_species','country','geo_location','clinical_status','isolation_source','sample_type','tissue')
if (any(!required %in% names(d))) stop('Missing columns: ', paste(setdiff(required,names(d)), collapse=', '))
clean <- function(x) {
  x[is.na(x)] <- ''
  x <- trimws(gsub('[[:space:]]+', ' ', x))
  x[tolower(x) %in% c('', 'na','n/a','none','null','unknown','not provided','not collected','not applicable','missing','not available','unspecified','-')] <- ''
  x
}
for (nm in required) d[[nm]] <- clean(d[[nm]])
if (!nrow(d)) stop('Empty input')
if (any(!nzchar(d$accession))) stop('Missing accession')
if (anyDuplicated(d$accession)) stop('Duplicate accession: resolve duplicates in step 06 before summarizing')
if (any(!nzchar(d$species))) stop('Missing GTDB species: check step 06')
# Excel-style letters (A..Z, AA..); A is reserved for exact matches.
letter <- function(n) {
  s <- ''
  while (n > 0L) { n <- n-1L; s <- paste0(LETTERS[n %% 26L + 1L], s); n <- n %/% 26L }
  s
}
parts <- function(x) {
  y <- clean(strsplit(x, ' | ', fixed=TRUE)[[1]])
  y[nzchar(y) & !duplicated(tolower(y))]
}
# Reconstruct complete locations for the summary only; accept both old and V2 tables.
geo <- function(country, locality) {
  co <- parts(country); loc <- parts(locality)
  if (!length(co)) return(loc)
  if (!length(loc)) return(co)
  if (length(co) != 1L) return(paste0('Country: ',paste(co,collapse=' / '),'; locality: ',paste(loc,collapse=' / ')))
  vapply(loc,function(z) {
    if (tolower(z)==tolower(co) || startsWith(tolower(z),paste0(tolower(co),':'))) z
    else paste0(co, ': ',z)
  }, character(1), USE.NAMES=FALSE)
}
# Prefer biological source/tissue, fall back to sample type when unavailable.
sample_values <- function(source,tissue,type) {
  z <- c(parts(source),parts(tissue)); z <- z[!duplicated(tolower(z))]
  if (!length(z)) z <- parts(type)
  z
}
# Count each distinct value once per genome; a genome may have multiple values.
frequency <- function(values) {
  values <- lapply(values,function(v) {v <- v[nzchar(v)]; if (!length(v)) 'Unknown' else v[!duplicated(tolower(v))]})
  v <- unlist(values,use.names=FALSE); k <- tolower(v)
  tab <- table(k); labels <- v[match(names(tab),k)]
  ord <- order(-as.integer(tab),tolower(labels),method='radix')
  paste0(labels[ord],' (',as.integer(tab)[ord],')',collapse='; ')
}
g <- Map(geo,d$country,d$geo_location)
s <- Map(sample_values,d$isolation_source,d$tissue,d$sample_type)
cstatus <- lapply(d$clinical_status,function(x) if(nzchar(x)) x else 'Unknown')
result <- list(); j <- 0L
for (sp in sort(unique(d$species),method='radix')) {
  ids <- which(d$species == sp)
  nc <- unique(d$ncbi_species[ids]) # first occurrence order within input
  other <- nc[nc != sp]
  # Missing NCBI is a separate group, like any other unmatched value.
  groups <- c(if (sp %in% nc) sp, other)
  for (ncb in groups) {
    ix <- ids[d$ncbi_species[ids] == ncb]
    suffix <- if (ncb == sp) 'A' else letter(match(ncb,other)+1L)
    j <- j+1L
    result[[j]] <- data.frame(species=paste(sp,suffix), counts=length(ix),
      ncbi_species=if(nzchar(ncb)) ncb else 'Unknown',
      geo_loc=frequency(g[ix]), clinical=frequency(cstatus[ix]), sample=frequency(s[ix]),
      stringsAsFactors=FALSE)
  }
}
out <- do.call(rbind,result)
stopifnot(sum(out$counts)==nrow(d), !anyDuplicated(out$species))
dir.create(dirname(args[2]),recursive=TRUE,showWarnings=FALSE)
write.table(out,args[2],sep='\t',quote=TRUE,row.names=FALSE,na='')
cat('Genomes:',nrow(d),'\nGTDB species:',length(unique(d$species)), '\nSpecies pairs:',nrow(out),'\nSum of counts:',sum(out$counts),'\nOutput:',args[2],'\n')
