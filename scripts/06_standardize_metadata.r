#!/usr/bin/env Rscript
# Usage: Rscript 06_standardize_metadata.R input.tsv output.tsv [mapping_dir]
# Example
# Rscript scripts/06_standardize_metadata.R     data/processed/Serratia_merged_metadata.tsv     data/processed/Serratia_standardized_metadata.tsv
args <- commandArgs(trailingOnly=TRUE)
if(length(args)<2 || length(args)>3) stop('Usage: Rscript 06_standardize_metadata.R input.tsv output.tsv [mapping_dir]')
script <- sub('^--file=', '', grep('^--file=',commandArgs(),value=TRUE)[1])
mapdir <- if(length(args)>=3) args[3] else file.path(dirname(dirname(normalizePath(script))), 'config', 'metadata')
readtsv <- function(p) read.delim(p, sep='\t', quote='"', comment.char='', check.names=FALSE, colClasses='character', na.strings=NULL)
writeout <- function(x,p) write.table(x,p,sep='\t',row.names=FALSE,quote=TRUE,na='')
d <- readtsv(args[1]); n <- nrow(d)
if(n == 0) stop('Input has no rows')
required_input <- c('BioSample','gtdb_accession','gtdb_gtdb_taxonomy')
if(any(!required_input %in% names(d))) stop('Missing input columns: ', paste(setdiff(required_input,names(d)),collapse=', '))
original_cols <- names(d)
if(anyDuplicated(d$gtdb_accession) || any(is.na(d$gtdb_accession) | trimws(d$gtdb_accession)=='')) stop('Genome accession is missing or duplicated; resolve merge before standardization')
get <- function(c) if(c %in% names(d)) ifelse(is.na(d[[c]]),'',trimws(d[[c]])) else rep('',n)
clean <- function(x) {x[is.na(x)] <- ''; k <- tolower(trimws(x)); x[k %in% c('','na','n/a','none','unknown','missing','not provided','not applicable','not available','other','unspecified','null','not collected') | grepl('^not (provided|available|collected)\\b',k,perl=TRUE)] <- ''; x}
fields <- list(geography=c('biosample_geo_loc_name','gtdb_ncbi_country'),host='biosample_host',isolation_source=c('biosample_isolation_source','gtdb_ncbi_isolation_source'),host_disease='biosample_host_disease',sample_type='biosample_sample_type',environment=c('biosample_env_broad_scale','biosample_env_local_scale','biosample_env_medium','biosample_environment'),tissue='biosample_host_tissue_sampled')
audit <- list(); clinical <- list(); ecology <- list()
combine <- function(m) apply(m,1,function(v) paste(unique(v[!is.na(v)&v!='']),collapse=' | '))
for(f in names(fields)) {
 m <- readtsv(file.path(mapdir,paste0(f,'_mapping.tsv')))
 vals <- countries <- list()
 for(c in fields[[f]]) {
  raw <- get(c); idx <- match(raw,trimws(m$raw)); mapped <- !is.na(idx); v <- clean(raw); v[mapped] <- m$standard[idx[mapped]]; v[is.na(v) | clean(raw)==''] <- ''
  status <- ifelse(raw=='','missing',ifelse(mapped,m$status[idx],'unmapped')); status[clean(raw)==''] <- 'missing'
  d[[paste0('std_',c)]] <- v; vals[[c]] <- v
  co <- m$country[idx]; co[is.na(co) | clean(raw)==''] <- ''; countries[[c]] <- co
  ec <- m$ecology[idx]; ec[is.na(ec) | clean(raw)==''] <- ''; ecology[[c]] <- ec
  cl <- m$clinical[idx]; cl[is.na(cl) | clean(raw)==''] <- ''; clinical[[c]] <- cl
  a <- data.frame(column=c,raw=raw,standard=v,status=status,stringsAsFactors=FALSE); audit[[c]] <- aggregate(list(count=rep(1,n)),a,sum)
 }
 d[[paste0('std_',f)]] <- combine(do.call(cbind,vals))
 if(f=='geography') {
  d$std_country <- combine(do.call(cbind,countries)); d$std_country_conflict <- as.integer(grepl(' | ',d$std_country,fixed=TRUE));
  a <- vals[[1]]; b <- vals[[2]]; same <- countries[[1]]!='' & countries[[1]]==countries[[2]]; takeb <- same & a==countries[[1]] & b!=countries[[2]]; takea <- same & b==countries[[2]] & a!=countries[[1]]; d$std_geography[takeb] <- b[takeb]; d$std_geography[takea] <- a[takea]
 }
 if(f=='isolation_source') d$std_isolation_source_disagreement <- as.integer(vals[[1]]!='' & vals[[2]]!='' & tolower(vals[[1]])!=tolower(vals[[2]]))
}
cl <- combine(do.call(cbind,clinical)); d$std_clinical_status <- ifelse(cl=='','Unknown',ifelse(grepl(' | ',cl,fixed=TRUE),'Conflict',cl))
ec <- combine(do.call(cbind,ecology)); d$std_ecology <- ifelse(ec=='','Unknown',ec)
parsecoord <- function(s) {
 s <- trimws(s); r <- regmatches(s,regexec('^([0-9]+(?:\\.[0-9]+)?)\\s*([NS])\\s+([0-9]+(?:\\.[0-9]+)?)\\s*([EW])$',toupper(s),perl=TRUE))[[1]]
 if(length(r)) z <- c(as.numeric(r[2])*ifelse(r[3]=='S',-1,1),as.numeric(r[4])*ifelse(r[5]=='W',-1,1)) else {
 r <- regmatches(s,regexec('^([+-]?[0-9]+(?:\\.[0-9]+)?)\\s*[,/]\\s*([+-]?[0-9]+(?:\\.[0-9]+)?)$',s,perl=TRUE))[[1]]; z <- if(length(r)) as.numeric(r[2:3]) else c(NA_real_,NA_real_) }
 if(anyNA(z)||abs(z[1])>90||abs(z[2])>180) c(NA_real_,NA_real_) else z
}
a <- t(vapply(get('biosample_lat_lon'),parsecoord,numeric(2))); b <- t(vapply(get('gtdb_ncbi_lat_lon'),parsecoord,numeric(2))); both <- !is.na(a[,1]) & !is.na(b[,1]); conflict <- both & (abs(a[,1]-b[,1])>0.01 | abs(a[,2]-b[,2])>0.01); useb <- is.na(a[,1]); a[useb,] <- b[useb,]; a[conflict,] <- NA_real_; d$std_latitude <- a[,1]; d$std_longitude <- a[,2]; d$std_coordinate_conflict <- as.integer(conflict)
date <- clean(get('biosample_collection_date')); pat <- '^[0-9]{4}(-[0-9]{2}(-[0-9]{2})?)?$'; valid <- grepl(pat,date); range <- grepl('^[0-9]{4}(-[0-9]{2}(-[0-9]{2})?)?/[0-9]{4}(-[0-9]{2}(-[0-9]{2})?)?$',date); d$std_collection_date <- ifelse(valid|range,date,''); d$std_collection_year_start <- ifelse(valid|range,substr(date,1,4),''); d$std_collection_year_end <- ifelse(range,substr(sub('^.*/','',date),1,4),d$std_collection_year_start); d$std_date_precision <- ifelse(date=='','Unknown',ifelse(range,'Range',ifelse(!valid,'Review',ifelse(nchar(date)==4,'Year',ifelse(nchar(date)==7,'Month','Day')))))
sp <- function(x) vapply(strsplit(x,';',fixed=TRUE),function(v) {v <- trimws(v); z <- v[startsWith(v,'s__')]; if(length(z)) sub('^s__','',z[1]) else ''},character(1))
d$std_species <- sp(get('gtdb_gtdb_taxonomy')); d$std_ncbi_species <- sp(get('gtdb_ncbi_taxonomy')); d$std_taxonomy_disagreement <- as.integer(d$std_species!='' & d$std_ncbi_species!='' & d$std_species!=d$std_ncbi_species)

# Final analysis column selection and clinical refinement
cols <- c(BioSample='BioSample', accession='gtdb_accession', species='std_species', ncbi_species='std_ncbi_species', taxonomy_disagreement='std_taxonomy_disagreement', bioproject='gtdb_ncbi_bioproject', representative='gtdb_gtdb_genome_representative', country='std_country', geo_location='std_geography', latitude='std_latitude', longitude='std_longitude', collection_date='std_collection_date', collection_year_start='std_collection_year_start', collection_year_end='std_collection_year_end', date_precision='std_date_precision', host='std_host', host_disease='std_host_disease', isolation_source='std_isolation_source', ecology='std_ecology', clinical_status='std_clinical_status', environment='std_environment', sample_type='std_sample_type', tissue='std_tissue', country_conflict='std_country_conflict', coordinate_conflict='std_coordinate_conflict', isolation_source_disagreement='std_isolation_source_disagreement', genome_size='gtdb_genome_size', gc_percentage='gtdb_gc_percentage', completeness='gtdb_checkm2_completeness', contamination='gtdb_checkm2_contamination', assembly_level='gtdb_ncbi_assembly_level')
required <- c('BioSample','gtdb_accession','std_species','std_clinical_status','std_isolation_source')
if(any(!required %in% names(d))) stop('Missing required columns: ',paste(setdiff(required,names(d)),collapse=', '))
for(c in setdiff(unname(cols),names(d))) d[[c]] <- rep(NA_character_,n)
out <- d[unname(cols)]; names(out) <- names(cols)
# Keep country and subnational locality separate; retain unparsed values.
for(i in seq_len(nrow(out))) {
  co <- out$country[i]; loc <- out$geo_location[i]
  if(is.na(loc) || loc=='') next
  if(is.na(co) || co=='') next
  countries <- trimws(strsplit(co,' | ',fixed=TRUE)[[1]])
  parts <- trimws(strsplit(loc,' | ',fixed=TRUE)[[1]])
  for(j in seq_along(parts)) for(k in countries) {
    if(tolower(parts[j])==tolower(k)) parts[j] <- ''
    else if(startsWith(tolower(parts[j]),paste0(tolower(k),':'))) parts[j] <- trimws(substring(parts[j],nchar(k)+2L))
  }
  out$geo_location[i] <- paste(unique(parts[nzchar(parts)]),collapse=' | ')
}
# Mark likely clinical specimens separately; do not infer confirmed patient origin.
s <- tolower(trimws(out$isolation_source)); s[is.na(s)] <- ''
probable <- grepl('^(blood|blood culture|sputum|urine|wound swab|endotracheal aspirate|tracheal aspirate|bronchoalveolar lavage|transtracheal aspiration)$',s)
out$clinical_evidence <- ifelse(out$clinical_status=='Clinical','Explicit clinical classification',ifelse(out$clinical_status=='Non-clinical','Existing non-clinical classification',ifelse(out$clinical_status=='Unknown' & probable,'Clinical specimen type; origin unconfirmed','Insufficient or conflicting evidence')))
out$clinical_status[out$clinical_status=='Unknown' & probable] <- 'Probable clinical'
# Collapse punctuation-only duplicates in combined isolation sources.
key <- function(s) gsub('[^[:alnum:]]','',tolower(s))
for(i in which(!is.na(out$isolation_source))) { v <- strsplit(out$isolation_source[i],' | ',fixed=TRUE)[[1]]; out$isolation_source[i] <- paste(v[!duplicated(key(v))],collapse=' | '); if(length(unique(key(v)))==1 && length(v)>1 && 'isolation_source_disagreement' %in% names(out)) out$isolation_source_disagreement[i] <- '0' }

# Write only the compact analysis table; preserve original input untouched.
prefix <- sub('\\.tsv$', '', args[2], ignore.case=TRUE)
if(normalizePath(args[1],mustWork=TRUE)==normalizePath(args[2],mustWork=FALSE)) stop('Input and output must differ')
dir.create(dirname(args[2]), recursive=TRUE, showWarnings=FALSE)
audit_dir <- file.path(dirname(args[2]), 'metadata_audit')
dir.create(audit_dir, recursive=TRUE, showWarnings=FALSE)
audit_prefix <- file.path(audit_dir, basename(prefix))
report <- do.call(rbind,audit)
writeout(report,paste0(audit_prefix,'_audit.tsv'))
writeout(report[report$status %in% c('review','unmapped'),],paste0(audit_prefix,'_review.tsv'))
conflict <- out$country_conflict==1 | out$coordinate_conflict==1 | out$isolation_source_disagreement==1 | out$clinical_status=='Conflict' | out$date_precision=='Review'
conflict[is.na(conflict)] <- FALSE
raw_cols <- intersect(unique(c(unlist(fields), 'biosample_lat_lon','gtdb_ncbi_lat_lon','biosample_collection_date')), original_cols)
writeout(cbind(out[conflict,,drop=FALSE],d[conflict,raw_cols,drop=FALSE]),paste0(audit_prefix,'_conflicts.tsv'))
tax <- which(out$taxonomy_disagreement==1)
writeout(out[tax,c('BioSample','accession','species','ncbi_species'),drop=FALSE],paste0(audit_prefix,'_taxonomy_disagreements.tsv'))
writeout(data.frame(column=original_cols,action=ifelse(original_cols %in% unname(cols),'retained_renamed',ifelse(original_cols %in% c(raw_cols,'gtdb_gtdb_taxonomy','gtdb_ncbi_taxonomy'),'used_for_standardization','excluded'))),paste0(audit_prefix,'_column_report.tsv'))
stopifnot(nrow(out)==n,ncol(out)==32)
writeout(out,args[2])
print(table(out$clinical_status,useNA='ifany'))
cat('Rows:',nrow(out),'Columns:',ncol(out),'\nOutput:',args[2],'\nAudit:',audit_dir,'\n')
