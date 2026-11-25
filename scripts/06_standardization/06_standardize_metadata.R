#!/usr/bin/env Rscript
# Usage: Rscript 06_standardize_metadata.R input.tsv output.tsv [mapping_dir]
args <- commandArgs(trailingOnly=TRUE)
if(length(args)<2) stop('Usage: Rscript 06_standardize_metadata.R input.tsv output.tsv [mapping_dir]')
script <- sub('^--file=', '', grep('^--file=',commandArgs(),value=TRUE)[1])
mapdir <- if(length(args)>=3) args[3] else file.path(dirname(normalizePath(script)), 'mappings')
readtsv <- function(p) read.delim(p, sep='\t', quote='"', comment.char='', check.names=FALSE, colClasses='character', na.strings=NULL)
writeout <- function(x,p) write.table(x,p,sep='\t',row.names=FALSE,quote=TRUE,na='')
d <- readtsv(args[1]); n <- nrow(d)
get <- function(c) if(c %in% names(d)) trimws(d[[c]]) else rep('',n)
clean <- function(x) {x[is.na(x)] <- ''; k <- tolower(trimws(x)); x[k %in% c('','na','n/a','none','unknown','missing','not provided','not applicable','not available','other','unspecified','null','not collected') | grepl('^not (provided|available|collected)\\b',k,perl=TRUE)] <- ''; x}
fields <- list(geography=c('biosample_geo_loc_name','gtdb_ncbi_country'),host='biosample_host',isolation_source=c('biosample_isolation_source','gtdb_ncbi_isolation_source'),host_disease='biosample_host_disease',sample_type='biosample_sample_type',environment=c('biosample_env_broad_scale','biosample_env_local_scale','biosample_env_medium','biosample_environment'),tissue='biosample_host_tissue_sampled')
audit <- list(); clinical <- list(); ecology <- list()
combine <- function(m) apply(m,1,function(v) paste(unique(v[!is.na(v)&v!='']),collapse=' | '))
for(f in names(fields)) {
 m <- readtsv(file.path(mapdir,paste0(f,'_mapping.tsv')))
 vals <- countries <- list()
 for(c in fields[[f]]) {
  raw <- get(c); idx <- match(raw,trimws(m$raw)); mapped <- !is.na(idx); v <- clean(raw); v[mapped] <- m$standard[idx[mapped]]; v[is.na(v)] <- ''
  status <- ifelse(raw=='','missing',ifelse(mapped,m$status[idx],'unmapped')); status[clean(raw)==''] <- 'missing'
  d[[paste0('std_',c)]] <- v; vals[[c]] <- v
  co <- m$country[idx]; co[is.na(co)] <- ''; countries[[c]] <- co
  ec <- m$ecology[idx]; ec[is.na(ec)] <- ''; ecology[[c]] <- ec
  cl <- m$clinical[idx]; cl[is.na(cl)] <- ''; clinical[[c]] <- cl
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
prefix <- sub('\\.tsv$','',args[2],ignore.case=TRUE); dir.create(dirname(args[2]),recursive=TRUE,showWarnings=FALSE); writeout(d,args[2]); report <- do.call(rbind,audit); writeout(report,paste0(prefix,'_audit.tsv')); writeout(report[report$status %in% c('review','unmapped'),],paste0(prefix,'_review.tsv')); writeout(d[d$std_country_conflict==1 | d$std_coordinate_conflict==1 | d$std_isolation_source_disagreement==1 | d$std_taxonomy_disagreement==1 | d$std_clinical_status=='Conflict' | d$std_date_precision=='Review',c(intersect(c('BioSample','gtdb_accession'),names(d)),grep('^std_',names(d),value=TRUE))],paste0(prefix,'_conflicts.tsv')); cat('Rows:',n,'\nOutput:',args[2],'\n')
