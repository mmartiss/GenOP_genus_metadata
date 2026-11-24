library(dplyr)

df <- read.delim(
  "data/processed/Serratia_merged_metadata.tsv",
  header = TRUE
)

filtered <- df %>%
  select(biosample_geo_loc_name, gtdb_ncbi_country)

new_df <- read.delim(
  "data/processed/Serratia_clean_metadata.tsv",
  header = TRUE
)
