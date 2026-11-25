# GenOP_genus_metadata
Collects genus metadata

### Needed
Rscript (apt install r-base)
ape and data.table (script -e 'install.packages(c("ape", "data.table"), repos="https://cloud.r-project.org")')


## Filter by genus
GTDB genus Serratia	2947 \
NCBI genus Serratia,	224 (But GTDB genus is different) \
zgrep 'g__Serratia;'	3196 \

## SRA
BioSamples processed: 2947 \
BioSamples with SRA: 1401 \
BioSamples without SRA: 1546 \
Errors: 0 \
Unique SRA runs: 1713 

## Metadata merge
GTDB genomes: 2947 \
GTDB BioSample IDs: 2947 \
Matched BioSample records: 1401 \
BioSample IDs with SRA runs: 1401 

## Clean metadata columns coverage (+ what columns removed and whats left)
| column | kept | non_missing | coverage_pct | unique_non_missing |
| --- | --- | ---: | ---: | ---: |
| BioSample | TRUE | 2947 | 100 | 2947 |
| gtdb_accession | TRUE | 2947 | 100 | 2947 |
| gtdb_checkm2_completeness | TRUE | 2947 | 100 | 120 |
| gtdb_checkm2_contamination | TRUE | 2947 | 100 | 215 |
| gtdb_checkm_completeness | TRUE | 2947 | 100 | 251 |
| gtdb_checkm_contamination | TRUE | 2947 | 100 | 193 |
| gtdb_contig_count | TRUE | 2947 | 100 | 422 |
| gtdb_gc_percentage | TRUE | 2947 | 100 | 2930 |
| gtdb_genome_size | TRUE | 2947 | 100 | 2912 |
| gtdb_gtdb_representative | TRUE | 2947 | 100 | 2 |
| gtdb_gtdb_taxonomy | TRUE | 2947 | 100 | 33 |
| gtdb_n50_contigs | TRUE | 2947 | 100 | 2774 |
| gtdb_n50_scaffolds | TRUE | 2947 | 100 | 2779 |
| gtdb_ncbi_assembly_level | TRUE | 2947 | 100 | 4 |
| gtdb_ncbi_bioproject | TRUE | 2947 | 100 | 135 |
| gtdb_ncbi_biosample | TRUE | 2947 | 100 | 2947 |
| gtdb_ncbi_date | TRUE | 2947 | 100 | 603 |
| gtdb_ncbi_genbank_assembly_accession | TRUE | 2947 | 100 | 2947 |
| gtdb_ncbi_organism_name | TRUE | 2947 | 100 | 288 |
| gtdb_ncbi_species_taxid | TRUE | 2947 | 100 | 248 |
| gtdb_ncbi_taxid | TRUE | 2947 | 100 | 290 |
| gtdb_ncbi_taxonomy | TRUE | 2947 | 100 | 26 |
| gtdb_scaffold_count | TRUE | 2947 | 100 | 409 |
| sra_run_count | TRUE | 2947 | 100 | 7 |
| gtdb_ncbi_strain_identifiers | TRUE | 2945 | 99.9 | 2852 |
| biosample_organism | TRUE | 1401 | 47.5 | 66 |
| biosample_taxid | TRUE | 1401 | 47.5 | 66 |
| biosample_metadata_source | TRUE | 1401 | 47.5 | 2 |
| biosample_has_location | TRUE | 1401 | 47.5 | 2 |
| biosample_has_ecology | TRUE | 1401 | 47.5 | 2 |
| sra_total_bases | TRUE | 1401 | 47.5 | 1399 |
| sra_total_spots | TRUE | 1401 | 47.5 | 1346 |
| sra_platforms | TRUE | 1401 | 47.5 | 10 |
| sra_library_strategies | TRUE | 1401 | 47.5 | 6 |
| sra_library_layouts | TRUE | 1401 | 47.5 | 3 |
| biosample_title | TRUE | 1397 | 47.4 | 574 |
| gtdb_ncbi_country | TRUE | 1314 | 44.6 | 267 |
| gtdb_ncbi_isolation_source | TRUE | 1225 | 41.6 | 311 |
| biosample_strain | TRUE | 1148 | 39 | 1109 |
| biosample_geo_loc_name | TRUE | 1014 | 34.4 | 194 |
| biosample_collection_date | TRUE | 1014 | 34.4 | 356 |
| biosample_host | TRUE | 914 | 31 | 58 |
| biosample_isolation_source | TRUE | 905 | 30.7 | 205 |
| biosample_package | TRUE | 818 | 27.8 | 12 |
| biosample_bioproject | TRUE | 661 | 22.4 | 156 |
| gtdb_ncbi_isolate | TRUE | 419 | 14.2 | 412 |
| gtdb_ncbi_lat_lon | TRUE | 412 | 14 | 113 |
| biosample_lat_lon | TRUE | 374 | 12.7 | 85 |
| biosample_isolate | TRUE | 309 | 10.5 | 302 |
| biosample_sample_type | TRUE | 305 | 10.3 | 40 |
| biosample_description | TRUE | 235 | 8 | 99 |
| biosample_host_disease | TRUE | 219 | 7.4 | 39 |
| biosample_env_broad_scale | TRUE | 59 | 2 | 21 |
| biosample_host_tissue_sampled | TRUE | 48 | 1.6 | 15 |
| biosample_env_local_scale | TRUE | 20 | 0.7 | 16 |
| biosample_env_medium | TRUE | 20 | 0.7 | 12 |
| biosample_sample_name | TRUE | 14 | 0.5 | 14 |
| biosample_environment | TRUE | 4 | 0.1 | 2 |
| gtdb_ncbi_assembly_type | TRUE | 0 | 0 | 0 |
| biosample_tissue | TRUE | 0 | 0 | 0 |
| biosample_body_site | TRUE | 0 | 0 | 0 |
| gtdb_ambiguous_bases | FALSE | 2947 | 100 | 125 |
| gtdb_checkm2_model | FALSE | 2947 | 100 | 2 |
| gtdb_checkm_marker_count | FALSE | 2947 | 100 | 3 |
| gtdb_checkm_marker_lineage | FALSE | 2947 | 100 | 2 |
| gtdb_checkm_marker_set_count | FALSE | 2947 | 100 | 3 |
| gtdb_checkm_strain_heterogeneity | FALSE | 2947 | 100 | 63 |
| gtdb_coding_bases | FALSE | 2947 | 100 | 2927 |
| gtdb_coding_density | FALSE | 2947 | 100 | 2936 |
| gtdb_gc_count | FALSE | 2947 | 100 | 2916 |
| gtdb_gtdb_genome_representative | FALSE | 2947 | 100 | 33 |
| gtdb_gtdb_type_designation_ncbi_taxa | FALSE | 2947 | 100 | 2 |
| gtdb_gtdb_type_species_of_genus | FALSE | 2947 | 100 | 2 |
| gtdb_l50_contigs | FALSE | 2947 | 100 | 109 |
| gtdb_l50_scaffolds | FALSE | 2947 | 100 | 104 |
| gtdb_longest_contig | FALSE | 2947 | 100 | 2794 |
| gtdb_longest_scaffold | FALSE | 2947 | 100 | 2812 |
| gtdb_lsu_23s_count | FALSE | 2947 | 100 | 12 |
| gtdb_lsu_5s_count | FALSE | 2947 | 100 | 15 |
| gtdb_mean_contig_length | FALSE | 2947 | 100 | 2904 |
| gtdb_mean_scaffold_length | FALSE | 2947 | 100 | 2904 |
| gtdb_mimag_high_quality | FALSE | 2947 | 100 | 2 |
| gtdb_mimag_low_quality | FALSE | 2947 | 100 | 1 |
| gtdb_mimag_medium_quality | FALSE | 2947 | 100 | 2 |
| gtdb_ncbi_genome_representation | FALSE | 2947 | 100 | 2 |
| gtdb_ncbi_molecule_count | FALSE | 2947 | 100 | 11 |
| gtdb_ncbi_seq_rel_date | FALSE | 2947 | 100 | 603 |
| gtdb_ncbi_spanned_gaps | FALSE | 2947 | 100 | 99 |
| gtdb_ncbi_taxonomy_unfiltered | FALSE | 2947 | 100 | 288 |
| gtdb_ncbi_total_gap_length | FALSE | 2947 | 100 | 674 |
| gtdb_ncbi_total_length | FALSE | 2947 | 100 | 2912 |
| gtdb_ncbi_ungapped_length | FALSE | 2947 | 100 | 2915 |
| gtdb_ncbi_unspanned_gaps | FALSE | 2947 | 100 | 1 |
| gtdb_protein_count | FALSE | 2947 | 100 | 994 |
| gtdb_ssu_count | FALSE | 2947 | 100 | 13 |
| gtdb_total_gap_length | FALSE | 2947 | 100 | 674 |
| gtdb_trna_aa_count | FALSE | 2947 | 100 | 11 |
| gtdb_trna_count | FALSE | 2947 | 100 | 77 |
| gtdb_trna_selenocysteine_count | FALSE | 2947 | 100 | 3 |
| gtdb_ncbi_submitter | FALSE | 2945 | 99.9 | 473 |
| gtdb_ncbi_assembly_name | FALSE | 2917 | 99 | 2917 |
| gtdb_ncbi_ncrna_count | FALSE | 2866 | 97.3 | 25 |
| gtdb_ncbi_rrna_count | FALSE | 2866 | 97.3 | 42 |
| gtdb_ncbi_ssu_count | FALSE | 2866 | 97.3 | 18 |
| gtdb_ncbi_translation_table | FALSE | 2866 | 97.3 | 1 |
| gtdb_ncbi_trna_count | FALSE | 2866 | 97.3 | 77 |
| gtdb_lsu_23s_contig_len | FALSE | 2784 | 94.5 | 1750 |
| gtdb_lsu_23s_length | FALSE | 2784 | 94.5 | 369 |
| gtdb_lsu_23s_query_id | FALSE | 2784 | 94.5 | 2784 |
| gtdb_lsu_5s_contig_len | FALSE | 2784 | 94.5 | 2457 |
| gtdb_lsu_5s_length | FALSE | 2784 | 94.5 | 21 |
| gtdb_lsu_5s_query_id | FALSE | 2784 | 94.5 | 2784 |
| gtdb_lsu_silva_23s_blast_align_len | FALSE | 2784 | 94.5 | 351 |
| gtdb_lsu_silva_23s_blast_bitscore | FALSE | 2784 | 94.5 | 455 |
| gtdb_lsu_silva_23s_blast_evalue | FALSE | 2784 | 94.5 | 1 |
| gtdb_lsu_silva_23s_blast_perc_identity | FALSE | 2784 | 94.5 | 179 |
| gtdb_lsu_silva_23s_blast_subject_id | FALSE | 2784 | 94.5 | 377 |
| gtdb_lsu_silva_23s_taxonomy | FALSE | 2784 | 94.5 | 79 |
| gtdb_ssu_contig_len | FALSE | 2760 | 93.7 | 1488 |
| gtdb_ssu_length | FALSE | 2760 | 93.7 | 215 |
| gtdb_ssu_query_id | FALSE | 2760 | 93.7 | 2760 |
| gtdb_ssu_silva_blast_align_len | FALSE | 2760 | 93.7 | 203 |
| gtdb_ssu_silva_blast_bitscore | FALSE | 2760 | 93.7 | 270 |
| gtdb_ssu_silva_blast_evalue | FALSE | 2760 | 93.7 | 2 |
| gtdb_ssu_silva_blast_perc_identity | FALSE | 2760 | 93.7 | 145 |
| gtdb_ssu_silva_blast_subject_id | FALSE | 2760 | 93.7 | 142 |
| gtdb_ssu_silva_taxonomy | FALSE | 2760 | 93.7 | 58 |
| gtdb_ncbi_contig_count | FALSE | 2623 | 89 | 422 |
| gtdb_ncbi_contig_n50 | FALSE | 2623 | 89 | 2470 |
| gtdb_ncbi_wgs_master | FALSE | 2589 | 87.9 | 2589 |
| biosample_submission_date | FALSE | 1401 | 47.5 | 1315 |
| biosample_publication_date | FALSE | 1401 | 47.5 | 353 |
| biosample_last_update | FALSE | 1401 | 47.5 | 1396 |
| gtdb_ncbi_scaffold_count | FALSE | 1233 | 41.8 | 251 |
| gtdb_ncbi_scaffold_l50 | FALSE | 1233 | 41.8 | 63 |
| gtdb_ncbi_scaffold_n50 | FALSE | 1233 | 41.8 | 1178 |
| gtdb_ncbi_scaffold_n75 | FALSE | 1233 | 41.8 | 1173 |
| gtdb_ncbi_scaffold_n90 | FALSE | 1233 | 41.8 | 1152 |
| gtdb_ncbi_protein_count | FALSE | 453 | 15.4 | 319 |
| gtdb_ssu_gg_blast_align_len | FALSE | 354 | 12 | 43 |
| gtdb_ssu_gg_blast_bitscore | FALSE | 354 | 12 | 75 |
| gtdb_ssu_gg_blast_evalue | FALSE | 354 | 12 | 3 |
| gtdb_ssu_gg_blast_perc_identity | FALSE | 354 | 12 | 52 |
| gtdb_ssu_gg_blast_subject_id | FALSE | 354 | 12 | 13 |
| gtdb_ssu_gg_taxonomy | FALSE | 354 | 12 | 8 |
| gtdb_ncbi_genome_category | FALSE | 107 | 3.6 | 1 |
| gtdb_gtdb_type_designation_ncbi_taxa_sources | FALSE | 59 | 2 | 1 |
| gtdb_ncbi_type_material_designation | FALSE | 32 | 1.1 | 2 |
| gtdb_ncbi_refseq_category | FALSE | 17 | 0.6 | 1 |

## Summary example
| Species | Counts | NCBI species | Geo_loc | Clinical | Sample |
| --- | ---: | --- | --- | --- | --- |
| S. marcescens A | 3 | S. marcescens | USA (2); China (1) | Yes (3) | Blood (2); Urine (1) |
| S. marcescens B | 2 | S. ficaria | Japan (1); China (1) | No (2) | Water (1); Soil (1) |
| S. marcescens C | 2 | S. nevei | USA (1); Germany (1) | Yes (2) | Sputum (1); Blood (1) |
| S. marcescens D | 1 | Unknown | Spain (1) | Unknown (1) | Soil (1) |
| S. nevei A | 1 | S. nevei | Germany (1) | No (1) | Soil (1) |
| S. nevei B | 2 | S. marcescens | France (1); Italy (1) | Yes (1); Unknown (1) | Blood (1); Unknown (1) |
| S. ficaria A | 1 | S. ficaria | France (1) | No (1) | Plant (1) |
