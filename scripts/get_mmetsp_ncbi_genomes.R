# Import libraries
library(tidyverse)
theme_set(theme_bw())
library("sf")
library("rnaturalearth")
library("rnaturalearthdata")
library(taxonomizr)

# Read raw data on MMETSP
mmetsp_raw_data <- read_tsv('data/sample-attr.tab.txt')

# Unpack attributes.
mmetsp_wider <- mmetsp_raw_data %>%
  rename(sample_name_main = sample_name) %>%
  pivot_wider(id_cols = c('sample_id', 'sample_name_main'), 
  names_from = "attr_type", 
  values_from = "attr_value", 
  names_repair = "unique")

# For duplicate lat/long, keep first latitude, and second longitude. Keep first sample name. 
# For other columns with duplicate entries, paste together (comma separated).
# keep first latitude, second longitude.
mmetsp_wider$longitude <- lapply(mmetsp_wider$longitude, function(x) if(length(x) > 1) x[[2]] else x)
mmetsp_wider$latitude <- lapply(mmetsp_wider$latitude, function(x) if(length(x) > 1) x[[1]] else x)

# keep first sample name
mmetsp_wider$sample_name <- lapply(mmetsp_wider$sample_name, function(x) if(length(x) > 1) x[[1]] else x)

# paste (comma separated) other duplicate entries.
mmetsp_wider$source_mat_id <- lapply(mmetsp_wider$source_mat_id, function(x) paste(x, collapse = ","))
mmetsp_wider$primary_citation <- lapply(mmetsp_wider$primary_citation, function(x) paste(x, collapse = ","))
mmetsp_wider$other_experimental_metadata_available <- lapply(mmetsp_wider$other_experimental_metadata_available, function(x) paste(x, collapse = ","))
mmetsp_wider$other_environmental_metadata_available <- lapply(mmetsp_wider$other_environmental_metadata_available, function(x) paste(x, collapse = ","))
mmetsp_wider$country <- lapply(mmetsp_wider$country, function(x) paste(x, collapse = ","))
mmetsp_wider$additional_citations <- lapply(mmetsp_wider$additional_citations, function(x) paste(x, collapse = ","))

# Unnest wide data, now that lists are only of length one.
mmetsp_unnested <- mmetsp_wider %>%
  unnest()

# Select and unnest taxon info.
mmetsp_taxon <- mmetsp_wider %>%
  select(sample_id, sample_name_main, taxon_id, phylum, class, order, genus, species, strain, fastq_file, latitude, longitude) %>%
  unnest() %>%
  mutate(
    genus_species_strain = gsub(" ", "_", paste(genus, species, strain, sep = "_"))
  )

# Select only barebones
mmetsp_select <- mmetsp_taxon %>%
  select('sample_id', 'sample_name_main', 'taxon_id', 'genus_species_strain', 'fastq_file')

# Specify column names.
colNames <- "assembly_accession, bioproject, biosample, wgs_master, refseq_category, taxid, species_taxid, organism_name, infraspecific_name, infraspecific_name2, isolateversion_status, assembly_level, release_type, genome_rep, seq_rel_date, asm_name, submitter, gbrs_paired_asm, paired_asm_comp, ftp_path, excluded_from_refseq, relation_to_type_material"
colNamesVec <- unlist(str_split(colNames, ", "))

# Read in metadata for all Genbank genomes.
# genome filename ending
suffix <- "_genomic.fna.gz"

genbank <- read_tsv('../data/assembly_summary_genbank.txt',
                    comment = "#",
                    col_names = colNamesVec) %>%
  # remove entries missing ftp genome path
  drop_na(ftp_path) %>%
  mutate(taxid = as.character(taxid),
         species_taxid = as.character(species_taxid)) %>%
  rowwise() %>%
  mutate(
    # add full filename and ftp path of genbank genomes
    genome_filename = paste(tail(str_split(ftp_path, '/')[[1]], 1), suffix, sep = ""),
    genome_ftp_path = paste(ftp_path, genome_filename, sep = "/")
  )


# Join MMETSP and Genbank data by taxon id.

# join mmetsp & genbank
join_taxid <- inner_join(mmetsp_unnested, genbank, by = c('taxon_id' = 'taxid')) %>%
  # remove duplicates
  distinct(sample_id, .keep_all = TRUE) %>%
  distinct(taxon_id, .keep_all = TRUE) %>%
  rowwise() %>%
  mutate(
    # add full filename and ftp path of genbank genomes
    genome_filename = paste(tail(str_split(ftp_path, '/')[[1]], 1), suffix, sep = ""),
    genome_ftp_path = paste(ftp_path, genome_filename, sep = "/"),
    strain = str_replace(strain, " ", "-"),
    new_fastq_ftp = paste("/iplant/home/shared/imicrobe/projects/104/samples", sample_id, 
                          paste(sample_name_main, ".fastq.tar", sep = ""), 
                          sep = "/"),
    genome_name = str_split(genome_filename, ".fna.gz")[[1]][1]
  ) 

# Write taxid joined 
# CSV of all info:
write_csv(join_taxid, '../data/mmetsp_ncbi_genome_info.csv')