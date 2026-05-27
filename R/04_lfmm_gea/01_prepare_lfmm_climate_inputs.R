
rm(list = ls())
gc()

# -----------------------------------------------------------------------------
# Paths
# -----------------------------------------------------------------------------

oak_dat_dir <- "/path/to_oak/dat"
beech_dat_dir <- "/path/to_beech/dat"
climate_dir <- "/path/to_climate_data"

# -----------------------------------------------------------------------------
# Helpers
# -----------------------------------------------------------------------------

recode_oak_provenance <- function(provenance) {
  provenance <- gsub("FR", "FRA_", provenance)
  provenance <- gsub("DE", "GER_", provenance)
  provenance <- gsub("PL", "POL_", provenance)
  provenance <- gsub("UK", "GBR_", provenance)
  provenance <- gsub("DK", "DEN_", provenance)

  provenance
}

# -----------------------------------------------------------------------------
# Oak LFMM climate data
# -----------------------------------------------------------------------------

setwd(oak_dat_dir)

oak <- read.csv("oak_dbh.csv")

oak_sample_info <- read.table(
  "common_ring_sample_info_no_reps.csv",
  sep = ",",
  header = TRUE
)

oak_vcf_samples <- read.table(
  "samples_populations_no_no277_reps_missing_removed_only_petrea_minDP3_meanDP10_maxDP80_NA0.90.txt",
  sep = "\t",
  header = FALSE
)
colnames(oak_vcf_samples)[1] <- "Sample"

oak_vcf_samples$RowNumber <- as.numeric(row.names(oak_vcf_samples))

oak_sample_info_reduced <- merge(
  oak_sample_info,
  oak_vcf_samples,
  by.x = "SampleID",
  by.y = "Sample",
  sort = FALSE
)
oak_sample_info_reduced <- oak_sample_info_reduced[order(oak_sample_info_reduced$RowNumber), ]
colnames(oak_sample_info_reduced)[2] <- "Treeid"

oak$Prov <- recode_oak_provenance(oak$Prov)

oak$merging_col <- paste(oak$Site, oak$Block, oak$Prov, oak$Tree, sep = "_")
oak_sample_info_reduced$merging_col <- paste(
  oak_sample_info_reduced$Site,
  oak_sample_info_reduced$Block,
  oak_sample_info_reduced$Prov3,
  oak_sample_info_reduced$Tree,
  sep = "_"
)

oak_merged <- merge(oak_sample_info_reduced, oak, by = "merging_col")

oak_climate <- read.csv(file.path(climate_dir, "oak_clim_for_lfmm_081224.csv"))

oak_sample_ids <- oak_merged[, c("SampleID", "Prov3")]

oak_lfmm_climate <- merge(
  oak_sample_ids,
  oak_climate,
  by.x = "Prov3",
  by.y = "prov"
)

write.csv(
  oak_lfmm_climate,
  file = file.path(climate_dir, "prepped_oak_lfmm_climate_data_081224.csv")
)

# -----------------------------------------------------------------------------
# Beech LFMM climate data
# -----------------------------------------------------------------------------

setwd(beech_dat_dir)

beech_dbh <- read.csv("beech_dbh.csv")
beech_tree_info <- beech_dbh[, c("Treeid", "Site", "Block", "Prov", "Tree")]

beech_sample_info <- read.table(
  file.path(
    beech_dat_dir,
    "samples_populations_fagus_no_high_miss_ind_and_incorr_removed_HWE0_0001_minDP3_meanDP10_maxDP80_NA0.90_reduced_thin1kb.MAF0.05.txt"
  ),
  sep = ",",
  header = FALSE,
  col.names = "SampleID"
)

beech_sample_info$RowNumber <- as.numeric(row.names(beech_sample_info))

split_ids <- strsplit(as.character(beech_sample_info$SampleID), "-")
beech_sample_info <- cbind(beech_sample_info, do.call(rbind, split_ids))

colnames(beech_sample_info)[1] <- "SampleID_vcf"
colnames(beech_sample_info)[2] <- "RowNumber"
colnames(beech_sample_info)[3] <- "sample_nr_vcf"
colnames(beech_sample_info)[4] <- "Treeid"

beech_sample_info <- beech_sample_info[, c("SampleID_vcf", "RowNumber", "sample_nr_vcf", "Treeid")]

beech_sample_info_merged <- merge(beech_sample_info, beech_tree_info, by = "Treeid")
beech_sample_info_merged <- beech_sample_info_merged[order(beech_sample_info_merged$RowNumber), ]

beech_climate <- read.csv(file.path(climate_dir, "beech_clim_for_lfmm_091224.csv"))
beech_climate$prov[beech_climate$prov == "UK018"] <- "GB018"

beech_sample_ids <- beech_sample_info_merged[, c("SampleID_vcf", "Prov")]

beech_lfmm_climate <- merge(
  beech_sample_ids,
  beech_climate,
  by.x = "Prov",
  by.y = "prov"
)

write.csv(
  beech_lfmm_climate,
  file = file.path(climate_dir, "prepped_beech_evis_climate_data_correct_091224.csv")
)
