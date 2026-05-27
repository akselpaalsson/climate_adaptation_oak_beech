rm(list = ls())
gc()

library(dplyr)
library(lme4)

 
setwd("/path/to/working_directory")

# -----------------------------------------------------------------------------
# Helper functions
# -----------------------------------------------------------------------------

write_gemma <- function(x, file) {
  write.table(x, file = file, quote = FALSE, row.names = FALSE)
}

remove_colnames <- function(x) {
  colnames(x) <- NULL
  x
}

recode_provenance <- function(x) {
  x <- as.character(x)
  x <- gsub("FR", "FRA_", x)
  x <- gsub("DE", "GER_", x)
  x <- gsub("PL", "POL_", x)
  x <- gsub("UK", "GBR_", x)
  x <- gsub("DK", "DEN_", x)
  x
}

site_to_numeric <- function(x) {
  as.integer(dplyr::recode(as.character(x), "DK" = "1", "FR" = "2", "PL" = "3", "UK" = "4"))
}

block_to_numeric <- function(x) {
  as.integer(dplyr::recode(
    as.character(x),
    "DK_1" = "1", "DK_2" = "2", "DK_3" = "3", "DK_4" = "4", "DK_5" = "5", "DK_6" = "6",
    "FR_1" = "7", "FR_2" = "8", "FR_3" = "9", "FR_4" = "10", "FR_5" = "11",
    "PL_1" = "12", "PL_2" = "13", "PL_3" = "14", "PL_4" = "15", "PL_5" = "16", "PL_6" = "17", "PL_7" = "18",
    "UK_1" = "19", "UK_2" = "20"
  ))
}

read_vcf_sample_order <- function(file) {
  x <- read.table(file, sep = "\t", header = FALSE)
  colnames(x)[1] <- "SampleID"
  x
}

read_sample_info_reduced <- function() {
  sample_info <- read.table("common_ring_sample_info_no_reps.csv", sep = ",", header = TRUE)
  sample_vcf <- read.table("samples_populations_no_no277_reps_missing_removed_only_petrea_minDP3_meanDP10_maxDP80_NA0.90.txt", sep = "\t", header = FALSE)
  colnames(sample_vcf)[1] <- "Sample"
  sample_vcf$RowNumber <- as.numeric(row.names(sample_vcf))
  
  sample_info_reduced <- merge(sample_info, sample_vcf, by.x = "SampleID", by.y = "Sample", sort = FALSE)
  sample_info_reduced <- sample_info_reduced[order(sample_info_reduced$RowNumber), ]
  colnames(sample_info_reduced)[2] <- "Treeid"
  sample_info_reduced$merging_col <- paste(sample_info_reduced$Site, sample_info_reduced$Block, sample_info_reduced$Prov3, sample_info_reduced$Tree, sep = "_")
  sample_info_reduced
}

make_gemma_tables <- function(data, sample_order, trait_cols, cov_cols = c("CI_BAsum", "Site.x")) {
  data$Block_for_model <- paste(data$Site.x, data$Block.x, sep = "_")
  
  cov_df <- data[, c("SampleID", cov_cols)]
  rand_df <- data[, c("SampleID", "Block_for_model")]
  trait_df <- data[, c("SampleID", trait_cols)]
  
  idx <- match(sample_order$SampleID, cov_df$SampleID)
  
  cov_df_sorted <- cov_df[idx, ]
  rand_df_sorted <- rand_df[idx, ]
  trait_df_sorted <- trait_df[idx, ]
  
  cov_df_sorted$Site.x <- site_to_numeric(cov_df_sorted$Site.x)
  cov_df_sorted <- cov_df_sorted %>%
    mutate(First_Column = 1, .before = 1) %>%
    dplyr::select(-SampleID)
  
  rand_df_sorted$Block_for_model <- block_to_numeric(rand_df_sorted$Block_for_model)
  rand_df_sorted <- rand_df_sorted %>%
    mutate(First_Column = 1, .before = 1) %>%
    dplyr::select(-SampleID)
  
  trait_df_sorted <- trait_df_sorted %>% dplyr::select(-SampleID)
  
  list(
    order = idx,
    cov = remove_colnames(cov_df_sorted),
    cov_no_competition = remove_colnames(cov_df_sorted %>% dplyr::select(-CI_BAsum)),
    rand = remove_colnames(rand_df_sorted),
    trait = remove_colnames(trait_df_sorted)
  )
}

write_pc_covariates <- function(data, sample_order, output_prefix, pca_file = "pca_gemma_dbh_height_CI_BAsum_Site.x_covs_maf0.05.txt") {
  pca <- read.table(pca_file, sep = "\t", header = TRUE) %>% dplyr::select(-FID)
  colnames(pca)[1] <- "SampleID"
  
  data_pca <- merge(data, pca, by = "SampleID")
  
  for (n_pcs in 1:10) {
    pc_cols <- paste0("PC", seq_len(n_pcs))
    cov_df <- data_pca[, c("SampleID", pc_cols)]
    idx <- match(sample_order$SampleID, cov_df$SampleID)
    cov_df_sorted <- cov_df[idx, ] %>% dplyr::select(-SampleID)
    cov_df_sorted <- remove_colnames(cov_df_sorted)
    
    write_gemma(
      cov_df_sorted,
      sprintf("%s%d.txt", output_prefix, n_pcs)
    )
  }
}

model_residuals <- function(formula, data, mixed = FALSE) {
  if (mixed) {
    residuals(lmer(formula, data = data))
  } else {
    residuals(lm(formula, data = data))
  }
}

# -----------------------------------------------------------------------------
# Oak DBH and height residual traits
# -----------------------------------------------------------------------------

oak <- read.csv("oak_dbh.csv")
sample.info.reduced <- read_sample_info_reduced()

oak$Prov <- recode_provenance(oak$Prov)
oak$merging_col <- paste(oak$Site, oak$Block, oak$Prov, oak$Tree, sep = "_")

merged_test1 <- merge(sample.info.reduced, oak, by = "merging_col")
merged_test1$Block_for_model <- paste(merged_test1$Site.x, merged_test1$Block.x, sep = "_")

sample_gwas_vcf <- read_vcf_sample_order("samples_gemma_dbh_height_CI_BAsum_Site.X_covs_populations_no_no277_reps_missing_removed_only_petrea_minDP3_meanDP10_maxDP80_NA0.90.recode.vcf.txt")
gemma_dbh_height <- make_gemma_tables(merged_test1, sample_gwas_vcf, trait_cols = c("DBH_cm", "Height_m"))

residual_traits_dbh_height <- data.frame(
  SampleID = merged_test1$SampleID,
  residuals_dbh_vs_site = model_residuals(DBH_cm ~ Site.x, merged_test1),
  residuals_dbh_vs_site_and_ci_basum = model_residuals(DBH_cm ~ Site.x + CI_BAsum, merged_test1),
  residuals_dbh_vs_site_and_ci_basum_and_block = model_residuals(DBH_cm ~ Site.x + CI_BAsum + (1 | Block_for_model), merged_test1, mixed = TRUE),
  residuals_dbh_vs_site_and_block = model_residuals(DBH_cm ~ Site.x + (1 | Block_for_model), merged_test1, mixed = TRUE),
  residuals_dbh_vs_site_and_height = model_residuals(DBH_cm ~ Site.x + Height_m, merged_test1),
  residuals_dbh_vs_site_and_ci_basum_and_height = model_residuals(DBH_cm ~ Site.x + CI_BAsum + Height_m, merged_test1),
  residuals_dbh_vs_site_and_ci_basum_and_height_and_block = model_residuals(DBH_cm ~ Site.x + CI_BAsum + Height_m + (1 | Block_for_model), merged_test1, mixed = TRUE),
  residuals_dbh_vs_site_and_height_block = model_residuals(DBH_cm ~ Site.x + Height_m + (1 | Block_for_model), merged_test1, mixed = TRUE),
  residuals_height_vs_site = model_residuals(Height_m ~ Site.x, merged_test1),
  residuals_height_vs_site_and_ci_basum = model_residuals(Height_m ~ Site.x + CI_BAsum, merged_test1),
  residuals_height_vs_site_and_ci_basum_and_block = model_residuals(Height_m ~ Site.x + CI_BAsum + (1 | Block_for_model), merged_test1, mixed = TRUE),
  residuals_height_vs_site_and_block = model_residuals(Height_m ~ Site.x + (1 | Block_for_model), merged_test1, mixed = TRUE),
  residuals_dbh_vs_site_and_dbh = model_residuals(Height_m ~ Site.x + DBH_cm, merged_test1),
  residuals_height_vs_site_and_ci_basum_and_dbh = model_residuals(Height_m ~ Site.x + CI_BAsum + DBH_cm, merged_test1),
  residuals_height_vs_site_and_ci_basum_and_dbh_and_block = model_residuals(Height_m ~ Site.x + CI_BAsum + DBH_cm + (1 | Block_for_model), merged_test1, mixed = TRUE),
  residuals_height_vs_site_and_dbh_block = model_residuals(Height_m ~ Site.x + DBH_cm + (1 | Block_for_model), merged_test1, mixed = TRUE)
)

residual_traits_df_sorted <- residual_traits_dbh_height[gemma_dbh_height$order, ] %>% dplyr::select(-SampleID)
residual_traits_df_sorted <- remove_colnames(residual_traits_df_sorted)
write_gemma(residual_traits_df_sorted, "final_gemma_dbh_height_residual_trait_all_variations.txt")

write_pc_covariates(
  merged_test1,
  sample_gwas_vcf,
  output_prefix = "final_gemma_dbh_height_residual_cov_file_pc"
)

merged_test1 <- merged_test1 %>%
  group_by(Site.x) %>%
  mutate(DBH_cm_standardized = scale(DBH_cm)) %>%
  ungroup()

residual_traits_dbh_stand <- data.frame(
  SampleID = merged_test1$SampleID,
  residuals_dbh_stand_vs_competition_and_block = model_residuals(DBH_cm_standardized ~ CI_BAsum + (1 | Block_for_model), merged_test1, mixed = TRUE),
  residuals_dbh_stand_vs_competition = model_residuals(DBH_cm_standardized ~ CI_BAsum, merged_test1)
)

residual_traits_dbh_stand_sorted <- residual_traits_dbh_stand[gemma_dbh_height$order, ] %>% dplyr::select(-SampleID)
residual_traits_dbh_stand_sorted <- remove_colnames(as.data.frame(residual_traits_dbh_stand_sorted))
write_gemma(residual_traits_dbh_stand_sorted, "final_gemma_dbh_stand_by_site_residual_trait_all_variations.txt")

# -----------------------------------------------------------------------------
# Oak SLA residual traits
# -----------------------------------------------------------------------------

sla_df <- merged_test1 %>% filter(!is.na(SLA_cm2_mg))
sla_df$Block_for_model <- paste(sla_df$Site.x, sla_df$Block.x, sep = "_")

samples_for_gemma_SLA_CI_BAsum_Site.X_covs <- data.frame(sla_df = sla_df$SampleID)

sample_gwas_vcf <- read_vcf_sample_order("samples_gemma_sla_CI_BAsum_Site.X_covs_populations_no_no277_reps_missing_removed_only_petrea_minDP3_meanDP10_maxDP80_NA0.90.txt")

cov_df <- sla_df[, c("SampleID", "CI_BAsum", "Site.x", "DBH_cm")]
rand_df <- sla_df[, c("SampleID", "Block_for_model")]
trait_df <- sla_df[, c("SampleID", "SLA_cm2_mg")]
order <- match(sample_gwas_vcf$SampleID, cov_df$SampleID)

cov_df_sorted <- cov_df[order, ]
rand_df_sorted <- rand_df[order, ]
trait_df_sorted <- trait_df[order, ]

cov_df_sorted$Site.x <- site_to_numeric(cov_df_sorted$Site.x)
cov_df_sorted <- cov_df_sorted %>%
  mutate(First_Column = 1, .before = 1) %>%
  dplyr::select(-SampleID)

cov_df_sorted_no_competition <- cov_df_sorted %>% dplyr::select(-CI_BAsum)
cov_df_sorted_no_dbh <- cov_df_sorted %>% dplyr::select(-DBH_cm)
cov_df_sorted_no_dbh_no_competition <- cov_df_sorted_no_competition %>% dplyr::select(-DBH_cm)

trait_df_sorted <- trait_df_sorted %>% dplyr::select(-SampleID)

rand_df_sorted$Block_for_model <- block_to_numeric(rand_df_sorted$Block_for_model)
rand_df_sorted <- rand_df_sorted %>%
  mutate(First_Column = 1, .before = 1) %>%
  dplyr::select(-SampleID)

trait_df_sorted <- remove_colnames(trait_df_sorted)
cov_df_sorted <- remove_colnames(cov_df_sorted)
cov_df_sorted_no_competition <- remove_colnames(cov_df_sorted_no_competition)
cov_df_sorted_no_dbh <- remove_colnames(cov_df_sorted_no_dbh)
cov_df_sorted_no_dbh_no_competition <- remove_colnames(cov_df_sorted_no_dbh_no_competition)
rand_df_sorted <- remove_colnames(rand_df_sorted)

residual_traits_sla <- data.frame(
  SampleID = sla_df$SampleID,
  residuals_sla_vs_site = model_residuals(SLA_cm2_mg ~ Site.x, sla_df),
  residuals_sla_vs_site_and_ci_basum = model_residuals(SLA_cm2_mg ~ Site.x + CI_BAsum, sla_df),
  residuals_sla_vs_site_and_ci_basum_and_block = model_residuals(SLA_cm2_mg ~ Site.x + CI_BAsum + (1 | Block_for_model), sla_df, mixed = TRUE),
  residuals_sla_vs_site_and_block = model_residuals(SLA_cm2_mg ~ Site.x + (1 | Block_for_model), sla_df, mixed = TRUE),
  residuals_sla_vs_site_and_dbh = model_residuals(SLA_cm2_mg ~ Site.x + DBH_cm, sla_df),
  residuals_sla_vs_site_and_ci_basum_and_dbh = model_residuals(SLA_cm2_mg ~ Site.x + CI_BAsum + DBH_cm, sla_df),
  residuals_sla_vs_site_and_ci_basum_and_dbh_and_block = model_residuals(SLA_cm2_mg ~ Site.x + CI_BAsum + DBH_cm + (1 | Block_for_model), sla_df, mixed = TRUE),
  residuals_sla_vs_site_and_dbh_block = model_residuals(SLA_cm2_mg ~ Site.x + DBH_cm + (1 | Block_for_model), sla_df, mixed = TRUE)
)

residual_traits_sla_df_sorted <- residual_traits_sla[, c(
  "SampleID",
  "residuals_sla_vs_site_and_ci_basum",
  "residuals_sla_vs_site_and_ci_basum_and_block",
  "residuals_sla_vs_site_and_block",
  "residuals_sla_vs_site_and_dbh",
  "residuals_sla_vs_site_and_ci_basum_and_dbh",
  "residuals_sla_vs_site_and_ci_basum_and_dbh_and_block",
  "residuals_sla_vs_site_and_dbh_block"
)]
residual_traits_sla_df_sorted <- residual_traits_sla_df_sorted[order, ] %>% dplyr::select(-SampleID)
residual_traits_sla_df_sorted <- remove_colnames(residual_traits_sla_df_sorted)

write_pc_covariates(
  sla_df,
  sample_gwas_vcf,
  output_prefix = "final_gemma_sla_height_residual_cov_file_pc"
)

write_gemma(residual_traits_sla_df_sorted, "final_gemma_sla_height_residual_trait_all_variations.txt")
write_gemma(cov_df_sorted, "gemma_sla_cov_ci_sum_and_site.txt")
write_gemma(cov_df_sorted_no_dbh_no_competition, "gemma_sla_cov_site.txt")
write_gemma(cov_df_sorted_no_competition, "gemma_sla_cov_site_dbh.txt")
write_gemma(cov_df_sorted_no_dbh, "gemma_sla_cov_site_competition.txt")
write_gemma(trait_df_sorted, "gemma_sla_traits.txt")
write_gemma(rand_df_sorted, "gemma_sla_block_as_random.txt")
write_gemma(residual_traits_sla_df_sorted, "gemma_sla_residuals_as_traits_all_variations.txt")

# -----------------------------------------------------------------------------
# Oak basal area increment and cumulative growth traits
# -----------------------------------------------------------------------------

oak_bai <- read.csv("oak_bai.csv")
oak_bai$Treeid2 <- substr(oak_bai$Treeid, nchar(oak_bai$Treeid), nchar(oak_bai$Treeid))
oak_bai$Prov <- recode_provenance(oak_bai$Prov)
oak_bai$merging_col <- paste(oak_bai$Site, oak_bai$Block, oak_bai$Prov, oak_bai$Treeid2, sep = "_")

oak_bai_sample_info_reduced <- merge(sample.info.reduced, oak_bai, by = "merging_col")

fit_lm_and_extract_bai <- function(data) {
  lm_model <- lm(BAI ~ year, data = data)
  data.frame(Estimate = coef(lm_model)[2], SampleID = unique(data$SampleID))
}

bai_slopes_as_trait <- oak_bai_sample_info_reduced %>%
  group_by(SampleID) %>%
  do(fit_lm_and_extract_bai(.))

sorted_bai_by_sample <- split(oak_bai_sample_info_reduced, oak_bai_sample_info_reduced$SampleID)
sorted_bai_by_sample <- lapply(sorted_bai_by_sample, function(df) df[order(df$year), ])

cumulative_bai_by_sample <- lapply(sorted_bai_by_sample, function(df) {
  df$size_per_year <- cumsum(df$BAI)
  df
})

fit_lm_and_extract_growth <- function(df) {
  lm_model <- lm(size_per_year ~ year, data = df)
  data.frame(Estimate = coef(lm_model)[2], SampleID = unique(df$SampleID))
}

slopes_year_size_df <- bind_rows(lapply(cumulative_bai_by_sample, fit_lm_and_extract_growth))

samples_for_gemma_bai_traits <- data.frame(bai_slopes_as_trait = bai_slopes_as_trait$SampleID)
write_gemma(samples_for_gemma_bai_traits, "samples_for_gemma_bai_traits.txt")

bai_and_growth_slopes_as_trait <- merge(bai_slopes_as_trait, slopes_year_size_df, by = "SampleID")
colnames(bai_and_growth_slopes_as_trait)[2] <- "bai"
colnames(bai_and_growth_slopes_as_trait)[3] <- "growth_year"

bai_and_growth_slopes_as_trait_df <- merge(bai_and_growth_slopes_as_trait, merged_test1, by = "SampleID")
bai_and_growth_slopes_as_trait_df$Block_for_model <- paste(bai_and_growth_slopes_as_trait_df$Site.x, bai_and_growth_slopes_as_trait_df$Block.x, sep = "_")

sample_gwas_vcf <- read_vcf_sample_order("samples_gemma_bai_traits_populations_no_no277_reps_missing_removed_only_petrea_minDP3_meanDP10_maxDP80_NA0.90.txt")
gemma_bai <- make_gemma_tables(bai_and_growth_slopes_as_trait_df, sample_gwas_vcf, trait_cols = c("bai", "growth_year"))

residual_traits_bai <- data.frame(
  SampleID = bai_and_growth_slopes_as_trait_df$SampleID,
  residuals_bai_vs_site = model_residuals(bai ~ Site.x, bai_and_growth_slopes_as_trait_df),
  residuals_bai_vs_site_and_ci_basum = model_residuals(bai ~ Site.x + CI_BAsum, bai_and_growth_slopes_as_trait_df),
  residuals_bai_vs_site_and_ci_basum_and_block = model_residuals(bai ~ Site.x + CI_BAsum + (1 | Block_for_model), bai_and_growth_slopes_as_trait_df, mixed = TRUE),
  residuals_bai_vs_site_and_block = model_residuals(bai ~ Site.x + (1 | Block_for_model), bai_and_growth_slopes_as_trait_df, mixed = TRUE),
  residuals_growth_year_vs_site = model_residuals(growth_year ~ Site.x, bai_and_growth_slopes_as_trait_df),
  residuals_growth_year_vs_site_and_ci_basum = model_residuals(growth_year ~ Site.x + CI_BAsum, bai_and_growth_slopes_as_trait_df),
  residuals_growth_year_vs_site_and_ci_basum_and_block = model_residuals(growth_year ~ Site.x + CI_BAsum + (1 | Block_for_model), bai_and_growth_slopes_as_trait_df, mixed = TRUE),
  residuals_growth_year_vs_site_and_block = model_residuals(growth_year ~ Site.x + (1 | Block_for_model), bai_and_growth_slopes_as_trait_df, mixed = TRUE)
)

residual_traits_bai_sorted <- residual_traits_bai[gemma_bai$order, ] %>% dplyr::select(-SampleID)
residual_traits_bai_sorted <- remove_colnames(residual_traits_bai_sorted)

write_gemma(gemma_bai$cov, "gemma_bai_traits_cov_ci_sum_and_site.txt")
write_gemma(gemma_bai$cov_no_competition, "gemma_bai_traits_cov_site.txt")
write_gemma(gemma_bai$trait, "gemma_bai_traits.txt")
write_gemma(gemma_bai$rand, "gemma_bai_block_as_random.txt")
write_gemma(residual_traits_bai_sorted, "gemma_bai_residual_traits.txt")

# -----------------------------------------------------------------------------
# Oak drought-response traits: resilience, resistance, recovery
# -----------------------------------------------------------------------------

oak_rs <- read.csv("oak_rs-2.csv")
sample.info.reduced <- read_sample_info_reduced()

oak_rs$Prov <- recode_provenance(oak_rs$Prov)
oak_rs$merging_col <- paste(oak_rs$Site, oak_rs$Block, oak_rs$Prov, oak_rs$tree, sep = "_")

oak_rs_rs <- subset(oak_rs, var == "rs")
oak_rs_rc <- subset(oak_rs, var == "rc")
oak_rs_rt <- subset(oak_rs, var == "rt")

merged_oak_rs <- merge(sample.info.reduced, oak_rs_rs, by = "merging_col")
merged_oak_rc <- merge(sample.info.reduced, oak_rs_rc, by = "merging_col")
merged_oak_rt <- merge(sample.info.reduced, oak_rs_rt, by = "merging_col")

samples_for_gemma_rs <- data.frame(merged_oak_rs = merged_oak_rs$SampleID)
samples_for_gemma_rc <- data.frame(merged_oak_rc = merged_oak_rc$SampleID)
samples_for_gemma_rt <- data.frame(merged_oak_rt = merged_oak_rt$SampleID)

write_gemma(samples_for_gemma_rs, "samples_for_merged_oak_rs.txt")
write_gemma(samples_for_gemma_rc, "samples_for_merged_oak_rc.txt")
write_gemma(samples_for_gemma_rt, "samples_for_merged_oak_rt.txt")

prepare_response_trait <- function(data, sample_file) {
  sample_order <- read_vcf_sample_order(sample_file)
  gemma <- make_gemma_tables(data, sample_order, trait_cols = "value")
  
  data$Block_for_model <- paste(data$Site.x, data$Block.x, sep = "_")
  
  list(
    sample_order = sample_order,
    gemma = gemma,
    residuals = data.frame(
      SampleID = data$SampleID,
      residuals_vs_site = model_residuals(value ~ Site.x, data),
      residuals_site_and_ci_basum = model_residuals(value ~ Site.x + CI_BAsum, data),
      residuals_vs_site_and_ci_basum_and_block = model_residuals(value ~ Site.x + CI_BAsum + (1 | Block_for_model), data, mixed = TRUE),
      residuals_vs_site_and_block = model_residuals(value ~ Site.x + (1 | Block_for_model), data, mixed = TRUE)
    )
  )
}

rs_prepared <- prepare_response_trait(
  merged_oak_rs,
  "samples_gemma_rs_populations_no_no277_reps_missing_removed_only_petrea_minDP3_meanDP10_maxDP80_NA0.90.MAF0.01.recode.vcf.txt"
)

rt_prepared <- prepare_response_trait(
  merged_oak_rt,
  "samples_gemma_rt_populations_no_no277_reps_missing_removed_only_petrea_minDP3_meanDP10_maxDP80_NA0.90.MAF0.01.recode.vcf.txt"
)

rc_prepared <- prepare_response_trait(
  merged_oak_rc,
  "samples_gemma_rc_populations_no_no277_reps_missing_removed_only_petrea_minDP3_meanDP10_maxDP80_NA0.90.MAF0.01.recode.vcf.txt"
)

residual_traits_rc <- rc_prepared$residuals
colnames(residual_traits_rc) <- c(
  "SampleID",
  "residuals_rc_vs_site",
  "residuals_rc_site_and_ci_basum",
  "residuals_rc_vs_site_and_ci_basum_and_block",
  "residuals_rc_vs_site_and_block"
)

residual_traits_rc_sorted <- residual_traits_rc[rc_prepared$gemma$order, ] %>% dplyr::select(-SampleID)
residual_traits_rc_sorted <- remove_colnames(residual_traits_rc_sorted)

write_gemma(rc_prepared$gemma$trait, "gemma_rc_trait.txt")
write_gemma(rc_prepared$gemma$cov, "gemma_rc_cov_competition_and_site.txt")
write_gemma(rc_prepared$gemma$cov_no_competition, "gemma_rc_cov_site.txt")
write_gemma(residual_traits_rc_sorted, "gemma_rc_residual_traits_all_combinations.txt")









# -----------------------------------------------------------------------------
# Oak climate-correlation response traits: precipitation and temperature
# -----------------------------------------------------------------------------

make_corr_trait_data <- function(sample_info, oak_corr, climate_variable) {
  corr_df <- subset(oak_corr, var0me == climate_variable)
  
  corr_df <- corr_df %>%
    mutate(tree_nr = as.numeric(sub(".*([0-9]+)$", "\\1", tree.id)))
  
  corr_df$Prov <- recode_provenance(corr_df$Prov)
  corr_df$merging_col <- paste(corr_df$Site, corr_df$Block, corr_df$Prov, corr_df$tree_nr, sep = "_")
  
  merged_df <- merge(sample_info, corr_df, by = "merging_col")
  merged_df <- merged_df[!duplicated(merged_df), ]
  merged_df$Block_for_model <- paste(merged_df$Site.x, merged_df$Block.x, sep = "_")
  
  merged_df
}

make_corr_gemma_tables <- function(data, sample_order) {
  cov_df <- data[, c("SampleID", "CI_Basum", "Site.x")]
  trait_df <- data[, c("SampleID", "SPR", "SUM")]
  
  idx <- match(sample_order$SampleID, cov_df$SampleID)
  
  cov_df_sorted <- cov_df[idx, ]
  trait_df_sorted <- trait_df[idx, ]
  
  cov_df_sorted$Site.x <- site_to_numeric(cov_df_sorted$Site.x)
  cov_df_sorted <- cov_df_sorted %>%
    mutate(First_Column = 1, .before = 1) %>%
    dplyr::select(-SampleID)
  
  cov_df_sorted_no_competition <- cov_df_sorted %>%
    dplyr::select(-CI_Basum)
  
  trait_df_sorted <- trait_df_sorted %>%
    dplyr::select(-SampleID)
  
  list(
    order = idx,
    trait = remove_colnames(trait_df_sorted),
    cov = remove_colnames(cov_df_sorted),
    cov_no_competition = remove_colnames(cov_df_sorted_no_competition)
  )
}

make_corr_residual_traits <- function(data, sample_order) {
  residual_traits <- data.frame(
    SampleID = data$SampleID,
    residuals_spr_vs_site = model_residuals(SPR ~ Site.x, data),
    residuals_spr_vs_site_and_ci_basum = model_residuals(SPR ~ Site.x + CI_Basum, data),
    residuals_sum_vs_site = model_residuals(SUM ~ Site.x, data),
    residuals_sum_vs_site_and_ci_basum = model_residuals(SUM ~ Site.x + CI_Basum, data),
    residuals_spr_vs_site_and_ci_basum_block = model_residuals(SPR ~ Site.x + CI_Basum + (1 | Block_for_model), data, mixed = TRUE),
    residuals_sum_vs_site_and_ci_basum_block = model_residuals(SUM ~ Site.x + CI_Basum + (1 | Block_for_model), data, mixed = TRUE)
  )
  
  residual_traits_sorted <- residual_traits[sample_order, ] %>%
    dplyr::select(-SampleID)
  
  remove_colnames(residual_traits_sorted)
}

make_pc_suffix <- function(n_pcs) {
  if (n_pcs == 1) {
    "pc1"
  } else {
    paste0("pc1_pc", n_pcs)
  }
}

write_corr_pc_covariates <- function(data, sample_order, pca_file, output_prefix) {
  pca <- read.table(pca_file, sep = "\t", header = TRUE) %>%
    dplyr::select(-FID)
  
  colnames(pca)[1] <- "SampleID"
  
  data_pca <- merge(data, pca, by = "SampleID")
  
  for (n_pcs in 1:10) {
    pc_cols <- paste0("PC", seq_len(n_pcs))
    
    cov_df <- data_pca[, c("SampleID", pc_cols)]
    idx <- match(sample_order$SampleID, cov_df$SampleID)
    
    cov_df_sorted <- cov_df[idx, ] %>%
      dplyr::select(-SampleID)
    
    cov_df_sorted <- remove_colnames(cov_df_sorted)
    
    write_gemma(
      cov_df_sorted,
      paste0(output_prefix, make_pc_suffix(n_pcs), ".txt")
    )
  }
}

write_corr_site_pc_covariates <- function(data, sample_order, pca_file, output_prefix) {
  pca <- read.table(pca_file, sep = "\t", header = TRUE) %>%
    dplyr::select(-FID)
  
  colnames(pca)[1] <- "SampleID"
  
  data_pca <- merge(data, pca, by = "SampleID")
  
  for (n_pcs in 1:10) {
    pc_cols <- paste0("PC", seq_len(n_pcs))
    
    cov_df <- data_pca[, c("SampleID", "Site.x", pc_cols)]
    idx <- match(sample_order$SampleID, cov_df$SampleID)
    
    cov_df_sorted <- cov_df[idx, ]
    cov_df_sorted$Site.x <- site_to_numeric(cov_df_sorted$Site.x)
    
    cov_df_sorted <- cov_df_sorted %>%
      mutate(First_Column = 1, .before = 1) %>%
      dplyr::select(-SampleID)
    
    cov_df_sorted <- remove_colnames(cov_df_sorted)
    
    write_gemma(
      cov_df_sorted,
      paste0(output_prefix, make_pc_suffix(n_pcs), ".txt")
    )
  }
}

export_plink_pca <- function(eigenvec_file, output_file) {
  pca <- read.table(eigenvec_file, header = FALSE)
  
  colnames(pca) <- c(
    "FID",
    "IID",
    paste0("PC", seq_len(ncol(pca) - 2))
  )
  
  write.table(
    pca,
    output_file,
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )
}

sample_info <- read_sample_info_reduced()
oak_corr <- read.csv("oak_corr.csv")


# -----------------------------------------------------------------------------
# Precipitation mean: SPR and SUM correlation-response traits
# -----------------------------------------------------------------------------

prec_mean <- make_corr_trait_data(
  sample_info = sample_info,
  oak_corr = oak_corr,
  climate_variable = "prec.mean"
)

samples_for_gemma_prec_mean <- data.frame(merged_test1 = prec_mean$SampleID)
write_gemma(samples_for_gemma_prec_mean, "samples_for_gemma_prec_mean.txt")

sample_gwas_vcf <- read.table(
  "gemma_prec_mean_populations_no_no277_reps_missing_removed_only_petrea_minDP3_meanDP10_maxDP80_NA0.90.MAF0.05.recode.vcf.txt",
  sep = "\t",
  header = FALSE
)
colnames(sample_gwas_vcf)[1] <- "SampleID"

prec_mean_gemma <- make_corr_gemma_tables(
  data = prec_mean,
  sample_order = sample_gwas_vcf
)

prec_mean_residual_traits <- make_corr_residual_traits(
  data = prec_mean,
  sample_order = prec_mean_gemma$order
)

write_gemma(
  prec_mean_gemma$trait,
  "gemma_prec_mean_spr_and_sum_trait.txt"
)

write_gemma(
  prec_mean_gemma$cov,
  "gemma_prec_mean_spr_and_sum_cov_competition_and_site.txt"
)

write_gemma(
  prec_mean_gemma$cov_no_competition,
  "gemma_prec_mean_spr_and_sum_cov_site.txt"
)

write_gemma(
  prec_mean_residual_traits,
  "gemma_prec_mean_spr_and_sum_residual_traits_all_combinations.txt"
)

write_gemma(
  prec_mean_gemma$cov,
  "gemma_prec_mean_spr_and_sum_cov_ci_sum_and_site.txt"
)

write_gemma(
  prec_mean_gemma$cov_no_competition,
  "gemma_prec_mean_spr_and_sum_cov_site.txt"
)

write_corr_pc_covariates(
  data = prec_mean,
  sample_order = sample_gwas_vcf,
  pca_file = "pca_gemma_prec_mean_maf0.05.txt",
  output_prefix = "gemma_prec_mean_spr_and_sum_"
)


# -----------------------------------------------------------------------------
# Temperature mean: SPR and SUM correlation-response traits
# -----------------------------------------------------------------------------

tmean_mean <- make_corr_trait_data(
  sample_info = sample_info,
  oak_corr = oak_corr,
  climate_variable = "tmea.mean"
)

samples_for_gemma_tmean_mean <- data.frame(merged_test1 = tmean_mean$SampleID)
write_gemma(samples_for_gemma_tmean_mean, "samples_for_gemma_tmean_mean.txt")

export_plink_pca(
  eigenvec_file = "PCA_gemma_tmean_mean_populations_no_no277_reps_missing_removed_only_petrea_minDP3_meanDP10_maxDP80_NA0.90.MAF0.05.eigenvec",
  output_file = "pca_gemma_tmean_mean_maf0.05.txt"
)

sample_gwas_vcf <- read.table(
  "gemma_tmean_mean_populations_no_no277_reps_missing_removed_only_petrea_minDP3_meanDP10_maxDP80_NA0.90.MAF0.05.recode.vcf.txt",
  sep = "\t",
  header = FALSE
)
colnames(sample_gwas_vcf)[1] <- "SampleID"

tmean_mean_gemma <- make_corr_gemma_tables(
  data = tmean_mean,
  sample_order = sample_gwas_vcf
)

tmean_mean_residual_traits <- make_corr_residual_traits(
  data = tmean_mean,
  sample_order = tmean_mean_gemma$order
)

write_gemma(
  tmean_mean_gemma$trait,
  "gemma_tmean_mean_spr_and_sum_trait.txt"
)

write_gemma(
  tmean_mean_gemma$cov,
  "gemma_tmean_mean_spr_and_sum_cov_competition_and_site.txt"
)

write_gemma(
  tmean_mean_gemma$cov_no_competition,
  "gemma_tmean_mean_spr_and_sum_cov_site.txt"
)

write_gemma(
  tmean_mean_residual_traits,
  "gemma_tmean_mean_spr_and_sum_residual_traits_all_combinations.txt"
)

write_gemma(
  tmean_mean_gemma$cov,
  "gemma_tmean_mean_spr_and_sum_cov_ci_sum_and_site.txt"
)

write_gemma(
  tmean_mean_gemma$cov_no_competition,
  "gemma_tmean_mean_spr_and_sum_cov_site.txt"
)

write_corr_pc_covariates(
  data = tmean_mean,
  sample_order = sample_gwas_vcf,
  pca_file = "pca_gemma_tmean_mean_maf0.05.txt",
  output_prefix = "gemma_tmean_mean_spr_and_sum_"
)

write_corr_site_pc_covariates(
  data = tmean_mean,
  sample_order = sample_gwas_vcf,
  pca_file = "pca_gemma_tmean_mean_maf0.05.txt",
  output_prefix = "gemma_tmean_mean_spr_and_sum_cov_site_"
)
