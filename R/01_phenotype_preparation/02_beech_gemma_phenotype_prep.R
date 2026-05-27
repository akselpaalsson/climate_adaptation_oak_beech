
rm(list = ls())
gc()

library(dplyr)
library(lme4)

setwd("/path/to/working_directory")

# -----------------------------------------------------------------------------
# Helper functions
# -----------------------------------------------------------------------------

write_gemma <- function(x, file, sep = " ") {
  write.table(x, file = file, quote = FALSE, row.names = FALSE, sep = sep)
}

remove_colnames <- function(x) {
  colnames(x) <- NULL
  x
}

read_beech_vcf_samples <- function(file = "samples_populations_fagus_no_high_miss_ind_and_incorr_removed_minDP3_meanDP10_maxDP80_NA0.90_reduced.MAF0.05.txt") {
  samples <- read.table(file, sep = "\t", header = FALSE)
  
  split_ids <- strsplit(as.character(samples$V1), "-")
  samples <- cbind(samples, do.call(rbind, split_ids))
  
  colnames(samples)[1:5] <- c(
    "SampleID_vcf",
    "RowNumber",
    "Treeid",
    "seq_id",
    "vcf_extra"
  )
  
  samples
}

read_sample_order <- function(file) {
  sample_order <- read.table(file, sep = "\t", header = FALSE)
  colnames(sample_order)[1] <- "SampleID_vcf"
  sample_order
}

model_residuals <- function(formula, data, mixed = FALSE) {
  if (mixed) {
    residuals(lmer(formula, data = data))
  } else {
    residuals(lm(formula, data = data))
  }
}

make_residual_traits <- function(data, response, prefix) {
  data$Block_for_model <- data$Blocknew
  
  residual_traits <- data.frame(
    SampleID_vcf = data$SampleID_vcf
  )
  
  residual_traits[[paste0("residuals_", prefix, "_vs_site_and_ci_basum_and_block")]] <-
    model_residuals(
      as.formula(paste0(response, " ~ Site + CI_BAsum + (1 | Block_for_model)")),
      data,
      mixed = TRUE
    )
  
  residual_traits[[paste0("residuals_", prefix, "_vs_site_and_block")]] <-
    model_residuals(
      as.formula(paste0(response, " ~ Site + (1 | Block_for_model)")),
      data,
      mixed = TRUE
    )
  
  residual_traits[[paste0("residuals_", prefix, "_vs_site")]] <-
    model_residuals(
      as.formula(paste0(response, " ~ Site")),
      data
    )
  
  residual_traits[[paste0("residuals_", prefix, "_vs_site_and_ci_basum")]] <-
    model_residuals(
      as.formula(paste0(response, " ~ Site + CI_BAsum")),
      data
    )
  
  residual_traits
}

sort_residual_traits <- function(residual_traits, sample_order_file) {
  sample_order <- read_sample_order(sample_order_file)
  idx <- match(sample_order$SampleID_vcf, residual_traits$SampleID_vcf)
  
  residual_traits_sorted <- residual_traits[idx, ] %>%
    select(-SampleID_vcf)
  
  remove_colnames(residual_traits_sorted)
}

export_pca <- function(eigenvec_file, output_file) {
  pca <- read.table(eigenvec_file, header = FALSE)
  
  colnames(pca)[1:ncol(pca)] <- c(
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
  
  pca
}

prepare_pca <- function(pca) {
  pca <- pca %>% select(-FID)
  colnames(pca)[1] <- "SampleID_vcf"
  pca
}

pc_suffix <- function(n_pcs) {
  paste0("pc", n_pcs)
}

write_pc_covariates <- function(data, pca, sample_order_file, output_prefix, max_pcs = 10) {
  sample_order <- read_sample_order(sample_order_file)
  pca <- prepare_pca(pca)
  data_pca <- merge(data, pca, by = "SampleID_vcf")
  
  for (n_pcs in seq_len(max_pcs)) {
    pc_cols <- paste0("PC", seq_len(n_pcs))
    cov_df <- data_pca[, c("SampleID_vcf", pc_cols)]
    
    idx <- match(sample_order$SampleID_vcf, cov_df$SampleID_vcf)
    
    cov_df_sorted <- cov_df[idx, ] %>%
      select(-SampleID_vcf)
    
    cov_df_sorted <- remove_colnames(cov_df_sorted)
    
    write_gemma(
      cov_df_sorted,
      paste0(output_prefix, pc_suffix(n_pcs), ".txt")
    )
  }
}

fit_bai_slope <- function(data) {
  lm_model <- lm(BAI_mm2 ~ year, data = data)
  
  data.frame(
    Estimate = coef(lm_model)[2],
    SampleID_vcf = unique(data$SampleID_vcf)
  )
}


# -----------------------------------------------------------------------------
# DBH, height, and SLA
# -----------------------------------------------------------------------------

samples_in_vcf <- read_beech_vcf_samples()
beech <- read.csv("beech_dbh.csv", sep = ",")

merged_traits <- merge(samples_in_vcf, beech, by = "Treeid")
merged_traits <- merged_traits[!is.na(merged_traits$CI_BAsum), ]

dbh_df <- merged_traits[!is.na(merged_traits$DBH_cm), ]
height_df <- merged_traits[!is.na(merged_traits$Height_m), ]
sla_df <- merged_traits[!is.na(merged_traits$SLA_cm2_mg), ]

beech_samples_dbh_gwas <- as.data.frame(dbh_df$SampleID_vcf)
write_gemma(beech_samples_dbh_gwas, "beech_samples_dbh_gwas.txt")

beech_samples_height_gwas <- as.data.frame(height_df$SampleID_vcf)
write_gemma(beech_samples_height_gwas, "beech_samples_height_gwas.txt")

dbh_residual_traits <- make_residual_traits(
  data = dbh_df,
  response = "DBH_cm",
  prefix = "dbh"
)

height_residual_traits <- make_residual_traits(
  data = height_df,
  response = "Height_m",
  prefix = "height"
)

sla_residual_traits <- make_residual_traits(
  data = sla_df,
  response = "SLA_cm2_mg",
  prefix = "sla"
)

write_gemma(
  sort_residual_traits(dbh_residual_traits, "samples_beech_dbh_gwas.txt"),
  "beech_dbh_residual_traits1.txt"
)

write_gemma(
  sort_residual_traits(height_residual_traits, "samples_beech_height_gwas.txt"),
  "beech_height_residual_traits.txt"
)

write_gemma(
  sort_residual_traits(sla_residual_traits, "samples_beech_sla_gwas.txt"),
  "beech_sla_residual_traits.txt"
)

pca_dbh <- export_pca(
  eigenvec_file = "PCA_beech_dbh_gwas.eigenvec",
  output_file = "pca_gemma_beech_dbh.txt"
)

pca_height <- export_pca(
  eigenvec_file = "PCA_beech_height_gwas.eigenvec",
  output_file = "pca_gemma_beech_height.txt"
)

pca_sla <- export_pca(
  eigenvec_file = "PCA_beech_sla_gwas.eigenvec",
  output_file = "pca_gemma_beech_sla.txt"
)

write_pc_covariates(
  data = dbh_df,
  pca = pca_dbh,
  sample_order_file = "samples_beech_dbh_gwas.txt",
  output_prefix = "beech_dbh_residual_cov_file_",
  max_pcs = 10
)

write_pc_covariates(
  data = height_df,
  pca = pca_height,
  sample_order_file = "samples_beech_height_gwas.txt",
  output_prefix = "beech_height_residual_cov_file_",
  max_pcs = 10
)

write_pc_covariates(
  data = sla_df,
  pca = pca_sla,
  sample_order_file = "samples_beech_sla_gwas.txt",
  output_prefix = "beech_sla_residual_cov_file_",
  max_pcs = 10
)


# -----------------------------------------------------------------------------
# BAI slope
# -----------------------------------------------------------------------------

beech_bai <- read.csv("beech_bai.csv")
samples_in_vcf <- read_beech_vcf_samples()

merged_bai <- merge(samples_in_vcf, beech_bai, by = "Treeid")

bai_slopes_as_trait <- merged_bai %>%
  group_by(SampleID_vcf) %>%
  do(fit_bai_slope(.))

merged_and_bai_slopes <- merge(bai_slopes_as_trait, merged_bai, by = "SampleID_vcf") %>%
  distinct(SampleID_vcf, .keep_all = TRUE)

beech_samples_gwas_bai <- as.data.frame(merged_and_bai_slopes$SampleID_vcf)
write_gemma(beech_samples_gwas_bai, "beech_samples_gwas_bai.txt")

bai_df <- merged_and_bai_slopes

bai_residual_traits <- make_residual_traits(
  data = bai_df,
  response = "Estimate",
  prefix = "bai"
)

write_gemma(
  sort_residual_traits(bai_residual_traits, "samples_beech_bai_gwas.txt"),
  "beech_bai_residual_traits.txt"
)

pca_bai <- export_pca(
  eigenvec_file = "PCA_beech_bai_gwas.eigenvec",
  output_file = "pca_gemma_beech_bai.txt"
)

write_pc_covariates(
  data = bai_df,
  pca = pca_bai,
  sample_order_file = "samples_beech_bai_gwas.txt",
  output_prefix = "beech_bai_residual_cov_file_",
  max_pcs = 20
)


# -----------------------------------------------------------------------------
# Resilience, resistance, and recovery traits
# -----------------------------------------------------------------------------

beech_rs <- read.csv("beech_rs.csv")
samples_in_vcf <- read_beech_vcf_samples()

merged_rs <- merge(samples_in_vcf, beech_rs, by = "Treeid")

beech_rs_rs <- subset(merged_rs, var == "rs")
beech_rs_rt <- subset(merged_rs, var == "rt")
beech_rs_rc <- subset(merged_rs, var == "rc")

beech_samples_gwas_rs <- as.data.frame(beech_rs_rs$SampleID_vcf)
write_gemma(beech_samples_gwas_rs, "beech_samples_gwas_rs.txt")

beech_samples_gwas_rt <- as.data.frame(beech_rs_rt$SampleID_vcf)
write_gemma(beech_samples_gwas_rt, "beech_samples_gwas_rt.txt")

beech_samples_gwas_rc <- as.data.frame(beech_rs_rc$SampleID_vcf)
write_gemma(beech_samples_gwas_rc, "beech_samples_gwas_rc.txt")

rs_trait_map <- list(
  list(
    data = beech_rs_rs,
    response = "value",
    prefix = "rs",
    sample_order = "samples_beech_rs_gwas.txt",
    residual_output = "beech_rs_residual_traits.txt",
    eigenvec = "PCA_beech_rs_gwas.eigenvec",
    pca_output = "pca_gemma_beech_rs.txt",
    pc_output_prefix = "beech_rs_residual_cov_file_"
  ),
  list(
    data = beech_rs_rt,
    response = "value",
    prefix = "rt",
    sample_order = "samples_beech_rt_gwas.txt",
    residual_output = "beech_rt_residual_traits.txt",
    eigenvec = "PCA_beech_rt_gwas.eigenvec",
    pca_output = "pca_gemma_beech_rt.txt",
    pc_output_prefix = "beech_rt_residual_cov_file_"
  ),
  list(
    data = beech_rs_rc,
    response = "value",
    prefix = "rc",
    sample_order = "samples_beech_rc_gwas.txt",
    residual_output = "beech_rc_residual_traits.txt",
    eigenvec = "PCA_beech_rc_gwas.eigenvec",
    pca_output = "pca_gemma_beech_rc.txt",
    pc_output_prefix = "beech_rc_residual_cov_file_"
  )
)

for (trait in rs_trait_map) {
  trait_data <- trait$data
  
  residual_traits <- make_residual_traits(
    data = trait_data,
    response = trait$response,
    prefix = trait$prefix
  )
  
  write_gemma(
    sort_residual_traits(residual_traits, trait$sample_order),
    trait$residual_output
  )
  
  pca <- export_pca(
    eigenvec_file = trait$eigenvec,
    output_file = trait$pca_output
  )
  
  write_pc_covariates(
    data = trait_data,
    pca = pca,
    sample_order_file = trait$sample_order,
    output_prefix = trait$pc_output_prefix,
    max_pcs = 10
  )
}


# -----------------------------------------------------------------------------
# Climate-correlation traits
# -----------------------------------------------------------------------------

beech_corr <- read.csv("beech_corr-2.csv")
samples_in_vcf <- read_beech_vcf_samples()

merged_corr <- merge(samples_in_vcf, beech_corr, by = "Treeid")

prec_mean <- subset(merged_corr, varname == "prec.mean")
tmean_mean <- subset(merged_corr, varname == "tmea.mean")

beech_samples_gwas_prec_mean_corr <- as.data.frame(prec_mean$SampleID_vcf)
write_gemma(beech_samples_gwas_prec_mean_corr, "beech_samples_gwas_prec_mean_corr.txt")

beech_samples_gwas_tmean_mean_corr <- as.data.frame(tmean_mean$SampleID_vcf)
write_gemma(beech_samples_gwas_tmean_mean_corr, "beech_samples_gwas_tmean_mean_corr.txt")

corr_trait_map <- list(
  list(
    data = prec_mean,
    sample_order = "samples_beech_prec_mean_corr_gwas.txt",
    eigenvec = "PCA_beech_corr_prec_mean_gwas.eigenvec",
    pca_output = "pca_gemma_beech_corr_prec_mean.txt",
    outputs = list(
      list(
        response = "SPR",
        residual_prefix = "corr",
        residual_output = "beech_corr_prec_mean_spr_residual_traits.txt",
        pc_output_prefix = "beech_corr_prec_mean_spr_residual_cov_file_"
      ),
      list(
        response = "SUM",
        residual_prefix = "corr",
        residual_output = "beech_corr_prec_mean_sum_residual_traits.txt",
        pc_output_prefix = "beech_corr_prec_mean_sum_residual_cov_file_"
      )
    )
  ),
  list(
    data = tmean_mean,
    sample_order = "samples_beech_tmean_mean_corr_gwas.txt",
    eigenvec = "PCA_beech_corr_tmean_mean_gwas.eigenvec",
    pca_output = "pca_gemma_beech_corr_tmean_mean.txt",
    outputs = list(
      list(
        response = "SPR",
        residual_prefix = "corr",
        residual_output = "beech_corr_tmean_mean_spr_residual_traits.txt",
        pc_output_prefix = "beech_corr_tmean_mean_spr_residual_cov_file_"
      ),
      list(
        response = "SUM",
        residual_prefix = "corr",
        residual_output = "beech_corr_tmean_mean_sum_residual_traits.txt",
        pc_output_prefix = "beech_corr_tmean_mean_sum_residual_cov_file_"
      )
    )
  )
)

for (corr_trait in corr_trait_map) {
  corr_data <- corr_trait$data
  corr_data$Block_for_model <- corr_data$Blocknew
  
  pca <- export_pca(
    eigenvec_file = corr_trait$eigenvec,
    output_file = corr_trait$pca_output
  )
  
  for (output in corr_trait$outputs) {
    residual_traits <- make_residual_traits(
      data = corr_data,
      response = output$response,
      prefix = output$residual_prefix
    )
    
    write_gemma(
      sort_residual_traits(residual_traits, corr_trait$sample_order),
      output$residual_output
    )
    
    write_pc_covariates(
      data = corr_data,
      pca = pca,
      sample_order_file = corr_trait$sample_order,
      output_prefix = output$pc_output_prefix,
      max_pcs = 10
    )
  }
}
