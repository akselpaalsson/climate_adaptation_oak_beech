
rm(list = ls())
gc()

library(dplyr)
library(qqman)

# Change these paths if the repository is moved.
gwas_traits_dir <- "path/to/wd1"
corr_recheck_dir <- "path/to/wd2"

bh_fdr_threshold <- function(pvals, fdr = 0.05) {
  pvals <- sort(pvals, decreasing = FALSE)
  m <- length(pvals)
  cutoffs <- (seq_len(m) / m) * fdr
  passed <- which(pvals <= cutoffs)
  k <- max(c(passed, 0))

  ((0:m) / m * fdr)[k + 1]
}

beech_chr_map <- c(
  "Bhaga_1" = 1,
  "Bhaga_2" = 2,
  "Bhaga_3" = 3,
  "Bhaga_4" = 4,
  "Bhaga_5" = 5,
  "Bhaga_6" = 6,
  "Bhaga_7" = 7,
  "Bhaga_8" = 8,
  "Bhaga_9" = 9,
  "Bhaga_10" = 10,
  "Bhaga_11" = 11,
  "Bhaga_12" = 12,
  "Bhaga_Mitochondria_Circular" = 14,
  "Bhaga_Unplaced_1180" = 14,
  "Bhaga_Unplaced_1192" = 14,
  "Bhaga_Unplaced_1597" = 14,
  "Bhaga_Unplaced_1952" = 14,
  "Bhaga_Unplaced_2098" = 14,
  "Bhaga_Unplaced_2332" = 14,
  "Bhaga_Unplaced_255" = 14,
  "Bhaga_Unplaced_2901" = 14,
  "Bhaga_Unplaced_317" = 14,
  "Bhaga_Unplaced_3606" = 14,
  "Bhaga_Unplaced_4025" = 14,
  "Bhaga_Unplaced_565" = 14,
  "Bhaga_Unplaced_621" = 14
)

independent_snp_count <- function(ld_file) {
  ld_data <- read.table(ld_file, sep = "", header = TRUE)

  ld_data$chr_bp_A <- paste(ld_data$CHR_A, ld_data$BP_A, sep = "_")
  ld_data$chr_bp_B <- paste(ld_data$CHR_B, ld_data$BP_B, sep = "_")

  filtered_ld_data <- ld_data %>%
    filter(R2 < 0.6)

  length(unique(c(filtered_ld_data$chr_bp_A, filtered_ld_data$chr_bp_B)))
}

read_gemma_result <- function(assoc_file) {
  result <- read.table(assoc_file, sep = "", header = TRUE)

  colnames(result)[1] <- "CHR"
  colnames(result)[2] <- "SNP"
  colnames(result)[3] <- "BP"
  colnames(result)[13] <- "P"

  result$CHR <- beech_chr_map[result$CHR]
  result <- subset(result, CHR != "14")

  result$CHR <- as.numeric(result$CHR)
  result$BP <- as.numeric(result$BP)
  result$P <- as.numeric(result$P)

  result <- result[complete.cases(result$CHR), ]
  result$bonf <- p.adjust(result$P, method = "bonferroni")
  result$fdr <- p.adjust(result$P, method = "fdr")

  result
}

extract_selected_snps <- function(result, assoc_file, trait, covariates, num_independent_snps) {
  n <- nrow(result)
  fdr_threshold <- bh_fdr_threshold(result$P, fdr = 0.05)

  suggestive_hits <- subset(result, P < (1 / num_independent_snps))
  genome_wide_hits <- subset(result, P < (0.05 / n))
  fdr_hits <- subset(result, P < fdr_threshold)

  if (nrow(suggestive_hits) > 0) {
    suggestive_hits$significance_typ <- "suggestive"
  }

  if (nrow(genome_wide_hits) > 0) {
    genome_wide_hits$significance_typ <- "genome_wide"
  }

  if (nrow(fdr_hits) > 0) {
    fdr_hits$significance_typ <- "fdr"
  }

  selected_snps <- rbind(suggestive_hits, genome_wide_hits, fdr_hits)

  if (nrow(selected_snps) > 0) {
    selected_snps$file <- basename(assoc_file)
    selected_snps$trait <- trait
    selected_snps$covariates <- covariates
  }

  selected_snps
}

plot_gwas_result <- function(result, assoc_file, num_independent_snps) {
  chisq <- qchisq(1 - result$P, 1)
  lambda <- median(na.omit(chisq)) / qchisq(0.5, 1)
  qq_title <- paste("QQ, lambda:", lambda)

  qq_file <- paste0("plot_qq_", basename(assoc_file), ".png")
  manhattan_file <- paste0("plot_manhattan_", basename(assoc_file), ".png")

  png(file = qq_file)
  qq(result$P, main = qq_title)
  dev.off()

  png(file = manhattan_file)
  qqman::manhattan(
    result,
    main = "",
    ylim = c(0, 10),
    cex = 0.6,
    cex.axis = 0.9,
    col = c("blue4", "orange3"),
    suggestiveline = -log10(1 / num_independent_snps),
    genomewideline = -log10(0.05 / num_independent_snps)
  )
  dev.off()
}

run_selected_model <- function(model) {
  old_wd <- getwd()
  on.exit(setwd(old_wd), add = TRUE)

  ld_file <- file.path(model$ld_dir, model$ld_file)
  assoc_file <- file.path(model$assoc_dir, model$assoc_file)

  num_independent_snps <- independent_snp_count(ld_file)
  result <- read_gemma_result(assoc_file)

  setwd(model$assoc_dir)

  plot_gwas_result(
    result = result,
    assoc_file = model$assoc_file,
    num_independent_snps = num_independent_snps
  )

  extract_selected_snps(
    result = result,
    assoc_file = model$assoc_file,
    trait = model$trait,
    covariates = model$covariates,
    num_independent_snps = num_independent_snps
  )
}

write_selected_snps <- function(selected_snps_list, output_file) {
  non_empty_snps <- Filter(
    Negate(is.null),
    lapply(selected_snps_list, function(x) if (nrow(x) > 0) x else NULL)
  )

  selected_snps <- do.call(rbind, non_empty_snps)
  write.csv(selected_snps, file = output_file)
}

models <- list(
  list(
    id = "selected_snp5",
    output_group = "bai",
    ld_dir = gwas_traits_dir,
    assoc_dir = gwas_traits_dir,
    ld_file = "filtered_beech_bai_gwas.recoderesiduals_bai_vs_site_and_block_lmm_with_cov_pc17_thinned_relatedness.assoc.txt.recode.vcf.ld",
    assoc_file = "beech_bai_gwas.recoderesiduals_bai_vs_site_and_ci_basum_and_block_lmm_with_cov_pc17_thinned_relatedness.assoc.txt",
    trait = "bai",
    covariates = "vs_site_and_competition_and_block_cov_pc17"
  ),
  list(
    id = "selected_snp38",
    output_group = "dbh",
    ld_dir = gwas_traits_dir,
    assoc_dir = gwas_traits_dir,
    ld_file = "filtered_beech_dbh_gwas.recoderesiduals_dbh_vs_site_lmm_with_cov_pc1_thinned_relatedness.assoc.txt.recode.vcf.ld",
    assoc_file = "beech_dbh_gwas.recoderesiduals_dbh_vs_site_lmm_with_cov_pc4_thinned_relatedness.assoc.txt",
    trait = "dbh",
    covariates = "vs_site_cov_pc4"
  ),
  list(
    id = "selected_snp15",
    output_group = "height",
    ld_dir = gwas_traits_dir,
    assoc_dir = gwas_traits_dir,
    ld_file = "filtered_beech_height_gwas.recoderesiduals_height_vs_site_and_ci_basum_and_block_lmm_with_cov_pc2_thinned_relatedness.assoc.txt.recode.vcf.ld",
    assoc_file = "beech_height_gwas.recoderesiduals_height_vs_site_and_ci_basum_and_block_lmm_with_cov_pc1_thinned_relatedness.assoc.txt",
    trait = "height",
    covariates = "vs_site_and_competition_and_block_cov_pc1"
  ),
  list(
    id = "selected_snp30",
    output_group = "sla",
    ld_dir = gwas_traits_dir,
    assoc_dir = gwas_traits_dir,
    ld_file = "filtered_beech_sla_gwas.recoderesiduals_sla_vs_site_and_ci_basum_lmm_with_cov_pc1_thinned_relatedness.assoc.txt.recode.vcf.ld",
    assoc_file = "beech_sla_gwas.recoderesiduals_sla_vs_site_and_ci_basum_lmm_with_cov_pc6_thinned_relatedness.assoc.txt",
    trait = "sla",
    covariates = "vs_site_and_competition_cov_pc6"
  ),
  list(
    id = "selected_snp9",
    output_group = "rt",
    ld_dir = gwas_traits_dir,
    assoc_dir = gwas_traits_dir,
    ld_file = "filtered_beech_rt_gwas.recoderesiduals_rt_vs_site_and_block_lmm_with_cov_pc3_thinned_relatedness.assoc.txt.recode.vcf.ld",
    assoc_file = "beech_rt_gwas.recoderesiduals_rt_vs_site_and_block_lmm_with_cov_pc5_thinned_relatedness.assoc.txt",
    trait = "rt",
    covariates = "vs_site_and_block_cov_pc5"
  ),
  list(
    id = "selected_snp36",
    output_group = "rs",
    ld_dir = gwas_traits_dir,
    assoc_dir = gwas_traits_dir,
    ld_file = "filtered_beech_rs_gwas.recoderesiduals_rs_vs_site_lmm_with_cov_pc1_thinned_relatedness.assoc.txt.recode.vcf.ld",
    assoc_file = "beech_rs_gwas.recoderesiduals_rs_vs_site_lmm_with_cov_pc2_thinned_relatedness.assoc.txt",
    trait = "rs",
    covariates = "vs_site_cov_pc2"
  ),
  list(
    id = "selected_snp31",
    output_group = "rc",
    ld_dir = gwas_traits_dir,
    assoc_dir = gwas_traits_dir,
    ld_file = "filtered_beech_rc_gwas.recoderesiduals_rc_vs_site_and_ci_basum_lmm_with_cov_pc1_thinned_relatedness.assoc.txt.recode.vcf.ld",
    assoc_file = "beech_rc_gwas.recoderesiduals_rc_vs_site_and_ci_basum_lmm_with_cov_pc7_thinned_relatedness.assoc.txt",
    trait = "rc",
    covariates = "vs_site_and_competition_cov_pc7"
  ),
  list(
    id = "selected_snp14",
    output_group = "corr",
    ld_dir = gwas_traits_dir,
    assoc_dir = corr_recheck_dir,
    ld_file = "filtered_beech_tmean_mean_corr_gwas.recoderesiduals_corr_tmean_mean_spr_vs_site_and_ci_baspr_and_block_lmm_with_thinned_relatedness.assoc.txt.recode.vcf.ld",
    assoc_file = "beech_tmean_mean_corr_141124_gwas.recoderesiduals_corr_tmean_mean_sum_vs_site_and_ci_basum_lmm_with_thinned_relatedness_with_cov_pc9_thinned_relatedness_181124.assoc.txt",
    trait = "corr_tmean_mean_sum",
    covariates = "vs_site_ci_basum_and_pc9"
  ),
  list(
    id = "selected_snp51",
    output_group = "corr",
    ld_dir = gwas_traits_dir,
    assoc_dir = corr_recheck_dir,
    ld_file = "filtered_beech_tmean_mean_corr_gwas.recoderesiduals_corr_tmean_mean_spr_vs_site_and_ci_baspr_and_block_lmm_with_thinned_relatedness.assoc.txt.recode.vcf.ld",
    assoc_file = "beech_tmean_mean_corr_141124_gwas.recoderesiduals_corr_tmean_mean_spr_vs_site_and_ci_basum_lmm_with_thinned_relatedness_with_cov_pc5_thinned_relatedness_181124.assoc.txt",
    trait = "corr_tmean_mean_spr",
    covariates = "vs_site_ci_basum_and_pc5"
  ),
  list(
    id = "selected_snp80",
    output_group = "corr",
    ld_dir = gwas_traits_dir,
    assoc_dir = corr_recheck_dir,
    ld_file = "filtered_beech_prec_mean_corr_gwas.recoderesiduals_corr_prec_mean_spr_vs_site_and_ci_baspr_and_block_lmm_with_thinned_relatedness.assoc.txt.recode.vcf.ld",
    assoc_file = "beech_prec_mean_corr_141124_gwas.recoderesiduals_corr_prec_mean_sum_vs_site_and_ci_basum_lmm_with_thinned_relatedness_with_cov_pc9_thinned_relatedness_181124.assoc.txt",
    trait = "corr_prec_mean_sum",
    covariates = "vs_site_ci_basum_and_pc9"
  ),
  list(
    id = "selected_snp129",
    output_group = "corr",
    ld_dir = gwas_traits_dir,
    assoc_dir = corr_recheck_dir,
    ld_file = "filtered_beech_prec_mean_corr_gwas.recoderesiduals_corr_prec_mean_spr_vs_site_and_ci_baspr_and_block_lmm_with_thinned_relatedness.assoc.txt.recode.vcf.ld",
    assoc_file = "beech_prec_mean_corr_141124_gwas.recoderesiduals_corr_prec_mean_spr_vs_site_and_ci_basum_and_block_lmm_with_thinned_relatedness_with_cov_pc4_thinned_relatedness_181124.assoc.txt",
    trait = "corr_prec_mean_spr",
    covariates = "vs_site_block_ci_basum_and_pc4"
  )
)

selected_snps_by_model <- lapply(models, run_selected_model)
names(selected_snps_by_model) <- vapply(models, `[[`, character(1), "id")

# Keep the original selected-SNP object names for transparency/reproducibility.
list2env(selected_snps_by_model, envir = .GlobalEnv)

write.csv(selected_snp5, file = file.path(gwas_traits_dir, "bai_beech_pc17.csv"))

output_files <- list(
  dbh = file.path(gwas_traits_dir, "beech_gwas_results_dbh_270624.csv"),
  height = file.path(gwas_traits_dir, "beech_gwas_results_height_270624.csv"),
  sla = file.path(gwas_traits_dir, "beech_gwas_results_sla_270624.csv"),
  rt = file.path(gwas_traits_dir, "beech_gwas_results_rt_270624.csv"),
  rs = file.path(gwas_traits_dir, "beech_gwas_results_rs_270624.csv"),
  rc = file.path(gwas_traits_dir, "beech_gwas_results_rc_270624.csv"),
  corr = file.path(corr_recheck_dir, "beech_gwas_results_corr_trait_191124.csv")
)

for (group_name in names(output_files)) {
  group_models <- vapply(models, `[[`, character(1), "output_group") == group_name
  write_selected_snps(
    selected_snps_list = selected_snps_by_model[group_models],
    output_file = output_files[[group_name]]
  )
}
