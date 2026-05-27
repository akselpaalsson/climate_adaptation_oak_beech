# Oak GWAS LMM: Manhattan plots, QQ plots, and selected SNP exports

rm(list = ls())
gc()

library(dplyr)
library(qqman)
library(stringr)

setwd("path/to/wd")

benjamini_hochberg_cutoff <- function(p_values, fdr = 0.05) {
  p_sorted <- sort(p_values, decreasing = FALSE)
  m <- length(p_sorted)
  cutoffs <- (seq_len(m) / m) * fdr
  k <- max(c(which(p_sorted <= cutoffs), 0))
  ((0:m) / m * fdr)[k + 1]
}

count_independent_snps <- function(ld_file, r2_threshold = 0.6) {
  ld_data <- read.table(ld_file, sep = "", header = TRUE)

  ld_data$chr_bp_A <- paste(ld_data$CHR_A, ld_data$BP_A, sep = "_")
  ld_data$chr_bp_B <- paste(ld_data$CHR_B, ld_data$BP_B, sep = "_")

  filtered_ld_data <- ld_data %>%
    filter(R2 < r2_threshold)

  length(unique(c(filtered_ld_data$chr_bp_A, filtered_ld_data$chr_bp_B)))
}

prepare_gemma_results <- function(assoc_file) {
  result <- read.table(assoc_file, sep = "", header = TRUE)

  colnames(result)[1] <- "CHR"
  colnames(result)[2] <- "SNP"
  colnames(result)[3] <- "BP"
  colnames(result)[13] <- "P"

  chromosome_map <- c(
    "OW028765.1" = 1, "OW028766.1" = 2, "OW028767.1" = 3,
    "OW028768.1" = 4, "OW028769.1" = 5, "OW028773.1" = 6,
    "OW028770.1" = 7, "OW028771.1" = 8, "OW028772.1" = 9,
    "OW028774.1" = 10, "OW028775.1" = 11, "OW028776.1" = 12
  )

  result$CHR <- chromosome_map[result$CHR]
  result$CHR <- as.numeric(result$CHR)
  result$BP <- as.numeric(result$BP)
  result$P <- as.numeric(result$P)

  result <- result[complete.cases(result$CHR), ]
  result$bonf <- p.adjust(result$P, method = "bonferroni")
  result$fdr <- p.adjust(result$P, method = "fdr")

  result
}

get_selected_snps <- function(result, assoc_file, trait, covariates, group, num_independent_snps) {
  n <- nrow(result)
  fdr_cutoff <- benjamini_hochberg_cutoff(result$P, fdr = 0.05)

  suggestive_hits <- subset(result, P < (1 / num_independent_snps))
  genome_wide_hits <- subset(result, P < (0.05 / n))
  fdr_hits <- subset(result, P < fdr_cutoff)

  if (nrow(suggestive_hits) > 0) {
    suggestive_hits$significance_typ <- "suggestive"
  }

  if (nrow(genome_wide_hits) > 0) {
    genome_wide_hits$significance_typ <- "genome_wide"
  }

  if (nrow(fdr_hits) > 0) {
    fdr_hits$significance_typ <- "fdr"
  }

  selected <- rbind(suggestive_hits, genome_wide_hits, fdr_hits)

  if (nrow(selected) > 0) {
    selected$file <- assoc_file
    selected$trait <- trait
    selected$covariates <- covariates
    selected$group <- group
  }

  selected
}

plot_gwas_results <- function(result, qq_file, manhattan_file, num_independent_snps) {
  chisq <- qchisq(1 - result$P, 1)
  lambda <- median(na.omit(chisq)) / qchisq(0.5, 1)

  png(file = qq_file)
  qq(result$P, main = paste("QQ, lambda: ", lambda, sep = ""))
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

run_gwas_model <- function(model) {
  num_independent_snps <- count_independent_snps(model$ld_file)
  result <- prepare_gemma_results(model$assoc_file)

  plot_gwas_results(
    result = result,
    qq_file = model$qq_file,
    manhattan_file = model$manhattan_file,
    num_independent_snps = num_independent_snps
  )

  get_selected_snps(
    result = result,
    assoc_file = model$assoc_file,
    trait = model$trait,
    covariates = model$covariates,
    group = model$group,
    num_independent_snps = num_independent_snps
  )
}

models <- list(
  list(
    trait = "prec_mean_spr",
    covariates = "vs_site",
    group = "full",
    ld_file = "filtered_gemma_prec_spr_populations_no_no277_reps_missing_removed_only_petrea_minDP3_meanDP10_maxDP80_NA0.90.MAF0.05.recode.vcf.ld",
    assoc_file = "gemma_prec_spr_populations_no_no277_reps_missing_removed_only_petrea_minDP3_meanDP10_maxDP80_NA0.90.MAF0.05.recoderesiduals_prec_spr_vs_site_lmm_with_thinned_relatedness.assoc.txt",
    qq_file = "plot_qq_gemma_prec_spr_populations_no_no277_reps_missing_removed_only_petrea_minDP3_meanDP10_maxDP80_NA0.90.MAF0.05.recoderesiduals_prec_spr_vs_site_lmm_with_thinned_relatedness.assoc.txt.png",
    manhattan_file = "plot_manhattan_gemma_prec_spr_populations_no_no277_reps_missing_removed_only_petrea_minDP3_meanDP10_maxDP80_NA0.90.MAF0.05.recoderesiduals_prec_spr_vs_site_lmm_with_thinned_relatedness.assoc.txt.png"
  ),
  list(
    trait = "prec_mean_sum",
    covariates = "vs_site",
    group = "full",
    ld_file = "filtered_gemma_prec_sum_populations_no_no277_reps_missing_removed_only_petrea_minDP3_meanDP10_maxDP80_NA0.90.MAF0.05.recode.vcf.ld",
    assoc_file = "gemma_prec_sum_populations_no_no277_reps_missing_removed_only_petrea_minDP3_meanDP10_maxDP80_NA0.90.MAF0.05.recoderesiduals_prec_sum_vs_site_lmm_with_thinned_relatedness.assoc.txt",
    qq_file = "plot_qq_gemma_prec_sum_populations_no_no277_reps_missing_removed_only_petrea_minDP3_meanDP10_maxDP80_NA0.90.MAF0.05.recoderesiduals_prec_sum_vs_site_lmm_with_thinned_relatedness.assoc.txt.png",
    manhattan_file = "plot_manhattan_gemma_prec_sum_populations_no_no277_reps_missing_removed_only_petrea_minDP3_meanDP10_maxDP80_NA0.90.MAF0.05.recoderesiduals_prec_sum_vs_site_lmm_with_thinned_relatedness.assoc.txt.png"
  ),
  list(
    trait = "tmean_spr",
    covariates = "vs_site_and_competition_and_block",
    group = "full",
    ld_file = "filtered_gemma_tmean_mean_populations_no_no277_reps_missing_removed_only_petrea_minDP3_meanDP10_maxDP80_NA0.90.MAF0.05.recode.vcf.ld",
    assoc_file = "gemma_tmean_spr_populations_no_no277_reps_missing_removed_only_petrea_minDP3_meanDP10_maxDP80_NA0.90.MAF0.05.recoderesiduals_tmean_spr_vs_site_and_block_and_competition_lmm_with_thinned_relatedness.assoc.txt",
    qq_file = "plot_qq_gemma_tmean_spr_populations_no_no277_reps_missing_removed_only_petrea_minDP3_meanDP10_maxDP80_NA0.90.MAF0.05.recoderesiduals_tmean_spr_vs_site_and_block_and_competition_lmm_with_thinned_relatedness.assoc.txt.png",
    manhattan_file = "plot_manhattan_gemma_tmean_spr_populations_no_no277_reps_missing_removed_only_petrea_minDP3_meanDP10_maxDP80_NA0.90.MAF0.05.recoderesiduals_tmean_spr_vs_site_and_block_and_competition_lmm_with_thinned_relatedness.assoc.txt.png"
  ),
  list(
    trait = "tmean_sum",
    covariates = "vs_site_and_competition_and_block",
    group = "full",
    ld_file = "filtered_gemma_tmean_sum_populations_no_no277_reps_missing_removed_only_petrea_minDP3_meanDP10_maxDP80_NA0.90.MAF0.05.recode.vcf.ld",
    assoc_file = "gemma_tmean_sum_populations_no_no277_reps_missing_removed_only_petrea_minDP3_meanDP10_maxDP80_NA0.90.MAF0.05.recoderesiduals_tmean_sum_vs_site_and_block_and_competition_lmm_with_thinned_relatedness.assoc.txt",
    qq_file = "plot_qq_gemma_tmean_sum_populations_no_no277_reps_missing_removed_only_petrea_minDP3_meanDP10_maxDP80_NA0.90.MAF0.05.recoderesiduals_tmean_sum_vs_site_and_block_and_competition_lmm_with_thinned_relatedness.assoc.txt.png",
    manhattan_file = "plot_manhattan_gemma_tmean_sum_populations_no_no277_reps_missing_removed_only_petrea_minDP3_meanDP10_maxDP80_NA0.90.MAF0.05.recoderesiduals_tmean_sum_vs_site_and_block_and_competition_lmm_with_thinned_relatedness.assoc.txt.png"
  ),
  list(
    trait = "rc",
    covariates = "vs_site_and_block",
    group = "full",
    ld_file = "filtered_gemma_rc_populations_no_no277_reps_missing_removed_only_petrea_minDP3_meanDP10_maxDP80_NA0.90.MAF0.05.recode.vcf.ld",
    assoc_file = "gemma_rc_populations_no_no277_reps_missing_removed_only_petrea_minDP3_meanDP10_maxDP80_NA0.90.MAF0.05.recoderesiduals_rc_vs_site_and_block_lmm_with_thinned_relatedness.assoc.txt",
    qq_file = "plot_qq_gemma_rc_populations_no_no277_reps_missing_removed_only_petrea_minDP3_meanDP10_maxDP80_NA0.90.MAF0.05.recoderesiduals_rc_vs_site_and_block_lmm_with_thinned_relatedness.assoc.txt.png",
    manhattan_file = "plot_manhattan_gemma_rc_populations_no_no277_reps_missing_removed_only_petrea_minDP3_meanDP10_maxDP80_NA0.90.MAF0.05.recoderesiduals_rc_vs_site_and_block_lmm_with_thinned_relatedness.assoc.txt.png"
  ),
  list(
    trait = "bai",
    covariates = "vs_site_and_competition_and_block",
    group = "full",
    ld_file = "filtered_gemma_bai_populations_no_no277_reps_missing_removed_only_petrea_minDP3_meanDP10_maxDP80_NA0.90.MAF0.05.recode.vcf.ld",
    assoc_file = "gemma_bai_populations_no_no277_reps_missing_removed_only_petrea_minDP3_meanDP10_maxDP80_NA0.90.MAF0.05.recoderesiduals_bai_vs_site_and_block_and_competition_lmm_with_thinned_relatedness.assoc.txt",
    qq_file = "plot_qq_gemma_bai_populations_no_no277_reps_missing_removed_only_petrea_minDP3_meanDP10_maxDP80_NA0.90.MAF0.05.recoderesiduals_bai_vs_site_and_block_and_competition_lmm_with_thinned_relatedness.assoc.txt.png",
    manhattan_file = "plot_manhattan_gemma_bai_populations_no_no277_reps_missing_removed_only_petrea_minDP3_meanDP10_maxDP80_NA0.90.MAF0.05.recoderesiduals_bai_vs_site_and_block_and_competition_lmm_with_thinned_relatedness.assoc.txt.png"
  ),
  list(
    trait = "rt",
    covariates = "vs_site_and_competition",
    group = "full",
    ld_file = "filtered_gemma_rt_populations_no_no277_reps_missing_removed_only_petrea_minDP3_meanDP10_maxDP80_NA0.90.MAF0.05.recode.vcf.ld",
    assoc_file = "gemma_rt_populations_no_no277_reps_missing_removed_only_petrea_minDP3_meanDP10_maxDP80_NA0.90.MAF0.05.recoderesiduals_rt_vs_site_and_competition_lmm_with_thinned_relatedness.assoc.txt",
    qq_file = "plot_qq_gemma_rt_populations_no_no277_reps_missing_removed_only_petrea_minDP3_meanDP10_maxDP80_NA0.90.MAF0.05.recoderesiduals_rt_vs_site_and_competition_lmm_with_thinned_relatedness.assoc.txt.png",
    manhattan_file = "plot_manhattan_gemma_rt_populations_no_no277_reps_missing_removed_only_petrea_minDP3_meanDP10_maxDP80_NA0.90.MAF0.05.recoderesiduals_rt_vs_site_and_competition_lmm_with_thinned_relatedness.assoc.txt.png"
  ),
  list(
    trait = "rs",
    covariates = "vs_site_and_competition_and_block",
    group = "full",
    ld_file = "filtered_gemma_rs_populations_no_no277_reps_missing_removed_only_petrea_minDP3_meanDP10_maxDP80_NA0.90.MAF0.05.recode.vcf.ld",
    assoc_file = "gemma_rs_populations_no_no277_reps_missing_removed_only_petrea_minDP3_meanDP10_maxDP80_NA0.90.MAF0.05.recoderesiduals_rs_vs_site_and_block_and_competition_lmm_with_thinned_relatedness.assoc.txt",
    qq_file = "plot_qq_gemma_rs_populations_no_no277_reps_missing_removed_only_petrea_minDP3_meanDP10_maxDP80_NA0.90.MAF0.05.recoderesiduals_rs_vs_site_and_block_and_competition_lmm_with_thinned_relatedness.assoc.txt.png",
    manhattan_file = "plot_manhattan_gemma_rs_populations_no_no277_reps_missing_removed_only_petrea_minDP3_meanDP10_maxDP80_NA0.90.MAF0.05.recoderesiduals_rs_vs_site_and_block_and_competition_lmm_with_thinned_relatedness.assoc.txt.png"
  ),
  list(
    trait = "sla",
    covariates = "vs_site_and_competition",
    group = "full",
    ld_file = "filtered_gemma_sla_populations_no_no277_reps_missing_removed_only_petrea_minDP3_meanDP10_maxDP80_NA0.90.MAF0.05.recode.vcf.ld",
    assoc_file = "gemma_sla_populations_no_no277_reps_missing_removed_only_petrea_minDP3_meanDP10_maxDP80_NA0.90.MAF0.05.recoderesiduals_sla_vs_site_and_competition_lmm_with_thinned_relatedness.assoc.txt",
    qq_file = "plot_qq_gemma_sla_populations_no_no277_reps_missing_removed_only_petrea_minDP3_meanDP10_maxDP80_NA0.90.MAF0.05.recoderesiduals_sla_vs_site_and_competition_lmm_with_thinned_relatedness.assoc.txt.png",
    manhattan_file = "plot_manhattan_gemma_sla_populations_no_no277_reps_missing_removed_only_petrea_minDP3_meanDP10_maxDP80_NA0.90.MAF0.05.recoderesiduals_sla_vs_site_and_competition_lmm_with_thinned_relatedness.assoc.txt.png"
  ),
  list(
    trait = "height",
    covariates = "vs_site",
    group = "full",
    ld_file = "filtered_gemma_dbh_populations_no_no277_reps_missing_removed_only_petrea_minDP3_meanDP10_maxDP80_NA0.90.MAF0.05.recode.vcf.ld",
    assoc_file = "gemma_height_populations_no_no277_reps_missing_removed_only_petrea_minDP3_meanDP10_maxDP80_NA0.90.MAF0.05.recoderesiduals_height_vs_site_lmm_with_thinned_relatedness.assoc.txt",
    qq_file = "plot_qq_gemma_height_populations_no_no277_reps_missing_removed_only_petrea_minDP3_meanDP10_maxDP80_NA0.90.MAF0.05.recoderesiduals_height_vs_site_lmm_with_thinned_relatedness.assoc.txt.png",
    manhattan_file = "plot_manhattan_gemma_height_populations_no_no277_reps_missing_removed_only_petrea_minDP3_meanDP10_maxDP80_NA0.90.MAF0.05.recoderesiduals_height_vs_site_lmm_with_thinned_relatedness.assoc.txt.png"
  ),
  list(
    trait = "dbh",
    covariates = "vs_site",
    group = "full",
    ld_file = "filtered_gemma_dbh_populations_no_no277_reps_missing_removed_only_petrea_minDP3_meanDP10_maxDP80_NA0.90.MAF0.05.recode.vcf.ld",
    assoc_file = "gemma_dbh_populations_no_no277_reps_missing_removed_only_petrea_minDP3_meanDP10_maxDP80_NA0.90.MAF0.05.recoderesiduals_dbh_vs_site_lmm_with_thinned_relatedness.assoc.txt",
    qq_file = "plot_qq_gemma_dbh_populations_no_no277_reps_missing_removed_only_petrea_minDP3_meanDP10_maxDP80_NA0.90.MAF0.05.recoderesiduals_dbh_vs_site_lmm_with_thinned_relatedness.assoc.txt.png",
    manhattan_file = "plot_manhattan_gemma_dbh_populations_no_no277_reps_missing_removed_only_petrea_minDP3_meanDP10_maxDP80_NA0.90.MAF0.05.recoderesiduals_dbh_vs_site_lmm_with_thinned_relatedness.assoc.txt.png"
  )
)

selected_snps <- lapply(models, run_gwas_model)
selected_snps <- selected_snps[vapply(selected_snps, nrow, integer(1)) > 0]
selected_snps_all <- do.call(rbind, selected_snps)
selected_snps_all$chr_bp <- paste(selected_snps_all$CHR, selected_snps_all$BP, sep = "_")

write.csv(
  selected_snps_all,
  file = "selected_snps_all_new_threshold_thinned_rel_mat_180624.csv"
)

# Export unique selected SNP positions for post-GWAS LD checking.
ld_chr_map <- c(
  "1" = "OW028765.1", "2" = "OW028766.1", "3" = "OW028767.1",
  "4" = "OW028768.1", "5" = "OW028768.1", "6" = "OW028768.1",
  "7" = "OW028770.1", "8" = "OW028771.1", "9" = "OW028772.1",
  "10" = "OW028774.1", "11" = "OW028775.1", "12" = "OW028776.1"
)

ld_check_snps <- selected_snps_all[, c("CHR", "BP", "chr_bp")]
ld_check_snps$chr <- ld_chr_map[ld_check_snps$CHR]
ld_check_snps <- ld_check_snps[order(ld_check_snps$chr_bp), ]
ld_check_snps <- ld_check_snps[!duplicated(ld_check_snps$chr_bp), ]
ld_check_snps <- ld_check_snps[, c("chr", "BP")]

write.table(
  ld_check_snps,
  file = "unique_snps_to_check_ld_thinned_rel_mat_180624.txt",
  quote = FALSE,
  row.names = FALSE
)

# Prepare unique LD-pruned hits for NCBI genome-coordinate annotation.
ld_pruned_hits <- read.csv("final_snp_ld_pruned_gwas_res_oak_200624.csv")

summary_table <- ld_pruned_hits %>%
  group_by(group, trait) %>%
  summarise(unique_entries = n_distinct(chr_bp), .groups = "drop")

unique_entries <- as.data.frame(unique(ld_pruned_hits$chr_bp))
colnames(unique_entries)[1] <- "chr_bp"

split_chr_bp <- str_split(unique_entries$chr_bp, "_")
unique_entries$CHR <- sapply(split_chr_bp, "[", 1)
unique_entries$BP <- as.numeric(sapply(split_chr_bp, "[", 2))
unique_entries <- subset(unique_entries, select = -chr_bp)

ncbi_chr_map <- c(
  "1" = "NC_065534.1", "2" = "NC_065535.1", "3" = "NC_065536.1",
  "4" = "NC_065537.1", "5" = "NC_065538.1", "6" = "NC_065539.1",
  "7" = "NC_065540.1", "8" = "NC_065541.1", "9" = "NC_065542.1",
  "10" = "NC_065543.1", "11" = "NC_065544.1", "12" = "NC_065545.1"
)

unique_entries$CHR <- ncbi_chr_map[unique_entries$CHR]
colnames(unique_entries) <- NULL

write.table(
  unique_entries,
  file = "gwas_hits_for_ncbi_genome_200624.txt",
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

