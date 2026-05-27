

rm(list = ls())
gc()

library(ggplot2)

setwd("/path/to/working_directory")

beech_chromosome_map <- c(
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









bslmm_models <- list(
  list(trait = "bai", param_file = "beech_bai_gwas_chain1.recoderesiduals_40m__bai_vs_site_and_ci_basum_and_block_bslmm_with_cov_pc17_thinned_relatedness_chain1.param.txt", y_max = 0.25, pip01_file = "beech_bai_pip01.csv", plot_file = "plot_bslmm_beech_bai_pip.png"),
  list(trait = "dbh", param_file = "beech_dbh_gwas_chain1.recoderesiduals_40m__dbh_vs_site_bslmm_with_cov_pc4_thinned_relatedness_chain1.param.txt", y_max = 0.25, pip01_file = "beech_dbh_pip01.csv", plot_file = "plot_bslmm_beech_dbh_pip.png"),
  list(trait = "height", param_file = "beech_height_gwas_chain1.recoderesiduals_40m__height_vs_site_and_ci_basum_and_block_bslmm_with_cov_pc1_thinned_relatedness_chain1.param.txt", y_max = 0.25, pip01_file = "beech_height_pip01.csv", plot_file = "plot_bslmm_beech_height_pip.png"),
  list(trait = "prec_spr", param_file = "beech_prec_mean_corr_141124_gwas_chain1_with_pcs.recoderesiduals_40m__beech_corr_prec_mean_spr_residual_traits_with_pcs_reprep_thinned_relatedness_chain1_with_pcs.param.txt", y_max = 0.25, pip01_file = "beech_prec_spr_pip01.csv", plot_file = "plot_bslmm_beech_prec_spr_pip.png"),
  list(trait = "prec_sum", param_file = "beech_prec_mean_corr_141124_gwas_chain1_with_pcs.recoderesiduals_40m__beech_corr_prec_mean_sum_residual_traits_with_pcs_reprep_thinned_relatedness_chain1_with_pcs.param.txt", y_max = 0.25, pip01_file = "beech_prec_sum_pip01.csv", plot_file = "plot_bslmm_beech_prec_sum_pip.png"),
  list(trait = "rc", param_file = "beech_rc_gwas_chain1.recoderesiduals_rc_vs_site_and_ci_basum_bslmm_with_cov_pc7_thinned_relatedness_chain1.param.txt", y_max = 0.25, pip01_file = "beech_rc_pip01.csv", plot_file = "plot_bslmm_beech_rc_pip.png"),
  list(trait = "rs", param_file = "beech_rs_gwas_chain1.recoderesiduals_rs_vs_site_bslmm_with_cov_pc2_thinned_relatedness_chain1.param.txt", y_max = 0.25, pip01_file = "beech_rs_pip01.csv", plot_file = "plot_bslmm_beech_rs_pip.png"),
  list(trait = "rt", param_file = "beech_rt_gwas_chain1.recoderesiduals_rt_vs_site_and_block_bslmm_with_cov_pc5_thinned_relatedness_chain1.param.txt", y_max = 0.25, pip01_file = "beech_rt_pip01.csv", plot_file = "plot_bslmm_beech_rt_pip.png"),
  list(trait = "sla", param_file = "beech_sla_gwas_chain1.recoderesiduals_40m__sla_vs_site_and_ci_basum_and_block_bslmm_with_cov_pc6_thinned_relatedness_chain1.param.txt", y_max = 0.25, pip01_file = "beech_sla_pip01.csv", plot_file = "plot_bslmm_beech_sla_pip.png"),
  list(trait = "tmean_spr", param_file = "beech_tmean_mean_corr_141124_gwas_chain1_with_pcs.recoderesiduals_40m__beech_corr_tmean_mean_spr_residual_traits_with_pcs_reprep_thinned_relatedness_chain1_with_pcs.param.txt", y_max = 0.25, pip01_file = "beech_tmean_spr_pip01.csv", plot_file = "plot_bslmm_beech_tmean_spr_pip.png"),
  list(trait = "tmean_sum", param_file = "beech_tmean_mean_corr_141124_gwas_chain2_with_pcs.recoderesiduals_40m__beech_corr_tmean_mean_sum_residual_traits_with_pcs_reprep_thinned_relatedness_chain2_with_pcs.param.txt", y_max = 0.25, pip01_file = "beech_tmean_sum_pip01.csv", plot_file = "plot_bslmm_beech_tmean_sum_pip.png")
)

prepare_bslmm_param <- function(param_file, chromosome_map) {
  result <- read.table(param_file, header = TRUE)

  colnames(result)[1] <- "CHR"
  colnames(result)[2] <- "SNP"
  colnames(result)[3] <- "BP"
  colnames(result)[7] <- "P"

  result$CHR <- chromosome_map[result$CHR]
  result <- subset(result, CHR != 14)

  result$CHR <- as.numeric(result$CHR)
  result$BP <- as.numeric(result$BP)
  result$P <- as.numeric(result$P)

  result <- result[complete.cases(result$CHR), ]
  result <- result[order(result$CHR, result$BP), ]

  chr_offsets <- c(0, cumsum(as.numeric(tapply(result$BP, result$CHR, max))))
  result$cumBP <- result$BP

  for (chr in unique(result$CHR)) {
    result$cumBP[result$CHR == chr] <- result$BP[result$CHR == chr] + chr_offsets[chr]
  }

  result
}

plot_bslmm_pip <- function(result, y_max, trait_label) {
  ggplot(result, aes(x = cumBP, y = P, color = factor(CHR))) +
    geom_point(alpha = 0.75, size = 1.2) +
    scale_color_manual(values = rep(c("darkblue", "darkred"), 22)) +
    scale_x_continuous(
      breaks = tapply(result$cumBP, result$CHR, mean),
      labels = unique(result$CHR)
    ) +
    coord_cartesian(ylim = c(0, y_max)) +
    geom_hline(yintercept = 0.01, linetype = "solid", color = "black", linewidth = 0.5) +
    labs(x = "Chromosome", y = "PIP", title = trait_label) +
    theme_minimal() +
    theme(
      legend.position = "none",
      axis.text.x = element_text(angle = 0, hjust = 1, face = "bold", size = 12),
      axis.text.y = element_text(face = "bold", size = 12),
      axis.title = element_text(face = "bold", size = 14),
      panel.grid = element_blank(),
      panel.border = element_blank(),
      axis.line = element_line(color = "black", linewidth = 0.7),
      axis.ticks = element_line(color = "black", linewidth = 0.7)
    )
}

process_bslmm_model <- function(model, chromosome_map) {
  result <- prepare_bslmm_param(model$param_file, chromosome_map)

  pip01 <- result[result$P >= 0.01, ]
  write.csv(pip01, file = model$pip01_file)

  plot <- plot_bslmm_pip(
    result = result,
    y_max = model$y_max,
    trait_label = paste("Beech", model$trait)
  )

  ggsave(
    filename = model$plot_file,
    plot = plot,
    width = 10,
    height = 5,
    dpi = 300
  )

  invisible(pip01)
}

pip01_results <- lapply(
  bslmm_models,
  process_bslmm_model,
  chromosome_map = beech_chromosome_map
)

names(pip01_results) <- vapply(bslmm_models, `[[`, character(1), "trait")
