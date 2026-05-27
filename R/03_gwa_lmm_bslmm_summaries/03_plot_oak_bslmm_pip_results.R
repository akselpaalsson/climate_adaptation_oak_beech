
rm(list = ls())
gc()

library(ggplot2)


setwd("/path/to/working_directory")

oak_chromosome_map <- c(
  "OW028765.1" = 1,
  "OW028766.1" = 2,
  "OW028767.1" = 3,
  "OW028768.1" = 4,
  "OW028769.1" = 5,
  "OW028773.1" = 6,
  "OW028770.1" = 7,
  "OW028771.1" = 8,
  "OW028772.1" = 9,
  "OW028774.1" = 10,
  "OW028775.1" = 11,
  "OW028776.1" = 12
)

bslmm_models <- list(
  list(trait = "dbh", param_file = "gemma_dbh_populations_no_no277_reps_missing_removed_only_petrea_minDP3_meanDP10_maxDP80_NA0.90.MAF0.05.recoderesiduals_dbh_vs_site_bslmm_with_thinned_relatedness_s_40m.param.txt", y_max = 0.08, pip01_file = "oak_dbh_pip01.csv", plot_file = "plot_bslmm_oak_dbh_pip.png"),
  list(trait = "height", param_file = "gemma_height_populations_no_no277_reps_missing_removed_only_petrea_minDP3_meanDP10_maxDP80_NA0.90.MAF0.05.recoderesiduals_height_vs_site_lmm_with_thinned_relatedness_s_40m.param.txt", y_max = 0.5, pip01_file = "oak_height_pip01.csv", plot_file = "plot_bslmm_oak_height_pip.png"),
  list(trait = "bai", param_file = "gemma_bai_populations_no_no277_reps_missing_removed_only_petrea_minDP3_meanDP10_maxDP80_NA0.90.MAF0.05.recoderesiduals_bai_vs_site_and_block_and_competition_lmm_with_thinned_relatedness_s_40m.param.txt", y_max = 0.75, pip01_file = "oak_bai_pip01.csv", plot_file = "plot_bslmm_oak_bai_pip.png"),
  list(trait = "sla", param_file = "gemma_sla_populations_no_no277_reps_missing_removed_only_petrea_minDP3_meanDP10_maxDP80_NA0.90.MAF0.05.recoderesiduals_sla_vs_site_and_competition_lmm_with_thinned_relatedness_s_40m.param.txt", y_max = 0.2, pip01_file = "oak_sla_pip01.csv", plot_file = "plot_bslmm_oak_sla_pip.png"),
  list(trait = "sprpre", param_file = "gemma_prec_spr_populations_no_no277_reps_missing_removed_only_petrea_minDP3_meanDP10_maxDP80_NA0.90.MAF0.05.recoderesiduals_prec_spr_vs_site_lmm_with_thinned_relatedness_s_40m.param.txt", y_max = 0.15, pip01_file = "oak_sprpre_pip01.csv", plot_file = "plot_bslmm_oak_sprpre_pip.png"),
  list(trait = "sumpre", param_file = "gemma_prec_sum_populations_no_no277_reps_missing_removed_only_petrea_minDP3_meanDP10_maxDP80_NA0.90.MAF0.05.recoderesiduals_prec_sum_vs_site_lmm_with_thinned_relatedness_s_40m.param.txt", y_max = 0.1, pip01_file = "oak_sumpre_pip01.csv", plot_file = "plot_bslmm_oak_sumpre_pip.png"),
  list(trait = "sprtmea", param_file = "gemma_tmean_spr_populations_no_no277_reps_missing_removed_only_petrea_minDP3_meanDP10_maxDP80_NA0.90.MAF0.05.recoderesiduals_tmean_spr_vs_site_and_block_and_competition_lmm_with_thinned_relatedness_s_40m.param.txt", y_max = 0.25, pip01_file = "oak_sprtmea_pip01.csv", plot_file = "plot_bslmm_oak_sprtmea_pip.png"),
  list(trait = "sumtmea", param_file = "gemma_tmean_sum_populations_no_no277_reps_missing_removed_only_petrea_minDP3_meanDP10_maxDP80_NA0.90.MAF0.05.recoderesiduals_tmean_sum_vs_site_and_block_and_competition_lmm_with_thinned_relatedness_s_40m.param.txt", y_max = 0.15, pip01_file = "oak_sumtmea_pip01.csv", plot_file = "plot_bslmm_oak_sumtmea_pip.png"),
  list(trait = "rc", param_file = "gemma_rc_populations_no_no277_reps_missing_removed_only_petrea_minDP3_meanDP10_maxDP80_NA0.90.MAF0.05.recoderesiduals_rc_vs_site_and_block_lmm_with_thinned_relatedness_s_40m.param.txt", y_max = 0.1, pip01_file = "oak_rc_pip01.csv", plot_file = "plot_bslmm_oak_rc_pip.png"),
  list(trait = "rs", param_file = "gemma_rs_populations_no_no277_reps_missing_removed_only_petrea_minDP3_meanDP10_maxDP80_NA0.90.MAF0.05.recoderesiduals_rs_vs_site_and_block_and_competition_lmm_with_thinned_relatedness_s_40m.param.txt", y_max = 0.2, pip01_file = "oak_rs_pip01.csv", plot_file = "plot_bslmm_oak_rs_pip.png"),
  list(trait = "rt", param_file = "gemma_rt_populations_no_no277_reps_missing_removed_only_petrea_minDP3_meanDP10_maxDP80_NA0.90.MAF0.05.recoderesiduals_rt_vs_site_and_competition_lmm_with_thinned_relatedness_s_40m.param.txt", y_max = 0.15, pip01_file = "oak_rt_pip01.csv", plot_file = "plot_bslmm_oak_rt_pip.png")
)

prepare_bslmm_param <- function(param_file, chromosome_map) {
  result <- read.table(param_file, sep = "\t", header = TRUE)

  colnames(result)[1] <- "CHR"
  colnames(result)[2] <- "SNP"
  colnames(result)[3] <- "BP"
  colnames(result)[7] <- "P"

  result$CHR <- chromosome_map[result$CHR]
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
    geom_point(alpha = 1, size = 2) +
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
    trait_label = paste("Oak", model$trait)
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
  chromosome_map = oak_chromosome_map
)

names(pip01_results) <- vapply(bslmm_models, `[[`, character(1), "trait")
