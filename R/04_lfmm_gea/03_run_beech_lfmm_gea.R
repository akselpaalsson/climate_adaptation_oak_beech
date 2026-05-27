
rm(list = ls())
gc()

library(LEA)
library(adegenet)
library(vcfR)
library(qvalue)
library(viridis)
library(scales)
library(tidyr)
library(dplyr)
library(ggplot2)

# -----------------------------------------------------------------------------
# Paths and input files
# -----------------------------------------------------------------------------

setwd("/path/to/working_directory")


dir_path <- file.path(test_dir, "res")
dat_path <- test_dir
env_path <- test_dir

dir.create(dir_path, showWarnings = FALSE, recursive = TRUE)

vcf_file <- file.path(
  dat_path,
  "populations_fagus_no_high_miss_ind_minDP3_meanDP10_maxDP80_NA0.90_reduced.MAF0.05_matching_evi_samples.recode.vcf"
)

lfmm_file <- file.path(
  dat_path,
  "impute_11_1beech_gen_dataset_match_eli_beech__proper_ped_MAF0.05.lfmm_imputed.lfmm"
)

env_file <- file.path(env_path, "prepped_beech_climate_data_091224.csv")

lfmm_raw_dir <- file.path(dir_path, "lfmm_raw_101224")
lfmm_dir <- file.path(dir_path, "lfmm")

univariate_summary_dir <- lfmm_raw_dir
multivariate_summary_dir <- dir_path

# -----------------------------------------------------------------------------
# Shared helpers
# -----------------------------------------------------------------------------

replace_na_with_median <- function(x) {
  x <- as.numeric(as.character(x))
  x[is.na(x)] <- median(x, na.rm = TRUE)
  x
}

write_candidate_table <- function(candidate, output_file) {
  if (NCOL(candidate) == 1) {
    write.table(
      t(candidate),
      output_file,
      sep = ",",
      row.names = FALSE,
      col.names = TRUE,
      quote = FALSE
    )
  } else if (NCOL(candidate) > 1) {
    write.table(
      candidate[order(as.numeric(candidate[, "pvalue"]), decreasing = FALSE), ],
      output_file,
      sep = ",",
      row.names = FALSE,
      col.names = TRUE,
      quote = FALSE
    )
  } else {
    write.table(
      candidate,
      output_file,
      sep = ",",
      row.names = TRUE,
      col.names = FALSE,
      quote = FALSE
    )
  }
}

write_pvalue_diagnostics <- function(pvalues, output_file) {
  pdf(output_file, width = 10, height = 6)
  par(mfrow = c(1, 2), mar = c(5, 5, 4, 1))
  hist(pvalues, col = "red", main = "P-value distribution")
  qqplot(
    rexp(length(pvalues), rate = log(10)),
    -log10(pvalues),
    xlab = "Expected quantile",
    pch = 19,
    cex = 0.4
  )
  abline(coef = c(0, 1))
  dev.off()
}

write_lfmm_manhattan_pdf <- function(pvalues, candidate, main, output_file) {
  pdf(output_file, width = 10, height = 6)
  par(mfrow = c(1, 1), mar = c(5, 5, 4, 1))
  plot(
    -log10(pvalues),
    main = main,
    cex.main = 1.2,
    xlab = "Locus",
    ylab = "-Log(P-value)",
    cex = 0.7,
    col = "grey"
  )
  points(candidate, -log10(pvalues)[candidate], pch = 19, cex = 0.7, col = "red")
  dev.off()
}

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

# -----------------------------------------------------------------------------
# Input preparation
# -----------------------------------------------------------------------------

setwd(dir_path)

vcf <- read.vcfR(vcf_file, verbose = FALSE)
snp_info <- getFIX(vcf)
snp_info <- data.frame(snp_info[, c("CHROM", "POS")])
snp_info$SNPid <- paste(snp_info$CHROM, snp_info$POS, sep = "_")

sample_ids <- colnames(vcf@gt)[-1]

gen_for_sample_order <- vcfR2genind(vcf)
gen_ind_names <- indNames(gen_for_sample_order)

env <- read.csv(env_file)
names(env)[names(env) == "X.1"] <- "X"
rownames(env) <- env$SampleID_vcf
env <- env[match(gen_ind_names, rownames(env)), ]
env <- subset(env, select = -c(X, Prov, SampleID_vcf))
env <- data.frame(apply(env, 2, replace_na_with_median))
rownames(env) <- gen_ind_names

X <- scale(env, center = TRUE, scale = TRUE)
cl <- kmeans(X, centers = 3)$cluster
Y_structure <- scaleGen(gen_for_sample_order, center = TRUE, scale = FALSE, NA.method = "mean")
gen_pc <- prcomp(Y_structure)
plot(gen_pc$x, pch = 19, col = viridis(4)[cl])

gen <- read.lfmm(lfmm_file)
rownames(gen) <- sample_ids
Y <- gen

# -----------------------------------------------------------------------------
# Univariate LFMM analyses
# -----------------------------------------------------------------------------

Ks <- 11
fdr_thres <- c(0.05, 0.01, 0.001)
fdr_output <- 0.05

dir.create(lfmm_raw_dir, recursive = TRUE, showWarnings = FALSE)

for (j in seq_len(NCOL(env))) {
  env_dir <- file.path(lfmm_raw_dir, paste0("env", j, "_", colnames(env)[j]))
  dir.create(env_dir, recursive = TRUE, showWarnings = FALSE)

  for (i in seq_len(Ks)) {
    dir.create(file.path(env_dir, paste0("K", i)), recursive = TRUE, showWarnings = FALSE)
  }
}

gif <- matrix(nrow = NCOL(env), ncol = Ks)
rownames(gif) <- colnames(env)
colnames(gif) <- paste0("K", seq_len(Ks))

nb_asso_q <- matrix(nrow = 1, ncol = length(fdr_thres))
rownames(nb_asso_q) <- "value"
colnames(nb_asso_q) <- c("q0.1", "q0.01", "q0.001")

nb_asso_k <- matrix(nrow = NCOL(env), ncol = Ks)
rownames(nb_asso_k) <- colnames(env)
colnames(nb_asso_k) <- paste0("K", seq_len(Ks))

for (j in seq_len(NCOL(env))) {
  for (i in seq_len(Ks)) {
    message(paste0("Env = ", j, " ", colnames(env)[j], " , K = ", i))

    run_dir <- file.path(lfmm_raw_dir, paste0("env", j, "_", colnames(env)[j]), paste0("K", i))
    setwd(run_dir)

    res <- matrix(nrow = NCOL(Y), ncol = 4)
    rownames(res) <- colnames(Y)
    colnames(res) <- c("SNPid", "zscore", "pvalue", "qvalue")

    mod_lfmm <- lfmm2(input = Y, env = env[, j], K = i, lambda = 1e-5)
    stats_lfmm <- lfmm2.test(
      object = mod_lfmm,
      input = Y,
      env = env[, j],
      full = FALSE,
      genomic.control = TRUE
    )

    gif[j, i] <- stats_lfmm$gif
    qv_lfmm <- qvalue::qvalue(as.vector(stats_lfmm$pvalues), fdr.level = fdr_output)

    res[, "SNPid"] <- snp_info$SNPid
    res[, "zscore"] <- stats_lfmm$zscores
    res[, "pvalue"] <- stats_lfmm$pvalues
    res[, "qvalue"] <- qv_lfmm$qvalues

    write.table(
      res,
      paste0("LFMM_AllResults_env", j, "_", colnames(env)[j], "_K", i, ".csv"),
      sep = ",",
      row.names = FALSE,
      col.names = TRUE,
      quote = FALSE
    )

    write_pvalue_diagnostics(
      stats_lfmm$pvalues,
      paste0("PvalueDistribution_env", j, "_", colnames(env)[j], "_K", i, ".pdf")
    )

    for (x in seq_along(fdr_thres)) {
      q <- fdr_thres[x]
      w <- which(sort(qv_lfmm$qvalues) <= q)
      candidate <- order(stats_lfmm$pvalues, decreasing = FALSE)[w]

      write_lfmm_manhattan_pdf(
        pvalues = stats_lfmm$pvalues,
        candidate = candidate,
        main = paste0("Manhattan plot | env", j, "_", colnames(env)[j], " | K", i, " | q=", q),
        output_file = paste0("ManhattanPlot_env", j, "_K", i, "_q", q, ".pdf")
      )

      if (q == fdr_output) {
        nb_asso_k[j, i] <- length(candidate)
      }

      nb_asso_q["value", x] <- length(candidate)
      candidate <- res[candidate, ]

      write_candidate_table(
        candidate,
        paste0("CandidatesOrdered_env", j, "_", colnames(env)[j], "_K", i, "_q", q, ".csv")
      )
    }

    write.table(
      nb_asso_q,
      paste0("AssociationsNb_env", j, "_", colnames(env)[j], "_K", i, "_q.csv"),
      sep = ",",
      row.names = FALSE,
      col.names = TRUE,
      quote = FALSE
    )
  }
}

gif_out <- cbind(rownames(gif), gif)
colnames(gif_out)[1] <- "env"
write.table(
  gif_out,
  file.path(univariate_summary_dir, paste0("GenomicInflationFactor_K1-", Ks, "_univariate.csv")),
  sep = ",",
  row.names = FALSE,
  col.names = TRUE,
  quote = FALSE
)

nb_asso_k_out <- cbind(rownames(nb_asso_k), nb_asso_k)
colnames(nb_asso_k_out)[1] <- "env"
write.table(
  nb_asso_k_out,
  file.path(univariate_summary_dir, paste0("NbSignifAssociations_K1-", Ks, "_q", fdr_output, "_univariate.csv")),
  sep = ",",
  row.names = FALSE,
  col.names = TRUE,
  quote = FALSE
)

# -----------------------------------------------------------------------------
# Plot significant univariate LFMM candidates against environmental variables
# -----------------------------------------------------------------------------

for (j in seq_len(NCOL(env))) {
  for (i in seq_len(Ks)) {
    candidate_file <- file.path(
      lfmm_raw_dir,
      paste0("env", j, "_", colnames(env)[j]),
      paste0("K", i),
      paste0("CandidatesOrdered_env", j, "_", colnames(env)[j], "_K", i, "_q", fdr_output, ".csv")
    )

    candidate_list <- read.table(candidate_file, row.names = "SNPid", header = TRUE, sep = ",")

    tmp_gen <- gen
    colnames(tmp_gen) <- snp_info$SNPid
    tmp_gen_t <- t(tmp_gen)
    rownames(tmp_gen_t) <- colnames(tmp_gen)
    colnames(tmp_gen_t) <- rownames(tmp_gen)

    candidate_gen <- merge(candidate_list, tmp_gen_t, by = "row.names")
    colnames(candidate_gen)[1] <- "SNPid"
    candidate_gen_ordered <- candidate_gen[order(candidate_gen$pvalue, decreasing = FALSE), ]

    if (NROW(candidate_list) > 0) {
      genotypes <- t(candidate_gen_ordered[, 5:NCOL(candidate_gen_ordered)])

      pdf(
        file.path(
          lfmm_raw_dir,
          paste0("env", j, "_", colnames(env)[j]),
          paste0("K", i),
          paste0("PlotOfSignifCandidates_K", i, "_q", fdr_output, ".pdf")
        ),
        width = 10,
        height = 10
      )

      par(mfrow = c(2, 2), mar = c(5, 5, 1, 1))

      for (p in seq_len(NROW(candidate_list))) {
        plot(
          env[, j],
          genotypes[, p] / 2,
          pch = 20,
          cex = 1.5,
          xlab = colnames(env)[j],
          ylab = "Genotype [-]",
          cex.lab = 1.2,
          ylim = c(0, 1),
          cex.main = 1.2,
          col = alpha("blue", 0.2)
        )

        legend(
          min(env[, j]),
          0.25,
          legend = c(
            paste(rownames(candidate_list)[p]),
            paste("z-score =", format(candidate_list$zscore[p], digits = 3)),
            paste("p-value =", format(candidate_list$pvalue[p], digits = 2)),
            paste("q-value =", format(candidate_list$qvalue[p], digits = 4))
          ),
          bty = "n",
          cex = 0.7
        )
      }

      dev.off()
    }
  }
}

# -----------------------------------------------------------------------------
# Inspect significant univariate associations at K = 11
# -----------------------------------------------------------------------------

full_list <- NULL

for (j in seq_len(NCOL(env))) {
  tmp_list <- read.table(
    file.path(
      lfmm_raw_dir,
      paste0("env", j, "_", colnames(env)[j]),
      paste0("K", Ks),
      paste0("CandidatesOrdered_env", j, "_", colnames(env)[j], "_K", Ks, "_q", fdr_output, ".csv")
    ),
    header = TRUE,
    sep = ","
  )

  full_list <- c(full_list, as.vector(unlist(tmp_list$SNPid)))
}

df_summary <- as.data.frame(table(full_list))

# -----------------------------------------------------------------------------
# Multivariate LFMM analyses
# -----------------------------------------------------------------------------

dir.create(file.path(lfmm_dir, "envAll"), recursive = TRUE, showWarnings = FALSE)

for (i in seq_len(Ks)) {
  dir.create(file.path(lfmm_dir, "envAll", paste0("K", i)), recursive = TRUE, showWarnings = FALSE)
}

gif <- matrix(nrow = 1, ncol = Ks)
rownames(gif) <- "envAll"
colnames(gif) <- paste0("K", seq_len(Ks))

nb_asso_q <- matrix(nrow = 1, ncol = length(fdr_thres))
rownames(nb_asso_q) <- "value"
colnames(nb_asso_q) <- c("q0.1", "q0.01", "q0.001")

nb_asso_k <- matrix(nrow = 1, ncol = Ks)
rownames(nb_asso_k) <- "envAll"
colnames(nb_asso_k) <- paste0("K", seq_len(Ks))

for (i in seq_len(Ks)) {
  message(paste0("Env All , K = ", i))

  run_dir <- file.path(lfmm_dir, "envAll", paste0("K", i))
  setwd(run_dir)

  res <- matrix(nrow = NCOL(Y), ncol = 5)
  rownames(res) <- colnames(Y)
  colnames(res) <- c("SNPid", "fscores", "pvalue", "adjRsquared", "qvalue")

  mod_lfmm <- lfmm2(input = Y, env = env, K = i, lambda = 1e-5, effect.sizes = TRUE)
  stats_lfmm <- lfmm2.test(
    object = mod_lfmm,
    input = Y,
    env = env,
    full = TRUE,
    genomic.control = TRUE
  )

  gif[, i] <- stats_lfmm$gif
  qv_lfmm <- qvalue::qvalue(as.vector(stats_lfmm$pvalues), fdr.level = fdr_output)

  res[, "SNPid"] <- snp_info$SNPid
  res[, "fscores"] <- stats_lfmm$fscores
  res[, "pvalue"] <- stats_lfmm$pvalues
  res[, "adjRsquared"] <- stats_lfmm$adj.r.squared
  res[, "qvalue"] <- qv_lfmm$qvalues

  write.table(
    res,
    paste0("LFMM_AllResults_envAll_K", i, ".csv"),
    sep = ",",
    row.names = FALSE,
    col.names = TRUE,
    quote = FALSE
  )

  write_pvalue_diagnostics(
    stats_lfmm$pvalues,
    paste0("PvalueDistribution_envAll_K", "_K", i, ".pdf")
  )

  for (x in seq_along(fdr_thres)) {
    q <- fdr_thres[x]
    w <- which(sort(qv_lfmm$qvalues) <= q)
    candidate <- order(stats_lfmm$pvalues, decreasing = FALSE)[w]

    write_lfmm_manhattan_pdf(
      pvalues = stats_lfmm$pvalues,
      candidate = candidate,
      main = paste0("Manhattan plot | envAll | K", i, " | q=", q),
      output_file = paste0("ManhattanPlot_envAll_K", i, "_q", q, ".pdf")
    )

    if (q == fdr_output) {
      nb_asso_k[1, i] <- length(candidate)
    }

    nb_asso_q["value", x] <- length(candidate)
    candidate <- res[candidate, ]

    if (NCOL(candidate) == 1) {
      write.table(
        t(candidate),
        paste0("CandidatesOrdered_envAll_K", i, "_q", q, ".csv"),
        sep = ",",
        row.names = FALSE,
        col.names = TRUE,
        quote = FALSE
      )
    } else if (NCOL(candidate) > 1) {
      write.table(
        candidate[order(as.numeric(candidate[, "pvalue"]), decreasing = FALSE), ],
        paste0("CandidatesOrdered_envAll", i, "_q", q, ".csv"),
        sep = ",",
        row.names = FALSE,
        col.names = TRUE,
        quote = FALSE
      )
    } else {
      write.table(
        candidate,
        paste0("CandidatesOrdered_envAll_K", i, "_q", q, ".csv"),
        sep = ",",
        row.names = TRUE,
        col.names = FALSE,
        quote = FALSE
      )
    }
  }

  write.table(
    nb_asso_q,
    paste0("AssociationsNb_envAll_K", i, "_q.csv"),
    sep = ",",
    row.names = FALSE,
    col.names = TRUE,
    quote = FALSE
  )
}

gif_out <- cbind(rownames(gif), gif)
colnames(gif_out)[1] <- "env"
write.table(
  gif_out,
  file.path(multivariate_summary_dir, paste0("GenomicInflationFactor_K1-", Ks, "_mulivariate.csv")),
  sep = ",",
  row.names = FALSE,
  col.names = TRUE,
  quote = FALSE
)

nb_asso_k_out <- cbind(rownames(nb_asso_k), nb_asso_k)
colnames(nb_asso_k_out)[1] <- "env"
write.table(
  nb_asso_k_out,
  file.path(multivariate_summary_dir, paste0("NbSignifAssociations_K1-", Ks, "_q", fdr_output, "_mulivariate.csv")),
  sep = ",",
  row.names = FALSE,
  col.names = TRUE,
  quote = FALSE
)

# -----------------------------------------------------------------------------
# Selected LFMM Manhattan plots from K = 11 output
# -----------------------------------------------------------------------------

prepare_lfmm_top1_manhattan_data <- function(result, chromosome_map) {
  result <- result %>%
    separate(
      SNPid,
      into = c("CHR", "BP"),
      sep = "_(?=[^_]+$)",
      extra = "merge",
      remove = FALSE
    )

  colnames(result)[1] <- "nonsense"
  colnames(result)[2] <- "CHR"
  colnames(result)[3] <- "BP"
  colnames(result)[4] <- "zscore"
  colnames(result)[5] <- "P_not_use"
  colnames(result)[6] <- "P"

  result$CHR <- chromosome_map[result$CHR]
  result$CHR <- as.numeric(result$CHR)
  result$BP <- as.numeric(result$BP)
  result$P <- as.numeric(result$P)

  result <- result[complete.cases(result$CHR), ]
  result <- result[order(result$CHR, result$BP), ]

  top_1_percent_threshold <- quantile(result$P_not_use, 0.01)
  top_1_percent <- result[result$P_not_use <= top_1_percent_threshold, ]

  result <- result %>%
    mutate(top_1_percent = ifelse(nonsense %in% top_1_percent$nonsense, "yes", "no"))

  chr_offsets <- c(0, cumsum(as.numeric(tapply(result$BP, result$CHR, max))))
  result$cumBP <- result$BP

  for (chr in unique(result$CHR)) {
    result$cumBP[result$CHR == chr] <- result$BP[result$CHR == chr] + chr_offsets[chr]
  }

  result <- result %>%
    mutate(plot_color = ifelse(top_1_percent == "yes", "black", paste0("CHR_", CHR)))

  result
}

plot_lfmm_top1_manhattan <- function(result) {
  chr_colors <- rep(c("darkblue", "darkred"), length.out = length(unique(result$CHR)))
  names(chr_colors) <- paste0("CHR_", unique(result$CHR))

  all_colors <- c(chr_colors, "black" = "black")

  ggplot(result, aes(x = cumBP, y = -log10(P_not_use), color = plot_color)) +
    geom_point(alpha = 0.75, size = 1.5) +
    scale_color_manual(values = all_colors) +
    theme_minimal() +
    labs(x = "Chromosome", y = "-log10(p)") +
    theme(legend.position = "none") +
    ylim(0, 10) +
    scale_x_continuous(
      breaks = tapply(result$cumBP, result$CHR, mean),
      labels = unique(result$CHR)
    ) +
    theme(
      axis.text.x = element_text(angle = 0, hjust = 1, face = "bold"),
      axis.text.y = element_text(face = "bold"),
      axis.title = element_text(face = "bold"),
      panel.grid = element_blank(),
      panel.border = element_rect(color = "black", fill = NA)
    )
}

selected_manhattan_inputs <- list(
  list(
    dir = file.path(lfmm_raw_dir, "env4_mean_tmea", "K11"),
    file = "LFMM_AllResults_env4_mean_tmea_K11.csv",
    output_file = "beech_lfmm_tmean_top1_identified_121224.csv"
  ),
  list(
    dir = file.path(lfmm_raw_dir, "env1_mean_prec", "K11"),
    file = "LFMM_AllResults_env1_mean_prec_K11.csv",
    output_file = "beech_lfmm_prec_top1_identified_121224.csv"
  )
)

selected_lfmm_manhattan_plots <- lapply(selected_manhattan_inputs, function(input) {
  setwd(input$dir)

  result <- read.csv(input$file)
  result <- prepare_lfmm_top1_manhattan_data(result, beech_chromosome_map)

  plot <- plot_lfmm_top1_manhattan(result)
  print(plot)

  write.csv(result, file = input$output_file)

  plot
})

