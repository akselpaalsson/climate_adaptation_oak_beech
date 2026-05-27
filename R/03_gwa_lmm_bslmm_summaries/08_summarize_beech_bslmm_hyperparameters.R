rm(list = ls())
gc()

library(data.table)

output_dir <- "/path/to/output_dir"

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

hyperparameters <- c("pi", "n_gamma", "rho", "h", "pve", "pge")

hyp_column_names <- c(
  "h",
  "pve",
  "rho",
  "pge",
  "pi",
  "n_gamma"
)

beech_hyper_models <- list(
  list(
    trait = "tmean_sum",
    species = "beech",
    input_dir = "path/to/hypers",
    hyp_files = c(
      "beech_tmean_mean_corr_gwas_chain1_with_pcs.recoderesiduals_40m__corr_tmean_mean_sum_vs_site_and_ci_basum_bslmm_with_cov_pc9_thinned_relatedness_chain1_with_pcs.hyp.txt",
      "beech_tmean_mean_corr_gwas_chain2_with_pcs.recoderesiduals_40m__corr_tmean_mean_sum_vs_site_and_ci_basum_bslmm_with_cov_pc9_thinned_relatedness_chain2_with_pcs.hyp.txt",
      "beech_tmean_mean_corr_gwas_chain3_with_pcs.recoderesiduals_40m__corr_tmean_mean_sum_vs_site_and_ci_basum_bslmm_with_cov_pc9_thinned_relatedness_chain3_with_pcs.hyp.txt",
      "beech_tmean_mean_corr_gwas_chain4_with_pcs.recoderesiduals_40m__corr_tmean_mean_sum_vs_site_and_ci_basum_bslmm_with_cov_pc9_thinned_relatedness_chain4_with_pcs.hyp.txt"
    ),
    output_file = "beech_tmean_sum_hypers_95CI.csv"
  ),
  list(
    trait = "tmean_spr",
    species = "beech",
    input_dir = "path/to/hypers",
    hyp_files = c(
      "beech_tmean_mean_corr_gwas_chain1_with_pcs.recoderesiduals_40m__corr_tmean_mean_spr_vs_site_and_block_bslmm_with_thinned_relatedness_chain1_with_pcs.hyp.txt",
      "beech_tmean_mean_corr_gwas_chain2_with_pcs.recoderesiduals_40m__corr_tmean_mean_spr_vs_site_and_block_bslmm_with_thinned_relatedness_chain2_with_pcs.hyp.txt",
      "beech_tmean_mean_corr_gwas_chain3_with_pcs.recoderesiduals_40m__corr_tmean_mean_spr_vs_site_and_block_bslmm_with_thinned_relatedness_chain3_with_pcs.hyp.txt",
      "beech_tmean_mean_corr_gwas_chain4_with_pcs.recoderesiduals_40m__corr_tmean_mean_spr_vs_site_and_block_bslmm_with_thinned_relatedness_chain4_with_pcs.hyp.txt"
    ),
    output_file = "beech_tmean_spr_hypers_95CI.csv"
  ),
  list(
    trait = "sla",
    species = "beech",
    input_dir = "path/to/hypers",
    hyp_files = c(
      "beech_sla_gwas_chain1_with_pcs.recoderesiduals_40m__sla_vs_site_and_ci_basum_and_block_bslmm_with_cov_pc6_thinned_relatedness_chain1_with_pcs.hyp.txt",
      "beech_sla_gwas_chain2_with_pcs.recoderesiduals_40m__sla_vs_site_and_ci_basum_and_block_bslmm_with_cov_pc6_thinned_relatedness_chain2_with_pcs.hyp.txt",
      "beech_sla_gwas_chain3_with_pcs.recoderesiduals_40m__sla_vs_site_and_ci_basum_and_block_bslmm_with_cov_pc6_thinned_relatedness_chain3_with_pcs.hyp.txt",
      "beech_sla_gwas_chain4_with_pcs.recoderesiduals_40m__sla_vs_site_and_ci_basum_and_block_bslmm_with_cov_pc6_thinned_relatedness_chain4_with_pcs.hyp.txt"
    ),
    output_file = "beech_sla_hypers_95CI.csv"
  ),
  list(
    trait = "rt",
    species = "beech",
    input_dir = "path/to/hypers",
    hyp_files = c(
      "beech_rt_gwas_chain1_with_pcs.recoderesiduals_40m__rt_vs_site_and_block_bslmm_with_cov_pc5_thinned_relatedness_chain1_with_pcs.hyp.txt",
      "beech_rt_gwas_chain2_with_pcs.recoderesiduals_40m__rt_vs_site_and_block_bslmm_with_cov_pc5_thinned_relatedness_chain2_with_pcs.hyp.txt",
      "beech_rt_gwas_chain3_with_pcs.recoderesiduals_40m__rt_vs_site_and_block_bslmm_with_cov_pc5_thinned_relatedness_chain3_with_pcs.hyp.txt",
      "beech_rt_gwas_chain4_with_pcs.recoderesiduals_40m__rt_vs_site_and_block_bslmm_with_cov_pc5_thinned_relatedness_chain4_with_pcs.hyp.txt"
    ),
    output_file = "beech_rt_hypers_95CI.csv"
  ),
  list(
    trait = "rs",
    species = "beech",
    input_dir = "path/to/hypers",
    hyp_files = c(
      "beech_rs_gwas_chain1_with_pcs.recoderesiduals_40m__rs_vs_site_bslmm_with_cov_pc2_thinned_relatedness_chain1_with_pcs.hyp.txt",
      "beech_rs_gwas_chain1_with_pcs.recoderesiduals_40m__rs_vs_site_bslmm_with_cov_pc2_thinned_relatedness_chain1_with_pcs.hyp.txt",
      "beech_rs_gwas_chain1_with_pcs.recoderesiduals_40m__rs_vs_site_bslmm_with_cov_pc2_thinned_relatedness_chain1_with_pcs.hyp.txt",
      "beech_rs_gwas_chain1_with_pcs.recoderesiduals_40m__rs_vs_site_bslmm_with_cov_pc2_thinned_relatedness_chain1_with_pcs.hyp.txt"
    ),
    output_file = "beech_rs_hypers_95CI.csv"
  ),
  list(
    trait = "rc",
    species = "beech",
    input_dir = "path/to/hypers",
    hyp_files = c(
      "beech_rc_gwas_chain1_with_pcs.recoderesiduals_40m__rc_vs_site_and_ci_basum_bslmm_with_cov_pc7_thinned_relatedness_chain1_with_pcs.hyp.txt",
      "beech_rc_gwas_chain2_with_pcs.recoderesiduals_40m__rc_vs_site_and_ci_basum_bslmm_with_cov_pc7_thinned_relatedness_chain2_with_pcs.hyp.txt",
      "beech_rc_gwas_chain3_with_pcs.recoderesiduals_40m__rc_vs_site_and_ci_basum_bslmm_with_cov_pc7_thinned_relatedness_chain3_with_pcs.hyp.txt",
      "beech_rc_gwas_chain4_with_pcs.recoderesiduals_40m__rc_vs_site_and_ci_basum_bslmm_with_cov_pc7_thinned_relatedness_chain4_with_pcs.hyp.txt"
    ),
    output_file = "beech_rc_hypers_95CI.csv"
  ),
  list(
    trait = "prec_sum",
    species = "beech",
    input_dir = "path/to/hypers",
    hyp_files = c(
      "beech_prec_mean_corr_gwas_chain1_with_pcs.recoderesiduals_40m__corr_prec_mean_sum_vs_site_and_block_bslmm_with_cov_pc7_thinned_relatedness_chain1_with_pcs.hyp.txt",
      "beech_prec_mean_corr_gwas_chain2_with_pcs.recoderesiduals_40m__corr_prec_mean_sum_vs_site_and_block_bslmm_with_cov_pc7_thinned_relatedness_chain2_with_pcs.hyp.txt",
      "beech_prec_mean_corr_gwas_chain3_with_pcs.recoderesiduals_40m__corr_prec_mean_sum_vs_site_and_block_bslmm_with_cov_pc7_thinned_relatedness_chain3_with_pcs.hyp.txt",
      "beech_prec_mean_corr_gwas_chain4_with_pcs.recoderesiduals_40m__corr_prec_mean_sum_vs_site_and_block_bslmm_with_cov_pc7_thinned_relatedness_chain4_with_pcs.hyp.txt"
    ),
    output_file = "beech_prec_sum_hypers_95CI.csv"
  ),
  list(
    trait = "prec_spr",
    species = "beech",
    input_dir = "path/to/hypers",
    hyp_files = c(
      "beech_prec_mean_corr_gwas_chain1_with_pcs.recoderesiduals_40m__corr_prec_mean_spr_vs_site_and_block_bslmm_with_cov_pc1_thinned_relatedness_chain1_with_pcs.hyp.txt",
      "beech_prec_mean_corr_gwas_chain2_with_pcs.recoderesiduals_40m__corr_prec_mean_spr_vs_site_and_block_bslmm_with_cov_pc1_thinned_relatedness_chain2_with_pcs.hyp.txt",
      "beech_prec_mean_corr_gwas_chain3_with_pcs.recoderesiduals_40m__corr_prec_mean_spr_vs_site_and_block_bslmm_with_cov_pc1_thinned_relatedness_chain3_with_pcs.hyp.txt",
      "beech_prec_mean_corr_gwas_chain4_with_pcs.recoderesiduals_40m__corr_prec_mean_spr_vs_site_and_block_bslmm_with_cov_pc1_thinned_relatedness_chain4_with_pcs.hyp.txt"
    ),
    output_file = "beech_prec_spr_hypers_95CI.csv"
  ),
  list(
    trait = "height",
    species = "beech",
    input_dir = "path/to/hypers",
    hyp_files = c(
      "beech_height_gwas_chain1.recoderesiduals_40m__height_vs_site_and_ci_basum_and_block_bslmm_with_cov_pc1_thinned_relatedness_chain1.hyp.txt",
      "beech_height_gwas_chain2.recoderesiduals_40m__height_vs_site_and_ci_basum_and_block_bslmm_with_cov_pc1_thinned_relatedness_chain2.hyp.txt",
      "beech_height_gwas_chain3.recoderesiduals_40m__height_vs_site_and_ci_basum_and_block_bslmm_with_cov_pc1_thinned_relatedness_chain3.hyp.txt",
      "beech_height_gwas_chain4.recoderesiduals_40m__height_vs_site_and_ci_basum_and_block_bslmm_with_cov_pc1_thinned_relatedness_chain4.hyp.txt"
    ),
    output_file = "beech_height_hypers_95CI.csv"
  ),
  list(
    trait = "dbh",
    species = "beech",
    input_dir = "path/to/hypers",
    hyp_files = c(
      "beech_dbh_gwas_chain1.recoderesiduals_40m__dbh_vs_site_bslmm_with_cov_pc4_thinned_relatedness_chain1.hyp.txt",
      "beech_dbh_gwas_chain2.recoderesiduals_40m__dbh_vs_site_bslmm_with_cov_pc4_thinned_relatedness_chain2.hyp.txt",
      "beech_dbh_gwas_chain3.recoderesiduals_40m__dbh_vs_site_bslmm_with_cov_pc4_thinned_relatedness_chain3.hyp.txt",
      "beech_dbh_gwas_chain4.recoderesiduals_40m__dbh_vs_site_bslmm_with_cov_pc4_thinned_relatedness_chain4.hyp.txt"
    ),
    output_file = "beech_dbh_hypers_95CI.csv"
  ),
  list(
    trait = "bai",
    species = "beech",
    input_dir = "path/to/hypers",
    hyp_files = c(
      "gemma_bai_populations_no_no277_reps_missing_removed_only_petrea_minDP3_meanDP10_maxDP80_NA0.90_chain1.MAF0.05.recoderesiduals_bai_vs_site_and_block_and_competition_lmm_with_thinned_relatedness_s_40m_chain1.hyp.txt",
      "gemma_bai_populations_no_no277_reps_missing_removed_only_petrea_minDP3_meanDP10_maxDP80_NA0.90_chain2.MAF0.05.recoderesiduals_bai_vs_site_and_block_and_competition_lmm_with_thinned_relatedness_s_40m_chain2.hyp.txt",
      "gemma_bai_populations_no_no277_reps_missing_removed_only_petrea_minDP3_meanDP10_maxDP80_NA0.90_chain3.MAF0.05.recoderesiduals_bai_vs_site_and_block_and_competition_lmm_with_thinned_relatedness_s_40m_chain3.hyp.txt",
      "gemma_bai_populations_no_no277_reps_missing_removed_only_petrea_minDP3_meanDP10_maxDP80_NA0.90_chain4.MAF0.05.recoderesiduals_bai_vs_site_and_block_and_competition_lmm_with_thinned_relatedness_s_40m_chain4.hyp.txt"
    ),
    output_file = "beech_bai_hypers_95CI.csv"
  )
)

read_hyp_chain <- function(input_dir, file_name) {
  file_path <- file.path(input_dir, file_name)
  
  if (!file.exists(file_path)) {
    stop("Missing .hyp file: ", file_path, call. = FALSE)
  }
  
  first_line <- readLines(file_path, n = 1, warn = FALSE)
  
  if (length(first_line) == 0 || !nzchar(trimws(first_line))) {
    stop("File has no readable header: ", file_path, call. = FALSE)
  }
  
  header_fields <- strsplit(trimws(first_line), "\\s+")[[1]]
  
  if (!all(hyp_column_names %in% header_fields)) {
    warning(
      "Header does not exactly match expected hyperparameter names in file: ",
      file_path,
      "\nObserved header: ",
      paste(header_fields, collapse = ", "),
      "\nExpected columns: ",
      paste(hyp_column_names, collapse = ", ")
    )
  }
  
  hyp <- fread(
    file_path,
    header = FALSE,
    skip = 1,
    fill = TRUE,
    data.table = TRUE
  )
  
  empty_cols <- names(hyp)[vapply(hyp, function(x) all(is.na(x) | trimws(as.character(x)) == ""), logical(1))]
  
  if (length(empty_cols) > 0) {
    hyp[, (empty_cols) := NULL]
  }
  
  if (ncol(hyp) == length(hyp_column_names) + 1) {
    message("Dropping leading extra/index column from: ", basename(file_path))
    hyp <- hyp[, -1, with = FALSE]
  }
  
  if (ncol(hyp) != length(hyp_column_names)) {
    stop(
      "Unexpected number of data columns in file: ",
      file_path,
      "\nObserved data columns after cleanup: ",
      ncol(hyp),
      "\nExpected data columns: ",
      length(hyp_column_names),
      "\nExpected names: ",
      paste(hyp_column_names, collapse = ", "),
      call. = FALSE
    )
  }
  
  setnames(hyp, hyp_column_names)
  
  for (column_name in hyp_column_names) {
    hyp[, (column_name) := as.numeric(get(column_name))]
  }
  
  hyp
}

pool_hyperparameter <- function(chains, hyperparameter) {
  values <- unlist(
    lapply(chains, function(chain) chain[[hyperparameter]]),
    use.names = FALSE
  )
  
  values <- as.numeric(values)
  values[is.finite(values)]
}

summarise_hyperparameter <- function(values, trait, species, hyperparameter) {
  if (length(values) == 0) {
    warning(
      "No finite values found for trait = ",
      trait,
      ", species = ",
      species,
      ", hyperparameter = ",
      hyperparameter
    )
    
    return(
      data.frame(
        Min. = NA_real_,
        `X2.5%` = NA_real_,
        Median = NA_real_,
        Mean = NA_real_,
        `X97.5%` = NA_real_,
        Max. = NA_real_,
        trait = trait,
        species = species,
        hyper = hyperparameter,
        check.names = FALSE
      )
    )
  }
  
  quantiles <- quantile(
    values,
    probs = c(0.025, 0.50, 0.975),
    names = FALSE,
    na.rm = TRUE
  )
  
  data.frame(
    Min. = min(values, na.rm = TRUE),
    `X2.5%` = quantiles[1],
    Median = quantiles[2],
    Mean = mean(values, na.rm = TRUE),
    `X97.5%` = quantiles[3],
    Max. = max(values, na.rm = TRUE),
    trait = trait,
    species = species,
    hyper = hyperparameter,
    check.names = FALSE
  )
}

summarise_trait <- function(model) {
  message("Processing ", model$species, " / ", model$trait)
  
  chains <- lapply(
    model$hyp_files,
    function(file_name) {
      read_hyp_chain(
        input_dir = model$input_dir,
        file_name = file_name
      )
    }
  )
  
  summary_table <- do.call(
    rbind,
    lapply(hyperparameters, function(hyperparameter) {
      values <- pool_hyperparameter(
        chains = chains,
        hyperparameter = hyperparameter
      )
      
      summarise_hyperparameter(
        values = values,
        trait = model$trait,
        species = model$species,
        hyperparameter = hyperparameter
      )
    })
  )
  
  output_path <- file.path(output_dir, model$output_file)
  
  write.csv(
    summary_table,
    file = output_path,
    row.names = FALSE
  )
  
  message("Wrote: ", output_path)
  
  summary_table
}

beech_hyper_summaries <- lapply(beech_hyper_models, summarise_trait)

names(beech_hyper_summaries) <- vapply(
  beech_hyper_models,
  `[[`,
  character(1),
  "trait"
)

beech_hyper_summary_all <- do.call(rbind, beech_hyper_summaries)

write.csv(
  beech_hyper_summary_all,
  file = file.path(output_dir, "beech_all_traits_hypers_95CI.csv"),
  row.names = FALSE
)

beech_hyper_summary_all
