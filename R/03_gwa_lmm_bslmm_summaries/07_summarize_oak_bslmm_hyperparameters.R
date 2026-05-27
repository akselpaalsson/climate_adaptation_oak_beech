
rm(list = ls())
gc()

library(data.table)

output_dir <- "path/to/output_dir/"

hyperparameters <- c("pi", "n_gamma", "rho", "h", "pve", "pge")

hyp_column_names <- c(
  "h",
  "pve",
  "rho",
  "pge",
  "pi",
  "n_gamma",
  "nonsense"
)

oak_hyper_models <- list(
  list(
    trait = "tmean_sum",
    species = "oak",
    input_dir = "/path/to/hypers",
    hyp_files = c(
      "gemma_tmean_sum_populations_no_no277_reps_missing_removed_only_petrea_minDP3_meanDP10_maxDP80_NA0.90_chain1.MAF0.05.recoderesiduals_tmean_sum_vs_site_and_block_and_competition_lmm_with_thinned_relatedness_s_40m_chain1.hyp.txt",
      "gemma_tmean_sum_populations_no_no277_reps_missing_removed_only_petrea_minDP3_meanDP10_maxDP80_NA0.90_chain2.MAF0.05.recoderesiduals_tmean_sum_vs_site_and_block_and_competition_lmm_with_thinned_relatedness_s_40m_chain2.hyp.txt",
      "gemma_tmean_sum_populations_no_no277_reps_missing_removed_only_petrea_minDP3_meanDP10_maxDP80_NA0.90_chain3.MAF0.05.recoderesiduals_tmean_sum_vs_site_and_block_and_competition_lmm_with_thinned_relatedness_s_40m_chain3.hyp.txt",
      "gemma_tmean_sum_populations_no_no277_reps_missing_removed_only_petrea_minDP3_meanDP10_maxDP80_NA0.90_chain4.MAF0.05.recoderesiduals_tmean_sum_vs_site_and_block_and_competition_lmm_with_thinned_relatedness_s_40m_chain4.hyp.txt"
    ),
    output_file = "oak_tmean_sum_hypers_95CI.csv"
  ),
  list(
    trait = "tmean_spr",
    species = "oak",
    input_dir = "/path/to/hypers",
    hyp_files = c(
      "gemma_tmean_spr_populations_no_no277_reps_missing_removed_only_petrea_minDP3_meanDP10_maxDP80_NA0.90_chain1.MAF0.05.recoderesiduals_tmean_spr_vs_site_and_block_and_competition_lmm_with_thinned_relatedness_s_40m_chain1.hyp.txt",
      "gemma_tmean_spr_populations_no_no277_reps_missing_removed_only_petrea_minDP3_meanDP10_maxDP80_NA0.90_chain2.MAF0.05.recoderesiduals_tmean_spr_vs_site_and_block_and_competition_lmm_with_thinned_relatedness_s_40m_chain2.hyp.txt",
      "gemma_tmean_spr_populations_no_no277_reps_missing_removed_only_petrea_minDP3_meanDP10_maxDP80_NA0.90_chain3.MAF0.05.recoderesiduals_tmean_spr_vs_site_and_block_and_competition_lmm_with_thinned_relatedness_s_40m_chain3.hyp.txt",
      "gemma_tmean_spr_populations_no_no277_reps_missing_removed_only_petrea_minDP3_meanDP10_maxDP80_NA0.90_chain4.MAF0.05.recoderesiduals_tmean_spr_vs_site_and_block_and_competition_lmm_with_thinned_relatedness_s_40m_chain4.hyp.txt"
    ),
    output_file = "oak_tmean_spr_hypers_95CI.csv"
  ),
  list(
    trait = "sla",
    species = "oak",
    input_dir = "/path/to/hypers",
    hyp_files = c(
      "gemma_sla_populations_no_no277_reps_missing_removed_only_petrea_minDP3_meanDP10_maxDP80_NA0.90_chain1.MAF0.05.recoderesiduals_sla_vs_site_and_competition_lmm_with_thinned_relatedness_s_40m_chain1.hyp.txt",
      "gemma_sla_populations_no_no277_reps_missing_removed_only_petrea_minDP3_meanDP10_maxDP80_NA0.90_chain2.MAF0.05.recoderesiduals_sla_vs_site_and_competition_lmm_with_thinned_relatedness_s_40m_chain2.hyp.txt",
      "gemma_sla_populations_no_no277_reps_missing_removed_only_petrea_minDP3_meanDP10_maxDP80_NA0.90_chain3.MAF0.05.recoderesiduals_sla_vs_site_and_competition_lmm_with_thinned_relatedness_s_40m_chain3.hyp.txt",
      "gemma_sla_populations_no_no277_reps_missing_removed_only_petrea_minDP3_meanDP10_maxDP80_NA0.90_chain4.MAF0.05.recoderesiduals_sla_vs_site_and_competition_lmm_with_thinned_relatedness_s_40m_chain4.hyp.txt"
    ),
    output_file = "oak_sla_hypers_95CI.csv"
  ),
  list(
    trait = "rt",
    species = "oak",
    input_dir = "/path/to/hypers",
    hyp_files = c(
      "gemma_rt_populations_no_no277_reps_missing_removed_only_petrea_minDP3_meanDP10_maxDP80_NA0.90_chain1.MAF0.05.recoderesiduals_rt_vs_site_and_competition_lmm_with_thinned_relatedness_s_40m_chain1.hyp.txt",
      "gemma_rt_populations_no_no277_reps_missing_removed_only_petrea_minDP3_meanDP10_maxDP80_NA0.90_chain2.MAF0.05.recoderesiduals_rt_vs_site_and_competition_lmm_with_thinned_relatedness_s_40m_chain2.hyp.txt",
      "gemma_rt_populations_no_no277_reps_missing_removed_only_petrea_minDP3_meanDP10_maxDP80_NA0.90_chain3.MAF0.05.recoderesiduals_rt_vs_site_and_competition_lmm_with_thinned_relatedness_s_40m_chain3.hyp.txt",
      "gemma_rt_populations_no_no277_reps_missing_removed_only_petrea_minDP3_meanDP10_maxDP80_NA0.90_chain4.MAF0.05.recoderesiduals_rt_vs_site_and_competition_lmm_with_thinned_relatedness_s_40m_chain4.hyp.txt"
    ),
    output_file = "oak_rt_hypers_95CI.csv"
  ),
  list(
    trait = "rs",
    species = "oak",
    input_dir = "/path/to/hypers",
    hyp_files = c(
      "gemma_rs_populations_no_no277_reps_missing_removed_only_petrea_minDP3_meanDP10_maxDP80_NA0.90_chain1.MAF0.05.recoderesiduals_rs_vs_site_and_block_and_competition_lmm_with_thinned_relatedness_s_40m_chain1.hyp.txt",
      "gemma_rs_populations_no_no277_reps_missing_removed_only_petrea_minDP3_meanDP10_maxDP80_NA0.90_chain2.MAF0.05.recoderesiduals_rs_vs_site_and_block_and_competition_lmm_with_thinned_relatedness_s_40m_chain2.hyp.txt",
      "gemma_rs_populations_no_no277_reps_missing_removed_only_petrea_minDP3_meanDP10_maxDP80_NA0.90_chain3.MAF0.05.recoderesiduals_rs_vs_site_and_block_and_competition_lmm_with_thinned_relatedness_s_40m_chain3.hyp.txt",
      "gemma_rs_populations_no_no277_reps_missing_removed_only_petrea_minDP3_meanDP10_maxDP80_NA0.90_chain4.MAF0.05.recoderesiduals_rs_vs_site_and_block_and_competition_lmm_with_thinned_relatedness_s_40m_chain4.hyp.txt"
    ),
    output_file = "oak_rs_hypers_95CI.csv"
  ),
  list(
    trait = "rc",
    species = "oak",
    input_dir = "/path/to/hypers",
    hyp_files = c(
      "gemma_rc_populations_no_no277_reps_missing_removed_only_petrea_minDP3_meanDP10_maxDP80_NA0.90_chain1.MAF0.05.recoderesiduals_rc_vs_site_and_block_lmm_with_thinned_relatedness_s_40m_chain1.hyp.txt",
      "gemma_rc_populations_no_no277_reps_missing_removed_only_petrea_minDP3_meanDP10_maxDP80_NA0.90_chain2.MAF0.05.recoderesiduals_rc_vs_site_and_block_lmm_with_thinned_relatedness_s_40m_chain2.hyp.txt",
      "gemma_rc_populations_no_no277_reps_missing_removed_only_petrea_minDP3_meanDP10_maxDP80_NA0.90_chain3.MAF0.05.recoderesiduals_rc_vs_site_and_block_lmm_with_thinned_relatedness_s_40m_chain3.hyp.txt",
      "gemma_rc_populations_no_no277_reps_missing_removed_only_petrea_minDP3_meanDP10_maxDP80_NA0.90_chain5.MAF0.05.recoderesiduals_rc_vs_site_and_block_lmm_with_thinned_relatedness_s_40m_chain5.hyp.txt"
    ),
    output_file = "oak_rc_hypers_95CI.csv"
  ),
  list(
    trait = "prec_sum",
    species = "oak",
    input_dir = "/path/to/hypers",
    hyp_files = c(
      "gemma_prec_sum_populations_no_no277_reps_missing_removed_only_petrea_minDP3_meanDP10_maxDP80_NA0.90_chain1.MAF0.05.recoderesiduals_prec_sum_vs_site_lmm_with_thinned_relatedness_s_40m_chain1.hyp.txt",
      "gemma_prec_sum_populations_no_no277_reps_missing_removed_only_petrea_minDP3_meanDP10_maxDP80_NA0.90_chain2.MAF0.05.recoderesiduals_prec_sum_vs_site_lmm_with_thinned_relatedness_s_40m_chain2.hyp.txt",
      "gemma_prec_sum_populations_no_no277_reps_missing_removed_only_petrea_minDP3_meanDP10_maxDP80_NA0.90_chain3.MAF0.05.recoderesiduals_prec_sum_vs_site_lmm_with_thinned_relatedness_s_40m_chain3.hyp.txt",
      "gemma_prec_sum_populations_no_no277_reps_missing_removed_only_petrea_minDP3_meanDP10_maxDP80_NA0.90_chain4.MAF0.05.recoderesiduals_prec_sum_vs_site_lmm_with_thinned_relatedness_s_40m_chain4.hyp.txt"
    ),
    output_file = "oak_prec_sum_hypers_95CI.csv"
  ),
  list(
    trait = "prec_spr",
    species = "oak",
    input_dir = "/path/to/hypers",
    hyp_files = c(
      "gemma_prec_spr_populations_no_no277_reps_missing_removed_only_petrea_minDP3_meanDP10_maxDP80_NA0.90_chain1.MAF0.05.recoderesiduals_prec_spr_vs_site_lmm_with_thinned_relatedness_s_40m_chain1.hyp.txt",
      "gemma_prec_spr_populations_no_no277_reps_missing_removed_only_petrea_minDP3_meanDP10_maxDP80_NA0.90_chain2.MAF0.05.recoderesiduals_prec_spr_vs_site_lmm_with_thinned_relatedness_s_40m_chain2.hyp.txt",
      "gemma_prec_spr_populations_no_no277_reps_missing_removed_only_petrea_minDP3_meanDP10_maxDP80_NA0.90_chain3.MAF0.05.recoderesiduals_prec_spr_vs_site_lmm_with_thinned_relatedness_s_40m_chain3.hyp.txt",
      "gemma_prec_spr_populations_no_no277_reps_missing_removed_only_petrea_minDP3_meanDP10_maxDP80_NA0.90_chain4.MAF0.05.recoderesiduals_prec_spr_vs_site_lmm_with_thinned_relatedness_s_40m_chain4.hyp.txt"
    ),
    output_file = "oak_prec_spr_hypers_95CI.csv"
  ),
  list(
    trait = "height",
    species = "oak",
    input_dir = "/path/to/hypers",
    hyp_files = c(
      "gemma_height_populations_no_no277_reps_missing_removed_only_petrea_minDP3_meanDP10_maxDP80_NA0.90_chain1.MAF0.05.recoderesiduals_height_vs_site_lmm_with_thinned_relatedness_s_40m_chain1.hyp.txt",
      "gemma_height_populations_no_no277_reps_missing_removed_only_petrea_minDP3_meanDP10_maxDP80_NA0.90_chain2.MAF0.05.recoderesiduals_height_vs_site_lmm_with_thinned_relatedness_s_40m_chain2.hyp.txt",
      "gemma_height_populations_no_no277_reps_missing_removed_only_petrea_minDP3_meanDP10_maxDP80_NA0.90_chain3.MAF0.05.recoderesiduals_height_vs_site_lmm_with_thinned_relatedness_s_40m_chain3.hyp.txt",
      "gemma_height_populations_no_no277_reps_missing_removed_only_petrea_minDP3_meanDP10_maxDP80_NA0.90_chain4.MAF0.05.recoderesiduals_height_vs_site_lmm_with_thinned_relatedness_s_40m_chain4.hyp.txt"
    ),
    output_file = "oak_height_hypers_95CI.csv"
  ),
  list(
    trait = "dbh",
    species = "oak",
    input_dir = "/path/to/hypers",
    hyp_files = c(
      "gemma_dbh_populations_no_no277_reps_missing_removed_only_petrea_minDP3_meanDP10_maxDP80_NA0.90_chain1.MAF0.05.recoderesiduals_dbh_vs_site_bslmm_with_thinned_relatedness_s_40m_chain1.hyp.txt",
      "gemma_dbh_populations_no_no277_reps_missing_removed_only_petrea_minDP3_meanDP10_maxDP80_NA0.90_chain2.MAF0.05.recoderesiduals_dbh_vs_site_bslmm_with_thinned_relatedness_s_40m_chain2.hyp.txt",
      "gemma_dbh_populations_no_no277_reps_missing_removed_only_petrea_minDP3_meanDP10_maxDP80_NA0.90_chain3.MAF0.05.recoderesiduals_dbh_vs_site_bslmm_with_thinned_relatedness_s_40m_chain3.hyp.txt",
      "gemma_dbh_populations_no_no277_reps_missing_removed_only_petrea_minDP3_meanDP10_maxDP80_NA0.90_chain4.MAF0.05.recoderesiduals_dbh_vs_site_bslmm_with_thinned_relatedness_s_40m_chain4.hyp.txt"
    ),
    output_file = "oak_dbh_hypers_95CI.csv"
  ),
  list(
    trait = "bai",
    species = "oak",
    input_dir = "/path/to/hypers",
    hyp_files = c(
      "gemma_bai_populations_no_no277_reps_missing_removed_only_petrea_minDP3_meanDP10_maxDP80_NA0.90_chain1.MAF0.05.recoderesiduals_bai_vs_site_and_block_and_competition_lmm_with_thinned_relatedness_s_40m_chain1.hyp.txt",
      "gemma_bai_populations_no_no277_reps_missing_removed_only_petrea_minDP3_meanDP10_maxDP80_NA0.90_chain2.MAF0.05.recoderesiduals_bai_vs_site_and_block_and_competition_lmm_with_thinned_relatedness_s_40m_chain2.hyp.txt",
      "gemma_bai_populations_no_no277_reps_missing_removed_only_petrea_minDP3_meanDP10_maxDP80_NA0.90_chain3.MAF0.05.recoderesiduals_bai_vs_site_and_block_and_competition_lmm_with_thinned_relatedness_s_40m_chain3.hyp.txt",
      "gemma_bai_populations_no_no277_reps_missing_removed_only_petrea_minDP3_meanDP10_maxDP80_NA0.90_chain4.MAF0.05.recoderesiduals_bai_vs_site_and_block_and_competition_lmm_with_thinned_relatedness_s_40m_chain4.hyp.txt"
    ),
    output_file = "oak_bai_hypers_95CI.csv"
  )
)

read_hyp_chain <- function(input_dir, file_name) {
  hyp <- fread(file.path(input_dir, file_name))
  colnames(hyp)[seq_along(hyp_column_names)] <- hyp_column_names
  as.data.frame(hyp)
}

pool_hyperparameter <- function(chains, hyperparameter) {
  values <- unlist(
    lapply(chains, function(chain) chain[[hyperparameter]]),
    use.names = FALSE
  )

  values[is.finite(values)]
}

summarise_hyperparameter <- function(values, trait, species, hyperparameter) {
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
  chains <- lapply(
    model$hyp_files,
    read_hyp_chain,
    input_dir = model$input_dir
  )

  summary_table <- do.call(
    rbind,
    lapply(hyperparameters, function(hyperparameter) {
      values <- pool_hyperparameter(chains, hyperparameter)

      summarise_hyperparameter(
        values = values,
        trait = model$trait,
        species = model$species,
        hyperparameter = hyperparameter
      )
    })
  )

  write.csv(
    summary_table,
    file = file.path(output_dir, model$output_file),
    row.names = FALSE
  )

  summary_table
}

oak_hyper_summaries <- lapply(oak_hyper_models, summarise_trait)
names(oak_hyper_summaries) <- vapply(oak_hyper_models, `[[`, character(1), "trait")
