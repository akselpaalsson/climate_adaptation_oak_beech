
rm(list = ls())
gc()

library(dplyr)
library(lme4)
library(lmerTest)
library(emmeans)
library(stringr)
library(tidyr)

# -----------------------------------------------------------------------------
# Paths
# -----------------------------------------------------------------------------
bslmm_dir <- "/path/to/bslmm_dir"
phenotype_dir <- "/path/to/phenotype_dir"
univar_dir <- "/path/to/univar_dir"
climate_distance_dir <- "/path/to/climate_data"
output_dir <- "/path/to/output_dir"
# -----------------------------------------------------------------------------
# Shared helpers
# -----------------------------------------------------------------------------

recode_oak_provenance <- function(provenance) {
  provenance <- gsub("FR", "FRA_", provenance)
  provenance <- gsub("DE", "GER_", provenance)
  provenance <- gsub("PL", "POL_", provenance)
  provenance <- gsub("UK", "GBR_", provenance)
  provenance <- gsub("DK", "DEN_", provenance)
  
  provenance
}

calculate_allele_freq <- function(genotypes) {
  genotypes <- genotypes[!is.na(genotypes)]
  
  total_alleles <- length(genotypes) * 2
  count_AA <- sum(genotypes == "AA") * 2
  count_AB <- sum(genotypes == "AB")
  count_BB <- sum(genotypes == "BB") * 2
  
  total_A <- count_AA + count_AB
  total_B <- count_BB + count_AB
  
  c(
    freq_A = total_A / total_alleles,
    freq_B = total_B / total_alleles
  )
}

first_existing_col <- function(data, candidates, required = TRUE, label = "column") {
  hit <- candidates[candidates %in% names(data)]
  
  if (length(hit) == 0) {
    if (required) {
      stop(
        "Could not find ", label, ". Tried: ",
        paste(candidates, collapse = ", "),
        "\nAvailable columns are: ", paste(names(data), collapse = ", "),
        call. = FALSE
      )
    }
    return(NULL)
  }
  
  hit[1]
}

add_analysis_columns <- function(data) {
  prov_col <- first_existing_col(
    data,
    c("Prov3", "Prov3.x", "Prov3.y", "Prov", "Prov.x", "Prov.y"),
    label = "provenance column"
  )
  
  site_col <- first_existing_col(
    data,
    c("Site.x", "Site.x.x", "Site", "Site.y", "Site.y.y"),
    label = "site column"
  )
  
  block_col <- first_existing_col(
    data,
    c("Block.x", "Block.x.x", "Block", "Block.y", "Block.y.y"),
    label = "block column"
  )
  
  ci_col <- first_existing_col(
    data,
    c(
      "CI_BAsum", "CI_BAsum.x", "CI_BAsum.y", "CI_BAsum.x.x", "CI_BAsum.y.y",
      "CI_Basum", "CI_Basum.x", "CI_Basum.y", "CI_Basum.x.x", "CI_Basum.y.y"
    ),
    label = "competition / CI_BAsum column"
  )
  
  data %>%
    mutate(
      Prov_analysis = .data[[prov_col]],
      Site_analysis = .data[[site_col]],
      Block_analysis = .data[[block_col]],
      CI_BAsum_analysis = .data[[ci_col]],
      merger = paste(Prov_analysis, Site_analysis, sep = "_")
    )
}

safe_mean <- function(x) {
  if (all(is.na(x))) return(NA_real_)
  mean(x, na.rm = TRUE)
}

get_interaction_row <- function(coefs, snp, distance_variable) {
  interaction_rows <- grep(":", coefs$type)
  interaction_rows <- interaction_rows[
    grepl(snp, coefs$type[interaction_rows], fixed = TRUE) &
      grepl(distance_variable, coefs$type[interaction_rows], fixed = TRUE)
  ]
  
  if (length(interaction_rows) == 0) {
    return(data.frame(type = NA_character_, Estimate = NA_real_, Std_Error = NA_real_, P_Value = NA_real_))
  }
  
  coefs[interaction_rows[1], , drop = FALSE]
}

make_chr_bp_snp_id <- function(snp_table) {
  snp_table$snp <- "snp"
  snp_table$snp_id <- paste(snp_table$snp, snp_table$chr_bp, sep = "_")
  
  as.character(snp_table$snp_id)
}

load_oak_bslmm_snp_table <- function() {
  snps <- read.csv(file.path(bslmm_dir, "oak_all_bslmm_snps_table_151024.csv"), sep = ",")
  
  snps <- snps %>%
    mutate(
      trait_real = ifelse(
        str_detect(trait, "oak_tmean|oak_prec"),
        str_extract(trait, "(?<=oak_)(\\w+_\\w+)"),
        str_extract(trait, "(?<=oak_)(\\w+)")
      )
    ) %>%
    mutate(trait_real = str_replace(trait_real, "_with.*", ""))
  
  snps
}

load_sample_info <- function() {
  sample_info <- read.table(
    file.path(phenotype_dir, "common_ring_sample_info_no_reps.csv"),
    sep = ",",
    header = TRUE
  )
  
  sample_vcf <- read.table(
    file.path(phenotype_dir, "samples_populations_no_no277_reps_missing_removed_only_petrea_minDP3_meanDP10_maxDP80_NA0.90.txt"),
    sep = "\t",
    header = FALSE
  )
  colnames(sample_vcf)[1] <- "Sample"
  
  sample_vcf$RowNumber <- as.numeric(row.names(sample_vcf))
  
  sample_info_reduced <- merge(
    sample_info,
    sample_vcf,
    by.x = "SampleID",
    by.y = "Sample",
    sort = FALSE
  )
  
  sample_info_reduced <- sample_info_reduced[order(sample_info_reduced$RowNumber), ]
  colnames(sample_info_reduced)[2] <- "Treeid"
  
  sample_info_reduced
}

load_univar_genotypes <- function() {
  read.csv(file.path(bslmm_dir, "extracted_oak_bslmm_genos_151024.csv"), sep = ",")
}

prepare_oak_dbh_sample_data <- function() {
  sample_info <- load_sample_info()
  
  oak <- read.csv(file.path(phenotype_dir, "dryad_oak_dbh.csv"))
  oak$Prov <- recode_oak_provenance(oak$Prov)
  
  oak$merging_col <- paste(oak$Site, oak$Block, oak$Prov, oak$Tree, sep = "_")
  sample_info$merging_col <- paste(
    sample_info$Site,
    sample_info$Block,
    sample_info$Prov3,
    sample_info$Tree,
    sep = "_"
  )
  
  merge(sample_info, oak, by = "merging_col")
}

prepare_corr_trait_data <- function(climate_variable, genotypes) {
  sample_info <- load_sample_info()
  
  oak_corr <- read.csv(file.path(phenotype_dir, "dryad_oak_corr.csv"))
  corr_trait <- subset(oak_corr, var0me == climate_variable)
  
  corr_trait <- corr_trait %>%
    mutate(tree_nr = as.numeric(sub(".*([0-9]+)$", "\\1", tree.id)))
  
  corr_trait$Prov <- recode_oak_provenance(corr_trait$Prov)
  
  corr_trait$merging_col <- paste(
    corr_trait$Site,
    corr_trait$Block,
    corr_trait$Prov,
    corr_trait$tree_nr,
    sep = "_"
  )
  
  sample_info$merging_col <- paste(
    sample_info$Site,
    sample_info$Block,
    sample_info$Prov3,
    sample_info$Tree,
    sep = "_"
  )
  
  merged <- merge(sample_info, corr_trait, by = "merging_col")
  merged <- merged[!duplicated(merged), ]
  
  merge(merged, genotypes, by = "SampleID")
}

prepare_growth_trait_data <- function(genotypes) {
  oak_sample_data <- prepare_oak_dbh_sample_data()
  merge(oak_sample_data, genotypes, by = "SampleID")
}

prepare_bai_trait_data <- function(genotypes) {
  sample_info <- load_sample_info()
  
  oak_bai <- read.csv(file.path(phenotype_dir, "dryad_oak_bai.csv"))
  oak_bai$Treeid2 <- substr(oak_bai$Treeid, nchar(oak_bai$Treeid), nchar(oak_bai$Treeid))
  oak_bai$Prov <- recode_oak_provenance(oak_bai$Prov)
  
  oak_bai$merging_col <- paste(
    oak_bai$Site,
    oak_bai$Block,
    oak_bai$Prov,
    oak_bai$Treeid2,
    sep = "_"
  )
  
  sample_info$merging_col <- paste(
    sample_info$Site,
    sample_info$Block,
    sample_info$Prov3,
    sample_info$Tree,
    sep = "_"
  )
  
  oak_bai_sample_info <- merge(sample_info, oak_bai, by = "merging_col")
  
  fit_lm_and_extract_estimate <- function(data) {
    lm_model <- lm(BAI ~ year, data = data)
    data.frame(
      Estimate = coef(lm_model)[2],
      SampleID = unique(data$SampleID)
    )
  }
  
  bai_slopes_as_trait <- oak_bai_sample_info %>%
    group_by(SampleID) %>%
    do(fit_lm_and_extract_estimate(.))
  
  bai_data <- merge(bai_slopes_as_trait, genotypes, by = "SampleID")
  oak_sample_data <- prepare_oak_dbh_sample_data()
  
  merge(oak_sample_data, bai_data, by = "SampleID")
}

prepare_resilience_trait_data <- function(resilience_variable, genotypes) {
  sample_info <- load_sample_info()
  sample_info$merging_col <- paste(
    sample_info$Site,
    sample_info$Block,
    sample_info$Prov3,
    sample_info$Tree,
    sep = "_"
  )
  
  oak_rs <- read.csv(file.path(phenotype_dir, "dryad_oak_rs-2.csv"))
  oak_rs$Prov <- recode_oak_provenance(oak_rs$Prov)
  oak_rs$merging_col <- paste(oak_rs$Site, oak_rs$Block, oak_rs$Prov, oak_rs$tree, sep = "_")
  
  oak_rs_trait <- subset(oak_rs, var == resilience_variable)
  resilience_data <- merge(sample_info, oak_rs_trait, by = "merging_col")
  resilience_data <- merge(resilience_data, genotypes, by = "SampleID")
  
  oak_sample_data <- prepare_oak_dbh_sample_data()
  
  merge(oak_sample_data, resilience_data, by = "SampleID")
}

load_distance_table <- function(file_name) {
  distance_table <- read.csv(file.path(climate_distance_dir, file_name))
  
  colnames(distance_table)[1] <- "Prov"
  colnames(distance_table)[2] <- "FR"
  colnames(distance_table)[3] <- "UK"
  colnames(distance_table)[4] <- "PL"
  colnames(distance_table)[5] <- "DK"
  
  distance_table <- distance_table %>%
    dplyr::select(-ID) %>%
    mutate(
      Prov = case_when(
        Prov == "DE263" ~ "GER_263",
        Prov == "PL033" ~ "POL_033",
        Prov == "DE031" ~ "GER_031",
        Prov == "DE087" ~ "GER_087",
        Prov == "DK260" ~ "DEN_260",
        Prov == "FR004" ~ "FRA_004",
        Prov == "FR050" ~ "FRA_050",
        Prov == "FR077" ~ "FRA_077",
        Prov == "UK026" ~ "GBR_026",
        TRUE ~ Prov
      )
    ) %>%
    pivot_longer(cols = -Prov, names_to = "country", values_to = "dist") %>%
    mutate(merger = paste(Prov, country, sep = "_")) %>%
    dplyr::select(Prov, dist, merger)
  
  distance_table
}

add_climate_distances <- function(selected_df) {
  selected_df <- add_analysis_columns(selected_df)
  
  dist_prec <- load_distance_table("distance_df_oak_std_prec_annual_051224.csv")
  dist_tmin <- load_distance_table("distance_df_oak_std_tmin_annual_051224.csv")
  dist_tmax <- load_distance_table("distance_df_oak_std_tmax_annual_051224.csv")
  dist_tmean <- load_distance_table("distance_df_oak_std_tmean_annual_051224.csv")
  
  selected_df_with_clima <- selected_df %>%
    left_join(dist_prec %>% dplyr::select(merger, dist_prec = dist), by = "merger") %>%
    left_join(dist_tmin %>% dplyr::select(merger, dist_tmin = dist), by = "merger") %>%
    left_join(dist_tmax %>% dplyr::select(merger, dist_tmax = dist), by = "merger") %>%
    left_join(dist_tmean %>% dplyr::select(merger, dist_tmean = dist), by = "merger")
  
  missing_distance <- sum(is.na(selected_df_with_clima$dist_prec) | is.na(selected_df_with_clima$dist_tmean))
  if (missing_distance > 0) {
    warning(
      missing_distance,
      " rows did not match climate-distance tables. Check Prov/Site keys in 'merger'.",
      call. = FALSE
    )
  }
  
  list(
    data = selected_df_with_clima,
    dist_tmean = dist_tmean
  )
}
get_trait_data <- function(config, genotypes) {
  if (config$data_type == "corr") {
    prepare_corr_trait_data(config$climate_variable, genotypes)
  } else if (config$data_type == "growth") {
    prepare_growth_trait_data(genotypes)
  } else if (config$data_type == "bai") {
    prepare_bai_trait_data(genotypes)
  } else if (config$data_type == "resilience") {
    prepare_resilience_trait_data(config$resilience_variable, genotypes)
  }
}

extract_emtrend_results <- function(snp, final_results_tmean, coefs_m1, analysis) {
  interaction_row <- get_interaction_row(coefs_m1, snp, analysis$interaction_distance)
  
  data.frame(
    snp_id = snp,
    dist.trend_AA_tmean = final_results_tmean[1, 3],
    dist.trend_AB_tmean = final_results_tmean[2, 3],
    dist.trend_SE_AA_tmean = final_results_tmean[1, 4],
    dist.trend_SE_AB_tmean = final_results_tmean[2, 4],
    dist.trend_p.value_AA_tmean = final_results_tmean[1, 7],
    dist.trend_p.value_AB_tmean = final_results_tmean[2, 7],
    dist.trend_lower.CL_AA_tmean = final_results_tmean[1, 9],
    dist.trend_lower.CL_AB_tmean = final_results_tmean[2, 9],
    dist.trend_upper.CL_AA_tmean = final_results_tmean[1, 10],
    dist.trend_upper.CL_AB_tmean = final_results_tmean[2, 10],
    contrast_tmean_id = final_results_tmean[3, 2],
    contrast_tmean_estimate = final_results_tmean[3, 8],
    contrast_tmean_estimate_lower_CL = final_results_tmean[3, 9],
    contrast_tmean_estimate_upper_CL = final_results_tmean[3, 10],
    contrast_tmean_p.value = final_results_tmean[3, 7],
    contrast_tmean_SE = final_results_tmean[3, 4],
    predictor = analysis$interaction_distance,
    coefs_m1_interaction_id_snp_dist_tmean = interaction_row$type[1],
    coefs_m1_interaction_snp_dist_tmean_estimate = interaction_row$Estimate[1],
    coefs_m1_interaction_snp_dist_tmean_SE = interaction_row$Std_Error[1],
    coefs_m1_interaction_snp_dist_tmean_p.value = interaction_row$P_Value[1]
  )
}
calculate_af_difference <- function(data, dist_tmean, snp) {
  selected_df_for_af <- add_analysis_columns(data)
  
  merged_for_af <- left_join(
    selected_df_for_af,
    dist_tmean %>% dplyr::select(merger, dist_tmean = dist),
    by = "merger"
  )
  
  filtered_df <- merged_for_af %>%
    filter(Prov_analysis %in% c("POL_033", "FRA_004"))
  
  if (nrow(filtered_df) == 0) return(NA_real_)
  
  allele_frequencies <- filtered_df %>%
    group_by(Prov_analysis) %>%
    summarise(
      freq_A = calculate_allele_freq(!!sym(snp))[1],
      freq_B = calculate_allele_freq(!!sym(snp))[2],
      .groups = "drop"
    )
  
  if (!all(c("POL_033", "FRA_004") %in% allele_frequencies$Prov_analysis)) {
    return(NA_real_)
  }
  
  allele_frequencies$freq_B[allele_frequencies$Prov_analysis == "POL_033"] -
    allele_frequencies$freq_B[allele_frequencies$Prov_analysis == "FRA_004"]
}
analyse_snp_trait <- function(data, snp, config, analysis) {
  if (!snp %in% names(data)) {
    warning("Skipping ", snp, ": SNP column not found in trait data.", call. = FALSE)
    return(NULL)
  }
  
  data <- add_analysis_columns(data)
  data[[snp]] <- as.factor(data[[snp]])
  
  data_no_bb <- data %>%
    filter(!!sym(snp) != "BB") %>%
    filter(!is.na(!!sym(snp))) %>%
    filter(!is.na(.data[[config$response]])) %>%
    filter(!is.na(CI_BAsum_analysis))
  
  if (!is.null(config$extra_non_na)) {
    data_no_bb <- data_no_bb %>%
      filter(!is.na(.data[[config$extra_non_na]]))
  }
  
  if (nrow(data_no_bb) == 0 || length(unique(data_no_bb[[snp]])) < 2) {
    warning("Skipping ", snp, ": not enough non-BB genotype classes after filtering.", call. = FALSE)
    return(NULL)
  }
  
  mean_AA <- safe_mean(data_no_bb[[config$response]][data_no_bb[[snp]] == "AA"])
  mean_AB <- safe_mean(data_no_bb[[config$response]][data_no_bb[[snp]] == "AB"])
  
  if (!is.na(mean_AA) && !is.na(mean_AB) && mean_AB > mean_AA) {
    data_no_bb[[snp]] <- ifelse(data_no_bb[[snp]] == "AA", "AB", "AA")
    data_no_bb[[snp]] <- factor(data_no_bb[[snp]])
  }
  
  response_sd <- sd(data_no_bb[[config$response]], na.rm = TRUE)
  if (is.na(response_sd) || response_sd == 0) {
    warning("Skipping ", snp, ": response has zero/NA standard deviation.", call. = FALSE)
    return(NULL)
  }
  
  data_no_bb <- data_no_bb %>%
    mutate(rt_standardized = (.data[[config$response]] - mean(.data[[config$response]], na.rm = TRUE)) / response_sd)
  
  climate_data <- add_climate_distances(data_no_bb)
  selected_df_with_clima <- climate_data$data %>%
    filter(!is.na(.data[[analysis$interaction_distance]])) %>%
    filter(!is.na(.data[[analysis$covariate_distance]]))
  
  selected_df_with_clima[[snp]] <- as.factor(selected_df_with_clima[[snp]])
  
  if (nrow(selected_df_with_clima) == 0 || length(unique(selected_df_with_clima[[snp]])) < 2) {
    warning("Skipping ", snp, ": not enough rows/genotype classes after climate-distance merge.", call. = FALSE)
    return(NULL)
  }
  
  model_formula_m1 <- as.formula(
    paste(
      "rt_standardized ~",
      paste0(
        snp,
        " * ",
        analysis$interaction_distance,
        " + ",
        analysis$covariate_distance,
        " + CI_BAsum_analysis + (1 | Site_analysis/Block_analysis)"
      )
    )
  )
  
  m1 <- tryCatch(
    lmer(model_formula_m1, data = selected_df_with_clima),
    error = function(e) {
      warning("Skipping ", snp, ": lmer failed: ", conditionMessage(e), call. = FALSE)
      NULL
    }
  )
  if (is.null(m1)) return(NULL)
  
  summary_m1 <- summary(m1)
  m1_coefficients_table <- summary_m1$coefficients
  coefs_m1 <- data.frame(
    type = rownames(m1_coefficients_table),
    Estimate = m1_coefficients_table[, "Estimate"],
    Std_Error = m1_coefficients_table[, "Std. Error"],
    P_Value = m1_coefficients_table[, "Pr(>|t|)"],
    row.names = NULL
  )
  
  emtrends_formula <- as.formula(paste("pairwise ~", snp))
  em_trends_tmean <- emtrends(
    m1,
    specs = emtrends_formula,
    var = analysis$interaction_distance,
    type = "response"
  )
  
  em_cont1_tmean <- test(em_trends_tmean)
  conf_intervals_tmean <- confint(em_trends_tmean)
  
  final_results_tmean <- data.frame(summary(em_cont1_tmean), conf_intervals_tmean)
  final_results_tmean <- final_results_tmean %>%
    dplyr::select(-contains(".1"))
  
  transformed_final_df <- extract_emtrend_results(
    snp = snp,
    final_results_tmean = final_results_tmean,
    coefs_m1 = coefs_m1,
    analysis = analysis
  )
  
  transformed_final_df$af_difference <- calculate_af_difference(
    data = data,
    dist_tmean = climate_data$dist_tmean,
    snp = snp
  )
  
  transformed_final_df
}
run_trait_analysis <- function(config, analysis, genotypes, oak_bslmm_all_snps) {
  snp_table <- subset(oak_bslmm_all_snps, trait_real == config$trait_filter)
  snp_list <- make_chr_bp_snp_id(snp_table)
  
  data <- get_trait_data(config, genotypes)
  
  final_results_all_snps <- data.frame()
  
  for (snp in snp_list) {
    transformed_final_df <- analyse_snp_trait(
      data = data,
      snp = snp,
      config = config,
      analysis = analysis
    )
    
    if (!is.null(transformed_final_df)) {
      final_results_all_snps <- rbind(final_results_all_snps, transformed_final_df)
    }
  }
  
  final_results_all_snps$trait <- config$output_trait
  
  write.csv(
    final_results_all_snps,
    file = file.path(output_dir, paste0(analysis$output_prefix, config$output_suffix, "_051224.csv"))
  )
  
  final_results_all_snps
}

# -----------------------------------------------------------------------------
# Analysis configuration
# -----------------------------------------------------------------------------

trait_configs <- list(
  list(
    trait_filter = "prec_mean_sum",
    output_trait = "sumpre",
    output_suffix = "sumpre",
    data_type = "corr",
    climate_variable = "prec.mean",
    response = "SUM",
    extra_non_na = "SPR",
    competition_covariate = "CI_BAsum_analysis",
    random_effect = "(1 | Site_analysis/Block_analysis)"
  ),
  list(
    trait_filter = "prec_mean_spr",
    output_trait = "sprpre",
    output_suffix = "sprpre",
    data_type = "corr",
    climate_variable = "prec.mean",
    response = "SPR",
    extra_non_na = "SPR",
    competition_covariate = "CI_BAsum_analysis",
    random_effect = "(1 | Site_analysis/Block_analysis)"
  ),
  list(
    trait_filter = "tmean_sum",
    output_trait = "sumtmea",
    output_suffix = "sumtmea",
    data_type = "corr",
    climate_variable = "tmea.mean",
    response = "SUM",
    extra_non_na = "SPR",
    competition_covariate = "CI_BAsum_analysis",
    random_effect = "(1 | Site_analysis/Block_analysis)"
  ),
  list(
    trait_filter = "tmean_spr",
    output_trait = "sprtmea",
    output_suffix = "sprtmea",
    data_type = "corr",
    climate_variable = "tmea.mean",
    response = "SPR",
    extra_non_na = "SPR",
    competition_covariate = "CI_BAsum_analysis",
    random_effect = "(1 | Site_analysis/Block_analysis)"
  ),
  list(
    trait_filter = "rt",
    output_trait = "rt",
    output_suffix = "rt",
    data_type = "resilience",
    resilience_variable = "rt",
    response = "value",
    extra_non_na = NULL,
    competition_covariate = "CI_BAsum_analysis",
    random_effect = "(1 | Site_analysis/Block_analysis)"
  ),
  list(
    trait_filter = "rs",
    output_trait = "rs",
    output_suffix = "rs",
    data_type = "resilience",
    resilience_variable = "rs",
    response = "value",
    extra_non_na = NULL,
    competition_covariate = "CI_BAsum_analysis",
    random_effect = "(1 | Site_analysis/Block_analysis)"
  ),
  list(
    trait_filter = "rc",
    output_trait = "rc",
    output_suffix = "rc",
    data_type = "resilience",
    resilience_variable = "rc",
    response = "value",
    extra_non_na = NULL,
    competition_covariate = "CI_BAsum_analysis",
    random_effect = "(1 | Site_analysis/Block_analysis)"
  ),
  list(
    trait_filter = "sla",
    output_trait = "sla",
    output_suffix = "sla",
    data_type = "growth",
    response = "SLA_cm2_mg",
    extra_non_na = "SLA_cm2_mg",
    competition_covariate = "CI_BAsum_analysis",
    random_effect = "(1 | Site_analysis/Block_analysis)"
  ),
  list(
    trait_filter = "bai",
    output_trait = "bai",
    output_suffix = "bai",
    data_type = "bai",
    response = "Estimate",
    extra_non_na = "Estimate",
    competition_covariate = "CI_BAsum_analysis",
    random_effect = "(1 | Site_analysis/Block_analysis)"
  ),
  list(
    trait_filter = "height",
    output_trait = "height",
    output_suffix = "height",
    data_type = "growth",
    response = "Height_m",
    extra_non_na = NULL,
    competition_covariate = "CI_BAsum_analysis",
    random_effect = "(1 | Site_analysis/Block_analysis)"
  ),
  list(
    trait_filter = "dbh",
    output_trait = "dbh",
    output_suffix = "dbh",
    data_type = "growth",
    response = "DBH_cm",
    extra_non_na = NULL,
    competition_covariate = "CI_BAsum_analysis",
    random_effect = "(1 | Site_analysis/Block_analysis)"
  )
)

analysis_configs <- list(
  list(
    name = "dist_prec_snp_tmean_as_cov",
    interaction_distance = "dist_prec",
    covariate_distance = "dist_tmean",
    output_prefix = "dist_prec_snp_tmean_as_cov_101224_oak_std_preds_snp_dist_tmean_dist_prec_as_int_results_all_"
  ),
  list(
    name = "dist_tmean_snp_prec_as_cov",
    interaction_distance = "dist_tmean",
    covariate_distance = "dist_prec",
    output_prefix = "dist_tmean_snp_prec_as_cov_101224_oak_std_preds_snp_dist_tmean_dist_prec_as_int_results_all_"
  )
)

# -----------------------------------------------------------------------------
# Run analyses
# -----------------------------------------------------------------------------

oak_bslmm_all_snps <- load_oak_bslmm_snp_table()
selected_univar_hits <- load_univar_genotypes()

for (analysis in analysis_configs) {
  for (config in trait_configs) {
    message("Running ", analysis$name, " for trait ", config$output_trait)
    
    run_trait_analysis(
      config = config,
      analysis = analysis,
      genotypes = selected_univar_hits,
      oak_bslmm_all_snps = oak_bslmm_all_snps
    )
  }
}
