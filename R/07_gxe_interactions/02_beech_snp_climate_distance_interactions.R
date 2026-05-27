
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

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# -----------------------------------------------------------------------------
# Helpers
# -----------------------------------------------------------------------------

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

make_chr_bp_snp_id <- function(snp_table, trait_name = NA_character_) {
  if (is.null(snp_table) || nrow(snp_table) == 0) {
    warning("No SNP rows found for trait: ", trait_name)
    return(character(0))
  }
  
  if (!"chr_bp" %in% colnames(snp_table)) {
    stop("SNP table does not contain required column 'chr_bp'. Available columns: ",
         paste(colnames(snp_table), collapse = ", "))
  }
  
  snp_table$snp <- "snp"
  snp_table$snp_id <- paste(snp_table$snp, snp_table$chr_bp, sep = "_")
  
  as.character(snp_table$snp_id)
}

read_vcf_sample_table <- function() {
  samples_in_vcf <- read.table(
    file.path(phenotype_dir, "samples_populations_fagus_no_high_miss_ind_and_incorr_removed_minDP3_meanDP10_maxDP80_NA0.90_reduced.MAF0.05.txt"),
    sep = "\t",
    header = FALSE
  )
  
  split_ids <- strsplit(as.character(samples_in_vcf$V1), "-")
  samples_in_vcf <- cbind(samples_in_vcf, do.call(rbind, split_ids))
  
  colnames(samples_in_vcf)[1] <- "SampleID_vcf"
  colnames(samples_in_vcf)[2] <- "RowNumber"
  colnames(samples_in_vcf)[3] <- "Treeid"
  colnames(samples_in_vcf)[4] <- "seq_id"
  colnames(samples_in_vcf)[5] <- "vcf_extra"
  
  samples_in_vcf
}

load_beech_bslmm_snp_table <- function(file_name, use_trait_real = TRUE) {
  snps <- read.csv(file.path(bslmm_dir, file_name), sep = ",")
  
  if (!"trait" %in% colnames(snps)) {
    stop("SNP table ", file_name, " does not contain required column 'trait'.")
  }
  
  # Always make robust trait labels.
  # This avoids zero-row matches when the file stores traits as:
  #   beech_prec_mean_spr
  #   beech_prec_mean_spr_with_...
  #   prec_mean_spr
  #   or older trait_real-style labels.
  snps <- snps %>%
    mutate(
      trait_raw = as.character(trait),
      trait_no_prefix = trait_raw %>%
        str_replace("^beech_", "") %>%
        str_replace("_with.*", ""),
      trait_real = case_when(
        str_detect(trait_raw, "beech_tmean|beech_prec") ~ str_extract(trait_raw, "(?<=beech_)(\\w+_\\w+_?\\w*)"),
        str_detect(trait_raw, "^beech_") ~ str_extract(trait_raw, "(?<=beech_)(\\w+)"),
        TRUE ~ trait_no_prefix
      ),
      trait_real = str_replace(trait_real, "_with.*", "")
    )
  
  snps
}

load_non_corr_genotypes <- function() {
  read.csv(file.path(bslmm_dir, "extracted_beech_bslmm_genos_151024.csv"), sep = ",")
}

load_corr_genotypes <- function() {
  read.csv(file.path(bslmm_dir, "extracted_univar_genos_only_corr_beech_041224.csv"), sep = ",")
}

prepare_dbh_height_sla_data <- function(genotypes) {
  samples_in_vcf <- read_vcf_sample_table()
  beech <- read.csv(file.path(phenotype_dir, "dryad_beech_dbh.csv"), sep = ",")
  
  merged <- merge(samples_in_vcf, beech, by = "Treeid")
  merged <- merged[!is.na(merged$CI_BAsum), ]
  
  merge(merged, genotypes, by.x = "SampleID_vcf", by.y = "SampleID")
}

prepare_bai_data <- function(genotypes) {
  samples_in_vcf <- read_vcf_sample_table()
  beech_bai <- read.csv(file.path(phenotype_dir, "dryad_beech_bai.csv"))
  
  merged <- merge(samples_in_vcf, beech_bai, by = "Treeid")
  
  fit_lm_and_extract_estimate <- function(data) {
    lm_model <- lm(BAI_mm2 ~ year, data = data)
    data.frame(
      Estimate = coef(lm_model)[2],
      SampleID_vcf = unique(data$SampleID_vcf)
    )
  }
  
  bai_slopes_as_trait <- merged %>%
    group_by(SampleID_vcf) %>%
    do(fit_lm_and_extract_estimate(.))
  
  merged_and_bai_slopes <- merge(bai_slopes_as_trait, merged, by = "SampleID_vcf")
  merged_and_bai_slopes <- merged_and_bai_slopes %>%
    distinct(SampleID_vcf, .keep_all = TRUE)
  
  merged_and_bai_slopes$Block_for_model <- merged_and_bai_slopes$Blocknew
  
  merge(merged_and_bai_slopes, genotypes, by.x = "SampleID_vcf", by.y = "SampleID")
}

prepare_resilience_data <- function(resilience_variable, genotypes) {
  samples_in_vcf <- read_vcf_sample_table()
  beech_rs <- read.csv(file.path(phenotype_dir, "beech_rs.csv"))
  
  merged <- merge(samples_in_vcf, beech_rs, by = "Treeid")
  trait_data <- subset(merged, var == resilience_variable)
  
  merge(trait_data, genotypes, by.x = "SampleID_vcf", by.y = "SampleID")
}

prepare_corr_trait_data <- function(climate_variable, genotypes) {
  samples_in_vcf <- read_vcf_sample_table()
  beech_corr <- read.csv(file.path(climate_distance_dir, "dryad_beech_corr_last.csv"))
  
  merged <- merge(samples_in_vcf, beech_corr, by = "Treeid")
  trait_data <- subset(merged, varname == climate_variable)
  
  merge(trait_data, genotypes, by.x = "SampleID_vcf", by.y = "SampleID")
}

load_distance_table <- function(file_name) {
  distance_table <- read.csv(file.path(climate_distance_dir, file_name))
  
  colnames(distance_table)[1] <- "Prov"
  colnames(distance_table)[2] <- "ES"
  colnames(distance_table)[3] <- "DE"
  colnames(distance_table)[4] <- "UK"
  colnames(distance_table)[5] <- "SE"
  
  distance_table %>%
    dplyr::select(-ID) %>%
    mutate(
      Prov = case_when(
        Prov == "UK018" ~ "GB018",
        TRUE ~ Prov
      )
    ) %>%
    pivot_longer(cols = -Prov, names_to = "country", values_to = "dist") %>%
    mutate(merger = paste(Prov, country, sep = "_")) %>%
    dplyr::select(Prov, dist, merger)
}

add_climate_distances <- function(selected_df) {
  selected_df$merger <- paste(selected_df$Prov, selected_df$Site, sep = "_")
  
  dist_prec <- load_distance_table("distance_df_beech_std_prec_annual_051224.csv")
  dist_tmin <- load_distance_table("distance_df_beech_std_tmin_annual_051224.csv")
  dist_tmax <- load_distance_table("distance_df_beech_std_tmax_annual_051224.csv")
  dist_tmean <- load_distance_table("distance_df_beech_std_tmean_annual_051224.csv")
  
  selected_df_with_clima <- merge(selected_df, dist_prec, by = "merger")
  colnames(selected_df_with_clima)[colnames(selected_df_with_clima) == "dist"] <- "dist_prec"
  
  selected_df_with_clima <- merge(selected_df_with_clima, dist_tmin, by = "merger")
  colnames(selected_df_with_clima)[colnames(selected_df_with_clima) == "dist"] <- "dist_tmin"
  
  selected_df_with_clima <- merge(selected_df_with_clima, dist_tmax, by = "merger")
  colnames(selected_df_with_clima)[colnames(selected_df_with_clima) == "dist"] <- "dist_tmax"
  
  selected_df_with_clima <- merge(selected_df_with_clima, dist_tmean, by = "merger")
  colnames(selected_df_with_clima)[colnames(selected_df_with_clima) == "dist"] <- "dist_tmean"
  
  list(
    data = selected_df_with_clima,
    dist_tmean = dist_tmean
  )
}

get_trait_data <- function(config, non_corr_genotypes, corr_genotypes) {
  if (config$data_type == "bai") {
    prepare_bai_data(non_corr_genotypes)
  } else if (config$data_type == "growth") {
    prepare_dbh_height_sla_data(non_corr_genotypes)
  } else if (config$data_type == "resilience") {
    prepare_resilience_data(config$resilience_variable, non_corr_genotypes)
  } else if (config$data_type == "corr") {
    prepare_corr_trait_data(config$climate_variable, corr_genotypes)
  }
}

extract_emtrend_results <- function(snp, final_results_tmean, coefs_m1) {
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
    predictor = "dist_tmean",
    coefs_m1_interaction_id_snp_dist_tmean = coefs_m1[6, 1],
    coefs_m1_interaction_snp_dist_tmean_estimate = coefs_m1[6, 2],
    coefs_m1_interaction_snp_dist_tmean_SE = coefs_m1[6, 3],
    coefs_m1_interaction_snp_dist_tmean_p.value = coefs_m1[6, 4]
  )
}

calculate_af_difference <- function(data, dist_tmean, snp) {
  selected_df_for_af <- data
  selected_df_for_af$merger <- paste(selected_df_for_af$Prov, selected_df_for_af$Site, sep = "_")
  
  merged_for_af <- merge(dist_tmean, selected_df_for_af, by = "merger")
  
  filtered_df <- merged_for_af %>%
    filter(Prov.x %in% c("CZ048", "FR008"))
  
  allele_frequencies <- filtered_df %>%
    group_by(Prov.x) %>%
    summarise(
      freq_A = calculate_allele_freq(!!sym(snp))[1],
      freq_B = calculate_allele_freq(!!sym(snp))[2],
      .groups = "drop"
    )
  
  af_difference <- allele_frequencies[1, 2] - allele_frequencies[2, 2]
  af_difference[1, 1]
}

analyse_snp_trait <- function(data, snp, config, analysis) {
  if (!snp %in% colnames(data)) {
    warning("Skipping SNP ", snp, " because it is not present in the genotype/trait data.")
    return(NULL)
  }
  
  data[[snp]] <- as.factor(data[[snp]])
  
  data_no_bb <- data %>%
    filter(!!sym(snp) != "BB") %>%
    filter(!is.na(!!sym(snp))) %>%
    filter(!is.na(.data[[config$response]]))
  
  mean_AA <- mean(data_no_bb[[config$response]][data_no_bb[[snp]] == "AA"])
  mean_AB <- mean(data_no_bb[[config$response]][data_no_bb[[snp]] == "AB"])
  
  if (mean_AB > mean_AA) {
    data_no_bb[[snp]] <- ifelse(data_no_bb[[snp]] == "AA", "AB", "AA")
  }
  
  data_no_bb <- data_no_bb %>%
    mutate(rt_standardized = (.data[[config$response]] - mean(.data[[config$response]], na.rm = TRUE)) / sd(.data[[config$response]], na.rm = TRUE))
  
  selected_df <- data_no_bb
  climate_data <- add_climate_distances(selected_df)
  selected_df_with_clima <- climate_data$data
  
  selected_df_with_clima[[snp]] <- as.factor(selected_df_with_clima[[snp]])
  
  model_formula_m1 <- as.formula(
    paste(
      "rt_standardized ~",
      paste0(
        snp,
        " * ",
        analysis$interaction_distance,
        " + ",
        analysis$covariate_distance,
        " + CI_BAsum + (1 | Site/Block)"
      )
    )
  )
  
  m1 <- lmer(model_formula_m1, data = selected_df_with_clima)
  
  summary_m1 <- summary(m1)
  m1_coefficients_table <- summary_m1$coefficients
  coefs_m1 <- data.frame(
    type = rownames(m1_coefficients_table),
    Estimate = m1_coefficients_table[, "Estimate"],
    Std_Error = m1_coefficients_table[, "Std. Error"],
    P_Value = m1_coefficients_table[, "Pr(>|t|)"]
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
    coefs_m1 = coefs_m1
  )
  
  transformed_final_df$af_difference <- calculate_af_difference(
    data = data,
    dist_tmean = climate_data$dist_tmean,
    snp = snp
  )
  
  transformed_final_df
}

match_snp_table_for_trait <- function(snp_table, config) {
  trait_filter <- config$trait_filter
  
  exact_match <- snp_table %>%
    filter(
      trait_raw == trait_filter |
        trait_no_prefix == trait_filter |
        trait_real == trait_filter
    )
  
  if (nrow(exact_match) > 0) {
    return(exact_match)
  }
  
  # Fallback: match trait_filter as a fixed substring of the raw trait column.
  # This catches values like "beech_prec_mean_spr_with_covariates".
  substring_match <- snp_table %>%
    filter(str_detect(trait_raw, fixed(trait_filter)))
  
  if (nrow(substring_match) > 0) {
    warning(
      "Used substring SNP-trait match for trait '", trait_filter,
      "'. Matched ", nrow(substring_match), " SNP rows."
    )
    return(substring_match)
  }
  
  available_traits <- sort(unique(c(
    snp_table$trait_raw,
    snp_table$trait_no_prefix,
    snp_table$trait_real
  )))
  available_traits <- available_traits[!is.na(available_traits)]
  
  warning(
    "No SNPs matched trait '", trait_filter, "' in file '", config$snp_table_file, "'. ",
    "First available trait labels are: ",
    paste(utils::head(available_traits, 25), collapse = ", ")
  )
  
  snp_table[0, , drop = FALSE]
}

run_trait_analysis <- function(config, analysis, non_corr_genotypes, corr_genotypes) {
  snp_table <- load_beech_bslmm_snp_table(
    file_name = config$snp_table_file,
    use_trait_real = config$use_trait_real
  )
  
  snp_table <- match_snp_table_for_trait(snp_table, config)
  
  snp_list <- make_chr_bp_snp_id(snp_table, trait_name = config$trait_filter)
  data <- get_trait_data(config, non_corr_genotypes, corr_genotypes)
  
  final_results_all_snps <- data.frame()
  
  if (length(snp_list) == 0) {
    warning("Skipping trait '", config$output_trait, "' because no SNPs were found after matching.")
    final_results_all_snps$trait <- config$output_trait
    
    write.csv(
      final_results_all_snps,
      file = file.path(output_dir, paste0(analysis$output_prefix, config$output_suffix, config$output_date_suffix, output_date, ".csv"))
    )
    
    return(final_results_all_snps)
  }
  
  message("  SNPs matched for trait ", config$output_trait, ": ", length(snp_list))
  
  for (snp in snp_list) {
    transformed_final_df <- tryCatch(
      analyse_snp_trait(
        data = data,
        snp = snp,
        config = config,
        analysis = analysis
      ),
      error = function(e) {
        warning("Skipping SNP ", snp, " for trait ", config$output_trait, ": ", conditionMessage(e))
        NULL
      }
    )
    
    if (!is.null(transformed_final_df) && nrow(transformed_final_df) > 0) {
      final_results_all_snps <- rbind(final_results_all_snps, transformed_final_df)
    }
  }
  
  final_results_all_snps$trait <- config$output_trait
  
  write.csv(
    final_results_all_snps,
    file = file.path(output_dir, paste0(analysis$output_prefix, config$output_suffix, config$output_date_suffix, output_date, ".csv"))
  )
  
  final_results_all_snps
}

# -----------------------------------------------------------------------------
# Analysis configuration
# -----------------------------------------------------------------------------

trait_configs <- list(
  list(
    trait_filter = "bai",
    output_trait = "bai",
    output_suffix = "bai",
    output_date_suffix = "_041224",
    snp_table_file = "beech_all_bslmm_snps_table_151024.csv",
    use_trait_real = TRUE,
    data_type = "bai",
    response = "Estimate"
  ),
  list(
    trait_filter = "dbh",
    output_trait = "dbh",
    output_suffix = "dbh",
    output_date_suffix = "_041224",
    snp_table_file = "beech_all_bslmm_snps_table_151024.csv",
    use_trait_real = TRUE,
    data_type = "growth",
    response = "DBH_cm"
  ),
  list(
    trait_filter = "height",
    output_trait = "height",
    output_suffix = "height",
    output_date_suffix = "_041224",
    snp_table_file = "beech_all_bslmm_snps_table_151024.csv",
    use_trait_real = TRUE,
    data_type = "growth",
    response = "Height_m"
  ),
  list(
    trait_filter = "prec_mean_spr",
    output_trait = "prec_mean_spr",
    output_suffix = "prec_mean_spr",
    output_date_suffix = "_051224",
    snp_table_file = "beech_all_bslmm_snps_table_261124_051224.csv",
    use_trait_real = FALSE,
    data_type = "corr",
    climate_variable = "prec.mean",
    response = "SPR"
  ),
  list(
    trait_filter = "prec_mean_sum",
    output_trait = "prec_mean_sum",
    output_suffix = "prec_mean_sum",
    output_date_suffix = "_051224",
    snp_table_file = "beech_all_bslmm_snps_table_261124_051224.csv",
    use_trait_real = FALSE,
    data_type = "corr",
    climate_variable = "prec.mean",
    response = "SUM"
  ),
  list(
    trait_filter = "rc",
    output_trait = "rc",
    output_suffix = "rc",
    output_date_suffix = "_041224",
    snp_table_file = "beech_all_bslmm_snps_table_151024.csv",
    use_trait_real = TRUE,
    data_type = "resilience",
    resilience_variable = "rc",
    response = "value"
  ),
  list(
    trait_filter = "rs",
    output_trait = "rs",
    output_suffix = "rs",
    output_date_suffix = "_041224",
    snp_table_file = "beech_all_bslmm_snps_table_151024.csv",
    use_trait_real = TRUE,
    data_type = "resilience",
    resilience_variable = "rs",
    response = "value"
  ),
  list(
    trait_filter = "rt",
    output_trait = "rt",
    output_suffix = "rt",
    output_date_suffix = "_041224",
    snp_table_file = "beech_all_bslmm_snps_table_151024.csv",
    use_trait_real = TRUE,
    data_type = "resilience",
    resilience_variable = "rt",
    response = "value"
  ),
  list(
    trait_filter = "sla",
    output_trait = "sla",
    output_suffix = "sla",
    output_date_suffix = "_041224",
    snp_table_file = "beech_all_bslmm_snps_table_151024.csv",
    use_trait_real = TRUE,
    data_type = "growth",
    response = "SLA_cm2_mg"
  ),
  list(
    trait_filter = "tmean_mean_spr",
    output_trait = "tmean_mean_spr",
    output_suffix = "tmean_mean_spr",
    output_date_suffix = "_051224",
    snp_table_file = "beech_all_bslmm_snps_table_261124_051224.csv",
    use_trait_real = FALSE,
    data_type = "corr",
    climate_variable = "tmea.mean",
    response = "SPR"
  ),
  list(
    trait_filter = "tmean_mean_sum",
    output_trait = "tmean_mean_sum",
    output_suffix = "tmean_mean_sum",
    output_date_suffix = "_051224",
    snp_table_file = "beech_all_bslmm_snps_table_261124_051224.csv",
    use_trait_real = FALSE,
    data_type = "corr",
    climate_variable = "tmea.mean",
    response = "SUM"
  )
)

analysis_configs <- list(
  list(
    name = "dist_prec_snp_tmean_as_cov",
    interaction_distance = "dist_prec",
    covariate_distance = "dist_tmean",
    output_prefix = "dist_prec_snp_tmean_as_cov_101224_beech_snp_dist_STD_tmean_dist_prec_as_int_results_all_"
  ),
  list(
    name = "dist_tmean_snp_prec_as_cov",
    interaction_distance = "dist_tmean",
    covariate_distance = "dist_prec",
    output_prefix = "dist_tmean_snp_prec_as_cov_101224_beech_snp_dist_STD_tmean_dist_prec_as_int_results_all_"
  )
)

# -----------------------------------------------------------------------------
# Run analyses
# -----------------------------------------------------------------------------

non_corr_genotypes <- load_non_corr_genotypes()
corr_genotypes <- load_corr_genotypes()

for (analysis in analysis_configs) {
  for (config in trait_configs) {
    message("Running ", analysis$name, " for trait ", config$output_trait)
    
    run_trait_analysis(
      config = config,
      analysis = analysis,
      non_corr_genotypes = non_corr_genotypes,
      corr_genotypes = corr_genotypes
    )
  }
}

