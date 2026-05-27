
rm(list = ls())
gc()

library(dplyr)
library(tidyr)
library(ggplot2)
library(UpSetR)
library(pheatmap)

# -----------------------------------------------------------------------------
# Paths
# -----------------------------------------------------------------------------

oak_snp_dir <- "/path/to_oak_snp"
oak_trait_dir <- "/path_to_oak_trait"

beech_snp_dir <- "/path/to_beech_snp"
beech_trait_dir <- "/path_to_beech_trait"

output_dir <- "/path/to/output_dir"

y_metric <- "overlap_count"

# -----------------------------------------------------------------------------
# Species configuration
# -----------------------------------------------------------------------------

species_configs <- list(
  oak = list(
    snp_file = file.path(oak_snp_dir, "oak_all_bslmm_snps_table_151024.csv"),
    snp_sep = ",",
    trait_file = file.path(oak_trait_dir, "merged_df_for_corr_mat_all_traits_300325.csv"),
    site_col = "Site.x",
    trait_cols = c(
      "SLA_cm2_mg", "DBH_cm", "Height_m", "Estimate",
      "Rc", "Rt", "RS", "sumpre", "sprpre", "SUM", "SPR"
    ),
    name_map = c(
      rt = "Rt",
      rc = "Rc",
      prec_mean_spr = "sprpre",
      dbh = "DBH_cm",
      rs = "RS",
      tmean_sum = "SUM",
      bai = "Estimate",
      prec_mean_sum = "sumpre",
      tmean_spr = "SPR",
      height = "Height_m",
      sla = "SLA_cm2_mg"
    ),
    recode_traits = NULL
  ),
  beech = list(
    snp_file = file.path(beech_snp_dir, "snps_for_pleio_beech_101224_proper.csv"),
    snp_sep = ";",
    trait_file = file.path(beech_trait_dir, "merged_and_bai_slopes_rc_rs_rt_prec_tmean.csv"),
    site_col = "Site.x",
    trait_cols = c(
      "SLA_cm2_mg", "DBH_cm", "Height_m", "BAI_mm2",
      "Rc", "Rt", "Rs", "sumpre", "sprpre", "sumtmea", "sprtmea"
    ),
    name_map = c(
      rt = "Rt",
      rc = "Rc",
      prec_mean_spr = "sprpre",
      dbh = "DBH_cm",
      rs = "Rs",
      tmean_mean_sum = "sumtmea",
      bai = "BAI_mm2",
      prec_mean_sum = "sumpre",
      tmean_mean_spr = "sprtmea",
      height = "Height_m",
      sla = "SLA_cm2_mg"
    ),
    recode_traits = c(
      corr_prec_mean_spr = "prec_mean_spr_with_pc_chain1_pip01",
      corr_prec_mean_sum = "prec_mean_sum_with_pc_chain1_pip01",
      corr_tmean_mean_spr = "tmean_mean_spr_with_pc_chain1_pip01"
    )
  )
)

# -----------------------------------------------------------------------------
# GWAS/BSLMM loci helpers
# -----------------------------------------------------------------------------

read_snp_table <- function(config, species) {
  snp_table <- read.csv(config$snp_file, sep = config$snp_sep)

  if (!is.null(config$recode_traits)) {
    snp_table <- snp_table %>%
      mutate(trait = recode(trait, !!!config$recode_traits))
  }

  snp_table <- snp_table %>%
    separate(trait, into = c("trait", "pip"), sep = "chain", extra = "merge", fill = "right") %>%
    mutate(
      pip = ifelse(!is.na(pip), paste0("chain", pip), NA),
      trait = gsub(paste0(species, "_|_with_pc_"), "", trait)
    )

  snp_table
}

make_trait_snp_list <- function(snp_table) {
  trait_snp_list <- split(snp_table$chr_bp, snp_table$trait)
  lapply(trait_snp_list, function(x) unique(na.omit(x)))
}

rename_snp_traits <- function(trait_snp_list, name_map) {
  old_names <- names(trait_snp_list)
  names(trait_snp_list) <- ifelse(
    old_names %in% names(name_map),
    name_map[old_names],
    old_names
  )

  trait_snp_list
}

write_upset_plot <- function(trait_snp_list, species) {
  png(
    filename = file.path(output_dir, paste0(species, "_gwa_loci_upset.png")),
    width = 2800,
    height = 1800,
    res = 300
  )

  upset(
    fromList(trait_snp_list),
    order.by = "degree",
    nintersects = NA,
    nsets = 200
  )

  dev.off()
}

# -----------------------------------------------------------------------------
# Phenotypic trait-correlation helpers
# -----------------------------------------------------------------------------

read_standardized_traits <- function(config) {
  trait_data <- read.csv(config$trait_file)

  selected <- trait_data[, c(config$site_col, config$trait_cols)]

  standardized <- selected %>%
    group_by(.data[[config$site_col]]) %>%
    mutate(across(all_of(config$trait_cols), scale)) %>%
    ungroup()

  standardized[, config$trait_cols]
}

write_trait_correlation_heatmap <- function(correlation_matrix, species) {
  write.csv(
    correlation_matrix,
    file = file.path(output_dir, paste0(species, "_phenotypic_trait_correlation_matrix.csv"))
  )

  pheatmap(
    correlation_matrix,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    display_numbers = TRUE,
    color = colorRampPalette(c("blue", "white", "red"))(100),
    filename = file.path(output_dir, paste0(species, "_phenotypic_trait_correlation_heatmap.png")),
    width = 8,
    height = 7
  )
}

# -----------------------------------------------------------------------------
# Locus-overlap vs phenotypic-correlation helpers
# -----------------------------------------------------------------------------

make_pairwise_overlap_table <- function(trait_snp_list, traits) {
  pairs <- combn(traits, 2, simplify = FALSE)

  bind_rows(lapply(pairs, function(pair) {
    trait1 <- pair[1]
    trait2 <- pair[2]

    snps1 <- trait_snp_list[[trait1]]
    snps2 <- trait_snp_list[[trait2]]

    tibble(
      trait1 = trait1,
      trait2 = trait2,
      n_snps_trait1 = length(snps1),
      n_snps_trait2 = length(snps2),
      overlap_count = length(intersect(snps1, snps2))
    )
  }))
}

make_plot_df <- function(trait_snp_list, traits_data) {
  common_traits <- intersect(names(trait_snp_list), colnames(traits_data))

  trait_snp_list <- trait_snp_list[common_traits]
  traits_data <- traits_data[, common_traits, drop = FALSE]

  cor_mat <- cor(traits_data, use = "complete.obs", method = "pearson")

  cor_df <- as.data.frame(as.table(cor_mat)) %>%
    rename(trait1 = Var1, trait2 = Var2, r_pheno = Freq) %>%
    mutate(
      trait1 = as.character(trait1),
      trait2 = as.character(trait2)
    ) %>%
    filter(trait1 < trait2)

  overlap_df <- make_pairwise_overlap_table(trait_snp_list, common_traits)

  left_join(cor_df, overlap_df, by = c("trait1", "trait2"))
}

plot_locus_overlap_vs_trait_correlation <- function(plot_df, species) {
  ggplot(plot_df, aes(x = r_pheno, y = overlap_count)) +
    geom_jitter(width = 0, height = 0, alpha = 1) +
    labs(
      x = "Phenotypic correlation (Pearson r)",
      y = "Number of shared loci",
      title = paste("GWA loci vs phenotypic trait correlation,", species)
    ) +
    theme_classic()
}

write_spearman_test <- function(plot_df, species) {
  test <- cor.test(
    plot_df$r_pheno,
    plot_df[[y_metric]],
    method = "spearman",
    exact = FALSE
  )

  sink(file.path(output_dir, paste0(species, "_trait_correlation_vs_locus_overlap_spearman_test.txt")))
  print(test)
  sink()

  test
}

# -----------------------------------------------------------------------------
# Run per species
# -----------------------------------------------------------------------------

run_species_analysis <- function(species, config) {
  snp_table <- read_snp_table(config, species)
  trait_snp_list <- make_trait_snp_list(snp_table)

  write_upset_plot(trait_snp_list, species)

  trait_snp_list <- rename_snp_traits(trait_snp_list, config$name_map)

  traits_data <- read_standardized_traits(config)
  correlation_matrix <- cor(traits_data, use = "complete.obs", method = "pearson")

  write_trait_correlation_heatmap(correlation_matrix, species)

  plot_df <- make_plot_df(trait_snp_list, traits_data)

  write.csv(
    plot_df,
    file = file.path(output_dir, paste0(species, "_trait_correlation_vs_locus_overlap_table.csv")),
    row.names = FALSE
  )

  correlation_plot <- plot_locus_overlap_vs_trait_correlation(plot_df, species)

  ggsave(
    filename = file.path(output_dir, paste0(species, "_trait_correlation_vs_locus_overlap.png")),
    plot = correlation_plot,
    width = 7,
    height = 5,
    dpi = 300
  )

  test <- write_spearman_test(plot_df, species)

  list(
    snp_table = snp_table,
    trait_snp_list = trait_snp_list,
    traits_data = traits_data,
    correlation_matrix = correlation_matrix,
    plot_df = plot_df,
    correlation_plot = correlation_plot,
    spearman_test = test
  )
}

oak_results <- run_species_analysis("oak", species_configs$oak)
beech_results <- run_species_analysis("beech", species_configs$beech)
