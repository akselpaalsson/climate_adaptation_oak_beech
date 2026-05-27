
rm(list = ls())
gc()

library(dplyr)

base_dir <- "path/to/base_dir/"

# -----------------------------------------------------------------------------
# Trait/PIP configuration
# -----------------------------------------------------------------------------

snp_sets <- list(
  list(
    directory = "dbh",
    prefix = "dbh_with_pc",
    pip = "01",
    representative_chain = 4,
    trait = "dbh_with_pc_chain4_pip01_pip01_subset"
  ),
  list(
    directory = "rc",
    prefix = "rc_with_pc",
    pip = "01",
    representative_chain = 1,
    trait = "rc_with_pc_chain1_pip01_pip01_subset"
  ),
  list(
    directory = "rt",
    prefix = "rt_with_pc",
    pip = "01",
    representative_chain = 1,
    trait = "rt_with_pc_chain1_pip01_pip01_subset"
  ),
  list(
    directory = "rs",
    prefix = "rs_with_pc",
    pip = "01",
    representative_chain = 1,
    trait = "rs_with_pc_chain1_pip01_pip01_subset"
  ),
  list(
    directory = "rs",
    prefix = "rs_with_pc",
    pip = "10",
    representative_chain = 1,
    trait = "rs_with_pc_chain1_pip10_pip10_subset"
  ),
  list(
    directory = "tmean_sum",
    prefix = "tmean_sum_with_pc",
    pip = "01",
    representative_chain = 3,
    trait = "tmean_sum_with_pc_chain3_pip01_pip01_subset"
  ),
  list(
    directory = "tmean_spr",
    prefix = "tmean_spr_with_pc",
    pip = "01",
    representative_chain = 4,
    trait = "tmean_spr_with_pc_chain4_pip01_pip01_subset"
  ),
  list(
    directory = "tmean_spr",
    prefix = "tmean_spr_with_pc",
    pip = "10",
    representative_chain = 4,
    trait = "tmean_spr_with_pc_chain4_pip10_pip10_subset"
  ),
  list(
    directory = "sla",
    prefix = "sla_with_pc",
    pip = "01",
    representative_chain = 4,
    trait = "sla_with_pc_chain4_pip01_pip01_subset"
  ),
  list(
    directory = "sla",
    prefix = "sla_with_pc",
    pip = "10",
    representative_chain = 4,
    trait = "sla_with_pc_chain4_pip10_pip10_subset"
  ),
  list(
    directory = "prec_mean_sum",
    prefix = "prec_mean_sum_with_pc",
    pip = "01",
    representative_chain = 1,
    trait = "prec_mean_sum_with_pc_chain1_pip01_pip01_subset"
  ),
  list(
    directory = "prec_mean_spr",
    prefix = "prec_mean_spr_with_pc",
    pip = "01",
    representative_chain = 2,
    trait = "prec_mean_spr_with_pc_chain2_pip01_pip01_subset"
  ),
  list(
    directory = "height",
    prefix = "height_with_pc",
    pip = "01",
    representative_chain = 2,
    trait = "height_with_pc_chain2_pip01_pip01_subset"
  ),
  list(
    directory = "height",
    prefix = "height_with_pc",
    pip = "10",
    representative_chain = 2,
    trait = "height_with_pc_chain2_pip10_pip10_subset"
  ),
  list(
    directory = "height",
    prefix = "height_with_pc",
    pip = "25",
    representative_chain = 2,
    trait = "height_with_pc_chain2_pip25_pip25_subset"
  ),
  list(
    directory = "bai",
    prefix = "bai_with_pc",
    pip = "01",
    representative_chain = 4,
    trait = "bai_with_pc_chain4_pip01_pip01_subset"
  ),
  list(
    directory = "bai",
    prefix = "bai_with_pc",
    pip = "10",
    representative_chain = 4,
    trait = "bai_with_pc_chain4_pip10_pip10_subset"
  ),
  list(
    directory = "bai",
    prefix = "bai_with_pc",
    pip = "25",
    representative_chain = 4,
    trait = "bai_with_pc_chain4_pip25_pip25_subset"
  ),
  list(
    directory = "bai",
    prefix = "bai_with_pc",
    pip = "50",
    representative_chain = 4,
    trait = "bai_with_pc_chain4_pip50_pip50_subset"
  )
)

# -----------------------------------------------------------------------------
# Helpers
# -----------------------------------------------------------------------------

read_chain_pip_table <- function(directory, prefix, chain, pip) {
  file_name <- paste0(prefix, "_chain", chain, "_pip", pip, ".csv")

  df <- read.csv(file.path(base_dir, directory, file_name))
  df$chr_bp <- paste(df$CHR, df$BP, sep = "_")

  df
}

find_common_chain_snps <- function(chain_tables) {
  Reduce(intersect, lapply(chain_tables, function(df) df$chr_bp))
}

make_common_snp_subset <- function(config) {
  chain_tables <- lapply(
    1:4,
    function(chain) {
      read_chain_pip_table(
        directory = config$directory,
        prefix = config$prefix,
        chain = chain,
        pip = config$pip
      )
    }
  )

  common_chr_bp <- find_common_chain_snps(chain_tables)

  subset_df <- chain_tables[[config$representative_chain]] %>%
    filter(chr_bp %in% common_chr_bp)

  subset_df$pip <- config$pip
  subset_df$trait <- config$trait

  subset_df
}

# -----------------------------------------------------------------------------
# Compile final SNP table
# -----------------------------------------------------------------------------

final_df <- bind_rows(lapply(snp_sets, make_common_snp_subset))

final_df$chr_bp <- paste(final_df$CHR, final_df$BP, sep = "_")

overlapping_entries <- final_df %>%
  group_by(chr_bp) %>%
  summarise(unique_traits = n_distinct(trait), .groups = "drop") %>%
  filter(unique_traits > 1)

overlapping_table <- final_df %>%
  filter(chr_bp %in% overlapping_entries$chr_bp)

write.csv(
  final_df,
  file = file.path(base_dir, "bai", "oak_all_bslmm_snps_table_151024.csv")
)
