
rm(list = ls())
gc()

library(dplyr)

# -----------------------------------------------------------------------------
# Paths
# -----------------------------------------------------------------------------

setwd("/path/to/working_directory")


read_selected_hits <- function(path) {
  read.csv(path, check.names = FALSE) %>%
    select(-any_of("X"))
}

beech_chromosome_names <- c(
  "1" = "Bhaga_1",
  "2" = "Bhaga_2",
  "3" = "Bhaga_3",
  "4" = "Bhaga_4",
  "5" = "Bhaga_5",
  "6" = "Bhaga_6",
  "7" = "Bhaga_7",
  "8" = "Bhaga_8",
  "9" = "Bhaga_9",
  "10" = "Bhaga_10",
  "11" = "Bhaga_11",
  "12" = "Bhaga_12"
)

selected_hit_files <- c(
  file.path(test_dir, "bai_beech_pc17_is_best_one.csv"),
  file.path(test_dir, "beech_gwas_results_dbh_270624.csv"),
  file.path(test_dir, "beech_gwas_results_height_270624.csv"),
  file.path(test_dir, "beech_gwas_results_sla_270624.csv"),
  file.path(test_dir, "beech_gwas_results_rt_270624.csv"),
  file.path(test_dir, "beech_gwas_results_rs_270624.csv"),
  file.path(test_dir, "beech_gwas_results_rc_270624.csv"),
  file.path(test_dir, "beech_gwas_results_corr_trait_191124.csv")
)

selected_hits <- bind_rows(lapply(selected_hit_files, read_selected_hits))

write.csv(
  selected_hits,
  file = file.path(test_dir, "all_snps_beech_results_280624.csv"),
  row.names = FALSE
)

unique_snp_positions <- selected_hits %>%
  mutate(chr_bp = paste(CHR, BP, sep = "_")) %>%
  select(CHR, BP, chr_bp) %>%
  arrange(chr_bp) %>%
  distinct(chr_bp, .keep_all = TRUE) %>%
  mutate(chr = beech_chromosome_names[as.character(CHR)]) %>%
  select(chr, BP)

write.table(
  unique_snp_positions,
  file = file.path(test_dir, "beech_unique_snps_to_check_ld_thinned_rel_mat_280624.txt"),
  quote = FALSE,
  row.names = FALSE
)
