
rm(list = ls())
gc()

library(dplyr)
library(tidyr)

# -----------------------------------------------------------------------------
# Paths
# -----------------------------------------------------------------------------

setwd("/path/to/working_directory")

# For this test version, all LFMM, GWAS/BSLMM, LD-marked, and vcftools LD files
# are expected directly in the shared test folder.
oak_lfmm_dir <- test_dir
beech_lfmm_prec_dir <- test_dir
beech_lfmm_tmean_dir <- test_dir
ld_check_dir <- test_dir

oak_bslmm_dir <- test_dir
beech_bslmm_dir <- test_dir
beech_corr_gwas_dir <- test_dir

# -----------------------------------------------------------------------------
# Chromosome mappings
# -----------------------------------------------------------------------------

oak_numeric_to_scaffold <- c(
  "1" = "OW028765.1",
  "2" = "OW028766.1",
  "3" = "OW028767.1",
  "4" = "OW028768.1",
  "5" = "OW028769.1",
  "6" = "OW028773.1",
  "7" = "OW028770.1",
  "8" = "OW028771.1",
  "9" = "OW028772.1",
  "10" = "OW028774.1",
  "11" = "OW028775.1",
  "12" = "OW028776.1"
)

beech_numeric_to_scaffold <- c(
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
  "12" = "Bhaga_12",
  "14" = "Bhaga_Mitochondria_Circular"
)

# -----------------------------------------------------------------------------
# Helpers
# -----------------------------------------------------------------------------

read_csv_from <- function(directory, file_name, sep = ",") {
  read.csv(file.path(directory, file_name), sep = sep)
}

write_table_to <- function(data, directory, file_name) {
  write.table(
    data,
    file = file.path(directory, file_name),
    row.names = FALSE,
    quote = FALSE
  )
}

make_lfmm_position_file <- function(tmean, prec, split_pattern, output_dir, output_file) {
  lfmm <- rbind(tmean, prec)

  lfmm <- lfmm %>%
    separate(
      nonsense,
      into = c("CHR", "BP"),
      sep = split_pattern,
      extra = "merge",
      remove = TRUE
    )

  top_lfmm <- subset(lfmm, top_1_percent == "yes")
  top_lfmm <- subset(top_lfmm, select = c(CHR, BP))
  top_lfmm_unique <- top_lfmm[!duplicated(top_lfmm[, c("CHR", "BP")]), ]
  top_lfmm_unique <- na.omit(top_lfmm_unique)

  write_table_to(top_lfmm_unique, output_dir, output_file)

  top_lfmm_unique
}

make_top1_lfmm_table <- function(tmean, prec, output_file) {
  prec$climate_variable <- "prec"
  tmean$climate_variable <- "tmean"

  lfmm <- rbind(prec, tmean)
  lfmm <- subset(lfmm, top_1_percent == "yes")

  write.csv(lfmm, file = file.path(ld_check_dir, output_file))

  lfmm
}

read_ld_marked_lfmm <- function(file_name) {
  file_path <- file.path(ld_check_dir, file_name)

  # The original script had a beech filename with a space before .csv.
  # Prefer the cleaned filename, but fall back to the original spelling if needed.
  if (!file.exists(file_path)) {
    alt_file_name <- gsub(" \.csv$", ".csv", file_name)
    alt_file_path <- file.path(ld_check_dir, alt_file_name)

    if (file.exists(alt_file_path)) {
      file_path <- alt_file_path
    } else {
      stop("LD-marked LFMM file not found: ", file_path,
           " or cleaned alternative: ", alt_file_path)
    }
  }

  lfmm <- read.csv(file_path, sep = ";", check.names = FALSE)

  if ("X" %in% names(lfmm)) {
    lfmm <- subset(lfmm, X != "remove")
  }

  lfmm
}

make_oak_gwas_positions <- function() {
  oak_snps <- read.csv(file.path(oak_bslmm_dir, "oak_all_bslmm_snps_table_151024.csv"))
  oak_snps$CHR <- oak_numeric_to_scaffold[as.character(oak_snps$CHR)]
  oak_snps <- oak_snps[, c("CHR", "BP")]

  oak_snps
}

make_beech_gwas_positions <- function() {
  beech_snps_no_corr <- read.csv(file.path(beech_bslmm_dir, "beech_all_bslmm_snps_table_151024.csv"))
  beech_snps_to_get_corr <- read.csv(
    file.path(beech_corr_gwas_dir, "snps_for_pleio_beech_101224_proper.csv"),
    sep = ";"
  )

  traits_to_keep <- c(
    "prec_mean_spr_with_pc_chain1_pip01",
    "corr_prec_mean_sum",
    "tmean_mean_spr_with_pc_chain1_pip01",
    "tmean_mean_sum_with_pc_chain1_pip01",
    "prec_mean_sum_with_pc_chain1_pip01",
    "corr_prec_mean_spr",
    "corr_tmean_mean_spr",
    "beech_prec_mean_spr_with_pc_chain1_pip01_subset"
  )

  traits_to_exclude <- c(
    "beech_tmean_mean_spr_with_pc_chain1_pip01_subset",
    "beech_tmean_mean_sum_with_pc_chain1_pip10_subset",
    "beech_tmean_mean_sum_with_pc_chain1_pip01_subset",
    "beech_prec_mean_spr_with_pc_chain1_pip01_subset"
  )

  beech_snps_filtered <- beech_snps_to_get_corr %>%
    filter(trait %in% traits_to_keep)

  beech_snps_no_corr_filtered <- beech_snps_no_corr %>%
    filter(!trait %in% traits_to_exclude)

  beech_snps_no_corr_filtered <- beech_snps_no_corr_filtered[, c("CHR", "BP")]
  beech_snps_filtered <- beech_snps_filtered[, c("CHR", "BP")]

  beech_snps_gwas <- rbind(beech_snps_no_corr_filtered, beech_snps_filtered)
  beech_snps_gwas$CHR <- beech_numeric_to_scaffold[as.character(beech_snps_gwas$CHR)]

  beech_snps_gwas
}

prepare_combined_oak_positions <- function(final_oak_lfmm, oak_snps, add_origin = FALSE) {
  lfmm <- final_oak_lfmm[, c("CHR_id", "BP")]
  colnames(lfmm)[1] <- "CHR"

  if (add_origin) {
    lfmm$origin <- "lfmm"
    oak_snps$origin <- "gwas"
  }

  rbind(lfmm, oak_snps)
}

prepare_combined_beech_positions <- function(final_beech_lfmm, beech_snps_gwas, add_origin = FALSE) {
  lfmm <- final_beech_lfmm[, c("nonsense_", "BP")]
  colnames(lfmm)[1] <- "CHR"

  if (add_origin) {
    lfmm$origin <- "lfmm"
    beech_snps_gwas$origin <- "gwas"
  }

  rbind(lfmm, beech_snps_gwas)
}

summarize_high_ld_pairs <- function(ld_table, snp_table) {
  required_ld_cols <- c("CHR", "POS1", "POS2", "R.2")
  missing_ld_cols <- setdiff(required_ld_cols, names(ld_table))
  if (length(missing_ld_cols) > 0) {
    stop("LD table is missing required column(s): ", paste(missing_ld_cols, collapse = ", "))
  }

  required_snp_cols <- c("CHR", "BP", "origin")
  missing_snp_cols <- setdiff(required_snp_cols, names(snp_table))
  if (length(missing_snp_cols) > 0) {
    stop("SNP table is missing required column(s): ", paste(missing_snp_cols, collapse = ", "))
  }

  snp_lookup <- snp_table %>%
    mutate(snp_id_col = paste(CHR, BP, sep = "_")) %>%
    select(snp_id_col, origin)

  ld_table <- ld_table %>%
    mutate(
      snpid_snp1 = paste(CHR, POS1, sep = "_"),
      snpid_snp2 = paste(CHR, POS2, sep = "_")
    )

  merged_snp1 <- merge(
    ld_table,
    snp_lookup,
    by.x = "snpid_snp1",
    by.y = "snp_id_col"
  ) %>%
    rename(snp1_origin = origin)

  merged_snp2 <- merge(
    merged_snp1,
    snp_lookup,
    by.x = "snpid_snp2",
    by.y = "snp_id_col"
  ) %>%
    rename(snp2_origin = origin)

  high_ld <- subset(merged_snp2, R.2 > 0.6)
  high_ld[high_ld$snp1_origin != high_ld$snp2_origin, ]
}

# -----------------------------------------------------------------------------
# 1. Prepare LFMM-only top-1% position files for vcftools
# -----------------------------------------------------------------------------

oak_lfmm_tmean <- read_csv_from(oak_lfmm_dir, "oak_lfmm_tmean_top1_identified_121224.csv")
oak_lfmm_prec <- read_csv_from(oak_lfmm_dir, "oak_lfmm_prec_top1_identified_121224.csv")

make_lfmm_position_file(
  tmean = oak_lfmm_tmean,
  prec = oak_lfmm_prec,
  split_pattern = "_",
  output_dir = oak_lfmm_dir,
  output_file = "oak_lfmm_tmean_and_prec_to_check_ld_121224.txt"
)

beech_lfmm_prec <- read_csv_from(beech_lfmm_prec_dir, "beech_lfmm_prec_top1_identified_121224.csv")
beech_lfmm_tmean <- read_csv_from(beech_lfmm_tmean_dir, "beech_lfmm_tmean_top1_identified_121224.csv")

make_lfmm_position_file(
  tmean = beech_lfmm_tmean,
  prec = beech_lfmm_prec,
  split_pattern = "_(?=[^_]+$)",
  output_dir = beech_lfmm_tmean_dir,
  output_file = "beech_lfmm_tmean_and_prec_to_check_ld_121224.txt"
)

# -----------------------------------------------------------------------------
# 2. Write LFMM top-1% tables used for manual LD marking
# -----------------------------------------------------------------------------

make_top1_lfmm_table(
  tmean = oak_lfmm_tmean,
  prec = oak_lfmm_prec,
  output_file = "oak_lfmm_top_1_percent_tmean_and_prec_121224.csv"
)

make_top1_lfmm_table(
  tmean = beech_lfmm_tmean,
  prec = beech_lfmm_prec,
  output_file = "beech_lfmm_top_1_percent_tmean_and_prec_121224.csv"
)

# -----------------------------------------------------------------------------
# 3. Prepare combined LFMM + GWAS position files for vcftools
#
# These inputs are expected to be manually LD-marked versions of the files above:
#   oak_lfmm_top_1_percent_tmean_and_prec_121224_ld_marked.csv
#   beech_lfmm_top_1_percent_tmean_and_prec_121224_ld_marked.csv
# -----------------------------------------------------------------------------

final_oak_lfmm <- read_ld_marked_lfmm("oak_lfmm_top_1_percent_tmean_and_prec_121224_ld_marked.csv")
final_beech_lfmm <- read_ld_marked_lfmm("beech_lfmm_top_1_percent_tmean_and_prec_121224_ld_marked.csv")

oak_snps <- make_oak_gwas_positions()
beech_snps_gwas <- make_beech_gwas_positions()

lfmm_and_gwas_oak <- prepare_combined_oak_positions(
  final_oak_lfmm = final_oak_lfmm,
  oak_snps = oak_snps,
  add_origin = FALSE
)

lfmm_and_gwas_beech <- prepare_combined_beech_positions(
  final_beech_lfmm = final_beech_lfmm,
  beech_snps_gwas = beech_snps_gwas,
  add_origin = FALSE
)

write_table_to(
  lfmm_and_gwas_oak,
  ld_check_dir,
  "lfmm_and_gwas_oak_131224_check_ld.txt"
)

write_table_to(
  lfmm_and_gwas_beech,
  ld_check_dir,
  "lfmm_and_gwas_beech_131224_check_ld.txt"
)

# -----------------------------------------------------------------------------
# 4. Summarize high-LD LFMM-GWAS pairs after vcftools has been run
# -----------------------------------------------------------------------------

ld_oak_unfiltered <- read.table(
  file.path(ld_check_dir, "lfmm_and_gwas_oak_131224_check_ld_test_filtered.geno.ld"),
  sep = "\t",
  header = TRUE
)

final_oak_lfmm <- read_ld_marked_lfmm("oak_lfmm_top_1_percent_tmean_and_prec_121224_ld_marked.csv")
oak_snps <- make_oak_gwas_positions()

lfmm_and_gwas_oak <- prepare_combined_oak_positions(
  final_oak_lfmm = final_oak_lfmm,
  oak_snps = oak_snps,
  add_origin = TRUE
)

oak_gwas_lfmm_high_ld_hits <- summarize_high_ld_pairs(
  ld_table = ld_oak_unfiltered,
  snp_table = lfmm_and_gwas_oak
)

ld_beech_unfiltered <- read.table(
  file.path(ld_check_dir, "lfmm_and_gwas_beech_131224_check_ld_test_filtered.geno.ld"),
  sep = "\t",
  header = TRUE
)

final_beech_lfmm <- read_ld_marked_lfmm("beech_lfmm_top_1_percent_tmean_and_prec_121224_ld_marked.csv")
beech_snps_gwas <- make_beech_gwas_positions()

lfmm_and_gwas_beech <- prepare_combined_beech_positions(
  final_beech_lfmm = final_beech_lfmm,
  beech_snps_gwas = beech_snps_gwas,
  add_origin = TRUE
)

beech_gwas_lfmm_high_ld_hits <- summarize_high_ld_pairs(
  ld_table = ld_beech_unfiltered,
  snp_table = lfmm_and_gwas_beech
)

write.csv(
  oak_gwas_lfmm_high_ld_hits,
  file = file.path(ld_check_dir, "oak_gwas_lfmm_high_ld_hits_131224.csv"),
  row.names = FALSE
)

write.csv(
  beech_gwas_lfmm_high_ld_hits,
  file = file.path(ld_check_dir, "beech_gwas_lfmm_high_ld_hits_131224.csv"),
  row.names = FALSE
)
