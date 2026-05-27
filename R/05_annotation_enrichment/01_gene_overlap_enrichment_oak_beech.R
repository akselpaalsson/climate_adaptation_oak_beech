
# -----------------------------------------------------------------------------

rm(list = ls())
gc()

library(dplyr)
library(VariantAnnotation)
library(GenomicRanges)
library(IRanges)
library(rtracklayer)

# -----------------------------------------------------------------------------
# Paths
# -----------------------------------------------------------------------------

setwd("/path/to/working_directory")

setwd(test_dir)

oak_work_dir <- test_dir
oak_gff_file <- file.path(test_dir, "GCF_932294415.1_dhQueRobu3.1_genomic.gff")
oak_vcf_file <- file.path(
  test_dir,
  "populations_no_no277_reps_missing_removed_only_petrea_minDP3_meanDP10_maxDP80_NA0.90_reduced.MAF0.05.recode.vcf"
)
oak_hit_file <- file.path(test_dir, "oak_all_bslmm_snps_table_151024.csv")

beech_work_dir <- test_dir
beech_gff_file <- file.path(test_dir, "Bhaga_genes.gff3")
beech_vcf_file <- file.path(
  test_dir,
  "populations_fagus_no_high_miss_ind_minDP3_meanDP10_maxDP80_NA0.90_reduced.MAF0.05.recode.vcf"
)
beech_hit_file <- file.path(test_dir, "snps_for_pleio_beech_101224_proper.csv")

output_dir <- test_dir

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

oak_scaffold_to_ncbi <- c(
  "OW028765.1" = "NC_065534.1",
  "OW028766.1" = "NC_065535.1",
  "OW028767.1" = "NC_065536.1",
  "OW028768.1" = "NC_065537.1",
  "OW028769.1" = "NC_065538.1",
  "OW028773.1" = "NC_065539.1",
  "OW028770.1" = "NC_065540.1",
  "OW028771.1" = "NC_065541.1",
  "OW028772.1" = "NC_065542.1",
  "OW028774.1" = "NC_065543.1",
  "OW028775.1" = "NC_065544.1",
  "OW028776.1" = "NC_065545.1"
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
  "12" = "Bhaga_12"
)

# -----------------------------------------------------------------------------
# Helpers
# -----------------------------------------------------------------------------

read_gene_annotations <- function(gff_file) {
  gff <- rtracklayer::import(gff_file)
  gff[gff$type %in% "gene"]
}

read_vcf_snps <- function(vcf_file, chromosome_map = NULL) {
  vcf <- readVcf(vcf_file)
  vcf_gr <- rowRanges(vcf)
  
  snps <- data.frame(
    chr = as.character(seqnames(vcf_gr)),
    pos = start(vcf_gr)
  )
  
  if (!is.null(chromosome_map)) {
    snps$chr <- chromosome_map[snps$chr]
  }
  
  snps <- snps[!is.na(snps$chr), ]
  
  snps
}

standardize_hit_coordinates <- function(hits) {
  # Some hit tables have CHR/BP columns, while others only have chr_bp.
  # This keeps the analysis unchanged but makes the input parsing robust.
  names(hits)[names(hits) == "chr"] <- "CHR"
  names(hits)[names(hits) == "pos"] <- "BP"
  
  if (!("CHR" %in% names(hits)) || !("BP" %in% names(hits))) {
    if (!("chr_bp" %in% names(hits))) {
      stop("Hit file must contain either CHR/BP columns or a chr_bp column.")
    }
    
    chr_bp_split <- strsplit(as.character(hits$chr_bp), "_")
    hits$CHR <- vapply(chr_bp_split, `[`, character(1), 1)
    hits$BP <- as.numeric(vapply(chr_bp_split, `[`, character(1), 2))
  }
  
  hits <- hits[!duplicated(paste(hits$CHR, hits$BP, sep = "_")), ]
  hits <- hits[, c("CHR", "BP")]
  colnames(hits) <- c("chr", "pos")
  
  hits
}

read_oak_hits <- function(hit_file) {
  hits <- read.csv(hit_file, check.names = FALSE)
  hits <- standardize_hit_coordinates(hits)
  
  hits$chr <- oak_numeric_to_scaffold[as.character(hits$chr)]
  hits$chr <- oak_scaffold_to_ncbi[hits$chr]
  
  hits <- hits[!is.na(hits$chr) & !is.na(hits$pos), ]
  
  hits
}

read_beech_hits <- function(hit_file) {
  hits <- read.csv(hit_file, sep = ";", check.names = FALSE)
  
  # If the file is actually comma-separated, reread it as comma-separated.
  if (ncol(hits) == 1 && grepl(",", names(hits)[1], fixed = TRUE)) {
    hits <- read.csv(hit_file, check.names = FALSE)
  }
  
  hits <- standardize_hit_coordinates(hits)
  
  hits <- subset(hits, chr != 14)
  hits$chr <- beech_numeric_to_scaffold[as.character(hits$chr)]
  
  hits <- hits[!is.na(hits$chr) & !is.na(hits$pos), ]
  
  hits
}

make_snp_ranges <- function(snps) {
  GRanges(
    seqnames = snps$chr,
    ranges = IRanges(start = snps$pos, width = 1)
  )
}

make_enrichment_result <- function(species, all_snps, gwas_hits, genes) {
  all_snps_gr <- make_snp_ranges(all_snps)
  gwas_hits_gr <- make_snp_ranges(gwas_hits)
  
  all_in_genes <- countOverlaps(all_snps_gr, genes) > 0
  hits_in_genes <- countOverlaps(gwas_hits_gr, genes) > 0
  
  table_data <- matrix(
    c(
      sum(hits_in_genes),
      sum(!hits_in_genes),
      sum(all_in_genes) - sum(hits_in_genes),
      sum(!all_in_genes) - sum(!hits_in_genes)
    ),
    nrow = 2,
    byrow = TRUE
  )
  
  rownames(table_data) <- c("GWAS_hits", "SNP_background_excluding_hits")
  colnames(table_data) <- c("in_genes", "not_in_genes")
  
  fisher_result <- fisher.test(table_data)
  
  summary_table <- data.frame(
    species = species,
    class = rownames(table_data),
    in_genes = table_data[, "in_genes"],
    not_in_genes = table_data[, "not_in_genes"],
    total = rowSums(table_data),
    percent_in_genes = table_data[, "in_genes"] / rowSums(table_data) * 100,
    fisher_p_value = fisher_result$p.value,
    odds_ratio = unname(fisher_result$estimate),
    fisher_conf_low = fisher_result$conf.int[1],
    fisher_conf_high = fisher_result$conf.int[2],
    row.names = NULL
  )
  
  list(
    table = summary_table,
    contingency_table = table_data,
    fisher = fisher_result
  )
}

write_enrichment_outputs <- function(result, table_file, fisher_file) {
  write.csv(
    result$table,
    file = file.path(output_dir, table_file),
    row.names = FALSE
  )
  
  sink(file.path(output_dir, fisher_file))
  print(result$contingency_table)
  cat("\n")
  print(result$fisher)
  sink()
}

# -----------------------------------------------------------------------------
# Oak enrichment test
# -----------------------------------------------------------------------------

oak_genes <- read_gene_annotations(oak_gff_file)
oak_all_snps <- read_vcf_snps(oak_vcf_file, chromosome_map = oak_scaffold_to_ncbi)
oak_hits <- read_oak_hits(oak_hit_file)

oak_result <- make_enrichment_result(
  species = "oak",
  all_snps = oak_all_snps,
  gwas_hits = oak_hits,
  genes = oak_genes
)

write_enrichment_outputs(
  oak_result,
  table_file = "oak_gene_overlap_enrichment_table.csv",
  fisher_file = "oak_gene_overlap_fisher_test.txt"
)

# -----------------------------------------------------------------------------
# Beech enrichment test
# -----------------------------------------------------------------------------

beech_genes <- read_gene_annotations(beech_gff_file)
beech_all_snps <- read_vcf_snps(beech_vcf_file)
beech_hits <- read_beech_hits(beech_hit_file)

beech_result <- make_enrichment_result(
  species = "beech",
  all_snps = beech_all_snps,
  gwas_hits = beech_hits,
  genes = beech_genes
)

write_enrichment_outputs(
  beech_result,
  table_file = "beech_gene_overlap_enrichment_table.csv",
  fisher_file = "beech_gene_overlap_fisher_test.txt"
)
