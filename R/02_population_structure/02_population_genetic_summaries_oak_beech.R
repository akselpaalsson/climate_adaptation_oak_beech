
################################################################################

rm(list = ls())


library(adegenet)
library(ade4)
library(vcfR)
library(hierfstat)
library(ggplot2)
library(MASS)
library(mefa)
library(ecodist)
library(dplyr)
library(tibble)
library(LEA)
library(viridis)
library(otuSummary)
library(stringr)
################################################################################
# 1. File paths
################################################################################

base_dir <- "path/to/base_dir/"

diversity_dir <- file.path(base_dir, "diversity_analyses")
snmf_script_dir <- file.path(base_dir, "scrp")
snmf_data_dir <- file.path(snmf_script_dir, "../data")
snmf_results_dir <- file.path(snmf_script_dir, "../res")

sample_info_file <- file.path(diversity_dir, "Common-Ring.sample.info.txt")
vcf_file <- file.path(diversity_dir, "populations_only_petrea.snps.vcf.gz.recode.vcf")
vcf_sample_file <- file.path(diversity_dir, "samples_populations_only_petrea.snps.vcf.gz.txt")

nei_fst_file <- file.path(diversity_dir, "nei_fst_oak_pops_only_petrea_050424.txt")


################################################################################
# 2. Import sample metadata and VCF
################################################################################

setwd(diversity_dir)

sample.info <- read.table(sample_info_file,sep = "\t",header = TRUE)

dat_vcf <- read.vcfR(vcf_file)
dat_genind <- vcfR2genind(dat_vcf)

sample.vcf <- read.table(vcf_sample_file,sep = "\t",header = FALSE)

colnames(sample.vcf)[1] <- "Sample"

sample.info.reduced <- merge(sample.info,sample.vcf,by.x = "SampleID",by.y = "Sample",sort = FALSE)

dat_genind$pop <- as.factor(sample.info.reduced$Prov2)

cols <- hcl.colors(length(levels(pop(dat_genind))),
palette = "viridis",alpha = 0.8,rev = FALSE,fixup = TRUE)


################################################################################
# 3. Basic population genetic statistics
################################################################################

dat_summary <- summary(dat_genind)

basic_stats_overall <- basic.stats(dat_genind)

inds_per_pop <- as.data.frame(dat_summary$n.by.pop)
inds_per_pop <- rownames_to_column(inds_per_pop, var = "pop")
colnames(inds_per_pop)[2] <- "nr"

populations <- unique(dat_genind$pop)

subsets <- lapply(populations, function(pop_i) {dat_genind[dat_genind$pop == pop_i, ]})

names(subsets) <- as.character(populations)

basic_stats_by_pop <- lapply(names(subsets), function(pop_i) {
  out <- as.data.frame(t(basic.stats(subsets[[pop_i]])$overall))
  out$pop <- pop_i
  out
})

basic_stats_overall_1 <- as.data.frame(t(basic_stats_overall$overall))
basic_stats_overall_1$pop <- "overall"

basic_pop_gen_stats <- do.call(rbind,
  c(basic_stats_by_pop, list(basic_stats_overall_1)))

basic_pop_gen_stats <- merge(inds_per_pop,
basic_pop_gen_stats,by = "pop",all = TRUE)


################################################################################
# 4. PCA
################################################################################

X <- scaleGen(dat_genind,center = TRUE,
scale = FALSE,NA.method = "mean")

pca1 <- dudi.pca(X,center = TRUE,
  scale = FALSE,scannf = FALSE,nf = 50)

summary(pca1)

cols <- hcl.colors(length(levels(pop(dat_genind))),
  palette = "viridis",alpha = 0.8,rev = FALSE,fixup = TRUE
)

s.class(
  pca1$li,
  pop(dat_genind),
  xax = 1,
  yax = 2,
  sub = "PCA 1-2",
  col = cols,
  clabel = 1
)
add.scatter.eig(
  pca1$eig[1:20],
  posi = "bottomright",
  nf = 3,
  xax = 1,
  yax = 2
)

s.class(
  pca1$li,
  pop(dat_genind),
  xax = 1,
  yax = 3,
  sub = "PCA 1-3",
  col = cols,
  clabel = 1
)
add.scatter.eig(
  pca1$eig[1:20],
  posi = "bottomright",
  nf = 3,
  xax = 1,
  yax = 3
)

s.class(
  pca1$li,
  pop(dat_genind),
  xax = 2,
  yax = 3,
  sub = "PCA 2-3",
  col = cols,
  clabel = 1
)
add.scatter.eig(
  pca1$eig[1:20],
  posi = "bottomright",
  nf = 3,
  xax = 2,
  yax = 3
)

barplot(
  pca1$eig[1:50],
  main = "PCA eigenvalues",
  col = heat.colors(50)
)


################################################################################
# 5. Pairwise Nei's FST
################################################################################

cats.matFst <- pairwise.neifst(dat_genind)

melted_data <- reshape2::melt(cats.matFst)

melted_data$Var1 <- as.character(melted_data$Var1)
melted_data$Var2 <- as.character(melted_data$Var2)

ggplot(melted_data, aes(Var1, Var2, fill = value)) +
  geom_tile(color = "white") +
  scale_fill_gradient(low = "white", high = "red") +
  theme_minimal() +
  labs(
    x = "Populations",
    y = "Populations",
    title = "Pairwise Fst Heatmap"
  )

colnames(melted_data) <- c("pop1", "pop2", "pwfst")

fst_clean <- na.omit(melted_data)

fit2 <- fitdistr(fst_clean$pwfst, densfun = "normal")
confint.default(fit2)


################################################################################
# 6. Isolation by distance
################################################################################

cats.matFst <- read.table(
  nei_fst_file,
  sep = ""
)

tablefst <- matrixConvert(cats.matFst)

subset_first_entries <- sample.info.reduced[
  !duplicated(sample.info.reduced$Prov2),
]

rownames(subset_first_entries) <- subset_first_entries$Prov

coords <- subset_first_entries[, c("Lon", "Lat")]

coords_matrix <- as.matrix(coords)

DM <- dist(coords_matrix)
distance_matrix <- as.matrix(DM)

rownames(distance_matrix) <- rownames(coords)

distance_df <- as.data.frame(distance_matrix)
table_distance_df <- matrixConvert(distance_df)

tablefst$sp2 <- gsub("X", "", tablefst$sp2)

tablefst$merging_col <- apply(
  tablefst[, c("sp1", "sp2")],
  1,
  function(x) paste(sort(x), collapse = "-")
)

table_distance_df$merging_col <- apply(
  table_distance_df[, c("sp1", "sp2")],
  1,
  function(x) paste(sort(x), collapse = "-")
)

tablefst <- tablefst %>%
  mutate(merging_col = trimws(merging_col))

ibd_data <- merge(
  tablefst,
  table_distance_df,
  by = "merging_col"
)

ibd_data$fst_clean <- ibd_data$dist.x / (1 - ibd_data$dist.x)

colnames(ibd_data)[7] <- "dist"
colnames(ibd_data)[8] <- "fst"

mantel_full <- ecodist::mantel(
  fst ~ dist,
  data = ibd_data,
  nperm = 1001
)

mantel_full

r_full <- mantel_full[1]
p1_full <- mantel_full[2]
p2_full <- mantel_full[3]
p3_full <- mantel_full[4]

plot_new <- ggplot(ibd_data, aes(x = dist, y = fst)) +
  geom_point() +
  labs(
    x = "Distance [°]",
    y = "FST [-]"
  ) +
  geom_smooth(
    method = "lm",
    se = FALSE,
    linetype = "dashed",
    color = "black"
  ) +
  annotate(
    "text",
    x = 18,
    y = 0.016,
    label = paste(
      "Mantel r =",
      round(r_full, 3),
      ", p =",
      round(p3_full, 3)
    )
  ) +
  theme_void() +
  theme(
    axis.line = element_line(color = "black"),
    axis.text.x = element_text(
      size = 12,
      color = "black",
      margin = margin(t = 5)
    ),
    axis.text.y = element_text(
      size = 12,
      color = "black",
      margin = margin(r = 5)
    ),
    axis.title = element_text(
      size = 14,
      color = "black"
    )
  ) +
  ggtitle("IBD")

plot_new


################################################################################
# 7. Import nucleotide diversity results from VCFtools
################################################################################

pi_files <- c(
  "4" = "populations_only_petrea_prov4_pi.windowed.pi",
  "26" = "populations_only_petrea_prov26_pi.windowed.pi",
  "31" = "populations_only_petrea_prov31_pi.windowed.pi",
  "33" = "populations_only_petrea_prov33_pi.windowed.pi",
  "50" = "populations_only_petrea_prov50_pi.windowed.pi",
  "77" = "populations_only_petrea_prov77_pi.windowed.pi",
  "87" = "populations_only_petrea_prov87_pi.windowed.pi",
  "260" = "populations_only_petrea_prov260_pi.windowed.pi",
  "263" = "populations_only_petrea_prov263_pi.windowed.pi"
)

pi_df <- do.call(
  rbind,
  lapply(names(pi_files), function(pop_i) {
    pi_data <- read.table(
      file.path(diversity_dir, pi_files[[pop_i]]),
      sep = "",
      header = TRUE
    )
    
    data.frame(
      pop = as.numeric(pop_i),
      pi = as.numeric(summary(pi_data$PI)[4])
    )
  })
)

baisc_pop_gen_stats_pi_added <- merge(
  basic_pop_gen_stats,
  pi_df,
  by = "pop",
  all = TRUE
)


################################################################################
# 8. SNMF analysis
################################################################################

setwd(snmf_results_dir)

dat.path <- "../data/"
dir.path <- "../res/"

ped_file <- paste0(
  dat.path,
  "populations_only_petrea.snps.vcf.gz.recode.vcf.snps.ped"
)

lfmm_file <- paste0(
  dat.path,
  "populations_only_petrea.snps.vcf.gz.recode.vcf.snps.ped.lfmm"
)

snmf_project_file <- paste0(
  dat.path,
  "populations_only_petrea.snps.vcf.gz.recode.vcf.snps.ped.snmfProject"
)

run_snmf <- FALSE

if (run_snmf) {
  LEA::ped2lfmm(
    input.file = ped_file,
    output.file = lfmm_file,
    force = TRUE
  )
  
  project.missing <- snmf(
    input.file = lfmm_file,
    K = 1:10,
    entropy = TRUE,
    repetitions = 10,
    CPU = 8,
    project = "new"
  )
}

proj.snmf <- load.snmfProject(snmf_project_file)

dir.create(
  file.path(snmf_results_dir, "snmf"),
  showWarnings = FALSE,
  recursive = TRUE
)

pdf(
  file.path(snmf_results_dir, "snmf", "SNMF_CrossEntropy_only_petrea.pdf"),
  width = 8.24 * 0.5,
  height = 8.24 * 0.5
)

plot(
  proj.snmf,
  cex = 1.2,
  col = "blue",
  pch = 19
)

dev.off()

K <- 5

sample.info <- read.table(
  file.path(snmf_data_dir, "Common-Ring.sample.info.txt"),
  sep = "\t",
  header = TRUE
)

sample.vcf <- read.table(
  vcf_sample_file,
  sep = "\t",
  header = FALSE
)

colnames(sample.vcf)[1] <- "Sample"

sample.info.reduced <- merge(
  sample.info,
  sample.vcf,
  by.x = "SampleID",
  by.y = "Sample",
  sort = FALSE
)

setwd(snmf_results_dir)

gc()

pdf(
  file.path(snmf_results_dir, "snmf", "SNMF_Qmatrix_K2_to_K10_only_petrea.pdf"),
  width = 10,
  height = 12
)

par(mfrow = c(3, 1))

for (clust in 2:10) {
  
  best <- which.min(cross.entropy(proj.snmf, K = clust))
  
  q_matrix_file <- paste0(
    dat.path,
    "populations_only_petrea.snps.vcf.gz.recode.vcf.snps.ped.snmf/K",
    clust,
    "/run",
    best,
    "/populations_only_petrea.snps.vcf.gz.recode.vcf.snps.ped_r",
    best,
    ".",
    clust,
    ".Q"
  )
  
  ap_best <- read.table(q_matrix_file)
  
  ap <- cbind(
    sample.info.reduced[1:7],
    ap_best
  )
  
  ap <- ap[order(ap$Prov2, ap$Site), ]
  
  write.table(
    ap,
    file.path(
      snmf_results_dir,
      "snmf",
      paste0("ap_DAPC_K_only_petrea", clust, ".txt")
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  
  ind.nb <- as.data.frame(table(ap$Prov2))
  colnames(ind.nb) <- c("Population", "Individuals")
  
  axis.position <- c()
  
  ap_1 <- as.matrix(t(ap))
  
  barplot(
    ap_1[8:(8 + clust - 1), ],
    main = paste0("Posterior assignment probailities K=", clust),
    col = viridis(clust),
    space = 0,
    border = NA,
    xaxt = "n",
    las = 3
  )
  
  for (j in seq_len(NROW(ind.nb))) {
    axis.position[j + 1] <- sum(axis.position[j], ind.nb[j, "Individuals"])
    
    abline(
      v = axis.position,
      col = "black",
      lty = 1,
      lwd = 1
    )
    
    text.position <- axis.position[j + 1] - ind.nb[j, "Individuals"] / 2
    
    text(
      text.position,
      -0.2,
      ind.nb[j, "Population"],
      srt = 90,
      cex = 1.0,
      xpd = TRUE,
      las = 1
    )
  }
}

dev.off()








###beech pop gen
rm(list = ls())
gc()


################################################################################
# 1. File paths
################################################################################

base_dir <- "path/to/base_dir"

analysis_dir <- file.path(base_dir, "pop_gen_beech")
results_dir <- file.path(base_dir, "results")
data_dir <- file.path(base_dir, "data")

vcf_prefix <- "populations_fagus_no_high_miss_ind_and_incorr_removed_HWE0_0001_minDP3_meanDP10_maxDP80_NA0.90_reduced_thin1kb.MAF0.05"

vcf_file <- file.path(
  analysis_dir,
  paste0(vcf_prefix, ".recode.vcf")
)

ped_file <- file.path(
  data_dir,
  paste0(vcf_prefix, ".snps.ped")
)

lfmm_file <- file.path(
  data_dir,
  paste0(
    "populations_fagus_no_high_miss_ind_HWE0_0001_minDP3_meanDP10_maxDP80_NA0.90_reduced_thin1kb.MAF0.05.snps.ped.lfmm"
  )
)

snmf_project_file <- file.path(
  data_dir,
  paste0(vcf_prefix, ".snps.snmfProject")
)

# Sample IDs are read directly from the VCF below.
# This avoids depending on a separate *.txt file that may not exist.

metadata_file <- file.path(data_dir, "beech_dbh.csv")
coords_file <- file.path(analysis_dir, "beech_prov_coords.csv")

snmf_dir <- file.path(results_dir, "snmf")

dir.create(snmf_dir, showWarnings = FALSE, recursive = TRUE)


################################################################################
# 2. Helper functions
################################################################################
read_beech_metadata <- function(sample_ids, metadata_file) {
  dbh_metadata <- read.csv(metadata_file, sep = ",", header = TRUE)
  
  new_sample.info <- dbh_metadata[, c("Treeid", "Site", "Block", "Prov", "Tree")]
  
  sample.info <- data.frame(
    SampleID = sample_ids,
    stringsAsFactors = FALSE
  )
  
  sample.info$RowNumber <- seq_len(nrow(sample.info))
  sample.info$RowNumber <- as.numeric(sample.info$RowNumber)
  
  split_ids <- strsplit(as.character(sample.info$SampleID), "-")
  sample.info <- cbind(sample.info, do.call(rbind, split_ids))
  
  colnames(sample.info)[1] <- "SampleID_vcf"
  colnames(sample.info)[2] <- "RowNumber"
  colnames(sample.info)[3] <- "sample_nr_vcf"
  colnames(sample.info)[4] <- "Treeid"
  
  sample.info <- sample.info[, c(
    "SampleID_vcf",
    "RowNumber",
    "sample_nr_vcf",
    "Treeid"
  )]
  
  sample_info_merged <- merge(
    sample.info,
    new_sample.info,
    by = "Treeid"
  )
  
  sample_info_merged[order(sample_info_merged$RowNumber), ]
}

standardize_comparison <- function(col) {
  sapply(col, function(x) {
    populations <- sort(unlist(strsplit(x, "-")))
    paste(populations, collapse = "-")
  })
}

replace_merging_col <- function(merging_col, replacements) {
  parts <- str_split(merging_col, "-")[[1]]
  replaced_parts <- sapply(parts, function(x) replacements[x])
  paste(replaced_parts, collapse = "-")
}


################################################################################
# 3. Import sample metadata and VCF
################################################################################

setwd(analysis_dir)

dat_vcf <- read.vcfR(vcf_file)

vcf_sample_ids <- colnames(dat_vcf@gt)[-1]

sample.info.reduced <- read_beech_metadata(
  sample_ids = vcf_sample_ids,
  metadata_file = metadata_file
)

dat_genind <- vcfR2genind(dat_vcf)
dat_genind$pop <- as.factor(sample.info.reduced$Prov)

cols <- hcl.colors(
  length(levels(pop(dat_genind))),
  palette = "viridis",
  alpha = 0.8,
  rev = FALSE,
  fixup = TRUE
)


################################################################################
# 4. Basic population genetic statistics
################################################################################

dat_summary <- summary(dat_genind)

basic_stats_overall <- basic.stats(dat_genind)

inds_per_pop <- as.data.frame(dat_summary$n.by.pop)
inds_per_pop <- rownames_to_column(inds_per_pop, var = "pop")
colnames(inds_per_pop)[2] <- "nr"

populations <- unique(dat_genind$pop)

subsets <- lapply(populations, function(pop_i) {
  dat_genind[dat_genind$pop == pop_i, ]
})

names(subsets) <- as.character(populations)

basic_stats_by_pop <- lapply(names(subsets), function(pop_i) {
  out <- as.data.frame(t(basic.stats(subsets[[pop_i]])$overall))
  out$pop <- pop_i
  out
})

basic_stats_overall_1 <- as.data.frame(t(basic_stats_overall$overall))
basic_stats_overall_1$pop <- "overall"

basic_pop_gen_stats <- do.call(
  rbind,
  c(basic_stats_by_pop, list(basic_stats_overall_1))
)

basic_pop_gen_stats <- merge(
  inds_per_pop,
  basic_pop_gen_stats,
  by = "pop",
  all = TRUE
)

write.csv(
  basic_pop_gen_stats,
  file = file.path(analysis_dir, "beech_basic_pop_gen_stats_250624.csv")
)

fst_results <- wc(dat_genind)
overall_fst <- fst_results$FST
overall_fst

boot_fst <- boot.ppfst(dat_genind, nboot = 1000)

ci <- boot_fst$ci
lower_limits <- boot_fst$ll
upper_limits <- boot_fst$ul

overall_lower_ci <- mean(lower_limits, na.rm = TRUE)
overall_upper_ci <- mean(upper_limits, na.rm = TRUE)


################################################################################
# 5. PCA
################################################################################

setwd(analysis_dir)

X <- scaleGen(
  dat_genind,
  center = TRUE,
  scale = FALSE,
  NA.method = "mean"
)

pca1 <- dudi.pca(
  X,
  center = TRUE,
  scale = FALSE,
  scannf = FALSE,
  nf = 50
)

summary(pca1)

cols <- hcl.colors(
  length(levels(pop(dat_genind))),
  palette = "viridis",
  alpha = 0.8,
  rev = FALSE,
  fixup = TRUE
)

s.class(
  pca1$li,
  pop(dat_genind),
  xax = 1,
  yax = 2,
  sub = "PCA 1-2",
  col = cols,
  clabel = 1
)
add.scatter.eig(
  pca1$eig[1:20],
  posi = "bottomright",
  nf = 3,
  xax = 1,
  yax = 2
)

s.class(
  pca1$li,
  pop(dat_genind),
  xax = 1,
  yax = 3,
  sub = "PCA 1-3",
  col = cols,
  clabel = 1
)
add.scatter.eig(
  pca1$eig[1:20],
  posi = "bottomright",
  nf = 3,
  xax = 1,
  yax = 3
)

s.class(
  pca1$li,
  pop(dat_genind),
  xax = 2,
  yax = 3,
  sub = "PCA 2-3",
  col = cols,
  clabel = 1
)
add.scatter.eig(
  pca1$eig[1:20],
  posi = "bottomright",
  nf = 3,
  xax = 2,
  yax = 3
)

barplot(
  pca1$eig[1:50],
  main = "PCA eigenvalues",
  col = heat.colors(50)
)


################################################################################
# 6. Pairwise Nei's FST
################################################################################

cats.matFst <- pairwise.neifst(dat_genind)


melted_data <- reshape2::melt(cats.matFst)

melted_data$Var1 <- as.character(melted_data$Var1)
melted_data$Var2 <- as.character(melted_data$Var2)

ggplot(melted_data, aes(Var1, Var2, fill = value)) +
  geom_tile(color = "white") +
  scale_fill_gradient(low = "white", high = "red") +
  theme_minimal() +
  labs(
    x = "Populations",
    y = "Populations",
    title = "Pairwise Fst Heatmap"
  )

colnames(melted_data) <- c("pop1", "pop2", "pwfst")

fst_clean <- na.omit(melted_data)

fit2 <- fitdistr(fst_clean$pwfst, densfun = "normal")
confint.default(fit2)


################################################################################
# 7. Isolation by distance
################################################################################

tablefst <- matrixConvert(cats.matFst)

beech_coords <- read.csv(coords_file)

coords <- beech_coords[, c("Lon", "Lat")]

DM <- dist(coords)
distance <- c(DM)
logD <- log(distance)

coords_matrix <- as.matrix(coords)

DM <- dist(coords_matrix)

distance_matrix <- as.matrix(DM)
rownames(distance_matrix) <- rownames(coords)

distance_df <- as.data.frame(distance_matrix)
table_distance_df <- matrixConvert(distance_df)

tablefst$sp2 <- gsub("X", "", tablefst$sp2)

tablefst$merging_col <- apply(
  tablefst[, c("sp1", "sp2")],
  1,
  function(x) paste(sort(x), collapse = "-")
)

table_distance_df$merging_col <- apply(
  table_distance_df[, c("sp1", "sp2")],
  1,
  function(x) paste(sort(x), collapse = "-")
)

tablefst <- tablefst %>%
  mutate(merging_col = trimws(merging_col))

replacements <- c(
  "1" = "FR002",
  "2" = "FR004",
  "3" = "FR008",
  "4" = "BE013",
  "5" = "NL014",
  "6" = "GB018",
  "7" = "SE023",
  "8" = "DE026",
  "9" = "DE031",
  "10" = "CZ048",
  "11" = "BG016"
)

table_distance_df <- table_distance_df %>%
  rowwise() %>%
  mutate(merging_col = replace_merging_col(merging_col, replacements)) %>%
  ungroup()

table_distance_df <- as.data.frame(table_distance_df)

table_distance_df2 <- table_distance_df %>%
  mutate(
    sp1 = replacements[as.character(sp1)],
    sp2 = replacements[as.character(sp2)]
  )

table_distance_df2$merging_col <- standardize_comparison(
  table_distance_df2$merging_col
)

tablefst$merging_col <- standardize_comparison(
  tablefst$merging_col
)

ibd_data <- merge(
  tablefst,
  table_distance_df2,
  by = "merging_col"
)

ibd_data$fst_clean <- ibd_data$dist.x / (1 - ibd_data$dist.x)

colnames(ibd_data)[7] <- "dist"
colnames(ibd_data)[8] <- "fst"

mantel_full <- ecodist::mantel(
  fst ~ dist,
  data = ibd_data,
  nperm = 1001
)

mantel_full

r_full <- mantel_full[1]
p1_full <- mantel_full[2]
p2_full <- mantel_full[3]
p3_full <- mantel_full[4]

plot_new <- ggplot(ibd_data, aes(x = dist, y = fst)) +
  geom_point() +
  labs(
    x = "Distance [°]",
    y = "FST [-]"
  ) +
  geom_smooth(
    method = "lm",
    se = FALSE,
    linetype = "dashed",
    color = "black"
  ) +
  annotate(
    "text",
    x = 18,
    y = 0.016,
    label = paste(
      "Mantel r =",
      round(r_full, 3),
      ", p =",
      round(p3_full, 3)
    )
  ) +
  theme_void() +
  theme(
    axis.line = element_line(color = "black"),
    axis.text.x = element_text(
      size = 12,
      color = "black",
      margin = margin(t = 5)
    ),
    axis.text.y = element_text(
      size = 12,
      color = "black",
      margin = margin(r = 5)
    ),
    axis.title = element_text(
      size = 14,
      color = "black"
    )
  ) +
  ggtitle("IBD")

plot_new


################################################################################
# 8. Import nucleotide diversity results from VCFtools
################################################################################

pi_files <- c(
  "FR002" = "populations_fagus_no_high_miss_ind_prov_fr2_pi.windowed.pi",
  "FR004" = "populations_fagus_no_high_miss_ind_prov_fr4_pi.windowed.pi",
  "FR008" = "populations_fagus_no_high_miss_ind_prov_fr8_pi.windowed.pi",
  "BE013" = "populations_fagus_no_high_miss_ind_prov_be13_pi.windowed.pi",
  "NL014" = "populations_fagus_no_high_miss_ind_prov_nl14_pi.windowed.pi",
  "GB018" = "populations_fagus_no_high_miss_ind_prov_gb18_pi.windowed.pi",
  "SE023" = "populations_fagus_no_high_miss_ind_prov_se23_pi.windowed.pi",
  "DE026" = "populations_fagus_no_high_miss_ind_prov_de26_pi.windowed.pi",
  "DE031" = "populations_fagus_no_high_miss_ind_prov_de31_pi.windowed.pi",
  "CZ048" = "populations_fagus_no_high_miss_ind_prov_cz48_pi.windowed.pi",
  "BG016" = "populations_fagus_no_high_miss_ind_prov_bg16_pi.windowed.pi",
  "overall" = "populations_fagus_no_high_miss_ind_overall_pi.windowed.pi"
)

pi_df <- do.call(
  rbind,
  lapply(names(pi_files), function(pop_i) {
    pi_data <- read.table(
      file.path(data_dir, pi_files[[pop_i]]),
      sep = "",
      header = TRUE
    )
    
    data.frame(
      pop = pop_i,
      pi = as.numeric(summary(pi_data$PI)[4])
    )
  })
)

baisc_pop_gen_stats_pi_added <- merge(
  basic_pop_gen_stats,
  pi_df,
  by = "pop",
  all = TRUE
)
################################################################################
# 9. SNMF analysis
################################################################################

setwd(results_dir)

run_snmf <- FALSE

if (run_snmf) {
  gc()
  
  LEA::ped2lfmm(
    input.file = ped_file,
    output.file = lfmm_file,
    force = TRUE
  )
  
  project.missing <- snmf(
    input.file = ped_file,
    K = 1:20,
    entropy = TRUE,
    repetitions = 10,
    CPU = 8,
    project = "new"
  )
}

proj.snmf <- load.snmfProject(snmf_project_file)

pdf(
  file.path(
    snmf_dir,
    paste0("SNMF_", vcf_prefix, ".pdf")
  ),
  width = 8.24 * 0.5,
  height = 8.24 * 0.5
)

plot(
  proj.snmf,
  cex = 1.2,
  col = "blue",
  pch = 19
)

dev.off()

K <- 17

gc()

pdf(
  file.path(snmf_dir, "SNMF_Qmatrix_K2_to_K20_fagus.pdf"),
  width = 10,
  height = 18
)

par(mfrow = c(3, 1))

for (clust in 2:20) {
  
  best <- which.min(cross.entropy(proj.snmf, K = clust))
  
  q_matrix_file <- file.path(
    data_dir,
    paste0(vcf_prefix, ".snps.snmf"),
    paste0("K", clust),
    paste0("run", best),
    paste0(vcf_prefix, ".snps_r", best, ".", clust, ".Q")
  )
  
  ap_best <- read.table(q_matrix_file)
  
  ap <- cbind(
    sample.info.reduced[1:7],
    ap_best
  )
  
  ap <- ap[order(ap$Prov, ap$Site), ]
  
  write.table(
    ap,
    file.path(snmf_dir, paste0("ap_DAPC_K", clust, ".txt")),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  
  ind.nb <- as.data.frame(table(ap$Prov))
  colnames(ind.nb) <- c("Population", "Individuals")
  
  axis.position <- c()
  
  ap_1 <- as.matrix(t(ap))
  
  barplot(
    ap_1[8:(8 + clust - 1), ],
    main = paste0("Posterior assignment probailities K=", clust),
    col = viridis(clust),
    space = 0,
    border = NA,
    xaxt = "n",
    las = 3
  )
  
  for (j in seq_len(NROW(ind.nb))) {
    axis.position[j + 1] <- sum(axis.position[j], ind.nb[j, "Individuals"])
    
    abline(
      v = axis.position,
      col = "black",
      lty = 1,
      lwd = 1
    )
    
    text.position <- axis.position[j + 1] - ind.nb[j, "Individuals"] / 2
    
    text(
      text.position,
      -0.2,
      ind.nb[j, "Population"],
      srt = 90,
      cex = 1.0,
      xpd = TRUE,
      las = 1
    )
  }
}

dev.off()


