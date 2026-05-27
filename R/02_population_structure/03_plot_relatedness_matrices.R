
setwd("/path/to/wd")

rel_mat_oak<-read.table("rel_mat_populations_only_petrea_dbh_height_df_for_rel_matrix_animal_model_minDP3_meanDP10_maxDP80_NA0.90_reduced_thin1kb.MAF0.05.recode.vcf.cXX.txt", sep="")
#install.packages("pheatmap")
library(pheatmap)
library(data.table)
# Convert the data.table to a matrix
str(rel_mat_oak)
class(rel_mat_oak)


#relatedness_matrix_oak <- fread(rel_mat_oak)
relatedness_matrix <- as.matrix(rel_mat_oak)

pheatmap(rel_mat_oak,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         main = "Genetic Relatedness Matrix Heatmap",
         color = colorRampPalette(c("white", "red"))(50))


pheatmap(rel_mat_oak,
         show_rownames = FALSE,
         show_colnames = FALSE,
         cluster_rows = TRUE,
         cluster_cols = TRUE)


pheatmap(rel_mat_oak,
         show_rownames = FALSE,
         fontsize_col = 4)



# threshold for half-sibs or closer
thr <- 0.25

# remove self-relatedness (diagonal)
diag(relatedness_matrix) <- NA

# logical matrix: TRUE if half-sib or closer
sib_mat <- relatedness_matrix >= thr

# number of sib / half-sib relationships per individual
n_sibs_per_ind <- rowSums(sib_mat, na.rm = TRUE)

n_sibs_per_ind

n_sib_pairs <- sum(sib_mat[upper.tri(sib_mat)], na.rm = TRUE)

setwd("path/to/wd")

rel_mat_beech<-read.table("beech_dbh_gwas_for_pca_chain2.recode.cXX.txt", sep="")
#install.packages("pheatmap")
library(pheatmap)
library(data.table)
# Convert the data.table to a matrix
str(rel_mat_beech)
class(rel_mat_beech)


relatedness_matrix <- fread(rel_mat_beech)
relatedness_matrix <- as.matrix(rel_mat_beech)

pheatmap(rel_mat_beech,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         main = "Genetic Relatedness Matrix Heatmap",
         color = colorRampPalette(c("white", "red"))(50))


pheatmap(rel_mat_oak,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         main = "Genetic Relatedness Matrix Heatmap",
         color = colorRampPalette(c("white", "red"))(50))


max_beech <- max(rel_mat_beech, na.rm = TRUE)
max_oak   <- max(rel_mat_oak, na.rm = TRUE)

min_beech <- min(rel_mat_beech, na.rm = TRUE)
min_oak   <- min(rel_mat_oak, na.rm = TRUE)

max_beech
max_oak


pheatmap(rel_mat_beech,
         show_rownames = FALSE,
         show_colnames = FALSE,
         cluster_rows = TRUE,
         cluster_cols = TRUE)


pheatmap(rel_mat_beech,
         show_rownames = FALSE,
         fontsize_col = 4)



# threshold for half-sibs or closer
thr <- 0.25

# remove self-relatedness (diagonal)
diag(relatedness_matrix) <- NA

# logical matrix: TRUE if half-sib or closer
sib_mat <- relatedness_matrix >= thr

# number of sib / half-sib relationships per individual
n_sibs_per_ind <- rowSums(sib_mat, na.rm = TRUE)

n_sibs_per_ind

n_sib_pairs <- sum(sib_mat[upper.tri(sib_mat)], na.rm = TRUE)



global_min <- min(min_beech, min_oak)
global_max <- max(max_beech, max_oak)
breaks <- seq(global_min, global_max, length.out = 51)
colors <- colorRampPalette(c("white", "red"))(50)


pheatmap(rel_mat_beech,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         main = "Genetic Relatedness Matrix – Beech",
         color = colors,
         breaks = breaks)

pheatmap(rel_mat_oak,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         main = "Genetic Relatedness Matrix – Oak",
         color = colors,
         breaks = breaks)

library(pheatmap)
library(gridExtra)

# global scale (as you already defined)
global_min <- min(min_beech, min_oak)
global_max <- max(max_beech, max_oak)

breaks <- seq(global_min, global_max, length.out = 51)
colors <- colorRampPalette(c("white", "red"))(50)

# Oak heatmap (TOP)
p_oak <- pheatmap(rel_mat_oak,
                  cluster_rows = TRUE,
                  cluster_cols = TRUE,
                  show_rownames = FALSE,
                  show_colnames = FALSE,
                  main = "Genetic Relatedness Matrix – Oak",
                  color = colors,
                  breaks = breaks,
                  legend = TRUE,
                  silent = TRUE)

# Beech heatmap (BOTTOM)
p_beech <- pheatmap(rel_mat_beech,
                    cluster_rows = TRUE,
                    cluster_cols = TRUE,
                    show_rownames = FALSE,
                    show_colnames = FALSE,
                    main = "Genetic Relatedness Matrix – Beech",
                    color = colors,
                    breaks = breaks,
                    legend = FALSE,   # avoid duplicate legend
                    silent = TRUE)

# Combine vertically (oak on top, beech bottom)
grid.arrange(p_oak$gtable,
             p_beech$gtable,
             ncol = 1,
             heights = c(1, 1))

############################################################################################################
############################################################################################################
############################################################################################################
############################################################################################################
############################################################################################################

## ============================================================
## Count pairwise relatedness > 0.25 (all comparisons)
## Works if rel_mat_* are data.frame / tibble / Matrix / etc.
## ============================================================

to_numeric_matrix <- function(x) {
  # convert common types to base numeric matrix
  if (inherits(x, "Matrix")) {
    x <- as.matrix(x)
  } else if (is.data.frame(x)) {
    x <- as.matrix(x)
  } else {
    x <- as.matrix(x)
  }
  
  # force numeric (handles characters/factors)
  suppressWarnings(storage.mode(x) <- "numeric")
  
  x
}

count_related_pairs_all <- function(K, threshold = 0.25) {
  K <- to_numeric_matrix(K)
  
  if (nrow(K) != ncol(K)) {
    stop("K must be square: nrow(K) must equal ncol(K).")
  }
  
  # remove diagonal
  diag(K) <- NA_real_
  
  # unique pairwise comparisons
  ut <- upper.tri(K)
  
  vals <- K[ut]
  vals <- vals[!is.na(vals)]
  
  list(
    threshold = threshold,
    total_pairs = length(vals),
    n_over = sum(vals > threshold),
    fraction_over = sum(vals > threshold) / length(vals),
    max_offdiag = max(vals, na.rm = TRUE)
  )
}

## =========================
## Run for both species
## =========================

res_beech <- count_related_pairs_all(rel_mat_beech, threshold = 0.25)
res_oak   <- count_related_pairs_all(rel_mat_oak,   threshold = 0.25)

## =========================
## Print results
## =========================

cat("\nBEECH\n")
cat("Pairs > 0.25:", res_beech$n_over, "out of", res_beech$total_pairs, "\n")
cat("Fraction:", round(res_beech$fraction_over, 6), "\n")
cat("Max off-diagonal:", round(res_beech$max_offdiag, 6), "\n")

cat("\nOAK\n")
cat("Pairs > 0.25:", res_oak$n_over, "out of", res_oak$total_pairs, "\n")
cat("Fraction:", round(res_oak$fraction_over, 6), "\n")
cat("Max off-diagonal:", round(res_oak$max_offdiag, 6), "\n")
