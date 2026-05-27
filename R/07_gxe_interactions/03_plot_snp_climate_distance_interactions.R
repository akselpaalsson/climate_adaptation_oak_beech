rm(list = ls())
gc()
library(dplyr)
library(tidyr)
library(ggplot2)
setwd("/path/to/working_directory")
bai_beech <- read.csv("dist_tmean_snp_prec_as_cov_101224_beech_snp_dist_STD_tmean_dist_prec_as_int_results_all_bai_041224_20260519.csv")
sla_beech <- read.csv("dist_tmean_snp_prec_as_cov_101224_beech_snp_dist_STD_tmean_dist_prec_as_int_results_all_sla_041224_20260519.csv")
height_beech <- read.csv("dist_tmean_snp_prec_as_cov_101224_beech_snp_dist_STD_tmean_dist_prec_as_int_results_all_height_041224_20260519.csv")
dbh_beech <- read.csv("dist_tmean_snp_prec_as_cov_101224_beech_snp_dist_STD_tmean_dist_prec_as_int_results_all_dbh_041224_20260519.csv")
rs_beech <- read.csv("dist_tmean_snp_prec_as_cov_101224_beech_snp_dist_STD_tmean_dist_prec_as_int_results_all_rs_041224_20260519.csv")
rc_beech <- read.csv("dist_tmean_snp_prec_as_cov_101224_beech_snp_dist_STD_tmean_dist_prec_as_int_results_all_rc_041224_20260519.csv")
rt_beech <- read.csv("dist_tmean_snp_prec_as_cov_101224_beech_snp_dist_STD_tmean_dist_prec_as_int_results_all_rt_041224_20260519.csv")
sprpre_beech <- read.csv("dist_tmean_snp_prec_as_cov_101224_beech_snp_dist_STD_tmean_dist_prec_as_int_results_all_prec_mean_spr_051224_20260519.csv")
sumpre_beech <- read.csv("dist_tmean_snp_prec_as_cov_101224_beech_snp_dist_STD_tmean_dist_prec_as_int_results_all_prec_mean_sum_051224_20260519.csv")
sumtmea_beech <- read.csv("dist_tmean_snp_prec_as_cov_101224_beech_snp_dist_STD_tmean_dist_prec_as_int_results_all_tmean_mean_sum_051224_20260519.csv")
sprtmea_beech <- read.csv("dist_tmean_snp_prec_as_cov_101224_beech_snp_dist_STD_tmean_dist_prec_as_int_results_all_tmean_mean_spr_051224_20260519.csv")
#
setwd("/path/to/working_directory")
bai_oak <- read.csv("dist_tmean_snp_prec_as_cov_101224_oak_std_preds_snp_dist_tmean_dist_prec_as_int_results_all_bai_051224.csv")
sla_oak <- read.csv("dist_tmean_snp_prec_as_cov_101224_oak_std_preds_snp_dist_tmean_dist_prec_as_int_results_all_sla_051224.csv")
height_oak <- read.csv("dist_tmean_snp_prec_as_cov_101224_oak_std_preds_snp_dist_tmean_dist_prec_as_int_results_all_height_051224.csv")
dbh_oak <- read.csv("dist_tmean_snp_prec_as_cov_101224_oak_std_preds_snp_dist_tmean_dist_prec_as_int_results_all_dbh_051224.csv")
rs_oak <- read.csv("dist_tmean_snp_prec_as_cov_101224_oak_std_preds_snp_dist_tmean_dist_prec_as_int_results_all_rs_051224.csv")
rc_oak <- read.csv("dist_tmean_snp_prec_as_cov_101224_oak_std_preds_snp_dist_tmean_dist_prec_as_int_results_all_rc_051224.csv")
rt_oak <- read.csv("dist_tmean_snp_prec_as_cov_101224_oak_std_preds_snp_dist_tmean_dist_prec_as_int_results_all_rt_051224.csv")
sprpre_oak <- read.csv("dist_tmean_snp_prec_as_cov_101224_oak_std_preds_snp_dist_tmean_dist_prec_as_int_results_all_sprpre_051224.csv")
sumpre_oak <- read.csv("dist_tmean_snp_prec_as_cov_101224_oak_std_preds_snp_dist_tmean_dist_prec_as_int_results_all_sumpre_051224.csv")
sumtmea_oak <- read.csv("dist_tmean_snp_prec_as_cov_101224_oak_std_preds_snp_dist_tmean_dist_prec_as_int_results_all_sumtmea_051224.csv")
sprtmea_oak <- read.csv("dist_tmean_snp_prec_as_cov_101224_oak_std_preds_snp_dist_tmean_dist_prec_as_int_results_all_sprtmea_051224.csv")


# Filter each DataFrame to include only unique snp_id entries
sumpre_oak <- distinct(sumpre_oak, snp_id, .keep_all = TRUE)
sumtmea_oak <- distinct(sumtmea_oak, snp_id, .keep_all = TRUE)
sprtmea_oak <- distinct(sprtmea_oak, snp_id, .keep_all = TRUE)
bai_oak <- distinct(bai_oak, snp_id, .keep_all = TRUE)
sla_oak <- distinct(sla_oak, snp_id, .keep_all = TRUE)
height_oak <- distinct(height_oak, snp_id, .keep_all = TRUE)
dbh_oak <- distinct(dbh_oak, snp_id, .keep_all = TRUE)
rs_oak <- distinct(rs_oak, snp_id, .keep_all = TRUE)
rc_oak <- distinct(rc_oak, snp_id, .keep_all = TRUE)
rt_oak <- distinct(rt_oak, snp_id, .keep_all = TRUE)
sprpre_oak <- distinct(sprpre_oak, snp_id, .keep_all = TRUE)
sumtmea_beech <- distinct(sumtmea_beech, snp_id, .keep_all = TRUE)
sprtmea_beech <- distinct(sprtmea_beech, snp_id, .keep_all = TRUE)
bai_beech <- distinct(bai_beech, snp_id, .keep_all = TRUE)
sla_beech <- distinct(sla_beech, snp_id, .keep_all = TRUE)
height_beech <- distinct(height_beech, snp_id, .keep_all = TRUE)
dbh_beech <- distinct(dbh_beech, snp_id, .keep_all = TRUE)
rs_beech <- distinct(rs_beech, snp_id, .keep_all = TRUE)
rc_beech <- distinct(rc_beech, snp_id, .keep_all = TRUE)
rt_beech <- distinct(rt_beech, snp_id, .keep_all = TRUE)
sprpre_beech <- distinct(sprpre_beech, snp_id, .keep_all = TRUE)

# Now, each DataFrame has only unique snp_id entries
sumpre_oak$species <- "oak";sumtmea_oak$species<- "oak";sprtmea_oak$species<- "oak";bai_oak$species<- "oak";sla_oak$species<- "oak"
height_oak$species<- "oak";dbh_oak$species<- "oak";rs_oak$species<- "oak";rc_oak$species<- "oak";rt_oak$species<- "oak"
sumpre_beech$species <- "beech";sumtmea_beech$species<- "beech";sprtmea_beech$species<- "beech";bai_beech$species<- "beech"
sla_beech$species<- "beech";height_beech$species<- "beech";dbh_beech$species<- "beech";rs_beech$species<- "beech"
rc_beech$species<- "beech";rt_beech$species<- "beech";sprpre_oak$species <- "oak";sprpre_beech$species <- "beech"

colnames(sumpre_beech)
colnames(sumpre_oak)
sumpre <- rbind(sumpre_oak, sumpre_beech)
sumtmea <- rbind(sumtmea_oak, sumtmea_beech)
bai <- rbind(bai_oak, bai_beech)
sla <- rbind(sla_oak, sla_beech)
height <- rbind(height_oak, height_beech)
dbh <-rbind(dbh_oak, dbh_beech)
rs <- rbind(rs_oak, rs_beech)
rt <- rbind(rt_oak, rt_beech)
rc <- rbind(rc_oak,rc_beech)
sprtmea <- rbind(sprtmea_oak, sprtmea_beech)
sprpre <- rbind(sprpre_oak, sprpre_beech)

invert_trait_effects <- function(df) {
  df$dist.trend_AA_tmean <- -df$dist.trend_AA_tmean
  df$dist.trend_AB_tmean <- -df$dist.trend_AB_tmean
  
  if (all(c("dist.trend_lower.CL_AA_tmean", "dist.trend_upper.CL_AA_tmean") %in% colnames(df))) {
    lower_aa <- df$dist.trend_lower.CL_AA_tmean
    upper_aa <- df$dist.trend_upper.CL_AA_tmean
    
    df$dist.trend_lower.CL_AA_tmean <- -upper_aa
    df$dist.trend_upper.CL_AA_tmean <- -lower_aa
  }
  
  if (all(c("dist.trend_lower.CL_AB_tmean", "dist.trend_upper.CL_AB_tmean") %in% colnames(df))) {
    lower_ab <- df$dist.trend_lower.CL_AB_tmean
    upper_ab <- df$dist.trend_upper.CL_AB_tmean
    
    df$dist.trend_lower.CL_AB_tmean <- -upper_ab
    df$dist.trend_upper.CL_AB_tmean <- -lower_ab
  }
  
  if ("contrast_tmean_estimate" %in% colnames(df)) {
    df$contrast_tmean_estimate <- -df$contrast_tmean_estimate
  }
  
  if (all(c("contrast_tmean_estimate_lower_CL", "contrast_tmean_estimate_upper_CL") %in% colnames(df))) {
    lower_contrast <- df$contrast_tmean_estimate_lower_CL
    upper_contrast <- df$contrast_tmean_estimate_upper_CL
    
    df$contrast_tmean_estimate_lower_CL <- -upper_contrast
    df$contrast_tmean_estimate_upper_CL <- -lower_contrast
  }
  
  if ("coefs_m1_interaction_snp_dist_tmean_estimate" %in% colnames(df)) {
    df$coefs_m1_interaction_snp_dist_tmean_estimate <- -df$coefs_m1_interaction_snp_dist_tmean_estimate
  }
  
  df
}

sumpre <- invert_trait_effects(sumpre)
sprpre <- invert_trait_effects(sprpre)
sumtmea <- invert_trait_effects(sumtmea)
sprtmea <- invert_trait_effects(sprtmea)

summary(sumpre$dist.trend_AA_tmean)
summary(sumpre$dist.trend_AB_tmean)
summary(sprpre$dist.trend_AA_tmean)
summary(sprpre$dist.trend_AB_tmean)
summary(sumtmea$dist.trend_AA_tmean)
summary(sumtmea$dist.trend_AB_tmean)
summary(sprtmea$dist.trend_AA_tmean)
summary(sprtmea$dist.trend_AB_tmean)
summary(bai$dist.trend_AA_tmean)
summary(bai$dist.trend_AB_tmean)
summary(dbh$dist.trend_AA_tmean)
summary(dbh$dist.trend_AB_tmean)
summary(sla$dist.trend_AA_tmean)
summary(sla$dist.trend_AB_tmean)
summary(height$dist.trend_AA_tmean)
summary(height$dist.trend_AB_tmean)
summary(rc$dist.trend_AA_tmean)
summary(rc$dist.trend_AB_tmean)
summary(rt$dist.trend_AA_tmean)
summary(rt$dist.trend_AB_tmean)
summary(rs$dist.trend_AA_tmean)
summary(rs$dist.trend_AB_tmean)



ggplot(dbh, aes(x = dist.trend_AB_tmean, y = dist.trend_AA_tmean)) +
  geom_point(aes(color = coefs_m1_interaction_snp_dist_tmean_p.value < 0.05,shape = factor(case_when(
    dist.trend_p.value_AA_tmean < 0.05 & dist.trend_p.value_AB_tmean < 0.05 & species == "oak" ~ "filled_oak",
    dist.trend_p.value_AA_tmean < 0.05 & dist.trend_p.value_AB_tmean < 0.05 & species == "beech" ~ "filled_beech",
    dist.trend_p.value_AA_tmean < 0.05 & dist.trend_p.value_AB_tmean >= 0.05 & species == "oak" ~ "hollow_oak",
    dist.trend_p.value_AA_tmean < 0.05 & dist.trend_p.value_AB_tmean >= 0.05 & species == "beech" ~ "hollow_beech",
    dist.trend_p.value_AA_tmean >= 0.05 & dist.trend_p.value_AB_tmean < 0.05 & species == "oak" ~ "cross_oak",
    dist.trend_p.value_AA_tmean >= 0.05 & dist.trend_p.value_AB_tmean < 0.05 & species == "beech" ~ "cross_beech",
    dist.trend_p.value_AA_tmean >= 0.05 & dist.trend_p.value_AB_tmean >= 0.05 & species == "beech" ~ "thick_outline_beech",
    dist.trend_p.value_AA_tmean >= 0.05 & dist.trend_p.value_AB_tmean >= 0.05 & species == "oak" ~ "thick_outline_oak")),
    fill = factor(case_when(dist.trend_p.value_AA_tmean < 0.05 & dist.trend_p.value_AB_tmean < 0.05 ~ "filled",
                            dist.trend_p.value_AA_tmean >= 0.05 & dist.trend_p.value_AB_tmean >= 0.05 ~ "not_filled",
                            TRUE ~ NA_character_)),size = 3, stroke = 1.5)) +
  geom_hline(yintercept = 0, linetype = "dashed") +geom_vline(xintercept = 0, linetype = "dashed") +
  xlim(-2.5, 2.5) +ylim(-2.5, 2.5) +labs(x = "beta_B", y = "beta_A", title = "dbh tmea") +
  scale_shape_manual(values = c("filled_oak" = 16, "filled_beech" = 18, "hollow_oak" = 1, "hollow_beech" = 5,
                                "cross_oak" = 10, "cross_beech" = 9, "thick_outline_oak" = 21, "thick_outline_beech" = 23)) +
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "blue")) +
  scale_fill_manual(values = c("filled" = "red", "not_filled" = "grey")) + # Set "not_filled" to grey
  theme_minimal() +theme(panel.grid.major = element_blank(),
                         panel.grid.minor = element_blank(),legend.position = "none",axis.text.x = element_text(size = 14),
                         axis.text.y = element_text(size = 14),panel.border = element_rect(color = "black", fill = NA, size = 1.5))


ggplot(height, aes(x = dist.trend_AB_tmean, y = dist.trend_AA_tmean)) +
  geom_point(aes(color = coefs_m1_interaction_snp_dist_tmean_p.value < 0.05,shape = factor(case_when(
    dist.trend_p.value_AA_tmean < 0.05 & dist.trend_p.value_AB_tmean < 0.05 & species == "oak" ~ "filled_oak",
    dist.trend_p.value_AA_tmean < 0.05 & dist.trend_p.value_AB_tmean < 0.05 & species == "beech" ~ "filled_beech",
    dist.trend_p.value_AA_tmean < 0.05 & dist.trend_p.value_AB_tmean >= 0.05 & species == "oak" ~ "hollow_oak",
    dist.trend_p.value_AA_tmean < 0.05 & dist.trend_p.value_AB_tmean >= 0.05 & species == "beech" ~ "hollow_beech",
    dist.trend_p.value_AA_tmean >= 0.05 & dist.trend_p.value_AB_tmean < 0.05 & species == "oak" ~ "cross_oak",
    dist.trend_p.value_AA_tmean >= 0.05 & dist.trend_p.value_AB_tmean < 0.05 & species == "beech" ~ "cross_beech",
    dist.trend_p.value_AA_tmean >= 0.05 & dist.trend_p.value_AB_tmean >= 0.05 & species == "beech" ~ "thick_outline_beech",
    dist.trend_p.value_AA_tmean >= 0.05 & dist.trend_p.value_AB_tmean >= 0.05 & species == "oak" ~ "thick_outline_oak")),
    fill = factor(case_when(dist.trend_p.value_AA_tmean < 0.05 & dist.trend_p.value_AB_tmean < 0.05 ~ "filled",
                            dist.trend_p.value_AA_tmean >= 0.05 & dist.trend_p.value_AB_tmean >= 0.05 ~ "not_filled",
                            TRUE ~ NA_character_)),size = 3, stroke = 1.5)) +
  geom_hline(yintercept = 0, linetype = "dashed") +geom_vline(xintercept = 0, linetype = "dashed") +
  xlim(-2.5, 2.5) +ylim(-2.5, 2.5) +labs(x = "beta_B", y = "beta_A", title = "height tmea") +
  scale_shape_manual(values = c("filled_oak" = 16, "filled_beech" = 18, "hollow_oak" = 1, "hollow_beech" = 5,
                                "cross_oak" = 10, "cross_beech" = 9, "thick_outline_oak" = 21, "thick_outline_beech" = 23)) +
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "blue")) +
  scale_fill_manual(values = c("filled" = "red", "not_filled" = "grey")) + # Set "not_filled" to grey
  theme_minimal() +theme(panel.grid.major = element_blank(),
                         panel.grid.minor = element_blank(),legend.position = "none",axis.text.x = element_text(size = 14),
                         axis.text.y = element_text(size = 14),panel.border = element_rect(color = "black", fill = NA, size = 1.5))




ggplot(bai, aes(x = dist.trend_AB_tmean, y = dist.trend_AA_tmean)) +
  geom_point(aes(color = coefs_m1_interaction_snp_dist_tmean_p.value < 0.05,shape = factor(case_when(
    dist.trend_p.value_AA_tmean < 0.05 & dist.trend_p.value_AB_tmean < 0.05 & species == "oak" ~ "filled_oak",
    dist.trend_p.value_AA_tmean < 0.05 & dist.trend_p.value_AB_tmean < 0.05 & species == "beech" ~ "filled_beech",
    dist.trend_p.value_AA_tmean < 0.05 & dist.trend_p.value_AB_tmean >= 0.05 & species == "oak" ~ "hollow_oak",
    dist.trend_p.value_AA_tmean < 0.05 & dist.trend_p.value_AB_tmean >= 0.05 & species == "beech" ~ "hollow_beech",
    dist.trend_p.value_AA_tmean >= 0.05 & dist.trend_p.value_AB_tmean < 0.05 & species == "oak" ~ "cross_oak",
    dist.trend_p.value_AA_tmean >= 0.05 & dist.trend_p.value_AB_tmean < 0.05 & species == "beech" ~ "cross_beech",
    dist.trend_p.value_AA_tmean >= 0.05 & dist.trend_p.value_AB_tmean >= 0.05 & species == "beech" ~ "thick_outline_beech",
    dist.trend_p.value_AA_tmean >= 0.05 & dist.trend_p.value_AB_tmean >= 0.05 & species == "oak" ~ "thick_outline_oak")),
    fill = factor(case_when(dist.trend_p.value_AA_tmean < 0.05 & dist.trend_p.value_AB_tmean < 0.05 ~ "filled",
                            dist.trend_p.value_AA_tmean >= 0.05 & dist.trend_p.value_AB_tmean >= 0.05 ~ "not_filled",
                            TRUE ~ NA_character_)),size = 3, stroke = 1.5)) +
  geom_hline(yintercept = 0, linetype = "dashed") +geom_vline(xintercept = 0, linetype = "dashed") +
  xlim(-2.5, 2.5) +ylim(-2.5, 2.5) +labs(x = "beta_B", y = "beta_A", title = "bai tmea") +
  scale_shape_manual(values = c("filled_oak" = 16, "filled_beech" = 18, "hollow_oak" = 1, "hollow_beech" = 5,
                                "cross_oak" = 10, "cross_beech" = 9, "thick_outline_oak" = 21, "thick_outline_beech" = 23)) +
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "blue")) +
  scale_fill_manual(values = c("filled" = "red", "not_filled" = "grey")) + # Set "not_filled" to grey
  theme_minimal() +theme(panel.grid.major = element_blank(),
                         panel.grid.minor = element_blank(),legend.position = "none",axis.text.x = element_text(size = 14),
                         axis.text.y = element_text(size = 14),panel.border = element_rect(color = "black", fill = NA, size = 1.5))



ggplot(sla, aes(x = dist.trend_AB_tmean, y = dist.trend_AA_tmean)) +
  geom_point(aes(color = coefs_m1_interaction_snp_dist_tmean_p.value < 0.05,shape = factor(case_when(
    dist.trend_p.value_AA_tmean < 0.05 & dist.trend_p.value_AB_tmean < 0.05 & species == "oak" ~ "filled_oak",
    dist.trend_p.value_AA_tmean < 0.05 & dist.trend_p.value_AB_tmean < 0.05 & species == "beech" ~ "filled_beech",
    dist.trend_p.value_AA_tmean < 0.05 & dist.trend_p.value_AB_tmean >= 0.05 & species == "oak" ~ "hollow_oak",
    dist.trend_p.value_AA_tmean < 0.05 & dist.trend_p.value_AB_tmean >= 0.05 & species == "beech" ~ "hollow_beech",
    dist.trend_p.value_AA_tmean >= 0.05 & dist.trend_p.value_AB_tmean < 0.05 & species == "oak" ~ "cross_oak",
    dist.trend_p.value_AA_tmean >= 0.05 & dist.trend_p.value_AB_tmean < 0.05 & species == "beech" ~ "cross_beech",
    dist.trend_p.value_AA_tmean >= 0.05 & dist.trend_p.value_AB_tmean >= 0.05 & species == "beech" ~ "thick_outline_beech",
    dist.trend_p.value_AA_tmean >= 0.05 & dist.trend_p.value_AB_tmean >= 0.05 & species == "oak" ~ "thick_outline_oak")),
    fill = factor(case_when(dist.trend_p.value_AA_tmean < 0.05 & dist.trend_p.value_AB_tmean < 0.05 ~ "filled",
                            dist.trend_p.value_AA_tmean >= 0.05 & dist.trend_p.value_AB_tmean >= 0.05 ~ "not_filled",
                            TRUE ~ NA_character_)),size = 3, stroke = 1.5)) +
  geom_hline(yintercept = 0, linetype = "dashed") +geom_vline(xintercept = 0, linetype = "dashed") +
  xlim(-2.5, 2.5) +ylim(-2.5, 2.5) +labs(x = "beta_B", y = "beta_A", title = "sla tmea") +
  scale_shape_manual(values = c("filled_oak" = 16, "filled_beech" = 18, "hollow_oak" = 1, "hollow_beech" = 5,
                                "cross_oak" = 10, "cross_beech" = 9, "thick_outline_oak" = 21, "thick_outline_beech" = 23)) +
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "blue")) +
  scale_fill_manual(values = c("filled" = "red", "not_filled" = "grey")) + # Set "not_filled" to grey
  theme_minimal() +theme(panel.grid.major = element_blank(),
                         panel.grid.minor = element_blank(),legend.position = "none",axis.text.x = element_text(size = 14),
                         axis.text.y = element_text(size = 14),panel.border = element_rect(color = "black", fill = NA, size = 1.5))




ggplot(rc, aes(x = dist.trend_AB_tmean, y = dist.trend_AA_tmean)) +
  geom_point(aes(color = coefs_m1_interaction_snp_dist_tmean_p.value < 0.05,shape = factor(case_when(
    dist.trend_p.value_AA_tmean < 0.05 & dist.trend_p.value_AB_tmean < 0.05 & species == "oak" ~ "filled_oak",
    dist.trend_p.value_AA_tmean < 0.05 & dist.trend_p.value_AB_tmean < 0.05 & species == "beech" ~ "filled_beech",
    dist.trend_p.value_AA_tmean < 0.05 & dist.trend_p.value_AB_tmean >= 0.05 & species == "oak" ~ "hollow_oak",
    dist.trend_p.value_AA_tmean < 0.05 & dist.trend_p.value_AB_tmean >= 0.05 & species == "beech" ~ "hollow_beech",
    dist.trend_p.value_AA_tmean >= 0.05 & dist.trend_p.value_AB_tmean < 0.05 & species == "oak" ~ "cross_oak",
    dist.trend_p.value_AA_tmean >= 0.05 & dist.trend_p.value_AB_tmean < 0.05 & species == "beech" ~ "cross_beech",
    dist.trend_p.value_AA_tmean >= 0.05 & dist.trend_p.value_AB_tmean >= 0.05 & species == "beech" ~ "thick_outline_beech",
    dist.trend_p.value_AA_tmean >= 0.05 & dist.trend_p.value_AB_tmean >= 0.05 & species == "oak" ~ "thick_outline_oak")),
    fill = factor(case_when(dist.trend_p.value_AA_tmean < 0.05 & dist.trend_p.value_AB_tmean < 0.05 ~ "filled",
                            dist.trend_p.value_AA_tmean >= 0.05 & dist.trend_p.value_AB_tmean >= 0.05 ~ "not_filled",
                            TRUE ~ NA_character_)),size = 3, stroke = 1.5)) +
  geom_hline(yintercept = 0, linetype = "dashed") +geom_vline(xintercept = 0, linetype = "dashed") +
  xlim(-2.7, 2.7) +ylim(-2.7, 2.7) +labs(x = "beta_B", y = "beta_A", title = "rc tmea") +
  scale_shape_manual(values = c("filled_oak" = 16, "filled_beech" = 18, "hollow_oak" = 1, "hollow_beech" = 5,
                                "cross_oak" = 10, "cross_beech" = 9, "thick_outline_oak" = 21, "thick_outline_beech" = 23)) +
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "blue")) +
  scale_fill_manual(values = c("filled" = "red", "not_filled" = "grey")) + # Set "not_filled" to grey
  theme_minimal() +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),legend.position = "none",
        axis.text.x = element_text(size = 14),axis.text.y = element_text(size = 14),
        axis.ticks = element_line(color = "black", linewidth = 0.8),axis.ticks.length = unit(4, "mm"),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1.5)
  )




ggplot(rt, aes(x = dist.trend_AB_tmean, y = dist.trend_AA_tmean)) +
  geom_point(aes(color = coefs_m1_interaction_snp_dist_tmean_p.value < 0.05,shape = factor(case_when(
    dist.trend_p.value_AA_tmean < 0.05 & dist.trend_p.value_AB_tmean < 0.05 & species == "oak" ~ "filled_oak",
    dist.trend_p.value_AA_tmean < 0.05 & dist.trend_p.value_AB_tmean < 0.05 & species == "beech" ~ "filled_beech",
    dist.trend_p.value_AA_tmean < 0.05 & dist.trend_p.value_AB_tmean >= 0.05 & species == "oak" ~ "hollow_oak",
    dist.trend_p.value_AA_tmean < 0.05 & dist.trend_p.value_AB_tmean >= 0.05 & species == "beech" ~ "hollow_beech",
    dist.trend_p.value_AA_tmean >= 0.05 & dist.trend_p.value_AB_tmean < 0.05 & species == "oak" ~ "cross_oak",
    dist.trend_p.value_AA_tmean >= 0.05 & dist.trend_p.value_AB_tmean < 0.05 & species == "beech" ~ "cross_beech",
    dist.trend_p.value_AA_tmean >= 0.05 & dist.trend_p.value_AB_tmean >= 0.05 & species == "beech" ~ "thick_outline_beech",
    dist.trend_p.value_AA_tmean >= 0.05 & dist.trend_p.value_AB_tmean >= 0.05 & species == "oak" ~ "thick_outline_oak")),
    fill = factor(case_when(dist.trend_p.value_AA_tmean < 0.05 & dist.trend_p.value_AB_tmean < 0.05 ~ "filled",
                            dist.trend_p.value_AA_tmean >= 0.05 & dist.trend_p.value_AB_tmean >= 0.05 ~ "not_filled",
                            TRUE ~ NA_character_)),size = 3, stroke = 1.5)) +
  geom_hline(yintercept = 0, linetype = "dashed") +geom_vline(xintercept = 0, linetype = "dashed") +
  xlim(-2.7, 2.7) +ylim(-2.7, 2.7) +labs(x = "beta_B", y = "beta_A", title = "rt tmea") +
  scale_shape_manual(values = c("filled_oak" = 16, "filled_beech" = 18, "hollow_oak" = 1, "hollow_beech" = 5,
                                "cross_oak" = 10, "cross_beech" = 9, "thick_outline_oak" = 21, "thick_outline_beech" = 23)) +
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "blue")) +
  scale_fill_manual(values = c("filled" = "red", "not_filled" = "grey")) + # Set "not_filled" to grey
  theme_minimal() +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),legend.position = "none",
        axis.text.x = element_text(size = 14),axis.text.y = element_text(size = 14),
        axis.ticks = element_line(color = "black", linewidth = 0.8),axis.ticks.length = unit(4, "mm"),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1.5)
  )


ggplot(rs, aes(x = dist.trend_AB_tmean, y = dist.trend_AA_tmean)) +
  geom_point(aes(
    color = coefs_m1_interaction_snp_dist_tmean_p.value < 0.05,
    shape = factor(case_when(
      dist.trend_p.value_AA_tmean < 0.05 & dist.trend_p.value_AB_tmean < 0.05 & species == "oak" ~ "filled_oak",
      dist.trend_p.value_AA_tmean < 0.05 & dist.trend_p.value_AB_tmean < 0.05 & species == "beech" ~ "filled_beech",
      dist.trend_p.value_AA_tmean < 0.05 & dist.trend_p.value_AB_tmean >= 0.05 & species == "oak" ~ "hollow_oak",
      dist.trend_p.value_AA_tmean < 0.05 & dist.trend_p.value_AB_tmean >= 0.05 & species == "beech" ~ "hollow_beech",
      dist.trend_p.value_AA_tmean >= 0.05 & dist.trend_p.value_AB_tmean < 0.05 & species == "oak" ~ "cross_oak",
      dist.trend_p.value_AA_tmean >= 0.05 & dist.trend_p.value_AB_tmean < 0.05 & species == "beech" ~ "cross_beech",
      dist.trend_p.value_AA_tmean >= 0.05 & dist.trend_p.value_AB_tmean >= 0.05 & species == "beech" ~ "thick_outline_beech",
      dist.trend_p.value_AA_tmean >= 0.05 & dist.trend_p.value_AB_tmean >= 0.05 & species == "oak" ~ "thick_outline_oak"
    )),
    fill = factor(case_when(
      dist.trend_p.value_AA_tmean < 0.05 & dist.trend_p.value_AB_tmean < 0.05 ~ "filled",
      dist.trend_p.value_AA_tmean >= 0.05 & dist.trend_p.value_AB_tmean >= 0.05 ~ "not_filled",
      TRUE ~ NA_character_
    ))
  ), size = 4.5, stroke = .7) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(x = "beta_B", y = "beta_A", title = "rs tmea square") +
  scale_shape_manual(values = c(
    "filled_oak" = 16, "filled_beech" = 18, "hollow_oak" = 1, "hollow_beech" = 5,
    "cross_oak" = 10, "cross_beech" = 9, "thick_outline_oak" = 21, "thick_outline_beech" = 23
  )) +
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "blue")) +
  scale_fill_manual(values = c("filled" = "red", "not_filled" = "grey")) +
  theme_minimal() +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),legend.position = "none",
        axis.text.x = element_text(size = 14),axis.text.y = element_text(size = 14),
        axis.ticks = element_line(color = "black", linewidth = 0.8),axis.ticks.length = unit(4, "mm"),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1.5)
  ) +
  coord_fixed(ratio = 1, xlim = c(-2.9, 2.9), ylim = c(-2.9, 2.9), expand = FALSE)




ggplot(rt, aes(x = dist.trend_AB_tmean, y = dist.trend_AA_tmean)) +
  geom_point(aes(
    color = coefs_m1_interaction_snp_dist_tmean_p.value < 0.05,
    shape = factor(case_when(
      dist.trend_p.value_AA_tmean < 0.05 & dist.trend_p.value_AB_tmean < 0.05 & species == "oak" ~ "filled_oak",
      dist.trend_p.value_AA_tmean < 0.05 & dist.trend_p.value_AB_tmean < 0.05 & species == "beech" ~ "filled_beech",
      dist.trend_p.value_AA_tmean < 0.05 & dist.trend_p.value_AB_tmean >= 0.05 & species == "oak" ~ "hollow_oak",
      dist.trend_p.value_AA_tmean < 0.05 & dist.trend_p.value_AB_tmean >= 0.05 & species == "beech" ~ "hollow_beech",
      dist.trend_p.value_AA_tmean >= 0.05 & dist.trend_p.value_AB_tmean < 0.05 & species == "oak" ~ "cross_oak",
      dist.trend_p.value_AA_tmean >= 0.05 & dist.trend_p.value_AB_tmean < 0.05 & species == "beech" ~ "cross_beech",
      dist.trend_p.value_AA_tmean >= 0.05 & dist.trend_p.value_AB_tmean >= 0.05 & species == "beech" ~ "thick_outline_beech",
      dist.trend_p.value_AA_tmean >= 0.05 & dist.trend_p.value_AB_tmean >= 0.05 & species == "oak" ~ "thick_outline_oak"
    )),
    fill = factor(case_when(
      dist.trend_p.value_AA_tmean < 0.05 & dist.trend_p.value_AB_tmean < 0.05 ~ "filled",
      dist.trend_p.value_AA_tmean >= 0.05 & dist.trend_p.value_AB_tmean >= 0.05 ~ "not_filled",
      TRUE ~ NA_character_
    ))
  ), size = 4.5, stroke = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(x = "beta_B", y = "beta_A", title = "rt tmea square") +
  scale_shape_manual(values = c(
    "filled_oak" = 16, "filled_beech" = 18, "hollow_oak" = 1, "hollow_beech" = 5,
    "cross_oak" = 10, "cross_beech" = 9, "thick_outline_oak" = 21, "thick_outline_beech" = 23
  )) +
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "blue")) +
  scale_fill_manual(values = c("filled" = "red", "not_filled" = "grey")) +
  theme_minimal() +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),legend.position = "none",
        axis.text.x = element_text(size = 14),axis.text.y = element_text(size = 14),
        axis.ticks = element_line(color = "black", linewidth = 0.8),axis.ticks.length = unit(4, "mm"),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1.5)
  ) +
  coord_fixed(ratio = 1, xlim = c(-2.9, 2.9), ylim = c(-2.9, 2.9), expand = FALSE)




ggplot(rc, aes(x = dist.trend_AB_tmean, y = dist.trend_AA_tmean)) +
  geom_point(aes(
    color = coefs_m1_interaction_snp_dist_tmean_p.value < 0.05,
    shape = factor(case_when(
      dist.trend_p.value_AA_tmean < 0.05 & dist.trend_p.value_AB_tmean < 0.05 & species == "oak" ~ "filled_oak",
      dist.trend_p.value_AA_tmean < 0.05 & dist.trend_p.value_AB_tmean < 0.05 & species == "beech" ~ "filled_beech",
      dist.trend_p.value_AA_tmean < 0.05 & dist.trend_p.value_AB_tmean >= 0.05 & species == "oak" ~ "hollow_oak",
      dist.trend_p.value_AA_tmean < 0.05 & dist.trend_p.value_AB_tmean >= 0.05 & species == "beech" ~ "hollow_beech",
      dist.trend_p.value_AA_tmean >= 0.05 & dist.trend_p.value_AB_tmean < 0.05 & species == "oak" ~ "cross_oak",
      dist.trend_p.value_AA_tmean >= 0.05 & dist.trend_p.value_AB_tmean < 0.05 & species == "beech" ~ "cross_beech",
      dist.trend_p.value_AA_tmean >= 0.05 & dist.trend_p.value_AB_tmean >= 0.05 & species == "beech" ~ "thick_outline_beech",
      dist.trend_p.value_AA_tmean >= 0.05 & dist.trend_p.value_AB_tmean >= 0.05 & species == "oak" ~ "thick_outline_oak"
    )),
    fill = factor(case_when(
      dist.trend_p.value_AA_tmean < 0.05 & dist.trend_p.value_AB_tmean < 0.05 ~ "filled",
      dist.trend_p.value_AA_tmean >= 0.05 & dist.trend_p.value_AB_tmean >= 0.05 ~ "not_filled",
      TRUE ~ NA_character_
    ))
  ), size = 4.5, stroke = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(x = "beta_B", y = "beta_A", title = "rc tmea square") +
  scale_shape_manual(values = c(
    "filled_oak" = 16, "filled_beech" = 18, "hollow_oak" = 1, "hollow_beech" = 5,
    "cross_oak" = 10, "cross_beech" = 9, "thick_outline_oak" = 21, "thick_outline_beech" = 23
  )) +
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "blue")) +
  scale_fill_manual(values = c("filled" = "red", "not_filled" = "grey")) +
  theme_minimal() +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),legend.position = "none",
        axis.text.x = element_text(size = 14),axis.text.y = element_text(size = 14),
        axis.ticks = element_line(color = "black", linewidth = 0.8),axis.ticks.length = unit(4, "mm"),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1.5)
  ) +
  coord_fixed(ratio = 1, xlim = c(-2.9, 2.9), ylim = c(-2.9, 2.9), expand = FALSE)










ggplot(sprtmea, aes(x = dist.trend_AB_tmean, y = dist.trend_AA_tmean)) +
  geom_point(aes(color = coefs_m1_interaction_snp_dist_tmean_p.value < 0.05,shape = factor(case_when(
    dist.trend_p.value_AA_tmean < 0.05 & dist.trend_p.value_AB_tmean < 0.05 & species == "oak" ~ "filled_oak",
    dist.trend_p.value_AA_tmean < 0.05 & dist.trend_p.value_AB_tmean < 0.05 & species == "beech" ~ "filled_beech",
    dist.trend_p.value_AA_tmean < 0.05 & dist.trend_p.value_AB_tmean >= 0.05 & species == "oak" ~ "hollow_oak",
    dist.trend_p.value_AA_tmean < 0.05 & dist.trend_p.value_AB_tmean >= 0.05 & species == "beech" ~ "hollow_beech",
    dist.trend_p.value_AA_tmean >= 0.05 & dist.trend_p.value_AB_tmean < 0.05 & species == "oak" ~ "cross_oak",
    dist.trend_p.value_AA_tmean >= 0.05 & dist.trend_p.value_AB_tmean < 0.05 & species == "beech" ~ "cross_beech",
    dist.trend_p.value_AA_tmean >= 0.05 & dist.trend_p.value_AB_tmean >= 0.05 & species == "beech" ~ "thick_outline_beech",
    dist.trend_p.value_AA_tmean >= 0.05 & dist.trend_p.value_AB_tmean >= 0.05 & species == "oak" ~ "thick_outline_oak")),
    fill = factor(case_when(dist.trend_p.value_AA_tmean < 0.05 & dist.trend_p.value_AB_tmean < 0.05 ~ "filled",
                            dist.trend_p.value_AA_tmean >= 0.05 & dist.trend_p.value_AB_tmean >= 0.05 ~ "not_filled",
                            TRUE ~ NA_character_)),size = 3, stroke = 1.5)) +
  geom_hline(yintercept = 0, linetype = "dashed") +geom_vline(xintercept = 0, linetype = "dashed") +
  xlim(-2.5, 2.5) +ylim(-2.5, 2.5) +labs(x = "beta_B", y = "beta_A", title = "sprtmea tmea") +
  scale_shape_manual(values = c("filled_oak" = 16, "filled_beech" = 18, "hollow_oak" = 1, "hollow_beech" = 5,
                                "cross_oak" = 10, "cross_beech" = 9, "thick_outline_oak" = 21, "thick_outline_beech" = 23)) +
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "blue")) +
  scale_fill_manual(values = c("filled" = "red", "not_filled" = "grey")) + # Set "not_filled" to grey
  theme_minimal() +theme(panel.grid.major = element_blank(),
                         panel.grid.minor = element_blank(),legend.position = "none",axis.text.x = element_text(size = 14),
                         axis.text.y = element_text(size = 14),panel.border = element_rect(color = "black", fill = NA, size = 1.5))


ggplot(sumtmea, aes(x = dist.trend_AB_tmean, y = dist.trend_AA_tmean)) +
  geom_point(aes(color = coefs_m1_interaction_snp_dist_tmean_p.value < 0.05,shape = factor(case_when(
    dist.trend_p.value_AA_tmean < 0.05 & dist.trend_p.value_AB_tmean < 0.05 & species == "oak" ~ "filled_oak",
    dist.trend_p.value_AA_tmean < 0.05 & dist.trend_p.value_AB_tmean < 0.05 & species == "beech" ~ "filled_beech",
    dist.trend_p.value_AA_tmean < 0.05 & dist.trend_p.value_AB_tmean >= 0.05 & species == "oak" ~ "hollow_oak",
    dist.trend_p.value_AA_tmean < 0.05 & dist.trend_p.value_AB_tmean >= 0.05 & species == "beech" ~ "hollow_beech",
    dist.trend_p.value_AA_tmean >= 0.05 & dist.trend_p.value_AB_tmean < 0.05 & species == "oak" ~ "cross_oak",
    dist.trend_p.value_AA_tmean >= 0.05 & dist.trend_p.value_AB_tmean < 0.05 & species == "beech" ~ "cross_beech",
    dist.trend_p.value_AA_tmean >= 0.05 & dist.trend_p.value_AB_tmean >= 0.05 & species == "beech" ~ "thick_outline_beech",
    dist.trend_p.value_AA_tmean >= 0.05 & dist.trend_p.value_AB_tmean >= 0.05 & species == "oak" ~ "thick_outline_oak")),
    fill = factor(case_when(dist.trend_p.value_AA_tmean < 0.05 & dist.trend_p.value_AB_tmean < 0.05 ~ "filled",
                            dist.trend_p.value_AA_tmean >= 0.05 & dist.trend_p.value_AB_tmean >= 0.05 ~ "not_filled",
                            TRUE ~ NA_character_)),size = 3, stroke = 1.5)) +
  geom_hline(yintercept = 0, linetype = "dashed") +geom_vline(xintercept = 0, linetype = "dashed") +
  xlim(-2.5, 2.5) +ylim(-2.5, 2.5) +labs(x = "beta_B", y = "beta_A", title = "sumtmea tmea") +
  scale_shape_manual(values = c("filled_oak" = 16, "filled_beech" = 18, "hollow_oak" = 1, "hollow_beech" = 5,
                                "cross_oak" = 10, "cross_beech" = 9, "thick_outline_oak" = 21, "thick_outline_beech" = 23)) +
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "blue")) +
  scale_fill_manual(values = c("filled" = "red", "not_filled" = "grey")) + # Set "not_filled" to grey
  theme_minimal() +theme(panel.grid.major = element_blank(),
                         panel.grid.minor = element_blank(),legend.position = "none",axis.text.x = element_text(size = 14),
                         axis.text.y = element_text(size = 14),panel.border = element_rect(color = "black", fill = NA, size = 1.5))



ggplot(sprpre, aes(x = dist.trend_AB_tmean, y = dist.trend_AA_tmean)) +
  geom_point(aes(color = coefs_m1_interaction_snp_dist_tmean_p.value < 0.05,shape = factor(case_when(
    dist.trend_p.value_AA_tmean < 0.05 & dist.trend_p.value_AB_tmean < 0.05 & species == "oak" ~ "filled_oak",
    dist.trend_p.value_AA_tmean < 0.05 & dist.trend_p.value_AB_tmean < 0.05 & species == "beech" ~ "filled_beech",
    dist.trend_p.value_AA_tmean < 0.05 & dist.trend_p.value_AB_tmean >= 0.05 & species == "oak" ~ "hollow_oak",
    dist.trend_p.value_AA_tmean < 0.05 & dist.trend_p.value_AB_tmean >= 0.05 & species == "beech" ~ "hollow_beech",
    dist.trend_p.value_AA_tmean >= 0.05 & dist.trend_p.value_AB_tmean < 0.05 & species == "oak" ~ "cross_oak",
    dist.trend_p.value_AA_tmean >= 0.05 & dist.trend_p.value_AB_tmean < 0.05 & species == "beech" ~ "cross_beech",
    dist.trend_p.value_AA_tmean >= 0.05 & dist.trend_p.value_AB_tmean >= 0.05 & species == "beech" ~ "thick_outline_beech",
    dist.trend_p.value_AA_tmean >= 0.05 & dist.trend_p.value_AB_tmean >= 0.05 & species == "oak" ~ "thick_outline_oak")),
    fill = factor(case_when(dist.trend_p.value_AA_tmean < 0.05 & dist.trend_p.value_AB_tmean < 0.05 ~ "filled",
                            dist.trend_p.value_AA_tmean >= 0.05 & dist.trend_p.value_AB_tmean >= 0.05 ~ "not_filled",
                            TRUE ~ NA_character_)),size = 3, stroke = 1.5)) +
  geom_hline(yintercept = 0, linetype = "dashed") +geom_vline(xintercept = 0, linetype = "dashed") +
  xlim(-2.5, 2.5) +ylim(-2.5, 2.5) +labs(x = "beta_B", y = "beta_A", title = "sprpre tmea") +
  scale_shape_manual(values = c("filled_oak" = 16, "filled_beech" = 18, "hollow_oak" = 1, "hollow_beech" = 5,
                                "cross_oak" = 10, "cross_beech" = 9, "thick_outline_oak" = 21, "thick_outline_beech" = 23)) +
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "blue")) +
  scale_fill_manual(values = c("filled" = "red", "not_filled" = "grey")) + # Set "not_filled" to grey
  theme_minimal() +theme(panel.grid.major = element_blank(),
                         panel.grid.minor = element_blank(),legend.position = "none",axis.text.x = element_text(size = 14),
                         axis.text.y = element_text(size = 14),panel.border = element_rect(color = "black", fill = NA, size = 1.5))


ggplot(sumpre, aes(x = dist.trend_AB_tmean, y = dist.trend_AA_tmean)) +
  geom_point(aes(color = coefs_m1_interaction_snp_dist_tmean_p.value < 0.05,shape = factor(case_when(
    dist.trend_p.value_AA_tmean < 0.05 & dist.trend_p.value_AB_tmean < 0.05 & species == "oak" ~ "filled_oak",
    dist.trend_p.value_AA_tmean < 0.05 & dist.trend_p.value_AB_tmean < 0.05 & species == "beech" ~ "filled_beech",
    dist.trend_p.value_AA_tmean < 0.05 & dist.trend_p.value_AB_tmean >= 0.05 & species == "oak" ~ "hollow_oak",
    dist.trend_p.value_AA_tmean < 0.05 & dist.trend_p.value_AB_tmean >= 0.05 & species == "beech" ~ "hollow_beech",
    dist.trend_p.value_AA_tmean >= 0.05 & dist.trend_p.value_AB_tmean < 0.05 & species == "oak" ~ "cross_oak",
    dist.trend_p.value_AA_tmean >= 0.05 & dist.trend_p.value_AB_tmean < 0.05 & species == "beech" ~ "cross_beech",
    dist.trend_p.value_AA_tmean >= 0.05 & dist.trend_p.value_AB_tmean >= 0.05 & species == "beech" ~ "thick_outline_beech",
    dist.trend_p.value_AA_tmean >= 0.05 & dist.trend_p.value_AB_tmean >= 0.05 & species == "oak" ~ "thick_outline_oak")),
    fill = factor(case_when(dist.trend_p.value_AA_tmean < 0.05 & dist.trend_p.value_AB_tmean < 0.05 ~ "filled",
                            dist.trend_p.value_AA_tmean >= 0.05 & dist.trend_p.value_AB_tmean >= 0.05 ~ "not_filled",
                            TRUE ~ NA_character_)),size = 3, stroke = 1.5)) +
  geom_hline(yintercept = 0, linetype = "dashed") +geom_vline(xintercept = 0, linetype = "dashed") +
  xlim(-2.5, 2.5) +ylim(-2.5, 2.5) +labs(x = "beta_B", y = "beta_A", title = "sumpre tmea") +
  scale_shape_manual(values = c("filled_oak" = 16, "filled_beech" = 18, "hollow_oak" = 1, "hollow_beech" = 5,
                                "cross_oak" = 10, "cross_beech" = 9, "thick_outline_oak" = 21, "thick_outline_beech" = 23)) +
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "blue")) +
  scale_fill_manual(values = c("filled" = "red", "not_filled" = "grey")) + # Set "not_filled" to grey
  theme_minimal() +theme(panel.grid.major = element_blank(),
                         panel.grid.minor = element_blank(),legend.position = "none",axis.text.x = element_text(size = 14),
                         axis.text.y = element_text(size = 14),panel.border = element_rect(color = "black", fill = NA, size = 1.5))



rm(list = ls())
gc()
# Load necessary libraries
library(dplyr)
library(tidyr)
library(ggplot2)
#setwd("/Users/kaffemobil/Desktop/wsl_from_kaffemobil/Aksel_061224/chelsa_data_redownloaded_to_match_elis_paper_020824/calculate_bio_clims_070824/")
#setwd("/Users/apalsson/Downloads/wsl_paper_work_01225 2/calculate_bio_clims_070824/")
setwd("/Users/apalsson/Desktop/WSL_DATA_FOR_PUBLICATION_120526/WSL_SCRIPT_TO_UPLOAD/raw_gxe_beech/testing_script/output_test/dryad")
bai_beech <- read.csv("dist_prec_snp_tmean_as_cov_101224_beech_snp_dist_STD_tmean_dist_prec_as_int_results_all_bai_041224_20260519.csv")
sla_beech <- read.csv("dist_prec_snp_tmean_as_cov_101224_beech_snp_dist_STD_tmean_dist_prec_as_int_results_all_sla_041224_20260519.csv")
height_beech <- read.csv("dist_prec_snp_tmean_as_cov_101224_beech_snp_dist_STD_tmean_dist_prec_as_int_results_all_height_041224_20260519.csv")
dbh_beech <- read.csv("dist_prec_snp_tmean_as_cov_101224_beech_snp_dist_STD_tmean_dist_prec_as_int_results_all_dbh_041224_20260519.csv")
rs_beech <- read.csv("dist_prec_snp_tmean_as_cov_101224_beech_snp_dist_STD_tmean_dist_prec_as_int_results_all_rs_041224_20260519.csv")
rc_beech <- read.csv("dist_prec_snp_tmean_as_cov_101224_beech_snp_dist_STD_tmean_dist_prec_as_int_results_all_rc_041224_20260519.csv")
rt_beech <- read.csv("dist_prec_snp_tmean_as_cov_101224_beech_snp_dist_STD_tmean_dist_prec_as_int_results_all_rt_041224_20260519.csv")
sprpre_beech <- read.csv("dist_prec_snp_tmean_as_cov_101224_beech_snp_dist_STD_tmean_dist_prec_as_int_results_all_prec_mean_spr_051224_20260519.csv")
sumpre_beech <- read.csv("dist_prec_snp_tmean_as_cov_101224_beech_snp_dist_STD_tmean_dist_prec_as_int_results_all_prec_mean_sum_051224_20260519.csv")
sumtmea_beech <- read.csv("dist_prec_snp_tmean_as_cov_101224_beech_snp_dist_STD_tmean_dist_prec_as_int_results_all_tmean_mean_sum_051224_20260519.csv")
sprtmea_beech <- read.csv("dist_prec_snp_tmean_as_cov_101224_beech_snp_dist_STD_tmean_dist_prec_as_int_results_all_tmean_mean_spr_051224_20260519.csv")
#
setwd("/Users/apalsson/Desktop/WSL_DATA_FOR_PUBLICATION_120526/WSL_SCRIPT_TO_UPLOAD/raw_gxe_oak/testing_script/output_test/dryad")
bai_oak <- read.csv("dist_prec_snp_tmean_as_cov_101224_oak_std_preds_snp_dist_tmean_dist_prec_as_int_results_all_bai_051224.csv")
sla_oak <- read.csv("dist_prec_snp_tmean_as_cov_101224_oak_std_preds_snp_dist_tmean_dist_prec_as_int_results_all_sla_051224.csv")
height_oak <- read.csv("dist_prec_snp_tmean_as_cov_101224_oak_std_preds_snp_dist_tmean_dist_prec_as_int_results_all_height_051224.csv")
dbh_oak <- read.csv("dist_prec_snp_tmean_as_cov_101224_oak_std_preds_snp_dist_tmean_dist_prec_as_int_results_all_dbh_051224.csv")
rs_oak <- read.csv("dist_prec_snp_tmean_as_cov_101224_oak_std_preds_snp_dist_tmean_dist_prec_as_int_results_all_rs_051224.csv")
rc_oak <- read.csv("dist_prec_snp_tmean_as_cov_101224_oak_std_preds_snp_dist_tmean_dist_prec_as_int_results_all_rc_051224.csv")
rt_oak <- read.csv("dist_prec_snp_tmean_as_cov_101224_oak_std_preds_snp_dist_tmean_dist_prec_as_int_results_all_rt_051224.csv")
sprpre_oak <- read.csv("dist_prec_snp_tmean_as_cov_101224_oak_std_preds_snp_dist_tmean_dist_prec_as_int_results_all_sprpre_051224.csv")
sumpre_oak <- read.csv("dist_prec_snp_tmean_as_cov_101224_oak_std_preds_snp_dist_tmean_dist_prec_as_int_results_all_sumpre_051224.csv")
sumtmea_oak <- read.csv("dist_prec_snp_tmean_as_cov_101224_oak_std_preds_snp_dist_tmean_dist_prec_as_int_results_all_sumtmea_051224.csv")
sprtmea_oak <- read.csv("dist_prec_snp_tmean_as_cov_101224_oak_std_preds_snp_dist_tmean_dist_prec_as_int_results_all_sprtmea_051224.csv")


#setwd("/Users/apalsson/Desktop/wsl/chelsa_data_redownloaded_to_match_elis_paper_020824/calculate_bio_clims_070824/")

# Filter each DataFrame to include only unique snp_id entries
sumpre_oak <- distinct(sumpre_oak, snp_id, .keep_all = TRUE)
sumtmea_oak <- distinct(sumtmea_oak, snp_id, .keep_all = TRUE)
sprtmea_oak <- distinct(sprtmea_oak, snp_id, .keep_all = TRUE)
bai_oak <- distinct(bai_oak, snp_id, .keep_all = TRUE)
sla_oak <- distinct(sla_oak, snp_id, .keep_all = TRUE)
height_oak <- distinct(height_oak, snp_id, .keep_all = TRUE)
dbh_oak <- distinct(dbh_oak, snp_id, .keep_all = TRUE)
rs_oak <- distinct(rs_oak, snp_id, .keep_all = TRUE)
rc_oak <- distinct(rc_oak, snp_id, .keep_all = TRUE)
rt_oak <- distinct(rt_oak, snp_id, .keep_all = TRUE)
sprpre_oak <- distinct(sprpre_oak, snp_id, .keep_all = TRUE)
sumtmea_beech <- distinct(sumtmea_beech, snp_id, .keep_all = TRUE)
sprtmea_beech <- distinct(sprtmea_beech, snp_id, .keep_all = TRUE)
bai_beech <- distinct(bai_beech, snp_id, .keep_all = TRUE)
sla_beech <- distinct(sla_beech, snp_id, .keep_all = TRUE)
height_beech <- distinct(height_beech, snp_id, .keep_all = TRUE)
dbh_beech <- distinct(dbh_beech, snp_id, .keep_all = TRUE)
rs_beech <- distinct(rs_beech, snp_id, .keep_all = TRUE)
rc_beech <- distinct(rc_beech, snp_id, .keep_all = TRUE)
rt_beech <- distinct(rt_beech, snp_id, .keep_all = TRUE)
sprpre_beech <- distinct(sprpre_beech, snp_id, .keep_all = TRUE)

# Now, each DataFrame has only unique snp_id entries
sumpre_oak$species <- "oak";sumtmea_oak$species<- "oak";sprtmea_oak$species<- "oak";bai_oak$species<- "oak";sla_oak$species<- "oak"
height_oak$species<- "oak";dbh_oak$species<- "oak";rs_oak$species<- "oak";rc_oak$species<- "oak";rt_oak$species<- "oak"
sumpre_beech$species <- "beech";sumtmea_beech$species<- "beech";sprtmea_beech$species<- "beech";bai_beech$species<- "beech"
sla_beech$species<- "beech";height_beech$species<- "beech";dbh_beech$species<- "beech";rs_beech$species<- "beech"
rc_beech$species<- "beech";rt_beech$species<- "beech";sprpre_oak$species <- "oak";sprpre_beech$species <- "beech"

colnames(sumpre_beech)
colnames(sumpre_oak)
sumpre <- rbind(sumpre_oak, sumpre_beech)
sumtmea <- rbind(sumtmea_oak, sumtmea_beech)
bai <- rbind(bai_oak, bai_beech)
sla <- rbind(sla_oak, sla_beech)
height <- rbind(height_oak, height_beech)
dbh <-rbind(dbh_oak, dbh_beech)
rs <- rbind(rs_oak, rs_beech)
rt <- rbind(rt_oak, rt_beech)
rc <- rbind(rc_oak,rc_beech)
sprtmea <- rbind(sprtmea_oak, sprtmea_beech)
sprpre <- rbind(sprpre_oak, sprpre_beech)

invert_trait_effects <- function(df) {
  df$dist.trend_AA_tmean <- -df$dist.trend_AA_tmean
  df$dist.trend_AB_tmean <- -df$dist.trend_AB_tmean
  
  if (all(c("dist.trend_lower.CL_AA_tmean", "dist.trend_upper.CL_AA_tmean") %in% colnames(df))) {
    lower_aa <- df$dist.trend_lower.CL_AA_tmean
    upper_aa <- df$dist.trend_upper.CL_AA_tmean
    
    df$dist.trend_lower.CL_AA_tmean <- -upper_aa
    df$dist.trend_upper.CL_AA_tmean <- -lower_aa
  }
  
  if (all(c("dist.trend_lower.CL_AB_tmean", "dist.trend_upper.CL_AB_tmean") %in% colnames(df))) {
    lower_ab <- df$dist.trend_lower.CL_AB_tmean
    upper_ab <- df$dist.trend_upper.CL_AB_tmean
    
    df$dist.trend_lower.CL_AB_tmean <- -upper_ab
    df$dist.trend_upper.CL_AB_tmean <- -lower_ab
  }
  
  if ("contrast_tmean_estimate" %in% colnames(df)) {
    df$contrast_tmean_estimate <- -df$contrast_tmean_estimate
  }
  
  if (all(c("contrast_tmean_estimate_lower_CL", "contrast_tmean_estimate_upper_CL") %in% colnames(df))) {
    lower_contrast <- df$contrast_tmean_estimate_lower_CL
    upper_contrast <- df$contrast_tmean_estimate_upper_CL
    
    df$contrast_tmean_estimate_lower_CL <- -upper_contrast
    df$contrast_tmean_estimate_upper_CL <- -lower_contrast
  }
  
  if ("coefs_m1_interaction_snp_dist_tmean_estimate" %in% colnames(df)) {
    df$coefs_m1_interaction_snp_dist_tmean_estimate <- -df$coefs_m1_interaction_snp_dist_tmean_estimate
  }
  
  df
}

sumpre <- invert_trait_effects(sumpre)
sprpre <- invert_trait_effects(sprpre)
sumtmea <- invert_trait_effects(sumtmea)
sprtmea <- invert_trait_effects(sprtmea)
summary(sumpre$dist.trend_AA_tmean)
summary(sumpre$dist.trend_AB_tmean)
summary(sprpre$dist.trend_AA_tmean)
summary(sprpre$dist.trend_AB_tmean)
summary(sumtmea$dist.trend_AA_tmean)
summary(sumtmea$dist.trend_AB_tmean)
summary(sprtmea$dist.trend_AA_tmean)
summary(sprtmea$dist.trend_AB_tmean)
summary(bai$dist.trend_AA_tmean)
summary(bai$dist.trend_AB_tmean)
summary(dbh$dist.trend_AA_tmean)
summary(dbh$dist.trend_AB_tmean)
summary(sla$dist.trend_AA_tmean)
summary(sla$dist.trend_AB_tmean)
summary(height$dist.trend_AA_tmean)
summary(height$dist.trend_AB_tmean)
summary(rc$dist.trend_AA_tmean)
summary(rc$dist.trend_AB_tmean)
summary(rt$dist.trend_AA_tmean)
summary(rt$dist.trend_AB_tmean)
summary(rs$dist.trend_AA_tmean)
summary(rs$dist.trend_AB_tmean)



ggplot(dbh, aes(x = dist.trend_AB_tmean, y = dist.trend_AA_tmean)) +
  geom_point(aes(color = coefs_m1_interaction_snp_dist_tmean_p.value < 0.05,shape = factor(case_when(
    dist.trend_p.value_AA_tmean < 0.05 & dist.trend_p.value_AB_tmean < 0.05 & species == "oak" ~ "filled_oak",
    dist.trend_p.value_AA_tmean < 0.05 & dist.trend_p.value_AB_tmean < 0.05 & species == "beech" ~ "filled_beech",
    dist.trend_p.value_AA_tmean < 0.05 & dist.trend_p.value_AB_tmean >= 0.05 & species == "oak" ~ "hollow_oak",
    dist.trend_p.value_AA_tmean < 0.05 & dist.trend_p.value_AB_tmean >= 0.05 & species == "beech" ~ "hollow_beech",
    dist.trend_p.value_AA_tmean >= 0.05 & dist.trend_p.value_AB_tmean < 0.05 & species == "oak" ~ "cross_oak",
    dist.trend_p.value_AA_tmean >= 0.05 & dist.trend_p.value_AB_tmean < 0.05 & species == "beech" ~ "cross_beech",
    dist.trend_p.value_AA_tmean >= 0.05 & dist.trend_p.value_AB_tmean >= 0.05 & species == "beech" ~ "thick_outline_beech",
    dist.trend_p.value_AA_tmean >= 0.05 & dist.trend_p.value_AB_tmean >= 0.05 & species == "oak" ~ "thick_outline_oak")),
    fill = factor(case_when(dist.trend_p.value_AA_tmean < 0.05 & dist.trend_p.value_AB_tmean < 0.05 ~ "filled",
                            dist.trend_p.value_AA_tmean >= 0.05 & dist.trend_p.value_AB_tmean >= 0.05 ~ "not_filled",
                            TRUE ~ NA_character_)),size = 3, stroke = 1.5)) +
  geom_hline(yintercept = 0, linetype = "dashed") +geom_vline(xintercept = 0, linetype = "dashed") +
  xlim(-2.5, 2.5) +ylim(-2.5, 2.5) +labs(x = "beta_B", y = "beta_A", title = "dbh prec") +
  scale_shape_manual(values = c("filled_oak" = 16, "filled_beech" = 18, "hollow_oak" = 1, "hollow_beech" = 5,
                                "cross_oak" = 10, "cross_beech" = 9, "thick_outline_oak" = 21, "thick_outline_beech" = 23)) +
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "blue")) +
  scale_fill_manual(values = c("filled" = "red", "not_filled" = "grey")) + # Set "not_filled" to grey
  theme_minimal() +theme(panel.grid.major = element_blank(),
                         panel.grid.minor = element_blank(),legend.position = "none",axis.text.x = element_text(size = 14),
                         axis.text.y = element_text(size = 14),panel.border = element_rect(color = "black", fill = NA, size = 1.5))


ggplot(height, aes(x = dist.trend_AB_tmean, y = dist.trend_AA_tmean)) +
  geom_point(aes(color = coefs_m1_interaction_snp_dist_tmean_p.value < 0.05,shape = factor(case_when(
    dist.trend_p.value_AA_tmean < 0.05 & dist.trend_p.value_AB_tmean < 0.05 & species == "oak" ~ "filled_oak",
    dist.trend_p.value_AA_tmean < 0.05 & dist.trend_p.value_AB_tmean < 0.05 & species == "beech" ~ "filled_beech",
    dist.trend_p.value_AA_tmean < 0.05 & dist.trend_p.value_AB_tmean >= 0.05 & species == "oak" ~ "hollow_oak",
    dist.trend_p.value_AA_tmean < 0.05 & dist.trend_p.value_AB_tmean >= 0.05 & species == "beech" ~ "hollow_beech",
    dist.trend_p.value_AA_tmean >= 0.05 & dist.trend_p.value_AB_tmean < 0.05 & species == "oak" ~ "cross_oak",
    dist.trend_p.value_AA_tmean >= 0.05 & dist.trend_p.value_AB_tmean < 0.05 & species == "beech" ~ "cross_beech",
    dist.trend_p.value_AA_tmean >= 0.05 & dist.trend_p.value_AB_tmean >= 0.05 & species == "beech" ~ "thick_outline_beech",
    dist.trend_p.value_AA_tmean >= 0.05 & dist.trend_p.value_AB_tmean >= 0.05 & species == "oak" ~ "thick_outline_oak")),
    fill = factor(case_when(dist.trend_p.value_AA_tmean < 0.05 & dist.trend_p.value_AB_tmean < 0.05 ~ "filled",
                            dist.trend_p.value_AA_tmean >= 0.05 & dist.trend_p.value_AB_tmean >= 0.05 ~ "not_filled",
                            TRUE ~ NA_character_)),size = 3, stroke = 1.5)) +
  geom_hline(yintercept = 0, linetype = "dashed") +geom_vline(xintercept = 0, linetype = "dashed") +
  xlim(-2.5, 2.5) +ylim(-2.5, 2.5) +labs(x = "beta_B", y = "beta_A", title = "height prec") +
  scale_shape_manual(values = c("filled_oak" = 16, "filled_beech" = 18, "hollow_oak" = 1, "hollow_beech" = 5,
                                "cross_oak" = 10, "cross_beech" = 9, "thick_outline_oak" = 21, "thick_outline_beech" = 23)) +
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "blue")) +
  scale_fill_manual(values = c("filled" = "red", "not_filled" = "grey")) + # Set "not_filled" to grey
  theme_minimal() +theme(panel.grid.major = element_blank(),
                         panel.grid.minor = element_blank(),legend.position = "none",axis.text.x = element_text(size = 14),
                         axis.text.y = element_text(size = 14),panel.border = element_rect(color = "black", fill = NA, size = 1.5))




ggplot(bai, aes(x = dist.trend_AB_tmean, y = dist.trend_AA_tmean)) +
  geom_point(aes(color = coefs_m1_interaction_snp_dist_tmean_p.value < 0.05,shape = factor(case_when(
    dist.trend_p.value_AA_tmean < 0.05 & dist.trend_p.value_AB_tmean < 0.05 & species == "oak" ~ "filled_oak",
    dist.trend_p.value_AA_tmean < 0.05 & dist.trend_p.value_AB_tmean < 0.05 & species == "beech" ~ "filled_beech",
    dist.trend_p.value_AA_tmean < 0.05 & dist.trend_p.value_AB_tmean >= 0.05 & species == "oak" ~ "hollow_oak",
    dist.trend_p.value_AA_tmean < 0.05 & dist.trend_p.value_AB_tmean >= 0.05 & species == "beech" ~ "hollow_beech",
    dist.trend_p.value_AA_tmean >= 0.05 & dist.trend_p.value_AB_tmean < 0.05 & species == "oak" ~ "cross_oak",
    dist.trend_p.value_AA_tmean >= 0.05 & dist.trend_p.value_AB_tmean < 0.05 & species == "beech" ~ "cross_beech",
    dist.trend_p.value_AA_tmean >= 0.05 & dist.trend_p.value_AB_tmean >= 0.05 & species == "beech" ~ "thick_outline_beech",
    dist.trend_p.value_AA_tmean >= 0.05 & dist.trend_p.value_AB_tmean >= 0.05 & species == "oak" ~ "thick_outline_oak")),
    fill = factor(case_when(dist.trend_p.value_AA_tmean < 0.05 & dist.trend_p.value_AB_tmean < 0.05 ~ "filled",
                            dist.trend_p.value_AA_tmean >= 0.05 & dist.trend_p.value_AB_tmean >= 0.05 ~ "not_filled",
                            TRUE ~ NA_character_)),size = 3, stroke = 1.5)) +
  geom_hline(yintercept = 0, linetype = "dashed") +geom_vline(xintercept = 0, linetype = "dashed") +
  xlim(-2.5, 2.5) +ylim(-2.5, 2.5) +labs(x = "beta_B", y = "beta_A", title = "bai prec") +
  scale_shape_manual(values = c("filled_oak" = 16, "filled_beech" = 18, "hollow_oak" = 1, "hollow_beech" = 5,
                                "cross_oak" = 10, "cross_beech" = 9, "thick_outline_oak" = 21, "thick_outline_beech" = 23)) +
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "blue")) +
  scale_fill_manual(values = c("filled" = "red", "not_filled" = "grey")) + # Set "not_filled" to grey
  theme_minimal() +theme(panel.grid.major = element_blank(),
                         panel.grid.minor = element_blank(),legend.position = "none",axis.text.x = element_text(size = 14),
                         axis.text.y = element_text(size = 14),panel.border = element_rect(color = "black", fill = NA, size = 1.5))



ggplot(sla, aes(x = dist.trend_AB_tmean, y = dist.trend_AA_tmean)) +
  geom_point(aes(color = coefs_m1_interaction_snp_dist_tmean_p.value < 0.05,shape = factor(case_when(
    dist.trend_p.value_AA_tmean < 0.05 & dist.trend_p.value_AB_tmean < 0.05 & species == "oak" ~ "filled_oak",
    dist.trend_p.value_AA_tmean < 0.05 & dist.trend_p.value_AB_tmean < 0.05 & species == "beech" ~ "filled_beech",
    dist.trend_p.value_AA_tmean < 0.05 & dist.trend_p.value_AB_tmean >= 0.05 & species == "oak" ~ "hollow_oak",
    dist.trend_p.value_AA_tmean < 0.05 & dist.trend_p.value_AB_tmean >= 0.05 & species == "beech" ~ "hollow_beech",
    dist.trend_p.value_AA_tmean >= 0.05 & dist.trend_p.value_AB_tmean < 0.05 & species == "oak" ~ "cross_oak",
    dist.trend_p.value_AA_tmean >= 0.05 & dist.trend_p.value_AB_tmean < 0.05 & species == "beech" ~ "cross_beech",
    dist.trend_p.value_AA_tmean >= 0.05 & dist.trend_p.value_AB_tmean >= 0.05 & species == "beech" ~ "thick_outline_beech",
    dist.trend_p.value_AA_tmean >= 0.05 & dist.trend_p.value_AB_tmean >= 0.05 & species == "oak" ~ "thick_outline_oak")),
    fill = factor(case_when(dist.trend_p.value_AA_tmean < 0.05 & dist.trend_p.value_AB_tmean < 0.05 ~ "filled",
                            dist.trend_p.value_AA_tmean >= 0.05 & dist.trend_p.value_AB_tmean >= 0.05 ~ "not_filled",
                            TRUE ~ NA_character_)),size = 3, stroke = 1.5)) +
  geom_hline(yintercept = 0, linetype = "dashed") +geom_vline(xintercept = 0, linetype = "dashed") +
  xlim(-2.5, 2.5) +ylim(-2.5, 2.5) +labs(x = "beta_B", y = "beta_A", title = "sla prec") +
  scale_shape_manual(values = c("filled_oak" = 16, "filled_beech" = 18, "hollow_oak" = 1, "hollow_beech" = 5,
                                "cross_oak" = 10, "cross_beech" = 9, "thick_outline_oak" = 21, "thick_outline_beech" = 23)) +
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "blue")) +
  scale_fill_manual(values = c("filled" = "red", "not_filled" = "grey")) + # Set "not_filled" to grey
  theme_minimal() +theme(panel.grid.major = element_blank(),
                         panel.grid.minor = element_blank(),legend.position = "none",axis.text.x = element_text(size = 14),
                         axis.text.y = element_text(size = 14),panel.border = element_rect(color = "black", fill = NA, size = 1.5))




ggplot(rc, aes(x = dist.trend_AB_tmean, y = dist.trend_AA_tmean)) +
  geom_point(aes(color = coefs_m1_interaction_snp_dist_tmean_p.value < 0.05,shape = factor(case_when(
    dist.trend_p.value_AA_tmean < 0.05 & dist.trend_p.value_AB_tmean < 0.05 & species == "oak" ~ "filled_oak",
    dist.trend_p.value_AA_tmean < 0.05 & dist.trend_p.value_AB_tmean < 0.05 & species == "beech" ~ "filled_beech",
    dist.trend_p.value_AA_tmean < 0.05 & dist.trend_p.value_AB_tmean >= 0.05 & species == "oak" ~ "hollow_oak",
    dist.trend_p.value_AA_tmean < 0.05 & dist.trend_p.value_AB_tmean >= 0.05 & species == "beech" ~ "hollow_beech",
    dist.trend_p.value_AA_tmean >= 0.05 & dist.trend_p.value_AB_tmean < 0.05 & species == "oak" ~ "cross_oak",
    dist.trend_p.value_AA_tmean >= 0.05 & dist.trend_p.value_AB_tmean < 0.05 & species == "beech" ~ "cross_beech",
    dist.trend_p.value_AA_tmean >= 0.05 & dist.trend_p.value_AB_tmean >= 0.05 & species == "beech" ~ "thick_outline_beech",
    dist.trend_p.value_AA_tmean >= 0.05 & dist.trend_p.value_AB_tmean >= 0.05 & species == "oak" ~ "thick_outline_oak")),
    fill = factor(case_when(dist.trend_p.value_AA_tmean < 0.05 & dist.trend_p.value_AB_tmean < 0.05 ~ "filled",
                            dist.trend_p.value_AA_tmean >= 0.05 & dist.trend_p.value_AB_tmean >= 0.05 ~ "not_filled",
                            TRUE ~ NA_character_)),size = 3, stroke = 1.5)) +
  geom_hline(yintercept = 0, linetype = "dashed") +geom_vline(xintercept = 0, linetype = "dashed") +
  xlim(-2.7, 2.7) +ylim(-2.7, 2.7) +labs(x = "beta_B", y = "beta_A", title = "rc prec") +
  scale_shape_manual(values = c("filled_oak" = 16, "filled_beech" = 18, "hollow_oak" = 1, "hollow_beech" = 5,
                                "cross_oak" = 10, "cross_beech" = 9, "thick_outline_oak" = 21, "thick_outline_beech" = 23)) +
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "blue")) +
  scale_fill_manual(values = c("filled" = "red", "not_filled" = "grey")) + # Set "not_filled" to grey
  theme_minimal() +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),legend.position = "none",
        axis.text.x = element_text(size = 14),axis.text.y = element_text(size = 14),
        axis.ticks = element_line(color = "black", linewidth = 0.8),axis.ticks.length = unit(4, "mm"),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1.5)
  )




ggplot(rt, aes(x = dist.trend_AB_tmean, y = dist.trend_AA_tmean)) +
  geom_point(aes(color = coefs_m1_interaction_snp_dist_tmean_p.value < 0.05,shape = factor(case_when(
    dist.trend_p.value_AA_tmean < 0.05 & dist.trend_p.value_AB_tmean < 0.05 & species == "oak" ~ "filled_oak",
    dist.trend_p.value_AA_tmean < 0.05 & dist.trend_p.value_AB_tmean < 0.05 & species == "beech" ~ "filled_beech",
    dist.trend_p.value_AA_tmean < 0.05 & dist.trend_p.value_AB_tmean >= 0.05 & species == "oak" ~ "hollow_oak",
    dist.trend_p.value_AA_tmean < 0.05 & dist.trend_p.value_AB_tmean >= 0.05 & species == "beech" ~ "hollow_beech",
    dist.trend_p.value_AA_tmean >= 0.05 & dist.trend_p.value_AB_tmean < 0.05 & species == "oak" ~ "cross_oak",
    dist.trend_p.value_AA_tmean >= 0.05 & dist.trend_p.value_AB_tmean < 0.05 & species == "beech" ~ "cross_beech",
    dist.trend_p.value_AA_tmean >= 0.05 & dist.trend_p.value_AB_tmean >= 0.05 & species == "beech" ~ "thick_outline_beech",
    dist.trend_p.value_AA_tmean >= 0.05 & dist.trend_p.value_AB_tmean >= 0.05 & species == "oak" ~ "thick_outline_oak")),
    fill = factor(case_when(dist.trend_p.value_AA_tmean < 0.05 & dist.trend_p.value_AB_tmean < 0.05 ~ "filled",
                            dist.trend_p.value_AA_tmean >= 0.05 & dist.trend_p.value_AB_tmean >= 0.05 ~ "not_filled",
                            TRUE ~ NA_character_)),size = 3, stroke = 1.5)) +
  geom_hline(yintercept = 0, linetype = "dashed") +geom_vline(xintercept = 0, linetype = "dashed") +
  xlim(-2.7, 2.7) +ylim(-2.7, 2.7) +labs(x = "beta_B", y = "beta_A", title = "rt prec") +
  scale_shape_manual(values = c("filled_oak" = 16, "filled_beech" = 18, "hollow_oak" = 1, "hollow_beech" = 5,
                                "cross_oak" = 10, "cross_beech" = 9, "thick_outline_oak" = 21, "thick_outline_beech" = 23)) +
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "blue")) +
  scale_fill_manual(values = c("filled" = "red", "not_filled" = "grey")) + # Set "not_filled" to grey
  theme_minimal() +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),legend.position = "none",
        axis.text.x = element_text(size = 14),axis.text.y = element_text(size = 14),
        axis.ticks = element_line(color = "black", linewidth = 0.8),axis.ticks.length = unit(4, "mm"),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1.5)
  )


ggplot(rs, aes(x = dist.trend_AB_tmean, y = dist.trend_AA_tmean)) +
  geom_point(aes(
    color = coefs_m1_interaction_snp_dist_tmean_p.value < 0.05,
    shape = factor(case_when(
      dist.trend_p.value_AA_tmean < 0.05 & dist.trend_p.value_AB_tmean < 0.05 & species == "oak" ~ "filled_oak",
      dist.trend_p.value_AA_tmean < 0.05 & dist.trend_p.value_AB_tmean < 0.05 & species == "beech" ~ "filled_beech",
      dist.trend_p.value_AA_tmean < 0.05 & dist.trend_p.value_AB_tmean >= 0.05 & species == "oak" ~ "hollow_oak",
      dist.trend_p.value_AA_tmean < 0.05 & dist.trend_p.value_AB_tmean >= 0.05 & species == "beech" ~ "hollow_beech",
      dist.trend_p.value_AA_tmean >= 0.05 & dist.trend_p.value_AB_tmean < 0.05 & species == "oak" ~ "cross_oak",
      dist.trend_p.value_AA_tmean >= 0.05 & dist.trend_p.value_AB_tmean < 0.05 & species == "beech" ~ "cross_beech",
      dist.trend_p.value_AA_tmean >= 0.05 & dist.trend_p.value_AB_tmean >= 0.05 & species == "beech" ~ "thick_outline_beech",
      dist.trend_p.value_AA_tmean >= 0.05 & dist.trend_p.value_AB_tmean >= 0.05 & species == "oak" ~ "thick_outline_oak"
    )),
    fill = factor(case_when(
      dist.trend_p.value_AA_tmean < 0.05 & dist.trend_p.value_AB_tmean < 0.05 ~ "filled",
      dist.trend_p.value_AA_tmean >= 0.05 & dist.trend_p.value_AB_tmean >= 0.05 ~ "not_filled",
      TRUE ~ NA_character_
    ))
  ), size = 4.5, stroke = .7) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(x = "beta_B", y = "beta_A", title = "rs prec square") +
  scale_shape_manual(values = c(
    "filled_oak" = 16, "filled_beech" = 18, "hollow_oak" = 1, "hollow_beech" = 5,
    "cross_oak" = 10, "cross_beech" = 9, "thick_outline_oak" = 21, "thick_outline_beech" = 23
  )) +
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "blue")) +
  scale_fill_manual(values = c("filled" = "red", "not_filled" = "grey")) +
  theme_minimal() +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),legend.position = "none",
        axis.text.x = element_text(size = 14),axis.text.y = element_text(size = 14),
        axis.ticks = element_line(color = "black", linewidth = 0.8),axis.ticks.length = unit(4, "mm"),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1.5)
  ) +
  coord_fixed(ratio = 1, xlim = c(-2.9, 2.9), ylim = c(-2.9, 2.9), expand = FALSE)




ggplot(rt, aes(x = dist.trend_AB_tmean, y = dist.trend_AA_tmean)) +
  geom_point(aes(
    color = coefs_m1_interaction_snp_dist_tmean_p.value < 0.05,
    shape = factor(case_when(
      dist.trend_p.value_AA_tmean < 0.05 & dist.trend_p.value_AB_tmean < 0.05 & species == "oak" ~ "filled_oak",
      dist.trend_p.value_AA_tmean < 0.05 & dist.trend_p.value_AB_tmean < 0.05 & species == "beech" ~ "filled_beech",
      dist.trend_p.value_AA_tmean < 0.05 & dist.trend_p.value_AB_tmean >= 0.05 & species == "oak" ~ "hollow_oak",
      dist.trend_p.value_AA_tmean < 0.05 & dist.trend_p.value_AB_tmean >= 0.05 & species == "beech" ~ "hollow_beech",
      dist.trend_p.value_AA_tmean >= 0.05 & dist.trend_p.value_AB_tmean < 0.05 & species == "oak" ~ "cross_oak",
      dist.trend_p.value_AA_tmean >= 0.05 & dist.trend_p.value_AB_tmean < 0.05 & species == "beech" ~ "cross_beech",
      dist.trend_p.value_AA_tmean >= 0.05 & dist.trend_p.value_AB_tmean >= 0.05 & species == "beech" ~ "thick_outline_beech",
      dist.trend_p.value_AA_tmean >= 0.05 & dist.trend_p.value_AB_tmean >= 0.05 & species == "oak" ~ "thick_outline_oak"
    )),
    fill = factor(case_when(
      dist.trend_p.value_AA_tmean < 0.05 & dist.trend_p.value_AB_tmean < 0.05 ~ "filled",
      dist.trend_p.value_AA_tmean >= 0.05 & dist.trend_p.value_AB_tmean >= 0.05 ~ "not_filled",
      TRUE ~ NA_character_
    ))
  ), size = 4.5, stroke = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(x = "beta_B", y = "beta_A", title = "rt prec square") +
  scale_shape_manual(values = c(
    "filled_oak" = 16, "filled_beech" = 18, "hollow_oak" = 1, "hollow_beech" = 5,
    "cross_oak" = 10, "cross_beech" = 9, "thick_outline_oak" = 21, "thick_outline_beech" = 23
  )) +
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "blue")) +
  scale_fill_manual(values = c("filled" = "red", "not_filled" = "grey")) +
  theme_minimal() +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),legend.position = "none",
        axis.text.x = element_text(size = 14),axis.text.y = element_text(size = 14),
        axis.ticks = element_line(color = "black", linewidth = 0.8),axis.ticks.length = unit(4, "mm"),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1.5)
  ) +
  coord_fixed(ratio = 1, xlim = c(-2.9, 2.9), ylim = c(-2.9, 2.9), expand = FALSE)




ggplot(rc, aes(x = dist.trend_AB_tmean, y = dist.trend_AA_tmean)) +
  geom_point(aes(
    color = coefs_m1_interaction_snp_dist_tmean_p.value < 0.05,
    shape = factor(case_when(
      dist.trend_p.value_AA_tmean < 0.05 & dist.trend_p.value_AB_tmean < 0.05 & species == "oak" ~ "filled_oak",
      dist.trend_p.value_AA_tmean < 0.05 & dist.trend_p.value_AB_tmean < 0.05 & species == "beech" ~ "filled_beech",
      dist.trend_p.value_AA_tmean < 0.05 & dist.trend_p.value_AB_tmean >= 0.05 & species == "oak" ~ "hollow_oak",
      dist.trend_p.value_AA_tmean < 0.05 & dist.trend_p.value_AB_tmean >= 0.05 & species == "beech" ~ "hollow_beech",
      dist.trend_p.value_AA_tmean >= 0.05 & dist.trend_p.value_AB_tmean < 0.05 & species == "oak" ~ "cross_oak",
      dist.trend_p.value_AA_tmean >= 0.05 & dist.trend_p.value_AB_tmean < 0.05 & species == "beech" ~ "cross_beech",
      dist.trend_p.value_AA_tmean >= 0.05 & dist.trend_p.value_AB_tmean >= 0.05 & species == "beech" ~ "thick_outline_beech",
      dist.trend_p.value_AA_tmean >= 0.05 & dist.trend_p.value_AB_tmean >= 0.05 & species == "oak" ~ "thick_outline_oak"
    )),
    fill = factor(case_when(
      dist.trend_p.value_AA_tmean < 0.05 & dist.trend_p.value_AB_tmean < 0.05 ~ "filled",
      dist.trend_p.value_AA_tmean >= 0.05 & dist.trend_p.value_AB_tmean >= 0.05 ~ "not_filled",
      TRUE ~ NA_character_
    ))
  ), size = 4.5, stroke = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(x = "beta_B", y = "beta_A", title = "rc prec square") +
  scale_shape_manual(values = c(
    "filled_oak" = 16, "filled_beech" = 18, "hollow_oak" = 1, "hollow_beech" = 5,
    "cross_oak" = 10, "cross_beech" = 9, "thick_outline_oak" = 21, "thick_outline_beech" = 23
  )) +
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "blue")) +
  scale_fill_manual(values = c("filled" = "red", "not_filled" = "grey")) +
  theme_minimal() +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),legend.position = "none",
        axis.text.x = element_text(size = 14),axis.text.y = element_text(size = 14),
        axis.ticks = element_line(color = "black", linewidth = 0.8),axis.ticks.length = unit(4, "mm"),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1.5)
  ) +
  coord_fixed(ratio = 1, xlim = c(-2.9, 2.9), ylim = c(-2.9, 2.9), expand = FALSE)










ggplot(sprtmea, aes(x = dist.trend_AB_tmean, y = dist.trend_AA_tmean)) +
  geom_point(aes(color = coefs_m1_interaction_snp_dist_tmean_p.value < 0.05,shape = factor(case_when(
    dist.trend_p.value_AA_tmean < 0.05 & dist.trend_p.value_AB_tmean < 0.05 & species == "oak" ~ "filled_oak",
    dist.trend_p.value_AA_tmean < 0.05 & dist.trend_p.value_AB_tmean < 0.05 & species == "beech" ~ "filled_beech",
    dist.trend_p.value_AA_tmean < 0.05 & dist.trend_p.value_AB_tmean >= 0.05 & species == "oak" ~ "hollow_oak",
    dist.trend_p.value_AA_tmean < 0.05 & dist.trend_p.value_AB_tmean >= 0.05 & species == "beech" ~ "hollow_beech",
    dist.trend_p.value_AA_tmean >= 0.05 & dist.trend_p.value_AB_tmean < 0.05 & species == "oak" ~ "cross_oak",
    dist.trend_p.value_AA_tmean >= 0.05 & dist.trend_p.value_AB_tmean < 0.05 & species == "beech" ~ "cross_beech",
    dist.trend_p.value_AA_tmean >= 0.05 & dist.trend_p.value_AB_tmean >= 0.05 & species == "beech" ~ "thick_outline_beech",
    dist.trend_p.value_AA_tmean >= 0.05 & dist.trend_p.value_AB_tmean >= 0.05 & species == "oak" ~ "thick_outline_oak")),
    fill = factor(case_when(dist.trend_p.value_AA_tmean < 0.05 & dist.trend_p.value_AB_tmean < 0.05 ~ "filled",
                            dist.trend_p.value_AA_tmean >= 0.05 & dist.trend_p.value_AB_tmean >= 0.05 ~ "not_filled",
                            TRUE ~ NA_character_)),size = 3, stroke = 1.5)) +
  geom_hline(yintercept = 0, linetype = "dashed") +geom_vline(xintercept = 0, linetype = "dashed") +
  xlim(-2.5, 2.5) +ylim(-2.5, 2.5) +labs(x = "beta_B", y = "beta_A", title = "sprtmea prec") +
  scale_shape_manual(values = c("filled_oak" = 16, "filled_beech" = 18, "hollow_oak" = 1, "hollow_beech" = 5,
                                "cross_oak" = 10, "cross_beech" = 9, "thick_outline_oak" = 21, "thick_outline_beech" = 23)) +
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "blue")) +
  scale_fill_manual(values = c("filled" = "red", "not_filled" = "grey")) + # Set "not_filled" to grey
  theme_minimal() +theme(panel.grid.major = element_blank(),
                         panel.grid.minor = element_blank(),legend.position = "none",axis.text.x = element_text(size = 14),
                         axis.text.y = element_text(size = 14),panel.border = element_rect(color = "black", fill = NA, size = 1.5))


ggplot(sumtmea, aes(x = dist.trend_AB_tmean, y = dist.trend_AA_tmean)) +
  geom_point(aes(color = coefs_m1_interaction_snp_dist_tmean_p.value < 0.05,shape = factor(case_when(
    dist.trend_p.value_AA_tmean < 0.05 & dist.trend_p.value_AB_tmean < 0.05 & species == "oak" ~ "filled_oak",
    dist.trend_p.value_AA_tmean < 0.05 & dist.trend_p.value_AB_tmean < 0.05 & species == "beech" ~ "filled_beech",
    dist.trend_p.value_AA_tmean < 0.05 & dist.trend_p.value_AB_tmean >= 0.05 & species == "oak" ~ "hollow_oak",
    dist.trend_p.value_AA_tmean < 0.05 & dist.trend_p.value_AB_tmean >= 0.05 & species == "beech" ~ "hollow_beech",
    dist.trend_p.value_AA_tmean >= 0.05 & dist.trend_p.value_AB_tmean < 0.05 & species == "oak" ~ "cross_oak",
    dist.trend_p.value_AA_tmean >= 0.05 & dist.trend_p.value_AB_tmean < 0.05 & species == "beech" ~ "cross_beech",
    dist.trend_p.value_AA_tmean >= 0.05 & dist.trend_p.value_AB_tmean >= 0.05 & species == "beech" ~ "thick_outline_beech",
    dist.trend_p.value_AA_tmean >= 0.05 & dist.trend_p.value_AB_tmean >= 0.05 & species == "oak" ~ "thick_outline_oak")),
    fill = factor(case_when(dist.trend_p.value_AA_tmean < 0.05 & dist.trend_p.value_AB_tmean < 0.05 ~ "filled",
                            dist.trend_p.value_AA_tmean >= 0.05 & dist.trend_p.value_AB_tmean >= 0.05 ~ "not_filled",
                            TRUE ~ NA_character_)),size = 3, stroke = 1.5)) +
  geom_hline(yintercept = 0, linetype = "dashed") +geom_vline(xintercept = 0, linetype = "dashed") +
  xlim(-2.5, 2.5) +ylim(-2.5, 2.5) +labs(x = "beta_B", y = "beta_A", title = "sumtmea prec") +
  scale_shape_manual(values = c("filled_oak" = 16, "filled_beech" = 18, "hollow_oak" = 1, "hollow_beech" = 5,
                                "cross_oak" = 10, "cross_beech" = 9, "thick_outline_oak" = 21, "thick_outline_beech" = 23)) +
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "blue")) +
  scale_fill_manual(values = c("filled" = "red", "not_filled" = "grey")) + # Set "not_filled" to grey
  theme_minimal() +theme(panel.grid.major = element_blank(),
                         panel.grid.minor = element_blank(),legend.position = "none",axis.text.x = element_text(size = 14),
                         axis.text.y = element_text(size = 14),panel.border = element_rect(color = "black", fill = NA, size = 1.5))



ggplot(sprpre, aes(x = dist.trend_AB_tmean, y = dist.trend_AA_tmean)) +
  geom_point(aes(color = coefs_m1_interaction_snp_dist_tmean_p.value < 0.05,shape = factor(case_when(
    dist.trend_p.value_AA_tmean < 0.05 & dist.trend_p.value_AB_tmean < 0.05 & species == "oak" ~ "filled_oak",
    dist.trend_p.value_AA_tmean < 0.05 & dist.trend_p.value_AB_tmean < 0.05 & species == "beech" ~ "filled_beech",
    dist.trend_p.value_AA_tmean < 0.05 & dist.trend_p.value_AB_tmean >= 0.05 & species == "oak" ~ "hollow_oak",
    dist.trend_p.value_AA_tmean < 0.05 & dist.trend_p.value_AB_tmean >= 0.05 & species == "beech" ~ "hollow_beech",
    dist.trend_p.value_AA_tmean >= 0.05 & dist.trend_p.value_AB_tmean < 0.05 & species == "oak" ~ "cross_oak",
    dist.trend_p.value_AA_tmean >= 0.05 & dist.trend_p.value_AB_tmean < 0.05 & species == "beech" ~ "cross_beech",
    dist.trend_p.value_AA_tmean >= 0.05 & dist.trend_p.value_AB_tmean >= 0.05 & species == "beech" ~ "thick_outline_beech",
    dist.trend_p.value_AA_tmean >= 0.05 & dist.trend_p.value_AB_tmean >= 0.05 & species == "oak" ~ "thick_outline_oak")),
    fill = factor(case_when(dist.trend_p.value_AA_tmean < 0.05 & dist.trend_p.value_AB_tmean < 0.05 ~ "filled",
                            dist.trend_p.value_AA_tmean >= 0.05 & dist.trend_p.value_AB_tmean >= 0.05 ~ "not_filled",
                            TRUE ~ NA_character_)),size = 3, stroke = 1.5)) +
  geom_hline(yintercept = 0, linetype = "dashed") +geom_vline(xintercept = 0, linetype = "dashed") +
  xlim(-2.5, 2.5) +ylim(-2.5, 2.5) +labs(x = "beta_B", y = "beta_A", title = "sprpre prec") +
  scale_shape_manual(values = c("filled_oak" = 16, "filled_beech" = 18, "hollow_oak" = 1, "hollow_beech" = 5,
                                "cross_oak" = 10, "cross_beech" = 9, "thick_outline_oak" = 21, "thick_outline_beech" = 23)) +
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "blue")) +
  scale_fill_manual(values = c("filled" = "red", "not_filled" = "grey")) + # Set "not_filled" to grey
  theme_minimal() +theme(panel.grid.major = element_blank(),
                         panel.grid.minor = element_blank(),legend.position = "none",axis.text.x = element_text(size = 14),
                         axis.text.y = element_text(size = 14),panel.border = element_rect(color = "black", fill = NA, size = 1.5))


ggplot(sumpre, aes(x = dist.trend_AB_tmean, y = dist.trend_AA_tmean)) +
  geom_point(aes(color = coefs_m1_interaction_snp_dist_tmean_p.value < 0.05,shape = factor(case_when(
    dist.trend_p.value_AA_tmean < 0.05 & dist.trend_p.value_AB_tmean < 0.05 & species == "oak" ~ "filled_oak",
    dist.trend_p.value_AA_tmean < 0.05 & dist.trend_p.value_AB_tmean < 0.05 & species == "beech" ~ "filled_beech",
    dist.trend_p.value_AA_tmean < 0.05 & dist.trend_p.value_AB_tmean >= 0.05 & species == "oak" ~ "hollow_oak",
    dist.trend_p.value_AA_tmean < 0.05 & dist.trend_p.value_AB_tmean >= 0.05 & species == "beech" ~ "hollow_beech",
    dist.trend_p.value_AA_tmean >= 0.05 & dist.trend_p.value_AB_tmean < 0.05 & species == "oak" ~ "cross_oak",
    dist.trend_p.value_AA_tmean >= 0.05 & dist.trend_p.value_AB_tmean < 0.05 & species == "beech" ~ "cross_beech",
    dist.trend_p.value_AA_tmean >= 0.05 & dist.trend_p.value_AB_tmean >= 0.05 & species == "beech" ~ "thick_outline_beech",
    dist.trend_p.value_AA_tmean >= 0.05 & dist.trend_p.value_AB_tmean >= 0.05 & species == "oak" ~ "thick_outline_oak")),
    fill = factor(case_when(dist.trend_p.value_AA_tmean < 0.05 & dist.trend_p.value_AB_tmean < 0.05 ~ "filled",
                            dist.trend_p.value_AA_tmean >= 0.05 & dist.trend_p.value_AB_tmean >= 0.05 ~ "not_filled",
                            TRUE ~ NA_character_)),size = 3, stroke = 1.5)) +
  geom_hline(yintercept = 0, linetype = "dashed") +geom_vline(xintercept = 0, linetype = "dashed") +
  xlim(-2.5, 2.5) +ylim(-2.5, 2.5) +labs(x = "beta_B", y = "beta_A", title = "sumpre prec") +
  scale_shape_manual(values = c("filled_oak" = 16, "filled_beech" = 18, "hollow_oak" = 1, "hollow_beech" = 5,
                                "cross_oak" = 10, "cross_beech" = 9, "thick_outline_oak" = 21, "thick_outline_beech" = 23)) +
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "blue")) +
  scale_fill_manual(values = c("filled" = "red", "not_filled" = "grey")) + # Set "not_filled" to grey
  theme_minimal() +theme(panel.grid.major = element_blank(),
                         panel.grid.minor = element_blank(),legend.position = "none",axis.text.x = element_text(size = 14),
                         axis.text.y = element_text(size = 14),panel.border = element_rect(color = "black", fill = NA, size = 1.5))



