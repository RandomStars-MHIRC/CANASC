# Example of application of the lambda evaluation
# Corresponding to Figure 2

# Torin Halvorson
# BC Children's Hospital Research Institute
# University of British Columbia

# Dependencies: SI_general_vectiorial and Lambda_evaluation scripts
# and following libraries
# 20250714

library(tidyverse)
library(dplyr)
library(glue)
library(ggpubr)

#-------------------------------------------------------------------------
#define functions for PREVENT analysis

#add a small epsilon to avoid zero
epsilon <- function(x) {x + 0.001}

#function to convert zero/negative values to 0.001
zero_to_nonzero <- function(x) {
  if_else (x <= 0, 0.001, x)
}

#an alternate function that makes values <0 equal to 1/2 of the minimum of the remaining values
pos_vec <- function(vec) {
  if_else(vec > 0, vec, NA)
}
zero_to_halfmin <- function(x) {
  if_else (x <= 0, (min(pos_vec(x), na.rm = TRUE))/2, x)
  }

#-------------------------------------------------------------------------------
#PREVENT data

PREVENT_data <- read.csv("20250328_PREVENT_RawData.csv")
PREVENT_data <- PREVENT_data %>%
  mutate(across(CD4_CD134pCD25p:CD8_CD137pCD69p, epsilon))

vars <- colnames(PREVENT_data[4:9])

#var <- "CD4_CD134pCD25p"

for (var in vars) {
  df <- PREVENT_data %>%
    dplyr::select(DonorID,Timepoint, Stim, .data[[var]]) %>%
    filter(Timepoint == "V3") %>%
    pivot_wider(names_from = Stim, values_from = .data[[var]]) %>%
    mutate(Difference = COVID_WT - DMSO,
           Ratio = COVID_WT/DMSO,
           BoxCox_L0.5 = SI(x1 = COVID_WT, x0 = DMSO)) %>%
    mutate(across(Difference:BoxCox_L0.5, zero_to_nonzero)) %>%
    mutate(across(COVID_WT:BoxCox_L0.5, log10))
  write.csv(df, glue("{var}_Log10_PREVENTData_WithBoxCox.csv"))
  
  corr_plot_Raw <- df %>%
    ggplot(aes(x = DMSO, y = COVID_WT)) +
    geom_point(size = 2, shape = 22, colour = "blue4", fill = "blue4") +
    theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    stat_cor(aes(x = DMSO, y = COVID_WT), method = "spearman") +
    stat_smooth(method = "lm", 
                formula = y ~ x, 
                geom = "smooth", color = "blue4") +
    labs(title = glue("{var}"), x = "log10(Unstimulated)", y = "log10(SARS-CoV-2 Spike)")
  ggsave(glue("{var}_PREVENT_Raw_CorrPlot.pdf"))
  corr_plot_Raw
  
  corr_plot_Subtract <- df %>%
    ggplot(aes(x = DMSO, y = Difference)) +
    geom_point(size = 2, shape = 22, colour = "blue3", fill = "blue3") +
    theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    stat_cor(aes(x = DMSO, y = Difference), method = "spearman") +
    stat_smooth(method = "lm", 
                formula = y ~ x, 
                geom = "smooth", color = "blue3") +
    labs(title = glue("{var}"), x = "log10(Unstimulated)", y = "log10(SARS-CoV-2 Spike - Unstimulated)")
  ggsave(glue("{var}_PREVENT_Subtract_CorrPlot.pdf"))
  corr_plot_Subtract
  
  corr_plot_Divide <- df %>%
    ggplot(aes(x = DMSO, y = Ratio)) +
    geom_point(size = 2, shape = 22, colour = "steelblue", fill = "steelblue") +
    theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    stat_cor(aes(x = DMSO, y = Ratio), method = "spearman") +
    stat_smooth(method = "lm", 
                formula = y ~ x, 
                geom = "smooth", color = "steelblue") +
    labs(title = glue("{var}"), x = "log10(Unstimulated)", y = "log10(SARS-CoV-2 Spike / Unstimulated)")
  ggsave(glue("{var}_PREVENT_Divide_CorrPlot.pdf"))
  corr_plot_Divide
  
  corr_plot_Bcx <- df %>%
    ggplot(aes(x = DMSO, y = BoxCox_L0.5)) +
    geom_point(size = 2, shape = 22, colour = "green4", fill = "green4") +
    theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    stat_cor(aes(x = DMSO, y = BoxCox_L0.5), method = "spearman") +
    stat_smooth(method = "lm", 
                formula = y ~ x, 
                geom = "smooth", color = "green4") +
    labs(title = glue("{var}"), x = "log10(Unstimulated)", y = "log10(Box-Cox(SARS-CoV-2 Spike) - Box-Cox(Unstimulated))")
  ggsave(glue("{var}_PREVENT_Bcx_CorrPlot.pdf"))
  corr_plot_Bcx
}

for (var in vars) {
  PREVENT_data_V3 <- PREVENT_data %>%
    filter(Timepoint == "V3")
  
  PREVENT_data_V3_stim <- PREVENT_data_V3 %>%
    filter(Stim == "COVID_WT")
  PREVENT_data_V3_unstim <- PREVENT_data_V3 %>%
    filter(Stim == "DMSO")
  
  x1 <- epsilon(PREVENT_data_V3_stim[[var]])
  x0 <- epsilon(PREVENT_data_V3_unstim[[var]])
  
  SI05 <- SI(x1, x0, L = 0.5, corrected = FALSE)
  
  # Figure 2D if we use prevent
  pdf(glue("Lambda_tests_Prevent_{var}_20250502.pdf"))
  sink(file=glue("Lambda_tests_Prevent_{var}_20250502.txt"))
  res <- test.lambda.lik(x1, x0)
  write.csv(res, "Lambda_tests_Prevent_CD25pCD134p_CD4p_20250502.csv")
  sink()
  dev.off()
  
  pdf(glue("Example_Fig2D_{var}_20250502.pdf"))
  plotres3(res)
  plot(log10(x0), log10(SI05), pch=15, xlab="Unstim (log10)", ylab="SI(L=0.5) (log10)")
  M1 <- lm(log10(SI05)~log10(x0))
  abline(M1$coef[1], M1$coef[2])
  tmp <- cor.test(SI05, x0, method="spearman")
  title(sub=paste("rho=", round(tmp$estimate,2), "  p=",  round(tmp$p.value,4)))
  dev.off()

}

#------------------------------------------------------------------------------------------------------
#CAN-ASC AIM Data
CANASC_data <- read.csv("20250509_CANASC_RingTestingRawData.csv") 
CANASC_data <- CANASC_data %>%
  dplyr::select(SampleID, Stimulant, Experiment, Site, Analysis, CD25pCD134p_CD4p:CD134pCD137p_CD4p, CD69pCD137p_CD8p:CD107apCD137p_CD8p) %>%
  filter(Stimulant == "CMV" | Stimulant == "unstimulated") %>%
  filter(Analysis == "Central") %>%
  mutate(across(CD25pCD134p_CD4p:CD107apCD137p_CD8p, epsilon))

vars <- colnames(CANASC_data[6:12])

#var <- "CD25pCD134p_CD4p"

for (var in vars) {
  df <- CANASC_data %>%
    dplyr::select(SampleID, Stimulant, Experiment, Site, Analysis, .data[[var]]) %>%
    pivot_wider(names_from = Stimulant, values_from = .data[[var]]) %>%
    filter(!is.na(unstimulated) & !is.na(CMV)) %>%
    mutate(Difference = CMV - unstimulated,
           Ratio = CMV/unstimulated,
           BoxCox_L0.5 = SI(x1 = CMV, x0 = unstimulated)) %>%
    mutate(across(Difference:BoxCox_L0.5, zero_to_nonzero)) %>%
    mutate(across(unstimulated:BoxCox_L0.5, log10))
  write.csv(df, glue("{var}_Log10_CANASCData_WithBoxCox.csv"))
  
  #get residuals for raw data after removing SampleID effect
  Mod4 <- lm(CMV~SampleID, data = df)
  resid4 <- resid(Mod4)
  cor.test(df$unstimulated, resid4, method = "s")
  summary(lm(resid4~unstimulated:Site, data = df))
  
  # get residuals for subtraction (Difference) after removing SampleID effect
  Mod1 <- lm(Difference~SampleID, data = df)
  resid1 <- resid(Mod1)
  cor.test(df$unstimulated, resid1, method="s")
  summary(lm(resid1~unstimulated:Site, data = df))
  
  # get residuals for ratio after removing SampleID effect
  Mod0 <- lm(Ratio~SampleID, data = df)
  resid0 <- resid(Mod0)
  cor.test(df$unstimulated, resid0, method="s")
  summary(lm(resid0~unstimulated:Site, data = df))
  
  # get residuals for SI05 after removing DonorID effect
  Mod3 <- lm(BoxCox_L0.5~SampleID, data = df)
  resid3 <- resid(Mod3)
  cor.test(df$unstimulated, resid3, method="s")
  summary(lm(resid3~unstimulated:Site, data = df))
  
  #make plots
  corr_plot_Raw <- df %>%
    ggplot(aes(x = unstimulated, y = resid4)) +
    geom_point(size = 2, shape = 22, colour = "blue4", fill = "blue4") +
    theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    stat_cor(aes(x = unstimulated, y = resid4), method = "spearman") +
    stat_smooth(method = "lm", 
                formula = y ~ x, 
                geom = "smooth", color = "blue4") +
    labs(title = glue("{var}"), x = "log10(Unstimulated)", y = "log10(CMV pp65)")
  ggsave(glue("{var}_CANASC_Raw_CorrPlot.pdf"))
  corr_plot_Raw
  
  corr_plot_Subtract <- df %>%
    ggplot(aes(x = unstimulated, y = resid1)) +
    geom_point(size = 2, shape = 22, colour = "blue3", fill = "blue3") +
    theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    stat_cor(aes(x = unstimulated, y = resid1), method = "spearman") +
    stat_smooth(method = "lm", 
                formula = y ~ x, 
                geom = "smooth", color = "blue3") +
    labs(title = glue("{var}"), x = "log10(Unstimulated)", y = "log10(CMV pp65 - Unstimulated)")
  ggsave(glue("{var}_CANASC_Subtract_CorrPlot.pdf"))
  corr_plot_Subtract
  
  corr_plot_Divide <- df %>%
    ggplot(aes(x = unstimulated, y = resid0)) +
    geom_point(size = 2, shape = 22, colour = "steelblue", fill = "steelblue") +
    theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    stat_cor(aes(x = unstimulated, y = resid0), method = "spearman") +
    stat_smooth(method = "lm", 
                formula = y ~ x, 
                geom = "smooth", color = "steelblue") +
    labs(title = glue("{var}"), x = "log10(Unstimulated)", y = "log10(CMV pp65 / Unstimulated)")
  ggsave(glue("{var}_CANASC_Divide_CorrPlot.pdf"))
  corr_plot_Divide
  
  corr_plot_Bcx <- df %>%
    ggplot(aes(x = unstimulated, y = resid3)) +
    geom_point(size = 2, shape = 22, colour = "green4", fill = "green4") +
    theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    stat_cor(aes(x = unstimulated, y = resid3), method = "spearman") +
    stat_smooth(method = "lm", 
                formula = y ~ x, 
                geom = "smooth", color = "green4") +
    labs(title = glue("{var}"), x = "log10(Unstimulated)", y = "log10(Box-Cox(CMV pp65) - Box-Cox(Unstimulated))")
  ggsave(glue("{var}_CANASC_Bcx_CorrPlot.pdf"))
  corr_plot_Bcx
  
}

for (var in vars) {
  DAT5 <- CANASC_data %>%
    dplyr::select(SampleID, Stimulant, Experiment, Site, Analysis, .data[[var]]) %>%
    pivot_wider(names_from = Stimulant, values_from = .data[[var]]) %>%
    filter(!is.na(unstimulated) & !is.na(CMV)) %>%
    pivot_longer(cols = c(unstimulated, CMV), names_to = "Stimulant", values_to = as.character(glue("{var}")))
  write.csv(DAT5, glue("{var}_CANASCData_WithBoxCox.csv"))
  
  
  DAT5.Stim <- DAT5 %>%
    filter(Stimulant == "CMV")
  DAT5.Unstim <- DAT5 %>%
    filter(Stimulant == "unstimulated")
  
  x1 <- epsilon(DAT5.Stim[[var]])
  x0 <- epsilon(DAT5.Unstim[[var]])
  
  SI05 <- SI(x1, x0, L = 0.5, corrected = FALSE)
  
  tmp.col <- rep("gray50", nrow(DAT5.Stim))
  tmp.col[DAT5.Stim$Site == "YHM" ] <- "gray50"
  tmp.col[DAT5.Stim$Site == "YUL" ] <- "dodgerblue"
  tmp.col[DAT5.Stim$Site == "YVR" ] <- "red"
  
  tmp.pch <- as.numeric(as.factor(DAT5.Stim$SampleID))
  
  # get residuals for SI05 after removing DonorID effect
  Mod3 <- lm(log10(SI05)~SampleID, data=DAT5.Stim)
  resid3 <- resid(Mod3)
  #plot(log10(x0[keep]), resid3, pch=tmp.pch[keep], col=tmp.col[keep])
  cor.test(x0, resid3, method="s")
  summary(lm(resid3~log10(x0):Site, data=DAT5.Stim))
  
  #optimal lambda calculations
  #This version is corrected for DonorID only
  pdf(glue("Lambda_tests_RingTesting_{var}_20250502.pdf"))
  sink(file=glue("Lambda_tests_RingTesting_{var}_20250502.txt"))
  res <- test.lambda.lik(x1,x0,f1="SampleID", DAT=DAT5.Stim)
  sink()
  dev.off()
  write.csv(res, glue("Lambda_tests_RingTesting_{var}_20250502.xlsx"))
  
  pdf(glue("Example_Fig2E_Ring_{var}_20250502.pdf"))
  plotres3(res)
  plot(log10(x0), resid3, pch=15, xlab="Unstim (log10)", ylab="SI(L=0.5) (log10)")
  M1 <- lm(resid3~log10(x0))
  abline(M1$coef[1], M1$coef[2])
  tmp <- cor.test(resid3, x0, method="spearman")
  title(sub=paste("rho=", round(tmp$estimate,2), "  p=",  round(tmp$p.value,4)))
  dev.off()
  
  
}