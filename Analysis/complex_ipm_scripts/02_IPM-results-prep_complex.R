###############################################################################X
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#################### Integrated Population Model (IPM): #######################X
#                 #---# PhD Dissertation: Chapter 1 #---#
#        Creating a Bayesian IPM to inform turkey management in PA
###                     *** Output organization ***                         ###X
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #X
# 
# Created by: Veronica A. Winter
# Last edited: 12/18/2024
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #X
###############################################################################X
# clean env
rm(list = ls())
gc()
###-----------------------------------------------------#X
# Set seed for reproducibility
set.seed(1235)
###-----------------------------------------------------#X
# Load necessary libraries
library(dplyr)
library(tidyr)
library(coda)
library(stringr)
library(purrr)

# source functions
source("Analysis/00_output-processing_funs.R")

# Load wmu area data
wmu_areas <- readRDS("Data/wmu_areas_km.rds")

###-----------------------------------------------------#X
# Organize data frame for IPM run with 4 WMUs for M, F, Ad, Juv

# Output includes:

# Recruitment:
# 'recruitment'

# Harvest rate:
# "female.h.ad.wmu", "female.h.juv.wmu",  
# "male.h.juv.wmu","male.h.ad.wmu"

# Survival: 
# "female.s.ad.wmu", "avg.juv.s.kf",
# "male.s.juv.wmu","male.s.ad.wmu"

# Abundance:
# "male.N.ad", "male.N.juv", 
# "female.N.ad", "female.N.juv"

###-----------------------------------------------------#X
# 12.28.2024 is the final version!
load("Data/Output/20241228_Complex_IPM_run.Rdata")

# `samples` is an MCMC array with dimensions [WMU, Year]
samples_df <- as.data.frame(ipm_run$chain1)

# Process each category using the function ----
# KF survival for..
# adult females
kf_female_surv_ad_df <- process_category_wmu(samples_df, "avg.ad.s.kf", "Female", "Adult")
# Juvenile females
kf_female_surv_juv_df <- process_category_wmu(samples_df, "avg.juv.s.kf", "Female", "Juvenile")

# DRM survival for..
# adult males
drm_male_surv_ad_df <- process_category(samples_df, "male.s.ad.wmu", "Male", "Adult")
# juvenile males
drm_male_surv_juv_df <- process_category(samples_df, "male.s.juv.wmu", "Male", "Juvenile")

# DRM harvest for..
# adult females
drm_female_harv_ad_df <- process_category(samples_df, "female.h.ad.wmu", "Female", "Adult")
# juv
drm_female_harv_juv_df <- process_category(samples_df, "female.h.juv.wmu", "Female", "Juvenile")
# adult males
drm_male_harv_ad_df <- process_category(samples_df, "male.h.ad.wmu", "Male", "Adult")
# juvenile males
drm_male_harv_juv_df <- process_category(samples_df, "male.h.juv.wmu", "Male", "Juvenile")

# PPB 
ppb <- process_category(samples_df, "aug31.ppb", "Female", "ppb")

ppb <- ppb %>% mutate(wmu = factor(wmu, levels = c(1, 2, 3, 4), 
                                     labels = c("2D", "3D", "4D", "5C")))
# HWB
hwb <- process_category(samples_df, "aug31.hwb", "Female", "hwb")
hwb <- hwb %>%  mutate(wmu = factor(wmu, levels = c(1, 2, 3, 4), 
                                    labels = c("2D", "3D", "4D", "5C")))
# recruitment
rec <- process_category(samples_df, "recruitment", "Female", "Adult")
rec <- rec %>%  mutate(wmu = factor(wmu, levels = c(1, 2, 3, 4), 
                                    labels = c("2D", "3D", "4D", "5C"))) %>% 
  left_join(wmu_areas, by = c("wmu" = "WMU_ID")) %>% 
  mutate(density_value = median_value / area_sq_km)

# Abundance
n_ad_female <- process_category(samples_df, "female.N.ad", "Female", "Adult")
n_juv_female <- process_category(samples_df, "female.N.juv", "Female", "Juvenile")
n_ad_male <- process_category(samples_df, "male.N.ad", "Male", "Adult")
n_juv_male <- process_category(samples_df, "male.N.juv", "Male", "Juvenile")

# Combine all data into one long format data frame ----
# Abundance
abundance_df <- bind_rows(
  n_ad_female, n_juv_female,
  n_ad_male, n_juv_male
) %>%
  mutate(wmu = factor(wmu, levels = c(1, 2, 3, 4), 
                      labels = c("2D", "3D", "4D", "5C"))) %>%
  left_join(wmu_areas, by = c("wmu" = "WMU_ID")) %>% 
  mutate(demographic_est = "Abundance",
         density_value = median_value / area_sq_km) #%>% 
  #filter(year != "5")

# DRM survival
drm_survival_df <- bind_rows(
  drm_male_surv_ad_df, drm_male_surv_juv_df
) %>% 
  mutate(wmu = factor(wmu, levels = c(1, 2, 3, 4), 
                      labels = c("2D", "3D", "4D", "5C")),
         demographic_est = "DRM_Survival")

# Known fate survival
kf_survival_df <- bind_rows(
  kf_female_surv_ad_df, kf_female_surv_juv_df
) %>% 
  mutate(wmu = factor(wmu, levels = c(1, 2, 3, 4),
                       labels = c("2D", "3D", "4D", "5C")),
          demographic_est = "KF_Survival")

# Combine kf and drm
# Step 1: Get the unique years in drm_survival_df
unique_years <- drm_survival_df %>% distinct(year) %>% pull(year)

# Step 2: Expand kf_survival_df to have one entry per year per WMU
expanded_kf_survival_df <- kf_survival_df %>%
  crossing(year = unique_years)

# Step 3: Bind the data frames together for plotting or analysis
combined_survival_df <- bind_rows(
  drm_survival_df,
  expanded_kf_survival_df
)


# DRM harvest rate 
drm_harvest_df <- bind_rows(
  drm_male_harv_ad_df, drm_male_harv_juv_df, 
   drm_female_harv_juv_df,
   drm_female_harv_ad_df
) %>% 
  mutate(wmu = factor(wmu, levels = c(1, 2, 3, 4), 
                      labels = c("2D", "3D", "4D", "5C")),
         demographic_est = "DRM_HarvestRate") 

# Save the  summary data frames to RDS files ----
date <- "20241228"
type = "Complex_"
saveRDS(kf_survival_df, paste0("Data/Output/", type, date,"_kf-survival_summary.rds"))
saveRDS(combined_survival_df, paste0("Data/Output/", type, date,"_comb-survival_summary.rds"))
saveRDS(drm_harvest_df, paste0("Data/Output/", type, date, "_harvest_summary.rds"))
saveRDS(abundance_df, paste0("Data/Output/", type, date, "_abundance_summary.rds"))
saveRDS(drm_survival_df, paste0("Data/Output/", type, date, "_drm_survival_summary.rds"))
saveRDS(ppb, paste0("Data/Output/", type, date, "_ppb_summary.rds"))
saveRDS(hwb, paste0("Data/Output/", type, date, "_hwb_summary.rds"))
saveRDS(rec, paste0("Data/Output/", type, date, "_rec_summary.rds"))
