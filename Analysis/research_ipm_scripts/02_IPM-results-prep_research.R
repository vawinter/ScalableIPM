###############################################################################X
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
############## Research Integrated Population Model (R_IPM): ##################X
#                     #---# PhD Dissertation: R_IPM #---#
#        Creating a Bayesian IPM to inform turkey management in PA
###                     *** Output organization ***                         ###X
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
wmu_areas <- readRDS("Data/wmu_areas_km.rds") %>% 
  mutate(WMU_ID = paste("WMU", WMU_ID, sep = " "))

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
#load("Data/Output/20250625_Parallel_O_vague_IPM_run.rdsR_IPM_run.Rdata")
combined_results <- readRDS("Data/Output/20251108_R_IPM_run23.rds")

# `samples` is an MCMC array with dimensions [WMU, Year]
samples_df <- as.data.frame(combined_results[[1]])

# Process each category using the function ----
# KF survival for..
# adult females
kf_female_surv_ad_df <- process_category_wmu(samples_df, "avg.ad.s.kf", "Female", "Adult")
# Juvenile females
kf_female_surv_juv_df <- process_category_wmu(samples_df, "avg.juv.s.kf", "Female", "Juvenile")
# Juvenile Males
kf_male_surv_juv_df <- process_category_wmu(samples_df, "juv.male.adj", "Male", "Juvenile")

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
                                     labels = c("WMU 2D", "WMU 3D", "WMU 4D", "WMU 5C")))
# HWB
hwb <- process_category(samples_df, "aug31.hwb", "Female", "hwb")
hwb <- hwb %>%  mutate(wmu = factor(wmu, levels = c(1, 2, 3, 4), 
                                    labels = c("WMU 2D", "WMU 3D", "WMU 4D", "WMU 5C")))
# recruitment
rec <- process_category(samples_df, "recruitment", "Female", "Adult")
rec <- rec %>%  mutate(wmu = factor(wmu, levels = c(1, 2, 3, 4), 
                                    labels = c("WMU 2D", "WMU 3D", "WMU 4D", "WMU 5C"))) %>% 
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
                      labels = c("WMU 2D", "WMU 3D", "WMU 4D", "WMU 5C"))) %>%
  left_join(wmu_areas, by = c("wmu" = "WMU_ID")) %>% 
  mutate(demographic_est = "Abundance",
         density_value = median_value / area_sq_km) 

# DRM survival
drm_survival_df <- bind_rows(
  drm_male_surv_ad_df, drm_male_surv_juv_df
) %>% 
  mutate(wmu = factor(wmu, levels = c(1, 2, 3), 
                      labels = c("WMU 2D", "WMU 3D", "WMU 4D")),
         demographic_est = "DRM_Survival")

# Known fate survival
kf_survival_df <- bind_rows(
  kf_female_surv_ad_df, kf_female_surv_juv_df, kf_male_surv_juv_df
) %>% 
  mutate(wmu = factor(wmu, levels = c(1, 2, 3, 4),
                       labels = c("WMU 2D", "WMU 3D", "WMU 4D", "WMU 5C")),
          demographic_est = "KF_Survival")

# Combine kf and drm
# Step 1: Get the unique years in drm_survival_df
unique_years <- drm_survival_df %>% distinct(year) %>% pull(year)

# Step 2: Expand kf_survival_df to have one entry per year per WMU
expanded_kf_survival_df <- kf_survival_df %>%
  filter(sex != "Male") %>% # dont need juvenile males in here for next df
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
  mutate(wmu = factor(wmu, levels = c(1, 2, 3), 
                      labels = c("WMU 2D", "WMU 3D", "WMU 4D")),
         demographic_est = "DRM_HarvestRate") 

# Save the  summary data frames to RDS files ----
folder_path <- "Data/Output/TEST/"
# Check if the folder exists, if not, create it
if (!dir.exists(folder_path)) {
  dir.create(folder_path)
  message("Folder created: ", folder_path)
} else {
  message("Folder already exists: ", folder_path)
}
type = "TEST/20251108_R_IPM_run23"
saveRDS(kf_survival_df, paste0("Data/Output/", type, "_kf-survival_summary.rds"))
saveRDS(combined_survival_df, paste0("Data/Output/", type, "_comb-survival_summary.rds"))
saveRDS(drm_harvest_df, paste0("Data/Output/", type,  "_harvest_summary.rds"))
saveRDS(abundance_df, paste0("Data/Output/", type,  "_abundance_summary.rds"))
saveRDS(drm_survival_df, paste0("Data/Output/", type,  "_drm_survival_summary.rds"))
saveRDS(ppb, paste0("Data/Output/", type,  "_ppb_summary.rds"))
saveRDS(hwb, paste0("Data/Output/", type,  "_hwb_summary.rds"))
saveRDS(rec, paste0("Data/Output/", type, "_rec_summary.rds"))
