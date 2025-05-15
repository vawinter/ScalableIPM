###############################################################################X
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#################### Integrated Population Model (IPM): #######################X
#                 #---# PhD Dissertation: Simple IPM #---#
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

# # for saving outputs
# date <- "20241218"
# type <- "Simple_"

# source functions
source("Analysis/00_output-processing_funs.R")

# Load wmu area data
wmu_areas <- readRDS("Data/wmu_km_areas_regions.rds")

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
load("Data/Simple_IPM_run.Rdata")
###-----------------------------------------------------#X
## Vague prior on female survival ----X
# A Note on Vague priors model:
# This model is the same as the Simple IPM with the only change being that female
# harvest rates and survival probability are constructed with a vague prior of
# Beta(1,1) as opposed to an informative prior
#
# load("Data/Outputs/Simple_IPM_run-vagueprior.Rdata")
###-----------------------------------------------------#X
# `samples` is an MCMC array with dimensions [WMU, Year]
samples_df <- as.data.frame(ipm_run$chain1)

# Process each category using the function ----
# KF survival for..
# adult females
kf_female_surv_ad_df <- process_category(samples_df, "avg.ad.s.kf", "Female", "Adult")
# Juvenile females
kf_female_surv_juv_df <- process_category(samples_df, "avg.juv.s.kf", "Female", "Juvenile")

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

ppb <- ppb %>% mutate(Region = factor(wmu, levels = c(1, 2, 3, 4, 5, 6, 7, 8, 9), 
                                   labels = c("Region 1",
                                              "Region 2", 
                                              "Region 3", 
                                              "Region 4",
                                              "Region 5",
                                              "Region 6",
                                              "Region 7",
                                              "Region 8",
                                              "Region 9")))
# HWB
hwb <- process_category(samples_df, "aug31.hwb", "Female", "hwb")
hwb <- hwb %>%  mutate(Region = factor(wmu, levels = c(1, 2, 3, 4, 5, 6, 7, 8, 9), 
                                              labels = c("Region 1",
                                                         "Region 2", 
                                                         "Region 3", 
                                                         "Region 4",
                                                         "Region 5",
                                                         "Region 6",
                                                         "Region 7",
                                                         "Region 8",
                                                         "Region 9")))
# recruitment
rec <- process_category(samples_df, "recruitment", "Female", "Adult")
rec <- rec %>%   mutate(Region = factor(wmu, levels = c(1, 2, 3, 4, 5, 6, 7, 8, 9), 
                                       labels = c("Region 1",
                                                  "Region 2", 
                                                  "Region 3", 
                                                  "Region 4",
                                                  "Region 5",
                                                  "Region 6",
                                                  "Region 7",
                                                  "Region 8",
                                                  "Region 9"))) %>% 
  left_join(wmu_areas, by = c("Region" = "WMU_ID")) %>% 
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
  mutate(Region = factor(wmu, levels = c(1, 2, 3, 4, 5, 6, 7, 8, 9), 
                        labels = c("Region 1",
                                   "Region 2", 
                                   "Region 3", 
                                   "Region 4",
                                   "Region 5",
                                   "Region 6",
                                   "Region 7",
                                   "Region 8",
                                   "Region 9"))) %>%
  left_join(wmu_areas, by = c("Region" = "WMU_ID")) %>% 
  mutate(demographic_est = "Abundance",
         density_value = median_value / area_sq_km)

# DRM survival
drm_survival_df <- bind_rows(
  drm_male_surv_ad_df, drm_male_surv_juv_df
) %>% 
  mutate(Region = factor(wmu, levels = c(1, 2, 3, 4, 5, 6, 7, 8, 9), 
                        labels = c("Region 1",
                                   "Region 2", 
                                   "Region 3", 
                                   "Region 4",
                                   "Region 5",
                                   "Region 6",
                                   "Region 7",
                                   "Region 8",
                                   "Region 9")),
         demographic_est = "DRM_Survival")

# Known fate survival
kf_survival_df <- bind_rows(
  kf_female_surv_ad_df, kf_female_surv_juv_df
) %>% 
  mutate(Region = factor(wmu, levels = c(1, 2, 3, 4, 5, 6, 7, 8, 9), 
                        labels = c("Region 1",
                                   "Region 2", 
                                   "Region 3", 
                                   "Region 4",
                                   "Region 5",
                                   "Region 6",
                                   "Region 7",
                                   "Region 8",
                                   "Region 9")),
         demographic_est = "KF_Survival")

# DRM harvest rate 
drm_harvest_df <- bind_rows(
  drm_male_harv_ad_df, drm_male_harv_juv_df, 
  drm_female_harv_juv_df,
  drm_female_harv_ad_df
) %>% 
  mutate(Region = factor(wmu, levels = c(1, 2, 3, 4, 5, 6, 7, 8, 9), 
                        labels = c("Region 1",
                                   "Region 2", 
                                   "Region 3", 
                                   "Region 4",
                                   "Region 5",
                                   "Region 6",
                                   "Region 7",
                                   "Region 8",
                                   "Region 9")),
         demographic_est = "DRM_HarvestRate") 


# DRM survival
drm_survival_df <- bind_rows(
  drm_male_surv_ad_df, drm_male_surv_juv_df
) %>% 
  mutate(Region = factor(wmu, levels = c(1:9), 
                      labels = c("Region 1",
                                 "Region 2", 
                                 "Region 3", 
                                 "Region 4",
                                 "Region 5",
                                 "Region 6",
                                 "Region 7",
                                 "Region 8",
                                 "Region 9")),
         demographic_est = "DRM_Survival")

# Known fate survival
kf_survival_df <- bind_rows(
  kf_female_surv_ad_df, kf_female_surv_juv_df
) %>% 
  mutate(Region = factor(wmu, levels = c(1:9), 
                      labels = c("Region 1",
                                 "Region 2", 
                                 "Region 3", 
                                 "Region 4",
                                 "Region 5",
                                 "Region 6",
                                 "Region 7",
                                 "Region 8",
                                 "Region 9")),
         demographic_est = "KF_Survival")

# Step 3: Bind the data frames together for plotting or analysis
combined_survival_df <- bind_rows(
  drm_survival_df,
  kf_survival_df
)


# Save the  summary data frames to RDS files ----
folder_path <- "Data/Output/"
# Check if the folder exists, if not, create it
if (!dir.exists(folder_path)) {
  dir.create(folder_path)
  message("Folder created: ", folder_path)
} else {
  message("Folder already exists: ", folder_path)
}

#type <- "Operational_vague_"
type <- "Operational"
saveRDS(kf_survival_df, paste0("Data/Output/",  type, "_kf-survival_summary.rds"))
saveRDS(combined_survival_df, paste0("Data/Output/", type, "_comb-survival_summary.rds"))
saveRDS(drm_harvest_df, paste0("Data/Output/",  type,  "_harvest_summary.rds"))
saveRDS(abundance_df, paste0("Data/Output/",  type,  "_abundance_summary.rds"))
saveRDS(drm_survival_df, paste0("Data/Output/",  type,  "_drm_survival_summary.rds"))
saveRDS(ppb, paste0("Data/Output/", type,  "_ppb_summary.rds"))
saveRDS(hwb, paste0("Data/Output/",  type,  "_hwb_summary.rds"))
saveRDS(rec, paste0("Data/Output/",  type,  "_rec_summary.rds"))

# Done!