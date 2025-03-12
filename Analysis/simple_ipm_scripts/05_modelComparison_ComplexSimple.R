################################################X
#      Normalize WAIC values for comparison    #X
## ========================================= ###X
# = Comparing model fit for full and simple    #X
##  ============ model outputs ============== ##X
#                Created by: VAW               #X
#                Date: 10/03/2024              #X
################################################X

# Clean environment
rm(list=ls())
gc()

##################################################################X
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#              ######### THINGS TO NOTE: #########
#
# * WAIC: Watanabe-Akaike information criteria
#  **also known as : widely applicable information criterion**
# 
# The WAIC is a measure of model fit that combines both the likelihood of 
# the data and the complexity of the model through pWAIC), so it is the 
# primary metric you'd use for comparison.
# 
# * Normalization:
# To normalize the WAIC, divide the WAIC by the number of data points in 
# each model, making the metric comparable across models with different 
# data sizes.
#
# * Formula for Normalized WAIC:
# Normalized WAIC=  WAIC/Number of data points
#
# ** This will give you a per-data-point WAIC that accounts for the   
# difference in data size.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
##################################################################X

# Prepare the data
# Simple model WAIC output (model 1) ----
# Field "WAIC":
#   [1] 69570.72
# Field "lppd":
#   [1] -34512.6
# Field "pWAIC":
#   [1] 272.7594

# # 1004 run
# Field "WAIC":
#   [1] 69550.27
# Field "lppd":
#   [1] -34495.77
# Field "pWAIC":
#   [1] 279.369

# Simple model
## Source data
source("Analysis/Scripts/00_data-formating_IPM_simple.R")

# Recruitment (HWB)
n_HWB <- length(HWB)

# Recruitment (PPB)
n_PHratio <- length(PHratio)

# DRM Male Data
n_male_y <- length(male.y)

# Abundance Data (Harvest)
n_harvest_ad_fall <- prod(dim(harvest.ad.fall))  # total elements in matrix
n_harvest_juv_fall <- prod(dim(harvest.juv.fall))
n_harvest_ad_spring <- prod(dim(harvest.ad.spring))
n_harvest_juv_spring <- prod(dim(harvest.juv.spring))

# First Year Harvest
n_th_year1_female_ad <- length(as.integer(round(harvest.ad.fall[1,])))
n_th_year1_female_juv <- length(as.integer(round(harvest.juv.fall[1,])))
n_th_year1_male_ad <- length(as.integer(round(harvest.ad.spring[1,])))
n_th_year1_male_juv <- length(as.integer(round(harvest.juv.spring[1,])))

# Total number of data points
total_data_points_simple <- n_HWB + n_PHratio + n_male_y + 
  n_harvest_ad_fall + n_harvest_juv_fall + 
  n_harvest_ad_spring + n_harvest_juv_spring +
  n_th_year1_female_ad
# 42177

# Full model (model 2) ----
# Field "WAIC":
#   [1] 68089.61
# Field "lppd":
#   [1] -33955.67
# Field "pWAIC":
#   [1] 89.13228

# clean env
rm(list = ls())

## Source data
source("Analysis/Scripts/00_data-formating_IPM.R")

# Full model
# Telemetry
n_telem <- length(status_matrix)

# Recruitment (HWB)
n_HWB <- length(HWB)

# Recruitment (PPB)
n_PHratio <- length(PHratio)

# DRM Male Data
n_male_y <- length(male.y)

# Abundance Data (Harvest)
n_harvest_ad_fall <- prod(dim(harvest.ad.fall))  # total elements in matrix
n_harvest_juv_fall <- prod(dim(harvest.juv.fall))
n_harvest_ad_spring <- prod(dim(harvest.ad.spring))
n_harvest_juv_spring <- prod(dim(harvest.juv.spring))

# First Year Harvest
n_th_year1_female_ad <- length(as.integer(round(harvest.ad.fall[1,])))
n_th_year1_female_juv <- length(as.integer(round(harvest.juv.fall[1,])))
n_th_year1_male_ad <- length(as.integer(round(harvest.ad.spring[1,])))
n_th_year1_male_juv <- length(as.integer(round(harvest.juv.spring[1,])))

# Total number of data points
total_data_points <- n_telem + n_HWB + n_PHratio + n_male_y + 
  n_harvest_ad_fall + n_harvest_juv_fall + 
  n_harvest_ad_spring + n_harvest_juv_spring +
  n_th_year1_female_ad

# Normalizing WAIC ----
# Model 1 (Simple)
waic_model1 <- 69550.27
n_data_model1 <- 42177  # Calculate this from your data
normalized_waic_model1 <- waic_model1 / n_data_model1

# Model 2 (Full)
waic_model2 <- 68089.61
n_data_model2 <- 35119  # Calculate this from your data
normalized_waic_model2 <- waic_model2 / n_data_model2

# Compare normalized WAIC
cat("Normalized WAIC for Model 1:", normalized_waic_model1, "\n") # ***
# Normalized WAIC for Model 1: 1.614378 
cat("Normalized WAIC for Model 2:", normalized_waic_model2, "\n")
# 1.938825 

# Done


