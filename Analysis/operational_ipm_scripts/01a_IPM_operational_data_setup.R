# Formatting data and model together 
# 03/2025

# Clean environment
rm(list = ls())
gc()

# Load necessary libraries
library(RODBC)
library(nimble)
library(parallel)
library(doParallel)
library(foreach)
library(coda)

# Set seed for reproducibility
set.seed(1235)

# Source required scripts
source("Analysis/Scripts/00_IPM_funs.R")
source("Analysis/Scripts/00_data-formating_IPM_simple.R")
source("models/ipm_simple.R")

##################################################X
# Estimate parameters in Nimble ----
##################################################X
### Nimble model set up ----
# Create n-1 variable
Nyears = male.n.occasions-1

# #---------------------------X
# Data
nimble.data <- list(
  
  ###-----------#X
  # Recruitment (HWB)
  HWB = HWB,
  hwb.doy.scale = hwb.doy.scale,
  hwb.doy.2 = hwb.doy.2,
  # hwb.Year2019 = hwb.Year2019,
  hwb.Year2020 = hwb.Year2020, 
  hwb.Year2021 = hwb.Year2021,
  hwb.Year2022 = hwb.Year2022,
  hwb.Year2023 = hwb.Year2023,
  hwb.aug31 = hwb.aug31,
  hwb.aug31.2 = hwb.aug31.2,
  
  ###-----------#X
  # Recruitment (PPB)
  PHratio = PHratio,
  ph.doy.scale = ph.doy.scale,
  ph.doy.2 = ph.doy.2,
  # ph.Year2019 = ph.Year2019,
  ph.Year2020 = ph.Year2020, 
  ph.Year2021 = ph.Year2021,
  ph.Year2022 = ph.Year2022,
  ph.Year2023 = ph.Year2023,
  ppb.aug31 = ppb.aug31,
  ppb.aug31.2 = ppb.aug31.2,
  
  ###-----------#X
  # DRM Male data
  male.y = male.y, 
  male.time.param = male.time.param,
  
  ###-----------#X
  # # Abundance (2020 - 2023)
  harvest.ad.fall = round(harvest.ad.fall[2:5,]),
  harvest.juv.fall = round(harvest.juv.fall[2:5,]),
  harvest.ad.spring  =  round(harvest.ad.spring[2:5,]),
  harvest.juv.spring = round(harvest.juv.spring[2:5,]),
  # First year harvest (2020)
  th.year1.female.ad = as.integer(round(harvest.ad.fall[2,])),
  th.year1.female.juv = as.integer(round(harvest.juv.fall[2,])),
  th.year1.male.ad = as.integer(round(harvest.ad.spring[2,])),
  th.year1.male.juv = as.integer(round(harvest.juv.spring[2,]))
)

# Constants
consts <- list(
  ###-----------#X
  # Recruitment (HWB)
  hwb.N = hwb.N,
  hwb.J = hwb.J,
  hwb.wmu = hwb.wmu,
  
  ###-----------#X
  # Recruitment (PPB)
  ph.N = ph.N,
  ph.J = ph.J,
  ph.wmu = ph.wmu,
  
  ###-----------#X
  # DRM Male constants
  male.f = male.f,
  male.I = male.I, # c(6866, 1094) 0, 1=juvenile
  male.II = male.II, # c(714, 1276) 0, 1=non-reward band
  male.nind = male.nind, # 1990
  male.n.occasions = male.n.occasions, # 5
  male.n.wmu = male.n.wmu, #10
  male.wmu = as.numeric(as.factor(male.wmu)),
  
  ###-----------#X
  female.n.wmu = 9, 
  female.n.occasions = female.n.occasions,
  ###-----------#X
  # Number of years we are estimating abundance for
  Nyears = Nyears
)


# Initial values
inits <- list(
  ###-----------#X
  # Recruitment (HWB)
  hwb.beta1 = 0,
  hwb.beta2 = 0,
  hwb.beta3 = 0,
  hwb.beta4 = 0,
  hwb.beta5 = 0,
  hwb.beta6 = 0,
  hwb.beta7 = 0,
  hwb.sigma = 1,
  hwb.u = rep(0, length(unique(hwb.wmu))),
  
  ###-----------#X
  # Recruitment (PPB)
  ph.beta1 = 0,
  ph.beta2 = 0,
  ph.beta3 = 0,
  ph.beta4 = 0,
  ph.beta5 = 0,
  ph.beta6 = 0,
  ph.beta7 = 0,
  ph.sigma.u = 1,
  ph.u = rep(0, length(unique(ph.wmu))),
  
  ###-----------#X
  # DRM Male inits
  male.z = male.z,
  male.juvenile.effect = rnorm(1, 0, 0.5),
  male.time.effect = rnorm((Nyears), 0, 0.5), 
  male.seber.recov = runif(2 ,0, 1),
  male.rrate.a = rnorm(1, 0.87, 0.039), 
  male.rrate.j = rnorm(1, 0.71, 0.072), 
  male.sigma = runif(1, 0, 10), 
  
  ###-----------#X
  # N1S1
  # Males
  male.N.ad.Survived = matrix(15000 , nrow = Nyears, ncol = male.n.wmu),
  male.N.juv.Survived = matrix(3000, nrow = Nyears, ncol = male.n.wmu),
  # Females
  female.N.ad.Survived = array(10000, dim = c(Nyears, 9)),
  female.N.juv.Survived = array(3000, dim = c(Nyears, 9)),
  ###-----------#X
  # Abundance inits
  # Set high initial N
  # suuuuuper sensitive to initial N
  # https://groups.google.com/g/nimble-users/c/fKx8OD7Qn9s
  # Males
  male.N.ad = matrix(100000 , nrow = Nyears, ncol = male.n.wmu),
  male.N.juv = matrix(90000, nrow = Nyears, ncol = male.n.wmu),
  # Females 
  female.N.ad = array(100000, dim = c(Nyears, 9)),
  female.N.juv = array(100000, dim = c(Nyears, 9))
)
save.image(file = "Data/Simple_IPM_setup-data/Simple_IPM_Nimble_data_setup.RData")
