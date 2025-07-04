###############################################################################X
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#################### Integrated Population Model (IPM): #######################X
#                        #---# Research IPM #---#
#        Creating a Bayesian IPM to inform turkey management in PA
###               ***  Data formatting and model fitting ***                ###X
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #X
###############################################################################X

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

# Import required scripts
source("Analysis/00_IPM_funs.R")
load("Data/Research_IPM_setup-data/Research_IPM_Nimble_data_setup.RData")
is.juv <- readRDS("Data/Research_IPM_setup-data/telem.juv.corrected.rds")
source("Models/research_ipm.R")

##################################################X
# Estimate parameters in Nimble ----
##################################################X
### Nimble model set up ----
# Create n-1 variable for years I am estimating abundance
Nyears = male.n.occasions-1

##---------------------------X
# Data
nimble.data <- list(
  # KF telemetered data 
  status = status_matrix,
  telem.juvenile = is.juv, 
  telem.wmu = telem.wmu, 
  
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
  #ph.Year2019 = ph.Year2019, 
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
  # Abundance (2020 - 2023)
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
  # KF telemetered constants 
  telem.nind = telem.nind,
  telem.first = telem.first,
  telem.last = telem.last,
  telem.year.start = telem.year.start,
  telem.year.end = telem.year.end,
  female.telem.wmu = 4,
  
  ###-----------#X
  # Recruitment (HWB)
  hwb.N = hwb.N,
  hwb.J = hwb.J,
  hwb.wmu = hwb.wmu,
  
  ###-----------#X
  # Recruitment (PPB)
  ph.N = ph.N,
  ph.J = ph.J, # make sure this is
  ph.wmu = ph.wmu,
  
  ###-----------#X
  # DRM Male constants
  male.f = male.f,
  male.I = male.I, # c(6866, 1094) 0, 1=juvenile
  male.II = male.II, # c(714, 1276) 0, 1=non-reward band
  male.nind = male.nind, 
  male.n.occasions = male.n.occasions, # 5
  male.n.wmu = male.n.wmu, 
  male.wmu = as.numeric(as.factor(male.wmu)),
  
  ###-----------#X
  female.n.wmu = female.n.wmu, # 4
  female.n.occasions = female.n.occasions,
  ###-----------#X
  # Number of years we are estimating abundance for
  Nyears = Nyears
)

# Initial values
inits <- list(
  # KF telemetered inits
  # Intercepts and coefficients
  telem.beta.int = rnorm(1, 0, 0.5),       
  telem.beta.age = rnorm(1, 0, 0.5),       
  telem.beta.wmu = rnorm(4, 0, 1),        
  telem.beta.month = rnorm(12, 0, 1),      
  
  # Variance parameters
  telem.sigma = runif(1, 0.1, 1),          
  telem.month.sigma = runif(1, 0, 2),
  # Survival probabilities 
  s.kf = array(runif(1, 0.2, 1),           
               dim = c(telem.nind, 4, 12)),
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
  # Males
  male.N.ad.Survived = matrix(3500 , nrow = Nyears, ncol = male.n.wmu),
  male.N.juv.Survived = matrix(1000, nrow = Nyears, ncol = male.n.wmu),
  # Females
  female.N.ad.Survived = array(300, dim = c(Nyears, female.n.wmu)),
  female.N.juv.Survived = array(300, dim = c(Nyears, female.n.wmu)),
  # Abundance inits
  # Set high initial N
  # suuuuuper sensitive to initial N
  # https://groups.google.com/g/nimble-users/c/fKx8OD7Qn9s
  # Males
  male.N.ad = matrix(80000 , nrow = Nyears, ncol = male.n.wmu),
  male.N.juv = matrix(90000, nrow = Nyears, ncol = male.n.wmu),
  # Females 
  female.N.ad = array(70000, dim = c(Nyears, female.n.wmu)),
  female.N.juv = array(70000, dim = c(Nyears, female.n.wmu))
)
save.image(file = "Data/Research_IPM_setup-data/R_IPM_Nimble_data_setup.RData")
