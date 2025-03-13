###############################################################################X
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
################### Full Integrated Population Model (IPM): ###################X
#                 #---# PhD Dissertation: Complex IPM  #---#
#        Creating a Bayesian IPM to inform turkey management in PA
###                       *** Real data run ***                             ###X
###                        *** Model 1/2 ***                                ###X
###                                                                         ###X
#            Creating a Bayesian integrated population model (IPM) 
#     to estimate recruitment, survival, and abundance per wildlife management 
#                unit (WMU), and age class over seasons/years
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #X
###                                                                         ###X
#    Modeling: Known fate (kf) for hen annual survival, Hen with brood (HWB) and 
#       Pouts per brood (PPB) recruitment models, Dead-recovery model (DRM) for 
#       harvest rate and back-transformed annual survival for males, per WMU, 
#       and age class, a Lincoln-Peterson abundance estimator for annual 
#       abundance per WMU
###                                                                         ###X
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #X
# 
# Created by: Veronica A. Winter
# Last edited: 01/01/2024
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
library(RODBC)
library(nimble)
###-----------------------------------------------------#X
# functions
source("Analysis/Scripts/00_IPM_funs.R")
# Prepared Data
source("Analysis/Scripts/00_data-formating_IPM.R")

##################################################################X
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#              ######### THINGS TO NOTE: #########
#
# * Variable names correspond to the associated model. For example:
# **'hwb' = proportion of hens with a brood
# **'ppb' = poults-per-brood e.g., poult-to-hen ratio
# **'kf' = known-fate survival model for females
# **'drm' = dead-recovery model
# **'N' = abundance
# 
# **'ad' & 'juv' = adult/juvenile age class
# **'M' & 'F' = male/female sex
#
# * 'fall' = fall hunt & 'spring' = spring hunt
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
##################################################################X
#         ######### Recent changes: 12/05/2024 #########
#
# * Changed known fate model to have month included
# * Ensured WMU orders matched
# * Changed ppb & hwb intercept to 2019, not 2023
#
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
  telem.juvenile = is_juvenile_matrix, # adult is intercept
  telem.wmu = telem.wmu, # 2D is intercept term
  
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
  telem.nind = telem.nind,# 405
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
  male.nind = male.nind, # 1990
  male.n.occasions = male.n.occasions, # 5
  male.n.wmu = male.n.wmu, #10
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

##################################################X
## Call Nimble model ----
##################################################X
# source model
source("models/ipm.R")

###-----------------------------------------------------#X
# Create the Nimble model
model <- nimbleModel(
  code = ipm,
  data = nimble.data,
  constants = consts,
  inits = inits
)

# Set MCMC specifications
ni <- 600000
nt <- 4
nb <- 400000
nc <- 2 # minimum 2 chains for rhat values

# Run MCMC 
start <- Sys.time()
ipm_run <- nimbleMCMC(
  model = model,
  monitors = c(
    ###--------------------------- Recruitment ----------------------#X
    "recruitment", "aug31.ppb", "aug31.hwb",
    ###------------------------- Male DRM ----------------------------#X
    "male.h.ad.wmu", "male.s.ad.wmu", "male.h.juv.wmu", "male.s.juv.wmu",
    ###------------------------- Female Harvest ----------------------#X
    "female.h.ad.wmu", "female.h.juv.wmu",  
    ###------------------------- Male abundance ----------------------#X
    "male.N.ad", "male.N.juv",
    ###------------------------ Female abundance ----------------------#X
    "female.N.ad", "female.N.juv",  
    ###------------------------- Female survival ----------------------#X
    "avg.ad.s.kf", "avg.juv.s.kf", "storage"
  ),
  niter = ni, 
  nburnin = nb, 
  thin = nt, 
  nchains = nc, 
  samplesAsCodaMCMC = TRUE,
  setSeed = 1235,
  WAIC = F
)
end <- Sys.time()
print(paste("Total time", (end - start)))

# Save output
folder_path <- "Data/"
# Check if the folder exists, if not, create it
if (!dir.exists(folder_path)) {
  dir.create(folder_path)
  message("Folder created: ", folder_path)
} else {
  message("Folder already exists: ", folder_path)
}
saveRDS(ipm_run, "Data/Complex_IPM_run.rds")
save.image(file = "Data/Complex_IPM_run.Rdata")
#############################################################X
#           Model diagnostics ----
#############################################################X
summary(ipm_run)
# Get a summary of the samples
MCMCvis::MCMCsummary(ipm_run,  params = c("avg.ad.s.kf", "avg.juv.s.kf"))

# Visualize the trace plots
MCMCvis::MCMCtrace(ipm_run, pdf = F, params = c("female.N.ad"), iter = ni)

#############################################################X
#       Extract posterior values check ---
#############################################################X
MCMCvis::MCMCpstr(ipm_run$samples, params = "male.h.ad.wmu", func = median)
MCMCvis::MCMCpstr(ipm_run$samples, params = "avg.ad.s.kf", func = median)
MCMCvis::MCMCpstr(ipm_run$samples, params = "N.lambda.ad.male", func = median)
MCMCvis::MCMCpstr(ipm_run$samples, params = "recruitment", func = median)
MCMCvis::MCMCpstr(ipm_run$samples, params = "male.N.ad", func = median)
MCMCvis::MCMCpstr(ipm_run$samples, params = "female.N.ad", func = mean)
MCMCvis::MCMCpstr(ipm_run$samples, params = "male.N.juv", func = median)

