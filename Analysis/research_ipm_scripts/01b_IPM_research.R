###############################################################################X
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
############## Research Integrated Population Model (R IPM): ##################X
###                                                                         ###X
#    Modeling: Known fate (kf) for hen annual survival, Hen with brood (HWB) and 
#       Pouts per brood (PPB) recruitment models, Dead-recovery model (DRM) for 
#       harvest rate and back-transformed annual survival for males, per WMU, 
#       and age class, a Lincoln-Peterson abundance estimator for annual 
#       abundance per WMU
###                                                                         ###X
###############################################################################X
# Parallelized Integrated Population Model (IPM) Script
# Based on original script by Veronica A. Winter
# Modified for parallel processing
# https://r-nimble.org/nimbleExamples/parallelizing_NIMBLE.html
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

# Set seed for reproducibility
set.seed(1235)
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

### Nimble model set up ----
# load("Data/Research_IPM_setup-data/R_IPM_Nimble_data_setup.RData") # original run
#load("Data/Research_IPM_setup-data/R_IPM_Nimble_data_setup_noKF24.RData") # 2024 no kf 24
#load("Data/Research_IPM_setup-data/R_IPM_Nimble_data_setup_24.RData") # all 2024
#load("Data/Research_IPM_setup-data/R_IPM_Nimble_data_setup_23updated22-24KF.RData") # 2023 Data with 22-24 kf from updated db
load("Data/Research_IPM_setup-data/TEST/R_IPM_Nimble_data_setup_23updated22-23KF_.RData")
#load("Data/Research_IPM_setup-data/TEST/R_IPM_24_kf23.RData") # all 2024 w. original kf23

source("Models/research_ipm_23_noabun.R")

# 10/29/2025 - on V:
# Running no abun model for 20-23 data with 'new' 22-23 KF data from most up to data '10-10' db
# 10/29/2025 - on Karen:
# Running no abun model for 20-23 data with 'new' 22-24* KF data from most up to data '10-10' db
#############################################################X
# Parallelization setup -----
#############################################################X
# Detect number of cores, leave one free for system processes
num_cores <- max(1, detectCores() - 1)

# Register parallel backend
cl <- makeCluster(num_cores)
registerDoParallel(cl)

# Export necessary objects to cluster
clusterExport(cl, c("nimble.data", "consts", "inits", "ipm"))

# Parallel function to run MCMC
run_parallel_mcmc <- function(seed_offset) {
  library(nimble)
  
  # Set unique seed for each chain
  set.seed(1235 + seed_offset)
  
  # Create the Nimble model
  model <- nimbleModel(
    code = ipm,
    data = nimble.data,
    constants = consts,
    inits = inits
  )
  
  # MCMC specifications
  ni <- 1200000  # Total iterations
  nt <- 10       # Thinning
  nb <- 450000   # Burn-in
  
  # Run MCMC 
  ipm_run <- nimbleMCMC(
    model = model,
    monitors = c(
      # Same monitoring parameters as original script
    #  "recruitment", 
      "aug31.ppb", "aug31.hwb",
      "male.h.ad.wmu", "male.s.ad.wmu", 
      "male.h.juv.wmu", "male.s.juv.wmu",
      "female.h.ad.wmu", "female.h.juv.wmu",  
      # "male.N.ad", "male.N.juv",
      # "female.N.ad", "female.N.juv",  
      "avg.ad.s.kf", "avg.juv.s.kf", "storage", 
      "juv.male.adj"
    ),
    niter = ni, 
    nburnin = nb, 
    thin = nt, 
    nchains = 1,  # One chain per core
    samplesAsCodaMCMC = TRUE,
    setSeed = 1235 + seed_offset,
    WAIC = FALSE
  )
  
  return(ipm_run)
}

# Parallel execution
start <- Sys.time()

# Run multiple chains in parallel
results <- foreach(
  i = 1:2,
  .packages = c("nimble")
) %dopar% {
  run_parallel_mcmc(i)
}

# Stop the cluster
stopCluster(cl)

end <- Sys.time()
print(paste("Total parallel processing time:", end - start))

#############################################################X
# Combine results -----
#############################################################X
combined_results <- coda::as.mcmc.list(lapply(results, coda::as.mcmc))
##------------------##X
# # if parallel need
# cl <- makeCluster(max(1, detectCores() - 1))
# registerDoParallel(cl)
# 
# # Parallel conversion of results to MCMC
# combined_results <- foreach(
#   result = results,
#   .packages = "coda"
# ) %dopar% {
#   as.mcmc(result)
# }

# Convert to mcmc.list
combined_results <- as.mcmc.list(combined_results)

# Stop the cluster
stopCluster(cl)
##------------------##X
# Save output
saveRDS(combined_results, "Data/Output/R_IPM_run23NoAbun23kf.rds")
save.image(file = "Data/Output/R_IPM_run23NoAbun23kf.Rdata")

#############################################################X
# Model diagnostics -----
#############################################################X
# You may need to adjust diagnostics for parallel results
library(MCMCvis)

# Summarize results
MCMCsummary(combined_results, params = c(
    "avg.ad.s.kf", "avg.juv.s.kf",
    "male.h.ad.wmu", "female.h.juv.wmu",
    "female.N.ad", "female.N.juv"
))


# Trace plots (adjust as needed for parallel processing)
MCMCvis::MCMCtrace(combined_results, pdf = FALSE, 
                   params = c("juv.male.adj"), iter = 1200000)

# Diagnostic checks
# Note: These may need modification for parallel results
lapply(combined_results, coda::geweke.diag)
lapply(combined_results, coda::effectiveSize)
lapply(combined_results, coda::raftery.diag)
lapply(combined_results, coda::heidel.diag)

# Done