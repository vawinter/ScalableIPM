###############################################################################X
# Parallelized V Integrated Population Model (IPM) Script
# Modified for parallel processing
# https://r-nimble.org/nimbleExamples/parallelizing_NIMBLE.html
###############################################################################X
# Note: this script is used for fitting either the O and V IPM. The only change is 
# the model script that is loaded in. Otherwise, data inputs are the same,
# the only change are the model priors.

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

### Nimble model set up ----
load("Data/Operational_IPM_setup-data/O_IPM_Nimble_data_setup24.RData")
source("Models/vague_ipm.R")

##################################################################X
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#              ######### THINGS TO NOTE: #########
#
# * Variable names correspond to the associated model. For example:
# **'hwb' = proportion of hens with a brood
# **'ppb' = poults-per-brood e.g., poult-to-hen ratio
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


#############################################################X
# Parallelization setup ----
#############################################################X
# Detect number of cores, leave one free for system processes
num_cores <- max(1, detectCores() - 1)

# Register parallel backend
cl <- makeCluster(num_cores)
registerDoParallel(cl)

# Export necessary objects to cluster
clusterExport(cl, c("nimble.data", "consts", "inits", "v_ipm"))

# Parallel function to run MCMC
run_parallel_mcmc <- function(seed_offset) {
  library(nimble)
  
  # Option to set unique seed for each chain
  set.seed(1235 + seed_offset)
  
  # Create the Nimble model
  model <- nimbleModel(
    code = v_ipm,
    data = nimble.data,
    constants = consts,
    inits = inits
  )
  
  # MCMC specifications
  ni <- 1400000  # Total iterations
  nt <- 10       # Thinning
  nb <- 550000   # Burn-in
  
  # Run MCMC 
  ipm_run <- nimbleMCMC(
    model = model,
    monitors = c(
      "recruitment", "aug31.ppb", "aug31.hwb",
      "male.h.ad.wmu", "male.s.ad.wmu", 
      "male.h.juv.wmu", "male.s.juv.wmu",
      "female.h.ad.wmu", "female.h.juv.wmu",  
      "male.N.ad", "male.N.juv",
      "female.N.ad", "female.N.juv",  
      "avg.ad.s.kf", "avg.juv.s.kf"
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
# 
# # Convert to mcmc.list
# combined_results <- as.mcmc.list(combined_results)

# Stop the cluster
stopCluster(cl)
##------------------##X
# Save output
saveRDS(combined_results, paste0("Data/Output/", format(Sys.Date(), "%Y%m%d"), "_V_IPM_run.rds"))
save.image(file = "Data/Output/V_IPM_run.Rdata")

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
MCMCvis::MCMCtrace(combined_results[[1]], pdf = FALSE, 
                   params = c("male.N.ad"), iter = 1200000)

# Diagnostic checks
# Note: These may need modification for parallel results
lapply(combined_results, coda::geweke.diag)
lapply(combined_results, coda::effectiveSize)
lapply(combined_results, coda::raftery.diag)
lapply(combined_results, coda::heidel.diag)

# Done