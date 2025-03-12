#############################################################X
# Simulating data to test abundance model modifications
# August 3, 2024
# Based on the essence of Gonnerman et al.
# males
#############################################################X

# Clear the workspace and free up memory
rm(list=ls())
gc()

# Call library
library(nimble)

# Set seed for reproducibility
set.seed(1235)  

# read in recruitment data from previous IPM run (commented out for this example)
rec_df <- readRDS("Data/Output/Recruitment.rds")
#load("simTest_abundance.RData")
#############################################################X
# Define the data and constants to be used in the NIMBLE model
male.n.wmu <- 10  # Number of Wildlife Management Units (WMUs) for males
true.occasions <- 10  # Number of time occasions (e.g., years) to simulate

# Simulate initial abundance for each WMU
th.year1.male.ad <- rpois(male.n.wmu, lambda = 1000)  # Initial abundance of adult males for each WMU
th.year1.male.juv <- rpois(male.n.wmu, lambda = 1000)  # Initial abundance of adult males for each WMU

# Simulate harvest rates for each WMU over time
male.h.ad.wmu <- matrix(runif(male.n.wmu * true.occasions, 
                              min = 0.003, max = 0.303), 
                        nrow = true.occasions, ncol = male.n.wmu)

male.h.juv.wmu <- matrix(runif(male.n.wmu * true.occasions, 
                              min = 0.003, max = 0.303), 
                        nrow = true.occasions, ncol = male.n.wmu)

# Simulate survival rates for each WMU over time
male.s.ad.wmu <- matrix(runif(male.n.wmu * true.occasions, 
                              min = 0.386, max = 0.834), 
                        nrow = true.occasions, ncol = male.n.wmu)

male.s.juv.wmu <- matrix(runif(male.n.wmu * true.occasions, 
                              min = 0.386, max = 0.834), 
                        nrow = true.occasions, ncol = male.n.wmu)
# male.s.ad.wmu <- 0.8

# Extract recruitment parameters from the recruitment data frame
# Assuming 'rec_df' is a predefined data frame with the necessary data
ppb_values <- subset(rec_df, type == "ppb")$Median  # Median values for PPB type
hwb_values <- subset(rec_df, type == "HWB")$Median  # Median values for HWB type

# Simulate recruitment parameters using the provided values
aug31.ppb <- matrix(sample(ppb_values, male.n.wmu * true.occasions, 
                           replace = TRUE), 
                    nrow = true.occasions, ncol = male.n.wmu)
aug31.hwb <- matrix(sample(hwb_values, male.n.wmu * true.occasions, 
                           replace = TRUE), 
                    nrow = true.occasions, ncol = male.n.wmu)

# Initialize harvest observations and other matrices
harvest.ad.spring <- matrix(0, nrow = true.occasions, ncol = male.n.wmu)
harvest.juv.spring <- matrix(0, nrow = true.occasions, ncol = male.n.wmu)
N.lambda.ad.male <- matrix(0, nrow = 1, ncol = male.n.wmu)
N.lambda.juv.male <- matrix(0, nrow = 1, ncol = male.n.wmu)
male.N.ad <- matrix(0, nrow = true.occasions, ncol = male.n.wmu)
male.N.juv <- matrix(0, nrow = true.occasions, ncol = male.n.wmu)
male.surv.ad.tot <- matrix(0, nrow = true.occasions, ncol = male.n.wmu)
male.surv.juv.tot <- matrix(0, nrow = true.occasions, ncol = male.n.wmu)
male.recruitment <- matrix(0, nrow = true.occasions, ncol = male.n.wmu)
delta.ad.male <- matrix(0, nrow = true.occasions, ncol = male.n.wmu)
delta.juv.male <- matrix(0, nrow = true.occasions, ncol = male.n.wmu)
#############################################################X
# Initialization Loop
for (u in 1:male.n.wmu) {
  # Initial abundance for each WMU for adult males
  N.lambda.ad.male[u] <- (th.year1.male.ad[u] + 1) / male.h.ad.wmu[1, u]
  male.N.ad[1, u] <- rpois(1, N.lambda.ad.male[u])
  
  # Initial abundance for each WMU for juvenile males
  N.lambda.juv.male[u] <- (th.year1.male.juv[u] + 1) / male.h.juv.wmu[1, u]
  male.N.juv[1, u] <- rpois(1, N.lambda.juv.male[u])
  
  # Loop over subsequent time occasions (excluding the last one)
  for (t in 1:(true.occasions - 1)) {
    # # Number of surviving adult males from time t to t+1
    # male.surv.ad.tot[t, u] <- rbinom(1, prob = male.s.ad.wmu[t, u], 
    #                                  size = round(male.N.ad[t, u]))
    # 
    # # Recruitment for the next time occasion
    # male.recruitment[t, u] <- ((male.N.ad[t, u] / 4 * aug31.hwb[t, u]) * aug31.ppb[t, u]) / 2
    # 
    # # Process model for the population count of adult males
    # delta.male[t, u] <- male.surv.ad.tot[t, u] + male.recruitment[t, u]
    # 
    # # Update the number of adult males for the next time occasion
    # male.N.ad[t + 1, u] <- rpois(1, delta.male[t, u])
    # 
    # Number of surviving adult males from time t-1 to t
    male.surv.ad.tot[t, u] <- rbinom(1, prob = male.s.ad.wmu[t, u], 
                                     size = round(male.N.ad[t, u]))
    
    # Number of surviving adult males from time t-1 to t
    male.surv.juv.tot[t, u] <- rbinom(1, prob = male.s.juv.wmu[t, u], 
                                     size = round(male.N.juv[t, u]))
    
    # Recruitment for the next time occasion
    male.recruitment[t, u] <- ((male.N.ad[t, u] / 4 * aug31.hwb[t, u]) * aug31.ppb[t, u]) / 2
    
    # Process model for the population count of adult males
    # Number of adult birds to survive and juveniles that transition into A from t to t+1
    delta.ad.male[t, u] <- male.surv.ad.tot[t, u] + male.surv.juv.tot[t, u]
    
    # Number of juvenile birds recruited into the population
    delta.juv.male[t, u] <- male.recruitment[t, u] 
    
    # Update the number of adult males for the next time occasion
    male.N.ad[t+1, u] <- rpois(1, delta.ad.male[t, u])
    male.N.juv[t+1, u] <- rpois(1, delta.juv.male[t, u])
  } # end t
} # end u

# Second loop to calculate harvest.ad.spring after initializing male.N.ad
for (u in 1:male.n.wmu) {
  for (t in 1:true.occasions) {
    # Total harvest observation for adult males in spring
    harvest.ad.spring[t, u] <- rbinom(n = 1, prob = male.h.ad.wmu[t, u], 
                                      size = round(male.N.ad[t, u]))
    
    # Total harvest observation for adult males in spring
    harvest.juv.spring[t, u] <- rbinom(n = 1, prob = male.h.juv.wmu[t, u], 
                                      size = round(male.N.juv[t, u]))
  } # end t
} # end u
#############################################################X
# Define the NIMBLE model code
#############################################################X
nimble_code <- nimbleCode({
  for (u in 1:male.n.wmu) {
    # Initial abundance for each WMU for adult males
    N.lambda.ad.male[u] <- (th.year1.male.ad[u] + 1) / male.h.ad.wmu[1, u]
    male.N.ad[1, u] ~ dpois(N.lambda.ad.male[u])
    
    # Initial abundance for each WMU for juvenile males
    N.lambda.juv.male[u] <- (th.year1.male.juv[u] + 1) / male.h.juv.wmu[1, u]
    male.N.juv[1, u] ~ dpois(N.lambda.juv.male[u])
    
    # Loop over time occasions (starting from the second occasion)
    for (t in 2:true.occasions) {
      # Number of surviving adult males from time t-1 to t
      male.surv.ad.tot[t-1, u] ~ dbin(prob = male.s.ad.wmu[t-1, u], 
                                      size = round(male.N.ad[t-1, u]))
      
      # Number of surviving adult males from time t-1 to t
      male.surv.juv.tot[t-1, u] ~ dbin(prob = male.s.juv.wmu[t-1, u], 
                                      size = round(male.N.juv[t-1, u]))
      
      # Recruitment for the next time occasion
      male.recruitment[t-1, u] <- ((male.N.ad[t-1, u] / 4 * aug31.hwb[t-1, u]) * aug31.ppb[t-1, u]) / 2
      
      # Process model for the population count of adult males
      # Number of adult birds to survive and juveniles that transition into A from t to t+1
      delta.ad.male[t-1, u] <- male.surv.ad.tot[t-1, u] + male.surv.juv.tot[t-1, u]
      
      # Number of juvenile birds recruited into the population
      delta.juv.male[t-1, u] <- male.recruitment[t-1, u] 
      
      # Update the number of adult males for the next time occasion
      male.N.ad[t, u] ~ dpois(delta.ad.male[t-1, u])
      male.N.juv[t, u] ~ dpois(delta.juv.male[t-1, u])
      
      # Total harvest observation for adult males in spring
      harvest.ad.spring[t, u] ~ dbin(prob = male.h.ad.wmu[t, u], 
                                     size = round(male.N.ad[t, u]))
      harvest.juv.spring[t, u] ~ dbin(prob = male.h.juv.wmu[t, u], 
                                     size = round(male.N.juv[t, u]))
    } # end t
  } # end u
})

#############################################################X
# Define the data and constants for the model
data <- list(
  aug31.ppb = aug31.ppb,
  aug31.hwb = aug31.hwb,
  
  harvest.ad.spring = harvest.ad.spring,
  male.h.ad.wmu = male.h.ad.wmu,
  male.s.ad.wmu = male.s.ad.wmu,
  th.year1.male.ad = th.year1.male.ad,
  
  harvest.juv.spring = harvest.juv.spring,
  male.h.juv.wmu = male.h.juv.wmu,
  male.s.juv.wmu = male.s.juv.wmu,
  th.year1.male.juv = th.year1.male.juv
)

constants <- list(
  male.n.wmu = male.n.wmu,
  true.occasions = true.occasions
)


# Define the initial values for the model parameters

# A note from Abbey Feuka: "Because the model removes individuals from the estimated N 
# at each time step, I needed to specify initial values for N that were 
# sufficiently high enough to not send it into the negatives when 
# animals are removed
# https://groups.google.com/g/nimble-users/c/fKx8OD7Qn9s

inits <- list(
  male.N.ad = matrix(50000, nrow = true.occasions, ncol = male.n.wmu),
  male.surv.ad.tot = matrix(100, nrow = true.occasions-1, ncol = male.n.wmu),
  male.N.juv = matrix(50000, nrow = true.occasions, ncol = male.n.wmu),
  male.surv.juv.tot = matrix(100, nrow = true.occasions-1, ncol = male.n.wmu)
)


# Build the NIMBLE model
nimble_model <- nimbleModel(
  code = nimble_code,
  data = data,
  constants = constants,
  inits = inits
)

# Run the MCMC
samples <- nimbleMCMC(
  model = nimble_model,
  monitors = c(
    "male.N.ad", "male.recruitment", "male.surv.ad.tot", "N.lambda.ad.male", "delta.ad.male",
    "male.N.juv", "male.surv.juv.tot", "N.lambda.juv.male", "delta.juv.male"
  ),
  niter = 250000, 
  nburnin = 40000, 
  nchains = 2, 
  thin = 5, 
  setSeed = TRUE,
  samplesAsCodaMCMC = TRUE
  )

# Ensure the samples object is correctly formatted for MCMCvis
class(samples) # should be mcmc.list

# Visualize the trace plots
MCMCvis::MCMCtrace(samples, pdf = FALSE, params = "delta.juv.male")

# Get a summary of the samples
MCMCvis::MCMCsummary(samples)

# Extract specific posterior values using median
posterior_median_surv <- MCMCvis::MCMCpstr(samples, params = "male.surv.ad.tot", func = median)
posterior_median_surv <- MCMCvis::MCMCpstr(samples, params = "male.surv.juv.tot", func = median)
posterior_median_delta.ad.male <- MCMCvis::MCMCpstr(samples, params = "delta.ad.male", func = median)
posterior_median_delta.juv.male <- MCMCvis::MCMCpstr(samples, params = "delta.juv.male", func = median)
posterior_median_N.lambda.ad.male <- MCMCvis::MCMCpstr(samples, params = "N.lambda.ad.male", func = median)
posterior_median_N.lambda.juv.male <- MCMCvis::MCMCpstr(samples, params = "N.lambda.juv.male", func = median)
posterior_median_male.recruitment <- MCMCvis::MCMCpstr(samples, params = "male.recruitment", func = median)

# Recovered sim values well!


# # Males
# for (u in 1:male.n.wmu) {
#   # Initial abundance for each WMU for males
#   N.lambda.male[u] <- (th.year1.male.ad[u] + 1) / male.h.ad.wmu[1, u]
#   male.N.ad[1, u] ~ dpois(N.lambda.male[u])
#   
#   # Loop over time occasions (starting from the second occasion)
#   for (t in 2:true.occasions) {
#     # Number of surviving adult males from time t-1 to t
#     male.surv.ad.tot[t-1, u] ~ dbin(prob = male.s.ad.wmu[t-1, u], 
#                                     size = round(male.N.ad[t-1, u]))
#     
#     # Recruitment for the next time occasion
#    # recruitment[t-1, u] <- ((male.N.ad[t-1, u] / 4 * aug31.hwb[t-1, u]) * aug31.ppb[t-1, u]) / 2
#     
#     # Process model for the population count of adult males
#     delta.male[t-1, u] <- male.surv.ad.tot[t-1, u] + recruitment[t-1, u]
#     
#     # Update the number of adult males for the next time occasion
#     male.N.ad[t, u] ~ dpois(delta.male[t-1, u])
#     
#     # Total harvest observation for adult males in spring
#     harvest.ad.spring[t, u] ~ dbin(prob = male.h.ad.wmu[t, u], 
#                                    size = round(male.N.ad[t, u]))
#   } # end t
# } # end u
# 