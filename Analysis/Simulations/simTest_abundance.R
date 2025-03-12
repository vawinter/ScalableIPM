#############################################################X
# Simulating data to test abundance model modifications
# August 3, 2024
# Based on the essence of Gonnerman et al.
# females and males
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
female.n.wmu <- 10  # Number of Wildlife Management Units (WMUs) for females
true.occasions <- 10  # Number of time occasions (e.g., years) to simulate

# Simulate initial abundance for each WMU
th.year1.female.ad <- rpois(female.n.wmu, lambda = 1000)  # Initial abundance of adult females for each WMU

# Simulate harvest rates for each WMU over time
female.h.ad.wmu <- matrix(runif(female.n.wmu * true.occasions, 
                                min = 0.014, max = 0.022), 
                          nrow = true.occasions, ncol = female.n.wmu)

# Simulate survival rates for each WMU over time
female.s.ad.wmu <- matrix(runif(female.n.wmu * true.occasions, 
                                min = 0.453, max = 0.663), 
                          nrow = true.occasions, ncol = female.n.wmu)

# Extract recruitment parameters from the recruitment data frame
# Assuming 'rec_df' is a predefined data frame with the necessary data
ppb_values <- subset(rec_df, type == "ppb")$Median  # Median values for PPB type
hwb_values <- subset(rec_df, type == "HWB")$Median  # Median values for HWB type

# Simulate recruitment parameters using the provided values
aug31.ppb <- matrix(sample(ppb_values, female.n.wmu * true.occasions, 
                           replace = TRUE), 
                    nrow = true.occasions, ncol = female.n.wmu)
aug31.hwb <- matrix(sample(hwb_values, female.n.wmu * true.occasions, 
                           replace = TRUE), 
                    nrow = true.occasions, ncol = female.n.wmu)

# Initialize harvest observations and other matrices
harvest.ad.fall <- matrix(0, nrow = true.occasions, ncol = female.n.wmu)
N.lambda.female <- matrix(0, nrow = 1, ncol = female.n.wmu)
female.N.ad <- matrix(0, nrow = true.occasions, ncol = female.n.wmu)
female.surv.ad.tot <- matrix(0, nrow = true.occasions, ncol = female.n.wmu)
female.recruitment <- matrix(0, nrow = true.occasions, ncol = female.n.wmu)
delta.female <- matrix(0, nrow = true.occasions, ncol = female.n.wmu)

#############################################################X
# Define the data and constants to be used in the NIMBLE model
male.n.wmu <- 10  # Number of Wildlife Management Units (WMUs) for males
true.occasions <- 10  # Number of time occasions (e.g., years) to simulate

# Simulate initial abundance for each WMU
th.year1.male.ad <- rpois(male.n.wmu, lambda = 1000)  # Initial abundance of adult males for each WMU

# Simulate harvest rates for each WMU over time
male.h.ad.wmu <- matrix(runif(male.n.wmu * true.occasions, 
                              min = 0.003, max = 0.303), 
                        nrow = true.occasions, ncol = male.n.wmu)

# Simulate survival rates for each WMU over time
male.s.ad.wmu <- matrix(runif(male.n.wmu * true.occasions, 
                              min = 0.386, max = 0.834), 
                        nrow = true.occasions, ncol = male.n.wmu)
# male.s.ad.wmu <- 0.8

# # Extract recruitment parameters from the recruitment data frame
# # Assuming 'rec_df' is a predefined data frame with the necessary data
# ppb_values <- subset(rec_df, type == "ppb")$Median  # Median values for PPB type
# hwb_values <- subset(rec_df, type == "HWB")$Median  # Median values for HWB type
# 
# # Simulate recruitment parameters using the provided values
# aug31.ppb <- matrix(sample(ppb_values, male.n.wmu * true.occasions, 
#                            replace = TRUE), 
#                     nrow = true.occasions, ncol = male.n.wmu)
# aug31.hwb <- matrix(sample(hwb_values, male.n.wmu * true.occasions, 
#                            replace = TRUE), 
#                     nrow = true.occasions, ncol = male.n.wmu)

# Initialize harvest observations and other matrices
harvest.ad.spring <- matrix(0, nrow = true.occasions, ncol = male.n.wmu)
N.lambda.male <- matrix(0, nrow = 1, ncol = male.n.wmu)
male.N.ad <- matrix(0, nrow = true.occasions, ncol = male.n.wmu)
male.surv.ad.tot <- matrix(0, nrow = true.occasions, ncol = male.n.wmu)
male.recruitment <- matrix(0, nrow = true.occasions, ncol = male.n.wmu)
delta.male <- matrix(0, nrow = true.occasions, ncol = male.n.wmu)
#############################################################X
# Initialization Loop
for (u in 1:female.n.wmu) {
  # Initial abundance for each WMU for females
  N.lambda.female[u] <- (1 + th.year1.female.ad[u]) / female.h.ad.wmu[1, u]
  female.N.ad[1, u] <- rpois(1, N.lambda.female[u])
  
  # Loop over subsequent time occasions (excluding the last one)
  for (t in 1:(true.occasions - 1)) {
    # Number of surviving adult females from time t to t+1
    female.surv.ad.tot[t, u] <- rbinom(1, prob = female.s.ad.wmu[t, u], 
                                       size = round(female.N.ad[t, u]))
    
    # Recruitment for the next time occasion
    female.recruitment[t, u] <- ((female.N.ad[t, u] / 4 * aug31.hwb[t, u]) * aug31.ppb[t, u]) / 2
    
    # Process model for the population count of adult females
    delta.female[t, u] <- female.surv.ad.tot[t, u] + female.recruitment[t, u]
    
    # Update the number of adult females for the next time occasion
    female.N.ad[t + 1, u] <- rpois(1, delta.female[t, u])
  } # end t
} # end u

# Second loop to calculate harvest.ad.spring after initializing female.N.ad
for (u in 1:female.n.wmu) {
  for (t in 1:true.occasions) {
    # Total harvest observation for adult females in spring
    harvest.ad.fall[t, u] <- rbinom(n = 1, prob = female.h.ad.wmu[t, u], 
                                    size = round(female.N.ad[t, u]))
  } # end t
} # end u

#############################################################X
# Initialization Loop
for (u in 1:male.n.wmu) {
  # Initial abundance for each WMU for males
  N.lambda.male[u] <- (1 + th.year1.male.ad[u]) / male.h.ad.wmu[1, u]
  male.N.ad[1, u] <- rpois(1, N.lambda.male[u])
  
  # Loop over subsequent time occasions (excluding the last one)
  for (t in 1:(true.occasions - 1)) {
    # Number of surviving adult males from time t to t+1
    male.surv.ad.tot[t, u] <- rbinom(1, prob = male.s.ad.wmu[t, u], 
                                     size = round(male.N.ad[t, u]))
    
    # Recruitment for the next time occasion
   # male.recruitment[t, u] <- ((male.N.ad[t, u] / 4 * aug31.hwb[t, u]) * aug31.ppb[t, u]) / 2
    
    # Process model for the population count of adult males
    delta.male[t, u] <- male.surv.ad.tot[t, u] + female.recruitment[t, u]
    
    # Update the number of adult males for the next time occasion
    male.N.ad[t + 1, u] <- rpois(1, delta.male[t, u])
  } # end t
} # end u

# Second loop to calculate harvest.ad.spring after initializing male.N.ad
for (u in 1:male.n.wmu) {
  for (t in 1:true.occasions) {
    # Total harvest observation for adult males in spring
    harvest.ad.spring[t, u] <- rbinom(n = 1, prob = male.h.ad.wmu[t, u], 
                                      size = round(male.N.ad[t, u]))
  } # end t
} # end u
#############################################################X
# Define the NIMBLE model code
#############################################################X
nimble_code <- nimbleCode({
  # Males
  for (u in 1:male.n.wmu) {
    # Initial abundance for each WMU for males
    N.lambda.male[u] <- (th.year1.male.ad[u] + 1) / male.h.ad.wmu[1, u]
    male.N.ad[1, u] ~ dpois(N.lambda.male[u])
    
    # Loop over time occasions (starting from the second occasion)
    for (t in 2:true.occasions) {
      # Number of surviving adult males from time t-1 to t
      male.surv.ad.tot[t-1, u] ~ dbin(prob = male.s.ad.wmu[t-1, u], 
                                      size = round(male.N.ad[t-1, u]))
      
      # Recruitment for the next time occasion
      # recruitment[t-1, u] <- ((male.N.ad[t-1, u] / 4 * aug31.hwb[t-1, u]) * aug31.ppb[t-1, u]) / 2
      
      # Process model for the population count of adult males
      delta.male[t-1, u] <- male.surv.ad.tot[t-1, u] + recruitment[t-1, u]
      
      # Update the number of adult males for the next time occasion
      male.N.ad[t, u] ~ dpois(delta.male[t-1, u])
      
      # Total harvest observation for adult males in spring
      harvest.ad.spring[t, u] ~ dbin(prob = male.h.ad.wmu[t, u], 
                                     size = round(male.N.ad[t, u]))
    } # end t
  } # end u
  
  #---------------------------------------------------------------------#X
  # Females
  for (u in 1:female.n.wmu) {
    # Initial abundance for each WMU for females
    N.lambda.female[u] <- (th.year1.female.ad[u] + 1) / female.h.ad.wmu[1, u]
    female.N.ad[1, u] ~ dpois(N.lambda.female[u])
    
    # Loop over time occasions (starting from the second occasion)
    for (t in 2:true.occasions) {
      # Number of surviving adult females from time t-1 to t
      female.surv.ad.tot[t-1, u] ~ dbin(prob = female.s.ad.wmu[t-1, u],
                                        size = round(female.N.ad[t-1, u]))
      
      # Recruitment for the next time occasion (assumes even sex ratio)
      recruitment[t-1, u] <- ((female.N.ad[t-1, u] * aug31.hwb[t-1, u]) * aug31.ppb[t-1, u]) / 2
      
      # Process model for the population count of adult females
      delta.female[t-1, u] <- female.surv.ad.tot[t-1, u] + recruitment[t-1, u]
      
      # Update the number of adult females for the next time occasion
      female.N.ad[t, u] ~ dpois(delta.female[t-1, u])
      
      # Total harvest observation for adult females in spring
      harvest.ad.fall[t, u] ~ dbin(prob = female.h.ad.wmu[t, u],
                                   size = round(female.N.ad[t, u]))
    } # end t
  } # end u
})

#############################################################X
# Define the data and constants for the model
data <- list(
  harvest.ad.fall = harvest.ad.fall,
  female.h.ad.wmu = female.h.ad.wmu,
  female.s.ad.wmu = female.s.ad.wmu,
  aug31.ppb = aug31.ppb,
  aug31.hwb = aug31.hwb,
  th.year1.female.ad = th.year1.female.ad,
  harvest.ad.spring = harvest.ad.spring,
  male.h.ad.wmu = male.h.ad.wmu,
  male.s.ad.wmu = male.s.ad.wmu,
  th.year1.male.ad = th.year1.male.ad
)

constants <- list(
  female.n.wmu = female.n.wmu,
  true.occasions = true.occasions,
  male.n.wmu = male.n.wmu
)


# Define the initial values for the model parameters

# A note from Abbey Feuka: "Because the model removes individuals from the estimated N 
# at each time step, I needed to specify initial values for N that were 
# sufficiently high enough to not send it into the negatives when 
# animals are removed
# https://groups.google.com/g/nimble-users/c/fKx8OD7Qn9s

inits <- list(
  female.N.ad = matrix(50000, nrow = true.occasions, ncol = female.n.wmu),
  female.surv.ad.tot = matrix(100, nrow = true.occasions-1, ncol = female.n.wmu),
  male.N.ad = matrix(50000, nrow = true.occasions, ncol = male.n.wmu),
  male.surv.ad.tot = matrix(100, nrow = true.occasions-1, ncol = male.n.wmu)
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
    "female.N.ad", "female.surv.ad.tot", 
    "N.lambda.female", "delta.female",     
    "male.N.ad", "recruitment", "male.surv.ad.tot", "N.lambda.male", 
    "delta.male"
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
MCMCvis::MCMCtrace(samples, pdf = FALSE, params = c("female.N.ad", "male.N.ad"))

# Get a summary of the samples
MCMCvis::MCMCsummary(samples)

# Extract specific posterior values using median
posterior_median_surv <- MCMCvis::MCMCpstr(samples, params = "female.surv.ad.tot", func = median)
posterior_median_delta <- MCMCvis::MCMCpstr(samples, params = "delta.female", func = median)

posterior_median_N_lambda <- MCMCvis::MCMCpstr(samples, params = "N.lambda.female", func = median)
posterior_median_female.N.ad <- MCMCvis::MCMCpstr(samples, params = "female.N.ad", func = median)
