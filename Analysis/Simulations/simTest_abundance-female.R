#############################################################X
# Simulating data to test abundance model modifications
# August 6, 2024
# Based on the essence of Gonnerman et al.
# females
#############################################################X

# Clear the workspace and free up memory
rm(list=ls())
gc()

# Call library
library(nimble)

# Set seed for reproducibility
set.seed(1235)  

# Example of reading in recruitment data from previous IPM run
rec_df <- readRDS("Data/Output/Recruitment.rds")

# note:
plot(density(rbeta(true.occasions * female.n.wmu, 3, 60)))
# [1] 0.003992896 0.102618377
range(rbeta(true.occasions * female.n.wmu, 3, 95))
# [1] 0.00731109 0.13370243
range(rbeta(true.occasions * female.n.wmu, 2, 97))
# [1] 0.004138058 0.059136682
range(rbeta(true.occasions * female.n.wmu, 5, 97))
# [1] 0.007256256 0.113716362
range(rbeta(true.occasions * female.n.wmu, 6, 97))
# [1] 0.01726918 0.13549650
range(rbeta(true.occasions * female.n.wmu, 7, 97))
# [1] 0.02330135 0.12955827
range(rbeta(true.occasions * female.n.wmu, 8, 97))
# [1] 0.02653213 0.14649130
#############################################################X
# Step 1: Define Constants and Initialize Matrices
# Define the number of WMUs and time occasions
female.n.wmu <- 10  # Number of WMUs
true.occasions <- 100  # Number of time occasions (e.g., years)

# Define initial abundance for adult females and poults
th.year1.female.ad <- rpois(female.n.wmu, lambda = 100)  # Adult females
th.year1.poult <- rpois(female.n.wmu, lambda = 100)  # Poults

female.seber.recov <- runif(1,0,1)

# Harvest rate and recruitment parameters
female.h.ad.wmu <- matrix(0.1, 
                          nrow = true.occasions, ncol = female.n.wmu)
female.s.ad.wmu <- 1 - (female.h.ad.wmu/female.seber.recov)

aug31.ppb <- matrix(runif(true.occasions * female.n.wmu, min = 0.8, max = 1.2),
                    nrow = true.occasions, ncol = female.n.wmu)
aug31.hwb <- matrix(runif(true.occasions * female.n.wmu, min = 1.5, max = 2.5),
                    nrow = true.occasions, ncol = female.n.wmu)

avg.juv.s.kf <- rep(0.7, female.n.wmu)  # Juvenile survival

# Initialize matrices for abundance, survival, and recruitment
female.N.ad <- matrix(0, nrow = true.occasions, ncol = female.n.wmu)
N.poult <- matrix(0, nrow = true.occasions, ncol = female.n.wmu)
female.surv.ad.tot <- matrix(0, nrow = true.occasions, ncol = female.n.wmu)
surv.poult.tot <- matrix(0, nrow = true.occasions, ncol = female.n.wmu)
recruitment <- matrix(0, nrow = true.occasions, ncol = female.n.wmu)

# Harvest data
harvest.ad.fall <- matrix(0, nrow = true.occasions, ncol = female.n.wmu)
harvest.poult.fall <- matrix(0, nrow = true.occasions, ncol = female.n.wmu)
#############################################################X
# Step 2: Simulate Data for Each Time Occasion
# Initial Time Step (t = 1): Simulate the initial abundance of adult females and poults.
for (u in 1:female.n.wmu) {
  # Initial abundance for adult females and poults
  N.lambda.ad.female <- th.year1.female.ad[u] / female.h.ad.wmu[1, u]
  female.N.ad[1, u] <- rpois(1, N.lambda.ad.female)
  
  N.lambda.poult <- th.year1.poult[u] / female.h.ad.wmu[1, u]
  N.poult[1, u] <- rpois(1, N.lambda.poult)
}


# Subsequent Time Steps (t = 2 to true.occasions): 
# Simulate the survival, recruitment, and harvest data for each time step.
for (t in 1:(true.occasions - 1)) {
  for (u in 1:female.n.wmu) {
    # Survival of adult females and poults
    female.surv.ad.tot[t, u] <- rbinom(1, size = female.N.ad[t, u], 
                                       prob = female.s.ad.wmu[t, u])
    surv.poult.tot[t, u] <- rbinom(1, size = N.poult[t, u], 
                                   prob = avg.juv.s.kf[u])
    
    # Recruitment
    recruitment[t, u] <- (female.N.ad[t, u] * aug31.hwb[t, u] * aug31.ppb[t, u])
    
    # Update population for next time period
    delta.ad.female <- female.surv.ad.tot[t, u] + (surv.poult.tot[t, u] / 2)
    delta.poult <- recruitment[t, u]
    
    female.N.ad[t + 1, u] <- rpois(1, delta.ad.female)
    N.poult[t + 1, u] <- rpois(1, delta.poult)
    
    # Harvest data
    harvest.ad.fall[t, u] <- rbinom(1, size = female.N.ad[t, u], 
                                    prob = female.h.ad.wmu[t, u])
    harvest.poult.fall[t, u] <- rbinom(1, size = N.poult[t, u], 
                                       prob = female.h.ad.wmu[t, u])
  }
}


# Now you have generated data for a more stable population
#############################################################X
# Define the NIMBLE model code
#############################################################X
nimble_code <- nimbleCode({
  #----------------------------------------------------------#
  # DRM: WMU survival rates
  #----------------------------------------------------------#
  # Proportion mortality due to hunting (Seber parameterization: (1-s)r)
  female.seber.recov ~ dunif(0,1) # Prior for proportion mortality due to hunting - adults
  
  for (t in 1:(female.true.occasions)) {
    for (u in 1:female.n.wmu) {

      # Option 2: Treating harvest rate as a latent parameter with a Beta prior.
      female.h.ad.wmu[t,u] ~ dbeta(shape1 = 3, shape2 = 30) # this was 3, 97 - but I expanded 
      # so the range was 0.02653213 0.14649130
      
      
      # Option 2: Calculate survival based on harvest rate and reporting rate
      female.s.ad.wmu[t, u] <- 1 - (female.h.ad.wmu[t, u]/female.seber.recov)
      
    }#
  }#t
  
  for (u in 1:female.n.wmu) {
    avg.juv.s.kf[u] ~ dunif(0.6, 0.8)
   
     # Initial abundance for adult females in each WMU
    N.lambda.ad.female[u] <- (th.year1.female.ad[u]) / female.h.ad.wmu[1, u]
    female.N.ad[1, u] ~ dpois(N.lambda.ad.female[u])
    
    # Initial abundance for juvenile females in each WMU
    N.lambda.poult[u] <- (th.year1.poult[u]) / female.h.ad.wmu[1, u]
    N.poult[1, u] ~ dpois(N.lambda.poult[u])
    
  } # end u
  
  # Loop over time occasions (starting from the second occasion)
  for (t in 1:(female.true.occasions-1)) {
    for (u in 1:female.n.wmu) {
      
      # Number of surviving adult females from time t-1 to t
      female.surv.ad.tot[t, u] ~ dbin(size = female.N.ad[t, u],
                                      prob = female.s.ad.wmu[t, u])
      
      
      # Number of surviving poults from time t-1 to t
      surv.poult.tot[t, u] ~ dbin(size = N.poult[t, u],
                                  prob = avg.juv.s.kf[u])
      
      # Recruitment for the next time occasion
      recruitment[t, u] <- ((female.N.ad[t, u] * aug31.hwb[t, u]) * aug31.ppb[t, u])
      
      # Process model for adult female population count
      delta.ad.female[t, u] <- female.surv.ad.tot[t, u] + (surv.poult.tot[t, u]/2)
      
      # Process model for poult count
      delta.poult[t, u] <- recruitment[t, u]
      
      # Update the number of adult females for the next time occasion
      female.N.ad[t+1, u] ~ dpois(delta.ad.female[t, u])
      
      # Update the number of poults for the next time occasion
      N.poult[t+1, u] ~ dpois(delta.poult[t, u])
      
      # Total harvest observation for adult females in fall
      harvest.ad.fall[t, u] ~ dbin(prob = female.h.ad.wmu[t, u],
                                   size = female.N.ad[t, u])
      
      # Total harvest observation for poults in fall
      harvest.poult.fall[t, u] ~ dbin(prob = female.h.ad.wmu[t, u],
                                      size = N.poult[t, u])
      
    } # end t
  } # end u
})



# Define the data and constants for the model
data_female <- list(
  harvest.ad.fall = harvest.ad.fall,
  harvest.poult.fall = harvest.poult.fall,
  aug31.ppb = aug31.ppb,
  aug31.hwb = aug31.hwb,
  th.year1.poult = harvest.poult.fall[1,],
  th.year1.female.ad = harvest.ad.fall[1,]
)

constants_female <- list(
  female.n.wmu = female.n.wmu,
  female.true.occasions = true.occasions
)

# Define the initial values for the model parameters
inits_female <- list(
  female.N.ad = array(50000, dim = c(true.occasions, female.n.wmu)),
  N.poult = array(50000, dim = c(true.occasions, female.n.wmu))
)

# Build the NIMBLE model
nimble_model_female <- nimbleModel(
  code = nimble_code,
  data = data_female,
  constants = constants_female,
  inits = inits_female
)

# Run the MCMC for females
samples_female <- nimbleMCMC(
  model = nimble_model_female,
  monitors = c("avg.juv.s.kf",
               # "N.poult",
               # "female.surv.ad.tot", 
               # "recruitment", 
               # "delta.ad.female",
               "female.h.ad.wmu",
               "female.s.ad.wmu"
               ),
  niter = 100000,
  nburnin = 60000,
  nchains = 1,
  thin = 1,
  setSeed = TRUE,
  samplesAsCodaMCMC = TRUE
)


#############################################################X
##         Look at parameter recovery
#############################################################X
# Visualize the trace plots
MCMCvis::MCMCtrace(samples_female, pdf = F, params = c("female.h.ad.wmu"))
#############################################################X
par(mfrow = c(3, 2))  # Set up the plotting window with 3 rows and 2 columns

for (t in 1:(true.occasions)) {
  for (i in 1:female.n.wmu) {
    # Construct the variable name dynamically for 'female.h.ad.wmu[t, i]'
    samp <- paste0("female.h.ad.wmu[", i, ", ", t, "]")
    
      # Plot trace of the sampled values
      matplot(samples_female[, samp], type = "l", 
              main = paste0("female.h.ad.wmu[", i, ", ", t, "]"), 
              ylab = "Sample values")
      
      # Add a reference line (e.g., rr[t, i] for the current time and WMU)
      abline(h = female.h.ad.wmu[i, t], col = "red", lwd = 2)
      
      # Plot density of the sampled values
      plot(density(samples_female[, samp]), main = "Density",
           xlab = "Beta Coefficient", ylab = "Density")
      
      # Add a vertical line at the true value
      abline(v = female.h.ad.wmu[i, t], col = "red", lwd = 2)
  }
}
################################################################
par(mfrow = c(3, 2))  # Set up the plotting window with 3 rows and 2 columns

for (t in 1:(true.occasions)) {
  for (i in 1:female.n.wmu) {
    # Construct the variable name dynamically for 'female.h.ad.wmu[t, i]'
    samp <- paste0("female.s.ad.wmu[", i, ", ", t, "]")
    
    # Plot trace of the sampled values
    matplot(samples_female[, samp], type = "l", 
            main = paste0("female.s.ad.wmu[", i, ", ", t, "]"), 
            ylab = "Sample values")
    
    # Add a reference line (e.g., rr[t, i] for the current time and WMU)
    abline(h = female.s.ad.wmu[i, t], col = "red", lwd = 2)
    
    # Plot density of the sampled values
    plot(density(samples_female[, samp]), main = "Density",
         xlab = "Beta Coefficient", ylab = "Density")
    
    # Add a vertical line at the true value
    abline(v = female.s.ad.wmu[i, t], col = "red", lwd = 2)
  }
}
###################################################################X
par(mfrow = c(3, 2))  # Set up the plotting window with 3 rows and 2 columns

for (t in 1:(true.occasions)) {
  for (i in 1:female.n.wmu) {
    # Construct the variable name dynamically for 'female.h.ad.wmu[t, i]'
    samp <- paste0("female.s.ad.wmu[", i, ", ", t, "]")
    
    # Plot trace of the sampled values
    matplot(samples_female[, samp], type = "l", 
            main = paste0("female.s.ad.wmu[", i, ", ", t, "]"), 
            ylab = "Sample values")
    
    # Add a reference line (e.g., rr[t, i] for the current time and WMU)
    abline(h = female.s.ad.wmu[i, t], col = "red", lwd = 2)
    
    # Plot density of the sampled values
    plot(density(samples_female[, samp]), main = "Density",
         xlab = "Beta Coefficient", ylab = "Density")
    
    # Add a vertical line at the true value
    abline(v = female.s.ad.wmu[i, t], col = "red", lwd = 2)
  }
}



