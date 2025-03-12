###############################################################################X
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#################### Integrated Population Model (IPM): #######################X
#                 #---# PhD Dissertation: Chapter 1 #---#
#        Creating a Bayesian IPM to inform turkey management in PA
###                   *** SIMULATION TEST ***                               ###X
###                                                                         ###X
#    Modeling: Creating a  known-fate model for telemetered hens by 
# wildlife management unit (WMU), and age class over months prior to
#                            fall hunting season.
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #X
# 
# Created by: Veronica A. Winter
# Last edited: 04/23/2024
#
# Added cumprod test 5/31/2024
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #X
###############################################################################X

# clean env
rm(list = ls())
gc()

# Source scripts of functions and data preparation
source("Analysis/Scripts/00_IPM_funs.R")

# Function to calculate survival probability using inverse cloglog in base r
icll <- function(x) {
  1 - exp(-exp(x))
}

# Set seed for reproducibility
set.seed(0219)

plot(density(icll(rnorm(10000, 0, 0.0001))))
plot(density(icll(rnorm(10000, 0, 0.001))))
plot(density(icll(rnorm(10000, 0, 0.01))))
plot(density(icll(rnorm(10000, 0, 0.1))))
plot(density(icll(rnorm(10000, 0, 1))))
plot(density(icll(rnorm(10000, 0, 10))))
plot(density(icll(rnorm(10000, 0, 100))))
dev.off()


#############################################################################X
# Simulating data: checking priors and model fit ----
#############################################################################X
## Parameters
nind <- 1000  # Number of individuals
nMonths <- 48  # Number of months
beta_age <- c(0.3, 0.5)  # Effects of being adult or juvenile
beta_wmu <- c(0.2, 0.3, 0.15, 0.5)  # Effect of different WMUs
survival_status <- matrix(NA, nrow = nind, ncol = nMonths)

# Simulate age groups and WMU assignments
is_adult <- matrix(sample(0:1, nind * nMonths, replace = TRUE), nrow = nind, ncol = nMonths)
is_juvenile <- 1 - is_adult  

# Fill in WMU matrix and include intercept
wmu<- matrix(0, nind*4, nrow = nind, ncol = 4)
wmu[1:125, 1] <- 1
wmu[125:250, 2] <- 1
wmu[251:375, 3] <- 1
wmu[, 4] <- 1 # intercept

# Initialize first and last matrices 
first <- sample(1:4, nind, replace = TRUE)
last <- NULL

# Initialize matrices to store survival data for each age class
survival_adult <- NULL
survival_juvenile <- NULL
cll_survival <-  NULL
survival_wmu <-  NULL

#############################################################################X
## Simulating data: simulating status and survival probabilities ----
#############################################################################X
for(i in 1:nind) {
  # setting so I can make sure that it gets 'overwritten'
  survival_status[i, first[i]] <- 9
  
  for (t in first[i]:nMonths) {
    # This is my estimation: intercepts for 'adult' and '4D'
    cll_survival[i] <- (beta_age[1] + beta_wmu[4]) + beta_age[2] * is_juvenile[i, t] + (beta_wmu[1:3]%*%wmu[i, 1:3])
    
    # this is creating my survival status, because they need to come from a DISTRIBUTION
    survival_status[i, t] <- rbern(n=1, prob = icll(cll_survival[i]))
    
    # Check if individual has 'died'
    if (survival_status[i, t] == 0) {
      last[i] <- t
      if (t < nMonths) {
        survival_status[i, (t+1):nMonths] <- NA
      }
      break  # Exit the loop for this individual
    }
    
    # Fill age-specific survival matrices
    if (is_adult[i, t] == 1) {
      survival_adult[t] <- cumprod(icll(cll_survival[i]))
    } else if (is_juvenile[i, t] == 1) {
      survival_juvenile[t] <- cumprod(icll(cll_survival[i]))
    }
    
  } # end t
} # end i

#############################################################################X
### Simulating data: Check model outputs ----
#############################################################################X
survival_adult[1:8]
survival_juvenile[1:8]

# Fill in NAs in 'last' with the max month
last <- ifelse(is.na(last),max(nMonths), last)
length(last) # Check dimensions
# last[1000] <- nMonths # add if dim < nind


# Note: dbern in nimble says 'this is a pdf for max likelihood' and
# dbern in Bayesian says 'this arises from a distribution' 

#############################################################################X
# Create Nimble model ----
#############################################################################X
survival <- nimbleCode({
  #----------------------------------------------------------#
  # Known-fate model: Priors on age and wmu
  #----------------------------------------------------------#
  # Intercept: combination of adult and wmu4 means
  beta.int[1] ~ dnorm(0, sd = 0.5)
  
  # Age
  beta.age[1] ~ dnorm(0, sd = 0.5)
  
  # Note: this is currently a fixed effect
  # to make a random effect, I would need to put a prior
  # on the sd: ex sigma.wmu ~ dinvgamma(0.01, 0.01)  # Weakly informative prior
  # # WMU random effects
  # for(u in 1:3){
  #   beta.wmu[u] ~ dnorm(0, sd = sqrt(sigma.wmu))
  # }
  
  # WMU
  for(u in 1:3){
    beta.wmu[u] ~ dnorm(0, sd = 1)
  }
  
  #----------------------------------------------------------#
  # Known fate model: Likelihood
  #----------------------------------------------------------#
  for (i in 1:nind) {
    # first to last encounter of each individual
    for(t in first[i]:last[i]){
      # Derived individual survival rate with inverse cloglog
      kf.survival[i, t] <-  beta.int[1] + beta.age[1] * is.juvenile[i, t] + inprod(beta.wmu[1:3], wmu[i, 1:3])
      
      # Likelihood for individual capture histories
      status[i, t] ~ dbern(prob = icloglog(kf.survival[i, t]))
      
      # Calculate cumulative survival probabilities
      if (t == first[i]) {
        telem.s.kf[i, t] <- icloglog(kf.survival[i, t])
      } else {
        telem.s.kf[i, t] <- telem.s.kf[i, t-1] * icloglog(kf.survival[i, t])
      }
      
    } # end t
  } # end i
})

#############################################################X
##           Model run set up ----
#############################################################X
# Data
test_dat <- list(
  status = survival_status, # 1 = alive, 0 = died, NA = not in study
  # is.adult = is_adult,
  is.juvenile = is_juvenile,
  wmu = wmu
)

# Constants
consts <- list(
  nind = nind,
  first = first,
  last = last
)

inits <- list(
  beta.int = beta_age[1],
  beta.age = beta_age[2],
  beta.wmu = beta_wmu[1:3]
)

# Create the Nimble model
model <- nimbleModel(
  code = survival,
  data = test_dat,
  constants = consts,
  inits = inits
)

# Check initialization issues
model$initializeInfo()
gc()

# Set MCMC specifications
burn <- 1000
iter <-  10000
t <- 1
chain <- 1

# Run MCMC 
survival_results <- nimbleMCMC(
  model = model,
  monitors = c("beta.wmu", 
               "beta.int",
               "beta.age",
               "telem.s.kf"
  ),
  niter = iter, 
  nburnin = burn, 
  thin = t, 
  nchains = chain, 
  setSeed = FALSE, 
  samplesAsCodaMCMC = TRUE,
  WAIC = TRUE)
beepr::beep(1)

# Print results
# print(survival_results)

#############################################################X
##           Model diagnostics ----
#############################################################X
### Look at parameter recovery ----
MCMCpstr(survival_results$samples, params = c("beta.int"))
# [1] 0.8391951
beta_age[1] + beta_wmu[4]
# [1] 0.8

MCMCpstr(survival_results$samples, params = c("beta.age"))
#   0.5152433

beta_age[2]
# [1] 0.50

MCMCpstr(survival_results$samples, params = c("beta.wmu"))
# [1] 0.2218685 0.4133110 0.2047299

beta_wmu[1:3]
# [1]  0.20 0.30 0.15 

MCMCpstr(survival_results$samples, params = c("telem.s.kf"))

### Check fit ----
#### Look at parameter recovery (density and chains) ----
MCMCtrace(survival_results$samples, pdf = F, iter = iter, params = c("beta.age"))
MCMCtrace(survival_results$samples, pdf = F, iter = iter, params = c("beta.wmu"))
MCMCtrace(survival_results$samples, pdf = F, iter = iter, params = c("beta.int"))

## find the effective sample size.
## If it is too low, increase niter above and re-run MCMC
coda::effectiveSize(survival_results)

##############################
## optional - Plot MCMC chains
##############################
# Intercept
par(mfrow = c(1,2))
matplot(survival_results$samples[,"beta.int[1]"],type="l", main = "Intercept")
abline(h = beta_age[1]+beta_wmu[4], col = "red", lwd = 2)
plot(density(survival_results$sample[,"beta.int[1]"]), main = "Density of Beta Coefficient",
     xlab = "Beta Coefficient", ylab = "Density")
abline(v = beta_age[1]+beta_wmu[4], col = "red", lwd = 2)

# beta age
par(mfrow = c(1,2))
matplot(survival_results$samples[,"beta.age[1]"],type="l", main = "Juvenile")
abline(h = beta_age[2], col = "red", lwd = 2)
plot(density(survival_results$sample[,"beta.age[1]"]), main = "Density of Beta Coefficient",
     xlab = "Beta Coefficient", ylab = "Density")
abline(v = beta_age[2], col = "red", lwd = 2)

# wmu[1]
par(mfrow = c(3,2))
matplot(survival_results$samples[,"beta.wmu[1]"],type="l", main = "WMU 1")
abline(h = beta_wmu[1], col = "red", lwd = 2)
plot(density(survival_results$sample[,"beta.wmu[1]"]), main = "Density of Beta Coefficient",
     xlab = "Beta Coefficient", ylab = "Density")
abline(v = beta_wmu[1], col = "red", lwd = 2)
# wmu[2]
#par(mfrow = c(1,2))
matplot(survival_results$samples[,"beta.wmu[2]"],type="l", main = "WMU 2")
abline(h = beta_wmu[2], col = "red", lwd = 2)
plot(density(survival_results$sample[,"beta.wmu[2]"]), main = "Density of Beta Coefficient",
     xlab = "Beta Coefficient", ylab = "Density")
abline(v = beta_wmu[2], col = "red", lwd = 2)
# wmu[3]
#par(mfrow = c(1,2))
matplot(survival_results$samples[,"beta.wmu[3]"],type="l", main = "WMU 3")
abline(h = beta_wmu[3], col = "red", lwd = 2)
plot(density(survival_results$sample[,"beta.wmu[3]"]), main = "Density of Beta Coefficient",
     xlab = "Beta Coefficient", ylab = "Density")
abline(v = beta_wmu[3], col = "red", lwd = 2)

##############################
## Get posterior means and CIs
##############################

## find posterior means
apply(survival_results$samples,2,mean)

## find posterior 95% Credible Intervals
apply(survival_results$samples,2,quantile,c(.025,.975))


# # Extract survival probabilities ---
saveRDS(survival_results, "Data/Output/20240423_known_fate-SIM.RDS")
