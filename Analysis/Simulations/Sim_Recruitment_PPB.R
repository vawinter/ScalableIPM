#----------------------------------------------------------#
# Simulating PPB (Poults Per Brood) data in base R
# Author: VAW
# Date: 1/6/2025
#----------------------------------------------------------#

# Set seed for reproducibility
set.seed(1235)

#----------------------------------------------------------#
# Parameters and Data Simulation Setup
#----------------------------------------------------------#

# Number of observations and groups
ph.N <- 1000  # Total number of observations
ph.J <- 10    # Number of groups (e.g., Wildlife Management Units)

#----------------------------------------------------------#
# Simulating Covariates
#----------------------------------------------------------#

# Simulate binary year indicators
ph.Year2019 <- rbinom(ph.N, 1, 0.2)  # Year 2019 indicator
ph.Year2020 <- rbinom(ph.N, 1, 0.2)  # Year 2020 indicator
ph.Year2021 <- rbinom(ph.N, 1, 0.2)  # Year 2021 indicator
ph.Year2022 <- rbinom(ph.N, 1, 0.2)  # Year 2022 indicator
ph.Year2023 <- rbinom(ph.N, 1, 0.2)  # Year 2023 indicator

# Simulate day-of-year (DOY) covariate and transformations
ph.doy.scale <- runif(ph.N, -1, 1)   # Scaled day-of-year
ph.doy.2 <- ph.doy.scale^2           # Quadratic term for DOY

# Assign each observation to a random group (e.g., WMU)
ph.wmu <- sample(1:ph.J, ph.N, replace = TRUE)

#----------------------------------------------------------#
# Simulating Priors
#----------------------------------------------------------#

# Coefficients for the regression model
ph.beta1 <- rnorm(1, 0, 0.01)  # Intercept
ph.beta2 <- rnorm(1, 0, 0.01)  # Year 2020 effect
ph.beta3 <- rnorm(1, 0, 0.01)  # Year 2021 effect
ph.beta4 <- rnorm(1, 0, 0.01)  # Year 2022 effect
ph.beta5 <- rnorm(1, 0, 0.01)  # Year 2023 effect
ph.beta6 <- rnorm(1, 0, 0.01)  # Linear DOY effect
ph.beta7 <- rnorm(1, 0, 0.01)  # Quadratic DOY effect

# Variance and random effect priors
ph.sigma <- rgamma(1, 1, 1)        # Gamma-distributed variance
ph.sigma.u <- runif(1, 0, 1)       # Random effect variance
ph.u <- rnorm(ph.J, 0, ph.sigma.u) # Random effects for each group

#----------------------------------------------------------#
# Calculating Mean and Gamma Parameters
#----------------------------------------------------------#

# Mean poults per brood for each observation
ph.mu <- exp(
  ph.beta1 +
    ph.beta2 * ph.Year2020 +
    ph.beta3 * ph.Year2021 +
    ph.beta4 * ph.Year2022 +
    ph.beta5 * ph.Year2023 +
    ph.beta6 * ph.doy.scale +
    ph.beta7 * ph.doy.2 +
    ph.u[ph.wmu]
)

# Gamma distribution parameters
ph.alpha <- (ph.mu^2) / (ph.sigma^2)  # Shape parameter
ph.theta <- (ph.sigma^2) / ph.mu      # Scale parameter

#----------------------------------------------------------#
# Simulating PPB Data
#----------------------------------------------------------#

# Simulate poults per brood ratio using Gamma distribution
PHratio <- rgamma(ph.N, shape = ph.alpha, scale = ph.theta)

#----------------------------------------------------------#
# Derived Variables for Aug 31 (DOY 244)
#----------------------------------------------------------#

# Aug 31 (DOY 244) DOY adjustment
ppb.aug31 <- (244 - ph.doy.scale) / ph.doy.scale
ppb.aug31.2 <- ppb.aug31^2

#----------------------------------------------------------#
# Summary and Validation
#----------------------------------------------------------#

# # Print summary of simulated PPB data
# cat("Summary of Simulated PPB Data\n")
# print(summary(PHratio))  # Summary of poults per brood ratios
# print(summary(ph.mu))    # Summary of mean poults per brood
# print(summary(ph.u))     # Summary of random effects (WMU-level)
# 
# # Verify derived variables
# cat("Sample Derived Values (Aug 31)\n")
# print(head(ppb.aug31))
# print(head(ppb.aug31.2))
