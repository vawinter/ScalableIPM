#----------------------------------------------------------#
# Simulating HWB data in base R
# Author: VAW
# Date: 1/6/2025
#----------------------------------------------------------#

# Set seed for reproducibility
set.seed(1235)

#----------------------------------------------------------#
# Parameters and Data Simulation Setup
#----------------------------------------------------------#

# Number of observations and groups
hwb.N <- 10000  # Total number of observations
hwb.J <- 10     # Number of groups (e.g., Wildlife Management Units)

#----------------------------------------------------------#
# Simulating Covariates
#----------------------------------------------------------#

# Simulate binary year indicators
hwb.Year2019 <- rbinom(hwb.N, 1, 0.2)  # Year 2019 indicator
hwb.Year2020 <- rbinom(hwb.N, 1, 0.2)  # Year 2020 indicator
hwb.Year2021 <- rbinom(hwb.N, 1, 0.2)  # Year 2021 indicator
hwb.Year2022 <- rbinom(hwb.N, 1, 0.2)  # Year 2022 indicator
hwb.Year2023 <- rbinom(hwb.N, 1, 0.2)  # Year 2023 indicator

# Simulate day-of-year (DOY) covariate and transformations
hwb.doy.scale <- runif(hwb.N, -1, 1)   # Scaled day-of-year
hwb.doy.2 <- hwb.doy.scale^2           # Quadratic term for DOY

# Assign each observation to a random group (e.g., WMU)
hwb.wmu <- sample(1:hwb.J, hwb.N, replace = TRUE)

#----------------------------------------------------------#
# Simulating Priors
#----------------------------------------------------------#

# Coefficients for the logistic regression model
hwb.beta1 <- rnorm(1, 0, 0.01)  # Intercept
hwb.beta2 <- rnorm(1, 0, 0.01)  # Year 2020 effect
hwb.beta3 <- rnorm(1, 0, 0.01)  # Year 2021 effect
hwb.beta4 <- rnorm(1, 0, 0.01)  # Year 2022 effect
hwb.beta5 <- rnorm(1, 0, 0.01)  # Year 2023 effect
hwb.beta6 <- rnorm(1, 0, 0.01)  # Linear DOY effect
hwb.beta7 <- rnorm(1, 0, 0.01)  # Quadratic DOY effect

# Random effect standard deviation
hwb.sigma <- runif(1, 0, 10)

# Random effects for each group (e.g., WMUs)
hwb.u <- rnorm(hwb.J, 0, hwb.sigma)

#----------------------------------------------------------#
# Calculating Mean Probability for HWB
#----------------------------------------------------------#

# Logistic regression model (logit scale converted to probability)
hwb.p <- expit(
  hwb.beta1 +
    hwb.beta2 * hwb.Year2020 +
    hwb.beta3 * hwb.Year2021 +
    hwb.beta4 * hwb.Year2022 +
    hwb.beta5 * hwb.Year2023 +
    hwb.beta6 * hwb.doy.scale +
    hwb.beta7 * hwb.doy.2 +
    hwb.u[hwb.wmu]
)

#----------------------------------------------------------#
# Simulating HWB Data
#----------------------------------------------------------#

# Simulate HWB observations using the Bernoulli distribution
HWB <- rbinom(hwb.N, 1, hwb.p)

#----------------------------------------------------------#
# Derived Variables for Aug 31 (DOY 244)
#----------------------------------------------------------#

# Aug 31 (DOY 244) DOY adjustment
hwb.aug31 <- (244 - hwb.doy.scale) / hwb.doy.scale
hwb.aug31.2 <- hwb.aug31^2

#----------------------------------------------------------#
# Summary and Validation
#----------------------------------------------------------#

# # Print summary of simulated HWB data
# cat("Summary of Simulated HWB Data\n")
# print(table(HWB))  # Count of 0s and 1s in HWB
# print(summary(hwb.p))  # Summary of HWB probabilities
# print(summary(hwb.u))  # Summary of random effects (WMU-level)
# 
# # Verify that derived variables align with expectations
# cat("Sample Derived Values (Aug 31)\n")
# print(head(hwb.aug31))
# print(head(hwb.aug31.2))

