# Define the nimble model
#----------------------------------------------------------#
hwb_code <- nimbleCode({
  ###########################################################X
  # Recruitment: Hen with brood (HWB) model
  ###########################################################X
  #----------------------------------------------------------#
  # HWB: Likelihood ----
  #----------------------------------------------------------#
  ###########################################################X
  for (i in 1:hwb.N) {
    # Hens with brood
    HWB[i] ~ dbern(hwb.p[i])
    logit(hwb.p[i]) <- hwb.beta1 * hwb.Year2019[i] + hwb.beta2 * hwb.Year2020[i] + hwb.beta3 * hwb.Year2021[i] + 
      hwb.beta4 * hwb.Year2022[i] + hwb.beta5 * hwb.Year2023[i] + hwb.beta6 * hwb.doy.scale[i] + 
      hwb.beta7 * hwb.doy.2[i] + hwb.u[hwb.wmu[i]]
  } # N
  
  #----------------------------------------------------------#
  # Priors on betas ----
  #----------------------------------------------------------#
  # year 2019 effect
  hwb.beta1 ~ dnorm(0, sd = 1)
  # year 2020 effect
  hwb.beta2 ~ dnorm(0, sd = 1)
  # year 2021 effect
  hwb.beta3 ~ dnorm(0, sd = 1)
  # year 2022 effect
  hwb.beta4 ~ dnorm(0, sd = 1)
  # year 2023 effect
  hwb.beta5 ~ dnorm(0, sd = 1)
  # DOY: make conformable for matrix multiplication
  # day of year
  hwb.beta6 ~ dnorm(0, sd = 1)
  # doy^2
  hwb.beta7 ~ dnorm(0, sd = 1) 
  # Random effect of wmu
  hwb.sigma ~ dunif(0, 10)
  for (j in 1:hwb.J) {
    hwb.u[j] ~ dnorm(0, sd = hwb.sigma)
  } # J
  
  
  # Loop over true occasions
  for (t in 1:5) {
    # Loop over management units
    for (u in 1:4) { 
      # Derived estimates for the number of hens on August 31 per WMU
      aug31.hwb[t,u] <- expit(hwb.beta1 * (t == 1) + (hwb.beta1 * (t == 1) + hwb.beta2 * (t == 2)) + 
                                (hwb.beta1 * (t == 1) + hwb.beta3 * (t == 3)) + (hwb.beta1 * (t == 1) + hwb.beta4 * (t == 4)) + 
                                (hwb.beta1 * (t == 1) + hwb.beta5 * (t == 5)) + (hwb.beta1 * (t == 1) + hwb.beta6 * hwb.aug31) + 
                                hwb.beta7 * hwb.aug31.2 + hwb.u[u])
    } # end u
  } # end t
  
})