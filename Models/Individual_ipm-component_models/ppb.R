# Define the nimble model
ppb_code <- nimbleCode({
  ###########################################################X
  # Recruitment: Poults per brood (PPB) model
  ###########################################################X
  #----------------------------------------------------------#
  # PPB: Likelihood ----
  #----------------------------------------------------------#
  ###########################################################X
  for (i in 1:ph.N) {
    # Note: For a gamma distribution, the variance is the shape times squared 
    # scale, which is equivalent to the mean (shape times scale) times the scale. 
    # So the variance scales with the mean (ph.mu).
    
    PHratio[i] ~ dgamma(shape = ph.alpha[i], scale = ph.theta[i])
    
    # Calculate the shape and scale parameter
    ph.alpha[i] <- ph.mu[i] * ph.disp
    ph.theta[i] <- 1/ph.disp
    
    # Calculate the mean for each observation
    log(ph.mu[i]) <- ph.beta1 * ph.Year2019[i] + ph.beta2 * ph.Year2020[i] + 
      ph.beta3 * ph.Year2021[i] + ph.beta4 * ph.Year2022[i] + 
      ph.beta5 * ph.Year2023[i] + ph.beta6 * ph.doy.scale[i] + 
      ph.beta7 * ph.doy.2[i] + ph.u[ph.wmu[i]]  
  }
  
  #----------------------------------------------------------#
  # Priors  ----
  #----------------------------------------------------------#
  # Using vague priors on betas
  # year 2019 effect
  ph.beta1 ~ dnorm(0, sd = 1)
  # year 2020 effect
  ph.beta2 ~ dnorm(0, sd = 1)
  # year 2021 effect
  ph.beta3 ~ dnorm(0, sd = 1)
  # year 2022 effect
  ph.beta4 ~ dnorm(0, sd = 1)
  # year 2023 effect
  ph.beta5 ~ dnorm(0, sd = 1)
  # day of year
  ph.beta6 ~ dnorm(0, sd = 1)
  # doy^2
  ph.beta7 ~ dnorm(0, sd = 1) 
  
  # Scaling parameter for Gamma distribution
  ph.disp ~ dunif(0, 10) 
  
  # Random effect of wmu
  ph.sigma.u ~ dunif(0, 1)
  for (j in 1:ph.J) {
    ph.u[j] ~ dnorm(0, sd = ph.sigma.u)
  }
  
  #----------------------------------------------------------#
  # Calculate Mean Over True Occasions ----
  #----------------------------------------------------------#
  # Loop over true occasions
  for (t in 1:female.n.occasions) {
    # Loop over management units
    for (u in 1:female.n.wmu) {
      # Derived estimates for the number of hens on August 31 per WMU and year
      aug31.ppb[t,u] <- exp(ph.beta1 * (t == 1) + ph.beta2 * (t == 2) + 
                              ph.beta3 * (t == 3) + ph.beta4 * (t == 4) + 
                              ph.beta5 * (t == 5) + ph.beta6 * aug31 + 
                              ph.beta7 * aug31.2 + ph.u[u])
    } # end u
  } # end t
  
})