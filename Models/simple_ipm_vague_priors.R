ipm_simple <- nimbleCode({
  ###########################################################X
  # Simple IPM for modeling statewide WMU group dynamics for PA ----
  ###########################################################X
  # Recruitment: Hen with brood (HWB) model ----
  ###########################################################X
  #----------------------------------------------------------#
  # HWB: Likelihood ----
  #----------------------------------------------------------#
  ###########################################################X

    for (i in 1:hwb.N) {
      # Hens with brood
      HWB[i] ~ dbern(hwb.p[i])
      logit(hwb.p[i]) <- hwb.beta1 + hwb.beta2 * hwb.Year2020[i] + hwb.beta3 * hwb.Year2021[i] +
        hwb.beta4 * hwb.Year2022[i] + hwb.beta5 * hwb.Year2023[i] + hwb.beta6 * hwb.doy.scale[i] +
        hwb.beta7 * hwb.doy.2[i] + hwb.u[hwb.wmu[i]]
    } # N
  
    #----------------------------------------------------------#
    # HWB: Priors ----
    #----------------------------------------------------------#
  
    # Intercept
    hwb.beta1 ~ dnorm(0, sd = 1)
    # year 2019 effect
    hwb.beta2 ~ dnorm(0, sd = 1)
    # year 2020 effect
    hwb.beta3 ~ dnorm(0, sd = 1)
    # year 2021 effect
    hwb.beta4 ~ dnorm(0, sd = 1)
    # year 2022 effect
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
  
    #----------------------------------------------------------#
    # Calculate the expected number of females on Aug 31 (Sept 1) ----
    #----------------------------------------------------------#
  
    # Loop over the years I am estimating abundance (4)
    for (t in 1:Nyears) {
      for (u in 1:female.n.wmu) {
        # Derived estimates for the number of hens on August 31 per WMU
        aug31.hwb[t,u] <- expit(
          hwb.beta1 + hwb.beta2 * (t == 1) +  hwb.beta3 * (t == 2) +
            hwb.beta4 * (t == 3) +  hwb.beta5 * (t == 4) +
            hwb.beta6 * hwb.aug31 + hwb.beta7 * hwb.aug31.2 +
            hwb.u[u]
        )
      } # end u
    } # end t
  # ###########################################################X
  # Recruitment: Poults per brood (PPB) model ----
  ###########################################################X
  #----------------------------------------------------------#
  # PPB: Likelihood ----
  #----------------------------------------------------------#
  ###########################################################X
  
  for (i in 1:ph.N) {
  
    #----------------------------------------------------------#
    # Note: For a gamma distribution, the variance is the shape times squared
    # scale, which is equivalent to the mean (shape times scale) times the scale.
    # So the variance scales with the mean (ph.mu).
    #----------------------------------------------------------#
  
    PHratio[i] ~ dgamma(shape = ph.alpha[i], scale = ph.theta[i])
  
    # Calculate the shape and scale parameter
    ph.alpha[i] <- ph.mu[i] * ph.disp
    ph.theta[i] <- 1/ph.disp
  
    # Calculate the mean for each observation
    log(ph.mu[i]) <- ph.beta1 + ph.beta2 * ph.Year2020[i] +
      ph.beta3 * ph.Year2021[i] + ph.beta4 * ph.Year2022[i] +
      ph.beta5 * ph.Year2023[i] + ph.beta6 * ph.doy.scale[i] +
      ph.beta7 * ph.doy.2[i] + ph.u[ph.wmu[i]]
  }

  #----------------------------------------------------------#
  # PPB: Priors  ----
  #----------------------------------------------------------#

  # Using vague priors on betas
  # Intercept
  ph.beta1 ~ dnorm(0, sd = 1)
  # year 2019 effect
  ph.beta2 ~ dnorm(0, sd = 1)
  # year 2020 effect
  ph.beta3 ~ dnorm(0, sd = 1)
  # year 2021 effect
  ph.beta4 ~ dnorm(0, sd = 1)
  # year 2022 effect
  ph.beta5 ~ dnorm(0, sd = 1)
  # DOY: make conformable for matrix multiplication
  # day of year
  ph.beta6 ~ dnorm(0, sd = 1)
  # doy^2
  ph.beta7 ~ dnorm(0, sd = 1)

  # Scaling parameter for Gamma distribution
  ph.disp ~ dunif(0, 1)

  # Random effect of wmu
  ph.sigma.u ~ dunif(0, 1)
  for (j in 1:ph.J) {
    ph.u[j] ~ dnorm(0, sd = ph.sigma.u)
  }

  #----------------------------------------------------------#
  # Calculate the expected ratio of poults on Aug 31 (Sept 1) ----
  #----------------------------------------------------------#

  # Loop over the years I am estimating abundance (4)
  for (t in 1:Nyears) {
    # Loop over management units
    for (u in 1:female.n.wmu) {
      # Derived estimates for the number of ppb on August 31 per WMU and year
      aug31.ppb[t,u] <- exp(
        ph.beta1 + ph.beta2 * (t == 1) + ph.beta3 * (t == 2) +
          ph.beta4 * (t == 3) + ph.beta5 * (t == 4) + ph.beta6 * ppb.aug31 +
          ph.beta7 * ppb.aug31.2 + ph.u[u]
      )
    } # end u
  } # end t

  ###########################################################X
  # Female survival  ----
  ###########################################################X
  # Estimate a survival rate for females based off full IPM 
  # known-fate model 
  for (t in 1:Nyears) {
  for (u in 1:female.n.wmu) {
    # avg.juv.s.kf[t, u] ~ dbeta(shape1 = 4, shape2 = 20)
    # avg.ad.s.kf[t, u] ~ dbeta(shape1 = 10, shape2 = 5)
    
    avg.juv.s.kf[t, u] ~ dbeta(shape1 = 1, shape2 = 1)
    avg.ad.s.kf[t, u] ~ dbeta(shape1 = 1, shape2 = 1)
 }
}

  ###########################################################X
  #----------------------------------------------------------#
  # Dead Recovery model for each males ----
  ###########################################################X
  # DRM: Female model based off males ----
  ###########################################################X
  #----------------------------------------------------------#
  # DRM: WMU harvest rates
  #----------------------------------------------------------#
  # Note: female.n.occasions-1 is done to match the notation for
  # my male DRM. See note in the male DRM likelihood for more info.
  #----------------------------------------------------------# 
  for (t in 1:(female.n.occasions-1)) {
    for (u in 1:female.n.wmu) {
      # For a vague prior
      female.h.ad.wmu[t,u] ~ dbeta(1, 1)
      female.h.juv.wmu[t,u] ~ dbeta(1, 1)
 
      # # Use a Beta prior on harvest rate
      # female.h.ad.wmu[t,u] ~ dbeta(shape1 = 2, shape2 = 50)
      # female.h.juv.wmu[t,u] ~ dbeta(shape1 = 2, shape2 = 50)
    }#u
  }#t

  ###########################################################X
  # DRM: Male model ----
  ###########################################################X
  # Dead Recovery model: Priors on age and wmu
  #----------------------------------------------------------#
  ###########################################################X

  # Reporting rates
  ### Informative priors for hunter reporting rates based on Diefenbach et al. (2012)
  male.rrate.j ~ dnorm(0.71, sd = 0.072) #  non-reward bands for juveniles
  male.rrate.a ~ dnorm(0.87, sd = 0.039) #  non-reward bands for adults

  #----------------------------------------------------------#
  # Proportion mortality due to hunting (Seber parameterization: (1-s)r)
  male.seber.recov[1] ~ dunif(0, 1)   # Prior for proportion mortality due to hunting - juveniles
  male.seber.recov[2] ~ dunif(0, 1)    # Prior for proportion of mortality due to hunting - adults

  #----------------------------------------------------------#

  # Prior on juvenile effect on survival
  male.juvenile.effect ~ dnorm(0, sd = 0.5)

  #----------------------------------------------------------#

  # Prior on time effect
  for (t in 1:(male.n.occasions-1)) {
    male.time.effect[t] ~ dnorm(0, sd = 0.5)
  }

  #----------------------------------------------------------#

  ## Random effect of wmu
  male.sigma ~ dunif(0,10)
  for (i in 1:10) {
    male.wmu.effect[i] ~ dnorm(0, sd = male.sigma)
  }

  #----------------------------------------------------------#

  # Survival and Recovery rate
  for (i in 1:male.nind) {
    for (t in male.f[i]:(male.n.occasions-1)) {

      #----------------------------------------------------------#
      # Note:
      # juvenile.effect is age (juv) coefficient, time.effect is time coefficient
      #----------------------------------------------------------#

      logit(male.s[i,t]) <- male.juvenile.effect*male.I[i,t] +
        # Hate this hard coding - need to find nimble work-around
        inprod(male.time.effect[1:4], male.time.param[t, 1:4]) +
        male.wmu.effect[male.wmu[i]]

      #----------------------------------------------------------#
      # Note:
      # the four terms in r[i,t] represent: juvenile non-reward bands,
      # juvenile reward bands, adult non-reward bands, adult reward band
      # the indicators I[] (juvenile=1) and II[] (reward = 0)
      # simply retain/cancel the appropriate terms
      #----------------------------------------------------------#
      # Recovery rate w. Seber parameterization: ((1-s)r)
      male.r[i,t] <- male.seber.recov[1] * male.I[i, t] * male.rrate.j * male.II[i] +
        male.seber.recov[1] * male.I[i, t] * (1 - male.II[i]) + male.seber.recov[2] * (1 - male.I[i, t]) * male.rrate.a * male.II[i] +
        male.seber.recov[2] * (1 - male.I[i, t]) * (1 - male.II[i])

      #----------------------------------------------------------#

    } #t
  } #i

  ###########################################################X
  # DRM: Derived Estimates
  ###########################################################X
  #----------------------------------------------------------#
  # DRM: Derived survival rate estimates
  #----------------------------------------------------------#

  for (t in 1:(male.n.occasions-1)) {
    logit(male.mean.s.ad[t]) <- male.time.effect[t]
    logit(male.mean.s.jv[t]) <- male.time.effect[t] + male.juvenile.effect

    #----------------------------------------------------------#
    # DRM: Derived harvest rate estimates
    #----------------------------------------------------------#

    male.mean.harv.jv[t] <- (1-male.mean.s.jv[t])*male.seber.recov[1]
    male.mean.harv.ad[t] <- (1-male.mean.s.ad[t])*male.seber.recov[2]
  }

  #----------------------------------------------------------#
  # DRM: WMU survival rates
  #----------------------------------------------------------#

  for (t in 1:(male.n.occasions-1)) {
    for (u in 1:male.n.wmu) {
      logit(male.s.ad.wmu[t,u]) <-  inprod(male.time.effect[1:4], male.time.param[t, 1:4]) + male.wmu.effect[u]
      logit(male.s.juv.wmu[t,u]) <- male.juvenile.effect +
        inprod(male.time.effect[1:4], male.time.param[t, 1:4]) + male.wmu.effect[u]
    }#u
  }#t

  #----------------------------------------------------------#
  # DRM: WMU harvest rates
  #----------------------------------------------------------#

  for (t in 1:(male.n.occasions-1)) {
    for (u in 1:male.n.wmu) {
      male.h.juv.wmu[t,u] <- (1-male.s.juv.wmu[t,u])*male.seber.recov[1]
      male.h.ad.wmu[t,u] <- (1-male.s.ad.wmu[t,u])*male.seber.recov[2]
    }# u
  }#t

  #----------------------------------------------------------#
  # DRM: Likelihood (Males) ----
  #----------------------------------------------------------#
  # Note: male.n.occasion here is 5, which is set so that my first year for
  # z is 1, when all my initial individuals were alive.
  # For more clarification, see Kery and Schaub:
  # "Bayesian Population Analysis using WinBUGS: A Hierarchical Perspective."
  #----------------------------------------------------------#

  for (i in 1:male.nind){
    # Define latent state at first capture
    male.z[i,male.f[i]] <- 1
    for (t in (male.f[i]+1):male.n.occasions){ # true.occasions + 1 b.c they're dead

      #------------------------#

      # State process
      male.z[i,t] ~ dbern(male.mu1[i,t])
      male.mu1[i,t] <- male.s[i,t-1] * male.z[i,t-1]
      #------------------------#

      # Observation process
      male.y[i,t] ~ dbern(male.mu2[i,t])
      male.mu2[i,t] <- male.r[i,t-1] * (male.z[i,t-1] - male.z[i,t])
    } #t
  } #i
  
  ##############################################################X
  # Abundance: Derived estimates via Lincoln-Peterson Est. ----
  ##############################################################X
  #
  # Abundance (N) is modeled during hunting season (spring, fall) annually, t,
  # for each WMU, u.
  #
  ##-------------------------------------------------------------------------#X
  # Estimating population abundance by age class and season for subsequent
  # occasions per WMU
  #---------------------------------------------------------------------#X
  # Male Abundance (May t -> May t+1)
  #---------------------------------------------------------------------#X

  for (u in 1:male.n.wmu) {
    # start: 2020 (year 1)
    # Initial abundance for adult males in each WMU using 2020 hr
    N.lambda.ad.male[u] <- (th.year1.male.ad[u]) / male.h.ad.wmu[1, u]
    male.N.ad[1, u] ~ dpois(N.lambda.ad.male[u])

    # Initial abundance for juvenile males in each WMU using 2020 hr
    N.lambda.juv.male[u] <- (th.year1.male.juv[u]) / male.h.juv.wmu[1, u]
    male.N.juv[1, u] ~ dpois(N.lambda.juv.male[u])
  } # end u

  # Loop over time occasions (2021-2023)
  for (t in 2:Nyears) {
    for (u in 1:male.n.wmu) {
      # start: 2021 (year 2)
      # Number of surviving adult males from time t-1 to t (N1S1)
      male.N.ad.Survived[t, u] ~ dbin(prob = male.s.ad.wmu[t-1, u],
                                      size = male.N.ad[t-1, u])

      # Number of surviving juvenile males from time t-1 to t (N1S1)
      male.N.juv.Survived[t, u] ~ dbin(prob = male.s.juv.wmu[t-1, u],
                                       size = male.N.juv[t-1, u])

      # Adults who survived, and juveniles who have aged up to adult and survived
      # Note: Poisson is to keep numbers whole for H/hr
      male.N.ad[t, u]  <- (male.N.ad.Survived[t, u] + male.N.juv.Survived[t, u])

      # New juveniles who entered the population
      # For males, these juveniles are from Nov (t-1) entering in May (t)
      male.N.juv[t, u] ~ dpois(recruitment[t-1, u])

      # Lincoln-Petersen: Adult males in spring
      harvest.ad.spring[t, u] ~ dbin(prob = male.h.ad.wmu[t, u],
                                     size = male.N.ad[t, u])

      # Lincoln-Petersen: Juvenile males in spring
      harvest.juv.spring[t, u] ~ dbin(prob = male.h.juv.wmu[t, u],
                                      size = male.N.juv[t, u])

    } # end t
  } # end u

  #---------------------------------------------------------------------#X
  # Female Abundance (November t -> November t+1)
  #---------------------------------------------------------------------#X
  
  for (u in 1:female.n.wmu) {
    # start: 2020
    # Initial abundance for adult females in each WMU using 2020 hr
    N.lambda.ad.female[u] <- (th.year1.female.ad[u]) / female.h.ad.wmu[1, u] 
    female.N.ad[1, u] ~ dpois(N.lambda.ad.female[u])
    
    # Initial abundance for juvenile females in each WMU using 2020 hr
    N.lambda.juv.female[u] <- (th.year1.female.juv[u]) / female.h.juv.wmu[1, u]
    female.N.juv[1, u] ~ dpois(N.lambda.juv.female[u])
    
    # Recruitment for the time occasion (2020)
    recruitment[1, u] <- ((female.N.ad[1, u] * aug31.hwb[1, u]) * aug31.ppb[1, u])/2
    
  } # end u
  
  # Loop over time occasions (2021-2023)
  for (t in 2:Nyears) { 
    for (u in 1:female.n.wmu) {
      # start: 2021
      # Recruitment for the time occasion (2021 - 2023)
      recruitment[t, u] <- ((female.N.ad[t, u] * aug31.hwb[t, u]) * aug31.ppb[t, u])/2
      
      # Number of surviving adult females from time t-1 to t (N1S1)
      female.N.ad.Survived[t, u] ~ dbin(size = female.N.ad[t-1, u],
                                        prob = avg.ad.s.kf[t-1, u])
      
      
      # Number of surviving adult females from time t-1 to t (N1S1)
      female.N.juv.Survived[t, u] ~ dbin(size = female.N.juv[t-1, u],
                                         prob = avg.juv.s.kf[t-1, u])
      
      # Adults who survived, and juveniles who have aged up to adult and survived
      # Note: Poisson is to keep numbers whole for H/hr
      female.N.ad[t, u]  <- (female.N.ad.Survived[t, u] + female.N.juv.Survived[t, u])
      
      # New juveniles who entered the population
      # For females, these juveniles are from Nov (t) entering in Nov (t)
      female.N.juv[t, u] ~ dpois(recruitment[t, u])
      
      
      # Lincoln-Petersen: Adult females in fall
      harvest.ad.fall[t, u] ~ dbin(prob = female.h.ad.wmu[t, u],
                                   size = female.N.ad[t, u])
      
      # Lincoln-Petersen: Juvenile females in fall
      harvest.juv.fall[t, u] ~ dbin(prob = female.h.juv.wmu[t, u],
                                    size = female.N.juv[t, u])
      
      
    } # end t
  } # end u
  
})
