ipm <- nimbleCode({
  ###########################################################X
  # Full IPM for modeling dynamics in 3 study site WMUs in PA ----
  ###########################################################X
  # Recruitment: Hen with brood (HWB) model ----
  ###########################################################X
  #----------------------------------------------------------#
  # HWB: Likelihood ----
  #----------------------------------------------------------#
  ###########################################################X
  # Note: intercept is 2019
  #----------------------------------------------------------#
  for (i in 1:hwb.N) {
    # Hens with brood
    HWB[i] ~ dbern(hwb.p[i])
    logit(hwb.p[i]) <- hwb.beta1  + hwb.beta2 * hwb.Year2020[i] + hwb.beta3 * hwb.Year2021[i] + 
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
  # Calculate the expected number of females on Aug 31 (Spet 1) ----
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
  ###########################################################X
  # Recruitment: Poults per brood (PPB) model ----
  ###########################################################X
  #----------------------------------------------------------#
  # PPB: Likelihood ----
  #----------------------------------------------------------#
  ###########################################################X
  
  for (k in 1:ph.N) {
    
    #----------------------------------------------------------#
    # Note: For a gamma distribution, the variance is the shape times squared 
    # scale, which is equivalent to the mean (shape times scale) times the scale. 
    # So the variance scales with the mean (ph.mu).
    #----------------------------------------------------------#
    
    PHratio[k] ~ dgamma(shape = ph.alpha[k], scale = ph.theta[k])
    
    # Calculate the shape and scale parameter
    ph.alpha[k] <- ph.mu[k] * ph.disp
    ph.theta[k] <- 1/ph.disp
    
    # Calculate the mean for each observation
    #----------------------------------------------------------#
    # Note: intercept is 2019
    #----------------------------------------------------------#
    log(ph.mu[k]) <- ph.beta1 + ph.beta2 * ph.Year2020[k] + 
      ph.beta3 * ph.Year2021[k] + ph.beta4 * ph.Year2022[k] + 
      ph.beta5 * ph.Year2023[k] + ph.beta6 * ph.doy.scale[k] + 
      ph.beta7 * ph.doy.2[k] + ph.u[ph.wmu[k]]  
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
  ph.disp ~ dunif(0, 10) 
  
  # Random effect of WMU
  ph.sigma.u ~ dunif(0, 1)
  for (j in 1:ph.J) { 
    ph.u[j] ~ dnorm(0, sd = ph.sigma.u)
  }
  
  #----------------------------------------------------------#
  # Calculate the expected ratio of poults on Aug 31 (Spet 1) ----
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
  # Known-fate model ----
  ###########################################################X
  ###########################################################X
  # Known-fate Model: Priors on Intercept, Age, and WMU ----
  ###########################################################X
  # Intercept
  telem.beta.int[1] ~ dnorm(0, sd = 0.5)
  
  # Age effect
  telem.beta.age[1] ~ dnorm(0, sd = 0.5)
  
  ###########################################################X
  # WMU random effect  ----
  ###########################################################X
  telem.sigma ~ dunif(0, 5)
  for (w in 1:female.telem.wmu) {
    telem.beta.wmu[w] ~ dnorm(0, sd = telem.sigma)
  }
  
  ###########################################################X
  # Month random effect  ----
  ###########################################################X
  telem.month.sigma ~ dunif(0, 5)
  for (m in 1:12) {
    telem.beta.month[m] ~ dnorm(0, sd = telem.month.sigma)
  }
  
  ###########################################################X
  # Known Fate Model: Likelihood ----
  ###########################################################X
  
  # Loop over all individuals in the dataset
  for (i in 1:telem.nind) {  
    # Loop over the range of years in which each individual was tracked
    for (y in telem.year.start[i]:telem.year.end[i]) {  
      
      # Loop over the months between the first and last capture for each individual
      for (m in telem.first[i]:telem.last[i]) {  
        #-------------------------------------------------------------------------X
        # Note:
        # Derived survival rate on the complementary log-log (cloglog) scale
        # s.kf[i, y, m]: Survival for individual i in year y and month m
        # telem.beta.int[1]: Intercept term: adults
        # telem.beta.age[1]: Effect of age (juvenile) on survival
        # telem.juvenile[i, y, m]: Indicator for juvenile status (1 if juvenile, 0 if adult)
        # inprod(telem.beta.wmu[1:4], telem.wmu[i, 1:4]): Captures spatial random effect for WMU
        # telem.beta.month[m]: Monthly survival effect to capture seasonal variation
        #-------------------------------------------------------------------------X
        
        s.kf[i, y, m] <-  telem.beta.int[1] +   # Adult
          telem.beta.age[1] * telem.juvenile[i, y, m] +  # Juvenile-specific effect
          inprod(telem.beta.wmu[1:female.telem.wmu], telem.wmu[i, 1:female.telem.wmu]) +  # WMU random effect
          telem.beta.month[m]   # Month random effect
        
        #-------------------------------------------------------------------------X
        # Likelihood for survival status (0 = dead, 1 = alive) using a Bernoulli distribution
        # The probability is the inverse complementary log-log of the survival 
        # icloglog() is appropriate when values are skewed low (towards 0) or high (towards 1)
        #-------------------------------------------------------------------------X
        
        status[i, y, m] ~ dbern(prob = icloglog(s.kf[i, y, m]))
      }
    }
  }
  
  ###########################################################X
  # Derived Quantities: Survival over month per age and WMU ----
  ###########################################################X
  
  # Initialize storage for monthly survival rates (Adult and Juvenile)
  for (m in 1:12) {
    # Adults, baseline WMU
    storage[1, 1, m] <- icloglog(telem.beta.int[1] + telem.beta.wmu[1] +  telem.beta.month[m])   
    # Juveniles, baseline WMU
    storage[1, 2, m] <- icloglog(telem.beta.int[1] + telem.beta.age[1] + telem.beta.wmu[1] +  telem.beta.month[m]) 
  }
  
  # For WMUs beyond the first
  for (j in 2:(female.telem.wmu)) {
    for (m in 1:12) {
      # Adults
      storage[j, 1, m] <- icloglog(telem.beta.int[1] + telem.beta.wmu[j] + telem.beta.month[m])
      # Juveniles
      storage[j, 2, m] <- icloglog(telem.beta.int[1] + telem.beta.age[1] + telem.beta.wmu[j] + telem.beta.month[m])
    }
  }
  #---------------------------------------------------------------X
  # Calculate survival for Adults, accounting for nesting season 
  #---------------------------------------------------------------X
  for (j in 1:female.telem.wmu) {
    # Jan to December survival
    avg.ad.s.kf[j] <- prod(storage[j, 1, 1:12])  
  }
  #---------------------------------------------------------------X
  # Calculate survival for Juvenile FEMALES transitioning to Adults in June
  #---------------------------------------------------------------X
  for (j in 1:female.telem.wmu) {
    # November and December survival using Jan as a proxy
    juvenile_part1[j] <- storage[j, 2, 1] * storage[j, 2, 1]     
    # January to May
    juvenile_part2[j] <- prod(storage[j, 2, 1:5])   
    # June to December (aging up to adult)
    adult_part[j] <- prod(storage[j, 1, 6:10])     
    
    # Cumulative survival for each entry month
    avg.juv.s.kf[j] <- (juvenile_part1[j] * juvenile_part2[j])* adult_part[j] 
  }
  
  #---------------------------------------------------------------X
  # Calculate survival for Juvenile MALES surviving until Spring harvest 
  #---------------------------------------------------------------X
  for (j in 1:male.n.wmu) { 
    # November and December survival using Jan as a proxy
    juvenile_male_1[j] <- storage[j, 2, 1] * storage[j, 2, 1]     
    # January to April
    juvenile_male_2[j] <- prod(storage[j, 2, 1:4])   

    # Cumulative survival until Spring harvest 
    juv.male.adj[j] <- (juvenile_male_1[j] * juvenile_male_2[j])
  }
  
  ###########################################################X
  #----------------------------------------------------------#
  # Dead Recovery model for each sex ----
  ###########################################################X
  # DRM: Female model ----
  ###########################################################X
  #----------------------------------------------------------#
  # DRM: WMU harvest rates
  #----------------------------------------------------------#
  # Note: female.n.occasions-1 is done to match the notation for
  # my male DRM. See note in the male DRM likelihood for more info.
  #----------------------------------------------------------#
  for (t in 1:(female.n.occasions-1)) {
    for (u in 1:female.n.wmu) {
      female.h.ad.wmu[t, u] ~ dbeta(shape1 = 2, shape2 = 50)
      female.h.juv.wmu[t,u] ~ dbeta(shape1 = 2, shape2 = 50)
    } #u
  } #t
  
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
  male.seber.recov[2] ~ dunif(0, 1)   # Prior for proportion of mortality due to hunting - adults
  
  #----------------------------------------------------------#
  
  # Prior on juvenile effect on survival
  male.juvenile.effect ~ dnorm(0, sd = 0.5)        
  
  #----------------------------------------------------------#
  
  # Prior on time effect     
  for (t in 1:(male.n.occasions-1)) {   
    male.time.effect[t] ~ dnorm(0, sd = 0.5)             
  }
  
  #----------------------------------------------------------#
  
  # Random effect of wmu
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
        inprod(male.time.effect[1:4], male.time.param[t, 1:4]) + 
        male.wmu.effect[male.wmu[i]] 
      
      #----------------------------------------------------------#
      # Note: 
      # the four terms in r[i,t] represent: juvenile non-reward bands, 
      # juvenile reward bands, adult non-reward bands, adult reward band
      # the indicators I[] (juvenile=1) and II[] (reward = 0) 
      # simply retain/cancel the appropriate terms
      #----------------------------------------------------------#
      # Reporting rate w. Seber parameterization: ((1-s)r)
      male.r[i,t] <- male.seber.recov[1] * male.I[i, t] * male.rrate.j * male.II[i] + 
        male.seber.recov[1] * male.I[i, t] * (1 - male.II[i]) + male.seber.recov[2] * (1 - male.I[i, t]) * male.rrate.a * male.II[i] + 
        male.seber.recov[2] * (1 - male.I[i, t]) * (1 - male.II[i])
   
    } #t
  } #i
  
  ###########################################################X
  # DRM: Derived Estimates
  ###########################################################X
  #----------------------------------------------------------#
  # DRM: WMU survival rates
  #----------------------------------------------------------#
  # Survival from 2020-2023
  for (t in 1:(male.n.occasions-1)) {
    for (u in 1:male.n.wmu) {
      # Juvenile
      logit(male.s.juv.wmu[t,u]) <- male.juvenile.effect + 
        inprod(male.time.effect[1:4], male.time.param[t, 1:4]) + male.wmu.effect[u] 
      
      # Adult
      logit(male.s.ad.wmu[t,u]) <- inprod(male.time.effect[1:4], male.time.param[t, 1:4]) + male.wmu.effect[u]
    } #u
  } #t
  
  #----------------------------------------------------------#
  # DRM: WMU harvest rates
  #----------------------------------------------------------#
  # Harvest rates from 2020-2023
  for (t in 1:(male.n.occasions-1)) {
    for (u in 1:male.n.wmu) {
      # Juvenile
      male.h.juv.wmu[t,u] <- (1-male.s.juv.wmu[t,u])*male.seber.recov[1] 
      
      # Adult
      male.h.ad.wmu[t,u] <- (1-male.s.ad.wmu[t,u])*male.seber.recov[2]
    } #u
  } #t
  
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
    for (t in (male.f[i]+1):male.n.occasions){ 
      
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
  ##--------------------------------------------------------------------#X
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
   
    # adjust recruitment for males 
     male.recruitment[1, u] <- recruitment[1, u] * juv.male.adj[u]
    # male.N.juv[1, u] ~ dpois(male.recruitment[1, u]) # I dont think this can work because we are using the previous years recruitment for males
  } # u

  # Loop over subsequent time occasions (2021-2023)
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
      #----------------------------------------------------------#
      # Note: Poisson is to keep numbers whole for H/hr
      #----------------------------------------------------------#
      male.N.ad[t, u]  <- (male.N.ad.Survived[t, u] + male.N.juv.Survived[t, u])

      # New juveniles who entered the population
      #----------------------------------------------------------#
      # Note: For males, these juveniles are from Nov (t-1) entering in May (t)
      #----------------------------------------------------------#
      # adjust recruitment for males 
      male.recruitment[t, u] <- recruitment[t, u] * juv.male.adj[u]
      male.N.juv[t, u] ~ dpois(male.recruitment[t-1, u])

      # Lincoln-Petersen: Adult males in spring
      harvest.ad.spring[t, u] ~ dbin(prob = male.h.ad.wmu[t, u],
                                     size = male.N.ad[t, u])

      # Lincoln-Petersen: Juvenile males in spring
      harvest.juv.spring[t, u] ~ dbin(prob = male.h.juv.wmu[t, u],
                                      size = male.N.juv[t, u])

    } # t
  } # u

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

  } # u

  # Loop over subsequent time occasions (2021-2023)
  for (t in 2:Nyears) {
    for (u in 1:female.n.wmu) {
      # start: 2021
      # Recruitment for the time occasion (2021 - 2023)
      recruitment[t, u] <- ((female.N.ad[t, u] * aug31.hwb[t, u]) * aug31.ppb[t, u])/2

      # Number of surviving adult females from time t-1 to t (N1S1)
      female.N.ad.Survived[t, u] ~ dbin(size = female.N.ad[t-1, u],
                                         prob = avg.ad.s.kf[u])


      # Number of surviving adult females from time t-1 to t (N1S1)
      female.N.juv.Survived[t, u] ~ dbin(size = female.N.juv[t-1, u],
                                          prob = avg.juv.s.kf[u])

      # Adults who survived, and juveniles who have aged up to adult and survived
      #----------------------------------------------------------#
      # Note: Poisson is to keep numbers whole for H/hr
      #----------------------------------------------------------#
      female.N.ad[t, u] <- (female.N.ad.Survived[t, u] + female.N.juv.Survived[t, u])

      # New juveniles who entered the population
      # For females, these juveniles are from Nov (t) entering in Nov (t), so no adj.
      # is needed
      female.N.juv[t, u] ~ dpois(recruitment[t, u])


      # Lincoln-Petersen: Adult females in fall
      harvest.ad.fall[t, u] ~ dbin(prob = female.h.ad.wmu[t, u],
                                   size = female.N.ad[t, u])

      # Lincoln-Petersen: Juvenile females in fall
      harvest.juv.fall[t, u] ~ dbin(prob = female.h.juv.wmu[t, u],
                                    size = female.N.juv[t, u])


    } # t
  } # u
})
