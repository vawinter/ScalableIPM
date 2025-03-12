  ###########################################################X
  # Recruitment: Hen with brood (HWB) model ----
  ###########################################################X
  #----------------------------------------------------------#
  # HWB: Likelihood ----
  #----------------------------------------------------------#
  ###########################################################X
  
  for (i in 1:hwb.N) {
    # Hens with brood
    HWB[i] <- rbinom(1, 1, prob = hwb.p[i])
    logit(hwb.p[i]) <- hwb.beta1  + hwb.beta2 * hwb.Year2020[i] + hwb.beta3 * hwb.Year2021[i] + 
      hwb.beta4 * hwb.Year2022[i] + hwb.beta5 * hwb.Year2023[i] + hwb.beta6 * hwb.doy.scale[i] + 
      hwb.beta7 * hwb.doy.2[i] + hwb.u[hwb.wmu[i]]
  } # N
  
  #----------------------------------------------------------#
  # HWB: Priors ----
  #----------------------------------------------------------#
  
  # Intercept
  hwb.beta1 <- rnorm(1, 0, sd = 1)
  # year 2019 effect
  hwb.beta2 <- rnorm(1, 0, sd = 1)
  # year 2020 effect
  hwb.beta3 <- rnorm(1, 0, sd = 1)
  # year 2021 effect
  hwb.beta4 <- rnorm(1, 0, sd = 1)
  # year 2022 effect
  hwb.beta5 <- rnorm(1, 0, sd = 1)
  # DOY: make conformable for matrix multiplication
  # day of year
  hwb.beta6 <- rnorm(1, 0, sd = 1)
  # doy^2
  hwb.beta7 <- rnorm(1, 0, sd = 1) 
  # Random effect of wmu
  hwb.sigma <- runif(1, 0, 10)
  for (j in 1:hwb.J) {
    hwb.u[j] <- rnorm(1, 0, sd = hwb.sigma)
  } # J
  
  #----------------------------------------------------------#
  # Calculate the expected number of females on Aug 31 (Spet 1) ----
  #----------------------------------------------------------#
  
  # Loop over true occasions
  for (t in 1:female.n.occasions) {
    for (u in 1:female.n.wmu+1) {
      # Derived estimates for the number of hens on August 31 per WMU
      aug31.hwb[t,u] <- expit(
        hwb.beta1 + hwb.beta2 * (t == 2) +  hwb.beta3 * (t == 3) + 
          hwb.beta4 * (t == 4) +  hwb.beta5 * (t == 5) + 
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
  
  for (i in 1:ph.N) {
    
    #----------------------------------------------------------#
    # Note: For a gamma distribution, the variance is the shape times squared 
    # scale, which is equivalent to the mean (shape times scale) times the scale. 
    # So the variance scales with the mean (ph.mu).
    #----------------------------------------------------------#
    
    PHratio[i] <- rgamma(1, shape = ph.alpha[i], scale = ph.theta[i])
    
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
  ph.beta1 <- rnorm(1, 0, sd = 1)
  # year 2019 effect
  ph.beta2 <- rnorm(1, 0, sd = 1)
  # year 2020 effect
  ph.beta3 <- rnorm(1, 0, sd = 1)
  # year 2021 effect
  ph.beta4 <- rnorm(1, 0, sd = 1)
  # year 2022 effect
  ph.beta5 <- rnorm(1, 0, sd = 1)
  # DOY: make conformable for matrix multiplication
  # day of year
  ph.beta6 <- rnorm(1, 0, sd = 1)
  # doy^2
  ph.beta7 <- rnorm(1, 0, sd = 1) 
  
  # Scaling parameter for Gamma distribution
  ph.disp <- runif(1, 0, 10) 
  
  # Random effect of wmu
  ph.sigma.u  <- runif(1, 0, 1)
  for (j in 1:ph.J) { # Note 12/4: change to match WMU with telemetry
    ph.u[j] <- rnorm(1, 0, sd = ph.sigma.u)
  }
  
  #----------------------------------------------------------#
  # Calculate the expected ratio of poults on Aug 31 (Spet 1) ----
  #----------------------------------------------------------#
  
  # Loop over true occasions
  for (t in 1:female.n.occasions) {
    # Loop over management units
    for (u in 1:female.n.wmu) {
      # Derived estimates for the number of ppb on August 31 per WMU and year
      aug31.ppb[t,u] <- exp(
        ph.beta1 + ph.beta2 * (t == 2) + ph.beta3 * (t == 3) + 
          ph.beta4 * (t == 4) + ph.beta5 * (t == 5) + ph.beta6 * ppb.aug31 + 
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
  telem.beta.int[1] <- rnorm(1, 0, sd = 0.5)
  
  # Age effect
  telem.beta.age[1] <- rnorm(1,0, sd = 0.5)
  
  # WMU random effect
  telem.sigma <- runif(1, 0, 2)
  for (w in 1:3) {
    telem.beta.wmu[w] <- rnorm(1, 0, sd = telem.sigma)
  }
  
  ###########################################################X
  # Month Effects  ----
  ###########################################################X
  telem.month.sigma <- runif(1, 0, 2)
  for (m in 1:12) {
    telem.beta.month[m] <- rnorm(1, 0, sd = telem.month.sigma)
  }
  
  ###########################################################X
  # Year Effects ----
  ###########################################################X
  # telem.sigma.year ~ dunif(0, 2)
  # for (y in 1:telem.nyears) {
  #   telem.beta.year[y] ~ dnorm(0, sd = telem.sigma.year)
  # }
  # 
  
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
        # note:
        # Derived survival rate on the complementary log-log (cloglog) scale
        # s.kf[i, y, m]: Survival for individual i in year y and month m
        # telem.beta.int[1]: Intercept term: adult in WMU 4D
        # telem.beta.age[1]: Effect of age (juvenile) on survival
        # telem.juvenile[i, y, m]: Indicator for juvenile status (1 if juvenile, 0 if adult)
        # inprod(telem.beta.wmu[1:3], telem.wmu[i, 1:3]): Captures spatial random effect for WMU
        # telem.beta.month[m]: Monthly survival effect to capture seasonal variation
        #-------------------------------------------------------------------------X
        
        s.kf[i, y, m] <-  telem.beta.int[1] +   # Adult and WMU 4D (intercept)
          telem.beta.age[1] * telem.juvenile[i, y, m] +  # Juvenile-specific effect
          inprod(telem.beta.wmu[1:3], telem.wmu[i, 1:3]) +  # WMU random effect
          telem.beta.month[m]  # Fixed effect for the month to capture seasonal variability
        
        #-------------------------------------------------------------------------X
        # Likelihood for survival status (0 = dead, 1 = alive) using a Bernoulli distribution
        # The probability is the inverse complementary log-log of the survival 
        # icloglog() is appropriate when values are skewed low (towards 0) or high (towards 1)
        #-------------------------------------------------------------------------X
        
        status[i, y, m] <- rnorm(1, 1, prob = icloglog(s.kf[i, y, m]))
      }
    }
  }
  
  ###########################################################X
  # Derived Quantities: Survival over m per age and WMU ----
  ###########################################################X
  
  # Initialize storage for monthly survival rates (Adult and Juvenile)
  for (m in 1:12) {
    storage[1, 1, m] <- icloglog(telem.beta.int[1] + telem.beta.month[m])  # Adult in WMU 1
    storage[1, 2, m] <- icloglog(telem.beta.int[1] + telem.beta.age[1] + telem.beta.month[m])  # Juvenile in WMU 1
  }
  
  # Loop over WMUs beyond the first
  for (j in 1:(female.telem.wmu - 1)) {
    for (m in 1:12) {
      storage[j + 1, 1, m] <- icloglog(telem.beta.int[1] + telem.beta.wmu[j] + telem.beta.month[m])  # Adult
      storage[j + 1, 2, m] <- icloglog(telem.beta.int[1] + telem.beta.age[1] + telem.beta.wmu[j] + telem.beta.month[m])  # Juvenile
    }
  }
  
  # Calculate cumulative survival for Adults (no age change needed)
  for (j in 1:female.telem.wmu) {
    adult_part1[j] <- prod(storage[j, 1, 1:5])  # Adult survival from Jan to May (Months 1 to 5)
    adult_part2[j] <- prod(storage[j, 1, 6:11])    # Adult survival from June to Dec (Months 6 to 11)
    avg.ad.s.kf[j] <- adult_part1[j] * adult_part2[j]  # Adult: product of all months
  }
  
  # Calculate cumulative survival for Juveniles transitioning to Adults in June
  for (j in 1:female.telem.wmu) {
    juvenile_part[j] <- prod(storage[j, 2, 1:5])  # Juvenile survival from Jan to May (Months 1 to 5)
    adult_part[j] <- prod(storage[j, 1, 6:11])    # Adult survival from June to Dec (Months 6 to 11)
    avg.juv.s.kf[j] <- juvenile_part[j] * adult_part[j]  # Total annual survival for a juvenile aging into an adult
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
  
  for (t in 1:(female.n.occasions-1)) {
    for (u in 1:female.n.wmu) {
      # Use a Beta prior on harvest rate 
      female.h.ad.wmu[t,u] <- rbeta(1, shape1 = 2, shape2 = 50)
      female.h.juv.wmu[t,u] <- rbeta(1, shape1 = 2, shape2 = 50)
    }#g
  }#t
  
  ###########################################################X
  # DRM: Male model ----
  ###########################################################X
  # Dead Recovery model: Priors on age and wmu 
  #----------------------------------------------------------#
  ###########################################################X
  
  # Reporting rates
  ### Informative priors for hunter reporting rates based on Diefenbach et al. (2012)
  male.rrate.j <- rnorm(1, 0.71, sd = 0.072) #  non-reward bands for juveniles
  male.rrate.a <- rnorm(1, 0.87, sd = 0.039) #  non-reward bands for adults
  
  #----------------------------------------------------------#
  # Proportion mortality due to hunting (Seber parameterization: (1-s)r)
  male.seber.recov[1] <- runif(1, 0, 1)   # Prior for proportion mortality due to hunting - juveniles
  male.seber.recov[2] <- runif(1, 0, 1)   # Prior for proportion of mortality due to hunting - adults
  
  #----------------------------------------------------------#
  
  # Prior on juvenile effect on survival
  male.juvenile.effect <- rnorm(1, 0, sd = 0.5)        
  
  #----------------------------------------------------------#
  
  # Prior on time effect     
  for (t in 1:(male.n.occasions-1)) {   
    male.time.effect[t] <- rnorm(1, 0, sd = 0.5)             
  }
  
  #----------------------------------------------------------#
  
  # Random effect of wmu
  male.sigma ~ dunif(0,10)
  for (i in 1:10) {
    male.wmu.effect[i] <- rnorm(1, 0, sd = male.sigma)
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
    # Juvenile
    logit(male.mean.s.juv[t]) <- male.time.effect[t] + male.juvenile.effect  
    # Adult
    logit(male.mean.s.ad[t]) <- male.time.effect[t]  
    #----------------------------------------------------------#
    # DRM: Derived harvest rate estimates
    #----------------------------------------------------------#
    # Juvenile
    male.mean.harv.jv[t] <- (1-male.mean.s.juv[t])*male.seber.recov[1]
    # Adult
    male.mean.harv.ad[t] <- (1-male.mean.s.ad[t])*male.seber.recov[2]
  } # t
  
  #----------------------------------------------------------#
  # DRM: WMU survival rates
  #----------------------------------------------------------#
  
  for (t in 1:(male.n.occasions-1)) {
    for (u in 1:male.n.wmu) {
      # Juvenile
      logit(male.s.juv.wmu[t,u]) <- male.juvenile.effect + 
        inprod(male.time.effect[1:4], male.time.param[t, 1:4]) + male.wmu.effect[u] 
      
      # Adult
      logit(male.s.ad.wmu[t,u]) <- inprod(male.time.effect[1:4], male.time.param[t, 1:4]) + male.wmu.effect[u]
    }#u
  }#t
  
  #----------------------------------------------------------#
  # DRM: WMU harvest rates
  #----------------------------------------------------------#
  
  for (t in 1:(male.n.occasions-1)) {
    for (u in 1:male.n.wmu) {
      # Juvenile
      male.h.juv.wmu[t,u] <- (1-male.s.juv.wmu[t,u])*male.seber.recov[1] 
      
      # Adult
      male.h.ad.wmu[t,u] <- (1-male.s.ad.wmu[t,u])*male.seber.recov[2]
    }#u
  }#t
  
  #----------------------------------------------------------#
  # DRM: Likelihood (Males) ----
  #----------------------------------------------------------#
  
  for (i in 1:male.nind){
    # Define latent state at first capture
    male.z[i,male.f[i]] <- 1 
    for (t in (male.f[i]+1):male.n.occasions){ 
      
      #------------------------#
      
      # State process
      male.z[i,t] <- rnorm(1, 1, male.mu1[i,t])
      male.mu1[i,t] <- male.s[i,t-1] * male.z[i,t-1]
      #------------------------#
      
      # Observation process
      male.y[i,t] <- rnorm(1, 1, male.mu2[i,t])
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
  # Male Abundance
  #---------------------------------------------------------------------#X
  
  for (u in 1:male.n.wmu) {
    
    # Initial abundance for adult males in each WMU
    N.lambda.ad.male[u] <- (th.year1.male.ad[u]) / male.h.ad.wmu[1, u]
    male.N.ad[1, u] <- rpois(1, N.lambda.ad.male[u])
    
    # Initial abundance for juvenile males in each WMU
    N.lambda.juv.male[u] <- (th.year1.male.juv[u]) / male.h.juv.wmu[1, u]
    male.N.juv[1, u] <- rpois(1, N.lambda.juv.male[u])
  } # end u
  
  # Loop over time occasions (starting from the second occasion)
  for (t in 1:(male.n.occasions-1)) {
    for (u in 1:male.n.wmu) {
      
      
      # Lincoln-Petersen: Adult males in spring
      harvest.ad.spring[t, u] <- rbinom(1, prob = male.h.ad.wmu[t, u],
                                     size = male.N.ad[t, u])
      
      # Lincoln-Petersen: Juvenile males in spring
      harvest.juv.spring[t, u] <- rbinom(1, prob = male.h.juv.wmu[t, u],
                                      size = male.N.juv[t, u])
      
      # Number of surviving adult males from time t-1 to t
      male.N.ad.post.harv[t, u] <- rbinom(1, prob = male.s.ad.wmu[t, u],
                                       size = male.N.ad[t, u])
      
      # Number of surviving juvenile males from time t-1 to t
      male.N.juv.post.harv[t, u] <- rbinom(1, prob = male.s.juv.wmu[t, u],
                                        size = male.N.juv[t, u])
      
      # Process model for adult male population count
      delta.ad.male[t, u] <- (male.N.ad.post.harv[t, u] + male.N.juv.post.harv[t, u])
      
      # Recruitment for juvenile males
      delta.juv.male[t, u] <- (recruitment[t, u])
      
      # Update the number of adult males for the next time occasion
      male.N.ad[t+1, u] <- rpois(1, delta.ad.male[t, u])
      
      # Update the number of juvenile males for the next time occasion
      male.N.juv[t+1, u] <- rpois(1, delta.juv.male[t, u])
      
    } # end t
  } # end u
  
  #---------------------------------------------------------------------#X
  # Female Abundance
  #---------------------------------------------------------------------#X
  
  for (u in 1:female.n.wmu) {
    
    # Initial abundance for adult females in each WMU
    N.lambda.ad.female[u] <- (th.year1.female.ad[u]) / female.h.ad.wmu[1, u]
    female.N.ad[1, u] <- rpois(1, N.lambda.ad.female[u])
    
    # Initial abundance for juvenile females in each WMU
    N.lambda.juv.female[u] <- (th.year1.female.juv[u]) / female.h.juv.wmu[1, u]
    female.N.juv[1, u] <- rpois(1, N.lambda.juv.female[u])
    
  } # end u
  
  # Loop over time occasions (starting from the second occasion)
  for (t in 1:(female.n.occasions-1)) {
    for (u in 1:female.n.wmu) {
      
      
      # Lincoln-Petersen: Adult females in fall
      harvest.ad.fall[t, u] <- rbinom(1, prob = female.h.ad.wmu[t, u],
                                   size = female.N.ad[t, u])
      
      # Lincoln-Petersen: Juvenile females in fall
      harvest.juv.fall[t, u] <- rbinom(1, prob = female.h.juv.wmu[t, u],
                                    size = female.N.juv[t, u])
      
      
      # Recruitment for the next time occasion
      recruitment[t, u] <- ((female.N.ad[t, u] * aug31.hwb[t, u]) * aug31.ppb[t, u])/2
      
      # Number of surviving adult females from time t-1 to t
      female.N.ad.post.harv[t, u] <- rbinom(1, size = female.N.ad[t, u],
                                         avg.ad.s.kf[u])
      
      
      # Number of surviving adult females from time t-1 to t
      female.N.juv.post.harv[t, u] <- rbinom(1, size = female.N.juv[t, u],
                                          prob = avg.juv.s.kf[u])
      
      # Process model for adult female population count
      delta.ad.female[t, u] <- female.N.ad.post.harv[t, u] + (female.N.juv.post.harv[t, u])
      
      # Process model for juveniles count
      delta.juv.female[t, u] <- (recruitment[t, u])
      
      # Update the number of adult females for the next time occasion
      female.N.ad[t+1, u]  <- rpois(1, delta.ad.female[t, u])
      
      # Update the number of juveniles for the next time occasion
      female.N.juv[t+1, u] <- rpois(1, delta.juv.female[t, u])
      
    } # end t
  } # end u
  