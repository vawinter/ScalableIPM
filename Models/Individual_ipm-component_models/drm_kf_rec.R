drm_kf_rec <- nimbleCode({
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
    logit(hwb.p[i]) <- hwb.beta1 * hwb.Year2019[i] + hwb.beta2 * hwb.Year2020[i] + hwb.beta3 * hwb.Year2021[i] + 
      hwb.beta4 * hwb.Year2022[i] + hwb.beta5 * hwb.Year2023[i] + hwb.beta6 * hwb.doy.scale[i] + 
      hwb.beta7 * hwb.doy.2[i] + hwb.u[hwb.wmu[i]]
  } # N
  
  #----------------------------------------------------------#
  # HWB: Priors ----
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
  ###########################################################X
  # Recruitment: Poults per brood (PPB) model ----
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
  # PPB: Priors  ----
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
  # DOY: make conformable for matrix multiplication
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
  ###########################################################X
  # Known-fate model ----
  ###########################################################X
  #----------------------------------------------------------#
  # Known-fate model: Priors on age and wmu ----
  #----------------------------------------------------------#
  ###########################################################X
  # Intercept
  telem.beta.int[1] ~ dnorm(0, sd = 0.5)
  
  # Age
  telem.beta.age[1] ~ dnorm(0, sd = 0.5) 
  
  
  # WMU
  for(w in 1:3){
    telem.beta.wmu[w] ~ dnorm(0, sd = 1)
  }
  
  #----------------------------------------------------------#
  # # Known-fate model: Likelihood ----
  #----------------------------------------------------------#
  for (i in 1:telem.nind) {
    # first to last encounter of each individual
    for (t in (telem.first[i]):telem.last[i]) {
      # Derived individual survival rate with inverse cloglog
      s.kf[i, t] <- telem.beta.int[1] + telem.beta.age[1] * telem.juvenile[i, t] +
        inprod(telem.beta.wmu[1:3], telem.wmu[i, 1:3])
      
      # Likelihood for individual capture histories
      status[i, t] ~ dbern(prob = icloglog(s.kf[i, t]))
      
    }
  }
  
  #----------------------------------------------------------#
  # Cumulative survival per age, wmu, and 'cohort' ----
  #----------------------------------------------------------#
  # Initialize storage vector
  # intercept = adult, wmu 1
  storage[1,1] <- icloglog(telem.beta.int[1]) # wmu 1 + adult
  storage[1,2] <- icloglog(telem.beta.int[1] + telem.beta.age[1]) # wmu 1 + juvenile
  
  # loop over wmus
  for (j in 1:(female.telem.wmu-1)) {
    # Adult in wmu + 1
    storage[j+1,1] <- icloglog(telem.beta.int[1] + telem.beta.wmu[j])
    # Juvenile in wmu + 1
    storage[j+1,2] <- icloglog(telem.beta.int[1] + telem.beta.age[1] + telem.beta.wmu[j])
  }
  
  # Calculate cumulative survival: adult
  female.s.kf[1:female.telem.wmu,1,1] <- storage[1:female.telem.wmu,1]^10 # Jan (Dec + 1)
  female.s.kf[1:female.telem.wmu,1,2] <- storage[1:female.telem.wmu,1]^9 # Feb (Jan + 1)
  female.s.kf[1:female.telem.wmu,1,3] <- storage[1:female.telem.wmu,1]^8 # March
  female.s.kf[1:female.telem.wmu,1,4] <- storage[1:female.telem.wmu,1]^7 # April
  
  # Calculate cumulative survival: juvenile (until June, then age up to adult)
  female.s.kf[1:female.telem.wmu,2,1] <- storage[1:female.telem.wmu,2]^5 * storage[1:female.telem.wmu,1]^5
  female.s.kf[1:female.telem.wmu,2,2] <- storage[1:female.telem.wmu,2]^4 * storage[1:female.telem.wmu,1]^5
  female.s.kf[1:female.telem.wmu,2,3] <- storage[1:female.telem.wmu,2]^3 * storage[1:female.telem.wmu,1]^5
  female.s.kf[1:female.telem.wmu,2,4] <- storage[1:female.telem.wmu,2]^2 * storage[1:female.telem.wmu,1]^5
  
  #----------------------------------------------------------#
  # Cumulative survival per age, wmu not above, and 'cohort' ----
  #----------------------------------------------------------#
  # Create prior on cumulative survival for remaining wmus
  # Adults
  female.s.kf.est[5, 1, 1] ~ dunif(0.7, 0.8) # alternative to cloglog, restricting values
  female.s.kf.est[5, 1, 2] ~ dunif(0.7, 0.8) 
  female.s.kf.est[5, 1, 3] ~ dunif(0.7, 0.8) 
  female.s.kf.est[5, 1, 4] ~ dunif(0.7, 0.8) 
  
  female.s.kf.est[6, 1, 1] ~ dunif(0.7, 0.8) 
  female.s.kf.est[6, 1, 2] ~ dunif(0.7, 0.8) 
  female.s.kf.est[6, 1, 3] ~ dunif(0.7, 0.8) 
  female.s.kf.est[6, 1, 4] ~ dunif(0.7, 0.8) 
  
  # Juvenile (until June, then age up to adult)
  female.s.kf.est[5, 2, 1] ~ dunif(0.7, 0.8) # alternative to cloglog, restricting values
  female.s.kf.est[5, 2, 2] ~ dunif(0.7, 0.8) 
  female.s.kf.est[5, 2, 3] ~ dunif(0.7, 0.8) 
  female.s.kf.est[5, 2, 4] ~ dunif(0.7, 0.8) 
  
  female.s.kf.est[6, 2, 1] ~ dunif(0.7, 0.8) 
  female.s.kf.est[6, 2, 2] ~ dunif(0.7, 0.8) 
  female.s.kf.est[6, 2, 3] ~ dunif(0.7, 0.8) 
  female.s.kf.est[6, 2, 4] ~ dunif(0.7, 0.8) 
  
  # cloglog option
  # icloglog(female.s.kf.est[i, 1, j]) ~ dnorm(0, sd = 0.5) # option 1
  
  #----------------------------------------------------------#
  # Creating array of adjusted survival consisting of derived (female.s.kf)
  # and estimated (female.s.kf.est) survival to use in DRM likelihood 
  # for females. 
  #
  # 'adj.female.s.kf' will have derived cumulative survival for wmu 1:4 and then
  # estimated cumulative survival for wmu 5:610, where we do not have this info.
  #----------------------------------------------------------#
  # Fill storage vector with KNOWN cumulative survival in telemetered wmus
  for (i in 1:female.telem.wmu) { # female.telem.wmu = 4 (wmus w. telemetered birds)
    for (j in 1:4) { # 1:4 = cohort group
      # Adult
      adj.female.s.kf[i, 1, j] <- female.s.kf[i, 1, j] # derived cumulative survival (adult)
      # Juvenile
      adj.female.s.kf[i, 2, j] <- female.s.kf[i, 2, j] # derived cumulative survival (juvenile)
    }
  }
  #----------------------------------------------------------#
  # Fill storage vector with prior on cumulative survival for remaining wmus
  for (i in 5:female.n.wmu) { # female.n.wmu = 10 (all wmus, looping to exclude first 4)
    for (j in 1:4) { # 1:4 = cohort group
      # Adult
      adj.female.s.kf[i, 1, j] <- female.s.kf.est[i, 1, j] # estimated cumulative survival (adult)
      # Juvenile
      adj.female.s.kf[i, 2, j] <- female.s.kf.est[i, 2, j] # estimated cumulative survival (juvenile)
    }
  }
  
  ###########################################################X
  #----------------------------------------------------------#
  # Dead Recovery model for each sex ----
  ###########################################################X
  # DRM: Female model ----
  ###########################################################X
  #----------------------------------------------------------#
  # Dead Recovery model: Priors on age and wmu 
  #----------------------------------------------------------#
  # Reporting rates:
  # Note:
  ## We don't have an informative prior like we use for males for females, so 
  ## we want a vague prior on reporting rates
  female.rrate.j ~ dunif(0,1)
  female.rrate.a ~ dunif(0,1)
  #----------------------------------------------------------#
  
  # Proportion mortality due to hunting (Seber parameterization: (1-s)r)
  female.seber.recov[1] ~ dunif(0,1) # Prior for proportion mortality due to hunting - juveniles
  female.seber.recov[2] ~ dunif(0,1) # Prior for proportion of mortality due to hunting - adults
  
  #----------------------------------------------------------#
  
  # Prior on juvenile effect on survival
  female.juvenile.effect ~ dnorm(0, sd = 0.5)        
  #----------------------------------------------------------#
  # Prior on time effect     
  for (t in 1:(true.occasions)) {   
    female.time.effect[t] ~ dnorm(0, sd = 0.5)             
  }
  #----------------------------------------------------------#
  ## Random effect of wmu
  female.sigma ~ dunif(0,10)
  for (i in 1:10) {
    female.wmu.effect[i] ~ dnorm(0, sd = female.sigma)
  }
  
  #----------------------------------------------------------#
  # Survival and Recovery rate
  for (i in 1:female.nind) {
    for (t in female.f[i]:(true.occasions)) {
      ### juvenile.effect is age (juv) coefficient, time.effect is time coefficient 
      logit(female.s[i,t]) <- female.juvenile.effect*female.I[i,t] + 
        inprod(female.time.effect[1:4], female.time.param[t, 1:4]) + 
        female.wmu.effect[female.wmu[i]] 
      ###########################################################X
      # Note:
      ### the four terms in r[i,t] represent: juvenile non-reward bands, 
      ### juvenile reward bands, adult non-reward bands, adult reward band
      ### the indicators I[] (juvenile=1) and II[] (reward = 0) 
      ### simply retain/cancel the appropriate terms
      ###########################################################X
      # Intermediate variables for r[i,t] rates 
      female.juvenile.hunting[i, t] <- female.seber.recov[1] * female.I[i, t] * female.rrate.j * female.II[i] + 
        female.seber.recov[1] * female.I[i, t] * (1 - female.II[i])
      female.adult.hunting[i, t] <- female.seber.recov[2] * (1 - female.I[i, t]) * female.rrate.a * female.II[i] + 
        female.seber.recov[2] * (1 - female.I[i, t]) * (1 - female.II[i])
      #----------------------------------------------------------#
      # Recovery rate w. Seber parameterization ((1-s)r)
      female.r[i,t] <- female.juvenile.hunting[i, t] + female.adult.hunting[i, t]
      #----------------------------------------------------------#
    } #t
  } #i
  
  ###########################################################X
  # DRM: Derived Estimates
  ###########################################################X
  #----------------------------------------------------------#
  # DRM: Derived survival rate estimates
  #----------------------------------------------------------#
  for (t in 1:(true.occasions)) {   
    logit(female.mean.s.ad[t]) <- female.time.effect[t]  
    logit(female.mean.s.jv[t]) <- female.time.effect[t] + female.juvenile.effect  
    
    
    #----------------------------------------------------------#
    # DRM: Derived harvest rate estimates
    #----------------------------------------------------------#
    female.mean.harv.jv[t] <- (1-female.mean.s.jv[t])*female.seber.recov[1]
    female.mean.harv.ad[t] <- (1-female.mean.s.ad[t])*female.seber.recov[2]
  }
  
  #----------------------------------------------------------#
  # DRM: WMU survival rates
  #----------------------------------------------------------#
  for (t in 1:(true.occasions)) {
    for (u in 1:female.n.wmu) {
      logit(female.s.ad.wmu[t,u]) <-  female.time.effect[t] + female.wmu.effect[u]  
      logit(female.s.juv.wmu[t,u]) <- female.juvenile.effect + 
        # I don't like this hard coding, need to find nimble work-around
        inprod(female.time.effect[1:4], female.time.param[t, 1:4]) +  
        female.wmu.effect[u] 
    }#g
  }#t
  
  #----------------------------------------------------------#
  # DRM: WMU harvest rates
  #----------------------------------------------------------#
  for (t in 1:(true.occasions)) {
    for (g in 1:female.n.wmu) {
      female.h.juv.wmu[t,g] <- (1-female.s.juv.wmu[t,g])*female.seber.recov[1] 
      female.h.ad.wmu[t,g] <- (1-female.s.ad.wmu[t,g])*female.seber.recov[2]
    }#g
  }#t
  
  #----------------------------------------------------------#
  # DRM: Likelihood (Females) ----
  #----------------------------------------------------------#
  for (i in 1:female.nind){
    # Define latent state at first capture
    female.z[i,female.f[i]] <- 1
    for (t in (female.f[i]+1):female.n.occasions){ # true.occasions + 1 b.c they're dead
      #------------------------#
      # State process
      female.z[i,t] ~ dbern(female.mu1[i,t])
      
      # For the first year the female enters the study (female.f.kf), multiply 
      # DRM survival probability (female.s) by known-fate survival 
      # probability (adj.female.s.kf)
      
      # Note: female.age is a vector of 1, 2=Juvenile for indexing purposes based on
      # age at capture. The juvenile to adult age change June 1 is captured in the derived calculation
      # for 'adj.female.s.kf' above (see 'Cumulative survival per age, wmu, and 'cohort'')
      
      female.mu1[i,t] <- female.f.kf[i,t-1]*adj.female.s.kf[female.wmu.index[i],
                                                            female.age[i], tagging.cohort[i]]*female.s[i,t-1]+
        (1-female.f.kf[i,t-1])*female.s[i,t-1] * female.z[i,t-1]
      #------------------------#
      # Observation process
      female.y[i,t] ~ dbern(female.mu2[i,t])
      female.mu2[i,t] <- female.r[i,t-1] * (female.z[i,t-1] - female.z[i,t])
    } #t
  } #i
  
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
  male.seber.recov[2] ~ dunif(0,1)    # Prior for proportion of mortality due to hunting - adults
  
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
    for (t in male.f[i]:(true.occasions)) {
      ###########################################################X
      # Note:
      ### juvenile.effect is age (juv) coefficient, time.effect is time coefficient 
      ###########################################################X
      logit(male.s[i,t]) <- male.juvenile.effect*male.I[i,t] + 
        # Hate this hard coding - need to find nimble work-around
        inprod(male.time.effect[1:4], male.time.param[t, 1:4]) + 
        male.wmu.effect[male.wmu[i]] 
      #----------------------------------------------------------#
      ###########################################################X
      # Note: 
      ### the four terms in r[i,t] represent: juvenile non-reward bands, 
      ### juvenile reward bands, adult non-reward bands, adult reward band
      ### the indicators I[] (juvenile=1) and II[] (reward = 0) 
      ### simply retain/cancel the appropriate terms
      ###########################################################X
      # Intermediate variables for r[i,t] rates 
      male.juvenile.hunting[i, t] <- male.seber.recov[1] * male.I[i, t] * male.rrate.j * male.II[i] + 
        male.seber.recov[1] * male.I[i, t] * (1 - male.II[i])
      male.adult.hunting[i, t] <- male.seber.recov[2] * (1 - male.I[i, t]) * male.rrate.a * male.II[i] + 
        male.seber.recov[2] * (1 - male.I[i, t]) * (1 - male.II[i])
      #----------------------------------------------------------#
      # Recovery rate w. Seber parameterization: ((1-s)r)
      male.r[i,t] <- male.juvenile.hunting[i, t] + male.adult.hunting[i, t]
      #----------------------------------------------------------#
    } #t
  } #i
  
  ###########################################################X
  # DRM: Derived Estimates
  ###########################################################X
  #----------------------------------------------------------#
  # DRM: Derived survival rate estimates
  #----------------------------------------------------------#
  for (t in 1:(true.occasions)) {   
    logit(male.mean.s.ad[t]) <- male.time.effect[t]  
    logit(male.mean.s.jv[t]) <- male.time.effect[t] + male.juvenile.effect  
    
    #----------------------------------------------------------#
    # DRM: Derivedharvest rate estimates
    #----------------------------------------------------------#
    male.mean.harv.jv[t] <- (1-male.mean.s.jv[t])*male.seber.recov[1]
    male.mean.harv.ad[t] <- (1-male.mean.s.ad[t])*male.seber.recov[2]
  }
  
  #----------------------------------------------------------#
  # DRM: WMU survival rates
  #----------------------------------------------------------#
  for (t in 1:(true.occasions)) {
    for (u in 1:male.n.wmu) {
      logit(male.s.ad.wmu[t,u]) <-  male.time.effect[t] + male.wmu.effect[u]  
      logit(male.s.juv.wmu[t,u]) <- male.juvenile.effect + 
        inprod(male.time.effect[1:4], male.time.param[t, 1:4]) + male.wmu.effect[u] 
    }#g
  }#t
  
  #----------------------------------------------------------#
  # DRM: WMU harvest rates
  #----------------------------------------------------------#
  for (t in 1:(true.occasions)) {
    for (g in 1:male.n.wmu) {
      male.h.juv.wmu[t,g] <- (1-male.s.juv.wmu[t,g])*male.seber.recov[1] 
      male.h.ad.wmu[t,g] <- (1-male.s.ad.wmu[t,g])*male.seber.recov[2]
    }#g
  }#t
  
  #----------------------------------------------------------#
  # DRM: Likelihood (Males) ----
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
  
  
})
