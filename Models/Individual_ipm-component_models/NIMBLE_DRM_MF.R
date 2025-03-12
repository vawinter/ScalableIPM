drm <- nimbleCode({
  ###########################################################X
  #----------------------------------------------------------#
  # Dead Recovery model for each sex
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
  # DRM: Likelihood (Males)
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
  
  ###########################################################X
  # DRM: Female model ----
  ###########################################################X
  #----------------------------------------------------------#
  # Dead Recovery model: Priors on age and wmu 
  #----------------------------------------------------------#
  # Reporting rates:
  # Note:
  ## We don't have an informative prior like we used for males for females, so 
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
  # DRM: Likelihood (Females)
  #----------------------------------------------------------#
  for (i in 1:female.nind){
    # Define latent state at first capture
    female.z[i,female.f[i]] <- 1
    for (t in (female.f[i]+1):female.n.occasions){ # true.occasions + 1 b.c they're dead
      #------------------------#
      # State process
      female.z[i,t] ~ dbern(female.mu1[i,t])
      female.mu1[i,t] <- female.s[i,t-1] * female.z[i,t-1]
      #------------------------#
      # Observation process
      female.y[i,t] ~ dbern(female.mu2[i,t])
      female.mu2[i,t] <- female.r[i,t-1] * (female.z[i,t-1] - female.z[i,t])
    } #t
  } #i

})