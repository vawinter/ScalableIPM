drm <- nimbleCode({
  #----------------------------------------------------------#
  # FEMALE DEAD RECOVERY MODEL ----
  #----------------------------------------------------------#
  #----------------------------------------------------------#
  # Dead Recovery model: Priors on age and wmu
  #----------------------------------------------------------#
  # # Priors and constraints
  for (i in 1:nind.female) {
    for (t in f.female[i]:(n.occasions.female-1)) {
      ### alpha is age (juv) coefficient, beta is time coefficient with beta[1] the mean level, gamma is WMU random effect
      logit(s.female[i,t]) <- alpha.female*I.female[i,t] + inprod(beta.time.female[1:(n.occasions.female-1)], time.female[t, 1:(n.occasions.female-1)]) + 
        gamma.wmu.female[wmu.female[i]]
     
       ### the four terms in r[i,t] represent: juvenile non-reward bands, juvenile reward bands, adult non-reward bands, adult reward band
      ### the indicators I[] (juvenile=1) and II[] (reward = 0) simply retain/cancel the appropriate terms
      r.female[i,t] <- rr.female[1]*I.female[i,t]*rrate.jenny*II.female[i] + rr.female[1]*I.female[i,t]*(1-II.female[i]) + 
        rr.female[2]*(1-I.female[i,t])*rrate.hen*II.female[i] + rr.female[2]*(1-I.female[i,t])*(1-II.female[i])
    } #t
  } #i
  
  ## Uninformative priors on hunter reporting rates 
  rrate.jenny ~ dnorm(0, sd = 1) #  non-reward bands for juveniles 
  rrate.hen ~ dnorm(0, sd = 1) #  non-reward bands for adults
  
  # Proportion mortality due to hunting
  rr.female[1] ~ dunif(0, 1)  # juveniles # rr*(1-S) = harvest rate r in seber param 
  rr.female[2] ~ dunif(0, 1)  # adults
  
  
  ## time effect
  for(t in 1:(n.occasions.female-1)){
    beta.time.female[t] ~ dnorm(0, sd = 0.5)
  }
  
  # WMU effect - random effect does not look great for females
  ## Random effect of wmu
  random.effect.wmu.female ~ dunif(0, 10)
  for(u in 1:n.wmu){
    gamma.wmu.female[u] ~ dnorm(0, sd = random.effect.wmu.female)
  }
  
  # Juvenile effect on survival
  alpha.female ~ dnorm(0, sd = 0.5)             
  ###########################################################X
  # DRM: Derived Estimates
  ###########################################################X
  #----------------------------------------------------------#
  # DRM: Derived STATEWIDE survival rates
  #----------------------------------------------------------#
  for (t in 1:(n.occasions.female-1)) {
    logit(s.ad.female[t]) <- inprod(beta.time.female[1:(n.occasions.female-1)], time.female[t, 1:(n.occasions.female-1)])
    logit(s.jv.female[t]) <- alpha.female + inprod(beta.time.female[1:(n.occasions.female-1)], time.female[t, 1:(n.occasions.female-1)]) # J = means parameterization
  }#t
  #----------------------------------------------------------#
  # DRM: Derived MEAN STATEWIDE harvest rate estimates
  #----------------------------------------------------------#
  logit(mean.s.ad.female) <- beta.time.female[1]
  logit(mean.s.jv.female) <- beta.time.female[1] + alpha.female
  
  #----------------------------------------------------------#
  # DRM: Derived STATEWIDE harvest rate estimates
  #----------------------------------------------------------#
  for (t in 1:(n.occasions.female-1)) {
    h.juv.female[t] <- (1-s.jv.female[t])*rr.female[1]
    h.ad.female[t] <- (1-s.ad.female[t])*rr.female[2]
  }#t
  #----------------------------------------------------------#
  # DRM: Derived MEAN STATEWIDE harvest rate estimates
  #----------------------------------------------------------#
  mean.h.jv.female <- (1-mean.s.jv.female)*rr.female[1]
  mean.h.ad.female <- (1-mean.s.ad.female)*rr.female[2]
  
  #----------------------------------------------------------#
  # DRM: WMU survival rates
  #----------------------------------------------------------#
  for (t in 1:(n.occasions.female-1)) {
    for (u in 1:n.wmu) {
      logit(s.ad.wmu.female[t,u]) <- inprod(beta.time.female[1:(n.occasions.female-1)], time.female[t, 1:(n.occasions.female-1)]) + gamma.wmu.female[u]
      logit(s.juv.wmu.female[t,u]) <- alpha.female + inprod(beta.time.female[1:(n.occasions.female-1)], time.female[t, 1:(n.occasions.female-1)]) + gamma.wmu.female[u]
    }#g
  }#t
  
  #----------------------------------------------------------#
  # DRM: WMU harvest rates
  #----------------------------------------------------------#
  for (t in 1:(n.occasions.female-1)) {
    for (g in 1:n.wmu) {
      h.juv.wmu.female[t,g] <- (1-s.juv.wmu.female[t,g])*rr.female[1]
      h.ad.wmu.female[t,g] <- (1-s.ad.wmu.female[t,g])*rr.female[2]
    }#g
  }#t
  #----------------------------------------------------------#
  # DRM: Likelihood
  #----------------------------------------------------------#
  for (i in 1:nind.female){
    # Define latent state at first capture
    z.female[i,f.female[i]] <- 1
    for (t in (f.female[i]+1):n.occasions.female){
      # State process
      mu1.female[i,t] <- s.female[i,t-1] * z.female[i,t-1]
      z.female[i,t] ~ dbern(mu1.female[i,t])
      # Observation process
      mu2.female[i,t] <- r.female[i,t-1] * (z.female[i,t-1] - z.female[i,t])
      y.female[i,t] ~ dbern(mu2.female[i,t])
    } #t
  } #i
  
  #----------------------------------------------------------#
  # MALE DEAD RECOVERY MODEL ----
  #----------------------------------------------------------#
  #----------------------------------------------------------#
  # Dead Recovery model: Priors on age and wmu
  #----------------------------------------------------------#
  # # Priors and constraints
  for (i in 1:nind) {
    for (t in f.male[i]:(n.occasions.male-1)) {
      ### alpha is age (juv) coefficient, beta is time coefficient with beta[1] the mean level, gamma is WMU random effect
      logit(s.male[i,t]) <- alpha.male*I.male[i,t] + inprod(beta.time.male[1:(n.occasions.male-1)], time.male[t, 1:(n.occasions.male-1)]) +
        gamma.male[wmu.male[i]]
      ### the four terms in r[i,t] represent: juvenile non-reward bands, juvenile reward bands, adult non-reward bands, adult reward band
      ### the indicators I[] (juvenile=1) and II[] (reward = 0) simply retain/cancel the appropriate terms
      r.male[i,t] <- rr.male[1]*I.male[i,t]*rrate.jake*II.male[i] + rr.male[1]*I.male[i,t]*(1-II.male[i]) +
        rr.male[2]*(1-I.male[i,t])*rrate.tom*II.male[i] + rr.male[2]*(1-I.male[i,t])*(1-II.male[i])
    } #t
  } #i

  ## Informative priors for hunter reporting rates
  rrate.jake ~ dnorm(0.71, sd = 0.072) #  non-reward bands for juveniles
  rrate.tom ~ dnorm(0.87, sd = 0.039) #  non-reward bands for adults

  # Proportion mortality due to hunting
  rr.male[1] ~ dunif(0, 1)  # juveniles # rr*(1-S) = harvest rate r in seber param
  rr.male[2] ~ dunif(0, 1)  # adults


  ## time effect
  for(t in 1:(n.occasions.male-1)){
    beta.time.male[t] ~ dnorm(0, sd = 0.5)
  }

  # WMU effect
  ## Random effect of wmu
  sigma.male ~ dunif(0, 10)
  for(u in 1:n.wmu){
    gamma.male[u] ~ dnorm(0, sd = sigma.male)
  }

  # Juvenile effect on survival
  alpha.male ~ dnorm(0, sd = 0.5)
  ###########################################################X
  # DRM: Derived Estimates
  ###########################################################X
  #----------------------------------------------------------#
  # DRM: Derived STATEWIDE survival rates
  #----------------------------------------------------------#
  for (t in 1:(n.occasions.male-1)) {
    logit(s.ad.male[t]) <- inprod(beta.time.male[1:(n.occasions.male-1)], time.male[t, 1:(n.occasions.male-1)])
    logit(s.jv.male[t]) <- alpha.male + inprod(beta.time.male[1:(n.occasions.male-1)], time.male[t, 1:(n.occasions.male-1)]) # J = means parameterization
  }#t
  #----------------------------------------------------------#
  # DRM: Derived MEAN STATEWIDE harvest rate estimates
  #----------------------------------------------------------#
  logit(mean.s.ad.male) <- beta.time.male[1]
  logit(mean.s.jv.male) <- beta.time.male[1] + alpha.male

  #----------------------------------------------------------#
  # DRM: Derived STATEWIDE harvest rate estimates
  #----------------------------------------------------------#
  for (t in 1:(n.occasions.male-1)) {
    h.juv.male[t] <- (1-s.jv.male[t])*rr.male[1]
    h.ad.male[t] <- (1-s.ad.male[t])*rr.male[2]
  }#t
  #----------------------------------------------------------#
  # DRM: Derived MEAN STATEWIDE harvest rate estimates
  #----------------------------------------------------------#
  mean.h.jv.male <- (1-mean.s.jv.male)*rr.male[1]
  mean.h.ad.male <- (1-mean.s.ad.male)*rr.male[2]

  #----------------------------------------------------------#
  # DRM: WMU survival rates
  #----------------------------------------------------------#
  for (t in 1:(n.occasions.male-1)) {
    for (u in 1:n.wmu) {
      logit(s.ad.wmu.male[t,u]) <- inprod(beta.time.male[1:(n.occasions.male-1)], time.male[t, 1:(n.occasions.male-1)]) + gamma.male[u]
      logit(s.juv.wmu.male[t,u]) <- alpha.male + inprod(beta.time.male[1:(n.occasions.male-1)], time.male[t, 1:(n.occasions.male-1)]) + gamma.male[u]
    }#g
  }#t

  #----------------------------------------------------------#
  # DRM: WMU harvest rates
  #----------------------------------------------------------#
  for (t in 1:(n.occasions.male-1)) {
    for (g in 1:n.wmu) {
      h.juv.wmu.male[t,g] <- (1-s.juv.wmu.male[t,g])*rr.male[1]
      h.ad.wmu.male[t,g] <- (1-s.ad.wmu.male[t,g])*rr.male[2]
    }#g
  }#t
  #----------------------------------------------------------#
  # DRM: Likelihood
  #----------------------------------------------------------#
  for (i in 1:nind.male){
    # Define latent state at first capture
    z.male[i,f.male[i]] <- 1
    for (t in (f.male[i]+1):n.occasions.male){
      # State process
      mu1.male[i,t] <- s.male[i,t-1] * z.male[i,t-1]
      z.male[i,t] ~ dbern(mu1.male[i,t])
      # Observation process
      mu2.male[i,t] <- r.male[i,t-1] * (z.male[i,t-1] - z.male[i,t])
      y.male[i,t] ~ dbern(mu2.male[i,t])
    } #t
  } #i
})
