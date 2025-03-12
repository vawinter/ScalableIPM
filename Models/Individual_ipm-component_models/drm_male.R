drm <- nimbleCode({
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
        gamma.wmu.male[wmu.male[i]]
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
  randome.effect.wmu.male ~ dunif(0, 10)
  for(u in 1:n.wmu){
    gamma.wmu.male[u] ~ dnorm(0, sd = randome.effect.wmu.male)
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
      logit(s.ad.wmu.male[t,u]) <- inprod(beta.time.male[1:(n.occasions.male-1)], time.male[t, 1:(n.occasions.male-1)]) + gamma.wmu.male[u]
      logit(s.juv.wmu.male[t,u]) <- alpha.male + inprod(beta.time.male[1:(n.occasions.male-1)], time.male[t, 1:(n.occasions.male-1)]) + gamma.wmu.male[u]
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
