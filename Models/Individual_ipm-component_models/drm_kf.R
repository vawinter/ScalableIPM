drm_kf <- nimbleCode({
  ###########################################################X
  # # Known-fate model ----
  # ###########################################################X
  # ###########################################################X
  # # Known-fate Model: Priors on Intercept, Age, and WMU ----
  # ###########################################################X
  # # Intercept
  # telem.beta.int[1] ~ dnorm(0, sd = 0.5)
  # 
  # # Age effect
  # telem.beta.age[1] ~ dnorm(0, sd = 0.5)
  # 
  # ###########################################################X
  # # WMU random effect  ----
  # ###########################################################X
  # telem.sigma ~ dunif(0, 5)
  # for (w in 1:female.telem.wmu) {
  #   telem.beta.wmu[w] ~ dnorm(0, sd = telem.sigma)
  # }
  # 
  # ###########################################################X
  # # Month random effect  ----
  # ###########################################################X
  # telem.month.sigma ~ dunif(0, 5)
  # for (m in 1:12) {
  #   telem.beta.month[m] ~ dnorm(0, sd = telem.month.sigma)
  # }
  # 
  # ###########################################################X
  # # Known Fate Model: Likelihood ----
  # ###########################################################X
  # 
  # # Loop over all individuals in the dataset
  # for (i in 1:telem.nind) {  
  #   # Loop over the range of years in which each individual was tracked
  #   for (y in telem.year.start[i]:telem.year.end[i]) {  
  #     
  #     # Loop over the months between the first and last capture for each individual
  #     for (m in telem.first[i]:telem.last[i]) {  
  #       #-------------------------------------------------------------------------X
  #       # Note:
  #       # Derived survival rate on the complementary log-log (cloglog) scale
  #       # s.kf[i, y, m]: Survival for individual i in year y and month m
  #       # telem.beta.int[1]: Intercept term: adults
  #       # telem.beta.age[1]: Effect of age (juvenile) on survival
  #       # telem.juvenile[i, y, m]: Indicator for juvenile status (1 if juvenile, 0 if adult)
  #       # inprod(telem.beta.wmu[1:4], telem.wmu[i, 1:4]): Captures spatial random effect for WMU
  #       # telem.beta.month[m]: Monthly survival effect to capture seasonal variation
  #       #-------------------------------------------------------------------------X
  #       
  #       s.kf[i, y, m] <-  telem.beta.int[1] +   # Adult
  #         telem.beta.age[1] * telem.juvenile[i, y, m] +  # Juvenile-specific effect
  #         inprod(telem.beta.wmu[1:female.telem.wmu], telem.wmu[i, 1:female.telem.wmu]) +  # WMU random effect
  #         telem.beta.month[m]   # Month random effect
  #       
  #       #-------------------------------------------------------------------------X
  #       # Likelihood for survival status (0 = dead, 1 = alive) using a Bernoulli distribution
  #       # The probability is the inverse complementary log-log of the survival 
  #       # icloglog() is appropriate when values are skewed low (towards 0) or high (towards 1)
  #       #-------------------------------------------------------------------------X
  #       
  #       status[i, y, m] ~ dbern(prob = icloglog(s.kf[i, y, m]))
  #     }
  #   }
  # }
  # 
  # ###########################################################X
  # # Derived Quantities: Survival over month per age and WMU ----
  # ###########################################################X
  # 
  # # Initialize storage for monthly survival rates (Adult and Juvenile)
  # for (m in 1:12) {
  #   # Adults, baseline WMU
  #   storage[1, 1, m] <- icloglog(telem.beta.int[1] + telem.beta.wmu[1] +  telem.beta.month[m])   
  #   # Juveniles, baseline WMU
  #   storage[1, 2, m] <- icloglog(telem.beta.int[1] + telem.beta.age[1] + telem.beta.wmu[1] +  telem.beta.month[m]) 
  # }
  # 
  # # For WMUs beyond the first
  # for (j in 2:(female.telem.wmu)) {
  #   for (m in 1:12) {
  #     # Adults
  #     storage[j, 1, m] <- icloglog(telem.beta.int[1] + telem.beta.wmu[j] + telem.beta.month[m])
  #     # Juveniles
  #     storage[j, 2, m] <- icloglog(telem.beta.int[1] + telem.beta.age[1] + telem.beta.wmu[j] + telem.beta.month[m])
  #   }
  # }
  # #---------------------------------------------------------------X
  # # Calculate survival for Adults, accounting for nesting season 
  # #---------------------------------------------------------------X
  # for (j in 1:female.telem.wmu) {
  #   # Jan to December survival
  #   avg.ad.s.kf[j] <- prod(storage[j, 1, 1:12])  
  # }
  # #---------------------------------------------------------------X
  # # Calculate survival for Juvenile FEMALES transitioning to Adults in June
  # #---------------------------------------------------------------X
  # for (j in 1:female.telem.wmu) {
  #   # November and December survival using Jan as a proxy
  #   juvenile_part1[j] <- storage[j, 2, 1] * storage[j, 2, 1]     
  #   # January to May
  #   juvenile_part2[j] <- prod(storage[j, 2, 1:5])   
  #   # June to December (aging up to adult)
  #   adult_part[j] <- prod(storage[j, 1, 6:10])     
  #   
  #   # Cumulative survival for each entry month
  #   avg.juv.s.kf[j] <- (juvenile_part1[j] * juvenile_part2[j])* adult_part[j] 
  # }
  # 
  # #---------------------------------------------------------------X
  # # Calculate survival for Juvenile MALES surviving until Spring harvest 
  # #---------------------------------------------------------------X
  # for (j in 1:male.n.wmu) { 
  #   # November and December survival using Jan as a proxy
  #   juvenile_male_1[j] <- storage[j, 2, 1] * storage[j, 2, 1]     
  #   # January to April
  #   juvenile_male_2[j] <- prod(storage[j, 2, 1:4])   
  #   
  #   # Cumulative survival until Spring harvest 
  #   juv.male.adj[j] <- (juvenile_male_1[j] * juvenile_male_2[j])
  # }
  # 
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

})
