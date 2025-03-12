survival <- nimbleCode({
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
  telem.sigma ~ dunif(0, 10)
  for (w in 1:4) {
    telem.beta.wmu[w] ~ dnorm(0, sd = telem.sigma)
  }
  
  ###########################################################X
  # Month random effect  ----
  ###########################################################X
  telem.month.sigma ~ dunif(0, 10)
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
        # telem.beta.int[1]: Intercept term: adult in WMU 4D
        # telem.beta.age[1]: Effect of age (juvenile) on survival
        # telem.juvenile[i, y, m]: Indicator for juvenile status (1 if juvenile, 0 if adult)
        # inprod(telem.beta.wmu[1:3], telem.wmu[i, 1:3]): Captures spatial random effect for WMU
        # telem.beta.month[m]: Monthly survival effect to capture seasonal variation
        # step(m - telem.first[i]): holds the first month at the intercept
        #-------------------------------------------------------------------------X
        
        s.kf[i, y, m] <-  telem.beta.int[1] +   # Adult, first month in study
          telem.beta.age[1] * telem.juvenile[i, y, m] +  # Juvenile-specific effect
          inprod(telem.beta.wmu[1:4], telem.wmu[i, 1:4]) +  # WMU random effect
          telem.beta.month[m]   # Month effects only beyond the first month
        
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
    storage[1, 1, m] <- icloglog(telem.beta.int[1] + telem.beta.wmu[1] +  telem.beta.month[m]) #* step(m - 1))  
    # Juveniles, baseline WMU
    storage[1, 2, m] <- icloglog(telem.beta.int[1] + telem.beta.age[1] + telem.beta.wmu[1] +  telem.beta.month[m]) # * step(m - 1)) 
  }
  
  # For WMUs beyond the first
  for (j in 1:(female.telem.wmu - 1)) {
    for (m in 1:12) {
      # Adults
      storage[j + 1, 1, m] <- icloglog(telem.beta.int[1] + telem.beta.wmu[j] + telem.beta.month[m])
      # Juveniles
      storage[j + 1, 2, m] <- icloglog(telem.beta.int[1] + telem.beta.age[1] + telem.beta.wmu[j] + telem.beta.month[m])
    }
  }
  #---------------------------------------------------------------X
  # Calculate survival for Adults, accounting for nesting season 
  #---------------------------------------------------------------X
  for (j in 1:female.telem.wmu) {
    # Jan to December
    avg.ad.s.kf[j] <- prod(storage[j, 1, 1:12])  
  }
  #---------------------------------------------------------------X
  # Calculate survival for Juveniles transitioning to Adults in June
  #---------------------------------------------------------------X
  for (j in 1:female.telem.wmu) {
    # November and December survival using Jan as a proxy
    juvenile_part1[j] <- storage[j, 2, 1] * storage[j, 2, 1]     
    # January to May
    juvenile_part2[j] <- prod(storage[j, 2, 1:5])   
    # June to November
    adult_part[j] <- prod(storage[j, 1, 6:12])     
    
    # Cumulative survival for each entry month
    avg.juv.s.kf[j] <- (juvenile_part1[j] * juvenile_part2[j])* adult_part[j] 
  }
  

})

