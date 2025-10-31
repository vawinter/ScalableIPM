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
  
  # Add to priors section
  telem.year.sigma ~ dunif(0, 5)
  for (y in 1:n_years) {
    telem.beta.year[y] ~ dnorm(0, sd = telem.year.sigma)
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
        # step(m - telem.first[i]): holds the first month at the intercept
        #-------------------------------------------------------------------------X
        
        s.kf[i, y, m] <-  telem.beta.int[1] +   # Adult
          telem.beta.age[1] * telem.juvenile[i, y, m] +  # Juvenile-specific effect
          inprod(telem.beta.wmu[1:female.telem.wmu], telem.wmu[i, 1:female.telem.wmu]) +  # WMU random effect
          telem.beta.month[m] +  # Month random effect
          telem.beta.year[y]  # NEW
        
        #-------------------------------------------------------------------------X
        # Likelihood for survival status (0 = dead, 1 = alive) using a Bernoulli distribution
        # The probability is the inverse complementary log-log of the survival 
        # icloglog() is appropriate when values are skewed low (towards 0) or high (towards 1)
        #-------------------------------------------------------------------------X
        
        status[i, y, m] ~ dbern(prob = icloglog(s.kf[i, y, m]))
      }
    }
  }
  
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
  # # Calculate survival for Juveniles transitioning to Adults in June
  # #---------------------------------------------------------------X
  # for (j in 1:female.telem.wmu) {
  #   # November and December survival using Jan as a proxy
  #   juvenile_part1[j] <- storage[j, 2, 1] * storage[j, 2, 1]     
  #   # January to May
  #   juvenile_part2[j] <- prod(storage[j, 2, 1:5])   
  #   # June to December (aging up to adult)
  #   adult_part[j] <- prod(storage[j, 1, 6:10])     # changed from 6:12.
  #   
  #   # Cumulative survival for each entry month
  #   avg.juv.s.kf[j] <- (juvenile_part1[j] * juvenile_part2[j])* adult_part[j] 
  # }
  # 
  # #---------------------------------------------------------------X
  # # Calculate survival for Juvenile MALES surviving until Spring harvest 
  # #---------------------------------------------------------------X
  # for (j in 1:male.n.wmu) { #exclude 5C
  #   # November and December survival using Jan as a proxy
  #   juvenile_male_1[j] <- storage[j, 2, 1] * storage[j, 2, 1]     
  #   # January to April
  #   juvenile_male_2[j] <- prod(storage[j, 2, 1:4])   
  #   
  #   # Cumulative survival until Spring harvest 
  #   juv.male.adj[j] <- (juvenile_male_1[j] * juvenile_male_2[j])
  # }
  ###########################################################X
  # Derived Quantities: Survival over month per age, WMU, and YEAR ----
  ###########################################################X
  
  # Now storage needs a year dimension: storage[wmu, age, month, year]
  for (y in 1:n_years) {
    for (m in 1:12) {
      # Adults, baseline WMU, year y
      storage[1, 1, m, y] <- icloglog(telem.beta.int[1] + 
                                        telem.beta.wmu[1] + 
                                        telem.beta.month[m] + 
                                        telem.beta.year[y])
      
      # Juveniles, baseline WMU, year y
      storage[1, 2, m, y] <- icloglog(telem.beta.int[1] + 
                                        telem.beta.age[1] + 
                                        telem.beta.wmu[1] + 
                                        telem.beta.month[m] + 
                                        telem.beta.year[y])
    }
    
    # For WMUs beyond the first
    for (j in 2:female.telem.wmu) {
      for (m in 1:12) {
        # Adults
        storage[j, 1, m, y] <- icloglog(telem.beta.int[1] + 
                                          telem.beta.wmu[j] + 
                                          telem.beta.month[m] + 
                                          telem.beta.year[y])
        # Juveniles
        storage[j, 2, m, y] <- icloglog(telem.beta.int[1] + 
                                          telem.beta.age[1] + 
                                          telem.beta.wmu[j] + 
                                          telem.beta.month[m] + 
                                          telem.beta.year[y])
      }
    }
  }
  
  #---------------------------------------------------------------X
  # Calculate ANNUAL survival by WMU and YEAR
  #---------------------------------------------------------------X
  
  # Adult survival by WMU and year
  for (j in 1:female.telem.wmu) {
    for (y in 1:n_years) {
      # Jan to December survival for year y
      avg.ad.s.kf[y, j] <- prod(storage[j, 1, 1:12, y])
    }
  }
  
  # Juvenile survival by WMU and year
  for (j in 1:female.telem.wmu) {
    for (y in 1:n_years) {
      # Nov-Dec (using Jan as proxy)
      juvenile_part1[y, j] <- storage[j, 2, 1, y] * storage[j, 2, 1, y]
      
      # Jan-May as juvenile
      juvenile_part2[y, j] <- prod(storage[j, 2, 1:5, y])
      
      # Jun-Oct as adult (aging up)
      adult_part[y, j] <- prod(storage[j, 1, 6:10, y])
      
      # Total juvenile survival for year y
      avg.juv.s.kf[y, j] <- (juvenile_part1[y, j] * juvenile_part2[y, j]) * adult_part[y, j]
    }
  }
  
  # Juvenile male survival to spring by WMU and year
  for (j in 1:male.n.wmu) {
    for (y in 1:n_years) {
      # Nov-Dec
      juvenile_male_1[y, j] <- storage[j, 2, 1, y] * storage[j, 2, 1, y]
      
      # Jan-Apr
      juvenile_male_2[y, j] <- prod(storage[j, 2, 1:4, y])
      
      # Cumulative to spring harvest
      juv.male.adj[y, j] <- (juvenile_male_1[y, j] * juvenile_male_2[y, j])
    }
  }

})

