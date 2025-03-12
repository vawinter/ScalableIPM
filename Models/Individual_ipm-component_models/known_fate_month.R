#############################################################################X
## Formatting Nimble model ----
#############################################################################X
survival <- nimbleCode({
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
  # # Known fate model: Likelihood ----
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
  # Monthly survival per age, wmu, and month ----
  #----------------------------------------------------------#
  # Initialize storage for WMU1 (Adult and Juvenile) across all months
  for (m in 1:num_months) {
    female.s.kf[1,1,m] <- icloglog(telem.beta.int[1])  # WMU1 + Adult
    female.s.kf[1,2,m] <- icloglog(telem.beta.int[1] + telem.beta.age[1])  # WMU1 + Juvenile
  }
  
  # Initialize storage as a 3D array: WMU x Age x Month
  for (j in 1:female.telem.wmu-1) {
    for (m in 1:num_months) {
    # Adult in wmu + 1
      female.s.kf[j+1,1, m] <- icloglog(telem.beta.int[1] + telem.beta.wmu[j])
    # Juvenile in wmu + 1
      female.s.kf[j+1,2, m] <- icloglog(telem.beta.int[1] + telem.beta.age[1] + telem.beta.wmu[j])
  }
}

})