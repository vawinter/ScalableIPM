###########################################################X
# Abundance: Derived estimates via Lincoln-Peterson Est. 
###########################################################X
# Estimating population abundance by age class and season for subsequent occasions
# per WMU
for (t in 1:true.occasions) {
  for(u in 1:male.n.wmu){ # female.wmu and male.wmu = same
    ###########################################################X
    # N: Male model ----
    ###########################################################X
    #-----------#
    # May 
    #-----------#
    male.N.ad.may[t, u] <- harvest.ad.spring[t, u]/(male.h.ad.wmu[t, u]) 
    male.N.juv.may[t, u] <- harvest.juv.spring[t, u]/(male.h.juv.wmu[t, u]) 
    #-----------#
    # June
    #-----------#
    male.N.ad.june[t, u] <-(male.N.ad.may[t, u] - harvest.ad.spring[t,u]) + (male.N.juv.may[t, u] - male.h.juv.wmu[t, u])
    #-----------#
    # Sept 
    #-----------#
    male.N.ad.sept[t, u] <-  male.N.ad.june[t, u] * (1 - male.s.ad.wmu[t, u] - male.h.ad.wmu[t, u])
    # Number of females on Aug 31 (Sept 1) * Number of hens w. brood (hwb) * number poults per brood (ppb)
    male.N.juv.sept[t, u] <- (female.N.ad.sept[t, u] * hwb.p[t] * ph.mu[t])/2 
    #-----------#
    # Nov 
    #-----------#
    male.N.ad.nov[t, u] <- male.N.ad.sept[t, u] # S = 1 for males 
    male.N.juv.nov[t, u] <- (male.N.juv.sept[t, u]) + (harvest.juv.fall[t, u]/male.h.juv.wmu[t, u])
    
    ###########################################################X
    # N: Feale model ----
    ###########################################################X
    #-----------#
    # May 
    #-----------#
    female.N.ad.may[t, u] <-  male.N.ad.may[t, u]*2
    female.N.juv.may[t, u] <- male.N.juv.may[t, u] 
    #-----------#
    # June
    #-----------#
    female.N.ad.june[t, u] <- female.N.ad.may[t, u] + female.N.juv.may[t, u] # S  = 1
    #-----------#
    # Sept 
    #-----------#
    female.N.ad.sept[t, u] <-  female.N.ad.june[t, u]
    # Number of females on Aug 31 (Sept 1) * Number of hens w. brood (hwb) * number poults per brood (ppb)
    female.N.juv.sept[t, u] <- (female.N.ad.sept[t, u] * hwb.p[t] * ph.mu[t])/2 
    #-----------#
    # Nov 
    #-----------#
    female.N.ad.nov[t, u] <-   (harvest.ad.fall[t, u]/female.h.ad.wmu[t, u]) 
    female.N.juv.nov[t, u] <- female.N.juv.sept[t, u] + (harvest.juv.fall[t, u]/female.h.juv.wmu[t, u])
  } # u
} # t
    