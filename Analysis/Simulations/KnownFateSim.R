##########################################################
#
# Simulate known fate data to better understand the effect of link function
#  on point estimates
#
# Duane R. Diefenbach, April 2024
##########################################################

library(R2jags)
library(dplyr)
library(mcmcplots)


sim <- 1000      # number of replicate datasets
smplsz <- 100  # number at risk each time period
rate <- 0.99   # monthly survival rate


######## vector filled with simulated data, assume fixed number at risk and survival rate
 MR <- rbinom(sim,smplsz,rate)
 

 jags.data <- list(x=MR,y=MR,z=MR, s=sim, n=smplsz)
 
 # Initial values
 inits <- function(){list(theta=rnorm(sim,0,1.3) )}
 
 ######## Call JAGS from R (BRT 4 min)
 # Parameters monitored
 parameters <- c("p", "ptheta", "plog")
 
 # MCMC settings
 ni <- 25000
 nt <- 1
 nb <- 7500
 nc <- 6
 
 ######## Call JAGS from R (BRT 4 min)
 m1 <- jags(jags.data, inits, parameters, "Binomjags.txt", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb)
 #summary(mr.ss.re)
 #print(m1, digits = 3)
 
print(paste("Mean p =", mean(m1$BUGSoutput$mean$p-rate)," SD = ",
  sd(m1$BUGSoutput$mean$p-rate)))
 mean(m1$BUGSoutput$mean$ptheta-rate)
  sd(m1$BUGSoutput$mean$ptheta-rate)
 mean(m1$BUGSoutput$mean$plog-rate)
  sd(m1$BUGSoutput$mean$plog-rate)
 
 
 mcmcplot(m1)
 