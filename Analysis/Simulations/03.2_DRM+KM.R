###############################################################################X
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#################### Integrated Population Model (IPM): #######################X
#                 #---# PhD Dissertation: Chapter 1 #---#
#        Creating a Bayesian IPM to inform turkey management in PA
###                *** SIMULATION TEST (Nimble, MF, kf) ***                 ###X
###                                                                         ###X
#    Known-fate Modeling: Creating a  known-fate model for telemetered hens by 
# wildlife management unit (WMU), and age class over months prior to
#                            fall hunting season.
#
#   DRM Modeling: Creating a  dead-recovery model for males and females by 
# wildlife management unit (WMU), and age class over seasons/years
###                                                                         ###X
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #X
# 
# Created by: Veronica A. Winter
# Last edited: 06/05/2024
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #X
###############################################################################X
# Clean env
rm(list = ls())
gc()

# Set options
options(max.print=1000000)

# Call functions
source("Analysis/Scripts/Simulations/Simulation dead recovery/00_funs.R")

# Libraries
library(nimble)
library(MCMCvis)
library(Rlab)

##################################################X
# Known-fate model (KF) Simulations ----
#############################################################################X
# Simulating data: checking priors and model fit ----
#############################################################################X
## Parameters
telem.nind <- 400  # Number of individuals
nMonths <- 48  # Number of months
telem.beta.age <- c(0.3, 0.5)  # Effects of being adult or juvenile
telem.beta.wmu <- c(0.2, 0.3, 0.15, 0.5)  # Effect of different WMUs
survival_status <- matrix(NA, nrow = telem.nind, ncol = nMonths)

# Simulate age groups and WMU assignments
is_adult <- matrix(sample(0:1, telem.nind * nMonths, replace = TRUE), nrow = telem.nind, ncol = nMonths)
telem.juvenile <- 1 - is_adult  

# Fill in WMU matrix and include intercept
telem.wmu <- matrix(0, telem.nind*4, nrow = telem.nind, ncol = 4)
telem.wmu[1:125, 1] <- 1
telem.wmu[125:250, 2] <- 1
telem.wmu[251:375, 3] <- 1
telem.wmu[, 4] <- 1 # intercept

# Initialize first and last matrices 
telem.first <- sample(1:4, telem.nind, replace = TRUE)
telem.last <- NULL

# Initialize matrices to store survival data for each age class
survival_adult <- NULL
survival_juvenile <- NULL
cll_survival <-  NULL
survival_wmu <-  NULL

#############################################################################X
## Simulating data: simulating status and survival probabilities ----
#############################################################################X
for(i in 1:telem.nind) {
  # setting so I can make sure that it gets 'overwritten'
  survival_status[i, telem.first[i]] <- 9
  
  for (t in telem.first[i]:nMonths) {
    # This is my estimation: intercepts for 'adult' and '4D'
    cll_survival[i] <- (telem.beta.age[1] + telem.beta.wmu[4]) + telem.beta.age[2] * telem.juvenile[i, t] + 
      (telem.beta.wmu[1:3]%*%telem.wmu[i, 1:3])
    
    # this is creating my survival status, because they need to come from a DISTRIBUTION
    survival_status[i, t] <- rbern(n=1, prob = icll(cll_survival[i]))
    
    # Check if individual has 'died'
    if (survival_status[i, t] == 0) {
      telem.last[i] <- t
      if (t < nMonths) {
        survival_status[i, (t+1):nMonths] <- NA
      }
      break  # Exit the loop for this individual
    }
    
    # Fill age-specific survival matrices
    if (is_adult[i, t] == 1) {
      survival_adult[t] <- cumprod(icll(cll_survival[i]))
    } else if (telem.juvenile[i, t] == 1) {
      survival_juvenile[t] <- cumprod(icll(cll_survival[i]))
    }
    
  } # end t
} # end i

#############################################################################X
### Simulating data: Check model outputs ----
#############################################################################X
survival_adult[1:8]
survival_juvenile[1:8]

# Fill in NAs in 'last' with the max month
telem.last <- ifelse(is.na(telem.last),max(nMonths), telem.last)
length(telem.last) # Check dimensions
# last[1000] <- nMonths # add if dim < nind


# Note: dbern in nimble says 'this is a pdf for max likelihood' and
# dbern in Bayesian says 'this arises from a distribution' 

##################################################X
# Dead-Recovery Model (DRM) Simulations ----
##################################################X
# Set simulation values ----
##################################################X
rrate.j <- 0.7 # reporting rate juveniles
rrate.a <- 0.9 # reporting rate adults
n.occasions <- 4 # number of occasions
S <- c(0.7, 0.5)  # survival c(juv, ad)
S_kf <- survival_adult[1] # known fate survival for females (assuming j and a are the same)
rr <- c(0.3, 0.5)  # seber recovery rate c(juv, ad)
marked <- 400 # number marked individuals
###-------------------#X
# Random effect for wmu
# raneff <- round(rnorm(10,0,0.05),3)
raneff <- c(-0.043, -0.006, -0.038,  0.056, -0.081,
            -0.094, -0.006,  0.003,  0.084, -0.057) # mean = -0.182, SD = 0.05667
########################################################X
# Function to simulate data for two different sexes ----
###-----------------------------------------------------#X
# Simulate M and F data separately w. same seed 
# that I have been using thus far
#########################################################X
# Males w. set.seed
male_data <- simulate_data(123, marked, n.occasions, S, rr, rrate.j, rrate.a, 
                           raneff, sex = "M")

# Females w. set.seed
female_data <- simulate_data(321, marked, n.occasions, S, rr, rrate.j, rrate.a, 
                             raneff, sex = "F", S_kf = S_kf)

### Check data ----
table(male_data$MR)
table(female_data$MR)
###-----------------------------------------------------#X
table(male_data$J)
table(female_data$J)
###-----------------------------------------------------#X
table(male_data$f)
table(female_data$f)
###-----------------------------------------------------#X
table(male_data$wmu)
table(female_data$wmu)
###-----------------------------------------------------#X
# F[i,t] indicator = known-fate (f.kf) 
# first capture occasion (matrix)
head(female_data$f.kf)
# remove column names
colnames(female_data$f.kf) <- NULL
###-----------------------------------------------------#X
##################################################X
### Rename variables ----
##################################################X
# Male data
male.y <- male_data$MR
male.I <- male_data$I
male.II <- male_data$II
male.wmu <- male_data$wmu
male.f <- male_data$f
male.time.param <- male_data$J
male.z <- male_data$known.state
male.rrate.j <- male_data$rrate.j
male.rrate.a <- male_data$rrate.a
male.n.occasions <- (dim(male_data$MR)[2])
male.nind <- dim(male_data$MR)[1]
###-----------------------------------------------------#X
# Female data
female.y <- female_data$MR
female.I <- female_data$I
female.II <- female_data$II
female.wmu <- female_data$wmu
female.f <- female_data$f # also tagging cohort
female.time.param <- female_data$J
female.z <- female_data$known.state
female.rrate.j <- female_data$rrate.j
female.rrate.a <- female_data$rrate.a
female.n.occasions <- (dim(female_data$MR)[2])
female.nind <- dim(female_data$MR)[1]
female.f.kf <- female_data$f.kf
##-----------------------------------------------------#X
# Create n.occasions-1 variable
true.occasions = n.occasions
##-----------------------------------------------------#X
# Apply function to get the first non-zero value in a row (Juvenile) 
# to each row of the matrix
female.age <- apply(female.I, 1, get_first_non_zero)
# NA = adult, other = Juvenile
# For indexing, I want 1 = adult, 2 = juvenile
female.age <- ifelse(is.na(female.age), 1, 2)
##################################################X
# Estimate parameters in Nimble ----
##################################################X
### Nimble model set up ----
# Data
nimble.data <- list(
  # KF telemetered data 
  status = survival_status,
  telem.juvenile = telem.juvenile, # adult is intercept
  telem.wmu = telem.wmu, # 4D is intercept term
  ###-----------#X
  # DRM Male data
  male.y = male.y, 
  male.time.param = male.time.param,
  ###-----------#X
  # DRM Female data
  female.y = female.y,
  female.time.param = female.time.param
)

# Constants
consts <- list(
  # KF telemetered constants 
  telem.nind = telem.nind,
  telem.first = telem.first,
  telem.last = telem.last,
  female.telem.wmu = 4, # 4 wmus but 1 is captured in intercept
  ###-----------#X
  # DRM Male constants
  male.f = male.f,
  male.I = male.I,
  male.II = male.II,
  male.nind = male.nind,
  male.n.occasions = male.n.occasions,
  male.n.wmu = 10,
  male.wmu = male.wmu,
  ###-----------#X
  # DRM Female constants
  female.f = female.f,
  female.I = female.I,
  female.II = female.II,
  female.nind = female.nind,
  female.n.occasions = female.n.occasions,
  female.n.wmu = 10,
  female.wmu = female.wmu,
  female.f.kf = female.f.kf,
  tagging.cohort = female.f, # tagging cohort and first capture are the same, but wanted to make its own index for clarity
  female.age = female.age, # age at first capture
  ###-----------#X
 # set.s.kf = survival_adult[1], # known fate survival for females (assuming j and a are the same)
  true.occasions = true.occasions
)

# Initial values
inits <- list(
  # KF telemetered inits
  telem.beta.int = rnorm(1, 0, 0.5),
  telem.beta.age = rnorm(1, 0, 0.5),
  telem.beta.wmu = rnorm(3, 0, 1),
  ###-----------#X
  # DRM Male inits
  male.z = male.z,
  male.juvenile.effect=rnorm(1,0,0.5),
  male.time.effect = rnorm((n.occasions), 0, 0.5), 
  male.seber.recov=runif(2,0,1),
  male.rrate.a = rnorm(1,0.87,0.039), 
  male.rrate.j = rnorm(1,0.71,0.072), 
  male.sigma=runif(1,0,10), 
  ###-----------#X
  # DRM Female inits
  female.z = female.z,
  female.juvenile.effect=rnorm(1,0,0.5),
  female.time.effect=rnorm((true.occasions),0,0.5),
  female.seber.recov=runif(2,0,1),
  female.rrate.a=runif(1,0,1),
  female.rrate.j=runif(1,0,1), 
  female.sigma=runif(1,0,10)
)

# ##################################################X
# ## Call Nimble model ----
# ##################################################
# # garbage clean up
# rm(female_data, male_data, cll_survival)
# gc()
# ###-----------------------------------------------------#X
# # source model
# source("models/drm_km.R")
# ###-----------------------------------------------------#X
# # Create the Nimble model
# model <- nimbleModel(
#   code = drm_kf,
#   data = nimble.data,
#   constants = consts,
#   inits = inits
# )
# ###-----------#X
# # Check initialization issues
# model$initializeInfo()
# gc()
# ###-----------#X
# # Set MCMC specifications
# ni <- 20000
# nt <- 1
# nb <- 7500
# nc <- 2 # minimum 2 chains for rhat values
# ###-----------#X
# # Run MCMC 
# c7_results <- nimbleMCMC(
#   model = model,
#   monitors = c("telem.beta.wmu", "telem.beta.int", "telem.beta.age", "female.s.kf",
#                "s.kf", "adj.female.s.kf",
#                ###------------------------------------------------------------#X
#                 "male.juvenile.effect","male.time.effect","male.wmu.effect",
#                 "male.rrate.a","male.rrate.j", "male.mean.s.ad","male.mean.s.jv",
#                 "male.mean.harv.ad", "male.mean.harv.jv","male.seber.recov",
#                 #"male.h.juv.wmu","male.h.ad.wmu", "male.s.juv.wmu","male.s.ad.wmu",
#                ###------------------------------------------------------------#X
#                "female.juvenile.effect","female.time.effect","female.wmu.effect",
#                "female.rrate.a","female.rrate.j", "female.mean.s.ad",
#                "female.mean.s.jv","female.mean.harv.ad","female.mean.harv.jv",
#                "female.seber.recov"#,"female.h.juv.wmu","female.h.ad.wmu",
#                #"female.s.juv.wmu","female.s.ad.wmu"
#                ),
#   niter = ni, 
#   nburnin = nb, 
#   thin = nt, 
#   nchains = nc, 
#   setSeed = FALSE, 
#   samplesAsCodaMCMC = TRUE,
#   WAIC = TRUE)
# # beepr::beep(1)
# gc()
# 
# #############################################################X
# #           Model diagnostics ----
# #############################################################X
# summary(c7_results$samples)
# MCMCvis::MCMCsummary(c7_results$samples) # to look at Rhats
# c7_results$WAIC
# MCMCvis::MCMCtrace(c7_results$samples, iter = ni)
# #############################################################X
# ## Known-fate output ----
# #############################################################X
# MCMCtrace(c7_results$samples$chain1, pdf = F, iter = ni, params = c("telem.beta.age"))
# MCMCtrace(c7_results$samples$chain1, pdf = F, iter = ni, params = c("telem.beta.int"))
# MCMCtrace(c7_results$samples$chain1, pdf = F, iter = ni, params = c("telem.beta.wmu"))
# 
# ### Extract probabilities ---
# MCMCpstr(c7_results$samples, params = c("adj.female.s.kf"))
# ### Extract probabilities ---
# MCMCpstr(c7_results$samples, params = c("female.s.kf"))
# MCMCtrace(c7_results$samples$chain1, pdf = F, iter = ni, params = c("female.s.kf"))
# 
# samplesSummary(c7_results$samples$chain1)
# #############################################################X
# # Intercept
# par(mfrow = c(1,2))
# matplot(c7_results$samples[,"telem.beta.int[1]"]$chain2,type="l", main = "Intercept")
# abline(h = telem.beta.age[1]+telem.beta.wmu[4], col = "darkmagenta", lwd = 2)
# plot(density(c7_results$sample[,"telem.beta.int[1]"]$chain2), main = "Density of Beta Coefficient",
#      xlab = "Beta Coefficient", ylab = "Density")
# abline(v = telem.beta.age[1]+telem.beta.wmu[4], col = "darkmagenta", lwd = 2)
# ###----------------------------------------------------#X
# # beta age
# par(mfrow = c(1,2))
# matplot(c7_results$samples[,"telem.beta.age[1]"]$chain1,type="l", main = "Juvenile")
# abline(h = telem.beta.age[2], col = "darkmagenta", lwd = 2)
# plot(density(c7_results$sample[,"telem.beta.age[1]"]$chain1), main = "Density of Beta Coefficient",
#      xlab = "Beta Coefficient", ylab = "Density")
# abline(v = telem.beta.age[2], col = "darkmagenta", lwd = 2)
# ###----------------------------------------------------#X
# # wmu[1]
# par(mfrow = c(3,2))
# matplot(c7_results$samples[,"telem.beta.wmu[1]"]$chain1,type="l", main = "WMU 1")
# abline(h = telem.beta.wmu[1], col = "darkmagenta", lwd = 2)
# plot(density(c7_results$sample[,"telem.beta.wmu[1]"]$chain1), main = "Density of Beta Coefficient",
#      xlab = "Beta Coefficient", ylab = "Density")
# abline(v = telem.beta.wmu[1], col = "darkmagenta", lwd = 2)
# ###----------------------------------------------------#X
# # wmu[2]
# #par(mfrow = c(1,2))
# matplot(c7_results$samples[,"telem.beta.wmu[2]"]$chain1,type="l", main = "WMU 2")
# abline(h = telem.beta.wmu[2], col = "darkmagenta", lwd = 2)
# plot(density(c7_results$sample[,"telem.beta.wmu[2]"]$chain1), main = "Density of Beta Coefficient",
#      xlab = "Beta Coefficient", ylab = "Density")
# abline(v = telem.beta.wmu[2], col = "darkmagenta", lwd = 2)
# ###----------------------------------------------------#X
# # wmu[3]
# #par(mfrow = c(1,2))
# matplot(c7_results$samples[,"telem.beta.wmu[3]"]$chain1,type="l", main = "WMU 3")
# abline(h = telem.beta.wmu[3], col = "darkmagenta", lwd = 2)
# plot(density(c7_results$sample[,"telem.beta.wmu[3]"]$chain1), main = "Density of Beta Coefficient",
#      xlab = "Beta Coefficient", ylab = "Density")
# abline(v = telem.beta.wmu[3], col = "darkmagenta", lwd = 2)
# 
# #############################################################X
# ## DRM output ----
# #############################################################X
# ##         Look at specific estimates
# #         Males = red, Females = blue
# #############################################################X
# # Male (red)
# #############################################################X
# # recovery
# par(mfrow = c(3,2))
# for(i in 1:length(rr)){
#   samp <- paste0("male.seber.recov[", i, "]")
#   matplot(c7_results$samples$chain1[,samp],type="l", main = paste0("male.seber.recov"))
#   abline(h = rr[i], col = "red", lwd = 2)
#   plot(density(c7_results$samples$chain1[,samp]), main = "Density",
#        xlab = "Beta Coefficient", ylab = "Density")
#   abline(v = rr[i], col = "red", lwd = 2)
# }
# ###----------------------------------------------------#X
# par(mfrow=c(1,2))
# # Survival c(juv, ad)
# matplot(c7_results$samples$chain1[,"male.mean.s.jv[1]"],type="l", main = "Male mean survival (J)")
# abline(h = S[1], col = "red", lwd = 2)
# plot(density(c7_results$samples$chain1[,"male.mean.s.jv[1]"]), main = "Density",
#      xlab = "Beta Coefficient", ylab = "Density")
# abline(v = S[1], col = "red", lwd = 2)
# ###----------------------------------------------------#X
# matplot(c7_results$samples$chain1[,"male.mean.s.ad[2]"],type="l", main = "Male mean survival (A)")
# abline(h = S[2], col = "red", lwd = 2)
# plot(density(c7_results$samples$chain1[,"male.mean.s.ad[2]"]), main = "Density",
#      xlab = "Beta Coefficient", ylab = "Density")
# abline(v = S[2], col = "red", lwd = 2)
# #############################################################X
# # Random effect on WMU
# par(mfrow = c(3,2))
# for(i in 1:length(raneff)){
#   samp <- paste0("male.wmu.effect[", i, "]")
#   matplot(c7_results$samples$chain1[,samp],type="l", main = paste0("Male random effect on wmu"))
#   abline(h = raneff[i], col = "red", lwd = 2)
#   plot(density(c7_results$samples$chain1[,samp]), main = "Density",
#        xlab = "Beta Coefficient", ylab = "Density")
#   abline(v = raneff[i], col = "red", lwd = 2)
# }
# #############################################################X
# # Female (blue)
# #############################################################X
# # recovery
# par(mfrow = c(3,2))
# for(i in 1:length(rr)){
#   samp <- paste0("female.seber.recov[", i, "]")
#   matplot(c7_results$samples$chain1[,samp],type="l", main = paste0("female.seber.recov"))
#   abline(h = rr[i], col = "blue", lwd = 2)
#   plot(density(c7_results$samples$chain1[,samp]), main = "Density",
#        xlab = "Beta Coefficient", ylab = "Density")
#   abline(v = rr[i], col = "blue", lwd = 2)
# }
# #############################################################X
# par(mfrow=c(1,2))
# # Survival c(juv, ad)
# ###----------------------------------------------------#X
# matplot(c7_results$samples$chain1[,"female.mean.s.jv[1]"],type="l", main = "Female mean survival (J)")
# abline(h = S[1], col = "blue", lwd = 2)
# plot(density(c7_results$samples$chain1[,"female.mean.s.jv[1]"]), main = "Density",
#      xlab = "Beta Coefficient", ylab = "Density")
# abline(v = S[1], col = "blue", lwd = 2)
# ###----------------------------------------------------#X
# matplot(c7_results$samples$chain1[,"female.mean.s.ad[1]"],type="l", main = "Female mean survival (A)")
# abline(h = S[2], col = "blue", lwd = 2)
# plot(density(c7_results$samples$chain1[,"female.mean.s.ad[1]"]), main = "Density",
#      xlab = "Beta Coefficient", ylab = "Density")
# abline(v = S[2], col = "blue", lwd = 2)
# #############################################################X
# # Random effect on WMU
# ###----------------------------------------------------#X
# par(mfrow = c(3,2))
# for(i in 1:length(raneff)){
#   samp <- paste0("female.wmu.effect[", i, "]")
#   matplot(c7_results$samples$chain1[,samp],type="l", main = paste0("Female random effect on wmu"))
#   abline(h = raneff[i], col = "blue", lwd = 2)
#   plot(density(c7_results$samples$chain1[,samp]), main = "Density",
#        xlab = "Beta Coefficient", ylab = "Density")
#   abline(v = raneff[i], col = "blue", lwd = 2)
# }
# 
