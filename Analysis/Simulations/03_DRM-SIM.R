###############################################################################X
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#################### Integrated Population Model (IPM): #######################X
#                 #---# PhD Dissertation: Chapter 1 #---#
#              a Bayesian IPM to inform turkey management in PA             ###X
###                *** SIMULATION TEST (Nimble, MF) ***                     ###X
###                                                                         ###X
#    Modeling: Creating a  dead-recovery model for males and females by 
# wildlife management unit (WMU), and age class over seasons/years
###                                                                         ###X
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #X
#  Code based off DRD simulations
# Created by: Veronica A. Winter
# Last edited: 05/28/2024
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #X
###############################################################################X
# Clean env
rm(list = ls())
gc()

# Set options
options(max.print=1000000)

# Call functions
source("Analysis/OneDrive_2024-05-23/Simulation dead recovery/00_funs.R")

# Libraries
library(nimble)
library(MCMCvis)

##################################################
# Set variables ----
##################################################X
# Note: these starting values recover better for F but not for M,
# I will leave the 0.7 and 0.9 now for consistency
# rrate.j <- 0.4 # reporting rate juveniles
# rrate.a <- 0.5 # reporting rate adults

rrate.j <- 0.7 # reporting rate juveniles
rrate.a <- 0.9 # reporting rate adults
n.occasions <- 4 # number of occasions
S <- c(0.7,0.5)  # survival c(juv, ad)
rr <- c(0.3,0.5)  # seber recovery rate c(juv, ad)
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
simulate_data <- function(seed, marked, n.occasions, S, rr, rrate.j, rrate.a, raneff) {
  set.seed(seed)
  
  ## Sim data for adults ----
  # indicator for age (juv = 1)
  I.ad <- matrix(0,nrow=marked*n.occasions,ncol=n.occasions)
  # indicator for reward (reward  = 0)
  II.ad <- sample(0:1, marked*n.occasions, replace = TRUE)
  # Create MR matrix for adults
  MR.ad <- DRM_simulate_data(n.occasions, marked, S, rr, I.ad, raneff, II.ad)
  ###-----------------------------------------------------#X
  ## Sim data for juveniles ----
  # indicator for age (juv = 1)
  I.jv <- matrix(c(rep(c(1,0,0,0), each=marked), rep(c(0,1,0,0), each=marked),
                   rep(c(0,0,1,0), each=marked), rep(c(0,0,0,1), each=marked)),
                 nrow=marked*n.occasions,ncol=n.occasions)
  # indicator for reward (reward  = 0)
  II.jv <- sample(0:1, marked*n.occasions, replace = TRUE)
  # Create MR matrix for Juveniles
  MR.jv <- DRM_simulate_data(n.occasions, marked, S, rr, I.jv, raneff, II.jv)
  ###-----------------------------------------------------#X
  # Merge separate matrices ---- 
  ### Adults and Juveniles, 100% reporting, 10 mgt units     
  MR <- rbind(MR.ad, MR.jv)
  I <- rbind(I.ad,I.jv)
  II <- c(II.ad,II.jv)
  wmu <- c(rep(1:10, marked * n.occasions / 10), rep(1:10, marked * n.occasions / 10))
  ###-----------------------------------------------------#X
  # Create other necessary variables for the model
  f <- apply(MR, 1, get.first)
  # Create design matrix for time coefficient using a mean parameterization
  j1 <- rep(1,dim(MR)[2]-1)
  j2 <- diag(dim(MR)[2]-2)
  j3 <- rep(-1,dim(MR)[2]-2)
  j23 <- rbind(j2,j3)
  J <- as.matrix(cbind(j1,j23)) ; colnames(J) <- NULL ; rownames(J) <- NULL
  ###-----------------------------------------------------#X
  known.state <- known.state.mr(MR)
  is.na(known.state) <- NA
  ###-----------------------------------------------------#X
  # Debugging: Print dimensions
  cat("Dimensions of MR:", dim(MR), "\n")
  cat("Dimensions of I:", dim(I), "\n")
  cat("Dimensions of II:", length(II), "\n")
  cat("Dimensions of wmu:", length(wmu), "\n")
  cat("Dimensions of J:", dim(J), "\n")
  cat("Dimensions of occasions:", (dim(MR)[2]), "\n")
  ###-----------------------------------------------------#X
  list(MR = MR, I = I, II = II, wmu = wmu, f = f, J = J, known.state = known.state,
       rrate.j = rrate.j, rrate.a = rrate.a, n.occasions = n.occasions)
}
##################################################X
# Simulate data separately for the sexes ----
##################################################X
# Males w. set.seed
male_data <- simulate_data(123, marked, n.occasions, S, rr, rrate.j, rrate.a, raneff)

# Females w. set.seed
female_data <- simulate_data(321, marked, n.occasions, S, rr, rrate.j, rrate.a, raneff)

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
female.f <- female_data$f
female.time.param <- female_data$J
female.z <- female_data$known.state
female.rrate.j <- female_data$rrate.j
female.rrate.a <- female_data$rrate.a
female.n.occasions <- (dim(female_data$MR)[2])
female.nind <- dim(female_data$MR)[1]
###-----------------------------------------------------#X
# Create n.occasions-1 variable
true.occasions = n.occasions

##################################################X
# Estimate parameters in Nimble ----
##################################################X
# source model
source("Analysis/Scripts/NIMBLE_DRM_MF.R")
###-----------------------------------------------------#X
### Nimble model set up ----
# Data
nimble.data <- list(
  # Male data
  male.y = male.y, 
  male.time.param = male.time.param,
  ###-----------#X
  # Female data
  female.y = female.y,
  female.time.param = female.time.param
  )

# Constants
consts <- list(
  # Male constants
  male.f = male.f,
  male.I = male.I,
  male.II = male.II,
  male.nind = male.nind,
  male.n.occasions = male.n.occasions,
  male.n.wmu = 10,
  male.wmu = male.wmu,
  ###-----------#X
  # Female constants
  female.f = female.f,
  female.I = female.I,
  female.II = female.II,
  female.nind = female.nind,
  female.n.occasions = female.n.occasions,
  female.n.wmu = 10,
  female.wmu = female.wmu
)

# Initial values
inits <- list(
  # Male inits
  male.z = male.z,
  male.juvenile.effect=rnorm(1,0,0.5),
  male.time.effect = rnorm((n.occasions), 0, 0.5), 
  male.seber.recov=runif(2,0,1),
  male.rrate.a = rnorm(1,0.87,0.039), 
  male.rrate.j = rnorm(1,0.71,0.072), 
  male.sigma=runif(1,0,10), 
  ###-----------#X
  # Female inits
  female.z = female.z,
  female.juvenile.effect=rnorm(1,0,0.5),
  female.time.effect=rnorm((true.occasions),0,0.5),
  female.seber.recov=runif(2,0,1),
  female.rrate.a=runif(1,0,1),
  female.rrate.j=runif(1,0,1), 
  female.sigma=runif(1,0,10) 
)

##################################################X
## Call Nimble model ----
##################################################X
# Create the Nimble model
model <- nimbleModel(
  code = drm,
  data = nimble.data,
  constants = consts,
  inits = inits
)

# Check initialization issues
model$initializeInfo()
gc()

# Set MCMC specifications
ni <- 50000
nt <- 1
nb <- 7500
nc <- 2 # minimum 2 chains for rhat values

# Run MCMC 
c5_drm_results <- nimbleMCMC(
  model = model,
  monitors = c("male.juvenile.effect","male.time.effect","male.wmu.effect",
               "male.rrate.a","male.rrate.j", "male.mean.s.ad","male.mean.s.jv",
               "male.mean.harv.ad", "male.mean.harv.jv","male.seber.recov",
               "male.h.juv.wmu","male.h.ad.wmu", "male.s.juv.wmu","male.s.ad.wmu",
               ###------------------------------------------------------------#X
               "female.juvenile.effect","female.time.effect","female.wmu.effect",
               "female.rrate.a","female.rrate.j", "female.mean.s.ad",
               "female.mean.s.jv","female.mean.harv.ad","female.mean.harv.jv",
               "female.seber.recov","female.h.juv.wmu","female.h.ad.wmu",
               "female.s.juv.wmu","female.s.ad.wmu"),
  niter = ni, 
  nburnin = nb, 
  thin = nt, 
  nchains = nc, 
  setSeed = FALSE, 
  samplesAsCodaMCMC = TRUE,
  WAIC = TRUE)
 beepr::beep(1)


#############################################################X
#           Model diagnostics ----
#############################################################X
summary(c5_drm_results$samples)
MCMCvis::MCMCsummary(c5_drm_results$samples) # to look at Rhats
c5_drm_results$WAIC
MCMCvis::MCMCtrace(c5_drm_results$samples, iter = ni)

#############################################################X
##         Look at specific estimates
#         Males = red, Females = blue
#############################################################X
# Male (red)
#############################################################X
# recovery
par(mfrow = c(3,2))
for(i in 1:length(rr)){
  samp <- paste0("male.seber.recov[", i, "]")
  matplot(c5_drm_results$samples$chain1[,samp],type="l", main = paste0("male.seber.recov"))
  abline(h = rr[i], col = "red", lwd = 2)
  plot(density(c5_drm_results$samples$chain1[,samp]), main = "Density",
       xlab = "Beta Coefficient", ylab = "Density")
  abline(v = rr[i], col = "red", lwd = 2)
}
###----------------------------------------------------#X
par(mfrow=c(1,2))
# Survival c(juv, ad)
matplot(c5_drm_results$samples$chain1[,"male.mean.s.jv[1]"],type="l", main = "Male mean survival (J)")
abline(h = S[1], col = "red", lwd = 2)
plot(density(c5_drm_results$samples$chain1[,"male.mean.s.jv[1]"]), main = "Density",
     xlab = "Beta Coefficient", ylab = "Density")
abline(v = S[1], col = "red", lwd = 2)
###----------------------------------------------------#X
matplot(c5_drm_results$samples$chain1[,"male.mean.s.ad[2]"],type="l", main = "Male mean survival (A)")
abline(h = S[2], col = "red", lwd = 2)
plot(density(c5_drm_results$samples$chain1[,"male.mean.s.ad[2]"]), main = "Density",
     xlab = "Beta Coefficient", ylab = "Density")
abline(v = S[2], col = "red", lwd = 2)
#############################################################X
# Random effect on WMU
par(mfrow = c(3,2))
for(i in 1:length(raneff)){
  samp <- paste0("male.wmu.effect[", i, "]")
  matplot(c5_drm_results$samples$chain1[,samp],type="l", main = paste0("Male random effect on wmu"))
  abline(h = raneff[i], col = "red", lwd = 2)
  plot(density(c5_drm_results$samples$chain1[,samp]), main = "Density",
       xlab = "Beta Coefficient", ylab = "Density")
  abline(v = raneff[i], col = "red", lwd = 2)
}
#############################################################X
# Female (blue)
#############################################################X
# recovery
par(mfrow = c(3,2))
for(i in 1:length(rr)){
  samp <- paste0("female.seber.recov[", i, "]")
  matplot(c5_drm_results$samples$chain1[,samp],type="l", main = paste0("female.seber.recov"))
  abline(h = rr[i], col = "blue", lwd = 2)
  plot(density(c5_drm_results$samples$chain1[,samp]), main = "Density",
       xlab = "Beta Coefficient", ylab = "Density")
  abline(v = rr[i], col = "blue", lwd = 2)
}
#############################################################X
par(mfrow=c(1,2))
# Survival c(juv, ad)
###----------------------------------------------------#X
matplot(c5_drm_results$samples$chain1[,"female.mean.s.jv[1]"],type="l", main = "Female mean survival (J)")
abline(h = S[1], col = "blue", lwd = 2)
plot(density(c5_drm_results$samples$chain1[,"female.mean.s.jv[1]"]), main = "Density",
     xlab = "Beta Coefficient", ylab = "Density")
abline(v = S[1], col = "blue", lwd = 2)
###----------------------------------------------------#X
matplot(c5_drm_results$samples$chain1[,"female.mean.s.ad[1]"],type="l", main = "Female mean survival (A)")
abline(h = S[2], col = "blue", lwd = 2)
plot(density(c5_drm_results$samples$chain1[,"female.mean.s.ad[1]"]), main = "Density",
     xlab = "Beta Coefficient", ylab = "Density")
abline(v = S[2], col = "blue", lwd = 2)
#############################################################X
# Random effect on WMU
###----------------------------------------------------#X
par(mfrow = c(3,2))
for(i in 1:length(raneff)){
  samp <- paste0("female.wmu.effect[", i, "]")
  matplot(c5_drm_results$samples$chain1[,samp],type="l", main = paste0("Female random effect on wmu"))
  abline(h = raneff[i], col = "blue", lwd = 2)
  plot(density(c5_drm_results$samples$chain1[,samp]), main = "Density",
       xlab = "Beta Coefficient", ylab = "Density")
  abline(v = raneff[i], col = "blue", lwd = 2)
}
