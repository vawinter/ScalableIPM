###############################################################################X
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#################### Integrated Population Model (IPM): #######################X
#                 #---# PhD Dissertation: Chapter 1 #---#
#        Creating a Bayesian IPM to inform turkey management in PA
###                *** SIMULATION TEST: Males ***                           ###X
###                                                                         ###X
#    Modeling: Creating a  dead-recovery model for males and females by 
# wildlife management unit (WMU), and age class over seasons/years
###                                                                         ###X
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #X
# 
# Created by: Veronica A. Winter
# Last edited: 04/24/2024
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #X
###############################################################################X

# clean env
rm(list = ls())
gc()

# Source scripts of functions and data preparation
source("Analysis/Scripts/00_IPM_funs.R")

set.seed(0123)
#set.seed(2)
#############################################################################X
# Simulating data: initial values ----
#############################################################################X
## Note: The function for simulating data was constructed following Kery & Schaub. 
# I wanted to use the same function for both sexes, since these models should be the same, 
# and therefore create the variables first to run the function then rename after. 

# Step 1: Simulate Parameters for MR
## Parameters
n.wmu <- 10  # Number of wildlife management units, example value
n.occasions <- 10  # Number of months
marked <- rep(500, n.occasions)  # Number of individuals
nind <- sum(marked)  # Total number of individuals
alpha <- 0.5
beta <- rep(0.3, n.occasions)
sigma <- runif(1, 0, 10)
gamma <- rnorm(n.wmu, 0, sigma)
rr <- c(0.8, 0.6)
rrate.j <- 0.7
rrate.a <- 0.9

# Correct initialization of I, II, J, wmu
is_adult <- sample(0:1, nind, replace = TRUE)  # Randomly assign adult/juvenile status
I <- matrix(is_adult, nrow = nind, ncol = n.occasions)
II <- sample(0:1, nind, replace = TRUE)
wmu <- sort(sample(1:n.wmu, nind, replace = TRUE))


# Create design matrix for time coefficient (beta) using a mean parameterization
j1 <- rep(1,n.occasions)
j2 <- diag(n.occasions-1)
j3 <- rep(-1,n.occasions-1)
j23 <- rbind(j2,j3)
J <- as.matrix(cbind(j1,j23)) ; colnames(J) <- NULL ; rownames(J) <- NULL
#
#       [,1] [,2] [,3] [,4] [,5]
# [1,]    1    1    0    0    0
# [2,]    1    0    1    0    0
# [3,]    1    0    0    1    0
# [4,]    1    0    0    0    1
# [5,]    1   -1   -1   -1   -1
#############################################################################X
## Simulating data: simulating mark recovery and survival probabilities ----
#############################################################################X
# Execute function
MR.male <- DRM_simulate_data(n.occasions, marked, alpha, beta, gamma, rr, rrate.j, rrate.a, I, II, J, wmu)
head(MR.male)
tail(MR.male)
table(MR.male)

# Step 2: Parameters fro DRM
# Initialize first  matrices 
f.male <- apply(MR.male, 1, get.first)
# Create a matrix of initial values for latent state z
z.male <- mr.init.z(MR.male)

# Rename variables
I.male  <- I # indicator of juvenile (juvenile = 1)
II.male <- II # indicator of reward band (reward = 0)
time.male <- J # time estiamte w. means parameterization 
beta.time.male <- beta # beta for time 
gamma.wmu.male <- gamma # estimate for wmu
randome.effect.wmu.male <- sigma # wmu random effect estimate
alpha.male <- alpha # juvenile estimate
rrate.tom <- rrate.a # reporting rate for adults
rrate.jake <- rrate.j# reporting rate for juveniles
rr.male <- rr # prior on proportion mortality due to hunting - adults and juveniles
wmu.male <- wmu # keeping it clear which variables belong to which sex

#############################################################################X
# Call Nimble model ----
#############################################################################X
source("models/drm_male.R")
#############################################################X
##           Model run set up ----
#############################################################X
# Data
test_dat <- list(
  y.male = MR.male,
  I.male = I.male,
  II.male = II.male,
  time.male = time.male,
  z.male = known.state.mr(MR.male)
)

# Constants
consts <- list(
  f.male=f.male,
  nind.male = dim(MR.male)[1],
  n.occasions.male = dim(MR.male)[2],
  wmu.male = wmu.male,
  n.wmu = n.wmu 
)

inits <- list(
  beta.time.male = beta.time.male,
  gamma.wmu.male = gamma.wmu.male,
  alpha.male = alpha.male,
  rrate.tom = rrate.tom,
  rrate.jake = rrate.jake,
  rr.male = rr.male,
  randome.effect.wmu.male = randome.effect.wmu.male
)

# Create the Nimble model
model <- nimbleModel(
  code = drm,
  data = test_dat,
  constants = consts,
  inits = inits
)

# Check initialization issues
model$initializeInfo()
gc()

# Set MCMC specifications
burn <- 1000
iter <-  40000
t <- 1
chain <- 1

# Run MCMC 
drm_results <- nimbleMCMC(
  model = model,
  monitors = c("rr.male",
               "beta.time.male",
               "gamma.wmu.male",
               "randome.effect.wmu.male",
               "alpha.male",
               "rrate.jake"),
  niter = iter, 
  nburnin = burn, 
  thin = t, 
  nchains = chain, 
  setSeed = FALSE, 
  samplesAsCodaMCMC = TRUE,
  WAIC = TRUE)
# beepr::beep(1)

#############################################################X
##           Model diagnostics ----
#############################################################X
### Look at parameter recovery & check fit ----
# Beta
par(mfrow = c(3,2))
for(i in 1:length(beta.time.male)){
  samp <- paste0("beta.time.male[", i, "]")
 matplot(drm_results$samples[,samp],type="l", main = paste0("Beta[", i, "]"))
  abline(h = beta.time.male[i], col = "red", lwd = 2)
  plot(density(drm_results$sample[,samp]), main = "Density",
       xlab = "Beta Coefficient", ylab = "Density")
  abline(v = beta.time.male[i], col = "red", lwd = 2)
  
}

# Gamma
par(mfrow = c(3,2))
for(i in 1:length(gamma.wmu.male)){
  samp <- paste0("gamma.wmu.male[", i, "]")
  matplot(drm_results$samples[,samp],type="l", main = paste0("Gamma[", i, "]"))
  abline(h = gamma.wmu.male[i], col = "purple", lwd = 2)
  plot(density(drm_results$sample[,samp]), main = "Density",
       xlab = "Beta Coefficient", ylab = "Density")
  abline(v = gamma.wmu.male[i], col = "purple", lwd = 2)
  
}

# Alpha
par(mfrow=c(3,2))
matplot(drm_results$samples[,"alpha.male"],type="l", main = "Alpha (effect of Juvenile)")
abline(h = alpha.male, col = "blue", lwd = 2)
plot(density(drm_results$sample[,"alpha.male"]), main = "Density",
     xlab = "Beta Coefficient", ylab = "Density")
abline(v = alpha.male, col = "blue", lwd = 2)

# rrate.jake
matplot(drm_results$samples[,"rrate.jake"],type="l", main = "rrate.j")
abline(h = rrate.jake, col = "red", lwd = 2)
plot(density(drm_results$sample[,"rrate.jake"]), main = "Density",
     xlab = "Beta Coefficient", ylab = "Density")
abline(v = rrate.jake, col = "red", lwd = 2)

# rr.male
matplot(drm_results$samples[,"rr.male[1]"],type="l", main = "rr.male[1]")
abline(h = rr.male[1], col = "darkorange", lwd = 2)
plot(density(drm_results$sample[,"rr.male[1]"]), main = "Density",
     xlab = "Beta Coefficient", ylab = "Density")
abline(v = rr.male[1], col = "darkorange", lwd = 2)

### Find the effective sample size ----
## If it is too low, increase niter above and re-run MCMC
coda::effectiveSize(drm_results$samples)
summary(drm_results$samples)




