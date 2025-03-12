###############################################################################X
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#################### Integrated Population Model (IPM): #######################X
#                 #---# PhD Dissertation: Chapter 1 #---#
#        Creating a Bayesian IPM to inform turkey management in PA
###                *** SIMULATION TEST: Females ***                         ###X
###                                                                         ###X
#    Modeling: Creating a  dead-recovery model for males and females by 
# wildlife management unit (WMU), and age class over seasons/years
###                                                                         ###X
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #X
# 
# Created by: Veronica A. Winter
# Last edited: 05/09/2024
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #X
###############################################################################X

# clean env
rm(list = ls())
gc()

# Source scripts of functions and data preparation
source("Analysis/Scripts/00_IPM_funs.R")
set.seed(0123)

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
marked <- rep(400, n.occasions)  # Number of individuals
nind <- sum(marked)  # Total number of individuals
alpha <- 0.5
beta <- rnorm(n.occasions, 0 , 0.5)
sigma <- runif(1, 0, 1)
gamma <- rnorm(n.wmu, 0, sigma)
rr <- c(0.5, 0.5)
rrate.j <- 0.7
rrate.a <- 0.3

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
MR.female <- DRM_simulate_data(n.occasions, marked, alpha, beta, gamma, rr, rrate.j, rrate.a, I, II, J, wmu)
head(MR.female)
tail(MR.female)
table(MR.female)

# Step 2: Parameters fro DRM
# Initialize first  matrices 
f.female <- apply(MR.female, 1, get.first)
# Create a matrix of initial values for latent state z
z.female <- mr.init.z(MR.female)

# Rename variables for to distinguish between sexes in model
I.female  <- I 
II.female  <- II
time.female <- J
beta.time.female <- beta
gamma.wmu.female <- gamma
random.effect.wmu.female <- sigma
alpha.female <- alpha
rrate.hen <- rrate.a
rrate.jenny <- rrate.j
rr.female <- rr
wmu.female <- wmu

#############################################################################X
# Call Nimble model ----
#############################################################################X
source("models/drm_female.R")
#############################################################X
##           Model run set up ----
#############################################################X
# Data
test_dat <- list(
  y.female = MR.female, 
  I.female = I.female, 
  II.female = II.female, 
  time.female = time.female, 
  z.female = known.state.mr(MR.female)
  )

# Constants
consts <- list(
  # f=f,
  # nind = dim(MR)[1], 
  f.female=f.female,
  nind.female = dim(MR.female)[1], 
  n.occasions.female = dim(MR.female)[2],
  n.wmu = n.wmu, 
  wmu.female = wmu.female
)

inits <- list(
  beta.time.female = beta.time.female,
  gamma.wmu.female = gamma.wmu.female,
  alpha.female = alpha.female,
  rrate.hen = rrate.hen,
  rrate.jenny = rrate.jenny,
  rr.female = rr.female,
  random.effect.wmu.female = random.effect.wmu.female
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
  monitors = c("rr.female",
               "beta.time.female",
               "gamma.wmu.female",
               "alpha.female",
               "rrate.jenny",
               "rr.female"),
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
for(i in 1:length(beta.time.female)){
  samp <- paste0("beta.time.female[", i, "]")
  matplot(drm_results$samples[,samp],type="l", main = paste0("beta.time.female[", i, "]"))
  abline(h = beta.time.female[i], col = "red", lwd = 2)
  plot(density(drm_results$sample[,samp]), main = "Density",
       xlab = "Beta Coefficient", ylab = "Density")
  abline(v = beta.time.female[i], col = "red", lwd = 2)
  
}

# Gamma
par(mfrow = c(3,2))
for(i in 1:length(gamma.wmu.female)){
  samp <- paste0("gamma.wmu.female[", i, "]")
  matplot(drm_results$samples[,samp],type="l", main = paste0("gamma.wmu.female[", i, "]"))
  abline(h = gamma.wmu.female[i], col = "purple", lwd = 2)
  plot(density(drm_results$sample[,samp]), main = "Density",
       xlab = "Beta Coefficient", ylab = "Density")
  abline(v = gamma.wmu.female[i], col = "purple", lwd = 2)
  
}

# Alpha
par(mfrow=c(3,2))
matplot(drm_results$samples[,"alpha.female"],type="l", main = "Alpha (effect of Juvenile)")
abline(h = alpha.female, col = "blue", lwd = 2)
plot(density(drm_results$sample[,"alpha.female"]), main = "Density",
     xlab = "Beta Coefficient", ylab = "Density")
abline(v = alpha.female, col = "blue", lwd = 2)

# rrate.jenny
matplot(drm_results$samples[,"rrate.jenny"],type="l", main = "rrate.jenny")
abline(h = rrate.jenny, col = "red", lwd = 2)
plot(density(drm_results$sample[,"rrate.jenny"]), main = "Density",
     xlab = "Beta Coefficient", ylab = "Density")
abline(v = rrate.jenny, col = "red", lwd = 2)

# rr.female
matplot(drm_results$samples[,"rr.female[1]"],type="l", main = "rr.female[1]")
abline(h = rr.female[1], col = "darkorange", lwd = 2)
plot(density(drm_results$sample[,"rr.female[1]"]), main = "Density",
     xlab = "Beta Coefficient", ylab = "Density")
abline(v = rr.female[1], col = "darkorange", lwd = 2)

### Find the effective sample size ----
## If it is too low, increase niter above and re-run MCMC
coda::effectiveSize(drm_results$samples)
summary(drm_results$samples)
