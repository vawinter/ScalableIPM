#----------------------------------------------------------#
# Simulating HWB data in base R
# VAW
# 1/6/2025
#-----------
# Set random seed for reproducibility
set.seed(1235)
#----------------------------------------------------------#
# Function to simulate MR (mark-recovery) data
#----------------------------------------------------------#
DRM_simulate_data <- function(marked, n.occasions, S, rr, rrate.j, rrate.a, raneff, I, II) {
  # Initialize matrix to store mark-recovery data
  MR <- matrix(NA, ncol = n.occasions + 1, nrow = marked * n.occasions)
  
  # Define marking occasions and WMU assignments
  mark.occ <- rep(1:n.occasions, each = marked)
  wmu <- rep(1:10, length.out = marked * n.occasions)
  
  # Simulate data for each individual
  for (i in 1:(marked * n.occasions)) {
    MR[i, mark.occ[i]] <- 1  # Mark individual at the initial occasion
    for (t in mark.occ[i]:n.occasions) {
      # Simulate survival probability
      s <- rbinom(1, 1, (S[2 - I[i, t]] + raneff[wmu[i]]))
      if (s == 1) next  # Individual survives to the next occasion
      
      # Simulate recovery probability
      r <- rr[1] * I[i, t] * rrate.j * II[i] + 
        rr[1] * I[i, t] * (1 - II[i]) +
        rr[2] * (1 - I[i, t]) * rrate.a * II[i] + 
        rr[2] * (1 - I[i, t]) * (1 - II[i])
      recovery_data <- rbinom(1, 1, r)
      
      # Update MR matrix based on recovery
      if (recovery_data == 0) {
        MR[i, t + 1] <- 0  # Not recovered
        break
      }
      if (recovery_data == 1) {
        MR[i, t + 1] <- 1  # Recovered
        break
      }
    }
  }
  
  # Replace NA values with 0 (no recovery)
  MR[is.na(MR)] <- 0
  return(MR)
}

#----------------------------------------------------------#
# Helper functions for state and design matrix generation
#----------------------------------------------------------#

# Function to find the first non-zero occurrence in a vector
get.first <- function(x) min(which(x != 0))

# Function to calculate known latent states
known.state.mr <- function(mr) {
  state <- matrix(NA, nrow = dim(mr)[1], ncol = dim(mr)[2])
  rec <- which(rowSums(mr) == 2)  # Rows with exactly two recoveries
  for (i in 1:length(rec)) {
    n1 <- min(which(mr[rec[i], ] == 1))
    n2 <- max(which(mr[rec[i], ] == 1))
    state[rec[i], n1:n2] <- 1
    state[rec[i], n1] <- NA
    state[rec[i], n2:dim(mr)[2]] <- 0
  }
  return(state)
}

# Function to create a design matrix for time coefficients
design_matrix <- function(n) {
  j1 <- rep(1, n - 1)
  j2 <- diag(n - 2)
  j3 <- rep(-1, n - 2)
  rbind(j1, cbind(j2, j3))
}

#----------------------------------------------------------#
# Simulation parameters
#----------------------------------------------------------#
rrate.j <- 0.7          # Reporting rate for juveniles
rrate.a <- 0.9          # Reporting rate for adults
n.occasions <- 4        # Number of occasions
S <- c(0.7, 0.5)        # Survival probabilities (juvenile, adult)
rr <- c(0.3, 0.5)       # Seber recovery rates (juvenile, adult)
marked <- 200           # Number of marked individuals
raneff <- c(-0.043, -0.006, -0.038, 0.056, -0.081, -0.094, -0.006, 0.003, 0.084, -0.057)

#----------------------------------------------------------#
# Simulate data for adults and juveniles
#----------------------------------------------------------#

# Adults: Indicator for juveniles (I) and reward (II)
I.ad <- matrix(0, nrow = marked * n.occasions, ncol = n.occasions)
II.ad <- sample(0:1, marked * n.occasions, replace = TRUE)
MR.ad <- DRM_simulate_data(marked, n.occasions, S, rr, rrate.j, rrate.a, raneff, I.ad, II.ad)

# Juveniles: Indicator for juveniles (I) and reward (II)
I.jv <- matrix(c(rep(c(1, 0, 0, 0), each = marked), 
                 rep(c(0, 1, 0, 0), each = marked), 
                 rep(c(0, 0, 1, 0), each = marked), 
                 rep(c(0, 0, 0, 1), each = marked)),
               nrow = marked * n.occasions, ncol = n.occasions)
II.jv <- sample(0:1, marked * n.occasions, replace = TRUE)
MR.jv <- DRM_simulate_data(marked, n.occasions, S, rr, rrate.j, rrate.a, raneff, I.jv, II.jv)

#----------------------------------------------------------#
# Combine data for adults and juveniles
#----------------------------------------------------------#
MR <- rbind(MR.ad, MR.jv)  # Combined mark-recovery matrix
I <- rbind(I.ad, I.jv)     # Combined juvenile indicators
II <- c(II.ad, II.jv)      # Combined reward indicators
wmu <- c(rep(1:10, length.out = nrow(MR)))

#----------------------------------------------------------#
# Generate additional variables for modeling
#----------------------------------------------------------#
f <- apply(MR, 1, get.first)         # First marking occasion
known.state <- known.state.mr(MR)    # Known latent states
J <- design_matrix(n.occasions + 1)  # Design matrix for time coefficients

#----------------------------------------------------------#
# Summary of simulated data
#----------------------------------------------------------#
cat("Summary of Simulated Data\n")
print(table(apply(MR, 1, sum)))  # Number of recoveries per individual
print(table(wmu))                # Distribution of WMU assignments

#----------------------------------------------------------#
# Prepare data for downstream modeling
#----------------------------------------------------------#
male.y <- MR                  # Mark-recovery matrix
male.I <- I                   # Juvenile indicators
male.II <- II                 # Reward indicators
male.wmu <- wmu               # WMU assignments
male.f <- f                   # First marking occasion
male.time.param <- J          # Design matrix
male.z <- known.state         # Known latent states
male.n.occasions <- ncol(MR)  # Number of occasions
male.nind <- nrow(MR)         # Number of individuals

