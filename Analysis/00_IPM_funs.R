###############################################################################X
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#################### Integrated Population Model (IPM): ########################
#                 #---# PhD Dissertation: Chapter 1 #---#
#        Creating a Bayesian IPM to inform turkey management in PA
###                                                                         ###X
# Script to create data preparation functions for both IPM models
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #X
# 
# Created by: Veronica A. Winter
# Created: January 2024
# Last edited: XX/XX/XXXX
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #X
###############################################################################X
# Load in libraries
library(jagsUI)
library(tidyr)
library(MCMCvis)
library(stringr)
library(dplyr)
library(ggplot2)
library(RODBC)
library(nimble)
library(Rlab)

# Functions
# From Key and Schaub
# Define function to create a matrix with information about known latent state z
known.state.mr <- function(mr){
  state <- matrix(NA, nrow = dim(mr)[1], ncol = dim(mr)[2])
  rec <- which(rowSums(mr)==2)
  for (i in 1:length(rec)){
    n1 <- min(which(mr[rec[i],]==1))
    n2 <- max(which(mr[rec[i],]==1))
    state[rec[i],n1:n2] <- 1
    state[rec[i],n1] <- NA
    state[rec[i],n2:dim(mr)[2]] <- 0
  }
  return(state)
}

# Define function to create a matrix of initial values for latent state z
mr.init.z <- function(mr){
  ch <- matrix(NA, nrow = dim(mr)[1], ncol = dim(mr)[2])
  rec <- which(rowSums(mr)==1)
  for (i in 1:length(rec)){
    n1 <- which(mr[rec[i],]==1)
    ch[rec[i],n1:dim(mr)[2]] <- 0
    ch[rec[i],n1] <- NA
  }
  return(ch)
}

# i. Create vector with occasion of marking 
get.first <- function(x) min(which(x!=0))

# cloglog fun
cloglog <- function(x, y){
  1-exp(-exp(t(x)))
}


hen_histories <- function(df, filter_sex = NA) {
  library(dplyr)
  library(lubridate)
  
  # Filter by sex if specified
  if (!is.na(filter_sex)) {
    df <- df %>% filter(sex == filter_sex)
  }
  
  # Parse dates and calculate fate and end date
  df <- df %>%
    mutate(
      begin = as.Date(paste0(captyr, "-", captmo, "-", captday)),
      censor = as.Date(paste0(cenyr, "-", cenmo, "-", cenday), format = "%Y-%m-%d"),
      mort = as.Date(paste0(estmortyr, "-", estmortmo, "-", estmortday), format = "%Y-%m-%d"),
      fate = case_when(
        !is.na(mort) ~ "Dead",
        !is.na(censor) ~ "Censor",
        TRUE ~ "Alive"
      ),
      end = case_when(
        fate == "Alive" ~ Sys.Date(),
        fate == "Censor" ~ censor,
        fate == "Dead" ~ mort
      )
    )
  
  # Transition juveniles to adults after year 1 (transition in June)
  df <- df %>%
    mutate(age = case_when(
      age == "J" & !is.na(end) & month(end) %in% 5:12 & year(end) > captyr ~ "A",
      age == "A" ~ "A",
      TRUE ~ "J"
    ))
  
  # Define study period
  start_year <- year(min(df$begin, na.rm = TRUE))
  end_year <- year(max(df$end, na.rm = TRUE))
  telem.nyears <- end_year - start_year + 1
  
  # Calculate start and stop months for each individual
  df <- df %>%
    mutate(
      start = month(begin) + 12 * (year(begin) - start_year),
      stop = month(end) + 12 * (year(end) - start_year),
      first = month(begin),
      last = month(end),
      year_start = year(begin) - start_year + 1,  # Relative start year
      year_end = year(end) - start_year + 1      # Relative end year
    )
  
  # Create is.juvenile variable
  df <- df %>%
    mutate(is.juvenile = ifelse(age == "J", 1, 0))
  
  # Add WMU indicator variables
  df <- df %>%
    mutate(
      is.2d = ifelse(studyarea == "2D", 1, 0),
      is.3d = ifelse(studyarea == "3D", 1, 0),
      is.5c = ifelse(studyarea == "5C", 1, 0),
      is.4d = ifelse(studyarea == "4D", 1, 0)  # Intercept WMU
    )
  
  # Initialize the telem_month array with dimensions: [individuals, years, months]
  telem.nind <- nrow(df)
  telem_month <- array(NA, dim = c(telem.nind, telem.nyears, 12))  # [individuals, years, months]
  
  # Fill the telem_month array: For each individual, year, and month
  for (i in 1:telem.nind) {
    for (y in 1:telem.nyears) {
      for (m in 1:12) {
        # Calculate the start and stop months relative to each year
        month_index <- (y - 1) * 12 + m
        # Check if the individual was alive in this month
        if (month_index >= df$start[i] & month_index <= df$stop[i]) {
          telem_month[i, y, m] <- m  # Store month number (1-12)
        }
      }
    }
  }
  
  # Initialize 3D arrays for juvenile status and observation data
  is_juvenile_matrix <- array(0, dim = c(telem.nind, telem.nyears, 12))
  is_data <- array(0, dim = c(telem.nind, telem.nyears, 12))
  
  # Fill the 3D arrays (is_juvenile_matrix and is_data)
  for (i in 1:telem.nind) {
    for (y in 1:telem.nyears) {
      for (m in 1:12) {
        # Calculate the start and stop months relative to each year
        month_index <- (y - 1) * 12 + m
        if (month_index >= df$start[i] && month_index <= df$stop[i]) {
          is_juvenile_matrix[i, y, m] <- df$is.juvenile[i]
          is_data[i, y, m] <- 1  # Mark the data as observed
        }
      }
    }
  }
  
  # Return the processed dataframe, is_juvenile_matrix, is_data, and telem_month
  return(list(
    hen = df,
    is_juvenile_matrix = is_juvenile_matrix,
    is_data = is_data,
    telem_month = telem_month
  ))
}


encounter_histories <- function(df, filter_sex = NA, start_year = NA) {
  library(dplyr)
  library(lubridate)
  
  # Filter by sex if specified
  if (!is.na(filter_sex)) {
    df <- df %>% filter(sex == filter_sex)
  }
  
  # Ensure required date columns are in the correct format
  # Calculate capture, censoring/mortality dates, and fate
  df$begin <- as.Date(paste0(df$captyr, "-", df$captmo, "-", df$captday))
  df$censor <- as.Date(paste0(df$cenyr, "-", df$cenmo, "-", df$cenday), format = "%Y-%m-%d")
  df$mort <- as.Date(paste0(df$estmortyr, "-", df$estmortmo, "-", df$estmortday), format = "%Y-%m-%d")
  df$fate <- ifelse(is.na(df$mort) & is.na(df$censor), "Alive", ifelse(!is.na(df$censor), "Censor", "Dead"))
  df$end <- case_when(df$fate == "Alive" ~ Sys.Date(), df$fate == "Censor" ~ df$censor, TRUE ~ df$mort)
  
  # Transition juveniles to adults after year 1
  # Make juveniles move to adult in June
  df <- df %>%
    mutate(age = case_when(
      age == "J" & !is.na(end) & month(end) %in% 5:12 & year(end) > captyr ~ "A",
      age == "A" ~ "A",
      TRUE ~ "J"
    )) 
  
  # Calculate study period
  if (is.na(start_year)) {
    start_year <- year(min(df$begin, na.rm = TRUE))
  }
  end_year <- year(max(df$end, na.rm = TRUE))
  
  # Calculate start and stop months for each animal
  df <- df %>%
    mutate(
      start = month(begin) + 12 * (year(begin) - start_year),
      stop = month(end) + 12 * (year(end) - start_year)
    )
  
  # Initialize 3D array for statuses [i, y, m]
  n_ind <- nrow(df)
  n_years <- end_year - start_year + 1
  STATUS <- array(NA, dim = c(n_ind, n_years, 12))
  
  # Fill the STATUS array
  for (i in 1:n_ind) {
    animal_dead <- FALSE
    for (y in 1:n_years) { # Year index
      for (m in 1:12) { # Month index
        month_index <- (y - 1) * 12 + m
        if (month_index >= df$start[i] && month_index <= df$stop[i] && !animal_dead) {
          STATUS[i, y, m] <- 1 # Alive
        }
        if (df$fate[i] == "Dead" && month_index == df$stop[i]) {
          STATUS[i, y, m] <- 0 # Mark death
          animal_dead <- TRUE
        }
        if (animal_dead && month_index > df$stop[i]) {
          STATUS[i, y, m] <- NA # After death
        }
      }
      if (animal_dead) break # Stop if animal is dead
    }
  }
  
  # Return the 3D STATUS array
  return(STATUS)
}


# Creating poult counts function for recruitment model ----
# Data organization from access database
poult_prep <- function(dbrood) {
  # Filter for specific count type
  broods <- dbrood %>% 
    filter(counttype == "4week")
  
  # Calculate number of unique WMUs and occasions
  n.wmu <- length(unique(broods$wmu_bc))
  n.occasions <- length(unique(broods$countyr))
  
  # Group by WMU and count year, then summarize
  poults <- dbrood %>% 
    dplyr::filter(counttype == "4week") %>% 
    rename(wmu = wmu_bc) %>% 
    group_by(wmu, countyr) %>% 
    dplyr::select(taghen_totalpoults) %>% 
    summarise(n = sum(taghen_totalpoults, na.rm = T), .groups = "drop") %>% 
    filter(!n == 0)
  
  # Aggregate by WMU, summing the total poults
  hp <- data.frame(
    model.matrix(countyr ~ wmu + n, data = poults)
  )
  
  p <- poults %>% 
    as.data.frame() %>% 
    mutate(wmu = case_when(wmu == "2D" ~ 1,
                           wmu == "3D" ~ 2,
                           wmu == "4d" ~ 3,
                           TRUE ~ 4)) %>% 
    relocate(wmu, .before = "countyr")
  
  # Create a two-dimensional array 'hp_array' based on unique 'wmu' and 'countyr' combinations
  unique_combinations <- unique(p[, c("wmu", "countyr")])
  hp_array <- matrix(0, nrow = nrow(unique_combinations), ncol = n.occasions)
  
  for (i in 1:nrow(unique_combinations)) {
    row_indices <- which(p$wmu == unique_combinations[i, "wmu"] & p$countyr == unique_combinations[i, "countyr"])
    hp_array[i, ] <- poults$n[row_indices]
  }
  
  
  # Return a list containing both data frames
  return(list(n_wmu = n.wmu, n_occasions = n.occasions, hp_array = hp_array, 
              hp = hp, p = p))
}

# Simulating recovery data
simul.mr <- function(S, R, marked){
  n.occasions <- dim(S)[2]
  MR <- matrix(NA, ncol = n.occasions+1, nrow = sum(marked))
  # Define vector with occasion of markings
  mark.occ <- rep(1:n.occasions, marked)
  # Fill CH matrix
  for(i in 1:sum(marked)){
    MR[i, mark.occ[i]] <- 1
    for(t in mark.occ[i]:n.occasions){
      # Bernoulli trial: has indiv. survived?
      sur <- rbinom(1, 1, S[i,t])
      if(sur == 1) next # if still alive, move on
      # Bernoulli trial: has it died and been recovered?
      rp <- rbinom(1, 1, R[i,t])
      if(rp==0){
        MR[i, t+1] <- 0
        break
      }
      if(rp ==1){
        MR[i, t+1] <- 1
        break
      }
    } # t
  } # i
  
  # Replace the NA in the file by 0
  MR[which(is.na(MR))] <- 0
  return(MR)
}

get.first <- function(x){
  min(which(x!=0))
} 


DRM_simulate_data<-function(n.occasions, marked, S, rr, I, raneff) {
  
  MR <- matrix(NA, ncol = n.occasions+1, nrow = marked*n.occasions)
  # Define vector with occasion of markings
  mark.occ <- rep(1:n.occasions, each=marked)
  wmu <- rep(c(1:10),marked*n.occasions/10)
  
  # Simulate data for each individual
  for (i in 1:(marked*n.occasions)) {
    MR[i, mark.occ[i]] <- 1
    for (t in mark.occ[i]:n.occasions) {
      # Calculate survival probability
      s <- rbinom(1, 1, (S[2-I[i,t]]+raneff[wmu[i]]) )
      if(s == 1) {next} # if still alive, move on
      # Bernoulli trial: has it died and been recovered?
      r <- rr[1] * I[i, t] + rr[2] * (1 - I[i, t])
      recovery_data <- rbinom(1, 1, r)
      if(recovery_data==0){
        MR[i, t+1] <- 0
        break
      }
      if(recovery_data ==1){
        MR[i, t+1] <- 1
        break
      }
    } # t
  } # i
  
  # Replace the NA in the file by 0
  MR[which(is.na(MR))] <- 0
  return(MR)
}

########################################################X
# Function to get the first non-zero value in a row
########################################################X
get_first_non_zero <- function(row) {
  return(row[row != 0][1])
}