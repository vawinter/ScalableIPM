###############################################################################X
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#################### Integrated Population Model (IPM): ########################
#                 #---# PhD Dissertation: Complex IPM #---#
#        Creating a Bayesian IPM to inform turkey management in PA
###                                                                         ###X
# Dead-recovery model for turkey banding data analysis of male gobblers, PA
#
# Using this script to format female data same as males (if we had enough F data)
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #X
###############################################################################X

# Clean env
rm(list = ls())
gc()

# Libraries
library(RODBC)
library(stringr)
library(dplyr)
library(lubridate)

# Bring in functions
source("00_IPM_funs.R")

# Data cleaning ----
#### Read in data and select data used for analysis
# # Release data
df2a <- read.csv("Data/Banding_harv_data/Band Data Export Report_20240909.csv")  # release data
df2a <- rename(df2a,bandid="band_id")
df2a$bandid <- as.character(df2a$bandid)
df2b <- read.csv("Data/Banding_harv_data/Reported Band Data Export Report_20240909_ForPSU.csv")  # recovery data
df2b <- rename(df2b, bandid="BandNumber")
df2b$bandid <- as.character(df2b$bandid)
##-----------------------------------------------------#X
# Data cleaning ----
df <- left_join(df2a[,c(1,2,6,7,11,17,22,35)], df2b[,c(1,6,7,15)], join_by(bandid))
df <- df[df$recapture!="Yes" & df$Turkey_Age!="U",] # remove recaptures and unknown age
df <- df[df$Turkey_Sex!="F",]  # remove females
df <- df[!df$EncounterReason %in% c("Alive", "Found band only", "Predation/Scavenged",
                                    "Retrieved by pet", "Trap related", "Unknown", "Vehicle"),]

df$capyr <- year(as.Date(df$Date_Capture, "%m/%d/%Y"))
df <- df[df$capyr != "2024", ]

### Group WMUs
# Table that defines grouping of WMUs
wmu <- c("1A","1B","2A","2B","2C","2D","2E","2F","2G","2H","3A","3B","3C","3D","4A","4B","4C","4D","4E","5A","5B","5C","5D")
grp <- c(1,1,2,3,2,2,2,4,4,4,4,4,5,6,7,7,8,7,8,9,9,10,10)
LookUp <- as.data.frame(cbind(wmu,grp))
df <- left_join(df,LookUp, by="wmu")

# filter groups I want
# df <- df[df$grp %in% c(2, 6, 7), ]
df <- df[df$wmu %in% c("2D", "3D", "4D"), ]

tst <- df %>% 
  group_by(wmu, capyr) %>% 
  summarise(num_banded = n()) 

table(df$wmu, df$capyr, df$EncounterReason, df$Turkey_Sex)


###Tabulate releases by Year and age
release <- as.data.frame(table(df$capyr,df$Turkey_Age))
release

### Set Year 1 to 2020 
df$rYr <- df$Year-2019 #recovery year 99=not recovered
df$rYr <- ifelse(is.na(df$rYr),99,df$rYr)
df$cYr <- df$capyr-2019 #release year

### Tabulate recoveries by year of release/recovery and age
df %>% 
  group_by(rYr,cYr,Turkey_Age, wmu) %>%
  summarise(n=n()) %>%
  mutate(total=sum(n))

### Create mark-recapture matrix MR
n.occasions <- max(df$cYr) # number release occasions
n.release <- dim(df)[1]
MR <- matrix(NA, ncol = n.occasions + 1, nrow = n.release)
# vector of occasion of marking
mark.occ <- df$cYr
recov.occ <- df$rYr
# Fill the CH matrix
for (i in 1:n.release) {
  MR[i,mark.occ[i]] <- 1 # write 1 at release occasion
  for (t in mark.occ[i]:n.occasions){
    if (t==recov.occ[i]) {MR[i,t+1] <- 1} # write 1 at recovery occasion
  } #t loop
} #i loop
MR[which(is.na(MR))] <- 0

### 8.2.2. Adapted model with constant parameters to allow for age-specific and time-specific rates
# Create vector with occasion of marking 
f <- apply(MR, 1, get.first)
##-----------------------------------------------------#X
# Create a known-fate (f.kf) matrix with 0, 1
f_factor <- factor(mark.occ)
# Create the 0-1 matrix using model.matrix
f.kf <- model.matrix(~ f_factor - 1)
##-----------------------------------------------------#X
# Create indicator array of which occasion an animal is juvenile
I <- array(0,c(dim(MR)[1],n.occasions))
for(i in 1:dim(MR)[1]){ 
  if(df$Turkey_Age[i]=="J"){
    I[i,f[i]] <- 1
  }
}
##-----------------------------------------------------#X
# Create indicator vector of whether animal has a reward band: 1=non-reward band
II <- ifelse(df$reward_band=="N",1,0)
##-----------------------------------------------------#X
# Create design matrix for age coefficient (beta) using a mean parameterization
j1 <- rep(1,dim(MR)[2]-1)
j2 <- diag(dim(MR)[2]-2)
j3 <- rep(-1,dim(MR)[2]-2)
j23 <- rbind(j2,j3)
J <- as.matrix(cbind(j1,j23)) ; colnames(J) <- NULL ; rownames(J) <- NULL
##-----------------------------------------------------#X
# Bundle data
known.state <- known.state.mr(MR)
is.na(known.state) <- NA

##-----------------------------------------------------#X
# Apply function to get the first non-zero value in a row (Juvenile) 
# to each row of the matrix
female.age <- apply(I, 1, get_first_non_zero)
# NA = adult, other = Juvenile
# For indexing, I want 1 = adult, 2 = juvenile
female.age <- ifelse(is.na(female.age), 1, 2)

##-----------------------------------------------------#X
female.wmu = df$wmu
female.wmu.index <- case_when(
  female.wmu == 2  ~ 1,
  female.wmu == 6  ~ 2,
  female.wmu == 7  ~ 3,
  female.wmu == 10 ~ 4)

##-----------------------------------------------------#X
# Create a list to store the variables
female_data <- list(
  female.f = f,
  female.I = I,
  female.II = II,
  female.time.param = J,
  female.nind = dim(MR)[1],
  female.n.occasions = (dim(MR)[2]),
  female.z = known.state,
  female.y = MR,
  female.n.wmu = length(unique(df$wmu)),
  female.wmu = df$wmu,
  female.age = female.age,
  female.f.kf = f.kf,
  female.wmu.index = female.wmu.index
)
##-----------------------------------------------------#X
# Define a directory where you want to save the RDS files
output_dir <- "Data/Research_IPM_setup-data/"

# Loop through the list and save each variable as a separate RDS file
for (var_name in names(female_data)) {
  variable <- female_data[[var_name]]
  file_name <- paste(output_dir, var_name, ".rds", sep = "")
  saveRDS(variable, file_name)
}
##################################################X
# Great, now for males ------------
##################################################X
rm(list = ls())
gc()

# Libraries
library(jagsUI)
library(MCMCvis)
library(RODBC)
library(stringr)
library(dplyr)
library(lubridate)
##-----------------------------------------------------#X
# Bring in functions
source("Analysis/Scripts/00_IPM_funs.R")
##-----------------------------------------------------#X
#### Read in data and select data used for analysis
# # Release data
df2a <- read.csv("Data/Banding_harv_data/Band Data Export Report_20240909.csv")  # release data
df2a <- rename(df2a,bandid="band_id")
df2a$bandid <- as.character(df2a$bandid)
df2b <- read.csv("Data/Banding_harv_data/Reported Band Data Export Report_20240909_ForPSU.csv")  # recovery data
df2b <- rename(df2b, bandid="BandNumber")
df2b$bandid <- as.character(df2b$bandid)
##-----------------------------------------------------#X
# Data cleaning ----
df <- left_join(df2a[,c(1,2,6,7,11,17,22,35)], df2b[,c(1,6,7,15)], join_by(bandid))
df <- df[df$recapture!="Yes" & df$Turkey_Age!="U",] # remove recaptures and unknown age
df <- df[df$Turkey_Sex!="F",]  # remove females
df <- df[!df$EncounterReason %in% c("Alive", "Found band only", "Predation/Scavenged",
                                    "Retrieved by pet", "Trap related", "Unknown", "Vehicle"),]

df$capyr <- year(as.Date(df$Date_Capture, "%m/%d/%Y"))
df <- df[df$capyr!=2024,]

### Group WMUs
# Table that defines grouping of WMUs
wmu <- c("1A","1B","2A","2B","2C","2D","2E","2F","2G","2H","3A","3B","3C","3D","4A","4B","4C","4D","4E","5A","5B","5C","5D")
grp <- c(1,1,2,3,2,2,2,4,4,4,4,4,5,6,7,7,8,7,8,9,9,10,10)
LookUp <- as.data.frame(cbind(wmu,grp))
df <- left_join(df,LookUp, by="wmu")

# filter groups I want
#df <- df[df$grp %in% c(2, 6, 7), ]
df <- df[df$wmu %in% c("2D", "3D", "4D", "5C"), ]
table(df$wmu)
##-----------------------------------------------------#X
###Tabulate releases by Year and age
release <- as.data.frame(table(df$capyr,df$Turkey_Age))
release
##-----------------------------------------------------#X
### Set Year 1 to 2020
df$rYr <- df$Year-2019 #recovery year 99=not recovered
df$rYr <- ifelse(is.na(df$rYr),99,df$rYr)
df$cYr <- df$capyr-2019 #release year
#-----------------------------------------------------#X
### Tabulate recoveries by year of release/recovery and age
df %>%
  group_by(rYr,cYr,Turkey_Age) %>%
  summarise(n=n()) %>%
  mutate(total=sum(n))
##-----------------------------------------------------#X
### Create mark-recapture matrix MR
n.occasions <- max(df$cYr) # number release occasions
n.release <- dim(df)[1]
MR <- matrix(NA, ncol = n.occasions + 1, nrow = n.release)
# vector of occasion of marking
mark.occ <- df$cYr
recov.occ <- df$rYr
##-----------------------------------------------------#X
# Fill the CH matrix
for (i in 1:n.release) {
  MR[i,mark.occ[i]] <- 1 # write 1 at release occasion
  for (t in mark.occ[i]:n.occasions){
    if (t==recov.occ[i]) {MR[i,t+1] <- 1} # write 1 at recovery occasion
  } #t loop
} #i loop
MR[which(is.na(MR))] <- 0
##-----------------------------------------------------#X
### 8.2.2. Adapted model with constant parameters to allow for age-specific and time-specific rates
# Create vector with occasion of marking 
f <- apply(MR, 1, get.first)
##-----------------------------------------------------#X
# Create indicator array of which occasion an animal is juvenile
I <- array(0,c(dim(MR)[1],n.occasions))
for(i in 1:dim(MR)[1]){ 
  if(df$Turkey_Age[i]=="J"){
    I[i,f[i]] <- 1
  }
}
##-----------------------------------------------------#X
# Create indicator vector of whether animal has a reward band: 1=non-reward band
II <- ifelse(df$reward_band=="N",1,0)

# Create design matrix for age coefficient (beta) using a mean parameterization
j1 <- rep(1,dim(MR)[2]-1)
j2 <- diag(dim(MR)[2]-2)
j3 <- rep(-1,dim(MR)[2]-2)
j23 <- rbind(j2,j3)
J <- as.matrix(cbind(j1,j23)) ; colnames(J) <- NULL ; rownames(J) <- NULL
##-----------------------------------------------------#X
# Bundle data
known.state <- known.state.mr(MR)
is.na(known.state) <- NA
##-----------------------------------------------------#X
# Create a list to store the variables
male_data <- list(
  male.f = f,
  male.I = I,
  male.II = II,
  male.time.param = J,
  male.nind = dim(MR)[1],
  male.n.occasions = (dim(MR)[2]),
  male.z = known.state,
  male.y = MR,
  male.n.wmu = length(unique(df$wmu)),
  male.wmu = df$wmu 
)
##-----------------------------------------------------#X
# Define a directory where you want to save the RDS files
output_dir <- "Data/Research_IPM_setup-data/"

# Loop through the list and save each variable as a separate RDS file
for (var_name in names(male_data)) {
  variable <- male_data[[var_name]]
  file_name <- paste(output_dir, var_name, ".rds", sep = "")
  saveRDS(variable, file_name)
}
##-----------------------------------------------------#X

#  Done!