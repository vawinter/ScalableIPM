###############################################################################X
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#################### Integrated Population Model (IPM): ########################
#                 #---# PhD Dissertation: Simple IPM #---#
#        Creating a Bayesian IPM to inform turkey management in PA
###                                                                         ###X
# Reading in data for fitting complex IPM
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #X
# 
# Created by: Veronica A. Winter
# Created: October 2023
# Last edited: XX/XX/XXXX
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
source("Analysis/Scripts/00_IPM_funs.R")

# Data cleaning ----
# This is mostly pulled from Duane's 'SpTurkeyHR2020-22RE_abundance.R' script

#### Read in data and select data used for analysis
df2a <- read.csv("Data/Banding_harv_data/Band Data Export Report_20240909.csv")  # release data
df2a <- rename(df2a,bandid="band_id")
df2a$bandid <- as.character(df2a$bandid)
df2b <- read.csv("Data/Banding_harv_data/Reported Band Data Export Report_20240909_ForPSU.csv")  # recovery data
df2b <- rename(df2b, bandid="BandNumber")
df2b$bandid <- as.character(df2b$bandid)

# Data cleaning ----
df <- left_join(df2a[,c(1,2,6,7,11,17,22,35)], df2b[,c(1,6,7,15)], join_by(bandid))
df <- df[df$recapture!="Yes" & df$Turkey_Age!="U",] # remove recaptures and unknown age
df <- df[df$Turkey_Sex!="M",]  # remove females
df$capyr <- year(as.Date(df$Date_Capture, "%m/%d/%Y"))
df <- df[df$capyr!=2024,]

### Group WMUs
# Table that defines grouping of WMUs
wmu <- c("1A","1B","2A","2B","2C","2D","2E","2F","2G","2H","3A","3B","3C","3D","4A","4B","4C","4D","4E","5A","5B","5C","5D")
grp <- c(1,1,2,3,2,2,2,4,4,4,4,4,5,6,7,7,8,7,8,9,9,10,10)
LookUp <- as.data.frame(cbind(wmu,grp))
df <- left_join(df,LookUp, by="wmu")

# For now, filter gropu 10
df <- df[df$grp!=10,]

### relate group WMUs to study areas
df$group.index <- case_when(
  substr(df$wmu, 1, 1) == "1" ~ 7,
  df$wmu %in% c("2A", "2C", "2D", "2E") ~ 2,
  df$wmu == "2B" ~ 2,
  df$wmu %in% c("2F", "2G", "2H", "3A", "3B") ~ 7,
  df$wmu == "3C" ~ 6,
  df$wmu == "3D" ~ 6,
  df$wmu %in% c("4A", "4B", "4D") ~ 7,
  df$wmu %in% c("4C", "4E") ~ 6,
 # df$wmu %in% c("5A", "5B") ~ 7,
  TRUE ~ 7
) 
##-----------------------------------------------------#X
###Tabulate releases by Year and age
release <- as.data.frame(table(df$capyr,df$Turkey_Age))
release
##-----------------------------------------------------#X
### Set Year 1 to 2020
df$rYr <- df$Year-2019 #recovery year 99=not recovered
df$rYr <- ifelse(is.na(df$rYr),99,df$rYr)
df$cYr <- df$capyr-2019 #release year

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

# for now...
#df <- df %>% filter(wmu %in% c(2, 6, 7))

##-----------------------------------------------------#X
# Create a list to store the variables
female_data <- list(
  female.n.occasions = (dim(MR)[2]),
  female.n.wmu = length(unique(df$grp)),
  female.wmu = df$grp,
  female.grp.index = df$group.index
)
##-----------------------------------------------------#X
# Define a directory where you want to save the RDS files
output_dir <- "Data/Simple_IPM_setup-data/"

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
# for now, remove 10
df <- df[df$grp!=10,]
table(df$grp)
##-----------------------------------------------------#X
###Tabulate releases by Year and age
release <- as.data.frame(table(df$capyr,df$Turkey_Age))
release
##-----------------------------------------------------#X
### Set Year 1 to 2020
df$rYr <- df$Year-2019 #recovery year 99=not recovered
df$rYr <- ifelse(is.na(df$rYr),99,df$rYr)
df$cYr <- df$capyr-2019 #release year
##-----------------------------------------------------#X
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
  male.n.wmu = length(unique(df$grp)),
  male.wmu = df$grp
)
##-----------------------------------------------------#X
# Define a directory where you want to save the RDS files
output_dir <- "Data/Simple_IPM_setup-data/"

# Loop through the list and save each variable as a separate RDS file
for (var_name in names(male_data)) {
  variable <- male_data[[var_name]]
  file_name <- paste(output_dir, var_name, ".rds", sep = "")
  saveRDS(variable, file_name)
}
##-----------------------------------------------------#X

