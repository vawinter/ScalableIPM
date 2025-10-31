###############################################################################X
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#################### Integrated Population Model (IPM): ########################
#                 #---# PhD Dissertation: R_IPM #---#
#        Creating a Bayesian IPM to inform turkey management in PA
###                                                                         ###X
#                Data formatting for dead-recovery model (DRM)
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
source("Analysis/00_IPM_funs.R")

########################################################X
#               Male data structuring
########################################################X
# Data cleaning ----
##-----------------------------------------------------#X
#Release data
df2a <- read.csv("../../PSUTurkey/turkey_IPM/SpringTurkeyHarvest2025/Band Data Export Report_20250711.csv")  # release data
df2a <- rename(df2a,bandid="Band.ID")
df2a$bandid <- as.character(df2a$bandid)
df2b <- read.csv("../../PSUTurkey/turkey_IPM/SpringTurkeyHarvest2025/Reported Band Data Export Report_20250722_ForPSU.csv")  # recovery data
df2b <- rename(df2b, bandid="BandNumber")
df2b$bandid <- as.character(df2b$bandid)

##-----------------------------------------------------#X
# Merge and structure data frame ------------
##-----------------------------------------------------#X
# Grab columns of interest from both DF
df <- left_join(df2a[,c("bandid", "Reward.Band","Turkey.Sex","Turkey.Age","Transmitter.ID",
                        "WMU", "Date.Capture", "Weight", "Recapture")], 
                df2b[,c("bandid","EncounterReason","Year")], 
                by="bandid")

# Remove recaptures and unknown age
df <- df[df$Recapture!="Yes" & df$Turkey.Age!="U",] 
df <- df[df$Turkey.Sex!="F",]  # remove females
df <- df[is.na(df$EncounterReason) | df$EncounterReason=="Harvest",]  # harvested birds only
df$capyr <- year(as.Date(df$Date.Capture, "%m/%d/%Y"))
df <- df[df$capyr!=2025,]

##-----------------------------------------------------#X
# Group WMUs ---
##-----------------------------------------------------#X
# Table that defines grouping of WMUs
WMU <- c("1A","1B","2A","2B","2C","2D","2E","2F","2G","2H","3A","3B","3C","3D","4A","4B","4C","4D","4E","5A","5B","5C","5D")
grp <- c(1,1,2,3,2,2,2,4,4,4,4,4,5,6,7,7,8,7,8,9,9,10,10)
LookUp <- as.data.frame(cbind(WMU,grp))
df <- left_join(df,LookUp, by="WMU")

# Tabulate releases by Year and age
release <- as.data.frame(table(df$capyr,df$Turkey.Age))
release

# Set Year 1 to 2020
df$rYr <- df$Year-2019 #recovery year 99=not recovered
df$rYr <- ifelse(is.na(df$rYr),99,df$rYr)
df$cYr <- df$capyr-2019 #release year

# Tabulate recoveries by year of release/recovery and age
df %>% 
  group_by(rYr,cYr,Turkey.Age, WMU) %>%
  summarise(n=n()) %>%
  mutate(total=sum(n))

# Filter WMUs that pertain to this IPM
df <- df[df$WMU %in% c("2D", "3D", "4D"), ]
table(df$WMU)

##-----------------------------------------------------#X
# Variable formatting -----
##-----------------------------------------------------#X
# Tabulate releases by Year and age
release <- as.data.frame(table(df$capyr,df$Turkey.Age))
release
##-----------------------------------------------------#X
### Set Year 1 to 2020
df$rYr <- df$Year-2019 #recovery year 99=not recovered
df$rYr <- ifelse(is.na(df$rYr),99,df$rYr)
df$cYr <- df$capyr-2019 #release year
#-----------------------------------------------------#X
### Tabulate recoveries by year of release/recovery and age
df %>%
  group_by(rYr,cYr,Turkey.Age) %>%
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
# Fill the CH matrix ----
##-----------------------------------------------------#X
for (i in 1:n.release) {
  MR[i,mark.occ[i]] <- 1 # write 1 at release occasion
  for (t in mark.occ[i]:n.occasions){
    if (t==recov.occ[i]) {MR[i,t+1] <- 1} # write 1 at recovery occasion
  } #t loop
} #i loop
MR[which(is.na(MR))] <- 0

##-----------------------------------------------------#X
# 8.2.2. Adapted model with constant parameters ---- 
#  to  allow for age-specific and time-specific rates
##-----------------------------------------------------#X
# Create vector with occasion of marking 
f <- apply(MR, 1, get.first)
##-----------------------------------------------------#X
# Create indicator array of which occasion an animal is juvenile
I <- array(0,c(dim(MR)[1],n.occasions))
for(i in 1:dim(MR)[1]){ 
  if(df$Turkey.Age[i]=="J"){
    I[i,f[i]] <- 1
  }
}
##-----------------------------------------------------#X
# Create indicator vector of whether animal has a reward band: 1=non-reward band
II <- ifelse(df$Reward.Band=="N",1,0)

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
  male.n.wmu = length(unique(df$WMU)),
  male.wmu = df$WMU 
)
##-----------------------------------------------------#X
# Define a directory where you want to save the RDS files
output_dir <- "Data/Research_IPM_setup-data/2024_DRM/"

# Loop through the list and save each variable as a separate RDS file
for (var_name in names(male_data)) {
  variable <- male_data[[var_name]]
  file_name <- paste(output_dir, var_name, ".rds", sep = "")
  saveRDS(variable, file_name)
}
#########################################################X
##-----------------------------------------------------#X
# Note: we did not have enough recovered females to do the DRM, so we only
# need to save the following variables:
##-----------------------------------------------------#X
# Female variables ----
female_data <- list(
  female.n.occasions = (dim(MR)[2]),
  female.n.wmu = length(unique(df$wmu))
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

# Done!

