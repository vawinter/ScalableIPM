########################################
#
# Dead-recovery model for turkey banding data analysis of male gobblers, Pennsylvania
#
# Survival and recovery rates are modeled by age. No random effects by WMU so assumes
#   harvest rates are the same across WMUs.
#
# July 2025
#########################################
rm(list=ls())
gc()
   options(max.print=1000000)
  # Sys.setenv(JAGS_HOME="C:/Program Files (x86)/JAGS/JAGS-4.3.1")
library(jagsUI)
library(MCMCvis)
library(RODBC)
library(stringr)
library(dplyr)
library(lubridate)
library(ggplot2)
library(data.table)


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



#### Read in data and select data used for analysis
df2a <- read.csv("../../PSUTurkey/turkey_IPM/SpringTurkeyHarvest2025/Band Data Export Report_20250711.csv")  # release data
        df2a <- rename(df2a,bandid="Band.ID")
        df2a$bandid <- as.character(df2a$bandid)
df2b <- read.csv("../../PSUTurkey/turkey_IPM/SpringTurkeyHarvest2025/Reported Band Data Export Report_20250722_ForPSU.csv")  # recovery data updated to include sex
        df2b <- rename(df2b, bandid="BandNumber")
        df2b$bandid <- as.character(df2b$bandid)
df <- left_join(df2a[,c(1,2,7,8,12,18,23,24,36)], df2b[,c(1,7,15)], by="bandid")
  df <- df[df$Recapture!="Yes" & df$Turkey.Age!="U",] # remove recaptures and unknown age
  df <- df[df$Turkey.Sex!="F",]  # remove females
  df <- df[is.na(df$EncounterReason) | df$EncounterReason=="Harvest",]  # harvested birds only
  df$capyr <- year(as.Date(df$Date.Capture, "%m/%d/%Y"))
  
### Group WMUs
    # Table that defines grouping of WMUs
      WMU <- c("1A","1B","2A","2B","2C","2D","2E","2F","2G","2H","3A","3B","3C","3D","4A","4B","4C","4D","4E","5A","5B","5C","5D")
	    grp <- c(1,1,2,3,2,2,2,4,4,4,4,4,5,6,7,7,8,7,8,9,9,10,10)
	    LookUp <- as.data.frame(cbind(WMU,grp))
	    
 df <- left_join(df,LookUp, by="WMU")
 df <- df[df$capyr!=2025,]
 # for now, remove 10
 df <- df[df$grp!=10,]
 table(df$grp)

###Tabulate releases by Year and age
 release <- as.data.frame(table(df$capyr,df$Turkey.Age))
 release

### Set Year 1 to 2020
    df$rYr <- df$Year-2019 #recovery year 99=not recovered
     df$rYr <- ifelse(is.na(df$rYr),99,df$rYr)
    df$cYr <- df$capyr-2019 #release year

### Tabulate recoveries by year of release/recovery and age
 df %>% group_by(rYr,cYr,Turkey.Age) %>%
          summarise(n=n()) %>%
          mutate(total=sum(n))
  
###Create mark-recapture matrix MR
    n.occasions <- max(df$cYr) # number release occasions
    n.release <- dim(df)[1]
    MR <- matrix(NA, ncol=n.occasions+1, nrow=n.release)
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
       get.first <- function(x) min(which(x!=0))
       f <- apply(MR, 1, get.first)
     # Create indicator array of which occasion an animal is juvenile
       I <- array(0,c(dim(MR)[1],n.occasions))
          for (i in 1:dim(MR)[1]) { if(df$Turkey.Age[i]=="J"){I[i,f[i]] <- 1}}
     # Create indicator vector of whether animal has a reward band: 1=non-reward band
       II <- ifelse(df$Reward.Band=="N",1,0)
     # Create design matrix for age coefficient (beta) using a mean parameterization
      j1 <- rep(1,dim(MR)[2]-1)
      j2 <- diag(dim(MR)[2]-2)
      j3 <- rep(-1,dim(MR)[2]-2)
      j23 <- rbind(j2,j3)
      J <- as.matrix(cbind(j1,j23)) ; colnames(J) <- NULL ; rownames(J) <- NULL

     
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
     output_dir <- "Data/Operational_IPM_setup-data/"
     
     # Loop through the list and save each variable as a separate RDS file
     for (var_name in names(male_data)) {
       variable <- male_data[[var_name]]
       file_name <- paste(output_dir, var_name, ".rds", sep = "")
       saveRDS(variable, file_name)
     }

# Female ----
     rm(list=ls())
     gc()
     options(max.print=1000000)
     # Sys.setenv(JAGS_HOME="C:/Program Files (x86)/JAGS/JAGS-4.3.1")
     library(jagsUI)
     library(MCMCvis)
     library(RODBC)
     library(stringr)
     library(dplyr)
     library(lubridate)
     library(ggplot2)
     library(data.table)
     
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
     
     #### Read in data and select data used for analysis
     df2a <- read.csv("../../PSUTurkey/turkey_IPM/SpringTurkeyHarvest2025/Band Data Export Report_20250711.csv")  # release data
     df2a <- rename(df2a,bandid="Band.ID")
     df2a$bandid <- as.character(df2a$bandid)
     df2b <- read.csv("../../PSUTurkey/turkey_IPM/SpringTurkeyHarvest2025/Reported Band Data Export Report_20250722_ForPSU.csv")  # recovery data updated to include sex
     df2b <- rename(df2b, bandid="BandNumber")
     df2b$bandid <- as.character(df2b$bandid)
     df <- left_join(df2a[,c(1,2,7,8,12,18,23,24,36)], df2b[,c(1,7,15)], by="bandid")
     df <- df[df$Recapture!="Yes" & df$Turkey.Age!="U",] # remove recaptures and unknown age
     df <- df[df$Turkey.Sex!="M",]  # remove females
     df <- df[is.na(df$EncounterReason) | df$EncounterReason=="Harvest",]  # harvested birds only
     df$capyr <- year(as.Date(df$Date.Capture, "%m/%d/%Y"))
     
     ### Group WMUs
     # Table that defines grouping of WMUs
     WMU <- c("1A","1B","2A","2B","2C","2D","2E","2F","2G","2H","3A","3B","3C","3D","4A","4B","4C","4D","4E","5A","5B","5C","5D")
     grp <- c(1,1,2,3,2,2,2,4,4,4,4,4,5,6,7,7,8,7,8,9,9,10,10)
     LookUp <- as.data.frame(cbind(WMU,grp))
     df <- left_join(df,LookUp, by="WMU")
     df <- df[df$capyr!=2025,]
     ### Group WMUs
     # For now, filter gropu 10
     df <- df[df$grp!=10,]
     
     ###Tabulate releases by Year and age
     release <- as.data.frame(table(df$capyr,df$Turkey.Age))
     release
     
     ### Set Year 1 to 2020
     df$rYr <- df$Year-2019 #recovery year 99=not recovered
     df$rYr <- ifelse(is.na(df$rYr),99,df$rYr)
     df$cYr <- df$capyr-2019 #release year
     
     ### Tabulate recoveries by year of release/recovery and age
     df %>% group_by(rYr,cYr,Turkey.Age) %>%
       summarise(n=n()) %>%
       mutate(total=sum(n))
     
     ###Create mark-recapture matrix MR
     n.occasions <- max(df$cYr) # number release occasions
     n.release <- dim(df)[1]
     MR <- matrix(NA, ncol=n.occasions+1, nrow=n.release)
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
     get.first <- function(x) min(which(x!=0))
     f <- apply(MR, 1, get.first)
     # Create indicator array of which occasion an animal is juvenile
     I <- array(0,c(dim(MR)[1],n.occasions))
     for (i in 1:dim(MR)[1]) { if(df$Turkey.Age[i]=="J"){I[i,f[i]] <- 1}}
     # Create indicator vector of whether animal has a reward band: 1=non-reward band
     II <- ifelse(df$Reward.Band=="N",1,0)
     # Create design matrix for age coefficient (beta) using a mean parameterization
     j1 <- rep(1,dim(MR)[2]-1)
     j2 <- diag(dim(MR)[2]-2)
     j3 <- rep(-1,dim(MR)[2]-2)
     j23 <- rbind(j2,j3)
     J <- as.matrix(cbind(j1,j23)) ; colnames(J) <- NULL ; rownames(J) <- NULL
     
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
     output_dir <- "Data/Operational_IPM_setup-data/"
     
     # Loop through the list and save each variable as a separate RDS file
     for (var_name in names(female_data)) {
       variable <- female_data[[var_name]]
       file_name <- paste(output_dir, var_name, ".rds", sep = "")
       saveRDS(variable, file_name)
     }
     