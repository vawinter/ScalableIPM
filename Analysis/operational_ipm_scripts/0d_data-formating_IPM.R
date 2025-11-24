###############################################################################X
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#################### Integrated Population Model (IPM): ########################
#                 #---# PhD Dissertation: Chapter 1 #---#
#        Creating a Bayesian IPM to inform turkey management in PA
###                                                                         ###X
# Duane model:
# Dead-recovery model for turkey banding data analysis of male gobblers, Pennsylvania
#
# Survival and recovery rates are modeled by age. No random effects by WMU so assumes
#   harvest rates are the same across WMUs.
#
# Using this script to format female data same as males
#
# Also organizing harvest data
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #X
# 
# Created by: Veronica A. Winter
# Created: October 2023
# Last edited: XX/XX/XXXX
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #X
###############################################################################X
rm(list = ls())
gc()

# Library
library(dplyr)

#############################################################X
# Load in DRM and known fate data ----
#############################################################X
dir <- "Data/Operational_IPM_setup-data/"
source("Analysis/00_IPM_funs.R")

# list files
x <- list.files(dir, full.names = T)[!list.files(dir) %in% c("Rdata")]
# name files and remove .rds
names <- list.files(dir, full.names = F)[!list.files(dir) %in% c("Rdata")] %>% 
  str_remove(pattern = ".rds")

# Load into env
myfiles = lapply(x, readRDS)
names(myfiles) <- names
list2env(myfiles,globalenv())

# remove what we don't need
rm(myfiles)
rm(names)

#############################################################X
#----------------------------------------------------------#
# location of Turkey database ----
#----------------------------------------------------------#
# access.dir <- "../../../TurkeyDatabase/"
# 
# # Read in tables from Access database
# db <- file.path(paste0(access.dir,"TurkeyDB.accdb"))
# ch <- odbcConnectAccess2007(db) # open channel to database
# dcap <- sqlFetch(ch,"captures", as.is = T)
# dcen <- sqlFetch(ch,"censors", as.is = T)
# dmor <- sqlFetch(ch,"mortalities", as.is = T)
# dtag <- sqlFetch(ch,"transmitter_status", as.is = T)
# close(ch) #close channel to database
# 
# # Survival ----
# # Grab necessary columns from dfs
# capt <- dcap[,c(1, 2, 3, 4,6:8,13:16,33,34)]
# # Save capt for cohort of who came into study when
# cohort <-  unique(dcap[,c(1,3)]) # only need id and month
# cens <- dcen[,c(1:4)]
# mor <- dmor[,c(1,12:14)]
# 
# d1 <- merge(capt, cens, by="bandid", all.x=TRUE, all.y=FALSE)
# df <- merge(d1, mor, by="bandid", all.x=TRUE, all.y=FALSE)
# 
# # Create encounter histories
# dat <- encounter_histories(df, filter_sex = "F")
# head(dat)
# tail(dat)