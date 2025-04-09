###############################################################################X
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#################### Integrated Population Model (IPM): ########################
#                 #---# PhD Dissertation: Complex IPM #---#
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

rm(list = ls())
gc()

# Library
library(dplyr)

#############################################################X
# Load in DRM and known fate data ----
#############################################################X
dir <- "Data/IPM_setup-data_full/"
source("Analysis/Scripts/00_IPM_funs.R")

# list files
x <- list.files(dir, full.names = T)[!list.files(dir) %in% c("Older_data", "test_wmu")]
# name files and remove .rds
names <- list.files(dir, full.names = F)[!list.files(dir) %in% c("Older_data", "test_wmu")] %>% 
  str_remove(pattern = ".rds")

# Load into env
myfiles = lapply(x, readRDS)
names(myfiles) <- names
list2env(myfiles,globalenv())

# remove what we don't need
rm(myfiles)
rm(names)



