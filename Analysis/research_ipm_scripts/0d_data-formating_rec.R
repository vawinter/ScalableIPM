###############################################################################X
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
############## Research Integrated Population Model (R_IPM): ##################X
#                     #---# PhD Dissertation: R_IPM #---#
#     Data formatting: Creating data for hens with brood (HWB) & 
#           poults per brood (PPB) per wildlife management unit (WMU) 
#                         to inform recruitment
###                                                                         ###X
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #X
# Based off frequentest (glmmTB) versions from Diefenbach et al. 2025
# Created by: Veronica A. Winter
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #X
###############################################################################X
# clean env
rm(list = ls())
gc()
###-----------------------------------------------------#X
# Set seed for reproducibility
set.seed(123)
###-----------------------------------------------------#X
# Load libraries
library(dplyr)
library(lubridate)
library(tidyr)
library(nimble)
library(ggplot2)

###-----------------------------------------------------#X
## Data organizing ----
###-----------------------------------------------------#X
# Read in data file
df.pa.old <- readRDS("../../PSUTurkey/turkey_IPM/Data/PoultHen/PA.rds")
pa2 <- read.csv("Data/TurkeySightingSurvey_2024.csv")

# Clean up new 2024 data
pa.new <- pa2 %>% 
  mutate(Date = as.Date(SightingDate, format = "%m/%d/%Y"),
         Year = year(Date),
         Month = month(Date),
         doy = yday(Date),
         State = "PA",
         Seen.Before = case_when(IsSameTurkey == "Yes" ~ "Y", 
                                 TRUE ~ "N")) 
# add in groups
WMU <- c("1A","1B","2A","2B","2C","2D","2E","2F","2G","2H","3A","3B","3C","3D","4A","4B","4C","4D","4E","5A","5B","5C","5D")
MU <- c(1,1,2,3,2,2,2,4,4,4,4,4,5,6,7,7,8,7,8,9,9,10,10)
LookUp <- as.data.frame(cbind(WMU,MU))
pa.new <- left_join(pa.new,LookUp, by="WMU")
head(pa.new)

# Get recruitment metrics
pa.df <- pa.new %>% 
  group_by(WMU) %>% 
  mutate(
    PHratio = case_when(
      Hens > 0 ~ round(Poults / Hens, 2), # keep two decimal places
      TRUE    ~ NA_real_
    ),
    HWB = if_else(!is.na(PHratio) & PHratio > 0, 1, 0),
    Total = Hens+Poults
  ) %>% 
  ungroup() %>% 
  as.data.frame() 

head(pa.df)

# Reorder columns
new_order <- c("Year", "Date", "State", "County", "Township", "Hens", "Poults", "Males", "Unknown",
               "Seen.Before", "Month", "MU", "Total", "PHratio", "HWB", "doy")

pa.df.new <- pa.df[, new_order]
head(pa.df.new)

# Combine with other dataset
df.pa <- df.pa.old %>% 
  bind_rows(pa.df)
table(df.pa$Year)
tail(df.pa)

############### Select State for  analysis
df1 <- df.pa
# filter too many unknown sex-age and too many poults/hen
df1 <- df1[df1$Unknown/df1$Total<0.25 & df1$PHratio<=16,] 
# filter too many hens without poults
df1 <- df1[!(df1$Hens>=8 & df1$Poults==0),] 
head(df1)
table(df1$Year)

# Filter any weird NAs
df1 <- df1 %>% 
  filter(!is.na(Year))

# Scale day of year
df1$doy.scale <- scale(df1$doy)
df1$doy.2 <- df1$doy.scale^2
df.pph <- df1
# Scale day of year (stand-alone variables)
df.pph$doy.scale <- scale(df.pph$doy)
df.pph$doy.2 <- df.pph$doy.scale^2

ppb.aug31 <- (244 - attr(scale(df.pph$doy), "scaled:center")) / attr(scale(df.pph$doy), "scaled:scale")
ppb.aug31.2 <- ppb.aug31^2

# Remove data before 2019 and no poults
df.all <- df1[df1$Year>2018 & df1$PHratio>0,]
df.all$Year <- as.factor(df.all$Year)
df.all$MU <- as.factor(df.all$MU) # this is by WMU

# Check for NA values
sum(is.na(df.all))

# Check for infinite values
hist(df.all$PHratio)
hist(df.all$doy.scale)
df.all.pph <-  df.all %>%
  arrange(Year)

# Creating 'dummy' variables for year in nimble model
df.all.pph$Year2019 <- ifelse(df.all.pph$Year == 2019, 1, 0)
df.all.pph$Year2020 <- ifelse(df.all.pph$Year == 2020, 1, 0)
df.all.pph$Year2021 <- ifelse(df.all.pph$Year == 2021, 1, 0)
df.all.pph$Year2022 <- ifelse(df.all.pph$Year == 2022, 1, 0)
df.all.pph$Year2023 <- ifelse(df.all.pph$Year == 2023, 1, 0)
df.all.pph$Year2024 <- ifelse(df.all.pph$Year == 2024, 1, 0)

# Create a list to store the variables
ppb_data <- list(
  PHratio = df.all.pph$PHratio,
  ph.doy.scale = as.vector(df.all.pph$doy.scale),
  ph.doy.2 = as.vector(df.all.pph$doy.2),
  ph.Year2019 = df.all.pph$Year2019,
  ph.Year2020 = df.all.pph$Year2020, 
  ph.Year2021 = df.all.pph$Year2021,
  ph.Year2022 = df.all.pph$Year2022,
  ph.Year2023 = df.all.pph$Year2023,
  ph.Year2024 = df.all.pph$Year2024,
  ph.N = nrow(df.all.pph),
  ph.J = length(unique(df.all.pph$MU)),
  ph.wmu = as.integer(df.all.pph$MU),
  ppb.aug31 = ppb.aug31,
  ppb.aug31.2 = ppb.aug31.2
)

# Define a directory where you want to save the RDS files
output_dir <- "Data/Research_IPM_setup-data/"

# Loop through the list and save each variable as a separate RDS file
for (var_name in names(ppb_data)) {
  variable <- ppb_data[[var_name]]
  file_name <- paste(output_dir, var_name, ".rds", sep = "")
  saveRDS(variable, file_name)
}


# HWB ----
# clean env
rm(list = ls())
gc()
###-----------------------------------------------------#X
# Set seed for reproducibility
set.seed(123)
###-----------------------------------------------------#X
# Load libraries
library(dplyr)
library(lubridate)
library(tidyr)
library(nimble)
library(ggplot2)
###-----------------------------------------------------#X
# Read in data file
df.pa.old <- readRDS("../../PSUTurkey/turkey_IPM/Data/PoultHen/PA.rds")
pa2 <- read.csv("Data/TurkeySightingSurvey_2024.csv")

# Clean up new 2024 data
pa.new <- pa2 %>% 
  mutate(Date = as.Date(SightingDate, format = "%m/%d/%Y"),
         Year = year(Date),
         Month = month(Date),
         doy = yday(Date),
         State = "PA") %>% 
  filter(!is.na(Year))

# add in groups
WMU <- c("1A","1B","2A","2B","2C","2D","2E","2F","2G","2H","3A","3B","3C","3D","4A","4B","4C","4D","4E","5A","5B","5C","5D")
MU <- c(1,1,2,3,2,2,2,4,4,4,4,4,5,6,7,7,8,7,8,9,9,10,10)
LookUp <- as.data.frame(cbind(WMU,MU))
pa.new <- left_join(pa.new,LookUp, by="WMU")

# Get recruitment metrics
pa.df <- pa.new %>% 
  group_by(WMU) %>% 
  mutate(
    PHratio = case_when(
      Hens > 0 ~ round(Poults / Hens, 2), # keep two decimal places
      TRUE    ~ NA_real_
    ),
    HWB = if_else(!is.na(PHratio) & PHratio > 0, 1, 0),
    Total = Hens+Poults
  ) %>% 
  ungroup() %>% 
  select(
    -c(WMU, County, SightingDate, IsSameTurkey, LandOwnership, Township)
  ) %>% 
  as.data.frame()

# Combine with other dataset
df.pa <- df.pa.old %>% 
  bind_rows(pa.df)
table(df.pa$Year)
###-----------------------------------------------------#X
## Data organizing ----
###-----------------------------------------------------#X
df.hwb <- df.pa
#filter too many unknown sex-age and too many poults/hen
df.hwb <- df.hwb[df.hwb$Unknown/df.hwb$Total < 0.25 & df.hwb$PHratio <= 16,] 
# Filter too many hens without poults
df.hwb <- df.hwb[!(df.hwb$Hens >= 8 & df.hwb$Poults == 0),] 
table(df.hwb$Year)

# Filter any weird NAs
df.hwb <- df.hwb %>% 
  filter(!is.na(Year))

# Create day of year variables
df.hwb$doy.scale <- as.vector(scale(df.hwb$doy))
df.hwb$doy.2 <- as.vector(df.hwb$doy.scale^2)

hwb.aug31 <- (244 - attr(scale(df.hwb$doy), "scaled:center")) / attr(scale(df.hwb$doy), "scaled:scale")
hwb.aug31.2 <- hwb.aug31^2

# Filter data for years after 2018
df.hwb <- df.hwb[df.hwb$Year > 2018,]
df.hwb$Year <- as.factor(df.hwb$Year)
df.hwb$MU <- as.integer(df.hwb$MU)

# Creating 'dummy' variables for year in nimble model
df.hwb$Year2019 <- ifelse(df.hwb$Year == 2019, 1, 0)
df.hwb$Year2020 <- ifelse(df.hwb$Year == 2020, 1, 0)
df.hwb$Year2021 <- ifelse(df.hwb$Year == 2021, 1, 0)
df.hwb$Year2022 <- ifelse(df.hwb$Year == 2022, 1, 0)
df.hwb$Year2023 <- ifelse(df.hwb$Year == 2023, 1, 0)
df.hwb$Year2024 <- ifelse(df.hwb$Year == 2024, 1, 0)

# Ensure it's an integer
df.hwb$Hens <- as.integer(df.hwb$Hens)

# Make each hen in a record a separate observation
df1.all.hwb <- uncount(df.hwb, Hens)
###-----------------------------------------------------#X
# Create a list to store the variables
hwb_data <- list(
  HWB = df1.all.hwb$HWB,
  hwb.Year2019 = df1.all.hwb$Year2019,
  hwb.Year2020 = df1.all.hwb$Year2020, # 2019 is intercept
  hwb.Year2021 = df1.all.hwb$Year2021,
  hwb.Year2022 = df1.all.hwb$Year2022,
  hwb.Year2023 = df1.all.hwb$Year2023,
  hwb.Year2024 = df1.all.hwb$Year2024,
  hwb.doy.scale = df1.all.hwb$doy.scale,
  hwb.doy.2 = df1.all.hwb$doy.2,
  hwb.N = nrow(df1.all.hwb),
  hwb.J = length(unique(df1.all.hwb$MU)),
  hwb.wmu = df1.all.hwb$MU,
  hwb.aug31 = hwb.aug31,
  hwb.aug31.2 = hwb.aug31.2
)

# Define a directory where you want to save the RDS files
output_dir <- "Data/Research_IPM_setup-data/"

# Loop through the list and save each variable as a separate RDS file
for (var_name in names(hwb_data)) {
  variable <- hwb_data[[var_name]]
  file_name <- paste(output_dir, var_name, ".rds", sep = "")
  saveRDS(variable, file_name)
}
