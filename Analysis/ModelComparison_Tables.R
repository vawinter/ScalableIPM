# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#################### Integrated Population Model (IPM): #######################X
#                 #---# PhD Dissertation: Complex IPM #---#
#        Creating a Bayesian IPM to inform turkey management in PA
###                       *** Real data run ***                             ###X
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #X
###                                                                         ###X
#                             Correlations
###                                                                         ###X
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #X
# 
# Created by: Veronica A. Winter
# Last edited: 10/01/2024
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #X
###############################################################################X
rm(list = ls())
gc()

# Load libraries
library(knitr)
library(dplyr)
library(kableExtra)

# Are we savig for datavis or manu?
selected_dir <- "Manuscript/"

# of data output
date <- "20250326"
type <- "Full_"

# Load data
abundance_df <- readRDS(paste0("Data/Output/", type, date, "_abundance_summary.rds"))

#############################
# Format harvest data
dat <- read.csv("../../PSUTurkey/turkey_IPM/Data/Banding_harv_data/FallSprHarvData_20240919.csv", header=T)
colnames(dat)
# If needed....
# dat$WMU.Group <- dat$ï..WMU.Group

#... by season
spring <- dat %>% 
  filter(season =="spring")

fall <- dat %>% 
  filter(season == "fall")

# get specific harvest
# Spring Adult Gobblers
wmu_spring_adult <- spring %>% 
  group_by(age, year, WMU.Group) %>% 
  reframe(wmu_count = sum(count)) %>% 
  filter(age == "gobbler",
         WMU.Group %in% c("2D", "3D", "4D")) %>% 
  dplyr:: select(year, WMU.Group, wmu_count, age) %>% 
  as.data.frame()

# Spring Juvenile Jakes
wmu_spring_juv <- spring %>% 
  group_by(age, year, WMU.Group) %>% 
  reframe(wmu_count = sum(count)) %>% 
  filter(age == "jake",    
         WMU.Group %in% c("2D", "3D", "4D")) %>% 
  dplyr::select(year, WMU.Group, wmu_count, age) %>% 
  as.data.frame()

# Fall Adult Hens
wmu_fall_adult <- fall %>% 
  filter(!is.na(count)) %>% 
  group_by(age, year, WMU.Group) %>% 
  reframe(wmu_count = sum(count))%>% 
  filter(age == "hen",
         WMU.Group %in% c("2D", "3D", "4D")) %>%  
  dplyr::select(year, WMU.Group, wmu_count, age) %>% 
  as.data.frame()

# Fall Juvenile (Non-Hens)
wmu_fall_juv <- fall %>% 
  filter(!is.na(count)) %>% 
  group_by(age, year, WMU.Group) %>% 
  reframe(wmu_count = sum(count))%>% 
  filter(!age == "hen",
         WMU.Group %in% c("2D", "3D", "4D")) %>% 
  dplyr::select(year, WMU.Group, wmu_count, age) %>% 
  as.data.frame()

# Combine data frames
harvest <- rbind(wmu_spring_adult, wmu_spring_juv, wmu_fall_adult, wmu_fall_juv)

# Format
Harvest <- harvest %>% 
  mutate(
    sex = case_when(
      age == "gobbler" ~ "Male",
      age == "jake"    ~ "Male",
      age == "hen"     ~ "Female",
      age == "jenny"   ~ "Female",
      TRUE             ~ NA_character_  # Handles any other cases
    ),
    age_class = case_when(
      age == "gobbler" ~ "Adult",
      age == "jake"    ~ "Juvenile",
      age == "hen"     ~ "Adult",
      age == "jenny"   ~ "Juvenile",
      TRUE             ~ NA_character_  # Handles any other cases
    )
  ) %>% 
  filter(!year %in% c(2019, 2024)) %>% 
  rename(harvest = wmu_count, wmu = WMU.Group) %>% 
  select(Year = year, wmu, harvest, age_class, sex)

# Merge data together
combined_df <- abundance_df %>% 
  rename(abundance = median_value) %>% 
  select(year, wmu, abundance, age_class, sex) %>% 
  mutate(Year = as.integer(case_when(
    year == 1 ~ "2020",   
    year == 2 ~ "2021",   
    year == 3 ~ "2022",  
    year == 4 ~ "2023",
    TRUE ~ NA_character_  # Default to NA if year is not 1, 2, 3, or 4
  ))) %>% left_join(Harvest) %>% 
  select(-c("year"))%>% 
  filter(Year != 2020)

# Calculate R² values for each group
correlations_groups_r2 <- combined_df %>%
  mutate(Sex_Age_Class = paste(sex, age_class)) %>%
  group_by(Sex_Age_Class, wmu,) %>%
  summarise(
    # Calculate Pearson correlation and square it to get R²
    r_squared = summary(lm(harvest ~ abundance))$r.squared,
    .groups = "drop"
  ) 

table_wmu <- kable(correlations_groups_r2,
                    booktabs = TRUE,   
                    format = "html", 
                    col.names = c("Sex & Age Class","WMU",  "R2")) %>%
  # Collapse rows for Sex/Age Class and WMU
  collapse_rows(columns = c(1, 2), latex_hline = "major", valign = "middle") %>%
  kable_styling(latex_options = c("hold_position", "repeat_header"), full_width = FALSE) 

save_kable(table_wmu, file = "harvestr2research_table.html")
#----------------------------------------------X
# of data output
date <- "20250326"
type <- "Simple_"

# Load data
abundance_df <- readRDS(paste0("Data/Output/", type, date, "_abundance_summary.rds"))

#############################
# Format harvest data
dat <- read.csv("../../PSUTurkey/turkey_IPM/Data/Banding_harv_data/FallSprHarvData_20240919.csv", header=T)
colnames(dat)
# If needed....
# dat$WMU.Group <- dat$ï..WMU.Group

#... by season
spring <- dat %>% 
  filter(season =="spring")

fall <- dat %>% 
  filter(season == "fall")

# get specific harvest
# Spring Adult Gobblers
wmu_spring_adult <- spring %>% 
  group_by(age, year, WMU.Group) %>% 
  reframe(wmu_count = sum(count)) %>% 
  filter(age == "gobbler",
         WMU.Group %in% c(1:9)) %>% 
  dplyr:: select(year, WMU.Group, wmu_count, age) %>% 
  as.data.frame()

# Spring Juvenile Jakes
wmu_spring_juv <- spring %>% 
  group_by(age, year, WMU.Group) %>% 
  reframe(wmu_count = sum(count)) %>% 
  filter(age == "jake",    
         WMU.Group %in% c(1:9)) %>% 
  dplyr::select(year, WMU.Group, wmu_count, age) %>% 
  as.data.frame()

# Fall Adult Hens
wmu_fall_adult <- fall %>% 
  filter(!is.na(count)) %>% 
  group_by(age, year, WMU.Group) %>% 
  reframe(wmu_count = sum(count))%>% 
  filter(age == "hen",
         WMU.Group %in% c(1:9)) %>%  
  dplyr::select(year, WMU.Group, wmu_count, age) %>% 
  as.data.frame()

# Fall Juvenile (Non-Hens)
wmu_fall_juv <- fall %>% 
  filter(!is.na(count)) %>% 
  group_by(age, year, WMU.Group) %>% 
  reframe(wmu_count = sum(count))%>% 
  filter(!age == "hen",
         WMU.Group %in% c(1:9)) %>%  
  dplyr::select(year, WMU.Group, wmu_count, age) %>% 
  as.data.frame()

# Combine data frames
harvest <- rbind(wmu_spring_adult, wmu_spring_juv, wmu_fall_adult, wmu_fall_juv)

# Format
Harvest <- harvest %>% 
  mutate(
    sex = case_when(
      age == "gobbler" ~ "Male",
      age == "jake"    ~ "Male",
      age == "hen"     ~ "Female",
      age == "jenny"   ~ "Female",
      TRUE             ~ NA_character_  # Handles any other cases
    ),
    age_class = case_when(
      age == "gobbler" ~ "Adult",
      age == "jake"    ~ "Juvenile",
      age == "hen"     ~ "Adult",
      age == "jenny"   ~ "Juvenile",
      TRUE             ~ NA_character_  # Handles any other cases
    )
  ) %>% 
  filter(!year %in% c(2019, 2024)) %>% 
  rename(harvest = wmu_count, Region = WMU.Group) %>% 
  mutate(Region = as.character(paste("Region", Region, sep = " "))) %>% 
  select(Year = year, Region, harvest, age_class, sex)

# Merge data together
combined_df <- abundance_df %>% 
  rename(abundance = median_value) %>% 
  select(year, Region, abundance, age_class, sex) %>% 
  mutate(Year = as.integer(case_when(
    year == 1 ~ "2020",   
    year == 2 ~ "2021",   
    year == 3 ~ "2022",  
    year == 4 ~ "2023",
    TRUE ~ NA_character_  # Default to NA if year is not 1, 2, 3, or 4
  ))) %>% left_join(Harvest) %>% 
  select(-c("year"))%>% 
  filter(Year != 2020)

# Calculate R² values for each group
correlations_groups_r2 <- combined_df %>%
  mutate(Sex_Age_Class = paste(sex, age_class)) %>%
#  group_by(Sex_Age_Class, Region) %>%
  summarise(
    # Calculate Pearson correlation and square it to get R²
    r_squared = summary(lm(harvest ~ abundance))$r.squared,
    .groups = "drop"
  ) 

table_reg <- kable(correlations_groups_r2,
                   booktabs = TRUE,   
                   format = "html", 
                   col.names = c("Sex & Age Class","Region",  "R2")) %>%
  # Collapse rows for Sex/Age Class and WMU
  collapse_rows(columns = c(1, 2), latex_hline = "major", valign = "middle") %>%
  kable_styling(latex_options = c("hold_position", "repeat_header"), full_width = FALSE) 

table_reg
