################################################X
#         Table formatting for manuscript      #X
## ========================================= ###X
# = Comparing model fit for full and simple    #X
##  ============ model outputs ============== ##X
#                Created by: VAW               #X
#                Date: 10/03/2024              #X
################################################X

rm(list = ls())
gc()

# Formatting table
library(knitr)
library(dplyr)
library(kableExtra)

# of data output
date <- "20241229"

# Load data 
kf_survival_df <- readRDS(paste0("Data/Output/Simple_", date,"_kf-survival_summary.rds"))
drm_harvest_df <- readRDS(paste0("Data/Output/Simple_", date, "_harvest_summary.rds"))
abundance_df <- readRDS(paste0("Data/Output/Simple_", date, "_abundance_summary.rds"))
drm_survival_df <- readRDS(paste0("Data/Output/Simple_", date, "_drm_survival_summary.rds"))
ppb <- readRDS(paste0("Data/Output/Simple_", date, "_ppb_summary.rds"))
hwb <- readRDS(paste0("Data/Output/Simple_", date, "_hwb_summary.rds"))
rec <- readRDS(paste0("Data/Output/Simple_", date, "_rec_summary.rds"))

# Load in WMU areas
wmu_areas <- readRDS("../turkey_IPM/Data/wmu_km_areas_w.groups.rds")

# Format harvest data ----
dat <- read.csv("Data/Banding_harv_data/FallSprHarvData_20240919.csv", header=T)
colnames(dat)
# If needed....
#dat$WMU.Group <- dat$ï..WMU.Group

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
  mutate(Region = paste("Region", WMU.Group, sep = " ")) %>% 
  left_join(wmu_areas, by = c("Region" = "WMU_ID")) %>% 
  filter(!year %in% c(2019, 2024)) 


## Structure harvest df ----
# Summarize data for total males, total females, and overall totals
harvest_df_totals <- Harvest %>%
  group_by(Region, year, sex) %>%
  summarise(
    median_value = sum(wmu_count),
    area_sq_km = first(area_sq_km),  # assuming area is the same across sex/age
    .groups = "drop"
  )

# Add 'Total' category by summing across sex and age_class 
harvest_df_overall <- harvest_df_totals %>%
  group_by(Region, year) %>%
  summarise(
    median_value = sum(median_value),
    area_sq_km = first(area_sq_km),  # assuming area is the same across sex/age
    .groups = "drop"
  ) %>%
  mutate(sex = "Total") %>%
  bind_rows(harvest_df_totals)  # Append the total rows to the original dataframe

# DRM harvest ----
drm_harvest_df2<- drm_harvest_df %>%
  mutate(Sex_Age_Class = paste(sex, age_class)) %>%
  dplyr::select(-c(sex, age_class, demographic_est)) %>%
  mutate(Year = case_when(year == 1 ~ "2020",
                          year == 2 ~ "2021",
                          year == 3 ~ "2022",
                          year == 4 ~ "2023")) %>%
  dplyr::select(-c(wmu, year)) %>%
  relocate(Sex_Age_Class, .before = median_value) %>%
  relocate(Region, .before = median_value) %>%
  relocate(Year, .before = median_value) %>% 
  rename(Median = median_value,
         "2.5% CI" = lower_ci,
         "97.5% CI" = upper_ci) %>%
  # Round numerical values to 3 decimal places
  mutate(across(c(Median, `2.5% CI`, `97.5% CI`), ~ round(., 3))) %>% 
  dplyr::select(-c(mean_value)) 

# Generate the table
drm_harvest_df_m <- drm_harvest_df2 %>% 
  filter(Sex_Age_Class %in% c("Male Adult", "Male Juvenile"))
table_hr_m <- kable(drm_harvest_df_m,
                  booktabs = TRUE,   
                  format = "latex", 
                  col.names = c("Sex & Age Class", "WMU Region", "Year","Median", "2.5% CI",  "97.5% CI"), 
                  caption = "Harvest rate with 95% credible intervals for each sex and age class across WMU groupss over four years.") %>%
  # Collapse rows for Sex/Age Class and WMU
  collapse_rows(columns = c(1, 2), latex_hline = "major", valign = "middle") %>%
  kable_styling(latex_options = c("hold_position", "repeat_header"), full_width = FALSE) 

table_hr_m

writeClipboard(table_hr_m)

ddrm_harvest_df_f <- drm_harvest_df2 %>% 
  filter(Sex_Age_Class %in% c("Female Adult", "Female Juvenile"))
table_hr_f <- kable(drm_harvest_df_m,
                    booktabs = TRUE,   
                    format = "latex", 
                    col.names = c("Sex & Age Class", "WMU Group", "Year","Median", "2.5% CI",  "97.5% CI"), 
                    caption = "Harvest rate with 95% credible intervals for each sex and age class across WMU groupss over four years.") %>%
  # Collapse rows for Sex/Age Class and WMU
  collapse_rows(columns = c(1, 2), latex_hline = "major", valign = "middle") %>%
  kable_styling(latex_options = c("hold_position", "repeat_header"), full_width = FALSE) 

table_hr_f

writeClipboard(table_hr_f)

# DRM survival ----
drm_survival_df2 <- drm_survival_df %>%
  mutate(Sex_Age_Class = paste(sex, age_class)) %>%
  dplyr::select(-c(sex, age_class, demographic_est)) %>%
   mutate(Year = case_when(year == 1 ~ "2020",
                   year == 2 ~ "2021",
                   year == 3 ~ "2022",
                   year == 4 ~ "2023")) %>%
  dplyr::select(-c(wmu, year)) %>%
  relocate(Sex_Age_Class, .before = median_value) %>%
  relocate(Region, .before = median_value) %>%
  relocate(Year, .before = median_value) %>% 
  rename(Median = median_value,
         "2.5% CI" = lower_ci,
         "97.5% CI" = upper_ci) %>%
  # Round numerical values to 3 decimal places
  mutate(across(c(Median, `2.5% CI`, `97.5% CI`), ~ round(., 3))) %>% 
  dplyr::select(-c(mean_value)) 

# Generate the table
drm_survival_df_a <- drm_survival_df2 %>% 
  filter(Sex_Age_Class %in% c("Male Adult"))

table_sp_a <- kable(drm_survival_df_a,
                    booktabs = TRUE,   
                    format = "latex", 
                    col.names = c("Sex & Age Class", "WMU Region", "Year","Median", "2.5% CI",  "97.5% CI"), 
                    caption = "Adult male survival probabilities with 95% credible intervals WMUs over four years.") %>%
  # Collapse rows for Sex/Age Class and WMU
  collapse_rows(columns = c(1, 2), latex_hline = "major", valign = "middle") %>%
  kable_styling(latex_options = c("hold_position", "repeat_header"), full_width = FALSE) 


writeClipboard(table_sp_a)
# Juv
drm_survival_df_j <- drm_survival_df2 %>% 
  filter(Sex_Age_Class %in% c("Male Juvenile"))

table_sp_j <- kable(drm_survival_df_j,
                    booktabs = TRUE,   
                    format = "latex", 
                    col.names = c("Sex & Age Class", "WMU Region", "Year","Median", "2.5% CI",  "97.5% CI"), 
                    caption = "Juvenile male survival probabilities with 95% credible intervals WMU regions over four years.") %>%
  # Collapse rows for Sex/Age Class and WMU
  collapse_rows(columns = c(1, 2), latex_hline = "major", valign = "middle") %>%
  kable_styling(latex_options = c("hold_position", "repeat_header"), full_width = FALSE) 


writeClipboard(table_sp_j)

# PPB ----
ppb2 <- ppb %>%
  select(-c(sex, age_class)) %>%
  mutate(Year = case_when(year == 1 ~ "2020",
                          year == 2 ~ "2021",
                          year == 3 ~ "2022",
                          year == 4 ~ "2023")) %>%
  select(-c(wmu, year)) %>%
  relocate(Region, .before = median_value) %>%
  relocate(Year, .before = median_value) %>% 
  rename(Median = median_value,
         "2.5% CI" = lower_ci,
         "97.5% CI" = upper_ci) %>%
  # Round numerical values to 3 decimal places
  mutate(across(c(Median, `2.5% CI`, `97.5% CI`), ~ round(., 3))) %>% 
  select(-c(mean_value)) 

# Generate the table
table_ppb <- kable(ppb2,
                  booktabs = TRUE,   
                  format = "latex", 
                  col.names = c("WMU Region", "Year","Median", "2.5% CI",  "97.5% CI"), 
                  caption = "Poult:Hen ratios with 95% credible intervals across WMU regions over four years.") %>%
  # Collapse rows for Sex/Age Class and WMU
  collapse_rows(columns = c(1, 2), latex_hline = "major", valign = "middle") %>%
  kable_styling(latex_options = c("hold_position", "repeat_header"), full_width = FALSE) 


writeClipboard(table_ppb)


# hwb
hwb2 <- hwb %>%
  select(-c(sex, age_class)) %>%
  mutate(Year = case_when(year == 1 ~ "2020",
                          year == 2 ~ "2021",
                          year == 3 ~ "2022",
                          year == 4 ~ "2023")) %>%
  select(-c(wmu, year)) %>%
  relocate(Region, .before = median_value) %>%
  relocate(Year, .before = median_value) %>% 
  rename(Median = median_value,
         "2.5% CI" = lower_ci,
         "97.5% CI" = upper_ci) %>%
  # Round numerical values to 3 decimal places
  mutate(across(c(Median, `2.5% CI`, `97.5% CI`), ~ round(., 3))) %>% 
  select(-c(mean_value)) 

# Generate the table
table_hwb <- kable(ppb2,
                   booktabs = TRUE,   
                   format = "latex", 
                   col.names = c("WMU Region", "Year","Median", "2.5% CI",  "97.5% CI"), 
                   caption = "Proportion of hens with a brood with 95% credible intervals across WMU regions over four years.") %>%
  # Collapse rows for Sex/Age Class and WMU
  collapse_rows(columns = c(1, 2), latex_hline = "major", valign = "middle") %>%
  kable_styling(latex_options = c("hold_position", "repeat_header"), full_width = FALSE) 


writeClipboard(table_hwb)


# rec
rec2 <- rec %>%
  select(-c(sex, age_class)) %>%
  mutate(Year = case_when(year == 1 ~ "2020",
                          year == 2 ~ "2021",
                          year == 3 ~ "2022",
                          year == 4 ~ "2023"),
         lower_ci  = lower_ci/area_sq_km,
         upper_ci  = upper_ci/area_sq_km) %>%
  select(-c(wmu, year, median_value, mean_value)) %>%
  relocate(Year, .before = density_value) %>% 
  relocate(group, .before = lower_ci ) %>% 
  relocate(Year, .before = lower_ci) %>% 
  relocate(density_value, .before = lower_ci) %>% 
  relocate(area_sq_km, .before= Year) %>% 
  rename(Density = density_value,
         "2.5% CI" = lower_ci,
         "97.5% CI" = upper_ci,
         "Area $km^{2}$" = area_sq_km) %>%
  # Round numerical values to 3 decimal places
  mutate(across(c(Density, `2.5% CI`, `97.5% CI`), ~ round(., 3))) 

# Generate the table
table_rec <- kable(rec2,
                   booktabs = TRUE,   
                   format = "latex", 
                   col.names = c("WMU Group", "Area $km^{2}$", "Year","Density", "2.5% CI",  "97.5% CI"), 
                   caption = "Recruitment density with 95% credible intervals across WMU groups over four years.") %>%
  # Collapse rows for Sex/Age Class and WMU
  collapse_rows(columns = c(1, 2), latex_hline = "major", valign = "middle") %>%
  kable_styling(latex_options = c("hold_position", "repeat_header"), full_width = FALSE) 


writeClipboard(table_rec)

# Abundance
abundance_df2<- abundance_df %>%
  dplyr::select(!area_sq_km) %>% 
  left_join(wmu_areas, by = c("Region" = "WMU_ID")) %>%  
  mutate(density_value = median_value / area_sq_km) %>% 
  mutate(Sex_Age_Class = paste(sex, age_class)) %>%
  select(-c(sex, age_class, demographic_est, mean_value)) %>%
  mutate(Year = case_when(year == 1 ~ "2020",
                          year == 2 ~ "2021",
                          year == 3 ~ "2022",
                          year == 4 ~ "2023"),
         lower_ci  = lower_ci/area_sq_km,
         upper_ci  = upper_ci/area_sq_km) %>%
  select(-c(wmu, year, median_value)) %>%
  relocate(Year, .before = density_value) %>% 
  relocate(Region, .before = lower_ci ) %>% 
  relocate(Year, .before = lower_ci) %>% 
  relocate(density_value, .before = lower_ci) %>% 
  relocate(area_sq_km, .before= Year) %>% 
  relocate(Sex_Age_Class, .before = Region) %>% 
  rename(Density = density_value,
         "2.5% CI" = lower_ci,
         "97.5% CI" = upper_ci,
         "`Area (km^2)`" = area_sq_km) %>%
  # Round numerical values to 3 decimal places
  mutate(across(c(Density, `2.5% CI`, `97.5% CI`), ~ round(., 3))) 

# Generate the table
abundance_am <- abundance_df2 %>% 
  filter(Sex_Age_Class %in% c("Male Adult"))
table_am <- kable(abundance_am,
                    booktabs = TRUE,   
                    format = "latex", 
                    col.names = c("Sex & Age Class","WMU Region", "`Area (km^2)`", "Year","Density", "2.5% CI",  "97.5% CI"), 
                    caption = "Adult male abundance density per abundance/$km^{2}$, with 95% credible intervals WMUs over four years.") %>%
  # Collapse rows for Sex/Age Class and WMU
  collapse_rows(columns = c(1, 2, 3), latex_hline = "major", valign = "middle") %>%
  kable_styling(latex_options = c("hold_position", "repeat_header"), full_width = FALSE) 


writeClipboard(table_am)

# Juv
abundance_jm <- abundance_df2 %>% 
  filter(Sex_Age_Class %in% c("Male Juvenile"))

table_jm <- kable(abundance_jm,
                  booktabs = TRUE,   
                  format = "latex", 
                  col.names = c("Sex & Age Class", "WMU Group", "`Area (km^2)`", "Year","Density", "2.5% CI",  "97.5% CI"), 
                  caption = "Juvenile male abundance density per abundance/$km^{2}$, with 95% credible intervals WMUs over four years.") %>%
  # Collapse rows for Sex/Age Class and WMU
  collapse_rows(columns = c(1, 2, 3), latex_hline = "major", valign = "middle") %>%
  kable_styling(latex_options = c("hold_position", "repeat_header"), full_width = FALSE) 


writeClipboard(table_jm)

# female
# Generate the table
abundance_af <- abundance_df2 %>% 
  filter(Sex_Age_Class %in% c("Female Adult")) 

table_af <- kable(abundance_af,
                  booktabs = TRUE,   
                  format = "latex", 
                  col.names = c("Sex & Age Class", "WMU Group", "`Area (km^2)`", "Year","Density", "2.5% CI",  "97.5% CI"), 
                  caption = "Adult female abundance density per abundance/$km^{2}$, with 95% credible intervals WMUs over four years.") %>%
  # Collapse rows for Sex/Age Class and WMU
  collapse_rows(columns = c(1, 2, 3), latex_hline = "major", valign = "middle") %>%
  kable_styling(latex_options = c("hold_position", "repeat_header"), full_width = FALSE) 


writeClipboard(table_af)

# Juv
abundance_jf <- abundance_df2 %>% 
  filter(Sex_Age_Class %in% c("Female Juvenile"))

table_jf <- kable(abundance_jf,
                  booktabs = TRUE,   
                  format = "latex", 
                  col.names = c("Sex & Age Class", "WMU Group", "`Area (km^2)`", "Year","Density", "2.5% CI",  "97.5% CI"), 
                  caption = "Juvenile female abundance density per abundance/$km^{2}$, with 95% credible intervals WMUs over four years.") %>%
  # Collapse rows for Sex/Age Class and WMU
  collapse_rows(columns = c(1, 2, 3), latex_hline = "major", valign = "middle") %>%
  kable_styling(latex_options = c("hold_position", "repeat_header"), full_width = FALSE) 


writeClipboard(table_jf)
