#####################################################X
#          Table formatting for manuscript          #X
## ============================================== ###X
# = Comparing model fit for R-IPM, O-IPM, and V_IPM #X
##  =============== model outputs ================ ##X
#####################################################X

rm(list = ls())
gc()

# Formatting table
library(knitr)
library(dplyr)
library(kableExtra)

# Load data 
kf_survival_df <- readRDS("Data/Output/O_24_kf-survival_summary.rds")
drm_harvest_df <- readRDS("Data/Output/O_24_harvest_summary.rds")
abundance_df <- readRDS("Data/Output/O_24_abundance_summary.rds")
drm_survival_df <- readRDS("Data/Output/O_24_drm_survival_summary.rds")
ppb <- readRDS("Data/Output/O_24_ppb_summary.rds")
hwb <- readRDS("Data/Output/O_24_hwb_summary.rds")
rec <- readRDS("Data/Output/O_24_rec_summary.rds")
combined_survival_df <- readRDS("Data/Output/O_24_comb-survival_summary.rds")
# Load in WMU areas
wmu_areas <- readRDS("Data/wmu_km_areas_regions.rds")

# DRM harvest ----
drm_harvest_df2<- drm_harvest_df %>%
  mutate(Sex_Age_Class = paste(sex, age_class),
         width = upper_ci - lower_ci) %>%
  dplyr::select(-c(sex, age_class, demographic_est)) %>%
  mutate(Year = case_when(year == 1 ~ "2020",
                          year == 2 ~ "2021",
                          year == 3 ~ "2022",
                          year == 4 ~ "2023",
                          year == 5 ~ "2024")) %>%
  dplyr::select(-c(wmu, year)) %>%
  relocate(Sex_Age_Class, .before = median_value) %>%
  relocate(Region, .before = median_value) %>%
  relocate(Year, .before = median_value) %>% 
  rename(Median = median_value,
         "2.5% CI" = lower_ci,
         "97.5% CI" = upper_ci) %>%
  # Round numerical values to 3 decimal places
  mutate(across(c(Median, `2.5% CI`, `97.5% CI`, width), ~ round(., 3))) %>% 
  dplyr::select(-c(mean_value)) 

# Summary statistics by sex-age class and WMU
harv_summary <- drm_harvest_df2 %>% 
  group_by(Sex_Age_Class) %>% 
  summarise(
    mean_harv = mean(Median),
    median_harv = median(Median),
    min_harv = min(Median),
    max_harv = max(Median),
    .groups = "drop"
  )


# Generate the table
table_hr_m <- kable(drm_harvest_df2,
                  booktabs = TRUE,   
                  format = "html", 
                  col.names = c("Sex & Age Class", "WMU Region", "Year","Median", "2.5% CI",  "97.5% CI", "CrI Width"), 
                  caption = "Harvest rate with 95% credible intervals for each sex and age class across WMU groupss over four years.") %>%
  # Collapse rows for Sex/Age Class and WMU
  collapse_rows(columns = c(1, 2), latex_hline = "major", valign = "middle") %>%
  kable_styling(latex_options = c("hold_position", "repeat_header"), full_width = FALSE) 

table_hr_m

save_kable(table_hr_m, file = "SubmissionMaterial/MajorRevisions/Appendix3Table1.html")

# DRM survival ----
drm_survival_df2 <- combined_survival_df %>%
  mutate(Sex_Age_Class = paste(sex, age_class),
         width = upper_ci - lower_ci) %>%
  dplyr::select(-c(sex, age_class, demographic_est)) %>%
   mutate(Year = case_when(year == 1 ~ "2020",
                   year == 2 ~ "2021",
                   year == 3 ~ "2022",
                   year == 4 ~ "2023",
                   year == 5 ~ "2024")) %>%
  dplyr::select(-c(wmu, year)) %>%
  relocate(Sex_Age_Class, .before = median_value) %>%
  relocate(Region, .before = median_value) %>%
  relocate(Year, .before = median_value) %>% 
  rename(Median = median_value,
         "2.5% CI" = lower_ci,
         "97.5% CI" = upper_ci) %>%
  # Round numerical values to 3 decimal places
  mutate(across(c(Median, `2.5% CI`, `97.5% CI`, width), ~ round(., 3))) %>% 
  dplyr::select(-c(mean_value)) 

# Generate the table
table_sp_a <- kable(drm_survival_df2,
                    booktabs = TRUE,   
                    format = "html", 
                    col.names = c("Sex & Age Class", "WMU Region", "Year","Median", "2.5% CI",  "97.5% CI", "CrI Width"), 
                    caption = "Adult male survival probabilities with 95% credible intervals WMUs over four years.") %>%
  # Collapse rows for Sex/Age Class and WMU
  collapse_rows(columns = c(1, 2), latex_hline = "major", valign = "middle") %>%
  kable_styling(latex_options = c("hold_position", "repeat_header"), full_width = FALSE) 


save_kable(table_sp_a, file = "SubmissionMaterial/MajorRevisions/Appendix3Table2.html")

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
                  format = "html", 
                  col.names = c("WMU Region", "Year","Median", "2.5% CI",  "97.5% CI"), 
                  caption = "Poult:Hen ratios with 95% credible intervals across WMU regions over four years.") %>%
  # Collapse rows for Sex/Age Class and WMU
  collapse_rows(columns = c(1, 2), latex_hline = "major", valign = "middle") %>%
  kable_styling(latex_options = c("hold_position", "repeat_header"), full_width = FALSE) 


writeClipboard(table_ppb)
save_kable(table_ppb, file = "SubmissionMaterial/MajorRevisions/Appendix3Table4.html")


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
                          year == 4 ~ "2023",
                          year == 5 ~ "2024"),
         lower_ci  = lower_ci/area_sq_km,
         upper_ci  = upper_ci/area_sq_km) %>%
  select(-c(wmu, year, median_value, mean_value)) %>%
  relocate(Year, .before = density_value) %>% 
  relocate(Region, .before = lower_ci ) %>% 
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
                   format = "html", 
                   col.names = c("Region", "Area $km^{2}$", "Year","Density", "2.5% CI",  "97.5% CI"), 
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
                          year == 4 ~ "2023",
                          year == 5 ~ "2024"),
         lower_ci  = lower_ci/area_sq_km,
         upper_ci  = upper_ci/area_sq_km,
         width = upper_ci - lower_ci) %>%
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
  mutate(across(c(Density, `2.5% CI`, `97.5% CI`, width), ~ round(., 3))) 

# Summary statistics by sex-age class and WMU
n_summary <- abundance_df2 %>% 
  group_by(Sex_Age_Class) %>% 
  summarise(
    mean_n = mean(Density),
    median_n = median(Density),
    min_n = min(Density),
    max_n = max(Density),
    .groups = "drop"
  )


# Generate the table
table_am <- kable(abundance_df2,
                    booktabs = TRUE,   
                    format = "html", 
                    col.names = c("Sex & Age Class","WMU Region", "`Area (km^2)`", "Year","Density", "2.5% CI",  "97.5% CI", "CrI Width"), 
                    caption = "Adult male abundance density per abundance/$km^{2}$, with 95% credible intervals WMUs over four years.") %>%
  # Collapse rows for Sex/Age Class and WMU
  collapse_rows(columns = c(1, 2, 3), latex_hline = "major", valign = "middle") %>%
  kable_styling(latex_options = c("hold_position", "repeat_header"), full_width = FALSE) 


save_kable(table_am, file = "SubmissionMaterial/MajorRevisions/Appendix3Table3.html")


