###############################################################################X
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#################### Integrated Population Model (IPM): ########################
#                 #---# PhD Dissertation: Complex IPM #---#
#        Creating a Bayesian IPM to inform turkey management in PA
###                                                                         ###X
# Formatting hunter harvest data the PGC received from game-take surveys
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #X
###############################################################################X


rm(list = ls())
gc()

# libraries
library(dplyr)

dat <- read.csv("Data/Banding_harv_data/FallSprHarvData_20240919.csv")
head(dat)

#... by season
spring <- dat %>% 
  filter(season =="spring")

fall <- dat %>% 
  filter(season == "fall")

# Spring Adult Gobblers
wmu_spring_adult <- spring %>% 
  group_by(age, year, WMU.Group) %>% 
  reframe(wmu_count = sum(count)) %>% 
  filter(age == "gobbler",
         WMU.Group %in% c("2D", "3D", "4D", "5C"),
         year %in% c("2019","2020", "2021", "2022", "2023")) %>% 
  select(year, WMU.Group, wmu_count) %>% 
  as.data.frame()

# Spring Juvenile Jakes
wmu_spring_juv <- spring %>% 
  group_by(age, year, WMU.Group) %>% 
  reframe(wmu_count = sum(count)) %>% 
  filter(age == "jake",    
         WMU.Group %in% c("2D", "3D", "4D", "5C"),
         year %in% c("2019","2020", "2021", "2022", "2023")) %>% 
  select(year, WMU.Group, wmu_count) %>% 
  as.data.frame()

# Fall Adult Hens
wmu_fall_adult <- fall %>% 
  filter(!is.na(count)) %>% 
  group_by(age, year, WMU.Group) %>% 
  reframe(wmu_count = sum(count))%>% 
  filter(age == "hen",
         WMU.Group %in% c("2D", "3D", "4D"),
         year %in% c("2019","2020", "2021", "2022", "2023")) %>%  
  select(year, WMU.Group, wmu_count) %>% 
  as.data.frame()

# Fall Juvenile (Non-Hens)
wmu_fall_juv <- fall %>% 
  filter(!is.na(count)) %>% 
  group_by(age, year, WMU.Group) %>% 
  reframe(wmu_count = sum(count))%>% 
  filter(!age == "hen",
         WMU.Group %in% c("2D", "3D", "4D"),
         year %in% c("2019","2020", "2021", "2022", "2023")) %>% 
  select(year, WMU.Group, wmu_count) %>% 
  as.data.frame()

# Inspect the range of year and WMU Group values
years <- unique(wmu_spring_juv$year)
wmu_groups <- unique(wmu_fall_juv$WMU.Group)
wmu_groups_m <- unique(wmu_spring_juv$WMU.Group)

# Initialize the matrices based on the ranges
harvest.ad.spring <- matrix(NA, nrow = length(years), ncol = length(wmu_groups_m))
harvest.juv.spring <- matrix(NA, nrow = length(years), ncol = length(wmu_groups_m))
harvest.ad.fall <- matrix(NA, nrow = length(years), ncol = length(wmu_groups))
harvest.juv.fall <- matrix(NA, nrow = length(years), ncol = length(wmu_groups))

# Create lookup tables for the year and WMU.Group indices
year_index <- match(wmu_spring_juv$year, years)
wmu_group_index <- match(wmu_spring_juv$WMU.Group, wmu_groups)

# Fill the matrices with appropriate harvest data
for (i in 1:nrow(wmu_spring_adult)) {
  row <- wmu_spring_adult[i, ]
  y_idx <- match(row$year, years)  # Get the row index
  w_idx <- match(row$WMU.Group, wmu_groups_m)  # Get the column index
  harvest.ad.spring[y_idx, w_idx] <- row$wmu_count
}

for (i in 1:nrow(wmu_spring_juv)) {
  row <- wmu_spring_juv[i, ]
  y_idx <- match(row$year, years)
  w_idx <- match(row$WMU.Group, wmu_groups_m)
  harvest.juv.spring[y_idx, w_idx] <- row$wmu_count
}

for (i in 1:nrow(wmu_fall_adult)) {
  row <- wmu_fall_adult[i, ]
  y_idx <- match(row$year, years)
  w_idx <- match(row$WMU.Group, wmu_groups)
  harvest.ad.fall[y_idx, w_idx] <- row$wmu_count
}

for (i in 1:nrow(wmu_fall_juv)) {
  row <- wmu_fall_juv[i, ]
  y_idx <- match(row$year, years)
  w_idx <- match(row$WMU.Group, wmu_groups)
  harvest.juv.fall[y_idx, w_idx] <- row$wmu_count
}



# ... save
saveRDS(harvest.ad.spring, "Data/Research_IPM_setup-data/harvest.ad.spring.rds")
saveRDS(harvest.juv.spring, "Data/Research_IPM_setup-data/harvest.juv.spring.rds")
saveRDS(harvest.ad.fall, "Data/Research_IPM_setup-data/harvest.ad.fall.rds")
saveRDS(harvest.juv.fall, "Data/Research_IPM_setup-data/harvest.juv.fall.rds")







