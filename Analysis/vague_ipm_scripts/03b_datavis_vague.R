# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#################### Integrated Population Model (IPM): #######################X
#                 #---# PhD Dissertation: Simple IPM #---#
#        Creating a Bayesian IPM to inform turkey management in PA
###                       *** Real data run ***                             ###X
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #X
###                                                                         ###X
#                             Data visualization
###                              VAGUE PRIORS                               ###X
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #X
###############################################################################X
rm(list = ls())
gc()

# Load libraries
library(dplyr)
library(ggplot2)
library(tidyr)
library(knitr)
library(ggridges)
library(reshape2)
library(patchwork)


# Load data 
kf_survival_df <- readRDS("Data/Output/V_kf-survival_summary.rds")
drm_harvest_df <- readRDS("Data/Output/V_harvest_summary.rds")
abundance_df <- readRDS("Data/Output/V_abundance_summary.rds")
drm_survival_df <- readRDS("Data/Output/V_drm_survival_summary.rds")
combined_survival_df <- readRDS("Data/Output/V_comb-survival_summary.rds")
ppb <- readRDS("Data/Output/V_ppb_summary.rds")
hwb <- readRDS("Data/Output/V_hwb_summary.rds")
rec <- readRDS("Data/Output/V_rec_summary.rds")

# Load in WMU areas
wmu_areas <- readRDS("Data/wmu_km_areas_regions.rds")

# Format harvest data ----
dat <- read.csv("Data/Banding_harv_data/FallSprHarvData_20240919.csv", header=T)
colnames(dat)
# If needed....
#dat$WMU.Group <- dat$X...WMU.Group

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


# Color palettes ----
colors <- c("#EB781B", "#CC5221", "#71250F", "#6F6534", "#365365")
colors2 <- c("#F6D9C0", "#C74442", "#71250E", "#6D383E", "#365365")

groups <- c("Region 1" = "#3c1518", 
            "Region 2" = "#f2f3ae", 
            "Region 3" = "#0f5c5d", 
            "Region 4" = "#69140e",
            "Region 5" = "#C74442", 
            "Region 6" = "#5f0f40",
            "Region 7" = "#043A70", 
            "Region 8" = "goldenrod",
            "Region 9" = "#f79d65", 
            "Region 10" = "#606c38")

combo_colors <- c(
  "Adult Female"    = "#7A370A",
  "Juvenile Female" = "#EB781B",
  "Adult Male"      = "#043A50",
  "Juvenile Male"   = "#5C7391"
)

sex_colors <- c(
  "Female"    = "#6D383E",
  "Male"      = "#416E7D",
  "Total"   = "#6F6534"
)

## Structure abundance df ----
# Summarize data for total males, total females, and overall totals
abundance_df_totals <- abundance_df %>%
  dplyr::select(!area_sq_km) %>% 
  left_join(wmu_areas, by = c("Region" = "WMU_ID")) %>% 
  group_by(Region, year, sex) %>%
  summarise(
    median_value = sum(median_value),
    lower_ci = sum(lower_ci),
    upper_ci = sum(upper_ci),
    area_sq_km = first(area_sq_km),  # assuming area is the same across sex/age
    .groups = "drop"
  )

# Add 'Total' category by summing across sex and age_class 
abundance_df_overall <- abundance_df %>%
  dplyr::select(!area_sq_km) %>% 
  left_join(wmu_areas, by = c("Region" = "WMU_ID")) %>% 
  group_by(Region, year) %>%
  summarise(
    median_value = sum(median_value),
    lower_ci = sum(lower_ci),
    upper_ci = sum(upper_ci),
    area_sq_km = first(area_sq_km),  # assuming area is the same across sex/age
    .groups = "drop"
  ) %>%
  mutate(sex = "Total") %>%
  bind_rows(abundance_df_totals)  # Append the total rows to the original dataframe

# Plots ----
## Harvest rate ----
harvest_plot <- ggplot(drm_harvest_df, aes(y = median_value, x = year, 
                                           shape = sex, color = interaction(age_class, sex, sep = " "))) +
  geom_point(size = 4, position = position_dodge(width = 0.5)) +
  geom_errorbar(aes(ymin = lower_ci, ymax = upper_ci), 
                width = 0.3, 
                position = position_dodge(width = 0.5), show_guide = F) +
  facet_wrap(~Region, ncol = 9) +
  
  # Assign shapes to sex: female and male
  scale_shape_manual(values = c("Female" = 16, "Male" = 17)) +
  
  annotate("segment", x = Inf, xend = Inf, y = -Inf, yend = Inf, color = "black", linetype = "dashed") +
  
  # Assign colors to age_class and sex combinations
  scale_color_manual(values = c(
    "Adult Female"    = "#7A370A",
    "Juvenile Female" = "#EB781B",
    "Adult Male"      = "#043A50",
    "Juvenile Male"   = "#5C7391"
  )) +
  
  
  scale_x_continuous(breaks = c(1, 2, 3, 4), labels = c("2020", "2021", "2022", "2023")) +  # Change x-axis labels
  scale_y_continuous(limits = c(0, 1)) +  # Set y-axis limits
  
  labs(x = "Year", y = "Harvest rate", shape = "", color = "") +  # Modify legend titles
  
  guides(color = guide_legend(override.aes = list(shape = c(16, 16, 17, 17))),  # Override shapes in color legend
         shape = guide_none()) +  # Hide the separate shape legend
  
  theme_classic() +
  theme(
    text = element_text(size = 16),           # Increase the base text size
    axis.title = element_text(size = 16),     # Increase axis title size
    axis.text = element_text(size = 16),      # Increase axis text (ticks) size
    strip.text = element_text(size = 16),     # Increase facet label text size
    legend.title = element_text(size = 16, face = "bold"),   # Customize legend title text
    legend.text = element_text(size = 16),    # Customize legend item text
    legend.position = "top",                # Position the legend
    legend.key = element_rect(fill = "white"),  # Customize legend key appearance
    legend.key.size = unit(1.5, "lines"),      # Adjust size of legend keys (shapes)
    axis.title.x = element_blank(),    # Remove x-axis title
    axis.text.x = element_blank()       # Remove x-axis text
  )

## Survival DRM ----
survival_plot <- ggplot(combined_survival_df, aes(y = median_value, x = year, 
                                                  shape = sex, color = interaction(age_class, sex, sep = " "))) +
  geom_point(size = 4, position = position_dodge(width = 0.5)) +
  geom_errorbar(aes(ymin = lower_ci, ymax = upper_ci), 
                width = 0.3, 
                position = position_dodge(width = 0.5), show_guide = F) +
  facet_wrap(~Region, ncol = 9) +
  
  # Use the same colors from the harvest plot for "adult male" and "juvenile male"
  # Assign colors to age_class and sex combinations
  scale_color_manual(values = c(
    "Adult Female"    = "#7A370A",
    "Juvenile Female" = "#EB781B",
    "Adult Male"      = "#043A50",
    "Juvenile Male"   = "#5C7391"
  )) +
  
  annotate("segment", x = Inf, xend = Inf, y = -Inf, yend = Inf, color = "black", linetype = "dashed") +
  
  # Shapes for juvenile and adult (keeping the same shapes as before)
  scale_shape_manual(values = c("Female" = 16, "Male" = 17)) +
  
  # Remove shape legend using guides()
  guides(shape = "none") +
  guides(color = "none") +
  
  # X-axis adjustments
  scale_x_continuous(breaks = c(1, 2, 3, 4), labels = c("2020", "2021", "2022", "2023")) +
  
  # Y-axis limits for survival probabilities
  scale_y_continuous(limits = c(0, 1)) +  
  
  # Update labels to reflect the shape and color scheme
  labs(x = "Year", y = "Survival probability", shape = "", color = "") +  
  
  # Keep theme consistent with the previous plots
  theme_classic() +
  theme(
    text = element_text(size = 16),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 16),
    strip.text = element_blank(),      # Remove WMU labels
    legend.title = element_text(size = 16, face = "bold"),
    legend.text = element_text(size = 16),
    legend.key = element_rect(fill = "white"),  # Customize legend key appearance
    legend.position = "none",
    legend.key.size = unit(1.5, "lines"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )


### Stack harvest/drm surv ----
combined_plot <- harvest_plot / survival_plot 

# Display the combined plot
print(combined_plot)

#### Save harvest/drm surv ----
ggsave("Dataviz/V_IPM_plot.png", plot = combined_plot,  width = 15, height = 12, dpi = 700)

# Abundance Plot ----
abundance_plot <- ggplot(abundance_df_overall, aes(y = median_value / area_sq_km, x = year, shape = sex, color = sex)) +
  geom_point(size = 6, position = position_dodge(width = 0.5)) +
  geom_errorbar(aes(ymin = lower_ci / area_sq_km, ymax = upper_ci / area_sq_km), 
                width = 0.3, position = position_dodge(width = 0.5), show_guide = F) +
  facet_wrap(~Region, ncol = 9) +  # Use custom labeller
  
  # Custom shapes and colors for sex
  scale_shape_manual(values = c("Female" = 16, "Male" = 17, "Total" = 8)) + 
  scale_color_manual(values = c(
    "Female"    = "#6D383E",
    "Male"      = "#416E7D",
    "Total"   = "#6F6534"
  )) +
  
  annotate("segment", x = Inf, xend = Inf, y = -Inf, yend = Inf, color = "black", linetype = "dashed") +
  scale_x_continuous(breaks = c(1, 2, 3, 4), labels = c("2020", "2021", "2022", "2023")) +
  
  labs(x = "Year", y = expression('Abundance/km'^2*''), shape = "", color = "") +  # Modify legend title
  # # Combine the shape and color legends
  guides(color = guide_legend(override.aes = list(shape = c(16, 17, 8))),  # Override shapes in the color legend
         shape = guide_none()) +  # Hide the separate shape legend
  ylim(0, 10) +
  theme_classic() +
  theme(
    text = element_text(size = 16),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 16),
    legend.title = element_text(size = 16, face = "bold"),
    legend.text = element_text(size = 16),
    legend.position = "top",
    legend.key = element_rect(fill = "white"),  # Customize legend key appearance
    legend.key.size = unit(1.5, "lines"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )


#### Save abundance ----
ggsave("Dataviz/V_abun_plot.png", plot = abundance_plot,   width = 16, height = 9, dpi = 700)

## Harvest density ----
harvest_den_plot <- ggplot(harvest_df_overall, aes(y = median_value/area_sq_km, x = year, shape = sex, color = sex))  +
  
  # Plot points with dodge position to avoid overlap
  geom_point(size = 4, position = position_dodge(width = 0.5)) +
  
  # Facet by WMU
  facet_wrap(~Region, ncol = 9) +
  
  # Custom shapes and colors for sex
  scale_shape_manual(values = c("Female" = 16, "Male" = 17, "Total" = 8)) + 
  scale_color_manual(values = c(
    "Female"    = "#6D383E",
    "Male"      = "#416E7D",
    "Total"   = "#6F6534"
  )) +
  
  ylim(0, 1) +
  # Combine the shape and color legends
  guides(color = guide_legend(override.aes = list(shape = c(16, 17, 8))),  # Override shapes in the color legend
         shape = guide_none()) +  # Hide the separate shape legend
  
  annotate("segment", x = Inf, xend = Inf, y = -Inf, yend = Inf, color = "black", linetype = "dashed") +
  
  # Modify axis labels and legend titles
  labs(x = "Year", y = expression('Harvest/km'^2*''), shape = "", color = "") +
  
  # Theme for consistent appearance
  theme_classic() +
  theme(
    text = element_text(size = 16),           # Increase the base text size
    axis.title = element_text(size = 16),     # Increase axis title size
    axis.text = element_text(size = 16),      # Increase axis text (ticks) size
    strip.text = element_text(size = 16),     # Increase facet label text size
    legend.title = element_text(size = 16, face = "bold"),   # Customize legend title text
    legend.text = element_text(size = 16),    # Customize legend item text
    legend.position = "top",                # Position the legend
    legend.key = element_rect(fill = "white"),  # Customize legend key appearance
    legend.key.size = unit(1.5, "lines") ,     # Adjust size of legend keys (shapes)
    axis.title.x = element_blank(),    # Remove x-axis title
    axis.text.x = element_blank()       # Remove x-axis text
  )

# Abundance Plot ----
abundance_plot <- ggplot(abundance_df_overall, aes(y = median_value / area_sq_km, x = year, shape = sex, color = sex)) +
  geom_point(size = 4, position = position_dodge(width = 0.5)) +
  geom_errorbar(aes(ymin = lower_ci / area_sq_km, ymax = upper_ci / area_sq_km), 
                width = 0.3, position = position_dodge(width = 0.5), show_guide = F) +
  facet_wrap(~Region, ncol = 9) +  # Use custom labeller
  
  # Custom shapes and colors for sex
  scale_shape_manual(values = c("Female" = 16, "Male" = 17, "Total" = 8)) + 
  scale_color_manual(values = c(
    "Female"    = "#6D383E",
    "Male"      = "#416E7D",
    "Total"   = "#6F6534"
  )) +
  
  annotate("segment", x = Inf, xend = Inf, y = -Inf, yend = Inf, color = "black", linetype = "dashed") +
  scale_x_continuous(breaks = c(1, 2, 3, 4, 5), labels = c("2019", "2020", "2021", "2022", "2023")) +
  
  labs(x = "Year", y = expression('Abundance/km'^2*''), shape = "", color = "") +  # Modify legend title
  
  # Combine the shape and color legends
  guides(color = guide_none(),  # Override shapes in the color legend
         shape = guide_none()) +  # Hide the separate shape legend
  ylim(0, 8) +
  theme_classic() +
  theme(
    text = element_text(size = 16),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 16),
    strip.text = element_blank(),      # Remove WMU labels
    legend.title = element_text(size = 16, face = "bold"),
    legend.text = element_text(size = 16),
    legend.position = "top",
    legend.key = element_rect(fill = "white"),  # Customize legend key appearance
    legend.key.size = unit(1.5, "lines"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )


### harvest/abundance ----
plot2 <- harvest_den_plot/abundance_plot

#### Save harvest/abundance ----
ggsave("Manuscript/Simple_abun-harv_plot.png", plot = plot2, width = 15, height = 12, dpi = 700)

## KF Survival Plot ----
kf_survival_plot <- ggplot(kf_survival_df, aes(y = median_value, x = year, 
                                               shape = sex, color = interaction(age_class, sex, sep = " "))) +
  geom_point(size = 4, position = position_dodge(width = 0.5)) +
  geom_errorbar(aes(ymin = lower_ci, ymax = upper_ci), 
                width = 0.3, 
                position = position_dodge(width = 0.5), show_guide = F) + 
  facet_wrap(~group, ncol = 3) +
  # Match colors for females as per the harvest plot
  scale_color_manual(values = c(
    "Adult Female"    = "#71250F",  # Color for adult female
    "Juvenile Female" = "#EB781B"  # Color for juvenile female
  )) +
  
  # Set shapes for age class (as you did before)
  scale_shape_manual(values = c("Female" = 16)) +
  annotate("segment", x = Inf, xend = Inf, y = -Inf, yend = Inf, color = "black", linetype = "dashed") +
  # Remove shape legend using guides()
  guides(shape = "none") +
  
  # Correct the x-axis for continuous year values
  scale_x_continuous(breaks = c(1, 2, 3, 4), labels = c("2020", "2021", "2022", "2023")) +
  
  # Y-axis limits
  scale_y_continuous(limits = c(0, 1)) +
  
  # Labels for axes and legend
  labs(x = "WMU group", y = "Survival probability", shape = "", color = "") +
  
  # Keep consistent theme
  theme_classic() +
  theme(
    text = element_text(size = 16),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 16),
    legend.title = element_text(size = 16, face = "bold"),
    legend.text = element_text(size = 16),
    legend.position = "top",
    legend.key = element_rect(fill = "white"),  # Customize legend key appearance
    legend.key.size = unit(1.5, "lines"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )


## Recruitment Plot ----
rec <- rec %>% 
  dplyr::select(!area_sq_km) %>% 
  left_join(wmu_areas, by = c("Region" = "WMU_ID")) %>% 
  mutate(density_value = median_value / area_sq_km)
recruitment_plot <- ggplot(rec, aes(y = density_value, x = year, color = as.factor(year))) +
  geom_point(size = 3, position = position_dodge(width = 0.5)) +
  geom_errorbar(aes(ymin = lower_ci / area_sq_km, ymax = upper_ci / area_sq_km), width = 0.3,
                position = position_dodge(width = 0.5), show_guide = F) +
  facet_wrap(~Region, ncol = 3) +
  
  # Y-axis with nice breaks
  scale_y_continuous(breaks = scales::pretty_breaks(n = 4)) +
  
  # Correct the x-axis for continuous year values
  scale_x_continuous(breaks = c(1, 2, 3, 4), labels = c("2020", "2021", "2022", "2023")) +
  
  # Axis labels and legend title
  labs(x = "Year", y = expression('Recruitment/km'^2*''), color = "") +
  
  # Add vertical lines at the right side of each facet
  annotate("segment", x = Inf, xend = Inf, y = -Inf, yend = Inf, color = "black", linetype = "dashed") +
  
  # Custom colors for years
  scale_color_manual(values = c("1" = "orange",
                                "2" = "#043A50",
                                "3" = "#C74442",
                                "4" = "plum4"), 
                     labels = c("2020", "2021", "2022", "2023")) +
  
  # Consistent theme
  theme_classic() +
  theme(
    text = element_text(size = 16),           # Increase base text size
    axis.title = element_text(size = 16),     # Increase axis title size
    axis.text = element_text(size = 16),      # Increase axis text (ticks) size
    strip.text = element_text(size = 16),     # Increase facet label text size
    legend.title = element_text(size = 16, face = "bold"),   # Customize legend title text
    legend.text = element_text(size = 16),    # Customize legend item text
    legend.position = "top",                  # Position the legend at the top
    legend.key = element_rect(fill = "white"),  # Customize legend key appearance
    legend.key.size = unit(1.5, "lines"),     # Adjust size of legend keys (shapes)
    axis.text.x = element_text(angle = 45, hjust = 1)  # Rotate x-axis labels
  )



#### Save the fem surv and rec ----
ggsave("Manuscript/V_kf_survival_plot.png", plot = kf_survival_plot,  width = 10, height = 6, dpi = 700)
ggsave("Manuscript/V_rec_plot.png", plot = recruitment_plot, width = 10, height = 8, dpi = 700)


