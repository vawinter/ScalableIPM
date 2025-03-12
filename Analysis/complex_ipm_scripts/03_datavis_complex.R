# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#################### Integrated Population Model (IPM): #######################X
#                 #---# PhD Dissertation: Complex IPM #---#
#        Creating a Bayesian IPM to inform turkey management in PA
###                       *** Real data run ***                             ###X
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #X
###                                                                         ###X
#                             Data visualization
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
library(dplyr)
library(ggplot2)
library(tidyr)
library(knitr)
library(ggridges)
library(reshape2)
library(patchwork)

# Are we savig for datavis or manu?
selected_dir <- "Manuscript/"

# of data output
date <- "20241228"
type <- "Full_"

# Load data 
kf_survival_df <- readRDS(paste0("Data/Output/", type, date,"_kf-survival_summary.rds"))
drm_harvest_df <- readRDS(paste0("Data/Output/", type, date, "_harvest_summary.rds")) 
abundance_df <- readRDS(paste0("Data/Output/", type, date, "_abundance_summary.rds"))
drm_survival_df <- readRDS(paste0("Data/Output/", type, date, "_drm_survival_summary.rds"))
combined_survival_df <- readRDS(paste0("Data/Output/", type, date,"_comb-survival_summary.rds"))
ppb <- readRDS(paste0("Data/Output/", type, date, "_ppb_summary.rds"))
hwb <- readRDS(paste0("Data/Output/", type, date, "_hwb_summary.rds"))
rec <- readRDS(paste0("Data/Output/", type, date, "_rec_summary.rds"))

# Harvest
# Load in WMU areas
wmu_areas <- readRDS("../turkey_IPM/Data/wmu_areas_km.rds")

# Format harvest data
dat <- read.csv("Data/Banding_harv_data/FallSprHarvData_20240919.csv", header=T)
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
  left_join(wmu_areas, by = c("WMU.Group" = "WMU_ID")) %>% 
  filter(!year %in% c(2019, 2024)) 


# Structure abundance df
# Summarize data for total males, total females, and overall totals
harvest_df_totals <- Harvest %>%
  group_by(WMU.Group, year, sex) %>%
  summarise(
    median_value = sum(wmu_count),
    area_sq_km = first(area_sq_km),  # assuming area is the same across sex/age
    .groups = "drop"
  )

# Add 'Total' category by summing across sex and age_class 
harvest_df_overall <- harvest_df_totals %>%
  group_by(WMU.Group, year) %>%
  summarise(
    median_value = sum(median_value),
    area_sq_km = first(area_sq_km),  # assuming area is the same across sex/age
    .groups = "drop"
  ) %>%
  mutate(sex = "Total") %>%
  bind_rows(harvest_df_totals)  # Append the total rows to the original dataframe


# Color palettes  ----
colors <- c("#EB781B", "#CC5221", "#71250F", "#6F6534", "#365365")
colors2 <- c("#F6D9C0", "#C74442", "#71250E", "#6D383E", "#365365")

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
ext_colors <- c(
  "2D" = "#CC5221", # indianred4
  "3D" = "#043A50", # goldenrod
  "4D" = "#C74442" # sienna3
)

# Filter too remove 5C survival estimate
combined_survival_df <- combined_survival_df %>% filter(wmu != "5C")

# Structure abundance df
# Summarize data for total males, total females, and overall totals
abundance_df_totals <- abundance_df %>%
  group_by(wmu, year, sex) %>%
  summarise(
    median_value = sum(median_value),
    lower_ci = sum(lower_ci),
    upper_ci = sum(upper_ci),
    area_sq_km = first(area_sq_km),  # assuming area is the same across sex/age
    .groups = "drop"
  )

# Add 'Total' category by summing across sex and age_class 
abundance_df_overall <- abundance_df %>%
  group_by(wmu, year) %>%
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
## Harvest rate Plot ----
#drm_harvest_df <- drm_harvest_df %>% filter(wmu == "3D")
harvest_plot <- ggplot(drm_harvest_df, aes(y = median_value, x = year, 
                                           shape = sex, color = interaction(age_class, sex, sep = " "))) +
  geom_point(size = 4, position = position_dodge(width = 0.5)) +
  geom_errorbar(aes(ymin = lower_ci, ymax = upper_ci), 
                width = 0.3, 
                position = position_dodge(width = 0.5), show.legend = F) +
  facet_wrap(~wmu, ncol = 3) +
  
  # Assign shapes to sex: female and male
  scale_shape_manual(values = c("Female" = 16, "Male" = 17)) +
  
  # Assign colors to age_class and sex combinations
  scale_color_manual(values = c(
    "Adult Female"    = "#7A370A",
    "Juvenile Female" = "#EB781B",
    "Adult Male"      = "#043A50",
    "Juvenile Male"   = "#5C7391"
    )) +

  scale_x_continuous(breaks = c(1, 2, 3, 4), labels = c("2020", "2021", "2022", "2023")) +  # Change x-axis labels
  scale_y_continuous(limits = c(0, 0.43)) +  # Set y-axis limits
  annotate("segment", x = Inf, xend = Inf, y = -Inf, yend = Inf, color = "black", linetype = "dashed") +
  labs(x = "Year", y = "Harvest rate", shape = "", color = "") +  # Modify legend titles
  
  guides(color = guide_legend(override.aes = list(shape = c(16, 16, 17, 17))),  # Override shapes in color legend
         shape = guide_none()) +  # Hide the separate shape legend

  theme_classic() +
  theme(
    text = element_text(size = 18),           # Increase the base text size
    axis.title = element_text(size = 18),     # Increase axis title size
    axis.text = element_text(size = 18),      # Increase axis text (ticks) size
    strip.text = element_text(size = 18),     # Increase facet label text size
    legend.title = element_text(size = 18, face = "bold"),   # Customize legend title text
    legend.text = element_text(18),    # Customize legend item text
    legend.position = "top",                # Position the legend
    legend.key = element_rect(fill = "white"),  # Customize legend key appearance
    legend.key.size = unit(1.5, "lines"),      # Adjust size of legend keys (shapes)
    axis.title.x = element_blank(),    # Remove x-axis title
    axis.text.x = element_blank()       # Remove x-axis text
  )

## DRM Survival Plot ----
#combined_survival_df <- combined_survival_df %>% filter(wmu == "3D")
survival_plot <- ggplot(combined_survival_df, aes(y = median_value, x = year, 
                                             shape = sex, color = interaction(age_class, sex, sep = " "))) +
  geom_point(size = 4, position = position_dodge(width = 0.5)) +
  geom_errorbar(aes(ymin = lower_ci, ymax = upper_ci), 
                width = 0.3, 
                position = position_dodge(width = 0.5), show.legend = F) +
  facet_wrap(~wmu, ncol = 3) +
  
  # Use the same colors from the harvest plot for "adult male" and "juvenile male"
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
    text = element_text(size = 18),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 18),
    strip.text = element_blank(),      # Remove WMU labels
    legend.title = element_text(size = 18, face = "bold"),
    legend.text = element_text(size = 20),
    legend.position = "none",
    legend.key = element_rect(fill = "white"),  # Customize legend key appearance
    legend.key.size = unit(1.5, "lines"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

## harv/den str ----
plot <- harvest_plot/survival_plot

### Save the harv/den ----
ggsave(paste0(selected_dir, "Full_IPM_plot.png"), plot = plot, width = 15, height = 10, dpi = 700)
ggsave("Manuscript/complex_3D-surv.png", plot = plot, width = 10, height = 10, dpi = 700)


## Abundance Plot ----
#abundance_df_overall <- abundance_df_overall %>% filter(wmu == "3D")
abundance_plot <- ggplot(abundance_df_overall, aes(y = median_value / area_sq_km, x = year, shape = sex, color = sex)) +
  geom_point(size = 4, position = position_dodge(width = 0.5)) +
  geom_errorbar(aes(ymin = lower_ci/area_sq_km, ymax = upper_ci/area_sq_km), 
                width = 0.3, position = position_dodge(width = 0.5), show.legend = F) +
  facet_wrap(~wmu, ncol = 3) +  # Use custom labeller
  
  # Custom shapes and colors for sex
  scale_shape_manual(values = c("Female" = 16, "Male" = 17, "Total" = 8)) + 
  scale_color_manual(values = c(
    "Female"    = "#6D383E",
    "Male"      = "#416E7D",
    "Total"   = "#6F6534"
  )) +

  scale_x_continuous(breaks = c(1, 2, 3, 4), labels = c("2020", "2021", "2022", "2023")) +
 # scale_y_continuous(breaks = c(1:8)) +

  labs(x = "Year", y = expression('Abundance/km'^2*''), shape = "", color = "") +  # Modify legend title
  
 # Combine the shape and color legends
  guides(color = guide_legend(override.aes = list(shape = c(16, 17, 8))),  # Override shapes in the color legend
         shape = guide_none()) +  # Hide the separate shape legend
 #  ylim(0, 3) +
  # # Combine the shape and color legends
  # guides(color = guide_none(),  # Override shapes in the color legend
  #        shape = guide_none()) +  # Hide the separate shape legend
  annotate("segment", x = Inf, xend = Inf, y = -Inf, yend = Inf, color = "black", linetype = "dashed") +
  theme_classic() +
  theme(
    text = element_text(size = 20),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 18),
 #   strip.text = element_blank(),      # Remove WMU labels
    legend.title = element_text(size = 18, face = "bold"),
    legend.text = element_text(size = 18),
    legend.position = "top",
    legend.key = element_blank(),  # Customize legend key appearance
    legend.key.size = unit(1.5, "lines"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )


### Save abudnance ----
ggsave(paste0(selected_dir, "Full_abun_plot.png"), plot = abundance_plot, width = 15, height = 10, dpi = 700)
ggsave("Manuscript/complex_3D.png", plot = abundance_plot, width = 10, height = 10, dpi = 700)

# age class abundance:
## Abundance Plot ----
abundance_plot <- ggplot(abundance_df, aes(y = median_value / area_sq_km, x = year, shape = sex, color = sex)) +
  geom_point(size = 4, position = position_dodge(width = 0.5)) +
  geom_errorbar(aes(ymin = lower_ci/area_sq_km, ymax = upper_ci/area_sq_km), 
                width = 0.3, position = position_dodge(width = 0.5), show_guide = F) +
  facet_wrap(~wmu, ncol = 3) +  # Use custom labeller
  
  # Custom shapes and colors for sex
 # Use the same colors from the harvest plot for "adult male" and "juvenile male"
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
  scale_x_continuous(breaks = c(1, 2, 3, 4), labels = c("2020", "2021", "2022", "2023")) +
  # scale_y_continuous(breaks = c(1:8)) +
  
  labs(x = "Year", y = "Abundance per sq. km", shape = "", color = "") +  # Modify legend title
  
  # Combine the shape and color legends
  guides(color = guide_legend(override.aes = list(shape = c(16, 17, 8))),  # Override shapes in the color legend
         shape = guide_none()) +  # Hide the separate shape legend
  #  ylim(0, 3) +
  # # Combine the shape and color legends
  # guides(color = guide_none(),  # Override shapes in the color legend
  #        shape = guide_none()) +  # Hide the separate shape legend
  annotate("segment", x = Inf, xend = Inf, y = -Inf, yend = Inf, color = "black", linetype = "dashed") +
  theme_classic() +
  theme(
    text = element_text(size = 16),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 16),
    #   strip.text = element_blank(),      # Remove WMU labels
    legend.title = element_text(size = 16, face = "bold"),
    legend.text = element_text(size = 16),
    legend.position = "top",
    legend.key = element_rect(fill = "white"),  # Customize legend key appearance
    legend.key.size = unit(1.5, "lines"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )


## Harvest density ----
harvest_den_plot <- ggplot(harvest_df_overall, aes(y = median_value/area_sq_km, x = year, shape = sex, color = sex))  +
  
  # Plot points with dodge position to avoid overlap
  geom_point(size = 4, position = position_dodge(width = 0.5)) +
  
  # Optional: Add error bars if needed (if you have CI data)
 # geom_errorbar(aes(ymin = lower_ci/area_sq_km, ymax = upper_ci/area_sq_km), width = 0.3, position = position_dodge(width = 0.5)) +
  annotate("segment", x = Inf, xend = Inf, y = -Inf, yend = Inf, color = "black", linetype = "dashed") +
  # Facet by WMU
  facet_wrap(~WMU.Group, ncol = 3) +
  
  # Custom shapes and colors for sex
  scale_shape_manual(values = c("Female" = 16, "Male" = 17, "Total" = 8)) + 
  scale_color_manual(values = c(
    "Female"    = "#6D383E",
    "Male"      = "#416E7D",
    "Total"   = "#6F6534"
  )) +

  # Combine the shape and color legends
  guides(color = guide_legend(override.aes = list(shape = c(16, 17, 8))),  # Override shapes in the color legend
         shape = guide_none()) +  # Hide the separate shape legend
  
  # Modify axis labels and legend titles
  labs(x = "Year", y = "Harvest per sq. km", shape = "", color = "") +
  
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

## Abundance Plot ----
abundance_plot <- ggplot(abundance_df_overall, aes(y = median_value / area_sq_km, x = year, shape = sex, color = sex)) +
  geom_point(size = 4, position = position_dodge(width = 0.5)) +
  geom_errorbar(aes(ymin = lower_ci/area_sq_km, ymax = upper_ci/area_sq_km), 
                width = 0.3, position = position_dodge(width = 0.5), show.legend = F) +
  facet_wrap(~wmu, ncol = 3) +  # Use custom labeller
  
  # Custom shapes and colors for sex
  scale_shape_manual(values = c("Female" = 16, "Male" = 17, "Total" = 8)) + 
  scale_color_manual(values = c(
    "Female"    = "#6D383E",
    "Male"      = "#416E7D",
    "Total"   = "#6F6534"
  )) +
  
  scale_x_continuous(breaks = c(1, 2, 3, 4), labels = c("2020", "2021", "2022", "2023")) +
 # scale_y_continuous(breaks = c(1:8)) +
  labs(x = "Year", y = "Abundance per sq. km", shape = "", color = "") +  # Modify legend title

  # Combine the shape and color legends
  guides(color = guide_none(),  # Override shapes in the color legend
         shape = guide_none()) +  # Hide the separate shape legend
  
  annotate("segment", x = Inf, xend = Inf, y = -Inf, yend = Inf, color = "black", linetype = "dashed") +

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

## harv/den str ----
plot2 <- harvest_den_plot/abundance_plot

### Save the harv/den ----
ggsave(paste0(selected_dir, "Full_abun-harv_plot.png"), plot = plot2, width = 15, height = 9, dpi = 700)

## KF Survival Plot ----
kf_survival_plot <- ggplot(kf_survival_df, aes(y = median_value, x = wmu, 
                                               shape = sex, color = interaction(age_class, sex, sep = " "))) +
  geom_point(size = 4, position = position_dodge(width = 0.5)) +
  geom_errorbar(aes(ymin = lower_ci, ymax = upper_ci), 
                width = 0.3, 
                position = position_dodge(width = 0.5), show_guide = F) + 
  
  # Match colors for females as per the harvest plot
  scale_color_manual(values = c(
    "Adult Female"    = "#71250F",  # Color for adult female
    "Juvenile Female" = "#EB781B"  # Color for juvenile female
  )) +

  # Set shapes for age class (as you did before)
  scale_shape_manual(values = c("Female" = 16)) +
  
  # Remove shape legend using guides()
  guides(shape = "none") +

  # Y-axis limits
  scale_y_continuous(limits = c(0, 1)) +
  
  # Labels for axes and legend
  labs(x = "WMU", y = "Survival probability", shape = "", color = "") +
  
  # Keep consistent theme
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
    legend.key.size = unit(1.5, "lines")
  )

## Recruitment Plot ----
recruitment_plot <- ggplot(rec, aes(y = density_value, x = year, shape = wmu, color = wmu)) +
  geom_point(size = 4, position = position_dodge(width = 0.5)) +
  geom_errorbar(aes(ymin = lower_ci/area_sq_km , ymax = upper_ci/area_sq_km ), width = 0.3,
                position = position_dodge(width = 0.5), show_guide = F) +
  
  # Custom shapes for WMUs
  scale_shape_manual(values = c("2D" = 18, "3D" = 20, "4D" = 15)) +
  
  # Custom colors for WMUs
  scale_color_manual(values = c(  "2D" = "orange", # indianred4
                                  "3D" = "#043A50", # goldenrod
                                  "4D" = "#C74442")) +
  
  # Y-axis with nice breaks
  scale_y_continuous(breaks = scales::pretty_breaks(n = 8)) +
  
  # X-axis labels for years
  scale_x_discrete(limits = c(1, 2, 3, 4), labels = c("2020", "2021", "2022", "2023")) +
  
  # Axis labels and legend title
  labs(x = "Year", y = "Recruitment per sq. km", shape = "", color = "") +
  
  
  # Consistent theme
  theme_classic() +
  theme(
    text = element_text(size = 16),           # Increase base text size
    axis.title = element_text(size = 16),     # Increase axis title size
    axis.text = element_text(size = 16),      # Increase axis text (ticks) size
    strip.text = element_text(size = 16),     # Increase facet label text size
    legend.title = element_text(size = 16, face = "bold"),   # Customize legend title text
    legend.text = element_text(size = 16),    # Customize legend item text
    legend.position = "top",                # Position the legend
    legend.key = element_rect(fill = "white"),  # Customize legend key appearance
    legend.key.size = unit(1.5, "lines"),      # Adjust size of legend keys (shapes)
    axis.text.x = element_text(angle = 45, hjust = 1)
  )


### Save the survival/rec ----
ggsave("Manuscript/Full_kf_survival_plot.png", plot = kf_survival_plot,  width = 8, height = 6, dpi = 700)
ggsave("Manuscript/Full_rec_plot.png", plot = recruitment_plot, width = 8, height = 7, dpi = 700)


# Summary tables ----
## Summarize harvest rate ----
harvest_rate_summary <- drm_harvest_df %>%
  group_by(sex, wmu, age_class) %>%
  filter(age_class == "Adult") %>% 
  summarise(
    lower_ci = round(min(lower_ci), 3),
   # min_median = round(min(median_value), 3),  # Minimum harvest rate
    median = round(median(median_value), 3),   # Maximum harvest rate
    upper_ci = round(min(upper_ci), 3)
  )

# Display the summary table
print(harvest_rate_summary)

## Summarize survival rate ----
survival_rate_summary <- drm_survival_df %>%
  group_by(sex, wmu,  age_class) %>%
  summarise(
    lower_ci = round(min(lower_ci), 3),
   # min_median = round(min(median_value), 3),  # Minimum harvest rate
    median = round(median(median_value), 3),   # Maximum harvest rate
    upper_ci = round(min(upper_ci), 3)
  )

# Display the summary table
print(survival_rate_summary)

## Summarize density ----
density_rate_summary <- abundance_df %>%
  filter(wmu == "3D") %>% 
  group_by(sex, wmu, age_class, year) %>%
  reframe(
    lower_ci = round(min(lower_ci/area_sq_km, na.rm = TRUE), 3),
  #  min_median = round(min(median_value/area_sq_km, na.rm = TRUE), 3),  # Minimum harvest rate
    median = round(median_value/area_sq_km, 3),   # Maximum harvest rate
    upper_ci = round(min(upper_ci/area_sq_km, na.rm = TRUE), 3)
  ) 
# Display the summary table
print(density_rate_summary)

## Summarize female survival ----
survival_rate_summary <- kf_survival_df %>%
  group_by(age_class, wmu, sex) %>%
  summarise(
    lower_ci = round(min(lower_ci, na.rm = TRUE), 3),
    median = round(median(median_value, na.rm = TRUE), 3),  # Minimum harvest rate
   # max_median = round(max(median_value, na.rm = TRUE), 3),   # Maximum harvest rate
    upper_ci = round(min(upper_ci, na.rm = TRUE), 3)
  )


# Display the summary table
print(survival_rate_summary)

## Summarize female survival ----
rec_rate_summary <- rec %>%
  group_by(age_class, wmu) %>%
  summarise(
    min_harvest_rate = round(min(median_value, na.rm = TRUE), 3),  # Minimum harvest rate
    max_harvest_rate = round(max(median_value, na.rm = TRUE), 3)   # Maximum harvest rate
  )

# Display the summary table
print(rec_rate_summary)

## Summarize ppb ----
rec_rate_summary <- ppb %>%
  group_by(sex, wmu) %>%
  filter(year != 1) %>% 
  summarise(
    mean = round(mean(median_value), 3),
    min_harvest_rate = round(min(median_value, na.rm = TRUE), 3),  # Minimum harvest rate
    max_harvest_rate = round(max(median_value, na.rm = TRUE), 3)   # Maximum harvest rate
  )

# Display the summary table
print(rec_rate_summary)

# testing to find a ratio between males:females
# Summarize density ----
ratio_summary <- abundance_df %>%
  group_by(sex, wmu, year) %>%
  summarise(
    lower_ci = round(min(lower_ci, na.rm = TRUE), 3),
    min_median = round(min(median_value/area_sq_miles, na.rm = TRUE), 3),  # Minimum harvest rate
    max_median = round(max(median_value/area_sq_miles, na.rm = TRUE), 3),   # Maximum harvest rate
    upper_ci = round(min(upper_ci, na.rm = TRUE), 3)
  ) %>%
  pivot_wider(
    names_from = sex,
    values_from = c(lower_ci, min_median, max_median, upper_ci),
    names_sep = "_"
  ) %>%
  mutate(
    male_female_ratio_min_median = min_median_Male / min_median_Female,
    male_female_ratio_max_median = max_median_Male / max_median_Female,
    # Calculate mean of min_median and max_median
    male_female_ratio_mean_median = (min_median_Male + max_median_Male) / 2 / ((min_median_Female + max_median_Female) / 2)  # Mean ratio
  ) 

# View the result
print(density_rate_summary)


# Calculate the average ratio across WMUs
average_ratios <- density_rate_summary %>%
  group_by(wmu) %>% 
  summarise(
    avgmedian_ratio = mean(male_female_ratio_min_median, na.rm = TRUE),
  )

# > average_ratios
# # A tibble: 3 × 2
# wmu   avgmedian_ratio
# <chr>           <dbl>
# 1 2D             3.05
# 2 3D               1.23
# 3 4D               1.37

# Calculate the mean of the 'avgmedian_ratio' column
mean_avgmedian_ratio <- mean(average_ratios$avgmedian_ratio)

# Print the result
print(round(mean_avgmedian_ratio, 2))
# [1] 1.88277





#################### Age class -----
# age class abundance:
## Abundance Plot ----
abundance_plot <- ggplot(abundance_df, aes(y = median_value / area_sq_km, x = year, shape = sex, color = interaction(age_class, sex, sep = " ")))  +
  geom_point(size = 4, position = position_dodge(width = 0.5)) +
  geom_errorbar(aes(ymin = lower_ci/area_sq_km, ymax = upper_ci/area_sq_km), 
                width = 0.3, position = position_dodge(width = 0.5), show_guide = F) +
  facet_wrap(~wmu, ncol = 4) +  # Use custom labeller
  
  # Shapes for juvenile and adult (keeping the same shapes as before)
  scale_shape_manual(values = c("Female" = 16, "Male" = 17)) +
  
  # Custom shapes and colors for sex
  # Use the same colors from the harvest plot for "adult male" and "juvenile male"
  scale_color_manual(values = c(
    "Adult Female"    = "#7A370A",
    "Juvenile Female" = "#EB781B",
    "Adult Male"      = "#043A50",
    "Juvenile Male"   = "#5C7391"
  )) +
  annotate("segment", x = Inf, xend = Inf, y = -Inf, yend = Inf, color = "black", linetype = "dashed") +
  
  # Remove shape legend using guides()
  guides(shape = "none") +
  guides(color = "none") +
  scale_x_continuous(breaks = c(1, 2, 3, 4), labels = c("2020", "2021", "2022", "2023")) +
  # scale_y_continuous(breaks = c(1:8)) +
  
  labs(x = "Year", y = "Abundance per sq. km", shape = "", color = "") +  # Modify legend title
  # Hide the separate shape legend
  annotate("segment", x = Inf, xend = Inf, y = -Inf, yend = Inf, color = "black", linetype = "dashed") +
  theme_classic() +
  theme(
    text = element_text(size = 16),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 16),
    #   strip.text = element_blank(),      # Remove WMU labels
    legend.title = element_text(size = 16, face = "bold"),
    legend.text = element_text(size = 16),
    legend.position = "top",
    legend.key = element_rect(fill = "white"),  # Customize legend key appearance
    legend.key.size = unit(1.5, "lines"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
abundance_plot
