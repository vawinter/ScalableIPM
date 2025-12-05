# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
############## Research Integrated Population Model (R_IPM): ##################X
#                     #---# PhD Dissertation: R_IPM #---#
#        Creating a Bayesian IPM to inform turkey management in PA
###                       *** Real data run ***                             ###X
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #X
###                                                                         ###X
#                             Data visualization
###                                                                         ###X
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
library(ggh4x)

#selected_dir <- "Datavis/"
selected_dir <- "SubmissionMaterial/MajorRevisions/PubFigs/"

type = "R24"
# Load data 
kf_survival_df <- readRDS(paste0("Data/Output/", type, "_kf-survival_summary.rds"))
drm_harvest_df <- readRDS(paste0("Data/Output/", type, "_harvest_summary.rds"))
comb_survival_df <- readRDS(paste0("Data/Output/", type, "_comb-survival_summary.rds"))
abundance_df <- readRDS(paste0("Data/Output/", type, "_abundance_summary.rds"))
drm_survival_df <- readRDS(paste0("Data/Output/", type, "_drm_survival_summary.rds"))
ppb <- readRDS(paste0("Data/Output/", type, "_ppb_summary.rds"))
hwb <- readRDS(paste0("Data/Output/", type, "_hwb_summary.rds"))
rec <- readRDS(paste0("Data/Output/", type, "_rec_summary.rds"))

# Harvest
# Load in WMU areas
wmu_areas <- readRDS("Data/wmu_areas_km.rds") %>% 
  mutate(WMU_ID = paste("WMU", WMU_ID, sep = " "))

# Format harvest data
dat <- read.csv("Data/FallSprHarvData.csv", header=T)
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
  mutate(WMU.Group = paste("WMU", WMU.Group, sep = " ")) %>% 
  dplyr:: select(year, WMU.Group, wmu_count, age) %>% 
  as.data.frame()

# Spring Juvenile Jakes
wmu_spring_juv <- spring %>% 
  group_by(age, year, WMU.Group) %>% 
  reframe(wmu_count = sum(count)) %>% 
  filter(age == "jake",    
         WMU.Group %in% c("2D", "3D", "4D")) %>% 
  mutate(WMU.Group = paste("WMU", WMU.Group, sep = " ")) %>% 
  dplyr::select(year, WMU.Group, wmu_count, age) %>% 
  as.data.frame()

# Fall Adult Hens
wmu_fall_adult <- fall %>% 
  filter(!is.na(count)) %>% 
  group_by(age, year, WMU.Group) %>% 
  reframe(wmu_count = sum(count))%>% 
  filter(age == "hen",
         WMU.Group %in% c("2D", "3D", "4D")) %>%  
  mutate(WMU.Group = paste("WMU", WMU.Group, sep = " ")) %>% 
  dplyr::select(year, WMU.Group, wmu_count, age) %>% 
  as.data.frame()

# Fall Juvenile (Non-Hens)
wmu_fall_juv <- fall %>% 
  filter(!is.na(count)) %>% 
  group_by(age, year, WMU.Group) %>% 
  reframe(wmu_count = sum(count))%>% 
  filter(!age == "hen",
         WMU.Group %in% c("2D", "3D", "4D")) %>% 
  mutate(WMU.Group = paste("WMU", WMU.Group, sep = " ")) %>% 
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
  filter(!year %in% c(2019), age != "ALL") %>% 
  filter(sex == "Female")


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

wmu_col <- c("WMU 2D" = "#0f4c5d", 
         "WMU 3D" = "#043A70", 
         "WMU 4D" = "goldenrod")

# Filter too remove 5C survival estimate
combined_survival_df <- comb_survival_df %>% filter(wmu != "WMU 5C")

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
## DRM Survival Plot ----
#combined_survival_df <- combined_survival_df %>% filter(wmu == "3D")
survival_plot <- ggplot(combined_survival_df, aes(y = median_value, x = year, 
                                             shape = sex, color = interaction(age_class, sex, sep = " "))) +
  geom_point(size = 4, position = position_dodge(width = 0.5)) +
  geom_errorbar(aes(ymin = lower_ci, ymax = upper_ci), 
                width = 0.3, 
                position = position_dodge(width = 0.5), show.legend = F) +
  facet_wrap2(~wmu, ncol = 3,
              strip = strip_themed(
                background_x = elem_list_rect(fill = c("WMU 2D" = "#0f4c5d",
                                                       "WMU 3D" = "#043A70",
                                                       "WMU 4D" = "goldenrod"))
              )
  ) +
  
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
  
  # Legend
  guides(color = guide_legend(override.aes = list(shape = c(16, 16, 17, 17))),  # Override shapes in color legend
         shape = guide_none()) +  # Hide the separate shape legend
  
  # X-axis adjustments
  scale_x_continuous(breaks = c(1, 2, 3, 4, 5), labels = c("2020", "2021", "2022", "2023", "2024")) +
  
  # Y-axis limits for survival probabilities
  scale_y_continuous(limits = c(0, 1)) +  
  
  # Update labels to reflect the shape and color scheme
  labs(x = "Year", y = "Survival probability", shape = "", color = "") +  
  
  # Keep theme consistent with the previous plots
  theme_classic() +
  theme(
    text = element_text(size = 18),           # Increase the base text size
    axis.title = element_text(size = 18),     # Increase axis title size
    axis.text.x = element_blank(),      # Increase axis text (ticks) size
    strip.text = element_text(size = 18, face = "bold", color = "white"),     # Increase facet label text size
    legend.title = element_text(size = 18, face = "bold"),   # Customize legend title text
    legend.position = "top",                # Position the legend
    legend.key = element_rect(fill = "white"),  # Customize legend key appearance
    legend.key.size = unit(1.5, "lines"),      # Adjust size of legend keys (shapes)
    axis.title.x = element_blank()    # Remove x-axis title
  )
## Harvest rate Plot ----
harvest_plot1 <- drm_harvest_df %>% 
  filter(sex == "Male") %>% 
  ggplot(aes(y = median_value, x = year, 
             shape = sex, color = interaction(age_class, sex, sep = " "))) +
  geom_point(size = 4, position = position_dodge(width = 0.5)) +
  geom_errorbar(aes(ymin = lower_ci, ymax = upper_ci), 
                width = 0.3, 
                position = position_dodge(width = 0.5), show.legend = F) +
  facet_wrap(~wmu, ncol = 3) +
  
  # Assign shapes to sex: female and male
  scale_shape_manual(values = c("Male" = 17)) +
  
  # Assign colors to age_class and sex combinations
  scale_color_manual(values = c(
    "Adult Male"      = "#043A50",
    "Juvenile Male"   = "#5C7391"
  )) +
  
  scale_x_continuous(breaks = c(1, 2, 3, 4, 5), labels = c("2020", "2021", "2022", "2023", "2024")) +  # Change x-axis labels
  scale_y_continuous(limits = c(0, 0.4)) +  # Set y-axis limits
  annotate("segment", x = Inf, xend = Inf, y = -Inf, yend = Inf, color = "black", linetype = "dashed") +
  labs(x = "", y = "Harvest rate", shape = "", color = "") +  # Modify legend titles
  
  guides(color = guide_none(),  # Override shapes in color legend
         shape = guide_none()) +  # Hide the separate shape legend
  
  theme_classic() +
  theme(
    text = element_text(size = 18),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 18),
    strip.text = element_blank(),      # Remove WMU labels
    legend.title = element_text(size = 18, face = "bold"),
    legend.text = element_text(size = 20),
    legend.position = "none",
    legend.key = element_rect(fill = "gray"),  # Customize legend key appearance
    legend.key.size = unit(1.5, "lines"),
    axis.text.x = element_blank()
  )
harvest_plot2 <- drm_harvest_df %>% 
  filter(sex == "Female") %>% 
  ggplot( aes(y = median_value, x = year, 
              shape = sex, color = interaction(age_class, sex, sep = " "))) +
  geom_point(size = 4, position = position_dodge(width = 0.5)) +
  geom_errorbar(aes(ymin = lower_ci, ymax = upper_ci), 
                width = 0.3, 
                position = position_dodge(width = 0.5), show.legend = F) +
  facet_wrap2(~wmu, ncol = 3,
              strip = strip_themed(
                background_x = elem_list_rect(fill = c("WMU 2D" = "#0f9c9d",
                                                       "WMU 3D" = "#044A70",
                                                       "WMU 4D" = "goldenrod"))
              )
  ) +
  
  # Assign shapes to sex: female and male
  scale_shape_manual(values = c("Female" = 16)) +
  
  # Assign colors to age_class and sex combinations
  scale_color_manual(values = c(
    "Adult Female"    = "#7A370A",
    "Juvenile Female" = "#EB781B"
  )) +
  
  scale_x_continuous(breaks = c(1, 2, 3, 4, 5), labels = c("2020", "2021", "2022", "2023", "2024")) +  # Change x-axis labels
  scale_y_continuous(limits = c(0, 0.03)) +  # Set y-axis limits
  annotate("segment", x = Inf, xend = Inf, y = -Inf, yend = Inf, color = "black", linetype = "dashed") +
  labs(x = "Year", y = "Harvest rate", shape = "", color = "") +  # Modify legend titles
  
  guides(color = guide_none(),  # Override shapes in color legend
         shape = guide_none()) +  # Hide the separate shape legend
  
  theme_classic() +
  theme(
    text = element_text(size = 18),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 18),
    strip.text = element_blank(),      # Remove WMU labels
    legend.title = element_text(size = 18, face = "bold", color = "white"),
    legend.text = element_text(size = 20),
    legend.position = "none",
    legend.key = element_rect(fill = "white"),  # Customize legend key appearance
    legend.key.size = unit(1.5, "lines"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

## harv/den str ----
plot <- survival_plot/harvest_plot1/harvest_plot2

### Save the harv/den ----
ggsave(paste0(selected_dir, type, "_IPM_plot.png"), plot = plot, width = 11, height = 13, dpi = 700)
ggsave(paste0(selected_dir, "Fig2.pdf"), plot = plot, device = "pdf",  bg="transparent",
       width = 11, height = 13, dpi = 700)


## Abundance Plot ----
#abundance_df_overall <- abundance_df_overall %>% filter(wmu == "3D")
abundance_plot <- ggplot(abundance_df_overall, aes(y = median_value/area_sq_km, x = year, shape = sex, color = sex)) +
  geom_point(size = 4, position = position_dodge(width = 0.5)) +
  geom_errorbar(aes(ymin = lower_ci/area_sq_km, ymax = upper_ci/area_sq_km), 
                width = 0.3, position = position_dodge(width = 0.5), show.legend = F) +
  facet_wrap2(~wmu, ncol = 3,
              strip = strip_themed(
                background_x = elem_list_rect(fill = c("WMU 2D" = "#0f4c5d",
                                                       "WMU 3D" = "#043A70",
                                                       "WMU 4D" = "goldenrod")))
              ) +
  
  # Custom shapes and colors for sex
  scale_shape_manual(values = c("Female" = 16, "Male" = 17, "Total" = 8)) + 
  scale_color_manual(values = c(
    "Female"    = "#6D383E",
    "Male"      = "#416E7D",
    "Total"   = "#6F6534"
  )) +

  scale_x_continuous(breaks = c(1, 2, 3, 4, 5), labels = c("2020", "2021", "2022", "2023", "2024")) +
 ylim(0, 10) +
  labs(x = "Year", y = expression('Abundance/km'^2*''), shape = "", color = "") +  # Modify legend title
  
 # Combine the shape and color legends
  guides(color = guide_legend(override.aes = list(shape = c(16, 17, 8))),  # Override shapes in the color legend
         shape = guide_none()) +  # Hide the separate shape legend

  annotate("segment", x = Inf, xend = Inf, y = -Inf, yend = Inf, color = "black", linetype = "dashed") +
  theme_classic() +
  theme(
    text = element_text(size = 20),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 18),
    strip.text = element_text(size = 18,
                               face = "bold", color = "white"),     
    legend.title = element_text(size = 18, face = "bold"),
    legend.text = element_text(size = 18),
    legend.position = "top",
    legend.key = element_blank(),  # Customize legend key appearance
    legend.key.size = unit(1.5, "lines"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )


### Save abudnance ----
ggsave(paste0(selected_dir,type, "_abun_plot.png"), plot = abundance_plot, width = 15, height = 10, dpi = 700)
ggsave(paste0(selected_dir, "Fig3.pdf"), plot = abundance_plot, device = "pdf",  bg="transparent",
       width = 15, height = 10, dpi = 700)
# age class abundance:
## Abundance Plot ----
abundance_plot <- ggplot(abundance_df_overall, aes(y = median_value/area_sq_km, x = year, shape = sex, color = sex)) +
  geom_point(size = 4, position = position_dodge(width = 0.5)) +
  geom_errorbar(aes(ymin = lower_ci/area_sq_km, ymax = upper_ci/area_sq_km), 
                width = 0.3, position = position_dodge(width = 0.5), show.legend = F) +
  facet_wrap2(~wmu, ncol = 3,
              strip = strip_themed(
                background_x = elem_list_rect(fill = c("WMU 2D" = "#0f4c5d",
                                                       "WMU 3D" = "#043A70",
                                                       "WMU 4D" = "goldenrod")))
  ) +
  
  # Custom shapes and colors for sex
  scale_shape_manual(values = c("Female" = 16, "Male" = 17, "Total" = 8)) + 
  scale_color_manual(values = c(
    "Female"    = "#6D383E",
    "Male"      = "#416E7D",
    "Total"   = "#6F6534"
  )) +
  
  scale_x_continuous(breaks = c(1, 2, 3, 4, 5), labels = c("2020", "2021", "2022", "2023", "2024")) +
  ylim(0, 10) +
  labs(x = "Year", y = expression('Abundance/km'^2*''), shape = "", color = "") +  # Modify legend title
  
  # Combine the shape and color legends
  guides(color = guide_legend(override.aes = list(shape = c(16, 17, 8))),  # Override shapes in the color legend
         shape = guide_none()) +  # Hide the separate shape legend
  
  annotate("segment", x = Inf, xend = Inf, y = -Inf, yend = Inf, color = "black", linetype = "dashed") +
  theme_classic() +
  theme(
    text = element_text(size = 16),
    axis.title = element_text(size = 16),
    axis.text = element_blank(),
    strip.text = element_text(size = 18, face = "bold", color = "white"),      # Remove WMU labels
    legend.title = element_text(size = 16, face = "bold"),
    legend.text = element_text(size = 16),
    legend.position = "top",
    legend.key = element_rect(fill = "white"),  # Customize legend key appearance
    legend.key.size = unit(1.5, "lines"),
    axis.title.x = element_blank()#,    # Remove x-axis title
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
    strip.text = element_blank(),     # Increase facet label text size
    legend.title = element_text(size = 16, face = "bold"),   # Customize legend title text
    legend.text = element_text(size = 16),    # Customize legend item text
    legend.position = "none",                # Position the legend
    legend.key = element_rect(fill = "white"),  # Customize legend key appearance
    legend.key.size = unit(1.5, "lines") ,     # Adjust size of legend keys (shapes)
   axis.text.x = element_text(angle = 45, hjust = 1)
  )


## harv/den str ----
plot2 <- abundance_plot/harvest_den_plot

### Save the harv/den ----
ggsave(paste0(selected_dir, type,"_harv_plot.png"), plot = harvest_den_plot, width = 15, height = 9, dpi = 700)

## KF Survival Plot ----
kf_survival_plot <- kf_survival_df %>% 
  filter(sex == "Female") %>% 
  ggplot(aes(y = median_value, x = wmu, 
                            shape = sex, color = interaction(age_class, sex, sep = " "))) +
  geom_point(size = 4, position = position_dodge(width = 0.5)) +
  geom_errorbar(aes(ymin = lower_ci, ymax = upper_ci), 
                width = 0.3, 
                position = position_dodge(width = 0.5), show.legend = F) + 
  
  # Match colors for females as per the harvest plot
  scale_color_manual(values = c(
    "Adult Female"    = "#71250F",  # Color for adult female
    "Juvenile Female" = "#EB781B"  # Color for juvenile female
  #  "Juvenile Male"   = "#5C7391"
  )) +

  # Set shapes for age class (as you did before)
  scale_shape_manual(values = c("Female" = 16,"Male" = 17)) +
  
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
recruitment_plot <- ggplot(rec, aes(y = median_value/area_sq_km, x = year, shape = wmu, color = wmu)) +
  geom_point(size = 4, position = position_dodge(width = 0.5)) +
  geom_errorbar(aes(ymin = lower_ci/area_sq_km , ymax = upper_ci/area_sq_km ), width = 0.3,
                position = position_dodge(width = 0.5), show.legend = F) +
  
  # Custom shapes for WMUs
  scale_shape_manual(values = c("WMU 2D" = 18, "WMU 3D" = 20, "WMU 4D" = 15)) +
  
  # Custom colors for WMUs
  scale_color_manual(values = c("WMU 2D" = "#0f4c5d", 
                                "WMU 3D" = "#043A70", 
                                "WMU 4D" = "goldenrod")) +
  # Y-axis with nice breaks
 # scale_y_continuous(breaks = scales::pretty_breaks(n = 8)) +
  
  # X-axis labels for years
  scale_x_discrete(limits = factor(c(1, 2, 3, 4, 5)), labels = c("2020", "2021", "2022", "2023", "2024")) +
  
  # Axis labels and legend title
  labs(x = "Year", y = "Recruitment per sq. km", shape = "", color = "") +
  
  ylim(0, 3) +
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
ggsave(paste0(selected_dir, type, "_kf_survival_plot.png"), plot = kf_survival_plot,  width = 8, height = 6, dpi = 700)
ggsave(paste0(selected_dir, type, "_rec_plot.png"), plot = recruitment_plot, width = 8, height = 7, dpi = 700)

#PPB and HWB

ppb_plot <-  ggplot(ppb, aes(y = median_value, x = year, shape = wmu, color = wmu)) +
  geom_point(size = 4, position = position_dodge(width = 0.5)) +
  geom_errorbar(aes(ymin = lower_ci , ymax = upper_ci ), width = 0.3,
                position = position_dodge(width = 0.5), show.legend = T) +
  
  # Custom shapes for WMUs
  scale_shape_manual(values = c("WMU 2D" = 18, "WMU 3D" = 20, "WMU 4D" = 15)) +
  
  # Custom colors for WMUs
  scale_color_manual(values = c("WMU 2D" = "#0f4c5d", 
                                "WMU 3D" = "#043A70", 
                                "WMU 4D" = "goldenrod")) +
  # X-axis labels for years
  scale_x_discrete(limits = factor(c(1, 2, 3, 4, 5)), labels = c("2020", "2021", "2022", "2023", "2024")) +
  
  # Axis labels and legend title
  labs(x = "Year", y = "Poults per brood", shape = "", color = "") +
  ylim(2.5, 4) +
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
ggsave(paste0(selected_dir, type, "_ppb_plot.png"), plot = ppb_plot,  width = 8, height = 6, dpi = 700) 

hwb_plot <-  ggplot(hwb, aes(y = median_value, x = year, shape = wmu, color = wmu)) +
  geom_point(size = 4, position = position_dodge(width = 0.5)) +
  geom_errorbar(aes(ymin = lower_ci , ymax = upper_ci ), width = 0.3,
                position = position_dodge(width = 0.5), show.legend = T) +
  
  # Custom shapes for WMUs
  scale_shape_manual(values = c("WMU 2D" = 18, "WMU 3D" = 20, "WMU 4D" = 15)) +
  
  # Custom colors for WMUs
  scale_color_manual(values = c("WMU 2D" = "#0f4c5d", 
                                "WMU 3D" = "#043A70", 
                                "WMU 4D" = "goldenrod")) +
  # X-axis labels for years
  scale_x_discrete(limits = factor(c(1, 2, 3, 4, 5)), labels = c("2020", "2021", "2022", "2023", "2024")) +
  
  # Axis labels and legend title
  labs(x = "Year", y = "Hens with a brood", shape = "", color = "") +
  ylim(0.5, 1) +
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

ggsave(paste0(selected_dir, type, "_hwb_plot.png"), plot = hwb_plot,  width = 8, height = 6, dpi = 700) 


drm_harvest_df_newprior <- readRDS("Data/Output/TEST/20251118_R_IPM_run24_harvest_summary.rds")%>% 
  filter(sex == "Female")

drm_harvest_df23 <- readRDS("Data/Output/R__harvest_summary.rds") %>% 
  filter(sex == "Female")

plot1 <-  ggplot(drm_harvest_df_OGprior, aes(y = median_value, x = year, 
                                             shape = sex, color = interaction(age_class, sex, sep = " "))) +
  geom_point(size = 4, position = position_dodge(width = 0.5)) +
  geom_errorbar(aes(ymin = lower_ci, ymax = upper_ci), 
                width = 0.3, 
                position = position_dodge(width = 0.5), show.legend = F) +
  facet_wrap(~wmu, ncol = 3) +
  
  # Assign shapes to sex: female and male
  scale_shape_manual(values = c("Female" = 16)) +
  
  # Assign colors to age_class and sex combinations
  scale_color_manual(values = c(
    "Adult Female"    = "#7A370A",
    "Juvenile Female" = "#EB781B"
  )) +
  
  scale_x_continuous(breaks = c(1, 2, 3, 4, 5), labels = c("2020", "2021", "2022", "2023", "2024")) +  # Change x-axis labels
  scale_y_continuous(limits = c(0, 0.05)) +  # Set y-axis limits
  annotate("segment", x = Inf, xend = Inf, y = -Inf, yend = Inf, color = "black", linetype = "dashed") +
  labs(x = "Year", y = "Harvest rate", shape = "", color = "", title = "Original prior") +  # Modify legend titles
  
  guides(color = guide_legend(override.aes = list(shape = c(16, 16))),  # Override shapes in color legend
         shape = guide_none()) +  # Hide the separate shape legend
  
  theme_classic() +
  theme(
    text = element_text(size = 18),           # Increase the base text size
    axis.title = element_text(size = 18),     # Increase axis title size
    # axis.text.x = element_blank(),      # Increase axis text (ticks) size
    strip.text = element_text(size = 18),     # Increase facet label text size
    legend.title = element_text(size = 18, face = "bold"),   # Customize legend title text
    legend.position = "top",                # Position the legend
    legend.key = element_rect(fill = "white"),  # Customize legend key appearance
    legend.key.size = unit(1.5, "lines")     # Adjust size of legend keys (shapes)
  )

plot2 <-  ggplot(drm_harvest_df_newprior, aes(y = median_value, x = year, 
                                              shape = sex, color = interaction(age_class, sex, sep = " "))) +
  geom_point(size = 4, position = position_dodge(width = 0.5)) +
  geom_errorbar(aes(ymin = lower_ci, ymax = upper_ci), 
                width = 0.3, 
                position = position_dodge(width = 0.5), show.legend = F) +
  facet_wrap(~wmu, ncol = 3) +
  
  # Assign shapes to sex: female and male
  scale_shape_manual(values = c("Female" = 16)) +
  
  # Assign colors to age_class and sex combinations
  scale_color_manual(values = c(
    "Adult Female"    = "#7A370A",
    "Juvenile Female" = "#EB781B"
  )) +
  
  scale_x_continuous(breaks = c(1, 2, 3, 4, 5), labels = c("2020", "2021", "2022", "2023", "2024")) +  # Change x-axis labels
  scale_y_continuous(limits = c(0, 0.05)) +  # Set y-axis limits
  annotate("segment", x = Inf, xend = Inf, y = -Inf, yend = Inf, color = "black", linetype = "dashed") +
  labs(x = "Year", y = "Harvest rate", shape = "", color = "", , title = "New prior") +  # Modify legend titles
  
  guides(color = guide_legend(override.aes = list(shape = c(16, 16))),  # Override shapes in color legend
         shape = guide_none()) +  # Hide the separate shape legend
  
  theme_classic() +
  theme(
    text = element_text(size = 18),           # Increase the base text size
    axis.title = element_text(size = 18),     # Increase axis title size
    #axis.text.x = element_blank(),      # Increase axis text (ticks) size
    strip.text = element_text(size = 18),     # Increase facet label text size
    legend.title = element_text(size = 18, face = "bold"),   # Customize legend title text
    legend.position = "top",                # Position the legend
    legend.key = element_rect(fill = "white"),  # Customize legend key appearance
    legend.key.size = unit(1.5, "lines")      # Adjust size of legend keys (shapes)
  )
plot3 <-  ggplot(drm_harvest_df, aes(y = median_value, x = year, 
                                     shape = sex, color = interaction(age_class, sex, sep = " "))) +
  geom_point(size = 4, position = position_dodge(width = 0.5)) +
  geom_errorbar(aes(ymin = lower_ci, ymax = upper_ci), 
                width = 0.3, 
                position = position_dodge(width = 0.5), show.legend = F) +
  facet_wrap(~wmu, ncol = 3) +
  
  # Assign shapes to sex: female and male
  scale_shape_manual(values = c("Female" = 16)) +
  
  # Assign colors to age_class and sex combinations
  scale_color_manual(values = c(
    "Adult Female"    = "#7A370A",
    "Juvenile Female" = "#EB781B"
  )) +
  
  scale_x_continuous(breaks = c(1, 2, 3, 4), labels = c("2020", "2021", "2022", "2023")) +  # Change x-axis labels
  scale_y_continuous(limits = c(0, 0.05)) +  # Set y-axis limits
  annotate("segment", x = Inf, xend = Inf, y = -Inf, yend = Inf, color = "black", linetype = "dashed") +
  labs(x = "Year", y = "Harvest rate", shape = "", color = "", , title = "2023") +  # Modify legend titles
  
  guides(color = guide_legend(override.aes = list(shape = c(16, 16))),  # Override shapes in color legend
         shape = guide_none()) +  # Hide the separate shape legend
  
  theme_classic() +
  theme(
    text = element_text(size = 18),           # Increase the base text size
    axis.title = element_text(size = 18),     # Increase axis title size
    #axis.text.x = element_blank(),      # Increase axis text (ticks) size
    strip.text = element_text(size = 18),     # Increase facet label text size
    legend.title = element_text(size = 18, face = "bold"),   # Customize legend title text
    legend.position = "top",                # Position the legend
    legend.key = element_rect(fill = "white"),  # Customize legend key appearance
    legend.key.size = unit(1.5, "lines")      # Adjust size of legend keys (shapes)
  )
plot1
plot2
plot3
