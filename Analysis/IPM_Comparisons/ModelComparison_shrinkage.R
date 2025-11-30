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
library(dplyr)
library(ggplot2)
library(kableExtra)

# Load in WMU areas
wmu_areas1 <- readRDS("Data/wmu_areas_km.rds")
wmu_areas2 <- readRDS("Data/wmu_km_areas_regions.rds")
wmu_areas <- rbind(wmu_areas1, wmu_areas2)

# Load data 
abundance_df_1 <- readRDS(paste0("Data/Output/R24_abundance_summary.rds")) %>% 
  mutate(Model = "Research")
abundance_df_2 <- readRDS(paste0("Data/Output/O_24_abundance_summary.rds")) %>% 
  mutate(wmu = as.character(wmu))%>% 
  mutate(Model = "Operational")
abundance_df_3 <- readRDS(paste0("Data/Output/V24_abundance_summary.rds")) %>% 
  mutate(wmu = as.character(wmu)) %>% 
  select(-c(area_sq_km, density_value)) %>% 
  left_join(wmu_areas2, by = c("Region" = "WMU_ID")) %>% 
  mutate(density_value = median_value/area_sq_km)%>% 
  mutate(Model = "Vague")

# Structure abundance df
# First, filter for only the regions we want (3D and Region 6)
filtered_abundance_df <- abundance_df_1 %>%
  filter(wmu %in% c("WMU 3D"))

filtered_abundance_df2 <- abundance_df_2 %>%
  filter(Region %in% c("Region 6")) %>% 
  select(-wmu) %>% 
  rename(wmu = Region)

filtered_abundance_df3 <- abundance_df_3 %>%
  filter(Region %in% c("Region 6")) %>% 
  select(-wmu) %>% 
  rename(wmu = Region)

wmu_areas <- readRDS("Data/wmu_km_areas_regions.rds") %>% 
  rename("Region" = "WMU_ID")

# Combine all models data
all_abundance_df <- filtered_abundance_df %>%
  bind_rows(filtered_abundance_df2) %>%
  bind_rows(filtered_abundance_df3) %>%
  group_by(wmu, year, sex, Model)

# Calculate annual variation metrics
annual_variation <- all_abundance_df %>%
  # Group by model, wmu/region, sex, and age_class to calculate means across years
  group_by(Model, wmu, sex, age_class) %>%
  summarize(
    # Calculate the mean abundance across all years
    mean_abundance = mean(median_value),
    # Calculate average absolute deviation from the mean
    avg_abs_deviation = mean(abs(median_value - mean_abundance)),
    # Calculate coefficient of variation (standardized measure)
    cv = sd(median_value) / mean_abundance,
    # Calculate the percent deviation
    deviation = avg_abs_deviation / mean_abundance,
    # Number of years in the data
    n_years = n(),
    .groups = "drop"
  )

# Print the annual variation metrics
print(annual_variation)

# Create a visualization of the percent deviation metric
ggplot(annual_variation, 
       aes(x = wmu, y = deviation*100, fill = Model)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  facet_grid(sex ~ age_class) +
  labs(title = "Annual Variation in Turkey Abundance Estimates",
       subtitle = "Measured as average percent deviation from the mean",
       x = "Wildlife Management Unit",
       y = "Percent Deviation (%)") +
  theme_minimal() +
  theme(
    legend.position = "top",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

# Save the plot
ggsave("Dataviz/annual_variation_by_wmu.png", width = 12, height = 8, dpi = 300)

# Table of results for the paper
annual_variation_table <- annual_variation %>%
  select(Model, wmu,  sex, age_class, deviation, cv) %>%
  arrange(Model, wmu, sex, age_class) %>%
  mutate(
    percent_deviation = round(deviation, 1),
    cv = round(cv, 3)
  )

# Print the table
print(annual_variation_table)

# Save the table
write.csv(annual_variation_table, "Dataviz/annual_variation_metrics.csv", row.names = FALSE)
#----------------------------
library(dplyr)
library(ggplot2)

# Focus on just WMU 3D (Research) and Region 6 (Operational)
# These should be equivalent geographic areas
comparison_data <- bind_rows(
  filtered_abundance_df %>% select(Model, wmu, year, sex, age_class, median_value),
  filtered_abundance_df2 %>% select(Model, wmu, year, sex, age_class, median_value),
  filtered_abundance_df3 %>% select(Model, wmu, year, sex, age_class, median_value)
)

# Create a new summary combining age classes
combined_age_summary <- comparison_data %>%
  # First sum the abundances across age classes for each year
  group_by(Model, wmu, year, sex) %>%
  summarize(
    total_abundance = sum(median_value),
    .groups = "drop"
  ) %>%
  # Then calculate variation metrics for these combined abundances
  group_by(Model, sex) %>%
  mutate(
    mean_abundance = mean(total_abundance),
    abs_deviation = abs(total_abundance - mean_abundance),
    percent_deviation = 100 * abs_deviation / mean_abundance
  ) %>%
  group_by(Model, sex) %>%
  summarize(
    mean_abundance = first(mean_abundance),
    avg_abs_deviation = mean(abs_deviation),
    avg_percent_deviation = mean(percent_deviation),
    cv = sd(total_abundance) / mean_abundance,
    .groups = "drop"
  )

# Create a more readable comparison
combined_comparison <- combined_age_summary %>%
  select(Model, sex, avg_percent_deviation, cv) %>%
  pivot_wider(
    names_from = Model,
    values_from = c(avg_percent_deviation, cv)
  )

# Print the results
print(combined_comparison)

# Create visualization for combined age classes
combined_plot <- ggplot(combined_age_summary, 
                        aes(x = sex, y = avg_percent_deviation, fill = Model)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  labs(title = "Annual Variation in WMU 3D/Region 6",
       subtitle = "Total males and females (age classes combined)",
       x = "Sex",
       y = "Average Percent Deviation (%)") +
  theme_minimal() +
  theme(legend.position = "top")

# Print and save the plot
print(combined_plot)
ggsave("Dataviz/combined_age_stochasticity.png", combined_plot, width = 8, height = 6, dpi = 300)

# Calculate mean abundances and deviations for each model/demographic group
shrinkage_metrics <- comparison_data %>%
  group_by(Model, sex, age_class) %>%
  mutate(
    # Mean abundance across all years for this model/demographic
    mean_abundance = mean(median_value),
    # Absolute deviation from mean for each year
    abs_deviation = abs(median_value - mean_abundance),
    # Percent deviation from mean
    percent_deviation = 100 * abs_deviation / mean_abundance
  ) %>%
  # Summarize to get average deviation metrics
  group_by(Model, sex, age_class) %>%
  summarize(
    mean_abundance = first(mean_abundance),
    avg_abs_deviation = mean(abs_deviation),
    avg_percent_deviation = mean(percent_deviation),
    max_percent_deviation = max(percent_deviation),
    cv = sd(median_value) / mean_abundance,
    .groups = "drop"
  ) %>%
  # Add a column identifying if this is the reference region (WMU 3D/Region 6)
  mutate(region_type = "WMU 3D/Region 6")

# Create a summary table
shrinkage_summary <- shrinkage_metrics %>%
  select(Model, sex, age_class, avg_percent_deviation, cv) %>%
  arrange(sex, age_class, Model) %>%
  pivot_wider(
    names_from = Model,
    values_from = c(avg_percent_deviation, cv)
  ) %>%
  mutate(
    # Calculate difference between models
    percent_diff_R_vs_O = avg_percent_deviation_Research - avg_percent_deviation_Operational,
    percent_diff_O_vs_V = avg_percent_deviation_Operational - avg_percent_deviation_Vague,
    # Calculate relative difference (as percentage)
    relative_diff_R_vs_O = 100 * percent_diff_R_vs_O / avg_percent_deviation_Research,
    relative_diff_O_vs_V = 100 * percent_diff_O_vs_V / avg_percent_deviation_Operational
  )

# Print the summary table
print(shrinkage_summary)

# Create a visualization
shrinkage_plot <- ggplot(shrinkage_metrics, 
                         aes(x = interaction(sex, age_class), y = avg_percent_deviation, fill = Model)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  labs(title = "Annual Variation in WMU 3D/Region 6",
       subtitle = "Measured as average percent deviation from the mean abundance",
       x = "Demographic Group",
       y = "Average Percent Deviation (%)") +
  theme_minimal() +
  theme(
    legend.position = "top",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

# Print the plot
print(shrinkage_plot)

# Save outputs
write.csv(shrinkage_metrics, "Dataviz/shrinkage_metrics_WMU3D_Region6.csv", row.names = FALSE)
ggsave("Dataviz/shrinkage_comparison_WMU3D_Region6.png", shrinkage_plot, width = 10, height = 6, dpi = 300)