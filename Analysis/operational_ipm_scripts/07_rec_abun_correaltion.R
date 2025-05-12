# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#################### Integrated Population Model (IPM): #######################X
#                 #---# PhD Dissertation: Simple IPM #---#
#        Creating a Bayesian IPM to inform turkey management in PA
###                       *** Real data run ***                             ###X
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #X
###                                                                         ###X
#                 Correlation between recruitment and abundance
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #X
# 
# Created by: Veronica A. Winter
# Last edited: 01/01/2024
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

# of data output
date <- "20250326"

# Load data 
abundance_df <- readRDS(paste0("../../PSUTurkey/turkey_IPM/Data/Output/Full_", date, "_abundance_summary.rds"))
rec <- readRDS(paste0("../../PSUTurkey/turkey_IPM/Data/Output/Full_", date, "_rec_summary.rds"))

abundance_df$group <- abundance_df$wmu
rec$group <- rec$wmu

# Data structuring
# Summarize data for total males, total females, and overall totals
abundance_df_totals <- abundance_df %>%
  group_by(group, year, sex) %>%
  summarise(
    median_value = sum(median_value),
    lower_ci = sum(lower_ci),
    upper_ci = sum(upper_ci),
    area_sq_km = first(area_sq_km),  # assuming area is the same across sex/age
    .groups = "drop"
  )

# Add 'Total' category by summing across sex and age_class 
abundance_df_overall <- abundance_df %>%
  group_by(group, year) %>%
  summarise(
    median_value = sum(median_value),
    lower_ci = sum(lower_ci),
    upper_ci = sum(upper_ci),
    area_sq_km = first(area_sq_km),  # assuming area is the same across sex/age
    .groups = "drop"
  ) %>%
  mutate(sex = "Total") %>%
  bind_rows(abundance_df_totals) %>%   # Append the total rows to the original dataframe
  mutate(abundance_density = median_value / area_sq_km) %>% 
  filter(sex == "Total")

# rename columns in rec
recruitment <- rec %>% 
  rename(
    rec_density = density_value,  # Rename density_value back to rec_density
    rec_lower_ci = lower_ci,      # Rename lower_ci to rec_lower_ci
    rec_upper_ci = upper_ci       # Rename upper_ci to rec_upper_ci
  ) %>% 
  dplyr::select(rec_density, group, year, rec_lower_ci, rec_upper_ci) %>% 
  mutate(Year = case_when(
    year == 1 ~ "2020",   # If year is 1, assign "2020"
    year == 2 ~ "2021",   # If year is 2, assign "2021"
    year == 3 ~ "2022",   # If year is 3, assign "2022"
    year == 4 ~ "2023",   # If year is 4, assign "2023"
    TRUE ~ NA_character_  # Default to NA if year is not 1, 2, 3, or 4
  )) %>% 
  select(-year)



# Merge data together
combined_df <- abundance_df_overall %>% 
  mutate(Year = case_when(
    year == 1 ~ "2020",   
    year == 2 ~ "2021",   
    year == 3 ~ "2022",  
    year == 4 ~ "2023",
    TRUE ~ NA_character_  # Default to NA if year is not 1, 2, 3, or 4
  )) %>% left_join(recruitment, by = c("Year", "group")) %>% 
  select(-c("year"))

head(combined_df)
colnames(combined_df)

# Scatter plot of recruitment vs. abundance_density
# Merge overall_trends and group_trends_correlations to add correlation values to the plot
plot_data <- combined_df 

#############################
# Calculate R² values for each group
#############################
# Calculate R² values for each group
correlations_groups_r2 <- combined_df %>%
  group_by(group) %>%
  mutate(rec_density_prev = lag(rec_density)) %>%  # create lagged rec_density
  filter(!is.na(rec_density_prev)) %>% 
  summarise(
    # Calculate Pearson correlation and square it to get R²
    r_squared = summary(lm(rec_density_prev ~ abundance_density))$r.squared,
    .groups = "drop"
  )


# View the R² values for each group
print(correlations_groups_r2)

# Merge R² values with the plot data
plot_data_r2 <- plot_data %>%
  left_join(correlations_groups_r2, by = "group")

# Create df_text with R² values for each group
df_text_r2 <- plot_data_r2 %>%
  group_by(group) %>%
  summarise(
    r2 = round(mean(r_squared, na.rm = TRUE), 3)  # Calculate the R² value for each group
  ) %>%
  mutate(
    x = 1.3,  # Position the text at the far right of each facet
    y = -Inf,  # Position the text at the bottom of each facet
    text = paste("R² = ", r2)  # Label with the R² value for each group
  )

# Scatter plot of recruitment vs. abundance_density with R² annotations
corr.plot_r2 <- ggplot(plot_data_r2, aes(x = rec_density, y = abundance_density, color = Year)) +
  geom_point(size = 3) +  # Scatter plot points
  geom_smooth(method = "lm", se = TRUE, color = "#6A87C1", fill = "lightblue") +  # Add regression line with CI
  labs(
    color = "",
    title = "",
    x = expression('Recruitment/km'^2*''),
    y = expression('Abundance/km'^2*'')
  ) +
  theme_classic() +
  facet_wrap(~group) +  # Facet by group
  annotate("segment", x = Inf, xend = Inf, y = -Inf, yend = Inf, color = "black", linetype = "dashed") + # Add dashed line
  scale_color_manual(values = c("2020" = "orange", "2021" = "#043A50", "2022" = "#C74442", "2023" = "plum4")) +  # Custom colors for each year
  geom_text(
    data = df_text_r2,  # Use df_text with R² annotations
    aes(x = x, y = y, label = text),
    color = "black",  # Set text color to black directly here
    fontface = "italic", hjust = 1, vjust = -1, size = 5, show.legend = FALSE, inherit.aes = FALSE
  ) +
  theme(
    text = element_text(size = 16),           # Increase base text size
    axis.title = element_text(size = 16),     # Increase axis title size
    axis.text = element_text(size = 16),      # Increase axis text (ticks) size
    strip.text = element_text(size = 16),     # Increase facet label text size
    legend.title = element_text(size = 16, face = "bold"),   # Customize legend title text
    legend.text = element_text(size = 16),    # Customize legend item text
    legend.position = "top",                # Position the legend
    legend.key = element_rect(fill = "white")  # Customize legend key appearance
  )

# Print the plot with R² annotations
print(corr.plot_r2)
combined_df <- combined_df %>% 
  mutate(rec_density_prev = lag(rec_density)) %>%  # create lagged rec_density
  filter(!is.na(rec_density_prev))  
summary(lm(combined_df$rec_density_prev ~ combined_df$abundance_density))$r.squared
#  0.39 simple
# 0.54 full
#----X

# Save plot
ggsave("Dataviz/rec_abun_cor_full.png", plot = corr.plot_r2, width = 10, height = 8, dpi = 700)

