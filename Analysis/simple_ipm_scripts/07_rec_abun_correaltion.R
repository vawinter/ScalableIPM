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
date <- "2021229"

# Load data 
abundance_df <- readRDS(paste0("Data/Output/Simple_", date, "_abundance_summary.rds"))
rec <- readRDS(paste0("Data/Output/Simple_", date, "_rec_summary.rds"))

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
  select(-year) %>% 
  filter(Year != 2020)


# Merge data together
combined_df <- abundance_df_overall %>% 
  mutate(Year = case_when(
    year == 1 ~ "2020",   
    year == 2 ~ "2021",   
    year == 3 ~ "2022",  
    year == 4 ~ "2023",
    TRUE ~ NA_character_  # Default to NA if year is not 1, 2, 3, or 4
  )) %>% left_join(recruitment, by = c("Year", "group")) %>% 
  select(-c("year"))%>% 
  filter(Year != 2020)

head(combined_df)
colnames(combined_df)

# Calculate the Pearson correlation between recruitment and abundance density
correlation <- cor(combined_df$rec_density, combined_df$abundance_density, method = "pearson")

# Print the result
print(correlation)

# Calculate the Spearman correlation
correlation_spearman <- cor(combined_df$rec_density, combined_df$abundance_density, method = "spearman")

# Calculate the Pearson and Spearman correlations across groups
correlations_groups <- combined_df %>%
  group_by(group) %>%
  summarise(
    # Pearson correlation
    pearson_correlation = cor(rec_density, abundance_density, method = "pearson"),
    
    # Spearman correlation
    spearman_correlation = cor(rec_density, abundance_density, method = "spearman"),
    
    # Perform the correlation test and extract p-value
    cor_test_result = list(cor.test(rec_density, abundance_density)),
    .groups = "drop"
  ) %>%
  mutate(
    # Extract p-value from the correlation test
    p_value = sapply(cor_test_result, function(x) x$p.value)
  )

# View the correlation results
print(correlations_groups)


# Scatter plot of recruitment vs. abundance_density
# Merge overall_trends and group_trends_correlations to add correlation values to the plot
plot_data <- combined_df %>%
  left_join(correlations_groups, by = "group")

# Create df_text with correlation values for each group (similar to your example)
df_text <- plot_data %>%
  group_by(group) %>%
  summarise(
    spear_corr = round(mean(spearman_correlation, na.rm = TRUE), 3)  # Calculate the Pearson's correlation for each group
  ) %>%
  mutate(
    x = 1.3,  # Position the text at the far right of each facet
    y = -Inf,  # Position the text at the bottom of each facet
    text = paste("Spearman's = ", spear_corr)  # Label the Pearson's correlation for each group
  )

# Scatter plot of recruitment vs. abundance_density with correlation annotations
corr.plot <- ggplot(plot_data, aes(x = rec_density, y = abundance_density, color = Year)) +
  geom_point(size = 3) +  # Scatter plot points
  geom_smooth(method = "lm", se = TRUE, color = "#6A87C1", fill = "lightblue") +  # Add regression line with CI
  labs(
    color = "",
    title = "",
    x = expression('Recruitment Density (birds/km'^2*')'),
    y = expression('Abundance Density (birds/km'^2*')')
  ) +
  theme_classic() +
  facet_wrap(~group) +  # Facet by group
  annotate("segment", x = Inf, xend = Inf, y = -Inf, yend = Inf, color = "black", linetype = "dashed") + # add dashed line
  scale_color_manual(values = c("2020" = "orange", "2021" = "#043A50", "2022" = "#C74442", "2023" = "plum4")) +  # Custom colors for each year
  geom_text(
    data = df_text,  # Use df_text with correlation annotations
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

# Print the plot
print(corr.plot)
#############################
# Calculate R² values for each group
correlations_groups_r2 <- combined_df %>%
  group_by(group) %>%
  summarise(
    # Calculate Pearson correlation and square it to get R²
    r_squared = summary(lm(rec_density ~ abundance_density))$r.squared,
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
summary(lm(combined_df$rec_density ~ combined_df$abundance_density))$r.squared
#  0.8090616
#----X


# Scatter plot of recruitment vs. abundance_density with correlation annotations
ggplot(combined_df, aes(x = rec_density, y = abundance_density, color = Year)) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", se = TRUE, color = "#6A87C1", fill = "lightblue") +  # Add regression line with CI
  labs(
    color = "",
    title = "Recruitment vs. Abundance Density by Year",
    x = expression('Recruitment Density (birds/km'^2*')'),
    y = expression('Abundance Density (birds/km'^2*')')
  ) +
  theme_classic() +
  scale_color_manual(values = c("2020" = "orange", "2021" = "#043A50", "2022" = "#C74442", "2023" = "plum4")) +
  facet_wrap(~group) +  # Facet by group
  theme(
    text = element_text(size = 16),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 16),
    strip.text = element_text(size = 16),
    legend.position = "top",
    legend.key = element_rect(fill = "white")
  )



# Perform correlation test (this also gives the p-value)
cor_test_result <- cor.test(combined_df$rec_density, combined_df$abundance_density)

# Print the result
print(cor_test_result)

# Save plot
ggsave("Manuscript/rec_abun_cor_wmu.png", plot = corr.plot_r2, width = 10, height = 8, dpi = 700)

### Realtive change ----
# Calculate relative change in abundance and recruitment starting from 2020
combined_df3 <- combined_df %>%
  arrange(group, Year) %>%  # Ensure the data is ordered by group and Year
  filter(sex == "Total") %>% 
  group_by(group) %>%
  mutate(
    # Calculate relative change in abundance density (as percentage)
    relative_change_abundance = (abundance_density / lag(abundance_density) - 1),
    
    # Calculate relative change in recruitment density (as percentage)
    relative_change_recruitment = (rec_density / lag(rec_density) - 1)
  ) %>%
  ungroup()

# View the resulting data frame with relative changes
head(combined_df3)

# Calculate Pearson and Spearman correlations for each group and year starting from 2020
correlations <- combined_df3 %>%
  group_by(group) %>%
  summarise(
    # Pearson correlation
    pearson_correlation = cor(relative_change_recruitment, relative_change_abundance, method = "pearson"),
    
    # Spearman correlation
    spearman_correlation = cor(relative_change_recruitment, relative_change_abundance, method = "spearman"),
    
    # Perform the correlation test and extract p-value
    cor_test_result = list(cor.test(relative_change_recruitment, relative_change_abundance)),
    .groups = "drop"
  ) %>%
  mutate(
    # Extract p-value from the correlation test
    p_value = sapply(cor_test_result, function(x) x$p.value)
  )

# Create a data frame to store the correlation annotations
df_text <- data.frame(
  group = unique(combined_df$group),  # Use unique group names from your plot data
  x = rep(1, length(unique(combined_df3$group))),  # Position text at the leftmost (or bottom) edge
  y = rep(1, length(unique(combined_df3$group))),  # Position text at the top edge (or bottom)
  text = paste("Pearson's = ", round(correlations$pearson_correlation, 3))  # Pearson's correlation for each group
)

# Scatter plot of relative change of recruitment vs. abundance_density
corr.plot <- ggplot(combined_df3, aes(x = relative_change_recruitment, y = relative_change_abundance, color = Year)) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", se = TRUE, color = "#6A87C1", fill = "lightblue") +  # Add regression line with CI
  labs(
    color = "",
    title = "Relative Change in Recruitment vs. Abundance Density",
    x = expression('Relative Change in Recruitment Density (%)'),
    y = expression('Relative Change in Abundance Density (%)')
  ) +
  theme_classic() +
  facet_wrap(~group) +  # Facet by group
  scale_color_manual(values = c("2020" = "orange", "2021" = "#043A50", "2022" = "#C74442", "2023" = "plum4")) +  # Custom colors for each year
  geom_text(
    data = df_text,  # Use df_text with correlation annotations
    aes(x = x, y = y, label = text),
    color = "black",  # Set text color to black
    fontface = "italic", hjust = 1, vjust = -0.9, size = 5, show.legend = FALSE, inherit.aes = FALSE
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

# Print the plot
print(corr.plot)
