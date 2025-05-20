rm(list = ls())
gc()

# Load libraries
library(knitr)
library(dplyr)
library(kableExtra)
library(ggplot2)
library(dplyr)
library(tidyr)

# Posteriors
harvest_df <- readRDS("Data/Output/Operational_harvest_summary.rds") %>% 
  filter(sex == "Female")


#-------------------------
# Function to calculate percent shifts and related metrics
calculate_shifts <- function(harvest_df, prior_params) {
  # Extract prior parameters for different demographic groups
  prior_shape1_ad <- prior_params$adult$shape1
  prior_shape2_ad <- prior_params$adult$shape2
  prior_shape1_juv <- prior_params$juvenile$shape1
  prior_shape2_juv <- prior_params$juvenile$shape2
  
  # Calculate metrics
  result_df <- harvest_df %>%
    mutate(
      prior_shape1 = ifelse(age_class == "Adult", prior_shape1_ad, prior_shape1_juv),
      prior_shape2 = ifelse(age_class == "Adult", prior_shape2_ad, prior_shape2_juv),
      prior_mean = prior_shape1 / (prior_shape1 + prior_shape2),
      prior_var = (prior_shape1 * prior_shape2) / 
        ((prior_shape1 + prior_shape2)^2 * (prior_shape1 + prior_shape2 + 1)),
      post_sd = (upper_ci - lower_ci) / 4,
      post_var = post_sd^2,
      var_ratio = prior_var / post_var,
      percent_shift = ((mean_value - prior_mean) / prior_mean) * 100,
      abs_percent_shift = abs(percent_shift)
    )
  
  return(result_df)
}

# Define prior parameters
prior_params <- list(
  adult = list(shape1 = 3, shape2 = 30),     # Adjust to your actual priors
  juvenile = list(shape1 = 3, shape2 = 30)   # Adjust to your actual priors
)

# Calculate metrics for all parameters
shift_metrics <- calculate_shifts(harvest_df, prior_params)

# Create summary by region and demographic group
region_summary <- shift_metrics %>%
  group_by(wmu, sex, age_class) %>%
  summarize(
    mean_percent_shift = mean(percent_shift),
    median_percent_shift = median(percent_shift),
    min_percent_shift = min(percent_shift),
    max_percent_shift = max(percent_shift),
    mean_abs_shift = mean(abs_percent_shift),
    mean_var_ratio = mean(var_ratio),
    .groups = 'drop'
  ) %>%
  arrange(wmu, sex, age_class) %>% 
  mutate(wmu = as.character(wmu))

# Create overall summary by demographic group
overall_summary <- shift_metrics %>%
  group_by(sex, age_class) %>%
  summarize(
    mean_percent_shift = mean(percent_shift),
    median_percent_shift = median(percent_shift),
    min_percent_shift = min(percent_shift),
    max_percent_shift = max(percent_shift),
    mean_abs_shift = mean(abs_percent_shift),
    mean_var_ratio = mean(var_ratio),
    .groups = 'drop'
  ) %>%
  mutate(wmu = "Average") %>%
  select(wmu, everything()) %>%
  arrange(sex, age_class)

# Create complete table
complete_table <- bind_rows(region_summary, overall_summary) %>%
  arrange(sex, age_class, wmu)


# Format for presentation
formatted_table <- complete_table %>%
  mutate(across(where(is.numeric), ~round(., 2))) %>% 
  mutate(Region = as.character(paste("Region", wmu, sep = " ")),
         Sex_Age_Class = paste(sex, age_class)) %>%
  select(Sex_Age_Class, Region, mean_percent_shift ,min_percent_shift ,max_percent_shift, mean_abs_shift)

# Print formatted table
print(formatted_table)

table_reg3 <- kable(formatted_table,
                    booktabs = TRUE,   
                    format = "html", 
                    col.names = c("Sex & Age Class","Region","Percent shift (mean)","Percent shift (min)","Percent shift (max)", "Absolute shift (mean)")) %>%
  # Collapse rows for Sex/Age Class and WMU
  collapse_rows(columns = c(1, 2), latex_hline = "major", valign = "middle") %>%
  kable_styling(latex_options = c("hold_position", "repeat_header"), full_width = FALSE) 

# Export to CSV
save_kable(table_reg3, file = "Dataviz/Appx7Tb1.html")

# vague model ----
# Posteriors
harvest_df <- readRDS("../../PSUTurkey/turkey_IPM/Data/Output/Simple_vague_20250108_harvest_summary.rds") %>% 
  filter(sex == "Female")
library(ggplot2)
library(dplyr)
library(tidyr)
# Define prior parameters
prior_params <- list(
  adult = list(shape1 = 1, shape2 = 1),     # Adjust to your actual priors
  juvenile = list(shape1 = 1, shape2 = 1)   # Adjust to your actual priors
)

# Calculate metrics for all parameters
shift_metrics <- calculate_shifts(harvest_df, prior_params)

# Create summary by region and demographic group
region_summary <- shift_metrics %>%
  group_by(wmu, sex, age_class) %>%
  summarize(
    mean_percent_shift = mean(percent_shift),
    median_percent_shift = median(percent_shift),
    min_percent_shift = min(percent_shift),
    max_percent_shift = max(percent_shift),
    mean_abs_shift = mean(abs_percent_shift),
    mean_var_ratio = mean(var_ratio),
    .groups = 'drop'
  ) %>%
  arrange(wmu, sex, age_class) %>% 
  mutate(wmu = as.character(wmu))

# Create overall summary by demographic group
overall_summary <- shift_metrics %>%
  group_by(sex, age_class) %>%
  summarize(
    mean_percent_shift = mean(percent_shift),
    median_percent_shift = median(percent_shift),
    min_percent_shift = min(percent_shift),
    max_percent_shift = max(percent_shift),
    mean_abs_shift = mean(abs_percent_shift),
    mean_var_ratio = mean(var_ratio),
    .groups = 'drop'
  ) %>%
  mutate(wmu = "Average") %>%
  select(wmu, everything()) %>%
  arrange(sex, age_class)

# Create complete table
complete_table <- bind_rows(region_summary, overall_summary) %>%
  arrange(sex, age_class, wmu)


# Format for presentation
formatted_table <- complete_table %>%
  mutate(across(where(is.numeric), ~round(., 2))) %>% 
  mutate(Region = as.character(paste("Region", wmu, sep = " ")),
         Sex_Age_Class = paste(sex, age_class)) %>%
  select(Sex_Age_Class, Region, mean_percent_shift ,min_percent_shift ,max_percent_shift, mean_abs_shift)

# Print formatted table
print(formatted_table)

table_reg4 <- kable(formatted_table,
                    booktabs = TRUE,   
                    format = "html", 
                    col.names = c("Sex & Age Class","Region","Percent shift (mean)","Percent shift (min)","Percent shift (max)", "Absolute shift (mean)")) %>%
  # Collapse rows for Sex/Age Class and WMU
  collapse_rows(columns = c(1, 2), latex_hline = "major", valign = "middle") %>%
  kable_styling(latex_options = c("hold_position", "repeat_header"), full_width = FALSE) 
save_kable(table_reg4, file = "Dataviz/Appx7Tb3.html")
