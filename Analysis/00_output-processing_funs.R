# Function to extract variables and WMU and Year from variable names
extract_wmu_year <- function(variable_name) {
  # Extract the indices from variable names
  matches <- str_match(variable_name, "\\[(\\d+)(?:,\\s*(\\d+))?\\]")
  list(
    wmu = as.integer(matches[3]),
    year = ifelse(is.na(matches[2]), NA_integer_, as.integer(matches[2]))
  )
}

extract_wmu <- function(variable_name) {
  # Extract the single index from variable names
  matches <- str_match(variable_name, "\\[(\\d+)\\]")
  
  # Return the extracted index as wmu
  list(
    wmu = as.integer(matches[2])
  )
}


# Function to calculate summary statistics including confidence intervals
calculate_summary_stats <- function(df) {
  df %>%
    group_by(wmu, year, sex, age_class) %>%
    summarize(
      mean_value = mean(value, na.rm = TRUE),
      median_value = median(value, na.rm = TRUE),
      lower_ci = quantile(value, probs = 0.025, na.rm = TRUE),
      upper_ci = quantile(value, probs = 0.975, na.rm = TRUE),
      .groups = 'drop'
    )
}

# Function to calculate summary statistics including confidence intervals
calculate_summary_stats_wmu <- function(df) {
  df %>%
    group_by(wmu, sex, age_class) %>%
    summarize(
      mean_value = mean(value, na.rm = TRUE),
      median_value = median(value, na.rm = TRUE),
      lower_ci = quantile(value, probs = 0.025, na.rm = TRUE),
      upper_ci = quantile(value, probs = 0.975, na.rm = TRUE),
      .groups = 'drop'
    )
}

# Updated function to process each category and include WMU, Year, and summary statistics
process_category <- function(df, prefix, sex, age_class) {
  df %>% 
    dplyr::select(starts_with(prefix)) %>% 
    pivot_longer(cols = everything(), names_to = "variable", values_to = "value") %>%
    mutate(
      sex = sex,
      age_class = age_class,
      wmu = map_int(variable, ~ extract_wmu_year(.)$wmu),
      year = map_int(variable, ~ extract_wmu_year(.)$year)
    ) %>%
    calculate_summary_stats()
} 
process_category_wmu <- function(df, prefix, sex, age_class) {
    df %>% 
      dplyr::select(starts_with(prefix)) %>% 
      pivot_longer(cols = everything(), names_to = "variable", values_to = "value") %>%
      mutate(
        sex = sex,
        age_class = age_class,
        wmu = map_int(variable, ~ extract_wmu(.)$wmu)
      ) %>%
    calculate_summary_stats_wmu()
}



