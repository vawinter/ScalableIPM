# Load necessary libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(dplyr)

# Assuming you've run your model and have results in 'survival_results'
# Extract storage parameter estimates
storage_params <- data.frame(survival_results[[1]])

# First, identify column names that start with "storage"
storage_cols <- grep("^storage", names(storage_params), value = TRUE)

# Create a function to extract indices from the parameter names
extract_indices <- function(param_name) {
  # Extract the indices from strings like "storage[1,1,1]"
  indices <- as.numeric(unlist(regmatches(
    param_name,
    gregexpr("\\d+", param_name)
  )))
  
  # Return WMU, age class, and month
  return(list(
    wmu = indices[1],
    age = indices[2],
    month = indices[3]
  ))
}

# Create a data frame for plotting
plot_data <- data.frame()

# Process each storage parameter
for (col in storage_cols) {
  # Extract indices
  indices <- extract_indices(col)
  
  # Skip if we couldn't extract three indices
  if (length(indices) < 3) next
  
  # Get the median survival value
  survival_value <- median(storage_params[[col]])
  survival_lower <- quantile(storage_params[[col]], 0.025)
  survival_upper <- quantile(storage_params[[col]], 0.975)
  
  # Add to our data frame
  plot_data <- rbind(plot_data, data.frame(
    wmu = indices$wmu,
    age_class = indices$age,
    month = indices$month,
    survival = survival_value,
    lower = survival_lower,
    upper = survival_upper,
    age_label = ifelse(indices$age == 1, "Adult", "Juvenile")
  ))
}


# Create the plot
ggplot(plot_data, aes(x = month, y = survival, group = age_label, color = age_label)) +
  geom_line(size = 1) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2, alpha = 0.7) +
  facet_wrap(~wmu) +
  scale_color_manual(values = c("Adult" = "darkred", "Juvenile" = "orange")) +
  labs(
    x = "Month",
    y = "Survival Probability",
    color = "Age Class"
  ) +
  ylim(0, 1) +
  xlim(1, 12) +
  theme_bw() +
  theme(legend.position = "top")


# Alternative: confidence ribbon instead of error bars
kf <- ggplot(plot_data, aes(x = month, y = survival, group = age_label, color = age_label, fill = age_label)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, color = NA) +
  geom_line(size = 1) +
  geom_point(size = 3) +
  facet_wrap(~wmu) +
  scale_color_manual(values = c("Adult" = "darkred", "Juvenile" = "orange")) +
  scale_fill_manual(values = c("Adult" = "darkred", "Juvenile" = "orange")) +
  labs(x = "Month", 
       y = "Survival Probability",
       color = "Age Class",
       fill = "Age Class") +
  theme_bw() +
  ylim(0, 1) +
  theme(legend.position = "top")
ggsave("../../Manuscripts/ScalableIPM/Datavis/20250326_kf_survival_plot_wmu_monthly_checkingIndex.png", plot = kf,  width = 8, height = 6, dpi = 700)
