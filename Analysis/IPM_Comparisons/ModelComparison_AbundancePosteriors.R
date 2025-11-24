# Comparing O-IPM and V-IPM posteriors

rm(list = ls())
gc()

# Libraries
library(ggridges)
library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)

# Load the data
operational_long <- readRDS("Data/Output/oipm_abundance.rds")
vague_long <- readRDS("Data/Output/vipm_abundance.rds")

wmu_areas <- readRDS("../../PSUTurkey/turkey_IPM/Data/wmu_km_areas_w.groups.rds") %>% 
  rename("Region" = "WMU_ID")

# Combine datasets
all_samples <- bind_rows(operational_long, vague_long)

# Extract demographic group, year, and region from parameter names
all_samples <- all_samples %>%
  mutate(Region = paste("Region", region, sep = " ")) %>% 
  left_join(wmu_areas) %>% 
  mutate(
    demographic_group = case_when(
      str_detect(parameter, "female.N.ad") ~ "Adult Female",
      str_detect(parameter, "female.N.juv") ~ "Juvenile Female",
      str_detect(parameter, "male.N.ad") ~ "Adult Male",
      str_detect(parameter, "male.N.juv") ~ "Juvenile Male",
      TRUE ~ NA_character_
    ),
    year = as.numeric(str_extract(parameter, "\\d+(?=\\.\\.)")) %>% 
      factor(labels = c("2020", "2021", "2022", "2023")),
    region = as.numeric(str_extract(parameter, "(?<=\\.\\.)(\\d+)(?=\\.)")),
    model = factor(model, levels = c("Operational", "Vague"))
  ) %>%
  # Order demographic group factor
  mutate(
    demographic_group = factor(demographic_group, 
                               levels = c("Adult Female", "Juvenile Female", 
                                          "Adult Male", "Juvenile Male")),
    density = abundance/area_sq_km
  )

# Create the density ridgeline plot for region 6
den_plot <- ggplot(all_samples, 
                       aes(x = density, 
                           y = demographic_group,
                           fill = model, 
                           color = model)) +
  # Use density ridges with stat="density" to better show the peaks
  geom_density_ridges(
    alpha = 0.7, 
    scale = 0.9, 
    rel_min_height = 0.01,
    bandwidth = 0.5,  # Adjust this to control smoothness
   # position = position_dodge(width = 0.3),
    quantile_lines = TRUE, 
    quantiles = c(0.25, 0.5, 0.75)
  ) +
  facet_wrap(~year) +
  # Set colors
  scale_fill_manual(values = c("Operational" = "#B34170", "Vague" = "#8856a7")) +
  scale_color_manual(values = c("Operational" = "#8f3359", "Vague" = "#6b4485")) +
  
  # Labels
  labs(x = expression('Abundance/km'^2*''),
       y = "") +
  
  # Theme
  theme_ridges() +
  theme(
    text = element_text(size = 14),
    legend.position = "top",
    strip.text = element_text(size = 16, face = "bold"),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 12),
    panel.spacing = unit(1, "lines"),
    legend.title = element_blank(),
    plot.title = element_text(size = 18, face = "bold"),
    plot.subtitle = element_text(size = 14)
  )

# Print the plot
print(region6_plot)

# Save the plot
ggsave("Dataviz/region6_posterior_density_comparison.png", region6_plot, width = 12, height = 8, dpi = 300)

# Create the density ridgeline plot faceted by region
den_plot <- ggplot(all_samples, 
                   aes(x = abundance, 
                       y = interaction(year, demographic_group),
                       fill = model, 
                       color = model)) +
  # Use density ridges with stat="density" to better show the peaks
  geom_density_ridges(
    alpha = 0.7, 
    scale = 0.9, 
    rel_min_height = 0.01,
    bandwidth = 0.5,  # Adjust this to control smoothness
    # position = position_dodge(width = 0.3),  # Removed as requested
    quantile_lines = TRUE, 
    quantiles = c(0.25, 0.5, 0.75)
  ) +

  # Set colors
  scale_fill_manual(values = c("Operational" = "#B34170", "Vague" = "#8856a7")) +
  scale_color_manual(values = c("Operational" = "#8f3359", "Vague" = "#6b4485")) +
  
  # Facet by region
  facet_wrap(~Region, scales = "free_x") +
  
  # Labels
  labs(x = "",
       y = "") +
  coord_cartesian(xlim = c(0, NA)) +  # Set lower limit to 0, let upper limit adjust per facet
  # Theme
  theme_ridges(grid = T,
               font_size = 12) +
  theme(
    text = element_text(size = 12),
    legend.position = "top",
    strip.text = element_text(size = 12),
   # axis.title = element_text(size = 16),
    axis.text = element_text(size = 12),
   # panel.spacing = unit(1, "lines"),
    legend.title = element_blank()
  )

# Print the plot
print(den_plot)

# Save the plot
ggsave("Dataviz/region_faceted_abundance.png", den_plot, width = 18, height = 16, dpi = 300)

# Save the plot
ggsave("Dataviz/all_regions_posterior_density_comparison.png", all_regions_plot, width = 16, height = 12, dpi = 300)