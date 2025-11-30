# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#################### Integrated Population Model (IPM): #######################X
#                 #---# PhD Dissertation: Simple IPM #---#
#        Creating a Bayesian IPM to inform turkey management in PA
###                       *** Real data run ***                             ###X
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #X
###                                                                         ###X
#                             Data visualization
###                           Mapping densities                             ###X
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #X
###############################################################################X

# Clear environment
rm(list=ls())
gc()

# Load necessary packages
library(sf)
library(tidyverse)
library(viridis)
library(lubridate)
library(rnaturalearth)
library(geosphere)
library(ggspatial)
library(maps)
library(showtext)

# Data read-in ----
# Add in states so we can include Pennsylvania boundaries
states <- st_as_sf(maps::map("state", plot = FALSE, fill = TRUE))

# Filter to Pennsylvania
pa <- states %>% filter(ID == "pennsylvania")

# Read in shapefile of WMUs
wmu_shapefile <- st_read("../../PSUTurkey/TurkeyProject/Data/PGC_BNDWildlifeManagementUnits2024.shp")

# Read in abundance data
abundance_df <- readRDS("Data/Output/O_24_abundance_summary.rds")

# Create WMU groups based on criteria
wmu_shapefile <- wmu_shapefile %>%
  mutate(WMU_Group = case_when(
    substr(WMU_ID, 1, 1) == "1" ~ "Region 1",
    WMU_ID %in% c("2A", "2C", "2D", "2E") ~ "Region 2",
    WMU_ID == "2B" ~ "Region 3",
    WMU_ID %in% c("2F", "2G", "2H", "3A", "3B") ~ "Region 4",
    WMU_ID == "3C" ~ "Region 5",
    WMU_ID == "3D" ~ "Region 6",
    WMU_ID %in% c("4A", "4B", "4D") ~ "Region 7",
    WMU_ID %in% c("4C", "4E") ~ "Region 8",
    WMU_ID %in% c("5A", "5B") ~ "Region 9",
    WMU_ID %in% c("5C", "5D") ~ "Region 10"
  )) %>% 
  group_by(WMU_Group) %>%
  summarize(geometry = st_union(geometry)) %>%
  ungroup()

# Summarize abundance data
abundance_df_totals <- abundance_df %>%
  group_by(Region, year, sex) %>%
  summarise(
    median_value = sum(median_value),
    lower_ci = sum(lower_ci),
    upper_ci = sum(upper_ci),
    area_sq_km = first(area_sq_km),
    .groups = "drop"
  )

# Rename "Group" to "Region" in the 'group' column
# abundance_df_totals <- abundance_df_totals %>%
#   mutate(group = str_replace(group, "Group", "Region"))

# Add 'Total' category by summing across sex and age_class
abundance_df_overall <- abundance_df %>%
#  mutate(group = str_replace(group, "Group", "Region")) %>% 
  group_by(Region, year) %>%
  summarise(
    median_value = sum(median_value),
    lower_ci = sum(lower_ci),
    upper_ci = sum(upper_ci),
    area_sq_km = first(area_sq_km),
    .groups = "drop"
  ) %>%
  mutate(sex = "Total") %>%
  bind_rows(abundance_df_totals) %>% 
  filter(sex == "Total") 

# Filter data for 2020 and 2023 and reshape it
abundance_df_filtered <- abundance_df_overall %>%
  filter(year %in% c(2, 5)) %>%
  select(Region, year, median_value)

# Reshape and calculate relative change over time
abundance_df_wide <- abundance_df_filtered %>%
  group_by(Region) %>% 
  pivot_wider(
    names_from = year,            
    values_from = median_value,   
    names_prefix = "year_"
  ) %>%
  mutate(
    # Calculate relative change: (year_4 - year_1) / year_1
    abundance_change_relative = (`year_5` - `year_2`) / `year_2`, 
    trend_color = case_when(
      abundance_change_relative > 0 ~ "increase",
      abundance_change_relative < 0 ~ "decrease",
      TRUE ~ "no change"
    )
  )

# Merge abundance data with shapefile
wmu_data <- wmu_shapefile %>%
  left_join(abundance_df_wide, by = c("WMU_Group" = "Region"))

# Filter out Region 10 for plotting
wmu_data_filtered <- wmu_data %>%
  filter(WMU_Group != "Region 10")

# Filter out Region 10 for plotting
wmu_10 <- wmu_shapefile %>%
  filter(WMU_Group == "Region 10")

# Plot relative abundance change by WMU
change_plot <- ggplot(wmu_data_filtered) +
  geom_sf(aes(fill = abundance_change_relative), color = "white", size = 0.1) +
  scale_fill_gradientn(
    colors = c(
      "#4C2E2F",  # Very strong decrease
      "#6E4546",  # Strong decrease
      "#8C5A5B",  # Moderate decrease
      "#B07B7E",  # Slight decrease
      "#F1E2CA",  # No change (this should map to 0)
      "#D4E6D0",  # Very light sage for minimal increase
      "#B8D4B0",  # Light sage green for small increase
      "#88A295",  # Muted sage green for slight increase
      "#628287",  # Balanced blue-green for good increase
      "#2B5E60",  # Dark blue-teal for exceptional increase
      "#1F4158",  # Nearly navy blue-teal for extraordinary increase
      "#183D3E"   # Deepest blue-teal for maximum increase
    ),
    # Map the 5th color (beige) to the value 0 within your range
    values = scales::rescale(c(-0.1, -0.075, -0.05, -0.025, 0, 0.1, 0.2, 0.4, 0.6, 0.8, 0.9, 1), 
                             to = c(0, 1)),
    limits = c(-0.1, 1),
    name = expression('ΔN'),
    guide = guide_colorbar(label = T)
  ) +
  geom_sf_pattern(data = wmu_10, fill = NA, color = "grey",
                  alpha = 0.8,
                  pattern = "stripe",
                  pattern_fill = "grey",
                  pattern_colour = NA,
                  pattern_density = 0.2,
                  pattern_spacing = 0.03) +
  theme_void() +
  geom_sf_text(data = wmu_data_filtered, aes(label = WMU_Group)) +
  geom_sf_text(data = wmu_10, aes(label = WMU_Group)) +
  theme(legend.position = "right") +
  xlab(NULL) +
  ylab(NULL) +
  theme(text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank())

change_plot

ggsave(file.path(paste0("Datavis/wmu_group_relative-change-1024.png")), plot = change_plot,
       device = "png", bg="transparent", dpi = 700, height = 5, width = 8)

# ggsave(file.path(paste0("Datavis/20241229-wmu_group_relative-change-10-21-23.png")), plot = change_plot,
#        device = "png", bg="transparent", dpi = 700, height = 5, width = 8)
# 

# For manu
# Assuming you have a variable in `wmu_data_filtered` like `abundance_change` and you want to plot triangles for significant increases

# Add a new variable for the shape based on `abundance_change` (e.g., for large changes)
wmu_data_filtered$shape <- ifelse(wmu_data_filtered$abundance_change > 0.1, 17,   # Triangle for positive change
                                  ifelse(wmu_data_filtered$abundance_change < -0.1, 24,  # Triangle down for negative change
                                         16))  # Circle for no or small change


##############################

# Summarize abundance data
abundance_df_totals <- abundance_df %>%
  group_by(Region, year, sex) %>%
  filter(year == "5") %>% 
  summarise(
    density_value = sum(density_value),
    lower_ci = sum(lower_ci),
    upper_ci = sum(upper_ci),
    area_sq_km = first(area_sq_km),
    .groups = "drop"
  )

# Add 'Total' category by summing across sex and age_class
abundance_df_overall2 <- abundance_df %>%
  #mutate(group = str_replace(group, "Group", "Region")) %>% 
  group_by(Region, year) %>%
  summarise(
    density_value = sum(density_value),
    lower_ci = sum(lower_ci),
    upper_ci = sum(upper_ci),
    area_sq_km = first(area_sq_km),
    .groups = "drop"
  ) %>%
  mutate(sex = "Total") %>%
  bind_rows(abundance_df_totals) %>% 
  filter(sex == "Total")

# Merge abundance data with shapefile
wmu_data2 <- wmu_shapefile %>%
  left_join(abundance_df_overall2, by = c("WMU_Group" = "Region"))

# Filter out Group 10 for plotting
wmu_data_filtered2 <- wmu_data2 %>%
  filter(WMU_Group != "Region 10")


den_plot <- ggplot(wmu_data_filtered2) +
  geom_sf(aes(fill = density_value), color = "white") +  # Fill based on density
  scale_fill_gradientn(
    colors = c(
      "white",
      "#F5F5DC",
      "#E8E8F3",
      "#DCEBF2", # Soft, very light warm blue (original low value)
      "#A6B9C9",  # Light, soft muted blue for low values
      "#8A9CB3",  # Muted light blue for moderate decrease
      "#7D8E9B",  # Slightly darker muted blue for slight decrease
      "#5D6E7D",  # Grayish blue for moderate increase
      "#3B4A56",  # Dark blue-teal for the highest value
      "#2A343C"   # Deeper, darker blue-teal for very high values
    ),
    name = expression('Abundance/km'^2*''),
    guide = guide_colorbar(label = T),
  values = scales::rescale(c(0, 5)),   # Set the rescaling range to match -0.4 to 0.4
   limits = c(0, 7)) +                       # Set legend limits from -0.4 to 0.4 
  geom_sf_pattern(data = wmu_10, fill = NA, color = "grey",
          alpha = 0.8,
          pattern = "stripe",
          pattern_fill = "grey",
          pattern_colour = NA,
          pattern_density = 0.2,
          pattern_spacing = 0.03) +
  theme_void() +
  geom_sf_text(data = wmu_10, aes(label = WMU_Group)) +
  geom_sf_text(data = wmu_data_filtered2, aes(label = WMU_Group)) +
  theme(legend.position = "right")
den_plot


ggsave(file.path(paste0("Datavis/wmu_group_den-10.png")), plot = den_plot,
       device = "png", bg="transparent", dpi = 700, height = 6, width = 9)

# Stack plots
# Add to each plot before combining
change_plot <- change_plot + 
  labs(tag = "A") +
  theme(plot.tag = element_text(size = 16, face = "bold", hjust = 0, vjust = 1))

den_plot <- den_plot + 
  labs(tag = "B") +
  theme(plot.tag = element_text(size = 16, face = "bold", hjust = 0, vjust = 1))

# Then combine however you're currently combining them
stack <- change_plot /den_plot

# If stacking vertically (one above the other)
ggsave(file.path("SubmissionMaterial/MajorRevisions/PubFigs/combined_plot.pdf"), 
       plot = stack,  # patchwork vertical stack
       device = "pdf", 
       bg = "transparent", 
       dpi = 700, 
       height = 9,  # Slightly less than 2x to account for shared elements
       width = 9)

# Line graph for density over time
# Add 'Total' category by summing across sex and age_class
#abundance_df <- readRDS("Data/Output/Simple_20241105_abundance_summary.rds")
## Structure abundance df ----
# Summarize data for total males, total females, and overall totals
abundance_df_totals <- abundance_df %>%
  group_by(Region, year, sex) %>%
  summarise(
    median_value = sum(median_value),
    lower_ci = sum(lower_ci),
    upper_ci = sum(upper_ci),
    area_sq_km = first(108759.7),  # assuming area is the same across sex/age
    .groups = "drop"
  )

# Add 'Total' category by summing across sex and age_class 
abundance_df_overall <- abundance_df %>%
  group_by(Region, year) %>%
  summarise(
    median_value = sum(median_value),
    lower_ci = sum(lower_ci),
    upper_ci = sum(upper_ci),
    area_sq_km = first(108759.7),  # assuming area is the same across sex/age
    .groups = "drop"
  ) %>%
  mutate(sex = "Total") %>%
  bind_rows(abundance_df_totals)  # Append the total rows to the original dataframe

# Aggregating to get statewide estimates
statewide_data <- abundance_df %>%
  group_by(year, sex, age_class) %>%  # Group by year, sex, and age_class for statewide summary
  summarize(
    mean_value = mean(mean_value, na.rm = TRUE),           # Average mean_value across groups
    median_value = median(median_value, na.rm = TRUE),     # Median of median_values across groups
    lower_ci = median(lower_ci/108759.7, na.rm = TRUE),    # Aggregate lower confidence interval
    upper_ci = median(upper_ci/108759.7, na.rm = TRUE),    # Aggregate upper confidence interval
    density_value = median_value/108759.7,                 # Average density values
    area_sq_km = sum(area_sq_km, na.rm = TRUE),            # Total area across groups
    .groups = "drop"  # Ensure that the groups are dropped after summarizing
  ) %>%
  ungroup()  # Ungroup after summarizing to avoid issues with further operations

# Pivoting data to show density values by year and group
statewide_density_table <- statewide_data %>%
  pivot_wider(
    names_from = c(sex, age_class),   # Create columns for each combination of sex and age_class
    values_from = density_value,      # Use density_value for each combination
    names_prefix = "Density_"         # Prefix to indicate the density
  ) %>% 
  distinct(year, Density_Female_Adult, Density_Female_Juvenile, Density_Male_Adult, Density_Male_Juvenile) %>%
  filter(
    !is.na(Density_Female_Adult) & 
      !is.na(Density_Female_Juvenile) & 
      !is.na(Density_Male_Adult) & 
      !is.na(Density_Male_Juvenile)
  ) %>% 
  


# View the table
print(statewide_data)

write.csv(statewide_data, "Output/statewide_density_data.csv", row.names = FALSE)

combo_colors <- c(
  "Adult Female"    = "#71250F",
  "Juvenile Female" = "#EB781B",
  "Adult Male"      = "#365365",
  "Juvenile Male"   = "#6F6534"
)


state <- ggplot(data = statewide_data, aes(
  x = year,
  y = density_value, 
  shape = sex,
  color = interaction(age_class, sex, sep = " "),
  fill = interaction(age_class, sex, sep = " ")  # Match fill with color
)) +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci), alpha = 0.2) +
  geom_point(size = 2) +
  scale_shape_manual(values = c("Female" = 16, "Male" = 17)) +
  scale_color_manual(values = combo_colors) +    # Apply custom colors to lines and points
  scale_fill_manual(values = combo_colors) +     # Apply custom colors to ribbons
  guides(color = guide_legend(override.aes = list(shape = c(16, 16, 17, 17))),  # Adjust color legend shapes
         shape = guide_none()) +  # Hide separate shape legend
  scale_x_continuous(breaks = c(1, 2, 3, 4, 5), labels = c("2020", "2021", "2022", "2023", "2024")) +  # Change x-axis labels
  #ylim(0, 0.10) +
  labs(
    title = 'Statewide Density (Birds/km²)',
    x = "",
    y = "",
    color = "",
    fill = ""
  ) +
  theme_bw() +
  theme(
    text = element_text(size = 20),           # Increase the base text size
    axis.title = element_text(size = 20),     # Increase axis title size
    axis.text = element_text(size = 20),      # Increase axis text (ticks) size
    strip.text = element_text(size = 20),     # Increase facet label text size
    legend.title = element_text(size = 20, face = "bold"),   # Customize legend title text
    legend.text = element_text(size = 20),    # Customize legend item text
    legend.position = "bottom",                # Position the legend
    legend.key = element_rect(fill = "white"),  # Customize legend key appearance
    legend.key.size = unit(1.5, "lines")      # Adjust size of legend keys (shapes)
  )
ggsave("Datavis/statewide_O-IPM_comparison.png", plot = state, width = 15, height = 10, dpi = 700)
