# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#################### Integrated Population Model (IPM): #######################X
#                 #---# PhD Dissertation: Simple IPM #---#
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
library(dplyr)
library(ggplot2)

# of data output
date <- "20250326"

# Load data 
abundance_df_1 <- readRDS(paste0("Data/Output/Full_", date, "_abundance_summary.rds"))
abundance_df_2 <- readRDS(paste0("Data/Output/Simple_", date, "_abundance_summary.rds")) %>% 
  mutate(wmu = as.character(wmu))

# Load in WMU areas
wmu_areas1 <- readRDS("../../PSUTurkey/turkey_IPM/Data/wmu_areas_km.rds")
wmu_areas2 <- readRDS("../../PSUTurkey/turkey_IPM/Data/wmu_km_areas_w.groups.rds")
wmu_areas <- rbind(wmu_areas1, wmu_areas2)

# Color palettes
colors <- c("#EB781B", "#CC5221", "#71250F", "#6F6534", "#365365")
colors2 <- c("#F6D9C0", "#C74442", "#71250E", "#6D383E", "#365365")

region_colors <- c(
  "3D" = "#043A50", # goldenrod
  "Region 6" = "#5f0f40"
)

# Combine the data frames as shown in your example
abundance_df <- left_join(abundance_df_1, abundance_df_2)

# Structure abundance df
# First, filter for only the regions we want (3D and Region 6)
filtered_abundance_df <- abundance_df_1 %>%
  filter(wmu %in% c("3D"))

filtered_abundance_df2 <- abundance_df_2 %>%
  filter(Region %in% c("Region 6")) %>% 
  select(-wmu) %>% 
  rename(wmu = Region)

# Summarize data by sex (combining age classes)
abundance_df_by_sex <- filtered_abundance_df %>%
  bind_rows(filtered_abundance_df2) %>% 
  group_by(wmu, year, sex) %>%
  summarise(
    median_value = sum(median_value),
    lower_ci = sum(lower_ci),
    upper_ci = sum(upper_ci),
    area_sq_km = first(area_sq_km),
    .groups = "drop"
  ) %>%
  # Ensure we only have Male and Female (no Total)
  filter(sex %in% c("Male", "Female"))

# Create the new plot with sex as facets and regions distinguished by color/shape
new_abundance_plot <- ggplot(abundance_df_by_sex, 
                             aes(y = median_value / area_sq_km, 
                                 x = year, 
                                 shape = wmu, 
                                 color = wmu)) +
  geom_point(size = 4, position = position_dodge(width = 0.5)) +
  geom_errorbar(aes(ymin = lower_ci/area_sq_km, ymax = upper_ci/area_sq_km), 
                width = 0.3, position = position_dodge(width = 0.5), show.legend = F) +
  
  # Use sex as facets
  facet_wrap(~sex, ncol = 2) +
  
  # Custom shapes and colors for regions
  scale_shape_manual(values = c("3D" = 17, "Region 6" = 16)) + 
  scale_color_manual(values = c("3D" = "#6A87C1", "Region 6" = "#F4A848")) +
  
  scale_x_continuous(breaks = c(1, 2, 3, 4), labels = c("2020", "2021", "2022", "2023")) +
  
  labs(x = "Year", y = expression('Abundance/km'^2*''), shape = "WMU", color = "WMU") +
  
  # Combine the shape and color legends
  guides(color = guide_legend(override.aes = list(shape = c(17, 16))),
         shape = guide_none()) +
  
  annotate("segment", x = Inf, xend = Inf, y = -Inf, yend = Inf, color = "black", linetype = "dashed") +
  theme_classic() +
  theme(
    text = element_text(size = 20),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 18),
    strip.text = element_text(size = 18, face = "bold"),
    legend.title = element_text(size = 18, face = "bold"),
    legend.text = element_text(size = 18),
    legend.position = "top",
    legend.key = element_blank(),
    legend.key.size = unit(1.5, "lines"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

# Print the plot
print(new_abundance_plot)

# If you want to save the plot
ggsave("Dataviz/Fig5.png", new_abundance_plot, width = 12, height = 6, dpi = 300)
ggsave("SubmissionMaterial/Fig5.pdf", new_abundance_plot, width = 12, height = 6, dpi = 300)
