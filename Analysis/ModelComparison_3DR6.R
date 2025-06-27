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
library(kableExtra)

# Load in WMU areas
wmu_areas1 <- readRDS("Data/wmu_areas_km.rds")
wmu_areas2 <- readRDS("Data/wmu_km_areas_regions.rds")
wmu_areas <- rbind(wmu_areas1, wmu_areas2)

# Load data 
abundance_df_1 <- readRDS(paste0("Data/Output/R__abundance_summary.rds")) %>% 
  mutate(Model = "Research")
abundance_df_2 <- readRDS(paste0("Data/Output/O_abundance_summary.rds")) %>% 
  mutate(wmu = as.character(wmu))%>% 
  mutate(Model = "Operational")
abundance_df_3 <- readRDS(paste0("Data/Output/V_abundance_summary.rds")) %>% 
  mutate(wmu = as.character(wmu)) %>% 
  select(-c(area_sq_km, density_value)) %>% 
  left_join(wmu_areas2, by = c("Region" = "WMU_ID")) %>% 
  mutate(density_value = median_value/area_sq_km)%>% 
  mutate(Model = "Vague")

# Color palettes
colors <- c("#EB781B", "#CC5221", "#71250F", "#6F6534", "#365365")
colors2 <- c("#F6D9C0", "#C74442", "#71250E", "#6D383E", "#365365")

region_colors <- c(
  "3D" = "#043A50", # goldenrod
  "Region 6" = "#5f0f40"
)

# Structure abundance df
# First, filter for only the regions we want (3D and Region 6)
filtered_abundance_df <- abundance_df_1 %>%
  filter(wmu %in% c("3D"))

filtered_abundance_df2 <- abundance_df_2 %>%
  filter(Region %in% c("Region 6")) %>% 
  select(-wmu) %>% 
  rename(wmu = Region)

filtered_abundance_df3 <- abundance_df_3 %>%
  filter(Region %in% c("Region 6")) %>% 
  select(-wmu) %>% 
  rename(wmu = Region)

# Summarize data by sex (combining age classes)
abundance_df_by_sex <- filtered_abundance_df %>%
  bind_rows(filtered_abundance_df2) %>%
  bind_rows(filtered_abundance_df3) %>%
  group_by(wmu, year, sex, Model) %>%
  summarise(
    median_value = sum(median_value),
    lower_ci = sum(lower_ci),
    upper_ci = sum(upper_ci),
    area_sq_km = first(area_sq_km),
    .groups = "drop"
  ) %>%
  # Ensure we only have Male and Female (no Total)
  filter(sex %in% c("Male", "Female")) %>%
  # Convert Model to factor with desired order
  mutate(Model = factor(Model, levels = c("Research", "Operational", "Vague")))

# Create the new plot with sex as facets and regions distinguished by color/shape
new_abundance_plot <- ggplot(abundance_df_by_sex, 
                             aes(y = median_value/area_sq_km, 
                                 x = year, 
                                 shape = Model, 
                                 color = Model)) +
  geom_point(size = 4, position = position_dodge(width = 0.5)) +
  geom_errorbar(aes(ymin = lower_ci/area_sq_km, ymax = upper_ci/area_sq_km), 
                width = 0.3, position = position_dodge(width = 0.5), show.legend = F) +
  
  # Use sex as facets
  facet_wrap(~sex, ncol = 2) +
  
  # Custom shapes and colors for regions
  scale_shape_manual(values = c("Research" = 17, "Operational" = 16, "Vague" = 15)) + 
  scale_color_manual(values = c("Research" = "#92b9d8", "Operational" = "#B34170", "Vague" = "#8856a7")) +
  
  scale_x_continuous(breaks = c(1, 2, 3, 4), labels = c("2020", "2021", "2022", "2023")) +
  
  labs(x = "Year", y = expression('Abundance/km'^2*''), shape = "", color = "") +
  
  # Combine the shape and color legends
  guides(color = guide_legend(override.aes = list(shape = c(17, 16, 15))),
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

# Table for appendix
abundance_tab <- filtered_abundance_df %>%
  bind_rows(filtered_abundance_df2) %>%
  bind_rows(filtered_abundance_df3) %>%
  mutate(CrIWidth = (upper_ci/area_sq_km) - (lower_ci/area_sq_km),
        lower_ci2 = lower_ci/area_sq_km,
        upper_ci2 = upper_ci/area_sq_km,
        SexAge = paste(age_class, sex, sep = " ")) %>% 
  mutate(Year = case_when(year == 1 ~ "2020",
                          year == 2 ~ "2021",
                          year == 3 ~ "2022",
                          year == 4 ~ "2023"),
         density_value = median_value/area_sq_km) %>% 
  # Round numerical values to 3 decimal places
  mutate(across(c(density_value, lower_ci2, upper_ci2, CrIWidth), ~ round(., 3))) %>% 
  dplyr::select(Model, wmu, SexAge, Year, density_value, lower_ci2, upper_ci2,
                CrIWidth) %>% 
  arrange(Model, SexAge)


table_reg <- kable(abundance_tab,
                   booktabs = TRUE,   
                   format = "html", 
                   col.names = c("Model","WMU/Region", "Sex & AGe Class", "Year", "Density", "2.5% CrI", "97.5% CrI", "CrI width")) %>%
  # Collapse rows for Sex/Age Class and WMU
  collapse_rows(columns = c(1, 2, 3), latex_hline = "major", valign = "middle") %>%
  kable_styling(latex_options = c("hold_position", "repeat_header"), full_width = FALSE) 

save_kable(table_reg, file = "Dataviz/Appx5Tb1.html")
