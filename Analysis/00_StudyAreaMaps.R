# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
################# Scalable Integrated Population Models : #####################X
#                   #---# PhD Dissertation: Study Area Map #---#
#        Creating a Bayesian IPM to inform turkey management in PA
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #X
###############################################################################X
rm(list = ls())
gc()

#Load packages----
library(lubridate) #for dates and times
library(dplyr) #for data.frame manipulation
library(ggplot2) # prefferd graphing package 
library(rnaturalearth) # for plotting maps (lines 25 - )
library("rnaturalearthdata")
library(geosphere)
library(ggspatial)
library(sf)
library(rgeos) 
library(maps)

# Data read in ----
# read in shapefile of WMUS
pa_wmu <- st_read("../../PSUTurkey/TurkeyProject/Data/PGC_BNDWildlifeManagementUnits2024.shp")
pa_public <- st_read("../../PSUTurkey/TurkeyProject/Data/State202412.shp")

# Formatting ----
# Create WMU groups based on criteria
groups <- pa_wmu %>%
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
  ))

# labels
group_lab <- c(paste("Region", 1:10, sep = " "))

# grab wmus for nest study
four <- pa_wmu %>% 
  filter(WMU_ID %in% c("2D","4D", "3D"))

# Load in global data
world <- ne_countries(scale = "medium", returnclass = "sf")
class(world)

# Add in states so I can label/input state boundaries
states <- st_as_sf(maps::map("state", plot = FALSE, fill = TRUE))
head(states)
class(states)

# get PA only
pa <- states %>% 
  filter(ID == "pennsylvania")

pa_groups <- groups %>% 
  group_by(WMU_Group) %>% 
  summarize(geometry = st_union(geometry)) %>%
  ungroup() %>% 
  mutate(WMU_Group = factor(WMU_Group,
                            levels = c(setdiff(WMU_Group, "Region 10"), "Region 10"))) |>
  arrange(WMU_Group)

# Maps ----
## Map of all wmus ----
grouped_wmu_map <- ggplot(data = pa) +
  geom_sf(fill = "white", color = "black") +
  
  # Region fills
  geom_sf(data = groups, aes(fill = WMU_Group), color = "black", alpha = 0.6) +
  
  # Telemetry outlines
  geom_sf(data = four, aes(color = "Telemetry"), fill = NA, linewidth = 2.5) +
  
  # ---- REGION COLORS ----
scale_fill_manual(
  values = c(
    "Region 1" = "#3c1518", 
    "Region 2" = "#0f5c5d", 
    "Region 3" = "#69140e",
    "Region 4" = "#C74442", 
    "Region 5" = "#5f0f40",
    "Region 6" = "#043A70", 
    "Region 7" = "goldenrod",
    "Region 8" = "#f79d65", 
    "Region 9" = "#606c38",
    "Region 10" = "#f2f3ae"
  ),
  breaks = c(
    "Region 1","Region 2","Region 3","Region 4","Region 5",
    "Region 6","Region 7","Region 8","Region 9","Region 10"
  ),
  name = ""
) +
  
  # ---- TELEMETRY LEGEND BELOW REGIONS ----
# We put telemetry second in the guides() call
scale_color_manual(
  name = "",
  values = c("Telemetry" = "black"),
  guide = guide_legend(
    override.aes = list(
      fill = NA,
      linewidth = 2,
      linetype = 1
))) +
  
  # ORDER LEGENDS: fill first, then color
  guides(
    fill = guide_legend(order = 1),
    color = guide_legend(order = 2)
  ) +
  
  # ---- SCALE BAR (top right) ----
ggspatial::annotation_scale(
  location = "tr",
  bar_cols = c("black", "white")
) +
  
  # ---- NORTH ARROW (bottom right) ----
ggspatial::annotation_north_arrow(
  location = "br",
  height = unit(1, "cm"),
  width = unit(1, "cm"),
  style = ggspatial::north_arrow_fancy_orienteering
) +
  
  theme_classic() +
  theme(
    text = element_text(size = 15, color = "black"),
    axis.ticks = element_blank(),
    legend.text = element_text(size = 12, color = "black")
  ) +
  geom_sf_text(data = pa_wmu, aes(label = WMU_ID)) +
  xlab(NULL) +
  ylab(NULL)
ggsave(file.path(paste0("SubmissionMaterial/MajorRevisions/PubFigs/study_area.png")), plot = grouped_wmu_map,
       device = "png", bg="transparent", dpi = 700, height = 5, width = 8)

ggsave(file.path(paste0("SubmissionMaterial/MajorRevisions/PubFigs/Fig1.pdf")), plot = grouped_wmu_map,
       device = "pdf", bg="transparent", dpi = 700, height = 5, width = 8)
