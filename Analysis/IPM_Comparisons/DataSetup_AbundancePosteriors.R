library(ggridges)
library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)

# Load the data
load("Data/Output/O_IPM_run24.Rdata")  # Operational model data

# Create data frames for the posterior samples
operational_samples <- as.data.frame(results) 

# Select only abundance columns
operational_abundance <- operational_samples %>%
  select(starts_with("female.N.ad"), starts_with("female.N.juv"), 
         starts_with("male.N.ad"), starts_with("male.N.juv"))

# Convert to long format
operational_long <- operational_abundance %>%
  pivot_longer(
    cols = everything(),
    names_to = "parameter",
    values_to = "abundance"
  ) %>%
  mutate(model = "Operational")

# Save
saveRDS(operational_long, "Data/Output/oipm24_abundance.rds")

# Vague IPM -----
# Load the data
load("Data/Output/V_IPM_run24.Rdata")  # Vague model data

# Create data frames for the posterior samples
vague_samples <- as.data.frame(results)   

# Select only abundance columns
vague_abundance <- vague_samples %>%
  select(starts_with("female.N.ad"), starts_with("female.N.juv"), 
         starts_with("male.N.ad"), starts_with("male.N.juv"))


# Convert to long format
vague_long <- vague_abundance %>%
  pivot_longer(
    cols = everything(),
    names_to = "parameter",
    values_to = "abundance"
  ) %>%
  mutate(model = "Vague")

# Save
saveRDS(vague_long, "Data/Output/vipm24_abundance.rds")

# Done