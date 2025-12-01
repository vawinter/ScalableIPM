# Appendix 7: Sensitivity Analysis of Female Harvest Rate Priors

library(dplyr)
library(ggplot2)
library(patchwork)
library(ggh4x)

# Table showing prior specifications
prior_comparison <- data.frame(
  Prior = c("Informed", "Wider"),
  Distribution = c("Beta(2, 50)", "Beta(3, 30)"),
  Mean = c(2/(2+50), 3/(3+30)),
  SD = c(sqrt((2*50)/((2+50)^2*(2+50+1))), 
         sqrt((3*30)/((3+30)^2*(3+30+1)))),
  `95% Interval` = c("[0.003, 0.107]", "[0.019, 0.241]")
)

# Plot showing both beta distributions
x <- seq(0, 0.5, length.out = 1000)
prior_densities <- data.frame(
  harvest_rate = rep(x, 2),
  density = c(dbeta(x, 2, 50), dbeta(x, 3, 30)),
  Prior = rep(c("Beta(2,50) - Informed", "Beta(3,30) - Wider"), each = 1000)
)

comp <- ggplot(prior_densities, aes(x = harvest_rate, y = density, color = Prior)) +
  geom_line(size = 1.2) +
  labs(x = "Female Harvest Rate", 
       y = "Density",
       title = "") +
  theme_classic() +
  theme(
    text = element_text(size = 16),           # Increase base text size
    axis.title = element_text(size = 16),     # Increase axis title size
    axis.text = element_text(size = 16),      # Increase axis text (ticks) size
    legend.title = element_blank(),   # Customize legend title text
    legend.text = element_text(size = 16),    # Customize legend item text
    legend.position = "top",                  # Position the legend at the top
    legend.key = element_rect(fill = "white"),  # Customize legend key appearance
    legend.key.size = unit(1.5, "lines")    # Adjust size of legend keys (shapes)
  )
ggsave("SubmissionMaterial/MajorRevisions/PubFigs/Appx7Fig1.png", plot = comp, device = "png",  bg="transparent",
       width = 10, height = 6, dpi = 700)


# Load posterior draws from both models
beta_2_50_results <- readRDS("Data/Output/R24_abundance_summary.rds")
beta_3_30_results <- readRDS("Data/Output/R24_beta330_abundance_summary.rds")
abundance_comparison <- bind_rows(
  beta_2_50_results %>% mutate(Prior = "Beta(2,50)"),
  beta_3_30_results %>% mutate(Prior = "Beta(3,30)")
) %>% 
  filter(sex == "Female") %>% 
  group_by(wmu, year, Prior) %>%  # Must include Prior here!
  summarise(
    median_value = sum(median_value),
    lower_ci = sum(lower_ci),
    upper_ci = sum(upper_ci),
    area_sq_km = first(area_sq_km),
    .groups = "drop"
  )

# Plot
abun <- ggplot(abundance_comparison, aes(x = year, y = median_value/area_sq_km, 
                                 shape = Prior, color = Prior)) +
  geom_point(size = 4, position = position_dodge(width = 0.3)) +
  geom_errorbar(aes(ymin = lower_ci/area_sq_km, ymax = upper_ci/area_sq_km), 
                width = 0.2, position = position_dodge(width = 0.3)) +
  facet_wrap2(~wmu, ncol = 3,
              strip = strip_themed(
                background_x = elem_list_rect(fill = c("WMU 2D" = "#0f4c5d",
                                                       "WMU 3D" = "#043A70",
                                                       "WMU 4D" = "goldenrod")))
  ) +
  scale_shape_manual(values = c("Beta(2,50)" = 16, "Beta(3,30)" = 17)) + 
  scale_color_manual(values = c("Beta(2,50)" = "#2C5F7C", "Beta(3,30)" = "#D97E3C")) +
  scale_x_continuous(breaks = c(1, 2, 3, 4, 5), 
                     labels = c("2020", "2021", "2022", "2023", "2024")) +
  ylim(0, 10) +
  labs(x = "Year", y = expression('Female Abundance/km'^2*''), 
       shape = "", color = "") +
  guides(color = guide_legend(override.aes = list(shape = c(16, 17))),
         shape = guide_none()) +
  annotate("segment", x = Inf, xend = Inf, y = -Inf, yend = Inf, 
           color = "black", linetype = "dashed") +
  theme_classic() +
  theme(
    text = element_text(size = 20),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 18),
    strip.text = element_text(size = 18, face = "bold", color = "white"),
    legend.title = element_text(size = 18, face = "bold"),
    legend.text = element_text(size = 18),
    legend.position = "top",
    legend.key = element_blank(),
    legend.key.size = unit(1.5, "lines"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

ggsave("SubmissionMaterial/MajorRevisions/PubFigs/Appx7Fig2.png", plot = abun, width = 15, height = 10, dpi = 700)

# Now, lets looks at the harvest rates
female_harvest_beta_2_50 <- readRDS("Data/Output/R24_harvest_summary.rds")
female_harvest_beta_3_30 <- readRDS("Data/Output/R24_beta330_harvest_summary.rds")
# Prepare harvest rate comparison data
harvest_comparison <- bind_rows(
  female_harvest_beta_2_50 %>% mutate(Prior = "Beta(2,50)"),
  female_harvest_beta_3_30 %>% mutate(Prior = "Beta(3,30)")
) %>% 
  filter(sex == "Female")  # Only female harvest rates
  

# Plot
harv1 <- harvest_comparison %>% 
  filter(age_class == "Adult") %>% 
  ggplot(aes(x = year, y = median_value, 
                               shape = Prior, color = Prior)) +
  geom_point(size = 4, position = position_dodge(width = 0.3)) +
  geom_errorbar(aes(ymin = lower_ci, ymax = upper_ci), 
                width = 0.2, position = position_dodge(width = 0.3)) +
  facet_wrap2(~wmu, ncol = 3,
              strip = strip_themed(
                background_x = elem_list_rect(fill = c("WMU 2D" = "#0f4c5d",
                                                       "WMU 3D" = "#043A70",
                                                       "WMU 4D" = "goldenrod")))
  ) +
  ylim(0, 0.05) +
  scale_shape_manual(values = c("Beta(2,50)" = 16, "Beta(3,30)" = 17)) + 
  scale_color_manual(values = c("Beta(2,50)" = "#2C5F7C", "Beta(3,30)" = "#D97E3C")) +
  scale_x_continuous(breaks = c(1, 2, 3, 4, 5), 
                     labels = c("2020", "2021", "2022", "2023", "2024")) +
  labs(x = "", y = expression('Adult Harvest Rate'), 
       shape = "", color = "") +
  guides(color = guide_legend(override.aes = list(shape = c(16, 17))),
         shape = guide_none()) +
  annotate("segment", x = Inf, xend = Inf, y = -Inf, yend = Inf, 
           color = "black", linetype = "dashed") +
  theme_classic() +
  theme(
    text = element_text(size = 20),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 18),
    strip.text = element_text(size = 18, face = "bold", color = "white"),
    legend.title = element_text(size = 18, face = "bold"),
    legend.text = element_text(size = 18),
    legend.position = "top",
    legend.key = element_blank(),
    legend.key.size = unit(1.5, "lines"),
    axis.text.x = element_blank()
  )

harv2 <- harvest_comparison %>% 
  filter(age_class != "Adult") %>% 
  ggplot(aes(x = year, y = median_value, 
             shape = Prior, color = Prior)) +
  geom_point(size = 4, position = position_dodge(width = 0.3)) +
  geom_errorbar(aes(ymin = lower_ci, ymax = upper_ci), 
                width = 0.2, position = position_dodge(width = 0.3)) +
  facet_wrap2(~wmu, ncol = 3,
              strip = strip_themed(
                background_x = elem_list_rect(fill = c("WMU 2D" = "#0f4c5d",
                                                       "WMU 3D" = "#043A70",
                                                       "WMU 4D" = "goldenrod")))
  ) +
  scale_shape_manual(values = c("Beta(2,50)" = 16, "Beta(3,30)" = 17)) + 
  scale_color_manual(values = c("Beta(2,50)" = "#2C5F7C", "Beta(3,30)" = "#D97E3C")) +
  scale_x_continuous(breaks = c(1, 2, 3, 4, 5), 
                     labels = c("2020", "2021", "2022", "2023", "2024")) +
  labs(x = "Year", y = expression('Juvenile Harvest Rate'), 
       shape = "", color = "") +
  guides(color = guide_legend(override.aes = list(shape = c(16, 17))),
         shape = guide_none()) +
  ylim(0, 0.05) +
  annotate("segment", x = Inf, xend = Inf, y = -Inf, yend = Inf, 
           color = "black", linetype = "dashed") +
  theme_classic() +
  theme(
    text = element_text(size = 20),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 18),
    strip.text = element_blank(),
    legend.title = element_text(size = 18, face = "bold"),
    legend.text = element_text(size = 18),
    legend.position = "none",
    legend.key = element_blank(),
    legend.key.size = unit(1.5, "lines"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
harvest <- harv1/harv2
ggsave("SubmissionMaterial/MajorRevisions/PubFigs/Appx7Fig3.png", plot = harvest, width = 15, height = 10, dpi = 700)
