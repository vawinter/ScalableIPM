rm(list=ls())
gc()

library(dplyr)
library(ggridges)
library(ggplot2)
library(coda)
library(kableExtra)

# Function to calculate CV and aggregate by base parameter
compute_cv <- function(data, model_name) {
  # Compute CV for each parameter
  cv_values <- apply(data, 2, function(x) {
    if (mean(x, na.rm = TRUE) != 0) {
      return(sd(x, na.rm = TRUE) / mean(x, na.rm = TRUE))
    } else {
      return(NA)
    }
  })
  
  # Convert to data frame
  cv_df <- data.frame(Parameter = names(cv_values), CV = cv_values, row.names = NULL)
  
  # Extract base parameter name (remove indices using regex)
  cv_df$Base_Parameter <- gsub("\\[.*\\]", "", cv_df$Parameter)
  
  # Compute mean CV for each base parameter
  cv_mean <- cv_df %>%
    group_by(Base_Parameter) %>%
    mutate(avgCV = mean(CV, na.rm = TRUE)) %>%
    # From the research model
    filter(Base_Parameter != "storage") %>%
    mutate(Model = model_name) %>% 
    ungroup()
  
  return(cv_mean)
}

# Load and process models
mod1 <- readRDS("Data/Output/R_IPM_run.rds")[[1]] %>% as.data.frame()
mod2 <- readRDS("Data/Output/20250623_Parallel_O_IPM_run.rds")[[1]] %>% as.data.frame()
vague <- readRDS("Data/Output/20250625_Parallel_O_vague_IPM_run.rds")[[1]] %>% as.data.frame()

# Process the posteriors
research_cv <- compute_cv(mod1, "Research")
operational_cv <- compute_cv(mod2, "Operational")
vague_cv <- compute_cv(vague, "Vague")

# Combine results
cv_df_all <- bind_rows(research_cv, operational_cv, vague_cv
                       ) %>% 
  mutate(Base_Parameter = recode(Base_Parameter,
                                 "male.N.juv" = "Juvenile Male Abundance",
                                 "female.N.juv" = "Juvenile Female Abundance",
                                 "male.N.ad" = "Adult Male Abundance",
                                 "female.N.ad" = "Adult Female Abundance",
                                 "avg.ad.s.kf" = "Adult Female Survival",
                                 "avg.juv.s.kf" = "Juvenile Female Survival",
                                 "female.h.ad.wmu" = "Adult Female Harvest Rate",
                                 "female.h.juv.wmu" = "Juvenile Female Harvest Rate",
                                 "male.h.ad.wmu" = "Adult Male Harvest Rate",
                                 "male.h.juv.wmu" = "Juvenile Male Harvest Rate",
                                 "male.s.ad.wmu" = "Adult Male Survival",
                                 "male.s.juv.wmu" = "Juvenile Male Survival",
                                 "recruitment" = "Recruitment"
  )) %>% 
  filter(!Base_Parameter %in% c("aug31.ppb", "aug31.hwb")) %>% 
  dplyr::select(Model, Base_Parameter, avgCV) %>% 
  mutate(avgCV = round(avgCV, 3)) %>% 
  distinct()


# Print result
print(cv_df_all)

saveRDS(cv_df_all, "Data/ModelComparison_CV.rds")

library(knitr)
library(kableExtra)
table_wmu <- kable(cv_df_all,
                   booktabs = TRUE,   
                   format = "html", 
                   col.names = c("Model", "Parameter","Average CV")) %>%
  # Collapse rows for Sex/Age Class and WMU
  collapse_rows(columns = c(1, 2, 3), latex_hline = "major", valign = "middle") %>%
  kable_styling(latex_options = c("hold_position", "repeat_header"), full_width = FALSE) 

save_kable(table_wmu, file = "Dataviz/Appx8Tb1.html")
##########################X
# Boxplots ----
##########################X
library(stringr)

box <- cv_df_all %>% 
  filter(Base_Parameter %in% c("Adult Female Abundance",
                               "Adult Female Survival",
                               "Adult Female Harvest Rate",
                               "Juvenile Female Abundance",
                               "Juvenile Female Survival",
                               "Juvenile Female Harvest Rate",
                               "Recruitment")) %>%  
  mutate(Model = factor(Model, levels = c("Research", "Operational", "Vague")),
         Base_Parameter = str_wrap(Base_Parameter, width = 10)) %>%  # Wrap at 10 characters
  ggplot(aes(x = Base_Parameter, y = CV, fill = Model)) +
  scale_fill_manual(values = c("Research" = "#b3cde3", "Operational" = "#B34170", "Vague" = "#8856a7")) +
  geom_boxplot(alpha = 0.7, outlier.shape = 16, outlier.size = 2) +
  labs(x = "", y = "CV", fill = "") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5),  # Remove angle since text is wrapped
        text = element_text(size = 14),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 14),
        strip.text = element_blank(),
        legend.title = element_text(size = 14, face = "bold"),
        legend.text = element_text(size = 14),
        legend.position = "top",
        legend.key = element_rect(fill = "white"),
        legend.key.size = unit(1.5, "lines"))

box

ggsave("Dataviz/Model_CV_rec.png", plot = box, width = 10, height = 6, dpi = 700)
ggsave("../../Manuscripts/ScalableIPM/Dataviz/Model_CV_rec.pdf", plot = box, width = 10, height = 6, dpi = 700)
#--------------X
male_box <- cv_df_all %>% 
  filter(Base_Parameter %in% c("Adult Male Abundance",
                               "Adult Male Survival",
                               "Adult Male Harvest Rate",
                               "Juvenile Male Abundance",
                               "Juvenile Male Survival",
                               "Juvenile Male Harvest Rate"))%>%  
  mutate(Model = factor(Model, levels = c("Research", "Operational", "Vague")),
         Base_Parameter = str_wrap(Base_Parameter, width = 10)) %>%  # Wrap at 10 characters
  ggplot(aes(x = Base_Parameter, y = CV, fill = Model)) +
  scale_fill_manual(values = c("Research" = "#b3cde3", "Operational" = "#B34170", "Vague" = "#8856a7")) +
  geom_boxplot(alpha = 0.7, outlier.shape = 16, outlier.size = 2) +
  labs(x = "", y = "CV", fill = "") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5),  # Remove angle since text is wrapped
        text = element_text(size = 14),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 14),
        strip.text = element_blank(),
        legend.title = element_text(size = 14, face = "bold"),
        legend.text = element_text(size = 14),
        legend.position = "top",
        legend.key = element_rect(fill = "white"),
        legend.key.size = unit(1.5, "lines"))  # Rotate axis labels for clarity

ggsave("Dataviz/appx6_fig1.png", plot = male_box, width = 10, height = 6, dpi = 700)




##########################X
# Density ----
##########################X
cv_df_all %>% 
  filter(Base_Parameter %in% c("Adult Female Abundance",
                               "Adult Female Survival",
                               "Adult Female Harvest Rate",
                               "Juvenile Female Abundance",
                               "Juvenile Female Survival",
                               "Juvenile Female Harvest Rate")) %>% 
ggplot(aes(x = CV, y = Base_Parameter, fill = Model)) +
  geom_density_ridges(alpha = 0.7) +  
  labs(x = "Coefficient of Variation (CV)", y = "Parameter", fill = "") +
  theme_minimal() +
  theme(legend.position = "top")
#--------------X
cv_df_all %>% 
  filter(Base_Parameter %in% c("Adult Male Abundance",
                               "Adult Male Survival",
                               "Adult Male Harvest Rate",
                               "Juvenile Male Abundance",
                               "Juvenile Male Survival",
                               "Juvenile Male Harvest Rate")) %>% 
  ggplot(aes(x = CV, y = Base_Parameter, fill = Model)) +
  geom_density_ridges(alpha = 0.7) +  
  labs(x = "Coefficient of Variation (CV)", y = "Parameter", fill = "") +
  theme_minimal() +
  theme(legend.position = "top")
##########################X
# Histograms ----
##########################X
cv_df_all %>% 
  filter(Base_Parameter %in% c("Adult Female Abundance",
                               "Adult Female Survival",
                               "Adult Female Harvest Rate",
                               "Juvenile Female Abundance",
                               "Juvenile Female Survival",
                               "Juvenile Female Harvest Rate")) %>% 
  ggplot(aes(x = CV, fill = Model)) +
  geom_histogram(binwidth = 0.1, alpha = 0.7, position = "identity") +
  labs(x = "Coefficient of Variation (CV)", y = "Frequency",
       fill = "Parameter") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  facet_wrap(~Base_Parameter)  # Facet by base parameter

cv_df_all %>% 
  filter(Base_Parameter %in% c("Adult Male Abundance",
                               "Adult Male Survival",
                               "Adult Male Harvest Rate",
                               "Juvenile Male Abundance",
                               "Juvenile Male Survival",
                               "Juvenile Male Harvest Rate")) %>% 
  ggplot(aes(x = CV, fill = Model)) +
  geom_histogram(binwidth = 0.1, alpha = 0.7, position = "identity") +
  labs(x = "Coefficient of Variation (CV)", y = "Frequency",
       fill = "Parameter") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  facet_wrap(~Base_Parameter)  # Facet by base parameter


