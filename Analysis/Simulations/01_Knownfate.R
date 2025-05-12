###############################################################################X
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#################### Integrated Population Model (IPM): #######################X
#                 #---# PhD Dissertation: Chapter 1 #---#
#        Creating a Bayesian IPM to inform turkey management in PA
###                                                                         ###X
#    Modeling: Creating a  known-fate model for telemetered hens by 
# wildlife management unit (WMU), and age class over months prior to
#                            fall hunting season.
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #X
# 
# Created by: Veronica A. Winter
# Last edited: 12/13/2024
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #X
###############################################################################X
# For simulation for this model, see 'Known_fate-sim.R'
#############################################################################X
#                       ----- Real data run -----
#############################################################################X
# clean env
rm(list = ls())
gc()
###-----------------------------------------------------#X
# Set seed for reproducibility
set.seed(1235)
###-----------------------------------------------------#X
# Load necessary libraries
library(RODBC)
library(nimble)
###-----------------------------------------------------#X
# source scripts
# source("Analysis/Scripts/00_IPM_funs.R")
load("Data/Research_IPM_run.Rdata")

#############################################################X
##           Source model ----
#############################################################X
# source model
source("Models/Individual_ipm-component_models/known_fate.R")
#############################################################X
##            Model run set up ----
#############################################################X
surv_data = list(
    # KF telemetered data 
    status = status_matrix,
    telem.juvenile = is_juvenile_matrix, # adult is intercept
    telem.wmu = telem.wmu # 2D is intercept term
)

# Constants
consts <- list(
  # KF telemetered constants 
  telem.nind = telem.nind,# 405
  telem.first = telem.first,
  telem.last = telem.last,
  telem.year.start = telem.year.start,
  telem.year.end = telem.year.end,
  female.telem.wmu = 4
)

inits <- list(
  # KF telemetered inits
  # Intercepts and coefficients
  telem.beta.int = rnorm(1, 0, 0.5),       
  telem.beta.age = rnorm(1, 0, 0.5),       
  telem.beta.wmu = rnorm(4, 0, 1),        
  telem.beta.month = rnorm(12, 0, 1),      
  
  # Variance parameters
  telem.sigma = runif(1, 0.1, 10),          
  telem.month.sigma = runif(1, 0, 10),
  # Survival probabilities 
  s.kf = array(runif(1, 0.2, 1),           
               dim = c(telem.nind, 4, 12))
)

# Create the Nimble model
model <- nimbleModel(
  code = survival,
  data = surv_data,
  constants = consts,
  inits = inits
)

# Check initialization issues
model$initializeInfo()
gc()

# Set MCMC specifications
burn <- 60000
iter <-  100000
t <- 10
chain <- 2

# Run MCMC 
survival_results <- nimbleMCMC(
  model = model,
  monitors = c("telem.beta.wmu", "telem.beta.int", "telem.beta.age",
               "s.kf", "telem.beta.month", "storage", "telem.sigma",
               "avg.ad.s.kf", "avg.juv.s.kf"
               ),
  niter = iter, 
  nburnin = burn, 
  thin = t, 
  nchains = chain, 
  setSeed = 1235, 
  samplesAsCodaMCMC = TRUE,
  WAIC = F)
 beepr::beep(1)
gc()

#############################################################X
##           Model diagnostics ----
#############################################################X
### Check rhats ----
MCMCvis::MCMCsummary(survival_results, params = c("avg.ad.s.kf", "storage"))
MCMCvis::MCMCsummary(survival_results, params = c("telem.beta.month", "telem.beta.wmu"))

### Check traceplot ----
MCMCtrace(survival_results, pdf = F, iter = iter, params = c("storage"))
MCMCtrace(survival_results, pdf = F, iter = iter, params = c("avg.ad.s.kf"))
MCMCtrace(survival_results, pdf = F, iter = iter, params = c("total.avg.ad.s.kf"))
MCMCtrace(survival_results, pdf = F, iter = iter, params = c("total.avg.juv.s.kf"))

#############################################################X
# Data viz ---
################################X
# Load necessary libraries
library(dplyr)
library(tidyr)
library(coda)
library(stringr)
library(purrr)
library(ggplot2)

# source
source("Analysis/00_output-processing_funs.R")

# Process each category using the function ----
samples_df <- as.data.frame(survival_results$chain1)

kf_female_surv_ad_df <- process_category_wmu(samples_df, "avg.ad.s.kf", "Female", "Adult")
#
kf_female_surv_juv_df <- process_category_wmu(samples_df, "avg.juv.s.kf", "Female", "Juvenile")


# Combine df
kf_survival_df <- bind_rows(
  kf_female_surv_ad_df, kf_female_surv_juv_df
) %>% 
  mutate(wmu = factor(wmu, levels = c(1, 2, 3, 4), 
                       labels = c("2D", "3D", "4D", "5C")),
          demographic_est = "KF_Survival")


# plot Female survival
kf_survival_plot <- ggplot(kf_survival_df, aes(y = median_value, x = wmu, 
                                               shape = sex, color = interaction(age_class, sex, sep = " "))) +
  geom_point(size = 4, position = position_dodge(width = 0.5)) +
  geom_errorbar(aes(ymin = lower_ci, ymax = upper_ci), 
                width = 0.3, 
                position = position_dodge(width = 0.5), show.legend = F) + 
  
  # Match colors for females as per the harvest plot
  scale_color_manual(values = c(
    "Adult Female"    = "#71250F",  # Color for adult female
    "Juvenile Female" = "#EB781B"  # Color for juvenile female
  )) +
  
  # Set shapes for age class (as you did before)
  scale_shape_manual(values = c("Female" = 16)) +
  
  # Remove shape legend using guides()
  guides(shape = "none") +
  
  # Y-axis limits
  scale_y_continuous(limits = c(0, 1)) +
  
  # Labels for axes and legend
  labs(x = "WMU", y = "Survival probability", shape = "", color = "") +
  
  # Keep consistent theme
  theme_classic() +
  theme(
    text = element_text(size = 16),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 16),
    strip.text = element_blank(),      # Remove WMU labels
    legend.title = element_text(size = 16, face = "bold"),
    legend.text = element_text(size = 16),
    legend.position = "top",
    legend.key = element_rect(fill = "white"),  # Customize legend key appearance
    legend.key.size = unit(1.5, "lines")
  )

ggsave("Datavis/20250326_kf_survival_plot_wmu_2.png", plot = kf_survival_plot,  width = 8, height = 6, dpi = 700)
#### NOTE: Writing regressions in Bayesian options: ----
# logit(s[i]) <- beta.age[1]*age.ad[i] + beta.age[2]*age.juv[i] ## mult vectors of age
# logit(s[i]) <- beta.age[1:2]%*%age[i,1:2]

library(dplyr)
library(knitr)
library(kableExtra)
library(flextable)

# Define the data
kf_survival_df2_f <- kf_survival_df %>%
  mutate(wmu = factor(wmu, levels = c(1, 2, 3, 4), 
                      labels = c("2D", "3D", "4D", "5C"))) %>% 
  mutate(Sex_Age_Class = paste(sex, age_class)) %>%
  filter(sex == "female") %>%
  dplyr::select(-c(sex, age_class, demographic_est)) %>%
  relocate(Sex_Age_Class, .before = median_value) %>%
  relocate(wmu, .before = median_value) %>%
  # relocate(Year, .before = median_value) %>%
  rename(Median = median_value,
         "2.5% CI" = lower_ci,
         "97.5% CI" = upper_ci,
         WMU = wmu) %>%
  # Round numerical values to 3 decimal places
  mutate(across(c(Median, `2.5% CI`, `97.5% CI`), ~ round(., 3))) %>%
  dplyr::select(-c(mean_value, jittered_wmu))

# Generate the table using flextable
ft <- flextable(kf_survival_df2_f)

# Add caption
ft <- set_caption(ft, "Survival probabilities with 95% credible intervals for each sex/age class across WMUs")

# Adjust column widths
ft <- autofit(ft)

# Optionally, style the table (bold headers, etc.)
ft <- bold(ft, part = "header") %>%
  align(align = "center", part = "all") %>%
  border_remove() %>%
  theme_vanilla()

ft