###############################################################################X
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
############## Research Integrated Population Model (R_IPM): ##################X
#                     #---# PhD Dissertation: R_IPM #---#
#        Creating a Bayesian IPM to inform turkey management in PA
###                                                                         ###X
#    Modeling: known-fate for telemetered hens by wildlife management unit (WMU), 
#                       and age class over months prior to
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
library(nimble)
###-----------------------------------------------------#X
# source scripts
source("Analysis/00_IPM_funs.R")
##-----------------------------------##x
# # list files [KF 24]
# dir <- "Data/Research_IPM_setup-data/kf_data_22-23/"
# K <- list.files(dir, pattern = "\\.rds$", full.names = TRUE)
# names <- sub("\\.rds$", "", basename(K))
# 
# myfiles <- lapply(K, readRDS)
# names(myfiles) <- names
# list2env(myfiles, globalenv())
# 
# rm(myfiles, names, K)
##-----------------------------------##x
nimble.data <- readRDS("Data/Research_IPM_setup-data/R_IPM_24_Data/R_IPM_Nimble_data_setup_24_nimble.data.rds")
inits <- readRDS("Data/Research_IPM_setup-data/R_IPM_24_Data/R_IPM_Nimble_data_setup_24_inits.rds")
consts <- readRDS("Data/Research_IPM_setup-data/R_IPM_24_Data/R_IPM_Nimble_data_setup_24_consts.rds")
#############################################################X
##           Source model ----
#############################################################X
# source model
source("models/Individual_ipm-component_models/known_fate.R")
#############################################################X
##            Model run set up ----
#############################################################X
# Prepare model data
surv_data = list(
  # KF telemetered data 
  status = nimble.data$status,
  telem.juvenile = nimble.data$telem.juvenile, # adult is intercept
  telem.wmu = nimble.data$telem.wmu
)

# Constants
consts <- list(
  # KF telemetered constants 
  telem.nind = consts$telem.nind,  # 405
  telem.first = consts$telem.first,
  telem.last = consts$telem.last,
  telem.year.start = consts$telem.year.start,
  telem.year.end = consts$telem.year.end,
  female.telem.wmu = 4,
  male.n.wmu = 3
)

# Initial values
inits <- list(
  # KF telemetered inits
  # Intercepts and coefficients
  telem.beta.int = rnorm(1, 0, 0.5),       
  telem.beta.age = rnorm(1, 0, 0.5),       
  telem.beta.wmu = rnorm(4, 0, 1),        
  telem.beta.month = rnorm(12, 0, 1),

  # Variance parameters
  telem.sigma = runif(1, 0.1, 1),          
  telem.month.sigma = runif(1, 0, 2),
  
  # Survival probabilities 
  s.kf = array(0, dim = c(consts$telem.nind, 4, 12))
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
  monitors = c(
    # Original parameters
    # "telem.beta.int", "telem.beta.age", "telem.beta.wmu", "telem.beta.month",
    # "telem.sigma", "telem.month.sigma", "storage", 
    # 
    # Derived parameters
    "avg.ad.s.kf", "avg.juv.s.kf", "juv.male.adj"
               ),
  niter = iter, 
  nburnin = burn, 
  thin = t, 
  nchains = chain, 
  samplesAsCodaMCMC = TRUE,
  WAIC = F)
 beepr::beep(1)
gc()

#############################################################X
##           Model diagnostics ----
#############################################################X
### Check rhats ----
MCMCvis::MCMCsummary(survival_results, params = c("avg.ad.s.kf", "avg.juv.s.kf", "juv.male.adj"))
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

avg.ad.s.kf <- process_category_wmu(samples_df, "avg.ad.s.kf", "Female", "Adult")
avg.juv.s.kf <- process_category_wmu(samples_df, "avg.juv.s.kf", "Female", "Juvenile")
avg.juv.m.s.kf <- process_category_wmu(samples_df, "juv.male.adj", "Male", "Juvenile")


# Combine dfs
kf_survival_df <- bind_rows(
  avg.ad.s.kf, avg.juv.s.kf,
) %>% 
  mutate(wmu = factor(wmu, levels = c(1, 2, 3, 4), 
                       labels = c("2D", "3D", "4D", "5C")))


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
  labs(x = "Year", y = "Survival probability", shape = "", color = "") +
  
  # Keep consistent theme
  theme_classic() +
  theme(
    text = element_text(size = 16),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 16),
#   strip.text = element_blank(),      # Remove WMU labels
    legend.title = element_text(size = 16, face = "bold"),
    legend.text = element_text(size = 16),
    legend.position = "top",
    legend.key = element_rect(fill = "white"),  # Customize legend key appearance
    legend.key.size = unit(1.5, "lines")
  )

# Tables:
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