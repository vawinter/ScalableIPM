###############################################################################X
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
############## Research Integrated Population Model (R_IPM): ##################X
#                     #---# PhD Dissertation: R_IPM #---#
#        Creating a Bayesian IPM to inform turkey management in PA
###                                                                         ###X
#               Data formatting for female known-fate model (KF)
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #X
###############################################################################X

rm(list = ls())
gc()

# Source scripts of functions and data preparation
source("Analysis/00_IPM_funs.R")
##-----------------------------------------------------#X
# location of Turkey database
access.dir <- "C:/Users/vaw5154/OneDrive - The Pennsylvania State University/PhD/PSUTurkey/PA Turkey Database/"
# # Read in tables from Access database
db <- file.path(paste0(access.dir,"TurkeyDB_10_10_2025.accdb"))
ch <- odbcConnectAccess2007(db) # open channel to database
dcap <- sqlFetch(ch,"captures", as.is = T)
dcen <- sqlFetch(ch,"censors", as.is = T)
dmor <- sqlFetch(ch,"mortalities", as.is = T)
dtag <- sqlFetch(ch,"transmitter_status", as.is = T)
close(ch) #close channel to database
##-----------------------------------------------------#X
# DF structuring ----
##-----------------------------------------------------#X
# Grab necessary columns from dfs
capt <- dcap[,c("bandid","captyr", "captmo", "captday", "studyarea",
                "lat","long","age","sex","weight","birdcondition", "recapture")]
cens <- dcen[,c("bandid","cenyr","cenmo","cenday")]
mor <- dmor[,c("bandid","estmortyr","estmortmo","estmortday","fate_m","subfate1_m")]

# Merge df's
d1 <- merge(capt, cens, by="bandid", all.x=TRUE, all.y=FALSE)
df <- merge(d1, mor, by="bandid", all.x=TRUE, all.y=FALSE)  

# Filter out rows where capyr is not 2024
df <- df[df$captyr != 2024, ]
df <- df[df$captyr != 2025, ]

# Filter for rows where sex is "F"
df <- df[df$sex == "F", ]
table(df$captyr)

# Define year we are evaluating until
stop_year <- 2024
##-----------------------------------------------------#X
# Create encounter histories ----
##-----------------------------------------------------#X
# w. function
status <- encounter_histories(df, filter_sex = "F", stop_year = 2024)
status_matrix <- status # created

##-----------------------------------------------------#X
# Create age variables, month of entry, etc ----
##-----------------------------------------------------#X
# Hen information (duration in study, location, age indicator)
# Filter and prep data
hen <- df %>%
  mutate(
    begin = as.Date(paste0(captyr, "-", captmo, "-", captday)),
    censor = as.Date(paste0(cenyr, "-", cenmo, "-", cenday), format = "%Y-%m-%d"),
    mort = as.Date(paste0(estmortyr, "-", estmortmo, "-", estmortday), format = "%Y-%m-%d"),
    fate = case_when(
      !is.na(mort) ~ "Dead",
      !is.na(censor) ~ "Censor",
      year(censor) > as.Date(paste0(stop_year, "-12-31")) ~ "Alive",
      year(mort) > as.Date(paste0(stop_year, "-12-31")) ~ "Alive",
      TRUE ~ "Alive"
    ),
    end = case_when(
      fate == "Alive" ~ as.Date(paste0(stop_year, "-12-31")),
      fate == "Censor" ~ censor,
      fate == "Dead" ~ mort
    )
  ) 

# Study period
start_year <- year(min(hen$begin))
end_year <- stop_year
telem.nyears <- end_year - start_year + 1

# Calculate month indices
hen_mo <- hen %>%
  mutate(
    start = month(begin) + 12 * (year(begin) - start_year),
    stop = month(end) + 12 * (year(end) - start_year),
    year_start = year(begin) - start_year + 1
  )

# Initialize arrays
telem.nind <- nrow(hen_mo)
is_juvenile_matrix <- array(0, dim = c(telem.nind, telem.nyears, 12))

# Fill arrays:
# Loop over each individual, study year, month (Jan=1, Dec=12)
for (i in 1:telem.nind) {
  for (y in 1:telem.nyears) {
    for (m in 1:12) {
      # Calculate absolute month index since study start
      month_index <- (y - 1) * 12 + m
      # Check if individual was being monitored during this month
      if (month_index >= hen_mo$start[i] && month_index <= hen_mo$stop[i]) {
        # Mark as juvenile if: captured as juvenile AND currently in capture year AND Jan-May
        is_juvenile_matrix[i, y, m] <- (hen_mo$age[i] == "J" && 
                                          y == hen_mo$year_start[i] && 
                                          m <= 5) # Juvenieles age to Audlt in June
      }
    } # end m
  } # end y
} # end i

# WMU indicator
telem.wmu = model.matrix(~ -1 + studyarea, data = df)

# When individual entered/left study (month)
hen_timing <- hen_mo %>%
  group_by(bandid) %>%
  mutate(
    year_start = year(begin) - start_year + 1,
    year_end = case_when(
      year(end) < stop_year ~ year(end) - start_year + 1,
      TRUE ~ stop_year - start_year + 1  # Default case
    ),
    first = month(begin),
    last = month(end),
    start = first + 12 * (year_start - 1),
    stop = last + 12 * (year_end - 1)
  ) %>%
  ungroup()

# Define variables for saving
telem.first <- hen_timing$first
telem.last <- hen_timing$last
telem.year.start = hen_timing$year_start
telem.year.end = hen_timing$year_end

##-----------------------------------------------------#X
# Save data ----
##-----------------------------------------------------#X
# Create a list to store the variables
kf_data <- list(
  status_matrix = status_matrix,
  is_juvenile_matrix2 = is_juvenile_matrix,
  telem.first = telem.first,
  telem.last = telem.last,
  telem.wmu = telem.wmu,
  telem.nind = telem.nind,
  telem.wmu = telem.wmu,
  telem.year.start = telem.year.start,
  telem.year.end = telem.year.end
)

# Define a directory where you want to save the RDS files
output_dir <- "Data/Research_IPM_setup-data/kf_data_22-23/"

# Loop through the list and save each variable as a separate RDS file
for (var_name in names(kf_data)) {
  variable <- kf_data[[var_name]]
  file_name <- paste(output_dir, var_name, ".rds", sep = "")
  saveRDS(variable, file_name)
}
# Done!
