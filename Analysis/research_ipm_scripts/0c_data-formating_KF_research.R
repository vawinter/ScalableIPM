###############################################################################X
# Data formatting Known fate ----
###############################################################################X
rm(list = ls())
gc()

# Source scripts of functions and data preparation
source("Analysis/00_IPM_funs.R")

# location of Turkey database
access.dir <- "TurkeyDatabase/"

# Read in tables from Access database
db <- file.path(paste0(access.dir,"TurkeyDB.accdb"))
ch <- odbcConnectAccess2007(db) # open channel to database
dcap <- sqlFetch(ch,"captures", as.is = T)
dcen <- sqlFetch(ch,"censors", as.is = T)
dmor <- sqlFetch(ch,"mortalities", as.is = T)
dtag <- sqlFetch(ch,"transmitter_status", as.is = T)
close(ch) #close channel to database

# Grab necessary columns from dfs
capt <- dcap[,c(1:4,6:8,13:16,34)]
cens <- dcen[,c(1:4)]
mor <- dmor[,c(1,12:14, 21, 22)]

# Merge df's
d1 <- merge(capt, cens, by="bandid", all.x=TRUE, all.y=FALSE)
df <- merge(d1, mor, by="bandid", all.x=TRUE, all.y=FALSE)  

# Filter out rows where capyr is not 2024
df <- df[df$captyr != 2024, ]
#df <- df[df$studyarea != "2D", ]

# Filter for rows where sex is "F"
df <- df[df$sex == "F", ]

# Create encounter histories w. function
status <- encounter_histories(df, filter_sex = "F")
dat <- hen_histories(df, filter_sex = "F")

## Format the data and initial values 
unique_bands <- unique(dat$hen$bandid)
telem.nind <- length(unique_bands)

## Format length of months vector
nMonths = max(dat$hen$last)

# Initialize matrices
status_matrix <- status
# is_adult_matrix <- matrix(dat$is.adult, nrow = nind, ncol = nMonths)
is_juvenile_matrix <- dat$is_juvenile_matrix

telem.first <- dat$hen$first
telem.last <- dat$hen$last
telem.wmu = model.matrix(~ -1 + is.2d + is.3d + is.4d + is.5c, data = dat$hen)

# Create a list to store the variables
kf_data <- list(
  status_matrix = status_matrix,
  is_juvenile_matrix2 = is_juvenile_matrix,
  telem.first = telem.first,
  telem.last = telem.last,
  telem.wmu = telem.wmu,
  telem.nind = telem.nind,
  telem.wmu = telem.wmu,
  telem.month = dat$telem_month,
  telem.year.start = dat$hen$year_start,
  telem.year.end = dat$hen$year_end
)

# Define a directory where you want to save the RDS files
output_dir <- "Data/Research_IPM_setup-data/"

# Loop through the list and save each variable as a separate RDS file
for (var_name in names(kf_data)) {
  variable <- kf_data[[var_name]]
  file_name <- paste(output_dir, var_name, ".rds", sep = "")
  saveRDS(variable, file_name)
}
