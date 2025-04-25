################       1. Get species list   ########################
library(googledrive)
library(googlesheets4)

# Authoritization
drive_auth()
gs4_auth(token = drive_token())

# Get data
ss <- drive_get("Special Status Species")

## E.g. sheet names
ss_meta <- gs4_get(ss) 
sheet_names <- ss_meta$sheets$name

dat <- read_sheet(ss, sheet = sheet_names[1])
head(dat)

# Get species' names
scientific_names_df <- dat[, c("Name (common)", "Name (Latin)")]

# Use Latin name for searching
scientific_names <- scientific_names_df$`Name (Latin)`

################       2. Get data from BIEN   ########################
# `BIEN_occurrence_county` Returns all occurrences records within a given state/province
install.packages("RBIEN")
library(BIEN)
library(stringr) # for "fuzzy matching" and filtering our target species 

country_vector<-c("United States","United States","United States")
state_vector<-c("California","California","California")
county_vector<-c("Santa Barbara","Ventura","San Luis Obispo")

# Get all records in the three counties
BIEN_occ<-BIEN_occurrence_county(country=country_vector, state = state_vector, county = county_vector)

# Filter to our target species
filtered_occ <- BIEN_occ %>%
  filter(
    sapply(scrubbed_species_binomial, function(x) {
      # Split each species name and get the first two words
      scrubbed_first_two <- str_split(x, " ", simplify = TRUE)[1:2]
      # Create a regular expression to match the first two words
      regex <- paste(scrubbed_first_two, collapse = " ")
      # Check if any name in scientific_names contains the same first two words
      any(str_detect(scientific_names, regex))
    }) &
      !is.na(latitude) & !is.na(longitude)  # Remove rows with NA in latitude or longitude
  )

write.csv(filtered_occ, "E:/path/to/your/folder/BIEN_filtered_occurrences.csv", row.names = FALSE)


