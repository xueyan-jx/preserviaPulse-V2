## Purpose of script: grab BIEN data
## Authors: GEOG 274
## Date: Spring, 2025
## Credits to: Yanni Zhan, Xue Yan, Dr. Diyang Cui

# ------------- Setting up --------------
# Install and load packages
library(here)
library(googledrive)
library(googlesheets4)
library(BIEN)
library(stringr) # for "fuzzy matching" and filtering our target species 
library(ggplot2)
library(sf)
library(terra)
library(tigris) # to download county boundary
library(dplyr)
library(CoordinateCleaner)

# Set a directory for data
here() # first check path

# Set a directory for data
occ_dir <- here("data/occurrences") # create a data folder with relative path
file.path(occ_dir) # double check its absolute path
if (!dir.exists(occ_dir)) dir.create(occ_dir) # if it's not there already, create it

# Authoritization to Google Drive
drive_auth()
gs4_auth(token = drive_token())

#------------- 1. Get species list  -------------
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


#-------------2. Get data from BIEN -------------
# Get boundaries for areas of interest
country_vector <- rep("United States", 3)
state_vector <- rep("California", 3)
county_vector <- c("Santa Barbara", "Ventura", "San Luis Obispo")

# Get all records in the three counties
BIEN_occ <- BIEN_occurrence_county(
  country = country_vector,
  state = state_vector,
  county = county_vector)

# Filter to our target species
# NOTE: This is a loose match that does not account for specific cases such as synonyms.
# Logic: if the first two words match, then assume two names are the same.
filtered_occ <- BIEN_occ %>% 
  # Remove rows with NA in latitude or longitude
  filter(!is.na(latitude) & !is.na(longitude)) %>% 
  # Split each species name and get the first two words
  # Check if any name in scientific_names contains the same first two words
  filter(tolower(word(scrubbed_species_binomial, 1, 2)) %in% 
           tolower(word(scientific_names, 1, 2)))


# Save out the downloaded data
raw_fname <- file.path(occ_dir, "BIEN_3counties_occurrences.csv")
write.csv(filtered_occ, raw_fname, row.names = FALSE)

# Upload data to google drive
# Get the target folder first to ensure it exists
speciesObs_folder <- drive_find(pattern = "speciesObs", type = "folder")
if(nrow(speciesObs_folder) > 0) {
  # Upload to Google Drive if folder found
  drive_upload(
    raw_fname,
    path = as_id(speciesObs_folder$id[1]),
    name = basename(raw_fname)
  )
} else {
  warning("Could not find 'specieObs' folder in Google Drive. File saved locally only.")
}

#-------------3. Data clean  -------------
BIEN_occ <- filtered_occ

# Or read the saved BIEN data directly
# BIEN_occ<-read.csv("replace to your path/BIEN_3counties_occurrences.csv")

# ------------- (1) Delete records before 1980----------
BIEN_occ <- BIEN_occ %>%
  mutate(date_collected = as.Date(date_collected)) %>%           
  filter(date_collected >= as.Date("1980-01-01") & date_collected <= as.Date("2025-05-17")) 

# ------------- (2) Other filter ----------
# CoordinateCleaner Wrapper
BIEN_occ <- BIEN_occ %>% mutate(record_id = row_number())
flag_cols_keep <- c('record_id', '.val', '.inst') #flag problems

flags_bien <- clean_coordinates(x = BIEN_occ, 
                                lon = "longitude", 
                                lat = "latitude", 
                                species = "scrubbed_species_binomial", 
                                tests = c("institutions", "centroids", "zeros"))[flag_cols_keep]

# remove rows with invalid coordinates
ids_bien_invalid <- flags_bien[flags_bien$`.val` == FALSE, ]$record_id
BIEN_occ_clean <- BIEN_occ %>% filter(!record_id %in% ids_bien_invalid)

# remove rows recorded in biodiversity institutions
ids_bien_inst <- flags_bien[flags_bien$`.inst` == FALSE, ]$record_id
BIEN_occ_clean <- BIEN_occ_clean %>% filter(!record_id %in% ids_bien_inst)

# Get datasource and summarize the information
unique(BIEN_occ_clean$datasource)

BIEN_occ_clean %>% group_by(datasource) %>%
  summarize(n=n())

#-------------4. Exclude records from GBIF -------------
BIEN_occ_subset <- BIEN_occ_clean %>% filter(datasource!="GBIF")


# ----------- 5. Keep only one record for each grid -----------
#---------(1) Create fishnet using climate data (CHELSA_bio1) ------
raster_template <- terra::rast("Change this to your path/CHELSA_bio1_1981-2010_V.2.1.tif")

# County boundary
target_counties <- ca_counties %>%
  filter(NAME %in% c("Santa Barbara", "Ventura", "San Luis Obispo")) # limit to counties of interest

# Keep sameCoordinate Reference Systems
counties <- st_transform(target_counties, sf::st_crs(raster_template))

# Get resolution and extent of the raster data
res_xy <- res(raster_template)
raster_bbox <- st_as_sfc(st_bbox(raster_template))

# Create full fishnet within the target extent
full_fishnet <- st_make_grid(raster_bbox, cellsize = res_xy, what = "polygons", square = TRUE)
full_fishnet_sf <- st_sf(full_fishnet)

# Crop and mask raster to the extent of the counties
fishnet_clipped <- st_intersection(full_fishnet_sf, counties)

# Add a unique ID column to the clipped fishnet polygons if it doesn't exist
if (!"grid_id" %in% colnames(fishnet_clipped)) {
  fishnet_clipped$grid_id <- 1:nrow(fishnet_clipped)
}

# ---------(2) Convert cleaned BIEN data to spatial data ------
# Convert to spatial data
BIEN_occ_subset_sf <- st_as_sf(BIEN_occ_subset, coords = c("longitude", "latitude"), crs = 4326)

# Reproject to match the fishnet CRS
BIEN_occ_subset_sf <- st_transform(BIEN_occ_subset_sf, st_crs(fishnet_clipped))

#### ---------c. Reduce records ------
# Spatial join: assign each point to a fishnet polygon
BIEN_occ_with_grid <- st_join(BIEN_occ_subset_sf, fishnet_clipped)

# Remove points that do not fall into any fishnet cell (i.e., those with NA grid_id)
BIEN_occ_with_grid <- BIEN_occ_with_grid %>% filter(!is.na(grid_id))

# Rename the column for species
BIEN_occ_with_grid <- BIEN_occ_with_grid %>%
  rename(species = scrubbed_species_binomial)

# For each species and grid cell, keep only one record
BIEN_occ_unique <- BIEN_occ_with_grid %>%
  group_by(species, grid_id) %>%
  slice(1) %>%  # Keep the first record in each group
  ungroup()

# save the cleaned dataframe
write_csv(BIEN_occ_unique,
          file.path(occ_dir, 'BIEN-final.csv'))

# Upload to google drive
# Get the target folder first to ensure it exists
speciesObs_folder <- drive_find(pattern = "speciesObs", type = "folder")
if(nrow(speciesObs_folder) > 0) {
  # Upload to Google Drive if folder found
  drive_upload(
    file.path(occ_dir, 'BIEN-final.csv'),
    path = as_id(speciesObs_folder$id[1]),
    name = 'BIEN-final.csv'
  )
} else {
  warning("Could not find 'speciesObs' folder in Google Drive. File saved locally only.")
}
