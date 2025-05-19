## Purpose of script: grab integrated_occurrences_dangermond_Portal data
## Authors: GEOG 274
## Date: Spring, 2025
## Credits to: Wenxin Yang, Jacqueline Vogel, Yanni Zhan, Xue Yan

# ------------- Setting up --------------
# Install and load packages

library(here)
library(sf)
library(terra) # read raster
library(dplyr)
library(rgbif)
library(googledrive)
library(googlesheets4)
library(tigris) # to download county boundary
library(tidyr)
library(ggplot2)
library(CoordinateCleaner) # add spatial filters
library(readr) # for using write_csv


# Set a directory for data
here() # first check path
occ_dir <- here("data/occurrences") # create a data folder with relative path
file.path(occ_dir) # double check its absolute path
if (!dir.exists(occ_dir)) dir.create(occ_dir) # if it's not there already, create it

# Authorization to Google Drive
drive_auth()
gs4_auth(token = drive_token())

# Get boundaries for areas of interest
ca_counties <- counties(state = "CA", cb = TRUE) # gather CA county boundary file
target_counties <- ca_counties %>%
  filter(NAME %in% c("Santa Barbara", "Ventura", "San Luis Obispo")) # limit to counties of interest

# ------------- Get occurrence info for a species list --------------
## ------------ 1. Get species list ------------
ss <- drive_get("Special Status Species")
all_names <- as.data.frame(matrix(ncol=2, nrow=0))
colnames(all_names) <- c('Name (Latin)', 'Name (common)')
for(taxon in c('Plants', 'Birds', 'Mammals', 'Herps', 'Inverterbrates')){
  dat <- read_sheet(ss, sheet=taxon) %>% dplyr::select(`Name (Latin)`, `Name (common)`)
  all_names <- rbind(all_names, dat)
}

all_names <- all_names %>% mutate(
  `Name (Latin)` = ifelse(`Name (Latin)`=='Vulpes vulpes ssp.', 'Vulpes vulpes', `Name (Latin)`)
) %>% unique()

scientific_names <- unique(all_names$`Name (Latin)`)

## ------------ 2. Get information from dangermond portal ------------
# download from google drive
dp_gdrive <- drive_find(pattern = "integrated_occurrences_dangermond.csv", type = "csv")
drive_download(dp_gdrive[1,]$id, path = here('data/occurrences/integrated_occurrences_dangermond.csv'))
DP_portal_occ <- read.csv("data/occurrences/integrated_occurrences_dangermond.csv",
                           quote = "\"",
                           stringsAsFactors = FALSE,
                           fileEncoding = "UTF-8")

# Xue's Comments: in this file, column "scientificName" includes species and subspecies
# Names in column "species" includes species (higher hierarchy) or has error (e.g., No. 41: Equisetum telmateia braunii is not the subspecies of Equisetum braunii)
# For matching the name with GBIF and BIEN, here we use "scientificName" and revise the column names
colnames(DP_portal_occ)[colnames(DP_portal_occ) == "species"] <- "species_name_with_error"
colnames(DP_portal_occ)[colnames(DP_portal_occ) == "scientificName"]=c("species")
# Exclude NA 
DP_portal_occ <- DP_portal_occ %>% 
  filter(!is.na(species))

#---------------3. Exclude records from GBIF and iNaturalist and non-plant record-------------
# change the syntax to include or exclude plants
#Plantae for plants, Animalia for animals
unique(DP_portal_occ$kingdom)
DP_portal_occ <- DP_portal_occ %>% 
  filter(source!="GBIF" & source!="iNaturalist" & kingdom =="Animalia")


# ------------- 4. Exclude records before 1980----------
DP_portal_occ <- DP_portal_occ %>% filter(year >=1980)
unique(DP_portal_occ$year)

# ------------ 5. Exclude species with no records -----------
# Function to remove taxonomic rank indicators such as ssp., var., subsp., etc.
# while keeping the genus, species, and any subspecies or variety names.
remove_rank_indicators_df <- function(name_vec) {
  cleaned_names <- sapply(name_vec, function(name) {
    # Convert the whole name to lowercase for consistent matching
    name_lower <- tolower(name)
    # Split by whitespace
    parts <- unlist(strsplit(name_lower, "\\s+"))
    # Define all rank indicators to remove (with and without dot)
    ranks <- c("ssp", "ssp.", "subsp", "subsp.", "var", "var.", "subvar", "subvar.", "f", "f.", "cv", "cv.")
    # Remove the rank indicators
    parts <- parts[!parts %in% ranks]
    
    # Capitalize first letter of the first word if exists
    if (length(parts) > 0) {
      parts[1] <- paste0(toupper(substring(parts[1], 1, 1)), substring(parts[1], 2))
    }
    
    # Join the remaining parts back together with space
    cleaned <- paste(parts, collapse = " ")
    return(cleaned)
  })
  
  # return data.frame
  return(data.frame(original_name = name_vec, cleaned_name = cleaned_names, stringsAsFactors = FALSE))
}

# Convert scientific_names to datapframe
scientific_names_df <- data.frame(scientificName = scientific_names, stringsAsFactors = FALSE)

# Clean the list of scientific names for comparison
clean_names <- remove_rank_indicators_df(scientific_names_df$scientificName)

# Find unique cleaned species names in the comparison list
unique_species <- unique(clean_names)

# Identify species in 'unique_species' that are missing from the occurrence data
missing_from_list <- setdiff(unique_species$cleaned_name, DP_portal_occ$species)
print(missing_from_list)

# Filter data to species with records
DP_portal_filtered <- DP_portal_occ %>%
  filter(tolower(species) %in% tolower(unique_species$cleaned_name))

# Find missing species
#this_info <- gbif_names[gbif_names$usageKey %in% missing_from_list,]
#this_info$`Handling person` = myname
#print(unique(this_info$canonicalName))

#missing_species <- rbind(missing_species, this_info)

# Update the Google spreadsheet
#write_sheet(missing_species, pt_info, sheet = "GBIF missing species")

# ------------- 6. Remove duplicates and rows with invalid coordinates----------
DP_portal_filtered <- DP_portal_filtered %>%
  filter(!is.na(longitude) & !is.na(latitude)) %>%   
  group_by(species, longitude, latitude, eventDate, year) %>%
  filter(row_number() == 1) %>%
  ungroup()

# ------------- 7. Other filter ----------
# CoordinateCleaner Wrapper
DP_portal_filtered <- DP_portal_filtered %>% mutate(record_id = row_number())
flag_cols_keep <- c('record_id', '.val', '.inst') #flag problems

flags_DP <- clean_coordinates(x = DP_portal_filtered, 
                                lon = "longitude", 
                                lat = "latitude", 
                                species = "species", 
                                tests = c("institutions", "centroids", "zeros"))[flag_cols_keep]

# remove rows with invalid coordinates
ids_DP_invalid <- flags_DP[flags_DP$`.val` == FALSE, ]$record_id
DP_occ_clean <- DP_portal_filtered %>% filter(!record_id %in% ids_DP_invalid)

# remove rows recorded in biodiversity institutions
ids_DP_inst <- flags_DP[flags_DP$`.inst` == FALSE, ]$record_id
DP_occ_clean <- DP_occ_clean %>% filter(!record_id %in% ids_DP_inst)

# Get datasource and summarize the information
unique(DP_occ_clean$source) 

DP_occ_clean %>% group_by(source) %>%
  summarize(n=n())

## Remove records with uncertainty >1000 meters
DP_occ_clean <- DP_occ_clean %>%
  filter(coordinateUncertaintyInMeters <=1000 | is.na(coordinateUncertaintyInMeters))

write.table(DP_occ_clean, here('data/occurrences/animals/DP_clean_animals_0519.csv'), sep= ';')

# ----------- 8. Keep only one record for each grid -----------
#---------(1) Create fishnet using climate data (CHELSA_bio1) ------
if(!file.exists(here('data/CHELSA_bio1_1981-2010_V.2.1.tif'))){
  tmplt <- drive_get("CHELSA_bio1_1981-2010_V.2.1.tif")
  drive_download(tmplt, path=here('data/CHELSA_bio1_1981-2010_V.2.1.tif'))
}

raster_template <- terra::rast(here("data/CHELSA_bio1_1981-2010_V.2.1.tif"))

# County boundary
target_counties <- target_counties

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

# ---------(2) Convert cleaned DP data to spatial data ------
# Convert to spatial data
DP_occ_clean_sf <- st_as_sf(DP_occ_clean, coords = c("longitude", "latitude"), crs = 4326)

# Reproject to match the fishnet CRS
DP_occ_clean_sf <- st_transform(DP_occ_clean_sf, st_crs(fishnet_clipped))

# ---------(3) Reduce records ------
# Spatial join: assign each point to a fishnet polygon
DP_occ_clean_with_grid <- st_join(DP_occ_clean_sf, fishnet_clipped)

# Remove points that do not fall into any fishnet cell (i.e., those with NA grid_id)
DP_occ_clean_with_grid <- DP_occ_clean_with_grid %>% filter(!is.na(grid_id))

# For each species and grid cell, keep only one record
DP_occ_unique <- DP_occ_clean_with_grid %>%
  group_by(species, grid_id) %>%
  slice(1) %>%  # Keep the first record in each group
  ungroup()

head(DP_occ_unique$geometry)
# save the cleaned dataframe
write_csv(DP_occ_unique,
          file.path(here(occ_dir, 'animals'), 'DP_portal_animals_cleaned.csv'))

# Upload to google drive
# Get the target folder first to ensure it exists
speciesObs_folder <- drive_find(pattern = "speciesObs", type = "folder")
if(nrow(speciesObs_folder) > 0) {
  # Upload to Google Drive if folder found
  drive_upload(
    file.path(occ_dir, 'DP-Portal-final.csv'),
    path = as_id(speciesObs_folder$id[1]),
    name = 'DP-Portal-final.csv'
  )
} else {
  warning("Could not find 'speciesObs' folder in Google Drive. File saved locally only.")
}
