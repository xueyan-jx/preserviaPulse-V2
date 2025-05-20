## Purpose of script: Merge records from GBIF, BIEN, DP_portal, and Calflora
## Authors: GEOG 274
## Date: Spring, 2025
## Credits to: Yanni Zhan, Xue Yan

#-------------Setting up -------------
library(here)
library(dplyr)
library(tidyr)  # Optional
library(ggplot2)
library(sf)
library(terra)
library(tigris) # to download county boundary
library(googledrive)
library(googlesheets4)
library(stringr)

# Set a directory for data
here() # first check path

# Set a directory for data
occ_dir <- here("data/occurrences") # create a data folder with relative path
file.path(occ_dir) # double check its absolute path
if (!dir.exists(occ_dir)) dir.create(occ_dir) # if it's not there already, create it

# Authoritization to Google Drive
drive_auth()
gs4_auth(token = drive_token())

#--------------0. Define a function to clean species name---------------
# Function to remove taxonomic rank indicators such as ssp., var., subsp., etc.
# while keeping the genus, species, and any subspecies or variety names
# Only keep the first three words in the name
# Note this is not a perfect function to clean names 
# Just works for "species" column. For "scientificName" doesn't work due to complex name
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
    
    # Keep at most the first 3 parts
    parts <- parts[1:min(3, length(parts))]
    
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

#------------- 1. Get species list from google sheet  -------------
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

#-------------2. Merge records from GBIF, BIEN and DP_portal -------------
# Read the final cleaned GBIF data and BIEN data
if(!file.exists(here('data/CHELSA_bio1_1981-2010_V.2.1.tif'))){
  tmplt <- drive_get("CHELSA_bio1_1981-2010_V.2.1.tif")
  drive_download(tmplt, path=here('data/CHELSA_bio1_1981-2010_V.2.1.tif'))
}

# List of filenames to check and download if missing
files_to_check <- c(
  "GBIF-plants-cleaned-0515.csv",
  "BIEN-final.csv",
  "DP_cal-final.csv"
)

# Loop through each file to check if it exists locally; if not, download from Google Drive
for (fname in files_to_check) {
  local_path <- here("data", fname)  # Construct the full local path to the file
  
  if (!file.exists(local_path)) {
    # Try to locate the file on Google Drive by its exact name
    file_on_drive <- drive_get(fname)
    
    if (nrow(file_on_drive) > 0) {
      message(paste("Downloading", fname, "from Google Drive..."))
      drive_download(file_on_drive, path = local_path)
    } else {
      warning(paste("File", fname, "not found on Google Drive."))
    }
  } else {
    message(paste("File", fname, "already exists locally."))
  }
}

# Read the three data files using `here` for consistent paths
GBIF_occ_final <- read_csv(here("data", "GBIF-plants-cleaned-0515.csv"))
BIEN_occ_final <- read_csv(here("data", "BIEN-final.csv"))
DP_cal_occ_final <- read_csv(here("data", "DP_cal-final.csv"))

# Check column names from the three datasets 
# and rename them so that they match each other
colnames(GBIF_occ_final)
colnames(BIEN_occ_final)
colnames(DP_cal_occ_final)
colnames(BIEN_occ_final)[1:2]=c("species","eventDate")

# Merge the three datasets
BIEN_occ_final$eventDate <- as.character(BIEN_occ_final$eventDate) #Make sure data types in eventDate are consistent
DP_cal_occ_final$eventDate <- as.character(DP_cal_occ_final$eventDate)
GBIF_BIEN_DP_Cal <- bind_rows(GBIF_occ_final, BIEN_occ_final, DP_cal_occ_final)

glimpse(GBIF_BIEN_DP_Cal)

#-------------3. Organize the merged data -------------
# Due to the "scientificName" column from GBIF includes some subspecies that we didn't recognize in column "species"
# Here we need further organize the merged data for the 6 containing subspecies (within the 20 species that DP most care about)
# The code in this part might need to revise if we want to take a look at all species in our list (over 200 ish)

## "scientificName" column
# Visual check scientificName needed
unique(GBIF_BIEN_DP_Cal$scientificName)[order(unique(GBIF_BIEN_DP_Cal$scientificName))]

# Filter based on scientificName
#scitic_lst <- c("Deinandra increscens subsp. villosa (Tanowitz) B.G.Baldwin",
#                "Astragalus nuttallii var. nuttallii",
#                "Malacothrix saxatilis var. saxatilis", 
#                "Lilium humboldtii subsp. ocellatum (Kellogg) Thorne",
#                "Ribes amarum var. hoffmannii Munz",
#                "Horkelia cuneata var. sericea (A.Gray) Ertter & Reveal")

# Rename "species" column if scientificName includes subspecies in scitic_lit
GBIF_BIEN_DP_Cal <- GBIF_BIEN_DP_Cal %>%
  mutate(species = case_when(
    scientificName == "Deinandra increscens subsp. villosa (Tanowitz) B.G.Baldwin" ~ "Deinandra increscens villosa",
    scientificName == "Astragalus nuttallii var. nuttallii"                         ~ "Astragalus nuttallii nuttallii",
    scientificName == "Malacothrix saxatilis var. saxatilis"                        ~ "Malacothrix saxatilis saxatilis",
    scientificName == "Lilium humboldtii subsp. ocellatum (Kellogg) Thorne"        ~ "Lilium humboldtii ocellatum",
    scientificName == "Ribes amarum var. hoffmannii Munz"                           ~ "Ribes amarum hoffmannii",
    scientificName == "Horkelia cuneata var. sericea (A.Gray) Ertter & Reveal"      ~ "Horkelia cuneata sericea",
    TRUE ~ species  # Keep original if no match
  ))


# View duplicates
GBIF_BIEN_DP_Cal %>%
  group_by(species, geometry) %>%
  filter(n() > 1)

# Remove duplicates
GBIF_BIEN_DP_Cal_unique <- GBIF_BIEN_DP_Cal %>%
  distinct(species, geometry, .keep_all = TRUE)

# Check the records number for species
GBIF_BIEN_DP_Cal_stat <- GBIF_BIEN_DP_Cal_unique %>%
  group_by(species) %>%
  summarize(n=n())

# Plot the results
ggplot(GBIF_BIEN_DP_Cal_stat, aes(x = reorder(species, -n), y = n)) +
  geom_bar(stat = "identity") +
  labs(
    title = "Species Occurrence Counts",
    x = "Species",
    y = "Number of Records"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    text = element_text(size = 14)
  )


# ----------- 4. Keep only one record for each grid -----------
#---------(1) Create fishnet using climate data (CHELSA_bio1) ------
if(!file.exists(here('data/CHELSA_bio1_1981-2010_V.2.1.tif'))){
  tmplt <- drive_get("CHELSA_bio1_1981-2010_V.2.1.tif")
  drive_download(tmplt, path=here('data/CHELSA_bio1_1981-2010_V.2.1.tif'))
}

raster_template <- terra::rast(here("data/CHELSA_bio1_1981-2010_V.2.1.tif"))

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

# ---------(2) Convert merged data to spatial data ------
# Extract coordinate
GBIF_BIEN_DP_Cal_unique <- GBIF_BIEN_DP_Cal_unique %>%
  mutate(
    x = as.numeric(str_extract(geometry, "(?<=c\\()[0-9.]+")),
    y = as.numeric(str_extract(geometry, "(?<=, )[0-9.]+"))
  )

# Convert to spatial data
GBIF_BIEN_DP_Cal_unique_sf <- st_as_sf(
  GBIF_BIEN_DP_Cal_unique,
  coords = c("x", "y"),
  crs = 2229 # same crs to fishnet_clipped
)

# ---------(3) Reduce records ------
# Spatial join: assign each point to a fishnet polygon
GBIF_BIEN_DP_Cal_grid <- st_join(GBIF_BIEN_DP_Cal_unique_sf, fishnet_clipped)

# Remove points that do not fall into any fishnet cell (i.e., those with NA grid_id)
GBIF_BIEN_DP_Cal_grid <- GBIF_BIEN_DP_Cal_grid %>% filter(!is.na(grid_id.x))

# For each species and grid cell, keep only one record
GBIF_BIEN_DP_Cal_final <- GBIF_BIEN_DP_Cal_grid %>%
  group_by(species, grid_id.x) %>%
  slice(1) %>%  # Keep the first record in each group
  ungroup()


# Clean species name
# Just found that some missed species names can be found in the column "scientificName"
# This part could be improved in the future once we found a better way to clean species names
# Current remove_rank_indicators_df function doesn't work
# clean_names_finaldata <- remove_rank_indicators_df(GBIF_BIEN_DP_final$scientificName)

# Add the cleaned name to GBIF_BIEN_DP_final
#GBIF_BIEN_DP_final <- GBIF_BIEN_DP_final %>%
#  mutate(cleaned_name = clean_names_finaldata$cleaned_name)

# save the cleaned dataframe
write_csv(GBIF_BIEN_DP_Cal_final,
          file.path(occ_dir, 'GBIF_BIEN_DP_Cal_final.csv'))

# Upload to google drive
# Get the target folder first to ensure it exists
speciesObs_folder <- drive_find(pattern = "speciesObs", type = "folder")
if(nrow(speciesObs_folder) > 0) {
  # Upload to Google Drive if folder found
  drive_upload(
    file.path(occ_dir, 'GBIF_BIEN_DP_Cal_final.csv'),
    path = as_id(speciesObs_folder$id[1]),
    name = 'GBIF_BIEN_DP_Cal_final.csv'
  )
} else {
  warning("Could not find 'speciesObs' folder in Google Drive. File saved locally only.")
}

# ----------- 4. Check final records number for species -----------
# ------------(1) Check final records number for all plants------
GBIF_BIEN_DP_Cal_final_stat <- GBIF_BIEN_DP_Cal_final %>%
  group_by(species) %>%
  summarize(n=n())
write.csv(GBIF_BIEN_DP_Cal_final_stat, 
          file.path(occ_dir, "Plant_occurrences_stat.csv"),
          row.names = FALSE)

# ------------(2) Check final records number for the first 20 species (Preserve most care)-----
# Convert scientific_names to dataframe
scientific_names_df <- data.frame(scientificName = scientific_names, stringsAsFactors = FALSE)

# Clean the list of scientific names for comparison
clean_names_list <- remove_rank_indicators_df(scientific_names_df$scientificName)

# Get the top 
top20_plants <- clean_names_list$cleaned_name[1:20]
GBIF_BIEN_DP_Cal_final_stat_top20 <- GBIF_BIEN_DP_Cal_final_stat %>%
  filter(species %in% top20_plants)

write.csv(GBIF_BIEN_DP_Cal_final_stat_top20, 
          file.path(occ_dir, "Plant_occurrences_stat_Top20.csv"),
          row.names = FALSE)

# Missing species
not_found <- setdiff(top20_plants, GBIF_BIEN_DP_Cal_final_stat$species)





