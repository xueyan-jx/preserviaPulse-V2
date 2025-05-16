## Purpose of script: Merge records from GBIF and BIEN
## Authors: GEOG 274
## Date: Spring, 2025
## Credits to: Yanni Zhan, Xue Yan

#-------------Setting up -------------
library(dplyr)
library(tidyr)  # Optional
library(ggplot2)
library(sf)
library(tigris) # to download county boundary
library(googledrive)
library(googlesheets4)
library(stringr)

#-------------1. Merge records from GBIF and BIEN -------------
# Read the final cleaned GBIF data and BIEN data
GBIF_occ_final <- read.csv("Replace to your path/Plant-GBIF0000607-250515123054153-final.csv",
                           quote = "\"",
                           stringsAsFactors = FALSE,
                           fileEncoding = "UTF-8")

BIEN_occ_final <- read.csv("E:/OneDrive/UCSB/Class/GEOG274/preserviaPulse/data/occurrences/BIEN0000607-250515123054153-final.csv",
                           quote = "\"",
                           stringsAsFactors = FALSE,
                           fileEncoding = "UTF-8")
# Check column names from the two datasets 
# and rename them so that they match each other
colnames(GBIF_occ_final)
colnames(BIEN_occ_final)
colnames(BIEN_occ_final)[1:2]=c("species","eventDate")

# Merge the two datasets
BIEN_occ_final$eventDate <- as.character(BIEN_occ_final$eventDate) #Make sure data types in eventDate are consistent
BIEN_GBIF <- bind_rows(BIEN_occ_final, GBIF_occ_final)

# View duplicates
BIEN_GBIF %>%
  group_by(species, geometry) %>%
  filter(n() > 1)

# Remove duplicates
BIEN_GBIF_unique <- BIEN_GBIF %>%
  distinct(species, geometry, .keep_all = TRUE)

# Check the records number for species
BIEN_GBIF_unique_stat <- BIEN_GBIF_unique %>%
  group_by(species) %>%
  summarize(n=n())

# Plot the results
ggplot(BIEN_GBIF_unique_stat, aes(x = reorder(species, -n), y = n)) +
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


# ----------- 2. Keep only one record for each grid -----------
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

# ---------(2) Convert merged data to spatial data ------
# Extract coordinate
BIEN_GBIF_unique <- BIEN_GBIF_unique %>%
  mutate(
    x = as.numeric(str_extract(geometry, "(?<=c\\()[0-9.]+")),
    y = as.numeric(str_extract(geometry, "(?<=, )[0-9.]+"))
  )

# Convert to spatial data
BIEN_GBIF_unique_sf <- st_as_sf(
  BIEN_GBIF_unique,
  coords = c("x", "y"),
  crs = 2229 # same crs to fishnet_clipped
)

# ---------(3) Reduce records ------
# Spatial join: assign each point to a fishnet polygon
BIEN_GBIF_with_grid <- st_join(BIEN_GBIF_unique_sf, fishnet_clipped)

# Remove points that do not fall into any fishnet cell (i.e., those with NA grid_id)
BIEN_GBIF_with_grid <- BIEN_GBIF_with_grid %>% filter(!is.na(grid_id))

# For each species and grid cell, keep only one record
BIEN_GBIF_final <- BIEN_GBIF_with_grid %>%
  group_by(species, grid_id) %>%
  slice(1) %>%  # Keep the first record in each group
  ungroup()

# save the cleaned dataframe
write_csv(BIEN_GBIF_final,
          file.path(occ_dir, paste0('BIEN_GBIF', download_id, "-final.csv")))

# Upload to google drive
# Get the target folder first to ensure it exists
speciesObs_folder <- drive_find(pattern = "speciesObs", type = "folder")
if(nrow(speciesObs_folder) > 0) {
  # Upload to Google Drive if folder found
  drive_upload(
    file.path(occ_dir, paste0('BIEN_GBIF', download_id, '-final.csv')),
    path = as_id(speciesObs_folder$id[1]),
    name = paste0('BIEN_GBIF', download_id, '-final.csv')
  )
} else {
  warning("Could not find 'speciesObs' folder in Google Drive. File saved locally only.")
}

# ----------- 3. Check final records number for species -----------
# Check final records number for all plants
BIEN_GBIF_unique_stat <- BIEN_GBIF_final %>%
  group_by(species) %>%
  summarize(n=n())
write.csv(BIEN_GBIF_unique_stat, 
          file.path(occ_dir, "Plant_occurrences_stat.csv"),
          row.names = FALSE)

# Check final records number for the first 20 species (Preserve most care)
top20_plants <- scientific_names[1:20]
BIEN_GBIF_unique_stat_top20 <- BIEN_GBIF_unique_stat %>%
  filter(species %in% top20_plants)

write.csv(BIEN_GBIF_unique_stat_top20, 
          file.path(occ_dir, "Plant_occurrences_stat_Top20.csv"),
          row.names = FALSE)

# Missing species
not_found <- setdiff(top20_plants, BIEN_GBIF_unique_stat$species)





