## Purpose of Script: Process and format CNDDB Rare Occurrence Data
## Authors: GEOG 274
## Date: Spring, 2025 
## Credits to: Olivia Ross, Wenxin Yang, Jacqueline Vogel, 
##             Yanni Zhan, Xue Yan, Dr. Lei Song

# packages 
librarian::shelf(sf, dplyr, here, terra, lubridate, googledrive, tigris, googlesheets4)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
# READING IN RARE OCCURENCE DATA & PROJECTING TO PROJECT CRS
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#  
# listing our occurrence data from CNDDB folder
occurrences <- list.files(here("storedlocation"), pattern = "\\.csv$", full.names = TRUE)

# creating empty list to store processed data
occurrence.list <- list()

# reading in occurrence df, projecting to proj crs, and combining in a list
for(df in occurrences) {
  
  occurrence.df <- read.csv(df)
  
  occurrence.sf <- st_as_sf(occurrence.df, coords = c("Longitude", "Latitude"), crs = 4326)
  
  occurrence.reproj <- st_transform(occurrence.sf, crs = 2229)
  
  occurrence.list[[df]] <- occurrence.reproj
  
}

# now binding all rows 
all <- do.call(rbind, occurrence.list)

crs(all)

# grabbing coordinate columns
all_coords <- st_coordinates(all)

all_added_coords <- cbind(all, all_coords)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
# DATA CLEANING

# adopting data cleaning steps and spatial filters from script 1

# 1. removed species before 1980
# 2. removed duplicates 
# 3. removed sea coordinates
# 4. checking accuracy of points and removing > 1000 meters

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
# 1. grabbing dates to make sure we are not including observations earlier than 1980
all_added_coords <- all_added_coords |>
  mutate(date_raw = gsub("AM", " AM", gsub("PM", " PM", LastUpdate)),
         parsed_date = mdy_hms(date_raw),
         Year = year(parsed_date)) |>
  dplyr::filter(Year > 1980)

# now adding 2 occurrences of the northern harrier, from eBird
nh <- data.frame(SciName = c("Circus hudsonius", "Circus hudsonius"),
                 ComName = c("northern harrier", "northern harrier"),
                 TaxonGroup = c("Birds", "Birds"),
                 Latitude = c(34.47103, 34.90286),
                 Longitude = c(-120.46648, -120.63444),
                 Year = c(2025, 2024))
# projecting df 
nh <- st_as_sf(nh, coords = c("Longitude", "Latitude"), crs = 4326)
nh <- st_transform(nh, 2229)

# grabbing coord columns and adding back to match other df 
nh_coords <- st_coordinates(nh)

nh <- cbind(nh, nh_coords)

# binding to rest of data 
all_added_coords <- bind_rows(all_added_coords, nh)
  
# 2. making sure there are no duplicates - produced df of 0 (no duplicates)
duplicates <- all_added_coords |>
  group_by(SciName, X, Y, Year) |>
  filter(n() > 1) |>
  ungroup()

# 3. clipping to target counties in case there are any sea coordinates
# reloading tri-counties
counties <- counties(state = "CA", cb = TRUE) 

target_counties <- counties |>
  filter(NAME %in% c("Santa Barbara", "Ventura", "San Luis Obispo"))

# transforming to proj crs
target_counties <- st_transform(target_counties, sf::st_crs(all_added_coords))

coords_counties <- st_intersection(all_added_coords, target_counties)

# 4. checking accuracy
unique(coords_counties$Accuracy)

# source for checking accuracy levels: 
# https://map.dfg.ca.gov/rarefind/view/RF_FieldDescriptions.htm

# dropping those with inaccuracy of > 1000 meters
all_added_coords <- coords_counties |>
  dplyr::filter(!Accuracy %in% c("Circular feature with a 1600 meter radius (1 mile)",
                              "Circular feature with a 1300 meter radius (4/5 mile)"))

# seeing which observations we dropped - 12
dropped <- coords_counties |>
  dplyr::filter(Accuracy %in% c("Circular feature with a 1600 meter radius (1 mile)",
                             "Circular feature with a 1300 meter radius (4/5 mile)"))

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
# FILTERING ON 1 KM GRID
# adopted from 01_GBIF_data_pulling script, corresponds to data cleaning step (5)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
# grabbing our Chelsea grid 
if(!file.exists(here('data/CHELSA_bio1_1981-2010_V.2.1.tif'))){
  tmplt <- drive_get("CHELSA_bio1_1981-2010_V.2.1.tif")
  drive_download(tmplt, path=here('data/CHELSA_bio1_1981-2010_V.2.1.tif'))
}

# creating a template
raster_template <- terra::rast(here("data/CHELSA_bio1_1981-2010_V.2.1.tif"))

# getting resolution and extent of the raster data
res_xy <- res(raster_template)
raster_bbox <- st_as_sfc(st_bbox(raster_template))

# creating a full fishnet within the target extent
full_fishnet <- st_make_grid(raster_bbox, cellsize = res_xy, what = "polygons", square = TRUE)
full_fishnet_sf <- st_sf(full_fishnet)

# cropping and masking raster to the extent of the counties
fishnet_clipped <- st_intersection(full_fishnet_sf, target_counties)

# Add a unique ID column to the clipped fishnet polygons if it doesn't exist
if (!"grid_id" %in% colnames(fishnet_clipped)) {
  fishnet_clipped$grid_id <- 1:nrow(fishnet_clipped)
}

# joining rare occurrence data with 1km fishnet 
occ_grid <- st_join(all_added_coords, fishnet_clipped)

# removing those that do not fall into a grid cell 
occ_with_grid <- occ_grid |> 
  dplyr::filter(!is.na(grid_id))

# for each species and grid cell, keep only one record
occ_unique <- occ_with_grid |>
  group_by(SciName, grid_id) |>
  slice(1) |> # keeping the first record in each group 
  ungroup()

# viewing the species that were dropped (0)
occ_dropped <- occ_grid |>
  dplyr::filter(is.na(grid_id))

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
# MATCHING DATA FRAME TO OTHER PLANT/ANIMAL DATA
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
# selecting columns we want and renaming to match other Animal-Plant merged df
all_rare <- occ_unique|>
  as.data.frame() |>
  dplyr::select(SciName, TaxonGroup, X, Y) |>
  dplyr::rename(species = SciName,
                taxon = TaxonGroup,
                x = X,
                y = Y) 

# renaming "Dicots" to "Plants" to match the plant/animal data
all_rare_clean <- all_rare |>
  mutate(taxon = case_when(taxon == "Dicots" ~ "plant", 
                           taxon == "Birds" ~ "birds",
                           taxon == "Mammals" ~ "mammals",
                            TRUE ~ taxon)) 

# now viewing counts by species 
rare_counts <- all_rare_clean |>
  group_by(species, taxon) |>
  count()

# saving df 
st_write(all_rare_clean, here("stored location"), append = FALSE)
