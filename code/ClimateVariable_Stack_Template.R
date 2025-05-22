# Load required libraries
library(terra)
library(sf)
library(dplyr)

# 
# Define paths
# 
data_dir <- ""      
output_dir <- ""       
shapefile_path <- ""  

# 
# Filter desired raster files, for example: bio1, bio3, bio5, etc.
# 
bio_vars <- c("bio1", "bio3", "bio5", "bio6", "bio15", "bio17", "bio18")
pattern <- paste0("^CHELSA_[a-z0-9]+_2011-2040_gfdl-esm4_ssp585_V\\.2\\.1\\.tif$")

bio_files <- list.files(
  path = data_dir,
  pattern = pattern,
  full.names = TRUE
)

# Stack the selected rasters
bio_stack <- rast(bio_files)

# 
# Load shapefile and buffer 3 counties
# 
counties_sf <- st_read(shapefile_path)

ca_counties <- counties_sf %>% 
  filter(STATEFP == "06")  # California

threecounty_buffer <- ca_counties %>%
  filter(NAME %in% c("Santa Barbara", "Ventura", "San Luis Obispo")) %>%
  st_transform(3310) %>%      # California Albers (meters) to better process buffer
  st_union() %>%
  st_buffer(dist = 5000) %>%  # 5 km buffer
  st_transform(4329)          # Final CRS: EPSG:2229

# 
# Crop and reproject rasters
# 

bio_stack_cropped <- crop(bio_stack, vect(threecounty_buffer))

target_crs <- "EPSG:2229"
bio_stack_projected <- project(bio_stack_cropped, target_crs, method = "bilinear")

# 
# Multi-layer raster
# 

output_file <- file.path(output_dir, "Projected__selectedbios_2040_2070_GFDL_SSP585_stack.tif")

writeRaster(bio_stack_projected, filename = output_file, overwrite = TRUE)


#plot(bio_stack_projected)


