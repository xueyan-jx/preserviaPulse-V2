## Purpose of script: Generate Slope, Aspect, TRI, Flow Accumulation from Elevation
## Authors: GEOG 274
## Date: Spring, 2025
## Credits to: Yifei Liu

# ------------- Setting up --------------
# Install and load packages
library(here)
library(terra)
library(sf)
library(dplyr)

# Set a directory for data
here() # first check path

# Define folder for boundary and elevation leyers
layer_dir <- here("data", "Elevation")
# if it's not there already, create it
if (!dir.exists(layer_dir)) dir.create(layer_dir, recursive = TRUE) 

#------------- 1. Boundary Data Preparation  -------------
# a) Load California county boundaries (download separately from: 
# https://data.ca.gov/dataset/ca-geographic-boundaries -> "CA County Boundaries")

# Unzip shapefile into a folder (e.g., "data/Elevation/ca_counties/")
county_shp_path <- file.path(layer_dir, "ca_counties", "CA_Counties.shp")

# Load counties using terra
ca_counties <- vect(county_shp_path)

# Filter only the three target counties
county_names <- c("Santa Barbara", "Ventura", "San Luis Obispo")
three_counties <- ca_counties[ca_counties$NAME %in% county_names, ]

# Plot check
plot(three_counties, col = "lightgreen", main = "Three Counties")
summary(three_counties)

# b) Merge three counties into one feature (dissolve boundary)
three_counties_union <- terra::aggregate(three_counties)

# Plot check
plot(three_counties_union, col = "lightblue", main = "Three Counties Union")

# c) Add 100m buffer around boundary to ensure no edge pixels are lost
three_counties_buffered <- buffer(three_counties_union, width = 100)

# Plot check
plot(three_counties_buffered, col = "pink", main = "Three Counties Buffered")

# d) Reproject to NAD83 / California State Plane Zone 5 (EPSG:2229)
three_counties_projected <- project(three_counties_buffered, "EPSG:2229")

# Plot check
plot(three_counties_projected, col = "orange", main = "Three Counties Projected")

# Save the processed boundary
writeVector(three_counties_projected, 
            file.path(layer_dir, "three_counties_projected.shp"), 
            overwrite = TRUE)


# ------------- 2. DEM Acquisition and Preprocessing --------------

# a) Load SRTM 30-meter DEM raster
# Download SRTM from USGS EarthExplorer or CGIAR-CSI SRTM sites
# Transform projected vector to WGS84 for EarthExplorer
three_counties_latlon <- project(three_counties_projected, "EPSG:4326")

# Get bounding box in lat/lon
ext <- ext(three_counties_latlon)
ext

# Define path to SRTM tiles
srtm_tile_dir <- here("data", "Elevation", "SRTM_raw")
srtm_tiles <- list.files(srtm_tile_dir, pattern = "\\.tif$", full.names = TRUE)

# Load and mosaic all tiles
srtm_list <- lapply(srtm_tiles, rast)
srtm_mosaic <- do.call(mosaic, srtm_list)

# Reproject to NAD83 / California State Plane Zone 5 (EPSG:2229)
srtm_projected <- project(srtm_mosaic, "EPSG:2229")

# Crop and mask using buffered tri-county boundary
srtm_cropped <- crop(srtm_projected, three_counties_projected)
srtm_clipped <- mask(srtm_cropped, three_counties_projected)

# Save final processed DEM
writeRaster(srtm_clipped,
            filename = file.path(layer_dir, "srtm_30m_clipped_projected.tif"),
            overwrite = TRUE)

# Plot check
plot(srtm_clipped, main = "Clipped & Projected DEM (Tri-County)")

# -------- 3a. Compute Slope --------

# Compute slope in degrees
slope_deg <- terrain(srtm_clipped, v = "slope", unit = "degrees")

# Plot check
plot(slope_deg, main = "Slope (degrees)", col = terrain.colors(100))

# Save output
writeRaster(slope_deg,
            filename = file.path(layer_dir, "slope_three_county_deg.tif"),
            overwrite = TRUE)

# -------- 3b. Compute Aspect --------

# Calculate aspect in degrees
aspect_deg <- terrain(srtm_clipped, v = "aspect", unit = "degrees")

# Plot check (note: aspect is circular — 0° and 360° = North)
plot(aspect_deg, main = "Aspect (degrees)", col = hcl.colors(100, "Zissou 1"))

# Save to file
writeRaster(aspect_deg,
            filename = file.path(layer_dir, "aspect_three_county_deg.tif"),
            overwrite = TRUE)

# -------- 3c. Compute Terrain Ruggedness Index (TRI) --------

# Calculate TRI from the clipped DEM
tri <- terrain(srtm_clipped, v = "TRI", unit='degrees')

# Plot check
rgb.palette <- colorRampPalette(c('lightblue',"blue","green","yellow","red", "red3"),
                                space = "rgb")
plot(tri, main = "Terrain Ruggedness Index (TRI)", col = rgb.palette(100))

# Save output raster
writeRaster(tri,
            filename = file.path(layer_dir, "TRI_three_county_deg.tif"),
            overwrite = TRUE)
