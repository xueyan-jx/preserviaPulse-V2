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

# Load raster template
template_raster <- rast("data/raster_template.tif")

# Load refined county boundry
refined_county_sf <- st_read("data/threecounties_refine.geojson")
refined_county_vect <- vect(refined_county_sf)
# Convert and reproject boundary to EPSG:2229
refined_county_vect_proj <- project(vect(refined_county_sf), crs(template_raster))


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

# ------------- 3.	Terrain Derivative Calculation --------------
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

# ------------- 4. Hydrological Analysis -----------------------------
# Compute Flow Direction
flow_dir <- terrain(srtm_clipped, v = "flowdir", unit='degrees')

# Plot check
plot(flow_dir, main = "Flow Direction", col = hcl.colors(100, "YlGnBu"))

# Compute Flow Accumulation
flow_accum <- flowAccumulation(flow_dir)

# Visualize with log transformation for better contrast
log_accum <- log1p(flow_accum)  # log(1 + x) to avoid log(0)
plot(log_accum, main = "Flow Accumulation (log scale)", col = hcl.colors(100, "YlGnBu"))

# Save output
writeRaster(flow_accum,
            filename = file.path(layer_dir, "flow_accum_three_county.tif"),
            overwrite = TRUE)

# ------------- 5. Distance to Coast -----------------------------
# Load your combined plant/animal data
anim_plant_df <- read.csv("data/Anim_Plant_merge.csv")

# Convert to sf object (assuming x and y are in State Plane Zone 5, EPSG:2229)
anim_plant_sf <- st_as_sf(anim_plant_df, coords = c("x", "y"), crs = 2229)

# Load coastline shapefile (adjust path)
coastline <- st_read("data/Elevation/ECU_clipped/county_ECU.shp")

# Reproject to match animal/plant points
coastline_2229 <- st_transform(coastline, crs = 2229)

# Compute distance from each observation to the coastline
dist_matrix <- st_distance(anim_plant_sf, coastline_2229)  # units: US survey feet

# Extract minimum distance for each row (across all coastal polygons/lines)
min_distances <- apply(dist_matrix, 1, min)

# Add as a new column (as numeric)
anim_plant_sf$distance_to_coast <- as.numeric(min_distances)

# Extract coordinates from geometry
coords <- st_coordinates(anim_plant_sf)

# Add x and y back to the dataframe
anim_plant_df_with_coords <- anim_plant_sf %>%
  st_drop_geometry() %>%
  mutate(x = coords[, 1],
         y = coords[, 2])

# Write to CSV with distance_to_coast, x, and y included
write.csv(anim_plant_df_with_coords, "data/Anim_Plant_merge_with_dist_to_coast.csv", row.names = FALSE)

# Convert to terra vector
dist_vect <- vect(anim_plant_sf)

# Rasterize using the template
dist_raster <- rasterize(dist_vect, template_raster, field = "distance_to_coast")

# Visualize
plot(dist_raster, main = "Distance to Coast (Rasterized)")

# Save output
writeRaster(dist_raster, "data/Elevation/dist_to_coast.tif", overwrite = TRUE)


# ------------- 6. Stack Environmental layers -----------------------------
# -------- 6a. Load and preprocess layers -----------------
# Load individual layers
slope <- rast(file.path(layer_dir, "slope_three_county_deg.tif"))
aspect <- rast(file.path(layer_dir, "aspect_three_county_deg.tif"))
flow_acc <- rast(file.path(layer_dir, "flow_accum_three_county.tif"))
dist_coast <- rast(file.path(layer_dir, "dist_to_coast.tif"))
solar <- rast(file.path(layer_dir, "AreaSol_srtm_30m_Clip3Counties.tif"))
# Set values < 691015 to NA (optional cleanup)
solar[solar < 691015] <- NA

climate_paths <- list(
  past = "projected_selectedbios_1980_2010_stack.tif",
  ssp126 = "projected_selectedbios_2040_2070_GFDL_SSP126_stack.tif",
  ssp370 = "projected_selectedbios_2040_2070_GFDL_SSP370_stack.tif",
  ssp585 = "Projected__selectedbios_2040_2070_GFDL_SSP585_stack.tif"
)

climate_rasters <- lapply(climate_paths, function(path) {
  rast(file.path(here("data", "ClimateVariables"), path))
})

# Rename all climate layer bands consistently
climate_names <- c("bio1", "bio15", "bio17", "bio18", "bio3", "bio5", "bio6")
names(climate_rasters$past) <- climate_names
names(climate_rasters$ssp126) <- climate_names
names(climate_rasters$ssp370) <- climate_names
names(climate_rasters$ssp585) <- climate_names

# -------- 6b. Resample & Mask All Layers -----------------
process_layer <- function(r) {
  r_res <- resample(mask(crop(r, refined_county_vect_proj), 
                         refined_county_vect_proj
                         ), 
                    template_raster
                    )
  mask(r_res, !is.na(template_raster))
}

# Apply processing to each environmental layer
slope_f <- process_layer(slope)
aspect_f <- process_layer(aspect)
flow_acc_f <- process_layer(flow_acc)
dist_coast_f <- process_layer(dist_raster)
solar_f <- process_layer(solar)
climate_f <- lapply(climate_rasters, process_layer)

# -------- 6d. Stack all environmental layers together -----------------
stack_and_save <- function(climate_layer, scenario_name) {
  env_stack <- c(slope_f, aspect_f, flow_acc_f, dist_coast_f, 
                 solar_f, climate_layer)
  names(env_stack) <- c("slope", "aspect", "flow_acc", 
                        "dist_coast", "solar", names(climate_layer))
  writeRaster(env_stack, paste0("data/Stack_Env/final_env_", 
                                scenario_name, ".tif"), 
              overwrite = TRUE)
}

# Save All Stacks
stack_and_save(climate_f$past, "1980_2010_stack")
stack_and_save(climate_f$ssp126, "ssp126_2040_2070_stack")
stack_and_save(climate_f$ssp370, "ssp370_2040_2070_stack")
stack_and_save(climate_f$ssp585, "ssp585_2040_2070_stack")

# Plot check
plot(rast(file.path("data/Stack_Env/final_env_ssp126_2040_2070_stack.tif")))

