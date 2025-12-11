## Purpose: Generate terrain, solar, and distance-to-coast layers,
##          and stack them with selected climate variables for
##          current (1980–2010) and optional future scenarios.
## Authors: GEOG 274
## Date: Spring, 2025
## Credits to: Yifei Liu

library(here)
library(terra)
library(sf)
library(dplyr)
library(tigris)

options(tigris_use_cache = TRUE)

# ---------------- 0. Paths & helpers ---------------- #

data_dir   <- here("data")
elev_dir   <- file.path(data_dir, "Elevation")
clim_dir   <- file.path(data_dir, "ClimateVariables")
stack_dir  <- file.path(data_dir, "Stack_Env")

if (!dir.exists(elev_dir))  dir.create(elev_dir,  recursive = TRUE)
if (!dir.exists(stack_dir)) dir.create(stack_dir, recursive = TRUE)

# Template raster grid (final modeling grid)
template_raster <- rast(file.path(data_dir, "raster_template.tif"))

# Refined tri-county boundary for clipping
refined_sf   <- st_read(file.path(data_dir, "threecounties_refine.geojson"))
refined_vect <- vect(refined_sf)
refined_proj <- project(refined_vect, crs(template_raster))

# Helper: align a raster to template and clip to boundary
align_and_clip <- function(r, template, boundary_vect) {
  if (crs(r) != crs(template)) {
    r <- project(r, crs(template), method = "near")
  }
  r <- crop(r, ext(template))
  r <- extend(r, ext(template))
  r <- resample(r, template, method = "near")
  r <- crop(r, boundary_vect)
  r <- mask(r, boundary_vect)
  r
}

#------------- 1. Boundary Data Preparation  -------------
# a) Load California county boundaries (download separately from: 
# https://data.ca.gov/dataset/ca-geographic-boundaries -> "CA County Boundaries")

# Unzip shapefile into a folder (e.g., "data/Elevation/ca_counties/")
county_shp_path <- file.path(elev_dir, "ca_counties", "CA_Counties.shp")

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
            file.path(elev_dir, "three_counties_projected.shp"), 
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
            filename = file.path(elev_dir, "srtm_30m_clipped_projected.tif"),
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
            filename = file.path(elev_dir, "slope_three_county_deg.tif"),
            overwrite = TRUE)

# -------- 3b. Compute Aspect --------

# Calculate aspect in degrees
aspect_deg <- terrain(srtm_clipped, v = "aspect", unit = "degrees")

# Plot check (note: aspect is circular — 0° and 360° = North)
plot(aspect_deg, main = "Aspect (degrees)", col = hcl.colors(100, "Zissou 1"))

# Save to file
writeRaster(aspect_deg,
            filename = file.path(elev_dir, "aspect_three_county_deg.tif"),
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
            filename = file.path(elev_dir, "TRI_three_county_deg.tif"),
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
            filename = file.path(elev_dir, "flow_accum_three_county.tif"),
            overwrite = TRUE)

# ------------- 5. Distance to Coast -----------------------------
# Query the coastlines from tigris package, transform, rasterize, and 
# calculate the distance.
## Note: unit parameter does not work in function distance. The default is m.
dist_to_cstl <- coastline() %>% 
  st_transform(crs(template_raster)) %>% 
  rasterize(template_raster) %>% 
  terra::distance(x = ., unit = "km", method = "geo") %>% 
  mask(template_raster)

# Deal with the bug in distance function
dist_to_cstl <- dist_to_cstl / 1000  # meters → km
names(dist_to_cstl) <- "dist_coast"

# Align and clip
dist_coast_f <- align_and_clip(dist_to_cstl, template_raster, refined_proj)

# Save output
writeRaster(
  dist_coast_f,
  file.path(elev_dir, "dist_coast_three_county_km.tif"),
  overwrite = TRUE
)

# ---------------- 6. Solar radiation ---------------- #

solar_raw <- rast(file.path(elev_dir, "AreaSol_srtm_30m_Clip3Counties.tif"))
solar_raw[solar_raw < 691015] <- NA  # optional cleanup threshold

solar_f <- align_and_clip(solar_raw, template_raster, refined_proj)
names(solar_f) <- "solar"

# ---------------- 7. Align terrain layers to template ---------------- #

slope_raw    <- rast(file.path(elev_dir, "slope_three_county_deg.tif"))
aspect_raw   <- rast(file.path(elev_dir, "aspect_three_county_deg.tif"))
flow_acc_raw <- rast(file.path(elev_dir, "flow_accum_three_county.tif"))

slope_f    <- align_and_clip(slope_raw,    template_raster, refined_proj)
aspect_f   <- align_and_clip(aspect_raw,   template_raster, refined_proj)
flow_acc_f <- align_and_clip(flow_acc_raw, template_raster, refined_proj)

names(slope_f)    <- "slope"
names(aspect_f)   <- "aspect"
names(flow_acc_f) <- "flow_acc"

# ------------- 8. Load selected climate stacks from script 11------------------
# -------- 8a. Load and preprocess layers -----------------
climate_paths <- c(
  past   = "projected_selectedbios_1980_2010_stack.tif",
  ssp126 = "projected_selectedbios_2040_2070_GFDL_SSP126_stack.tif",
  ssp370 = "projected_selectedbios_2040_2070_GFDL_SSP370_stack.tif",
  ssp585 = "projected_selectedbios_2040_2070_GFDL_SSP585_stack.tif"
)

climate_rasters <- list()
for (nm in names(climate_paths)) {
  f <- file.path(clim_dir, climate_paths[[nm]])
  if (file.exists(f)) {
    message("→ Loading climate stack: ", f)
    r <- rast(f)
    climate_rasters[[nm]] <- align_and_clip(r, template_raster, refined_proj)
  } else {
    message("⚠ Skipping ", nm, " climate stack (file not found): ", f)
  }
}

# -------- 8b. Stack env layers + climate for each scenario -----------------
stack_and_save <- function(climate_layer, scenario_name) {
  env_stack <- c(
    slope_f,
    aspect_f,
    flow_acc_f,
    dist_coast_f,
    solar_f,
    climate_layer
  )
  names(env_stack) <- c(
    "slope", "aspect", "flow_acc",
    "dist_coast", "solar",
    names(climate_layer)
  )
  out_file <- file.path(stack_dir, paste0("final_env_", scenario_name, ".tif"))
  writeRaster(env_stack, out_file, overwrite = TRUE)
  message("✓ Wrote environmental stack: ", out_file)
}

# Current stack (required)
if (!is.null(climate_rasters$past)) {
  stack_and_save(climate_rasters$past, "1980_2010_stack")
} else {
  stop("No baseline climate stack found (projected_selectedbios_1980_2010_stack.tif).")
}

# Optional future stacks (only if files exist)
if (!is.null(climate_rasters$ssp126)) {
  stack_and_save(climate_rasters$ssp126, "ssp126_2040_2070_stack")
}
if (!is.null(climate_rasters$ssp370)) {
  stack_and_save(climate_rasters$ssp370, "ssp370_2040_2070_stack")
}
if (!is.null(climate_rasters$ssp585)) {
  stack_and_save(climate_rasters$ssp585, "ssp585_2040_2070_stack")
}

message("✓ Finished building environmental stacks.")

