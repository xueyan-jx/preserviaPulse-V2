## Purpose of script: Generate background sample for all species using KDE to have bias correction
## Authors: GEOG 274
## Date: Spring, 2025
## Credits to: Yifei Liu, Zijun Li, Isabella Perez

# ------------- Setting up --------------
# Install and load packages
library(here)
library(terra)
library(sf)
library(dplyr)
library(spatialEco)

# Set a directory for data
here() # first check path

# Define folder for boundary and elevation leyers
layer_dir <- here("data")
# if it's not there already, create it
if (!dir.exists(layer_dir)) dir.create(layer_dir, recursive = TRUE)

# Load raster template
template_raster <- rast("data/raster_template.tif")

# Load refined county boundry
refined_county_sf <- st_read("data/threecounties_refine.geojson")
refined_county_vect <- vect(refined_county_sf)
# Convert and reproject boundary to EPSG:2229
refined_county_vect_proj <- project(vect(refined_county_sf), crs(template_raster))
# Convert SpatVector back to sf
refined_county_sf_proj <- st_as_sf(refined_county_vect_proj)


#------------- 1. Load Current Species Obs Data  -------------
# Load your combined plant/animal data
anim_plant_df <- read.csv("data/Anim_Plant_merge_final0613.csv")

# Convert to sf object (assuming x and y are in State Plane Zone 5, EPSG:2229)
anim_plant_sf <- st_as_sf(anim_plant_df, coords = c("x", "y"), crs = 2229)

# ------------- 2. Kernel Density Estimation -------------
# Unweighted KDE (spatial location only) with 1km resolution
pt.kde<- sf.kde(x = anim_plant_sf, ref = template_raster,
                standardize = TRUE, res = 1000)
plot(pt.kde, main="Unweighted kde")

# ------------- 3. Sample Background Points -------------
# Set seed for reproducibility
set.seed(274)

# Number of background points
n_bg <- 10000

# Sample background points based on KDE probabilities
bg_points <- spatSample(
  x = pt.kde,
  size = n_bg,
  method = "weights",       # KDE as weights
  as.points = TRUE,
  na.rm = TRUE
)

# Convert to sf for easy plotting and export
bg_sf <- st_as_sf(bg_points)

# Remove background points outside refined boundary (just in case)
bg_sf <- bg_sf[refined_county_sf_proj, ]


# Plot for inspection
plot(pt.kde, main = "KDE + Sampled Background Points")
plot(bg_sf$geometry, col = "red", pch = 16, cex = 0.3, add = TRUE)


# ------------- 4. Output -------------

# KDE raster layer save to tif
writeRaster(pt.kde, "data/background_kde_raster.tif",overwrite = TRUE)

# Extract coordinates from geometry
coords <- st_coordinates(bg_sf)

# Add x and y back to the dataframe
bg_sf_with_coords <- bg_sf %>%
  st_drop_geometry() %>%
  rename(kde_weight = z) %>%
  mutate(x = coords[, 1],
         y = coords[, 2],
         background = 1)

# Background point save to csv
write.csv(bg_sf_with_coords, "data/sampled_background_points.csv", row.names = FALSE)
