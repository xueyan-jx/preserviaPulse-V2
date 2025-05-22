---
  #title: "Climate Variables Conservation GIS 2025"
  #author: "Zijun Li"
  #date: "2025-05-12"
  #output: html_document
---
  
# Load necessary libraries, may only need terra for this task. 
library(terra)
library(dplyr)
library(sf)

# Define path, extent, and environmental variable names
data_dir <- "/LOCAL_PATH"
env_names <- c(paste0("bio", 1:19))

# Load shapefile (adjust path to your file)
counties_sf <- st_read("/PATH_TO_BOUNDARY")

# Get California counties
ca_counties <- counties_sf %>% 
  filter(STATEFP == "06")

threecounty_buffer <- ca_counties %>%
  filter(NAME %in% c("Santa Barbara", "Ventura", "San Luis Obispo")) %>%
  st_transform(3310) %>% # California Albers (units: meters)
  st_union() %>%
  st_buffer(dist = 5000) %>% # buffer 5 km
  st_transform(4326)

# Load and stack raster files
env_files <- list.files(
  path = data_dir, pattern = "^CHELSA_.*\\.tif$", full.names = TRUE)
env_stack <- rast(env_files)
env_stack_cropped <- crop(env_stack, threecounty_buffer)

# Define the target CRS (EPSG:2229 - NAD83 / California zone 5 (ftUS))
target_crs <- "EPSG:2229"

# Reproject the cropped raster stack to EPSG:2229
env_stack_projected <- project(env_stack_cropped, target_crs, method = "bilinear")

# Option 1: Export as a Single Multi-Layer Raster (stacked GeoTIFF)
# Define output path and filename
output_dir <-"/LOCAL_PATH_STORE"
output_file <- file.path(output_dir,"projected_env_stack.tif")

# Write the stack to a single multi-layer GeoTIFF
writeRaster(env_stack_projected, filename = output_file, overwrite = TRUE)

plot(env_stack_projected[[1:16]]) 
# Base R plotting windows (especially on Windows or macOS) sometimes default to 
# showing only 16 layers per multi-panel plot.

#plot individually in a loop
# plot(env_stack_projected[[1:19]])
# for (i in 1:nlyr(env_stack_projected)) {
#  +     plot(env_stack_projected[[i]], main = names(env_stack_projected)[i])
#  +     Sys.sleep(1)  # Pause 1 second between plots
#  + }

#Option 2: Export Each Layer as a Separate File
# Create an output subdirectory
output_dir <- file.path(output_dir, "projected_layers")
dir.create(output_dir, showWarnings = FALSE)

# Export each layer individually
for (i in 1:nlyr(env_stack_projected)) {
  layer_name <- names(env_stack_projected)[i]
  output_path <- file.path(output_dir, paste0(layer_name, ".tif"))
  writeRaster(env_stack_projected[[i]], filename = output_path, overwrite = TRUE)
}

#check resolution: res(raster)
#check projection: crs(raster)
#check raster extent: ext(raster)
#nlyr(raster_stack)
#were all layers loaded? names(raster_stack)
#are some layers blank or fully NA? is.na(raster_stack[[i]])


#---------------------------------------------------------------------------------------
#Do pearson correlation of 19 bio climate variables
#---------------------------------------------------------------------------------------

# Extract and clean values
env_values <- values(env_stack_projected, na.rm = TRUE)
env_values <- na.omit(env_values)
env_df <- as.data.frame(env_values)

# Correlation plot (full)
cor_matrix <- cor(env_df)
pdf(file.path(data_dir, "chelsa_hist_clim_corrplot.pdf"), width = 10, height = 10)
corrplot(cor_matrix, method = "color", type = "upper", tl.cex = 0.8,
         addCoef.col = "black", number.cex = 0.5, diag = FALSE)
dev.off()

# Use tidymodels' recipe to remove highly correlated variables, change the threshold here to compare.
rec <- recipe(~ ., data = env_df) %>%
  step_corr(all_predictors(), threshold = 0.9) 

prep_rec <- prep(rec)
filtered_df <- juice(prep_rec)

# Variables kept after correlation filter
selected_vars <- colnames(filtered_df)

# Subset raster stack
env_stack_selected <- env_stack_projected[[selected_vars]]

# Recompute and plot correlation matrix after filtering
env_selected_values <- values(env_stack_selected, na.rm = TRUE)
env_selected_values <- na.omit(env_selected_values)
cor_matrix <- cor(env_selected_values)

pdf(file.path(data_dir, "chelsa_hist_selected_clim_corrplot.pdf"), width = 6, height = 6)
corrplot(cor_matrix, method = "color", type = "upper", tl.cex = 0.8,
         addCoef.col = "black", number.cex = 0.5, diag = FALSE)
dev.off()
