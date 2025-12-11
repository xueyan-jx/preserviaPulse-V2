## Purpose: Extend the current environmental stack with the full
##          set of CHELSA baseline bioclimatic variables (bio1–bio19).
## Authors: GEOG 274
## Date: Spring, 2025

# Install and load packages
library(here)
library(terra)
library(sf)

# Set a directory for data
here() # first check path

# ---------------- 0. Paths & inputs ---------------- #
data_dir   <- here("data")
clim_dir   <- file.path(data_dir, "ClimateVariables")
stack_dir  <- file.path(data_dir, "Stack_Env")

template_raster <- rast(file.path(data_dir, "raster_template.tif"))

refined_sf   <- st_read(file.path(data_dir, "threecounties_refine.geojson"))
refined_vect <- vect(refined_sf)
refined_proj <- project(refined_vect, crs(template_raster))

# Helper: align and clip to template + boundary
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

# ---------------- 1. Load existing current env stack ---------------- #
env_sel_file <- file.path(stack_dir, "final_env_1980_2010_stack.tif")

if (!file.exists(env_sel_file)) {
  stop("Current env stack not found: ", env_sel_file,
       "\nRun 12_Envi_layer_elev_stack.R first.")
}

env_sel   <- rast(env_sel_file)
# Get the names of bands already in that stack
# (Assuming each layer has a meaningful name like "bio1", "bio2", …)
sel_names <- names(env_sel)

message("Current env stack layers:")
print(sel_names)

# ---------------- 2. Load CHELSA baseline rasters ---------------- #

# List all CHELSA files and derive their layer names
chelsa_dir <- file.path(clim_dir, "Chelsa_1981_2010_Climatologies_bio19")

chelsa_files <- list.files(
  path       = chelsa_dir,
  pattern    = "^CHELSA_bio[0-9]{1,2}_1981-2010_V\\.2\\.1\\.tif$",
  full.names = TRUE
)

if (length(chelsa_files) == 0) {
  stop("No CHELSA baseline rasters found in ", chelsa_dir)
}

chelsa_names <- gsub("^.*(bio[0-9]+)_.*$", "\\1", basename(chelsa_files))

# Identify which bios are missing from env_sel
missing_idx   <- which(!chelsa_names %in% sel_names)
missing_files <- chelsa_files[missing_idx]
missing_names <- chelsa_names[missing_idx]

if (length(missing_files) == 0) {
  message("No missing CHELSA variables. full_env_1980_2010_stack = final_env_1980_2010_stack.")
  out_file <- file.path(stack_dir, "full_env_1980_2010_stack.tif")
  writeRaster(env_sel, out_file, overwrite = TRUE)
  quit(save = "no")
}

message("Missing CHELSA variables to add:")
print(missing_names)

# ---------------- 3. Align & clip missing CHELSA layers ---------------- #
aligned_clipped_list <- lapply(seq_along(missing_files), function(i) {
  r <- rast(missing_files[i])
  r <- align_and_clip(r, template_raster, refined_proj)
  names(r) <- missing_names[i]
  r
})

# Combine all aligned rasters
new_layers   <- do.call(c, aligned_clipped_list)
full_current <- c(env_sel, new_layers)

# Write out the combined stack
out_file <- file.path(stack_dir, "full_env_1980_2010_stack.tif")
writeRaster(full_current, out_file, overwrite = TRUE)

# Check
message("✓ Wrote extended baseline stack: ", out_file)
message("Final layer names:")
print(names(full_current))
