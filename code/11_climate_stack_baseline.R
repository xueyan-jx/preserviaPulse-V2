## Purpose: Preprocess CHELSA baseline climate variables (1981–2010),
##          run correlation analysis, and create a selected variable stack.
## Author: Zijun Li, Xue Yan, Yifei Liu
## Date: "2025-05-12"
## Output: html_document

  
library(here)
library(terra)
library(sf)
library(dplyr)
library(corrplot)
library(recipes)

# ---------------- 0. Paths & inputs ---------------- #

data_dir  <- here("data")
clim_dir  <- file.path(data_dir, "ClimateVariables")
stack_dir <- file.path(data_dir, "Stack_Env")

if (!dir.exists(stack_dir)) dir.create(stack_dir, recursive = TRUE)

# Refined tri-county boundary (for cropping)
refined_boundary <- st_read(file.path(data_dir, "threecounties_refine.geojson"))
refined_vect     <- vect(refined_boundary)

# Target CRS for modeling
target_crs <- "EPSG:2229"  # NAD83 / California zone 5 (ftUS)

# ---------------- 1. Load CHELSA baseline rasters ---------------- #

chelsa_files <- list.files(
  path       = file.path(clim_dir, "Chelsa_1981_2010_Climatologies_bio19"),
  pattern    = "^CHELSA_bio[0-9]{1,2}_1981-2010_V\\.2\\.1\\.tif$",
  full.names = TRUE
)

if (length(chelsa_files) == 0) {
  stop("No CHELSA baseline rasters found in Chelsa_1981_2010_Climatologies_bio19/")
}

# Stack all baseline CHELSA bioclim variables
chelsa_stack <- rast(chelsa_files)

# Extract layer names like "bio1", "bio2", ..., "bio19"
chelsa_names <- gsub("^.*(bio[0-9]+)_.*$", "\\1", basename(chelsa_files))
names(chelsa_stack) <- chelsa_names

# ---------------- 2. Crop & project to study region ---------------- #

# Make sure boundary and CHELSA stack are in the same CRS for cropping
if (crs(refined_vect) != crs(chelsa_stack)) {
  refined_for_crop <- project(refined_vect, crs(chelsa_stack))
} else {
  refined_for_crop <- refined_vect
}

# 2a. Crop & mask CHELSA rasters to tri-county boundary (in CHELSA CRS)
chelsa_cropped <- crop(chelsa_stack, refined_for_crop)
chelsa_cropped <- mask(chelsa_cropped, refined_for_crop)

# 2b. Project cropped stack to modeling CRS (EPSG:2229)
chelsa_projected <- project(chelsa_cropped, target_crs, method = "bilinear")

# Save full projected stack (baseline climate)
full_hist_file <- file.path(clim_dir, "chelsa_hist_1981_2010_projected.tif")
writeRaster(chelsa_projected, full_hist_file, overwrite = TRUE)

# ---------------- 3. Correlation analysis & variable selection ---------------- #

# Extract values to a data frame (may be large; this is a full-grid correlation)
clim_vals <- values(chelsa_projected, na.rm = TRUE)
clim_vals <- na.omit(clim_vals)
clim_df   <- as.data.frame(clim_vals)

# Full correlation matrix
cor_matrix <- cor(clim_df)

# Plot full correlation matrix
pdf(file.path(clim_dir, "chelsa_hist_clim_corrplot.pdf"),
    width = 10, height = 10)
corrplot(cor_matrix, method = "color", type = "upper", tl.cex = 0.8,
         addCoef.col = "black", number.cex = 0.5, diag = FALSE)
dev.off()

# Use recipes to remove highly correlated variables (threshold can be tuned)
rec <- recipe(~ ., data = clim_df) %>%
  step_corr(all_predictors(), threshold = 0.9)

prep_rec     <- prep(rec)
filtered_df  <- juice(prep_rec)
selected_vars <- colnames(filtered_df)

message("Selected CHELSA variables (after correlation filter):")
print(selected_vars)

# Subset raster stack to selected variables
chelsa_selected <- chelsa_projected[[selected_vars]]

# Recompute correlation matrix for selected variables
sel_vals <- values(chelsa_selected, na.rm = TRUE)
sel_vals <- na.omit(sel_vals)
cor_matrix_sel <- cor(sel_vals)

pdf(file.path(clim_dir, "chelsa_hist_selected_clim_corrplot.pdf"),
    width = 6, height = 6)
corrplot(cor_matrix_sel, method = "color", type = "upper", tl.cex = 0.8,
         addCoef.col = "black", number.cex = 0.5, diag = FALSE)
dev.off()

# ---------------- 4. Save selected baseline climate stack ---------------- #

# This is the stack used by the environmental script and SDMs
sel_hist_file <- file.path(clim_dir, "projected_selectedbios_1980_2010_stack.tif")
writeRaster(chelsa_selected, sel_hist_file, overwrite = TRUE)

message("✓ Wrote baseline CHELSA stacks:")
message("  full:     ", full_hist_file)
message("  selected: ", sel_hist_file)

# Note:
# Future climate scenarios (e.g., CHELSA CMIP6 SSP126/370/585) can be
# processed with the same pattern and saved with analogous file names, e.g.:
#   projected_selectedbios_2040_2070_GFDL_SSP126_stack.tif
#   projected_selectedbios_2040_2070_GFDL_SSP370_stack.tif
#   projected_selectedbios_2040_2070_GFDL_SSP585_stack.tif