## Purpose of script: Uncertainty map
## Authors: GEOG 274
## Date: Spring, 2025
## Credits to: Yifei Liu, Wenxin Yang, Yanni Zhan, Xue Yan

rm(list=ls())

library(terra)
# library(raster)
library(here)
library(RColorBrewer)

# Denfine species name and file paths
species_name <- "Arctostaphylos purissima"  #"Taxidea taxus"  #"Haliaeetus leucocephalus" 
eva_dir       <- here("results", "plants_results", "evaluations")
jldp_path     <- here("visualization", "jldp_boundary.shp")
current_unc_path <- here("results", "current", paste0(species_name, "_current_uncertainty.tif"))
future_unc_path <- here(eva_dir, paste0(species_name, "_ssp370_uncertainty.tif"))
var_scenarios_path <- here(eva_dir, paste0(species_name, "_var_scenarios.tif"))


#    Read them all into a named list of SpatRaster
r_current_unc <- rast(current_unc_path)
r_future_unc  <- rast(future_unc_path)
r_var_scen    <- rast(var_scenarios_path)
rasters      <- list(r_current_unc, r_future_unc, r_var_scen)

# Read the JLD preserve boundary
jldp_vec <- vect(jldp_path)

#  Define concise titles for each subfigure
short_titles <- c(
  "Current Uncertainty",
  "Future Uncertainty (ssp370)",
  "Variance across scenarios"
)

# Compute a shared min/max for color scaling
all_min <- min(
  global(r_current_unc, "min", na.rm = TRUE),
  global(r_future_unc,  "min", na.rm = TRUE),
  global(r_var_scen,    "min", na.rm = TRUE)
)
all_max <- max(
  global(r_current_unc, "max", na.rm = TRUE),
  global(r_future_unc,  "max", na.rm = TRUE),
  global(r_var_scen,    "max", na.rm = TRUE)
)

# Build a red-yellow-blue color palette
palette_n   <- 100
cw_palette  <- colorRampPalette(rev(brewer.pal(11, "RdBu")))(palette_n)


# Define three fig regions for 1*3 layout
fig_regions <- list(
  c(0.00, 0.333, 0.00, 1.00),  # left third
  c(0.333, 0.667, 0.00, 1.00), # middle third
  c(0.667, 1.000, 0.00, 1.00)  # right third
)

# Open a PNG device (adjust filename, size, and resolution as desired)
png(
  filename = here("visualization", paste0(species_name, "_uncertainty_compare.png")),
  width    = 12,    # in inches
  height   = 5,     # in inches
  units    = "in",
  res      = 300    # dpi
)

#  Loop through each raster + its fig region, and plot main + inset
for (i in seq_along(rasters)) {
  r        <- rasters[[i]]
  title_i  <- short_titles[i]
  fr_main  <- fig_regions[[i]]  # main map region
  
  #  Plot the main map in its cell
  par(fig = fr_main, new = (i > 1), 
      mar = c(0.5, 0.5, 1, 0.5), 
      mgp = c(0, 0.5, 0), tck = -0.01)
  plot(
    r,
    col   = cw_palette,
    zlim  = c(all_min, all_max),
    main = title_i,
    axes = TRUE,
    box  = FALSE, 
    legend = TRUE
  )

  #  Overlay JLDP boundary from shapefile in red on the main map
  lines(jldp_vec, col = "red", lwd = 1)

  #  Compute inset coordinates slightly left of the top‐right of the main cell
  x1 <- fr_main[1]
  x2 <- fr_main[2]
  y1 <- fr_main[3]
  y2 <- fr_main[4]

  fr_inset <- c(
    x1 + 0.4 * (x2 - x1),
    x1 + 0.85 * (x2 - x1),
    y1 + 0.55 * (y2 - y1),
    y1 + 0.88* (y2 - y1)
  )
  
  #  Plot inset rectangle (cropped extent, unmasked) with its own red boundary
  par(fig = fr_inset, new = TRUE, mar = c(0, 0, 0, 0))
  r_crop <- crop(r, jldp_vec)    # only crop, do not mask
  plot(
    r_crop,
    col = cw_palette,
    zlim  = c(all_min, all_max),
    axes   = FALSE,
    box    = FALSE,
    legend = FALSE     # remove inset’s color bar
  )
  lines(jldp_vec, col = "red", lwd = 1.5)  # JLDP boundary in red
  # box(col = "black", lwd = 1)            # black rectangle around the inset
}


#  Reset plotting parameters
# par(fig = c(0, 1, 0, 1), new = FALSE)
dev.off()

