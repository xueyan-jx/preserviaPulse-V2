## Purpose of script: Visualize Plant Occurance Data
## Authors: GEOG 274
## Date: Spring, 2025
## Credits to: Yifei Liu, Xue Yan, Yanni Zhan

# ------------- Setting up --------------
# Install and load packages
library(here)
library(sf)
library(dplyr)
library(ggplot2)
library(ggspatial)
library(tigris)
library(viridis)

# Set a directory for data
here() # first check path

# Define folder for boundary and occurrence data
layer_dir <- here("data")
# if it's not there already, create it
if (!dir.exists(layer_dir)) dir.create(layer_dir, recursive = TRUE) 

#------------- 1. Load three counties boundary  -------------
# Gather CA county boundary file
ca_counties <- counties(state = "CA", cb = TRUE)

# Filter target three counties
target_counties <- ca_counties %>%
  filter(NAME %in% c("Santa Barbara", "Ventura", "San Luis Obispo")) 

# Reproject to NAD83 / California State Plane Zone 5 (EPSG:2229)
target_counties_projected <- st_transform(target_counties, "EPSG:2229")

# Convert to sf
target_counties_sf <- st_as_sf(target_counties_projected)

# ------------- 2. Load plant occurrence data --------------
# Top 6 species
top6_df <- read.csv("data/BIEN_GBIF0000607-250515123054153-final_top6_xy.csv")
top6_sf <- st_as_sf(top6_df, coords = c("x", "y"), crs = 2229)

# Final 17 species
final17_df <- read.csv("data/BIEN_GBIF0000607-250515123054153-final_17_xy.csv")
final17_sf <- st_as_sf(final17_df, coords = c("x", "y"), crs = 2229)

# ------------- 3. Plot Top 5/ Final 17 plant species --------------
# Confirm CRS is EPSG:2229
st_crs(target_counties_sf)  # Should return 2229
st_crs(top6_sf)         # Should also be 2229
st_crs(final17_sf)         # Should also be 2229

# Plot Top 5
ggplot() +
  geom_sf(data = target_counties_sf, fill = "gray95", color = "black") +
  geom_sf(data = top6_sf, aes(color = species), alpha = 0.7, size = 1) +
  labs(
    title = "Top 6 Most Frequent Plant Species",
    subtitle = "NAD83 / California State Plane Zone 5 (EPSG:2229)",
    color = "Species"
  ) +
  coord_sf(crs = 2229, default_crs = NA, datum = NA) +
  annotation_scale(
    location = "bl",        # bottom-left
    width_hint = 0.25,      # proportion of scale bar relative to width
    bar_cols = c("grey60", "white")
  ) +
  annotation_north_arrow(
    location = "tr",        # top-right
    which_north = "true",   # geographic north
    pad_x = unit(0.5, "cm"),
    pad_y = unit(0.5, "cm"),
    style = north_arrow_fancy_orienteering()
  ) +
  theme_bw()
ggsave("visualization/top6_species_map.png", width = 10, height = 8, dpi = 300)


# Plot Final 17
ggplot() +
  geom_sf(data = target_counties_sf, fill = "gray95", color = "black") +
  geom_sf(data = final17_sf, aes(color = species), alpha = 0.7, size = 1) +
  scale_color_viridis_d(option = "turbo") +
  labs(
    title = "Final 17 Plant Species",
    subtitle = "NAD83 / California State Plane Zone 5 (EPSG:2229)",
    color = "Species"
  ) +
  coord_sf(crs = 2229, default_crs = NA, datum = NA) +
  annotation_scale(
    location = "bl",        # bottom-left
    width_hint = 0.25,      # proportion of scale bar relative to width
    bar_cols = c("grey60", "white")
  ) +
  annotation_north_arrow(
    location = "tr",        # top-right
    which_north = "true",   # geographic north
    pad_x = unit(0.5, "cm"),
    pad_y = unit(0.5, "cm"),
    style = north_arrow_fancy_orienteering()
  ) +
  theme_bw()
ggsave("visualization/final17_species_map.png", width = 10, height = 8, dpi = 300)


