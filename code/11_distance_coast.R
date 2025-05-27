# ecu<-st_read("path to ECU shapefile")
# 
# ecu <- st_transform(ecu, crs="EPSG:2229")
# 
# variable_name_for_observations<-st_transform(variable_name_for_observations, crs="EPSG:2229") #change "variable_name_for_observations" to the variable name associated with the species observations
# 
# dist_to_coast <- sf::st_distance(variable_name_for_observations, ecu) #finds distance to coast (US survey foot)
# 
# min_distances <- apply(dist_to_coast, 1, min) #finds minimum distance to coast
# 
# variable_name_for_observations$distance_to_coast <- min_distances

# Load libraries
library(terra)
library(sf)
library(tigris)
library(here)

dat_dir <- here("data")

# Load a template to use
template <- rast(
  file.path(dat_dir, "Stack_Env", "final_env_1980_2010_stack.tif"))[[1]]

# Query the coastlines from tigris package, transform, rasterize, and 
# calculate the distance.
## Note: unit parameter does not work in function distance. The default is m.
dist_to_cstl <- coastline() %>% 
  st_transform(crs(template)) %>% 
  rasterize(template) %>% 
  terra::distance(x = ., unit = "km", method = "geo") %>% 
  mask(template)

## Deal with the bug in distance function
dist_to_cstl <- dist_to_cstl / 1000 # to km

fnames <- list.files(
  file.path(dat_dir, "Stack_Env"), full.names = TRUE)

# Replace the layer in the Env stacks
for(fname in fnames){
  lyrs <- rast(fname)
  lyrs$dist_coast <- dist_to_cstl
  
  msk <- sum(lyrs, na.rm = FALSE)
  lyrs <- mask(lyrs, msk)
  
  writeRaster(lyrs, gsub("Stack_Env", "Stack_Env_fix", fname))
}

# If everything is okay, remove old files and rename the fix folder
unlink(file.path(dat_dir, "Stack_Env"), recursive = TRUE)
file.rename(file.path(dat_dir, "Stack_Env_fix"), 
            file.path(dat_dir, "Stack_Env"))
