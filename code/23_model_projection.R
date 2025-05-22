library(SSDM)
#library(terra)
library(dplyr)
library(raster)
library(here)
library(spatialEco)
select <- dplyr::select

# ------------ 1. Get environmental and occurrence data ------------
# Environmental data
#Env <- raster::stack(here("data", "Stack_Env", "final_env_1980_2010_stack.tif"))


Env_files <- list.files(here("data", "Stack_Env"), 
                          pattern = "^final_env_ssp.*_stack\\.tif$", 
                          full.names = TRUE)

Env_normalized_list <- lapply(Env_files, function(f) {
  s <- raster::stack(f)
  
  s_std <- raster::stack(lapply(1:nlayers(s), function(i) {
    r <- s[[i]]
    mu <- cellStats(r, stat = 'mean', na.rm = TRUE)
    sigma <- cellStats(r, stat = 'sd', na.rm = TRUE)
    (r - mu) / sigma
  }))
  
  # Extract file name, e.g., "final_env_ssp126_2040_2070_stack.tif"
  base_name <- basename(f)
  
  # Extract scenario names
  ssp_tag <- sub(".*(ssp[0-9]+).*", "\\1", base_name)
  
  # Change layer names
  names(s_std) <- paste0(ssp_tag, "_", names(s), "_std")
  
  return(s_std)
})



#Env_normalized_list <- lapply(Env_files, function(f) {
#  s <- raster::stack(f)  # read raster stack
#  s_std <- raster::stack(lapply(1:nlayers(s), function(i)raster.Zscore(s[[i]])))  # normalized
#  names(s_std) <- paste0(names(s), "_std")  # change names
#  return(s_std)  # Return the normalized data
#})


occ <- read.csv(here("data", "occurrences", "Anim_Plant_merge.csv")) %>% 
  dplyr::rename("SPECIES" = species, 
         "LONGITUDE" = x, "LATITUDE" = y) %>% 
  select(-taxon) %>% 
  filter(SPECIES == "Agelaius tricolor")




# ------------ 2. Model setting ------------
mods <- ensemble_modelling(c('RF', "GAM", "MAXENT"), occ, Env, 
  Xcol = 'LONGITUDE', Ycol = 'LATITUDE', verbose = TRUE,
  cv = "LOO")

# ------------ 3. Function for loop ------------
model_projection <- function(mod, ){
  
  
  
}