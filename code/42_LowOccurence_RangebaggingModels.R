## Purpose of Script: Range bagging Models for Rare Occurrence Data
## Authors: GEOG 274
## Date: Spring, 2025 
## Credits to: Olivia Ross, Dr. Lei Song

# packages
librarian::shelf(SSDM, terra, raster, dplyr, sf, tigris, here, bssdm)

# installing bssdm from github
remotes::install_github("cebra-analytics/bssdm")
# Source: https://github.com/cebra-analytics/bssdm/tree/main/R

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
# READING IN OUR DATA FOR MODELING

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
# cleaned rare occurrences
rare_occurences <- read.csv(here("stored data from script 41"))

# reading in climate stacks 
clim <- rast(here("data/Stack_Env/final_env_1980_2010_stack.tif"))
s126 <- rast(here("data/Stack_Env/final_env_ssp126_2040_2070_stack.tif"))
s370 <- rast(here("data/Stack_Env/final_env_ssp370_2040_2070_stack.tif"))
s585 <- rast(here("data/Stack_Env/final_env_ssp585_2040_2070_stack.tif"))

# target boundaries for mapping purposes
counties <- counties(state = "CA", cb = TRUE) 

target_counties <- counties |>
  filter(NAME %in% c("Santa Barbara", "Ventura", "San Luis Obispo"))

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
# CLEANING DATA FOR OUR MODELS 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
# occurrence data needs to be in WGS 84 projection with lat and lon columns 
rare_clean <- rare_occurences |>
  rename(lon = x, lat = y) |>
  st_as_sf(coords = c("lon", "lat"), crs = 2229) |>
  st_transform(crs = 4326) 

# grabbing coordinates for their separate columns 
rare_coords <- st_coordinates(rare_clean)

rare_added_coords <- cbind(rare_clean, rare_coords)

# renaming columns and dropping geometry
rare_final <- rare_added_coords |>
  rename(lon = X, lat = Y) |>
  as.data.frame() |>
  dplyr::select(-geometry)

# making sure this is correct 
plot(target_counties$geometry, col = NA)
plot(rare_clean$geometry, add=TRUE)

# cleaning climate data - 
# requires any crs or WGS84
crs(clim)

# reprojecting climate data
clim.stacks <- list(s126 = s126, s370 = s370, s585 = s585, clim = clim)

names(clim.stacks)

clim.projected <- lapply(clim.stacks, function(x) {
  terra::project(x, "EPSG:4326")
})

# making sure names are preserved 
names(clim.projected) <- paste0(names(clim.stacks), ".wgs84")

# adding to the environment
list2env(clim.projected, envir = .GlobalEnv)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
# RUNNING THE RANGEBAG MODEL
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
# running rangebag model on all species to test 
rare.rangebag <- rangebag(clim.wgs84, rare_final, n_models = 100, 
                          n_dim = 2, 
                          sample_prop = 0.7,
                          limit_occur = FALSE)

# now running prediction model on all species to test 

current <- predict(rare.rangebag, clim.wgs84)

pred.126 <- predict(rare.rangebag, s126.wgs84)

pred.370 <- predict(rare.rangebag, s370.wgs84)

pred.585 <- predict(rare.rangebag, s585.wgs84)

plot(current)
plot(pred.126)
plot(pred.370)
plot(pred.585)

# now modeling each individual species

# seeing how many occurrences of each species we have 
rare_counts <- rare_final |>
  group_by(species, taxon) |>
  count()

View(rare_counts)

# creating a species df for each species to model 

# willow flycatcher
Empidonax_traillii_extimus <- rare_final |>
  filter(species == "Empidonax traillii extimus")

# peregrine falcon
Falco_peregrinus_anatum <- rare_final |>
  filter(species == "Falco peregrinus anatum")

# northern harrier
Circus_hudsonius <- rare_final |>
  filter(species == "Circus hudsonius")

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
# LOOPING THROUGH EACH SPECIES AND CREATING HABITAT SUITABILITY MAPS FOR EACH 
# CLIMATE SCENARIO
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
# modeling each species on a loop 
species <- list(Circus_hudsonius = Circus_hudsonius,
                Empidonax_traillii_extimus = Empidonax_traillii_extimus,
                Falco_peregrinus_anatum = Falco_peregrinus_anatum)

names(species)

# storing our outputs by climate scenario
output.folder.current <- here("data/low_occurence_current")
output.folder.126 <- here("data/low_occurence_ssp126")
output.folder.370 <- here("data/low_occurence_ssp370")
output.folder.585 <- here("data/low_occurence_ssp585")

for (i in seq_along(species)) {
 
  # reading in each df and grabbing the name  
  indiv <- species[[i]]
  
  name <- names(species)[i]
  
  # rangebag for each species 
  rangebagi <- rangebag(clim.wgs84, indiv, n_models = 100, 
                            n_dim = 2, 
                            sample_prop = 0.7,
                            limit_occur = FALSE)
 
  # habitat suitability predictions for each scenario 
  currenti <- predict(rangebagi, clim.wgs84)
  
  s126i <- predict(rangebagi, s126.wgs84)
  
  s370i <- predict(rangebagi, s370.wgs84)
  
  s585i <- predict(rangebagi, s585.wgs84)
  
  # creating our file paths and names for each scenario
  filename.current <- file.path(output.folder.current, paste0(name, "_current.tif"))
  filename.126 <- file.path(output.folder.126, paste0(name, "_ssp126.tif"))
  filename.370 <- file.path(output.folder.370, paste0(name, "_ssp370.tif"))
  filename.585 <- file.path(output.folder.585, paste0(name, "_ssp585.tif"))
  
  # writing our rasters for each scenario
  writeRaster(currenti, filename.current, overwrite = TRUE)
  writeRaster(s126i, filename.126, overwrite = TRUE)
  writeRaster(s370i, filename.370, overwrite = TRUE)
  writeRaster(s585i, filename.585, overwrite = TRUE)
  
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
# VIEWING RASTERS & CHECKING OUR WORK 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
# northern harrier
nh.current <- rast(here("data/low_occurence_current/Circus_hudsonius_current.tif"))
nh.s126 <- rast(here("data/low_occurence_ssp126/Circus_hudsonius_ssp126.tif"))
nh.s370 <- rast(here("data/low_occurence_ssp370/Circus_hudsonius_ssp370.tif"))
nh.s585 <- rast(here("data/low_occurence_ssp585/Circus_hudsonius_ssp585.tif"))

# plotting results
plot(nh.current, main = "Northern Harrier Habitat Suitability in Current Climate")
plot(nh.s126, main = "Northern Harrier Habitat Suitability in SSP 126 Climate Scenario")
plot(nh.s370, main = "Northern Harrier Habitat Suitability in SSP 370 Climate Scenario")
plot(nh.s585, main = "Northern Harrier Habitat Suitability in SSP 585 Climate Scenario")

plot(target_counties, col = NA, add = TRUE)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
# MODEL UNCERTAINTIES - LOADING DATA
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
library("pROC")

# adding our presence data to a list 
presences <- list(Empidonax_traillii_extimus = Empidonax_traillii_extimus,
                  Falco_peregrinus_anatum = Falco_peregrinus_anatum,
                  Circus_hudsonius = Circus_hudsonius)

species_names <- c("Empidonax_traillii_extimus", "Falco_peregrinus_anatum",
                   "Circus_hudsonius")

names(presences) <- species_names

# reading in sdms for each species
scenario_folders <- c("data/low_occurence_current", "data/low_occurence_ssp126", "data/low_occurence_ssp370", 
             "data/low_occurence_ssp585")

n_species <- length(presences)
n_scenarios <- length(scenario_folders)

# creating raster lists: 1 list per species, each containing 4 rasters (current + one per climate scenario)
raster_lists <- vector("list", n_species)
names(raster_lists) <- species_names  

# looping over each folder and assigning each SDM output to each species
for (sp in seq_along(species_names)) {
  
  # grabbing our species names 
  species <- species_names[sp]
  
  # creating our lists of species and each SDM scenario
  species_rasters <- vector("list", n_scenarios)
  
  # looping over each scenario folder
  for (sc in seq_along(scenario_folders)) {
    folder <- scenario_folders[sc]
    
    # matching raster file for species in this folder 
    files <- list.files(folder, pattern = "\\.tif$", full.names = TRUE)
    match_file <- grep(species, files, ignore.case = TRUE, value = TRUE)
    
    # loading the matched file 
    species_rasters[[sc]] <- rast(match_file)
  }

  # appending to raster list 
  raster_lists[[sp]] <- species_rasters
}


# checking our work 
sapply(raster_lists, length)
r <- raster_lists[["Falco_peregrinus_anatum"]][[1]]

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
# MODEL UNCERTAINTIES - AUC SCORES & TRUE SKILL STATISTIC
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
n_species <- length(presences)
n_scenarios <- length(scenario_folders)
species_names <- names(presences)

# constructing our AUC matrix and adding names for rows and columns
auc_matrix <- matrix(NA, nrow = n_species, ncol = n_scenarios,
                     dimnames = list(species_names, paste0("Scenario_", LETTERS[1:n_scenarios])))

# creating our AUC matrix
auc_matrix <- matrix(NA, nrow = n_species, ncol = n_scenarios,
                     dimnames = list(species_names, basename(scenario_folders)))

# creating our TSS matrix
tss_matrix <- matrix(NA, nrow = n_species, ncol = n_scenarios,
                     dimnames = list(species_names, basename(scenario_folders)))

# setting seed for background points
set.seed(142)

# looping over each species df and sdm scenario
for (i in seq_along(species_names)) {
  
  # reading in each species df 
  sp_name <- species_names[i]
  pres_df <- presences[[sp_name]]

# calculating an AUC score for each SDM scenario, per species     
  for (j in seq_len(n_scenarios)) {
    pred_raster <- raster_lists[[sp_name]][[j]]
    
    # extracting suitability values at presence points
    presence_vals <- extract(pred_raster, pres_df[, c("lon", "lat")])[, 2]
    presence_vals2 <- presence_vals[!is.na(presence_vals)]
    
    # sampling 1000 background points from raster
    bg_sample <- terra::spatSample(pred_raster, size = 1000, method = "random", values = TRUE)
    background_vals <- bg_sample[[1]]
    background_vals2 <- background_vals[!is.na(background_vals)]
    
    # Combining presence and background predictions
    all_preds <- c(presence_vals2, background_vals2)
    
    # labeling points, 1 for presence and 0 for absence (background)
    labels <- c(rep(1, length(presence_vals2)), rep(0, length(background_vals2)))
    
    stopifnot(length(all_preds) == length(labels))
    
    # calculating ROC and AUC
    roc_obj <- roc(response = labels, predictor = all_preds)
    auc_val <- auc(roc_obj)
    
    # adding to our matrix
    auc_matrix[i, j] <- auc_val
    
    # CALCULATING TRUE SKILL STATISTIC
    
    # finding the optimal threshold
    opt_thresh <- pROC::coords(roc_obj, x = "best", best.method = "youden", transpose = TRUE)[1]
    
    # Binarize predictions
    binary_preds <- as.numeric(all_preds >= opt_thresh)
    
    # Confusion matrix components
    TP <- sum(binary_preds == 1 & labels == 1)
    TN <- sum(binary_preds == 0 & labels == 0)
    FP <- sum(binary_preds == 1 & labels == 0)
    FN <- sum(binary_preds == 0 & labels == 1)
    
    # Sensitivity and Specificity
    sensitivity <- TP / (TP + FN)
    specificity <- TN / (TN + FP)
    
    # TSS
    tss_val <- sensitivity + specificity - 1
    tss_matrix[i, j] <- tss_val
  }
}

# printing the final tables
print(auc_matrix)
print(tss_matrix)