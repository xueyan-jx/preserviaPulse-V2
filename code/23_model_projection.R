## Purpose of script: Species Distribution Projection 
## Authors: GEOG 274
## Date: Spring, 2025
## Credits to: Wenxin Yang, Yanni Zhan, Xue Yan, Yifei Liu

library(SSDM)
#library(terra)
library(dplyr)
library(raster)
library(here)
library(spatialEco)
select <- dplyr::select

# ------------ Output setting ------------
# Output folder
#dir.create("results/projections", recursive = TRUE, showWarnings = FALSE)
#dir.create("results/evaluations", recursive = TRUE, showWarnings = FALSE)


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
  
  # Add attribute: scenario_name
  attr(s_std, "scenario_name") <- ssp_tag
  return(s_std)
})

# Occurrence data
## Read data
occ <- read.csv(here("data", "occurrences", "Anim_Plant_merge.csv")) 
occ_plant <- read.csv(here("data", "occurrences", "Plant_occurrences_Top19.csv")) # Only model the most important 19 species with more than 10 records

## Filter to target species
occ_filtered <- occ %>%
  filter(
    taxon == "birds" |
    taxon == "herps" |
    taxon == "mammals" |
      (taxon == "plant" & species %in% occ_plant$species)
  )

## Check if plant filter is correct
occ_filtered %>%
  filter(taxon == "plant") %>%
  summarise(n_species = n_distinct(species))

## Organize data to fit function ensemble_modelling
occ_organized <- occ_filtered %>%
  dplyr::rename("SPECIES" = species, 
                "LONGITUDE" = x, "LATITUDE" = y) %>% 
  select(-taxon)

## Get species list
species_list <- unique(occ_organized$SPECIES)

# ------------ 2. Model setting ------------
ESDM <- ensemble_modelling(c('RF', "GAM", "MAXENT"), occ, Env, 
  Xcol = 'LONGITUDE', Ycol = 'LATITUDE', verbose = TRUE,
  cv = "LOO")

# ------------ 3. Function for loop ------------
model_species_by_env_list <- function(ESDM, species_list, occ_organized, Env_normalized_list,
                                      projection_dir = here("results", "projections"),
                                      evaluation_dir = here("results", "evaluations")) {
  # Output directories
  dir.create(projection_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(evaluation_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Output variable
  all_algorithm_evaluation <- list()
  all_variable_importance <- list()
  
  for (sp in species_list) {
    cat("=== Modeling species:", sp, "===\n")
    occ_sp <- occ_organized %>% filter(SPECIES == sp) # Occurrence for each species
    
    for (Env in Env_normalized_list) {
  
      # -----Get scenario name-----
      scenario_name <- attr(Env, "scenario_name")
      base_name <- paste0(sp, "_", scenario_name)
      
      # ---- projection----
      proj <- project(ESDM, Env, update.projections = FALSE, SDM.projections = TRUE)
      
      # ---- save projection ----
      projection_path <- file.path(projection_dir, paste0(base_name, "_projection.tif"))
      writeRaster(proj@projection, filename = projection_path, format = "GTiff", overwrite = TRUE)
      
      # ---- save uncertainty ----
      uncertainty_path <- file.path(evaluation_dir, paste0(base_name, "_uncertainty.tif"))
      writeRaster(proj@uncertainty, filename = uncertainty_path, format = "GTiff", overwrite = TRUE)
      
      # ---- algorithm.evaluation ----
      alg_eval <- as.data.frame(proj@algorithm.evaluation)
      alg_eval$Species <- sp
      alg_eval$Scenario <- scenario_name
      all_algorithm_evaluation[[length(all_algorithm_evaluation) + 1]] <- alg_eval
      
      # ---- variable.importance ----
      var_imp <- as.data.frame(proj@variable.importance)
      var_imp$Species <- sp
      var_imp$Scenario <- scenario_name
      all_variable_importance[[length(all_variable_importance) + 1]] <- var_imp

    }
  }
  
  # ---- Save all forms ----
  write.csv(do.call(rbind, all_algorithm_evaluation), 
            file = file.path(evaluation_dir, "all_algorithm_evaluation.csv"), row.names = FALSE)
  
  write.csv(do.call(rbind, all_variable_importance), 
            file = file.path(evaluation_dir, "all_variable_importance.csv"), row.names = FALSE)
  
  # ---- Return all forms as a list----
  return(list(
    algorithm_evaluation = all_algorithm_evaluation,
    variable_importance = all_variable_importance
  ))
}


# ------------ 4. Run model for all target species ------------
Projection_result <- model_species_by_env_list(ESDM, 
                                               species_list, 
                                               occ_organized, 
                                               Env_normalized_list,
                                               projection_dir = here("results", "projections"),
                                               evaluation_dir = here("results", "evaluations"))



# -------------5. Uncertainty between different scenarios-------
# Define function
calculate_uncertainty <- function(species_name, projection_dir, output_dir) {

  ssp_list <- c("ssp126", "ssp370", "ssp585")
  r_files <- file.path(projection_dir, paste0(species_name, "_", ssp_list, "_projection.tif"))
  
  # Check if file exists
  if (!all(file.exists(r_files))) {
    warning(paste("Missing projections for species:", species_name))
    return(NULL)
  }
  
  # Read the three rasters and merge to stack
  r_stack <- stack(r_files)
  
  # Calculate mean and variance
  r_mean <- calc(r_stack, mean, na.rm = TRUE)
  r_sd <- calc(r_stack, sd, na.rm = TRUE)
  r_var <- r_sd^2
  
  # Create output dir
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Write to local
  writeRaster(r_mean, filename = file.path(output_dir, paste0(species_name, "_mean_scenarios.tif")),
              format = "GTiff", overwrite = TRUE)
  
  writeRaster(r_var, filename = file.path(output_dir, paste0(species_name, "_var_scenarios.tif")),
              format = "GTiff", overwrite = TRUE)
  
  return(list(species = species_name, mean = r_mean, variance = r_var))
}

# Apply function to all species
projection_dir <- "results/projections"
output_dir <- "results/evaluations"

uncertainty_results <- lapply(species_list, function(sp) {
  cat("Calculating uncertainty for:", sp, "\n")
  calculate_uncertainty(sp, projection_dir, output_dir)
})

