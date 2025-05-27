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


# ------------ 1. Get environmental and species list ------------
# Environmental data
## Reference data
ref_Env <- raster::stack(here("data", "Stack_Env", "final_env_1980_2010_stack.tif"))

## Projection data
Env_files <- list.files(here("data", "Stack_Env"), 
                          pattern = "^final_env_ssp.*_stack\\.tif$", 
                          full.names = TRUE)
## Normalization 
Env_normalized_list <- lapply(Env_files, function(f) {
  s <- raster::stack(f)
  
  s_std <- raster::stack(lapply(1:nlayers(s), function(i) {
    r <- s[[i]]  
    ref <- ref_Env[[i]] 
    
    mu <- cellStats(ref, stat = 'mean', na.rm = TRUE)
    sigma <- cellStats(ref, stat = 'sd', na.rm = TRUE)
    (r - mu) / sigma
  }))
  
  # Extract file name, e.g., "final_env_ssp126_2040_2070_stack.tif"
  base_name <- basename(f)
  
  # Extract scenario names
  ssp_tag <- sub(".*(ssp[0-9]+).*", "\\1", base_name)
  
  # Change layer names
  # names(s_std) <- paste0(ssp_tag, "_", names(s), "_std")
  
  # Add attribute: scenario_name
  attr(s_std, "scenario_name") <- ssp_tag
  return(s_std)
})

# Get species list (need to change)
species_list <- list("Agelaius tricolor")

# ------------ 2. Model setting ------------
model_dir <- here("data", "models")
model_files <- list.files(model_dir, pattern = "\\.rds$", full.names = TRUE)

# read as list
ESDM <- setNames(
  lapply(model_files, readRDS),
  tools::file_path_sans_ext(basename(model_files))
)

#SDM_ens <- readRDS(here("data", "Agelaiustricolor_ntree500.rds"))

# ------------ 3. Function for model loop ------------
model_species_by_env_list <- function(ESDM, species_list, Env_normalized_list,
                                      projection_dir = here("results", "projections"),
                                      evaluation_dir = here("results", "evaluations"),
                                      submodel_dir = here("results", "submodelsEva")) {
  # Output directories
  dir.create(projection_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(evaluation_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(submodel_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Output variables
  all_algorithm_evaluation <- list()
  all_variable_importance <- list()
  all_submodel_evaluation <- list()            
  all_submodel_variable_importance <- list()
  
  for (sp in species_list) {
    cat("=== Modeling species:", sp, "===\n")
    
    for (Env in Env_normalized_list) {
      
      # -----Get scenario name-----
      scenario_name <- attr(Env, "scenario_name")
      base_name <- paste0(sp, "_", scenario_name)
      
      # -----Get model for each species----
      sp_model <- ESDM[[sp]]
      
      # ---- projection----
      proj <- project(ESDM, Env, update.projections = FALSE, SDM.projections = TRUE)
      
      # ---- save projection ----
      projection_path <- file.path(projection_dir, paste0(base_name, "_projection.tif"))
      writeRaster(proj@projection, filename = projection_path, format = "GTiff", overwrite = TRUE)
      
      # ---- save uncertainty ----
      uncertainty_path <- file.path(evaluation_dir, paste0(base_name, "_uncertainty.tif"))
      writeRaster(proj@uncertainty, filename = uncertainty_path, format = "GTiff", overwrite = TRUE)
      
      # ---- save overall algorithm.evaluation ----
      alg_eval <- as.data.frame(proj@algorithm.evaluation)
      alg_eval$Species <- sp
      alg_eval$Scenario <- scenario_name
      all_algorithm_evaluation[[length(all_algorithm_evaluation) + 1]] <- alg_eval
      
      # ---- save overall variable.importance ----
      var_imp <- as.data.frame(proj@variable.importance)
      var_imp$Species <- sp
      var_imp$Scenario <- scenario_name
      all_variable_importance[[length(all_variable_importance) + 1]] <- var_imp
      
      # ---- save each submodel's projection, evaluation, and variable importance ----
      for (i in seq_along(proj@sdms)) {
        submodel <- proj@sdms[[i]]
        sub_name <- submodel@name
        sub_base <- paste0(base_name, "_", sub_name)
        
        # Save projection for each model
        sub_proj_path <- file.path(submodel_dir, paste0(sub_base, "_projection.tif"))
        writeRaster(submodel@projection, filename = sub_proj_path, format = "GTiff", overwrite = TRUE)
        
        # Save evaluation (CSV)
        sub_eval <- as.data.frame(submodel@evaluation)
        sub_eval$Species <- sp
        sub_eval$Scenario <- scenario_name
        sub_eval$Submodel <- sub_name 
        sub_eval_path <- file.path(submodel_dir, paste0(sub_base, "_evaluation.csv"))
        all_submodel_evaluation[[length(all_submodel_evaluation) + 1]] <- sub_eval
        
        # Save variable importance (CSV)
        sub_var_imp <- as.data.frame(submodel@variable.importance)
        sub_var_imp$Species <- sp
        sub_var_imp$Scenario <- scenario_name
        sub_var_imp$Submodel <- sub_name 
        sub_var_imp_path <- file.path(submodel_dir, paste0(sub_base, "_var_importance.csv"))
        all_submodel_variable_importance[[length(all_submodel_variable_importance) + 1]] <- sub_var_imp
      } 
    }    
  }      
  
  # ---- Save all forms ----
  write.csv(do.call(rbind, all_algorithm_evaluation), 
            file = file.path(evaluation_dir, "all_algorithm_evaluation.csv"), row.names = FALSE)
  
  write.csv(do.call(rbind, all_variable_importance), 
            file = file.path(evaluation_dir, "all_variable_importance.csv"), row.names = FALSE)
  
  write.csv(do.call(rbind, all_submodel_evaluation),             
            file = file.path(evaluation_dir, "all_submodel_evaluation.csv"), row.names = FALSE)
  
  write.csv(do.call(rbind, all_submodel_variable_importance),   
            file = file.path(evaluation_dir, "all_submodel_variable_importance.csv"), row.names = FALSE)
  
  # ---- Return all forms as a list----
  return(list(
    algorithm_evaluation = all_algorithm_evaluation,
    variable_importance = all_variable_importance,
    submodel_evaluation = all_submodel_evaluation,                         
    submodel_variable_importance = all_submodel_variable_importance 
  ))
} 


# ------------ 4. Run model for all target species ------------
Projection_result <- model_species_by_env_list(ESDM, 
                                               species_list, 
                                               Env_normalized_list,
                                               projection_dir = here("results", "projections"),
                                               evaluation_dir = here("results", "evaluations"),
                                               submodel_dir = here("results", "submodelsEva"))



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

