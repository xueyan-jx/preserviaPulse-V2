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
library(tryCatchLog)
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
occ <- read.csv(here("data", "occurrences", "Anim_Plant_merge.csv")) 
species_list <- unique(occ$species)

#species_list <- list("Agelaius tricolor")


# ------------ 2. Function for model loop ------------
Predict_species <- function(species_list, Env_normalized_list,
                            model_dir = here("data", "models"),
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
    
    model_path <- file.path(model_dir, paste0(sp, ".rds"))
    
    # Check if model exists
    if (!file.exists(model_path)) {
      warning(paste("Model file not found for", sp))
      next
    }
    
    # Get model for each species (skip and record the failed ones)
    sp_model <- tryCatchLog(readRDS(model_path), error = function(e) {
      warning(paste("Failed to read model for", sp, ":", e$message))
      return(NULL)
    },write.error.dump.file = TRUE)
    
    if (is.null(sp_model)) next
    
    for (Env in Env_normalized_list) {
      
      # -----Get scenario name-----
      scenario_name <- attr(Env, "scenario_name")
      base_name <- paste0(sp, "_", scenario_name)
      
      # ---- projection----
      proj <- tryCatch({
        project(sp_model, Env, update.projections = FALSE, SDM.projections = TRUE)
      }, error = function(e) {
        warning(paste("Projection failed for", sp, "in", scenario_name, ":", e$message))
        return(NULL)
      })
      
      if (is.null(proj)) next
      
      
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

      # ---- save each submodel's projection ----
      for (i in seq_along(proj@sdms)) {
        submodel <- proj@sdms[[i]]
        sub_name <- submodel@name
        sub_base <- paste0(base_name, "_", sub_name)
        
        # Save projection for each model
        sub_proj_path <- file.path(submodel_dir, paste0(sub_base, "_projection.tif"))
        writeRaster(submodel@projection, filename = sub_proj_path, format = "GTiff", overwrite = TRUE)
      } 
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
    variable_importance = all_variable_importance,
    submodel_evaluation = all_submodel_evaluation,                         
    submodel_variable_importance = all_submodel_variable_importance 
  ))
} 


# ------------ 3. Run model for all target species ------------
Projection_result <- Predict_species(species_list, 
                                     Env_normalized_list,
                                     model_dir = here("data", "models"),
                                     projection_dir = here("results", "projections"),
                                     evaluation_dir = here("results", "evaluations"),
                                     submodel_dir = here("results", "submodelsEva"))



# -------------4. Uncertainty between different scenarios-------
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


# Check results
r<- rast(here("results", "projections", "Agelaius tricolor_ssp370_projection.tif"))
plot(r)