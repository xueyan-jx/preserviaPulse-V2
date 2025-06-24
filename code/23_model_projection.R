## Purpose of script: Species Distribution Projection 
## Authors: GEOG 274
## Date: Spring, 2025
## Credits to: Wenxin Yang, Yanni Zhan, Xue Yan, Yifei Liu


#library(terra)
library(dplyr)
library(raster)
library(here)
library(spatialEco)
library(tryCatchLog)
library(mgcv)
library(dismo)
library(randomForest)
library(pROC)
library(caret)

select <- dplyr::select


# ------------ 1. Get environmental and species list ------------
# Environmental data
## Reference data
ref_Env <- raster::stack(here("data", "Stack_Env", "combined_stack_hist.tif"))

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


# ------------ 2. Function for model loop ------------
Predict_species <- function(Env_normalized_list,
                            model_dir = here("results", "models","mammal"),
                            projection_dir = here("results", "projections","mammal")) {
  
  # Get species list
  rds_files <- list.files(model_dir, pattern = "\\.RDS$", ignore.case = TRUE, full.names = FALSE)
  if (length(rds_files) == 0) {
    stop("No .rds model files found in the model directory.")
  }
  species_list <- sub("\\.RDS$", "", rds_files)
  
  # Output directories
  dir.create(projection_dir, recursive = TRUE, showWarnings = FALSE)

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
      
      proj_list <- list()
      
      # -----Get scenario name-----
      scenario_name <- attr(Env, "scenario_name")
      base_name <- paste0(sp, "_", scenario_name)
      cat("Projecting:", base_name, "\n")
      
      # ---- projection----
      # GAM
      if (!is.null(sp_model$gam)) {
        gam_proj <- tryCatch({
          predict(Env, sp_model$gam, type = "response") 
        }, error = function(e) {
          warning(paste("GAM projection failed for", sp, ":", e$message))
          NULL
        })
        proj_list[["gam"]] <- gam_proj
      }
      
      # Random Forest
      if (!is.null(sp_model$rf)) {
        rf_proj <- tryCatch({
          predict(Env, sp_model$rf, type = "prob", index = 2)
        }, error = function(e) {
          warning(paste("RF projection failed for", sp, ":", e$message))
          NULL
        })
        proj_list[["rf"]] <- rf_proj
      }
      
      # Maxent
      if (!is.null(sp_model$maxent)) {
        maxent_proj <- tryCatch({
          predict(Env, sp_model$maxent)
        }, error = function(e) {
          warning(paste("Maxent projection failed for", sp, ":", e$message))
          NULL
        })
        proj_list[["maxent"]] <- maxent_proj
      }
      
      # Stack output raster from different models
      valid_projs <- proj_list[!sapply(proj_list, is.null)]
      if (length(valid_projs) == 0) {
        warning(paste("No valid projections for", sp, "in", scenario_name))
        next
      }
      
      proj_stack <- stack(valid_projs)
      names(proj_stack) <- names(valid_projs)
      
      # Model ensemble
      proj_ensemble <- calc(proj_stack, fun = mean, na.rm = TRUE)
      
      # ---- save projection ----
      projection_path <- file.path(projection_dir, paste0(base_name, "_projection.tif"))
      writeRaster(proj_ensemble, filename = projection_path, format = "GTiff", overwrite = TRUE)
      
      } 
    }    
  }      


# ------------ 3. Run model for all target species ------------
Projection_result <- Predict_species(Env_normalized_list,
                                     model_dir = here("results", "models","mammal"),
                                     projection_dir = here("results", "projections","mammal"))



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
  
  # Calculate range: max - min
  r_max <- calc(r_stack, max, na.rm = TRUE)
  r_min <- calc(r_stack, min, na.rm = TRUE)
  r_range <- r_max - r_min
  
  # Create output dir
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Write to local
  writeRaster(r_mean, filename = file.path(output_dir, paste0(species_name, "_mean_scenarios.tif")),
              format = "GTiff", overwrite = TRUE)
  
  writeRaster(r_var, filename = file.path(output_dir, paste0(species_name, "_var_scenarios.tif")),
              format = "GTiff", overwrite = TRUE)
  
  writeRaster(r_range, filename = file.path(output_dir, paste0(species_name, "_range_scenarios.tif")),
              format = "GTiff", overwrite = TRUE)
  
  return(list(species = species_name, mean = r_mean, variance = r_var, Projrange = r_range))
}

# Apply function to all species
projection_dir <- "results/projections/mammal"
output_dir <- "results/evaluations/mammal"
species_list <- sub("\\.RDS$", "", list.files(here("results", "models","mammal"), pattern = "\\.RDS$", ignore.case = TRUE))
if (length(species_list) == 0) {
  stop("No .rds model files found in the model directory.")
}


uncertainty_results <- lapply(species_list, function(sp) {
  cat("Calculating uncertainty for:", sp, "\n")
  calculate_uncertainty(sp, projection_dir, output_dir)
})

# Check results
r <- rast(here("results", "evaluations", "plant","Abronia maritima_mean_scenarios.tif"))
plot(r)