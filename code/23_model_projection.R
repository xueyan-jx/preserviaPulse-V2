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
      binary_rasters <- list() 
      
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
        
        
        # GAM binary
        gam_bin <- (gam_proj > as.numeric(sp_model$gam_thresh))*1
        gam_bin <- mask(gam_bin, gam_proj)
        binary_rasters[["gam"]] <- gam_bin
        
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
        
        # RF binary
        rf_bin <- (rf_proj > as.numeric(sp_model$rf_thresh))*1
        rf_bin <- mask(rf_bin, rf_proj)
        binary_rasters[["rf"]] <- rf_bin
        
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
        
        # Maxent binary
        maxent_bin <- (maxent_proj > as.numeric(sp_model$maxent_thresh))*1
        maxent_bin <- mask(maxent_bin, maxent_proj)
        binary_rasters[["maxent"]] <- maxent_bin
        
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
      
      # ---- save binary projection ---
      if (length(binary_rasters) >= 2) {
        binary_stack <- stack(binary_rasters)
        binary_sum <- calc(binary_stack, sum)
        ensemble_binary <- (binary_sum >= 2)*1
        ensemble_binary <- mask(ensemble_binary, binary_sum)
        ensemble_binary_path <- file.path(projection_dir, paste0(base_name, "_ensemble_binary.tif"))
        writeRaster(ensemble_binary, filename = ensemble_binary_path, format = "GTiff", overwrite = TRUE)
      }
      
      rm(gam_proj, rf_proj, maxent_proj, gam_bin, rf_bin, maxent_bin)
      
      } 
    }    
  }      


# ------------ 3. Run model for all target species ------------
Projection_result <- Predict_species(Env_normalized_list,
                                     model_dir = here("results", "models","mammal"),
                                     projection_dir = here("results", "projections","mammal"))


# -------------4. Richness calculation ------------------------
calculate_species_richness <- function(projection_dir, 
                                       scenario = c("ssp126", "ssp370", "ssp585"), 
                                       output_dir = here("results", "projections","mammal")) {
  richness_list <- list()
  
  for (s in scenario){
    cat("Processing scenario:", s, "\n")
    
    # get ensemble binary files
    bin_files <- list.files(projection_dir, 
                            pattern = paste0("_", s, "_ensemble_binary\\.tif$"), 
                            full.names = TRUE)
    
    if (length(bin_files) == 0) {
      stop("No binary maps found for scenario: ", s)
    }
    
    # Read all binary files
    bin_stack <- stack(bin_files)
    
    # sum
    richness_raster <- calc(bin_stack, sum, na.rm = TRUE)
    richness_raster <- mask(richness_raster, bin_stack[[1]])
    
    output_path <- file.path(output_dir, paste0("species_richness_", s, ".tif"))
    
    writeRaster(richness_raster, 
                filename = output_path, 
                format = "GTiff", 
                datatype = "INT2U", 
                overwrite = TRUE)
    cat("Saved richness raster:", output_path, "\n")
    
    richness_list[[s]] <- richness_raster
    
  }
  return(richness_list)
}

calculate_species_richness(projection_dir = here("results", "projections","mammal"), 
                           scenario = c("ssp126", "ssp370", "ssp585"), 
                           output_dir = here("results", "projections","mammal"))
  
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