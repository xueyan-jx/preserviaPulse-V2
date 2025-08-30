## Purpose of script: 
## 1. Species Distribution Projection 
## 2. Richness calculation for each species group
## 3. Uncertainty between scenarios
## 4. Write variable importance from the models
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
library(readxl)
library(ggplot2)
library(tibble)
library(stringr)
library(writexl)
library(tools)

select <- dplyr::select


# ------------ 1. Get environmental data ------------
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
                            model_dir = here("results", "models","Rerun"),
                            projection_dir = here("results", "projections","Rerun")) {
  
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
        gam_bin <- (gam_proj > as.numeric(sp_model$gam_thresh))*1 # as.integer should also work
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
                                     model_dir = here("results", "models","Rerun"),
                                     projection_dir = here("results", "projections","Rerun"))


# -------------4. Richness calculation ------------------------
## -------------(1) Richness calculation ------------------------
calculate_species_richness <- function(projection_dir, 
                                       scenario = c("ssp126", "ssp370", "ssp585"), 
                                       output_dir = here("results", "projections","Rerun")) {
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

calculate_species_richness(projection_dir = here("results", "projections","Habitat","Herb","plant"), 
                           scenario = c("ssp126", "ssp370", "ssp585"), 
                           output_dir = here("results", "projections","Habitat","Herb","plant"))

# r <- terra::rast(here("results", "projections","Habitat","Herb","plant", "Plant_species_richness_ssp126.tif"))
# plot(r)  

## -------------(2) Current-future difference ------------------------
# r1 <- rast("raster1.tif")
# r2 <- rast("raster2.tif")
# 
# r_diff <- r1 - r2
# 
# writeRaster(r_diff, here("output", "raster_diff.tif"), overwrite = TRUE)



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

# -------------6. Write variable importance from the models-------
# Order variables
var_order <- c("bio1", "bio3", "bio5", "bio6", "bio15", "bio17", "bio18", "bio19", 
               "slope", "aspect", "flow_acc", "solar")


# Function for cleaning column names
clean_varname <- function(x) {
  x %>%
    str_replace_all("^s\\((.+)\\)$", "\\1") %>%      # delete s() in gam 
    str_replace_all("\\.contribution$", "") %>%     # delete .contribution in maxent
    str_trim()
}

# Function for processing varimp tables
process_varimp_unified <- function(df) {
  df2 <- df %>%
    select(2) %>%                         # Extract Importance 
    rownames_to_column(var = "Variable") %>%  # Convert row name to column
    rename(Importance = 2) %>%
    mutate(Variable = clean_varname(Variable)) %>%
    filter(Variable %in% var_order) %>%
    mutate(Variable = factor(Variable, levels = var_order)) %>%
    arrange(Variable)
  return(df2)
}

# Read models
rds_dir <- here("results", "models","Rerun")
rds_files <- list.files(rds_dir, pattern = "\\.RDS$", ignore.case = TRUE, full.names = FALSE)


all_data <- list()  

for (rds_file in rds_files) {
  rds_path <- file.path(rds_dir, rds_file)   
  Targetmodel <- readRDS(rds_path)
  species <- file_path_sans_ext(basename(rds_file))
  
  if (!is.null(Targetmodel[["rf_varimp"]])) {
    df_rf <- process_varimp_unified(Targetmodel[["rf_varimp"]]) %>%
      mutate(species = species, model = "rf")
    all_data[[length(all_data) + 1]] <- df_rf
  }
  
  if (!is.null(Targetmodel[["maxent_varimp"]])) {
    df_maxent <- process_varimp_unified(Targetmodel[["maxent_varimp"]]) %>%
      mutate(species = species, model = "maxent")
    all_data[[length(all_data) + 1]] <- df_maxent
  }
  
  if (!is.null(Targetmodel[["gam_varimp"]])) {
    df_gam <- process_varimp_unified(Targetmodel[["gam_varimp"]]) %>%
      mutate(species = species, model = "gam")
    all_data[[length(all_data) + 1]] <- df_gam
  }
}

if (length(all_data) > 0) {
  combined_df <- bind_rows(all_data) %>%
    select(species, model, Variable, Importance)
  
  out_file <- file.path(rds_dir, "all_species_varimp_combined.xlsx")
  write_xlsx(list(all_varimp = combined_df), out_file)
  message("Saved to：", out_file)
} else {
  message("No more data to be written")
}


# Plot
eval_dir <- here("results", "evaluations", "KeyHabitat","Shrub")
file_path <- file.path(eval_dir, "Shrub_variable_importance.xlsx")

allSpecies_df <- read_excel(file_path, sheet = "allSpecies")
bird_df <- read_excel(file_path, sheet = "Bird")
plant_df <- read_excel(file_path, sheet = "plant")
other_animal_df <- read_excel(file_path, sheet = "animal_on_the_ground")

# allSpecies_df$group <- "All species" 
# bird_df$group <- "Bird"
# plant_df$group <- "Plant"
# other_animal_df$group <- "Other animal"
# df_long <- bind_rows(allSpecies_df, bird_df, plant_df, other_animal_df)
# 
# df_long <- bind_rows(allSpecies_df, bird_df, plant_df, other_animal_df) %>%
#   distinct(species, model, Variable, .keep_all = TRUE)

plant_vars <- plant_df %>% select(species, model, Variable, Importance) %>% mutate(group_tag = "Plant")
bird_vars <- bird_df %>% select(species, model, Variable, Importance) %>% mutate(group_tag = "Bird")
other_animal_vars <- other_animal_df %>% select(species, model, Variable, Importance) %>% mutate(group_tag = "Other animal")

df_long <- bind_rows(plant_vars, bird_vars, other_animal_vars) %>%
  mutate(all_species_flag = "All species")


df_long <- df_long %>%
  group_by(species, model, group_tag) %>%
  mutate(
    Importance_rank = rank(-Importance, ties.method = "average"),
    Importance_norm = (Importance - min(Importance, na.rm = TRUE)) / 
      (max(Importance, na.rm = TRUE) - min(Importance, na.rm = TRUE))
  ) %>%
  ungroup()

all_species_df <- df_long %>%
  mutate(group_tag = "All species")

df_long <- bind_rows(df_long, all_species_df)

df_long %>% count(species, model)


summary_df <- df_long %>%
  group_by(group_tag, Variable) %>%
  summarise(
    median = median(Importance_rank, na.rm = TRUE),
    q1 = quantile(Importance_rank, 0.25, na.rm = TRUE),
    q3 = quantile(Importance_rank, 0.75, na.rm = TRUE),
    .groups = "drop"
  )


ggplot(summary_df, aes(x = Variable, y = median, color = group_tag)) +
  geom_errorbar(aes(ymin = q1, ymax = q3), 
                position = position_dodge(width = 0.6),
                width = 0.2, linewidth = 0.8)  +
  
  geom_point(position = position_dodge(width = 0.6),
             size = 2.5)   +
  labs(
    title = "Model-averaged variable importance across species groups",
    x = "Variable",
    y = "Importance (median with IQR)"
  ) +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_color_manual(
    name = "Species Group",
    values = c(
    "All species" = "#4a4e69",
    "Bird" = "#FFCA3A",
    "Plant" = "#8AC926",
    "Other animal" = "#1982C4"
  ))  +
  theme_minimal(base_size = 14) +
  theme(
    axis.line = element_line(color = "black", size = 0.6),
    axis.ticks = element_line(color = "black", linewidth = 0.4),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )


# across_species_summary <- per_species_mean %>%
#   group_by(Variable) %>%
#   summarise(
#     mean_of_means = mean(mean_importance, na.rm = TRUE),
#     sd_of_means = sd(mean_importance, na.rm = TRUE),
#     min_of_means = min(mean_importance, na.rm = TRUE),
#     max_of_means = max(mean_importance, na.rm = TRUE),
#     n_species = n(),
#     .groups = "drop"
#   ) %>%
#   arrange(desc(mean_of_means))

