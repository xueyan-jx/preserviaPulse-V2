## Purpose of script: 
## 1. Model uncertainty between scenarios
## 2. Write variable importance and AUC from the models
## 3. Plot percentage variable importance

## Authors: GEOG 274
## Date: Spring, 2025
## Credits to:  Xue Yan

library(tidyverse)   
library(readxl)
library(here)
library(scales)
library(raster)      
library(terra)       
library(tools)       


# -------------1. Uncertainty between different scenarios-------
## ----- (1) Define function----
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

## ----- (2) Apply function to all species ------
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
# r <- rast(here("results", "evaluations", "plant","Abronia maritima_mean_scenarios.tif"))
# plot(r)

# -------------2. Write variable importance and AUC from the models-------
# Order variables
# var_order <- c("bio1", "bio3", "bio5", "bio6", "bio14", "bio15", "bio16", "bio18")

var_order <- c("bio1", "bio3", "bio5", "bio6", "bio15", "bio17", "bio18", "bio19", 
               "slope", "aspect", "flow_acc", "solar")

## ---(1) Define functions ----
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


## ----(2) Read models and extract variable importance-----
rds_dir <- here("results", "models","Rerun","K_fold","bird")
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
  
  out_file <- file.path(rds_dir, "OtherAnimal_species_varimp_combined.xlsx")
  write_xlsx(list(all_varimp = combined_df), out_file)
  message("Saved to：", out_file)
} else {
  message("No more data to be written")
}

## -----(4) Read models and extract AUC ------
rds_dir <- here("results", "models","Rerun","K_fold","bird")
rds_files <- list.files(rds_dir, pattern = "\\.RDS$", ignore.case = TRUE, full.names = FALSE)

all_auc <- list()

for (rds_file in rds_files) {
  rds_path <- file.path(rds_dir, rds_file)   
  Targetmodel <- readRDS(rds_path)
  species <- file_path_sans_ext(basename(rds_file))
  
  if (!is.null(Targetmodel[["kfold_auc"]])) {
    auc_df <- Targetmodel[["kfold_auc"]] %>%
      summarise(
        gam_auc_mean = mean(GAM, na.rm = TRUE),
        rf_auc_mean = mean(RF, na.rm = TRUE),
        maxent_auc_mean = mean(MaxEnt, na.rm = TRUE)
      ) %>%
      mutate(species = species) %>%
      select(species, gam_auc_mean, rf_auc_mean, maxent_auc_mean)
  } else {
    auc_df <- tibble(
      species = species,
      gam_auc_mean = NA_real_,
      rf_auc_mean = NA_real_,
      maxent_auc_mean = NA_real_
    )
  }
  
  all_auc[[length(all_auc) + 1]] <- auc_df
}

if (length(all_auc) > 0) {
  combined_auc <- bind_rows(all_auc)
  
  out_file <- file.path(rds_dir, "bird_species_auc.xlsx")
  write_xlsx(combined_auc, out_file)
  message("Saved AUC results to: ", out_file)
} else {
  message("No AUC results found")
}


# ----3. Plot ----
##  ---- (1) Organize Species scientific name and common name ----
species_lookup <- tribble(
  ~Latin, ~Common,
  
  "Arctostaphylos purissima", "La Purisima manzanita",
  "Cirsium rhothophilum", "Surf thistle",
  "Deinandra increscens ssp. villosa", "Gaviota Tarplant",
  "Horkelia cuneata var. sericea", "Kellogg’s Horkelia",
  "Juglans californica", "Southern California black walnut",
  "Malacothrix saxatilis var. saxatilis", "Cliff aster",
  "Mucronea californica", "California spineflower",
  "Ribes amarum var. hoffmannii", "Hoffman's bitter gooseberry",
  "Eriodictyon capitatum", "Lompoc yerba santa",
  "Astragalus nuttallii var. nuttallii", "Nuttall’s milkvetch",
  "Senecio blochmaniae", "Dune (Blochman's) ragwort",
  "Phacelia hubbyi", "Hubby's phacelia",
  "Abronia maritima", "Red sand verbena",
  "Scrophularia atrata", "Black-flowered figwort",
  
  "Athene cunicularia", "Burrowing Owl",
  "Lanius ludovicianus", "Loggerhead shrike",
  "Setophaga petechia", "Yellow warbler",
  "Elanus leucurus", "White-tailed kites",
  "Pandion haliaetus", "Osprey",
  "Cepphus columba", "Pigeon guillemot",
  "Aquila chrysaetos", "Golden eagle",
  
  "Rana draytonii", "California red-legged frog",
  "Taxidea taxus", "American badger",
  "Thamnophis hammondii", "Two-striped gartersnake",
  "Actinemys marmorata", "Pacific pond turtle",
  "Danaus plexippus", "Monarch butterfly",
  "Puma concolor", "Mountain Lion"
)
##  ---- (2) Read all sheet ----
# Read the file which is the combination of exported variable importance tables for all types of species
# There are three sheets in the file: allSpecies, Bird, plant
sheets <- excel_sheets(here("results","evaluations","KeyHabitat","Varimp_all_species", "NewEnv_species_varimp_combined_shrub_herb.xlsx"))
sheets_to_read <- sheets[1:3]

df_list <- lapply(sheets_to_read, function(sheet) {
  read_excel(here("results","evaluations","KeyHabitat","Shrub" ,"NewEnv_species_varimp_combined_shrub.xlsx"), sheet = sheet) %>%
    mutate(group = sheet)  
})

# Merge all required lists
df_all <- bind_rows(df_list)

# Step 2. species × model normalization
df_scaled <- df_all %>%
  group_by(group, species, model) %>%
  mutate(
    # in case max == min 
    Importance_scaled = if (max(Importance, na.rm = TRUE) == min(Importance, na.rm = TRUE)) {
      0
    } else {
      (Importance - min(Importance, na.rm = TRUE)) /
        (max(Importance, na.rm = TRUE) - min(Importance, na.rm = TRUE))
    }
  ) %>%
  ungroup()

## ---- (3) convert to percentage（total = 100）----
df_pct <- df_scaled %>%
  group_by(group, species, model) %>%
  mutate(
    Importance_pct = if (sum(Importance_scaled, na.rm = TRUE) == 0) {
      0
    } else {
      Importance_scaled / sum(Importance_scaled, na.rm = TRUE) * 100
    }
  ) %>%
  ungroup()

## Test if the last step is correct
#df_pct <- df_pct %>%
#          group_by(group, species, model) %>%
#          summarise(total_pct = sum(Importance_pct, na.rm = TRUE)) %>%
#          ungroup()

## ---- (4) Average the results from the three models ----
df_avg <- df_pct %>%
  group_by(group, species, Variable) %>%
  summarise(Importance_mean = mean(Importance_pct, na.rm = TRUE), .groups = "drop")

## -----(5) Convert to percentage importance -----
## desired variable order

desired_order <- c("bio1", "bio3", "bio5", "bio6", "bio15", "bio17", "bio18", "bio19", 
                   "slope", "aspect", "flow_acc", "solar")

# desired_order <- c("bio14", "bio15", "bio16", "bio18",
#                   "bio1", "bio3", "bio5", "bio6") # Might want to change the variable names 


## organize data
df_avg1 <- df_avg %>%
  mutate(Variable = as.character(Variable)) %>%
  
  mutate(
    Variable = factor(Variable, levels = desired_order, ordered = TRUE)
  ) %>%
  
  group_by(species) %>%
  mutate(
    Importance_pct = Importance_mean / sum(Importance_mean, na.rm = TRUE) * 100
  ) %>%
  ungroup() %>%
  
  mutate(
    group = recode(group,
                   "plant" = "Plants",
                   "Bird" = "Birds",
                   "animal_on_the_ground" = "Other Animals") 
  ) %>%
  mutate(
    species_corrected = recode(
      species,
      "Astragalus nuttallii nuttallii" = "Astragalus nuttallii var. nuttallii",
      "Deinandra increscens villosa" = "Deinandra increscens ssp. villosa",
      "Horkelia cuneata sericea" = "Horkelia cuneata var. sericea",
      "Malacothrix saxatilis saxatilis" = "Malacothrix saxatilis var. saxatilis",
      "Ribes amarum hoffmannii" = "Ribes amarum var. hoffmannii"
      
    )
  ) %>%
  left_join(species_lookup, by = c("species_corrected" = "Latin")) %>%
  mutate(
    species_label = paste0(species_corrected, " (", Common, ")")
  )

## Test if step 5 is correct
df_avg11 <- df_avg1 %>%
  group_by(group, species_corrected) %>%
  summarise(total_pct = sum(Importance_pct, na.rm = TRUE)) %>%
  ungroup()


## ---------(6) IQR plot - Use percentage variable importance -------
summary_df <- df_avg1 %>%
  group_by(group, Variable) %>%
  summarise(
    median = median(Importance_pct, na.rm = TRUE),
    q1 = quantile(Importance_pct, 0.25, na.rm = TRUE),
    q3 = quantile(Importance_pct, 0.75, na.rm = TRUE),
    .groups = "drop"
  )

summary_all <- df_avg1 %>% 
  group_by(Variable) %>% 
  summarise(
    median = median(Importance_pct, na.rm = TRUE),
    q1 = quantile(Importance_pct, 0.25, na.rm = TRUE),
    q3 = quantile(Importance_pct, 0.75, na.rm = TRUE),
    .groups = "drop"
  ) %>% 
  mutate(group = "All species")

summary_df2 <- bind_rows(summary_df, summary_all)

ggplot(summary_df2, aes(x = Variable, y = median, color = group)) +
  geom_errorbar(
    aes(ymin = q1, ymax = q3),
    position = position_dodge(width = 0.6),
    width = 0.2, linewidth = 0.8
  ) +
  geom_point(
    position = position_dodge(width = 0.6),
    size = 2.5
  ) +
  labs(
    title = "Model-averaged variable importance across species groups",
    x = "Variable",
    y = "Importance contribution (%)"
  ) +
  scale_color_manual(
    name = "Species Group",
    values = c(
      "All species" = "#4a4e69",
      "Birds" = "#FFCA3A",
      "Plants" = "#8AC926",
      "Other Animals" = "#1982C4"
    )
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.line = element_line(color = "black", size = 0.6),
    axis.ticks = element_line(color = "black", linewidth = 0.4),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

## --------(7) Percentage variable importance plot for each species (Stacked barplot)-----
p_a <- ggplot(df_avg1, aes(x = species_label, y = Importance_pct, fill = Variable)) +
  geom_col(stat = "identity") +         # stack bar （default position="stack"）
  facet_wrap(~ group, scales = "free_y", ncol = 1) +  # each group will be one raw
  coord_flip() +                        # flip coordinate：species is（y）
  labs(x = NULL, y = "Percentage contribution (%)", fill = "Variable") +
  theme_bw(base_size = 13) +
  theme(
    axis.text.y = element_text(size = 10),   # adjust the size of species name 
    panel.grid.major.x = element_blank(),
    strip.background = element_rect(fill = "grey90", color = NA),
    strip.text = element_text(size = 13, face = "bold")
  ) +
  scale_fill_manual(
    values = c(
      "bio14" = "#440154",  
      "bio15" ="#3B528B",  
      "bio16" = "#21908C",  
      "bio18" = "#5DC863",  
      "bio1" = "#FDE725",  
      "bio3" = "#F8961E",  
      "bio5" = "#D94801", 
      "bio6" = "#8C2D04",
      breaks = desired_order 
    )
  )

