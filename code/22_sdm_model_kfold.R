## Purpose of script: K-fold evaluation for species distribution model
## Authors: GEOG 274
## Date: Summer, 2025
## Credits to: Yanni Zhan, Xue Yan, Wenxin Yang, Yifei Liu

library(mgcv)
library(dismo)
library(randomForest)
library(raster)
library(pROC)
library(dplyr)
library(caret)
library(doParallel)
library(tidyr)
library(rsample)

rm(list=ls())
#setwd("C:/Users/yzhan/Desktop/courses/geog_274")

setwd("C:/Users/xyan253/OneDrive/UCSB/Class/GEOG274/preserviaPulse-V2")
setwd("E:/OneDrive/UCSB/Class/GEOG274/preserviaPulse-V2")

# --- Inputs ---
##### Occ data
occ_pre <- read.csv('./data/occurrences/Anim_Plant_merge_final0711.csv')
occ_pre <- na.omit(occ_pre)
occ_pre$Presence <- 1


##### env_stack: RasterStack of predictors
Env <- stack("./data/env/final_env_1980_2010_stack.tif")
#Env <- Env[[c(1:3,5:12,19)]] 

#Env <- Env[[c(6,10,11,12,17,7,18,9)]] 

names(Env)
plot(Env)

# z score
means <- cellStats(Env, "mean", na.rm=TRUE)
sds   <- cellStats(Env, "sd",   na.rm=TRUE)

env_stack <- brick(lapply(seq_len(nlayers(Env)), function(i) {
  (Env[[i]] - means[i]) / sds[i]
}))

names(env_stack)
plot(env_stack)


##### species list
# plant list
taxon <- "plant"

# shrub
species_lst <- c("Arctostaphylos purissima", 
                 "Cirsium rhothophilum",
                 "Deinandra increscens villosa",
                 "Horkelia cuneata sericea",
                 "Juglans californica", 
                 "Malacothrix saxatilis saxatilis",
                 "Mucronea californica", 
                 "Ribes amarum hoffmannii",
                 "Eriodictyon capitatum",
                 "Astragalus nuttallii nuttallii",
                 "Senecio blochmaniae")

# herb
species_lst <- c(#"Arctostaphylos purissima",
  #"Astragalus nuttallii nuttallii",
  #"Deinandra increscens villosa",
  #"Horkelia cuneata sericea",
  "Phacelia hubbyi",
  "Abronia maritima", 
  #"Cirsium rhothophilum",
  "Scrophularia atrata" 
  #"Mucronea californica"
)

# all
species_lst <- c("Arctostaphylos purissima", 
                 "Cirsium rhothophilum",
                 "Deinandra increscens villosa",
                 "Horkelia cuneata sericea",
                 "Juglans californica", 
                 "Malacothrix saxatilis saxatilis",
                 "Mucronea californica", 
                 "Ribes amarum hoffmannii",
                 "Eriodictyon capitatum",
                 "Astragalus nuttallii nuttallii",
                 "Senecio blochmaniae",
                 "Phacelia hubbyi",
                 "Abronia maritima", 
                 "Scrophularia atrata")

taxon <- "bird"

# bird in shrub
species_lst <- c("Athene cunicularia",
                 "Lanius ludovicianus",
                 "Setophaga petechia",
                 "Elanus leucurus")

# bird in herb
species_lst <- c("Athene cunicularia",
                 "Elanus leucurus",
                 "Pandion haliaetus",
                 "Cepphus columba",
                 "Lanius ludovicianus", 
                 "Aquila chrysaetos")
# all
species_lst <- c("Athene cunicularia",
                 "Elanus leucurus",
                 "Pandion haliaetus",
                 "Cepphus columba",
                 "Lanius ludovicianus", 
                 "Aquila chrysaetos",
                 "Setophaga petechia")

taxon <- "OtherAnimal"
# other animals in shrub
species_lst <- c("Rana draytonii",
                 "Taxidea taxus",
                 "Thamnophis hammondii",
                 "Actinemys marmorata", # herps
                 "Danaus plexippus") # invert

# other animals in herb
species_lst <- c("Puma concolor",
                 "Rana draytonii",
                 "Taxidea taxus",
                 "Thamnophis hammondii",
                 "Actinemys marmorata", # herps
                 "Danaus plexippus") 


# check stats
species_lst %in% unique(occ_pre$species)

occ_pre %>% filter(species %in% species_lst) %>% group_by(species) %>% summarize(obs=n())


##### paralel
#cl = makeCluster(10)
#registerDoParallel(cl)

# --- Prediction storage for each fold ---
gam_preds <- list()
rf_preds <- list()
maxent_preds <- list()

#foreach( j=1:length(species_lst),.packages=c("raster","mgcv","dismo","randomForest","pROC","dplyr","caret","rsample") ) %dopar% {
for(j in 1:length(species_lst)){
  
  # extract target species
  #target_sp <- "Juglans californica"
  target_sp <- species_lst[j]
  print(target_sp)
  
  occ_pre_target <- subset(occ_pre, occ_pre$species == target_sp)
  head(occ_pre_target)
  
  
  ##### Abs data
  occ_abs <- read.csv("./data/sampled_background_points.csv")
  occ_abs <- na.omit(occ_abs)
  occ_abs <- occ_abs[,2:3]
  occ_abs$species <- target_sp
  occ_abs$Presence <- 0
  
  
  
  ##### Set k fold
  #set.seed(274)
  k <- 10
  #occ_pre_target$fold <- sample(1:k, nrow(occ_pre_target), replace = TRUE)
  #occ_abs$fold <- sample(1:k, nrow(occ_abs), replace = TRUE)
  
  # merge
  #occ_target <- rbind(occ_pre_target[,-2], occ_abs)
  #occ_data <- occ_target
  
  # Extract environmental variables
  env_occ <- raster::extract(env_stack, occ_pre_target[, c("x", "y")])
  env_occ <- na.omit(cbind(occ_pre_target[,-2], env_occ))
  
  set.seed(274)
  cv_splits_occ <- vfold_cv(env_occ, v = k)
  
  
  env_abs <- raster::extract(env_stack, occ_abs[, c("x", "y")])
  env_abs <- na.omit(cbind(occ_abs, env_abs))
  
  set.seed(274)
  cv_splits_abs <- vfold_cv(env_abs, v = k)
  
  
  
  # ---  GAM and RF formula ---
  var_names <- names(env_stack)
  gam_formula <- as.formula(paste("Presence ~", paste0("s(", var_names, ")", collapse = " + ")))
  rf_formula  <- as.formula(paste("as.factor(Presence) ~", paste(var_names, collapse = " + ")))
  
  
  
  # --- Storage ---
  auc_list <- data.frame(Fold = 1:k, GAM = NA, RF = NA, MaxEnt = NA)
  thresh_list <- data.frame(Fold = 1:k, GAM = NA, RF = NA, MaxEnt = NA)
  sensitivities <- data.frame(Fold = 1:k, GAM = NA, RF = NA, MaxEnt = NA)
  specificities <- data.frame(Fold = 1:k, GAM = NA, RF = NA, MaxEnt = NA)
  tss_list <- data.frame(Fold = 1:k, GAM = NA, RF = NA, MaxEnt = NA)
  kfold <- data.frame(Fold = 1:k, train_0 = NA, train_1 = NA, test_0 = NA, test_1 = NA)
  
  
  
  # --- k-fold loop ---
  # Note: Defining a helper function for robustness outside the main loop is best practice, 
  # but for simplicity, we'll keep the logic inline for MaxEnt and use the unified approach for all three.
  
  for (i in 1:k) {
    
    cat("Fold", i, "\n")
    
    # Access a single split:
    #split <- cv_splits$splits[[i]]
    #train_data <- analysis(split)
    #test_data <- assessment(split)
    
    
    # Extract fold i for presence and absence
    pres_split <- cv_splits_occ$splits[[i]]
    abs_split  <- cv_splits_abs$splits[[i]]
    
    # Training data
    train_pres <- analysis(pres_split)
    train_abs  <- analysis(abs_split)
    train_data <- bind_rows(train_pres, train_abs)
    
    # Testing data
    test_pres <- assessment(pres_split)
    test_abs  <- assessment(abs_split)
    test_data <- bind_rows(test_pres, test_abs)
    
    
    kfold$train_0[i] <- sum(train_data$Presence==0)
    kfold$train_1[i] <- sum(train_data$Presence==1)
    kfold$test_0[i] <- sum(test_data$Presence==0)
    kfold$test_1[i] <- sum(test_data$Presence==1)
    
    
    #train <- occ_data %>% filter(fold != i)
    #test  <- occ_data %>% filter(fold == i)
    
    # Extract environmental variables
    #env_train <- extract(env_stack, train[, c("x", "y")])
    #env_test  <- extract(env_stack, test[, c("x", "y")])
    
    #train_data <- na.omit(cbind(train, env_train))
    #test_data  <- na.omit(cbind(test, env_test))
    
    # --- Helper function for robust coord extraction ---
    # Function to get the optimal threshold coordinates and handle matrix output
    get_best_coords <- function(roc_obj) {
      best_coords <- coords(roc_obj, "best", 
                            ret = c("threshold", "specificity", "sensitivity"), 
                            transpose = FALSE)
      # Robustness check: if multiple optimal thresholds exist, coords returns a matrix.
      if (is.matrix(best_coords)) {
        # Take the results from the first optimal threshold found
        best_coords <- best_coords[1, ]
      }
      return(best_coords)
    }
    
    # --- GAM ---
    print("GAM")
    # train
    prNum <- sum(train_data$Presence == 1) # number of presence records
    bgNum <- sum(train_data$Presence == 0)
    wt <- ifelse(train_data$Presence == 1, 1, prNum / bgNum)
    
    gam_model <- gam(gam_formula, data = train_data, 
                     family = binomial(link = "logit"), weights = wt, method = "REML")
    gam_pred <- predict(env_stack, gam_model, type = "response")
    gam_preds[[i]] <- gam_pred
    
    # test
    pred_test1 <- predict(gam_model, newdata = test_data, type = "response")
    roc_obj1 <- roc(test_data$Presence, pred_test1)
    auc_list$GAM[i] <- auc(roc_obj1)
    
    # Robustly extract optimal threshold and metrics
    best_coords1 <- get_best_coords(roc_obj1)
    thresh_list$GAM[i] <- as.numeric(best_coords1["threshold"])
    sensitivities$GAM[i] <- as.numeric(best_coords1["sensitivity"])
    specificities$GAM[i] <- as.numeric(best_coords1["specificity"])
    tss_list$GAM[i] <- sensitivities$GAM[i] + specificities$GAM[i] - 1 # Calculate TSS
    
    #thresh_list$GAM[i] <- coords(roc_obj1, "best")[1]
    #sensitivities$GAM[i] <- as.numeric(coords(roc_obj1, "best")[3])
    #specificities$GAM[i] <- as.numeric(coords(roc_obj1, "best")[2])
    #tss_list$GAM[i] <- sensitivities$GAM[i] + specificities$GAM[i] - 1
    
    
    # --- RF ---
    print("RF")
    # train
    spsize <- c("0" = prNum, "1" = prNum)
    
    set.seed(274)
    rf_model <- randomForest(rf_formula, data = train_data, importance=TRUE, ntree = 1000, sampsize = spsize, replace = TRUE)
    rf_pred <- predict(env_stack, rf_model, type = "prob", index = 2)
    rf_preds[[i]] <- rf_pred
    
    # test
    pred_test2 <- predict(rf_model, newdata = test_data, type = "prob")[,2]
    roc_obj2 <- roc(test_data$Presence, pred_test2)
    auc_list$RF[i] <- auc(roc_obj2)
    
    # Robustly extract optimal threshold and metrics
    best_coords2 <- get_best_coords(roc_obj2)
    
    thresh_list$RF[i] <- as.numeric(best_coords2["threshold"])
    sensitivities$RF[i] <- as.numeric(best_coords2["sensitivity"])
    specificities$RF[i] <- as.numeric(best_coords2["specificity"])
    tss_list$RF[i] <- sensitivities$RF[i] + specificities$RF[i] - 1 # Calculate TSS
    
    #thresh_list$RF[i] <- coords(roc_obj2, "best")[1]
    #sensitivities$RF[i] <- as.numeric(coords(roc_obj2, "best")[3])
    #specificities$RF[i] <- as.numeric(coords(roc_obj2, "best")[2])
    #tss_list$RF[i] <- sensitivities$RF[i] + specificities$RF[i] - 1
    
    
    
    # --- MaxEnt ---
    print("MaxEnt")
    # train
    maxent_model <- maxent(x=train_data[,-(1:4)], p=train_data$Presence)
    maxent_pred <- predict(env_stack, maxent_model)
    maxent_preds[[i]] <- maxent_pred
    
    # test
    pred_test3 <- predict(maxent_model, test_data)
    roc_obj3 <- roc(test_data$Presence, pred_test3)
    auc_list$MaxEnt[i] <- auc(roc_obj3)
    #thresh_list$MaxEnt[i] <- coords(roc_obj3, "best")[1]
    #sensitivities$MaxEnt[i] <- as.numeric(coords(roc_obj3, "best")[3])
    #specificities$MaxEnt[i] <- as.numeric(coords(roc_obj3, "best")[2])
    #tss_list$MaxEnt[i] <- sensitivities$MaxEnt[i] + specificities$MaxEnt[i] - 1
    
    # Robustly extract optimal threshold and metrics (TSS Max)
    best_coords3 <- get_best_coords(roc_obj3)
    thresh_list$MaxEnt[i] <- as.numeric(best_coords3["threshold"])
    sensitivities$MaxEnt[i] <- as.numeric(best_coords3["sensitivity"])
    specificities$MaxEnt[i] <- as.numeric(best_coords3["specificity"])
    tss_list$MaxEnt[i] <- sensitivities$MaxEnt[i] + specificities$MaxEnt[i] - 1
    
  }
  
  #auc_list
  #thresh_list
  #sensitivities
  #specificities
  #tss list
  
  
  # Organize AUC and TSS data
  auc_long <- auc_list %>% 
    tidyr::pivot_longer(cols = -Fold, names_to = "Model", values_to = "AUC")
  tss_long <- tss_list %>% 
    tidyr::pivot_longer(cols = -Fold, names_to = "Model", values_to = "TSS")
  summary_metrics <- bind_rows(
    auc_long %>% group_by(Model) %>% summarise(Metric = "AUC", Mean = mean(AUC), SD = sd(AUC), .groups = 'drop'),
    tss_long %>% group_by(Model) %>% summarise(Metric = "TSS", Mean = mean(TSS), SD = sd(TSS), .groups = 'drop')
  )
  
  # Model uncertainty based on k-fold
  gam_stack <- stack(gam_preds)
  gam_mean <- mean(gam_stack)
  gam_sd <- calc(gam_stack, fun=sd)
  
  rf_stack <- stack(rf_preds)
  rf_mean <- mean(rf_stack)
  rf_sd <- calc(rf_stack, fun=sd)
  
  maxent_stack <- stack(maxent_preds)
  maxent_mean <- mean(maxent_stack)
  maxent_sd <- calc(maxent_stack, fun=sd)
  
  ensemble_stack_means <- stack(gam_mean, rf_mean, maxent_mean)
  ensemble_mean <- mean(ensemble_stack_means)
  ensemble_sd_model_diff <- calc(ensemble_stack_means, fun=sd)
  
  
  output_dir <- paste0("results/evaluations/K_fold_auc_tss/prediction_maps/", taxon, "/")
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  writeRaster(ensemble_mean, 
              filename = paste0(output_dir, target_sp, "_ensemble_mean.tif"), overwrite=TRUE)
  writeRaster(ensemble_sd_model_diff, 
              filename = paste0(output_dir, target_sp, "_ensemble_sd_model_diff.tif"), overwrite=TRUE)
  
  writeRaster(gam_sd, filename = paste0(output_dir, target_sp, "_GAM_sd.tif"), overwrite=TRUE)
  
  writeRaster(rf_sd, filename = paste0(output_dir, target_sp, "_RF_sd.tif"), overwrite=TRUE)
  
  writeRaster(maxent_sd, filename = paste0(output_dir, target_sp, "_MaxEnt_sd.tif"), overwrite=TRUE)
  
  
  
  # --- Output ---
  model_list <- list(
    
    kfold = kfold,
    kfold_auc = auc_list,
    kfold_thresh = thresh_list,
    kfold_sen = sensitivities,
    kfold_spec = specificities,
    kfold_tss = tss_list,
    summary_metrics = summary_metrics
    
  )
  
  #saveRDS(model_list, file = paste0("results/models/Rerun/K_fold/",taxon,"/",target_sp, "_kfold.RDS"))
  saveRDS(model_list, file = paste0("results/evaluations/K_fold_auc_tss/",taxon,"/",target_sp, "_kfold.RDS"))
  #rds_dir <- paste0("results/evaluations/K_fold_auc_tss_data/", taxon, "/")
  #dir.create(rds_dir, showWarnings = FALSE, recursive = TRUE)
  #saveRDS(model_list, file = paste0(rds_dir, target_sp, "_kfold.RDS"))
  
}








##### check results
species_lst[1]
test=readRDS(paste0("evaluation/",taxon,"/",species_lst[1], "_kfold.RDS"))
test0=readRDS(paste0("new/",taxon,"/",species_lst[1], ".RDS"))
#test$gam
#test$maxent
#test$rf

test$kfold
test$kfold_auc
test$kfold_thresh
test$kfold_sen
test$kfold_spec

test$gam_thresh
test$rf_thresh
test$maxent_thresh

plot(test$rf_prob)
plot(test$rf_binary)
plot(test$maxent_prob)
plot(test$maxent_binary)
plot(test$gam_prob)
plot(test$gam_binary)




##### create csv
kfold_files=list.files(path="./evaluation",pattern="*.kfold.RDS",recursive=T,full.names=T)
kfold_files

library(dplyr)
library(tidyr)
library(purrr)

results_list <- list()

for (file in kfold_files) {
  rds_data <- readRDS(file)
  
  # Extract species name (remove "_kfold*.RDS")
  species_name <- gsub("_kfold.*\\.RDS$", "", basename(file))
  
  # Extract taxon from folder name (e.g., "plant", "bird")
  taxon <- basename(dirname(file))
  
  # Pivot each metric to long format
  auc_long <- rds_data$kfold_auc %>%
    pivot_longer(cols = -Fold, names_to = "Model", values_to = "AUC")
  
  thresh_long <- rds_data$kfold_thresh %>%
    pivot_longer(cols = -Fold, names_to = "Model", values_to = "Threshold")
  
  sen_long <- rds_data$kfold_sen %>%
    pivot_longer(cols = -Fold, names_to = "Model", values_to = "Sensitivity")
  
  spec_long <- rds_data$kfold_spec %>%
    pivot_longer(cols = -Fold, names_to = "Model", values_to = "Specificity")
  
  # Combine all metrics into one data frame
  combined <- reduce(list(auc_long, thresh_long, sen_long, spec_long), full_join, by = c("Fold", "Model"))
  
  # Add Species and Taxon columns
  combined <- combined %>%
    mutate(Species = species_name, Taxon = taxon)
  
  results_list[[length(results_list) + 1]] <- combined
}

# Combine all into one final data frame
final_df <- bind_rows(results_list) %>%
  dplyr::select(Taxon, Species, Fold, Model, AUC, Threshold, Sensitivity, Specificity)

# check
final_df$Threshold[sapply(final_df$Threshold, function(x) length(x) > 1)]
final_df$Specificity[sapply(final_df$Specificity, function(x) length(x) > 1)]

final_df[sapply(final_df$Specificity, function(x) length(x) > 1), ]

# View or write to file
print(head(final_df))


#####
# Define a function to check and collapse
collapse_if_multiple <- function(x) {
  sapply(x, function(el) {
    if (is.null(el) || length(el) == 0) return(NA)
    if (is.data.frame(el) || is.matrix(el)) el <- as.vector(as.matrix(el))
    if (length(el) > 1) paste(as.numeric(el), collapse = ",") else as.numeric(el)
  })
}

final_df$AUC <- collapse_if_multiple(final_df$AUC)
final_df$Threshold <- collapse_if_multiple(final_df$Threshold)
final_df$Sensitivity <- collapse_if_multiple(final_df$Sensitivity)
final_df$Specificity <- collapse_if_multiple(final_df$Specificity)

write.csv(final_df, "./evaluation/combined_kfold_results.csv", row.names = FALSE)

#####
evaluation_csv=read.csv("./evaluation/combined_kfold_results2.csv")
colnames(evaluation_csv)

evaluation_csv <- evaluation_csv %>%
  mutate(TSS = as.numeric(Sensitivity) + as.numeric(Specificity) - 1)

write.csv(evaluation_csv, "./evaluation/combined_kfold_results2.csv", row.names = FALSE)


avg_results <- evaluation_csv %>%
  group_by(Taxon, Species, Model) %>%
  summarise(
    AUC = mean(as.numeric(AUC), na.rm = TRUE),
    Threshold = mean(as.numeric(Threshold), na.rm = TRUE),
    Sensitivity = mean(as.numeric(Sensitivity), na.rm = TRUE),
    Specificity = mean(as.numeric(Specificity), na.rm = TRUE),
    TSS = mean(TSS, na.rm = TRUE),
    .groups = "drop"
  )

avg_results2 <- evaluation_csv %>%
  group_by(Taxon, Species) %>%
  summarise(
    AUC = mean(as.numeric(AUC), na.rm = TRUE),
    Threshold = mean(as.numeric(Threshold), na.rm = TRUE),
    Sensitivity = mean(as.numeric(Sensitivity), na.rm = TRUE),
    Specificity = mean(as.numeric(Specificity), na.rm = TRUE),
    TSS = mean(TSS, na.rm = TRUE),
    .groups = "drop"
  )

library(openxlsx)

wb <- createWorkbook()
addWorksheet(wb, "GAM MaxEnt RF")
addWorksheet(wb, "Ens Model")

writeData(wb, sheet = "GAM MaxEnt RF", avg_results)
writeData(wb, sheet = "Ens Model", avg_results2)

saveWorkbook(wb, "./evaluation/model_summary.xlsx", overwrite = TRUE)



# SSDM
#
j=1
species_lst[j]
a=readRDS(paste0("results_plants/models/",species_lst[j],".RDS"))
plot(a@projection)
plot(a@uncertainty) #different colors

all(values(a@sdms[[2]]@projection)==values(maxent_pred_all))
identical(a@sdms[[2]]@projection,maxent_pred_all)
plot(a@sdms[[2]]@projection)

plot(a@sdms[[1]]@projection)

