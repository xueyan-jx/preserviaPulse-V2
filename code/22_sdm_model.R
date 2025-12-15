## Purpose of script: Species Distribution Model
## Authors: GEOG 274
## Date: Spring, 2025
## Credits to: Wenxin Yang, Yanni Zhan, Xue Yan, Yifei Liu

library(SSDM)
library(dplyr)
library(raster)


rm(list=ls())


# ------------ 1. Get environmental and species list ------------
##### Env data
Env <- stack("final_env_1980_2010_stack.tif")
means <- cellStats(Env, "mean", na.rm=TRUE)
sds   <- cellStats(Env, "sd",   na.rm=TRUE)

# z-score normalization
Env_z <- brick(lapply(seq_len(nlayers(Env)), function(i) {
  (Env[[i]] - means[i]) / sds[i]
}))

names(Env_z) <- names(Env)
plot(Env_z)



##### Occ data
occ_pre <- read.csv('Anim_Plant_merge.csv')
occ_pre <- na.omit(occ_pre)
occ_pre$Presence <- 1

# plant list
species_lst <- c("Juglans californica", 
                 "Abronia maritima", 
                 "Phacelia hubbyi", 
                 "Arctostaphylos purissima", 
                 "Suaeda taxifolia", 
                 "Mucronea californica", 
                 "Scrophularia atrata", 
                 "Astragalus nuttallii nuttallii",
                 "Deinandra increscens villosa",
                 "Senecio blochmaniae", 
                 "Erysimum suffrutescens", 
                 "Calochortus fimbriatus", 
                 "Dichondra occidentalis",
                 "Malacothrix saxatilis saxatilis",
                 "Delphinium umbraculorum",
                 "Lilium humboldtii ocellatum",
                 "Ribes amarum hoffmannii", 
                 "Cirsium rhothophilum",
                 "Eriodictyon capitatum")

species_lst %in% unique(occ_pre$species)


# ------------ 2. Run model for all target species ------------
##### Parallel
cl = makeCluster(19)
registerDoParallel(cl)

foreach( i=1:length(species_lst),.packages=c("raster","SSDM") ) %dopar% {
  #for(i in 1:length(species_lst)){
  
  print(i)
  
  target_sp = species_lst[i]
  occ_pre_target <- subset(occ_pre, occ_pre$species == target_sp)
  head(occ_pre_target)
  
  
  ##### Abs data
  occ_abs <- read.csv("sampled_background_points.csv")
  occ_abs <- na.omit(occ_abs)
  occ_abs <- occ_abs[,2:3]
  occ_abs$species <- target_sp
  occ_abs$Presence <- 0
  
  # merge
  occ_target <- rbind(occ_pre_target[,-2], occ_abs)
  
  
  
  ##### Sampling weight
  prNum <- sum(occ_target$Presence == 1) # number of presence records
  bgNum <- sum(occ_target$Presence == 0)
  wt <- ifelse(occ_target$Presence == 1, 1, prNum / bgNum)
  spsize <- c("0" = bgNum, "1" = prNum)
  
  
  
  ##### SDM building
  set.seed(1)
  
  SDM_ens <- ensemble_modelling(c('GAM','MAXENT','RF'),
                                occ_target, Env_z, Xcol = 'x', Ycol = 'y', Pcol="Presence",
                                cv = "holdout", cv.param = c(0.7, 10),
                                rep = 1, cores = 50, #parmode  = "algorithms",
                                uncertainty = TRUE,
                                SDM.projections = TRUE, # save individual SDM models
                                ensemble.metric = c("AUC"),
                                ensemble.thresh = c(0.75),
                                #axes.metric = "AUC", # access variable importance
                                final.fit.data = "all", # controls which data go into the final fits, "all" uses every record
                                bin.thresh = "SES", # thresholding rule, "SES": sensitivity–specificity equality
                                weight = TRUE, # models are weighted by their performance scores; otherwise each selected model contributes equally
                                gam.args = list(family = binomial(link = "logit"),
                                                weights = wt,
                                                method = "REML"),
                                cta.args = list(ntree = 500,
                                                sampsize = spsize,
                                                replace = TRUE))
  
  if (is.null(SDM_ens)) {
    
    message("No model for ", target_sp, " — skipping save.")
    next
    
  } else {
    
    save.esdm(SDM_ens,name=target_sp, path="results")
    saveRDS(SDM_ens,paste0("results/",target_sp,".RDS"))
    
    
  }
  
  
}

# check results
plot(SDM_ens@projection)
plot(SDM_ens@uncertainty)