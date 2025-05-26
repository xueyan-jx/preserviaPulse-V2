library(SSDM)
library(raster)
library(here)

occ_pre <- read.csv(here('data/occurrences/animals/animals-data-ready-occ-0519.csv'))
colnames(occ_pre)
occ_pre <- na.omit(occ_pre)
occ_pre$Presence <- 1

# select species
i <- 1
target_sp <- unique(occ_pre$species)[i]
print(target_sp)

occ_pre_target <- subset(occ_pre, occ_pre$species == target_sp)
head(occ_pre_target)


##### Abs data
occ_abs <- read.csv(here("data/occurrences/sampled_background_points.csv"))
occ_abs <- na.omit(occ_abs)
occ_abs <- occ_abs[,2:3]
occ_abs$species <- target_sp
occ_abs$Presence <- 0

# merge
occ_target <- rbind(occ_pre_target[,-2], occ_abs)


Env <- load_var(path = 'data/env', format = '.tif', verbose = FALSE)
Env


seq_len(length(occ_target$fold))

SDM <- modelling('GLM', occ_target, 
                 Env, Xcol = 'x', Ycol = 'y', Pcol = "Presence",
                 verbose = FALSE,
                 cv = "k-fold",
                 cv.param = c(10,1)
)

SDM
plot(SDM@projection, main = 'SDM\nwith GLM algorithm')



cv.param <- c(10, 1)

if(cv.param[1]<2){
  warning("less than 2 folds selected, automatic adjustment to k=5")
  cv.param[1] <- 5
}


k <- cv.param[1]
rep <- cv.param[2]

for (i in seq_len(length(rep))) {
  data <- occ_target
  data$fold <- 0
  for (p in 0:1) {
    datap <- data[which(data$PRESENCE == p), ]
    indices <- seq_len(length(datap$fold))
    fold <- 1
    while (length(indices) > 0) {
      j <- ifelse(length(indices)==1,indices,sample(indices, 1))
      datap$fold[j] <- fold
      indices <- indices[-which(indices == j)]
      if (fold != k) {
        fold <- fold + 1
      } else {
        fold <- 1
      }
    }
    data[which(data$PRESENCE == p), ] <- datap
  }
  for (j in 1:k) {
    eval.testdata <- data[which(data$fold == j), ]
    eval.testdata <- eval.testdata[-which(names(data) == "fold")]
    eval.traindata <- data[which(data$fold != j), ]
    eval.traindata <- eval.traindata[-which(names(data) == "fold")]
    eval.traindata$train <- TRUE
    # evalobj <- obj
    evalobj@data <- eval.traindata
    model <- get_model(evalobj, ...)
    predicted.values <- c(predict(model, eval.testdata))
    if(!is.null(metric)){
      threshval <- .optim.thresh(eval.testdata$Presence, predicted.values, thresh)
      threshval <- mean(threshval[[which(names(threshval) == metric)]])
      roweval <- dismo::evaluate(p=predicted.values[which(eval.testdata$Presence==1)], a=predicted.values[which(eval.testdata$Presence==0)], threshval)
    } else {
      roweval <- dismo::evaluate(p=predicted.values[which(eval.testdata$Presence==1)], a=predicted.values[which(eval.testdata$Presence==0)])
      threshval <- dismo::threshold(roweval,stat=bin.thresh)
      roweval <- dismo::evaluate(p=predicted.values[which(eval.testdata$Presence==1)], a=predicted.values[which(eval.testdata$Presence==0)], threshval)
    }
    caleval <- sdm::calibration(eval.testdata$Presence,predicted.values, nbin=20, weight=TRUE)
    
    evaldf <- data.frame(threshold=threshval, AUC=roweval@auc, omission.rate=roweval@MCR, sensitivity=roweval@TPR, specificity=roweval@TNR, prop.correct=roweval@CCR, Kappa=roweval@kappa, calibration=caleval@statistic)
    if (i == 1 && j == 1) {
      evaluation <- evaldf
    } else {
      evaluation <- rbind(evaluation, evaldf)
    }
  }
}
} else if (cv == "LOO") {
  data <- obj@data
  predicted.values <- c()
  for (j in seq_len(length(data[, 1]))) {
    eval.testdata <- data[j, ]
    eval.traindata <- data[-j, ]
    eval.traindata$train <- TRUE
    evalobj <- obj
    evalobj@data <- eval.traindata
    model <- get_model(evalobj, ...)
    predicted.values[j] <- c(predict(model, eval.testdata))
  }
  if(!is.null(metric)){
    threshval <- .optim.thresh(data$Presence, predicted.values, thresh)
    threshval <- mean(threshval[[which(names(threshval) == metric)]])
    roweval <- dismo::evaluate(p=predicted.values[which(data$Presence==1)], a=predicted.values[which(data$Presence==0)], threshval)
  } else {
    roweval <- dismo::evaluate(p=predicted.values[which(data$Presence==1)], a=predicted.values[which(data$Presence==0)])
    threshval <- dismo::threshold(roweval,stat=bin.thresh)
    roweval <- dismo::evaluate(p=predicted.values[which(data$Presence==1)], a=predicted.values[which(data$Presence==0)], threshval)
  }
  caleval <- sdm::calibration(data$Presence,predicted.values, nbin=20, weight=TRUE)
  
  evaluation <- data.frame(threshold=threshval, AUC=roweval@auc, omission.rate=roweval@MCR, sensitivity=roweval@TPR, specificity=roweval@TNR, prop.correct=roweval@CCR, Kappa=roweval@kappa, calibration=caleval@statistic)
}
obj@evaluation <- evaluation[1, ]
for (i in seq_len(length(evaluation))) {
  obj@evaluation[i] <- mean(evaluation[, i], na.rm = TRUE)
}

} else {
  # Continuous values of MEMs
  warning("Evaluation is not yet implemented for continuous data of MEMs !")
}

# assign train/test fractions for final model training

if(final.fit.data=='holdout'){
  obj@data <- data
}
if(is.numeric(final.fit.data)){
  if(final.fit.data>0 & final.fit.data<=1){
    data <- obj@data
    data$train <- FALSE
    for (p in 0:1) {
      datap <- data[which(data$Presence == p), ]
      datap$train[sample.int(length(datap$train), round(length(datap$train)*cv.param[3]))] <- TRUE
      data[which(data$Presence == p), ] <- datap
    }
    obj@data <- data
  } else{
    warning("Training fraction needs to be between 0 and 1 for sampling, assuming 1")
    final.fit.data <- 'all'}
}
if(final.fit.data=='all'){
  data$train <- TRUE
  obj@data <- data
}


return(obj)
})

#' @rdname evaluate
#' @export
setMethod("evaluate", "MAXENT.SDM", function(obj, cv, cv.param, final.fit.data='all', bin.thresh = 'SES', metric = NULL, thresh = 1001, Env, ...) {
  # Parameters
  text.cv.param <- character()
  for (i in seq_len(length(cv.param))) {
    text.cv.param <- paste0(text.cv.param, "|", cv.param[i])
  }
  obj@parameters$cv <- cv
  obj@parameters$cv.param <- text.cv.param
  
  if(!is.null(metric)){
    warning("Argument 'metric' is deprecated and will be removed in future versions. Please consider using 'bin.thresh' instead.")
    obj@parameters$metric <- metric
  } else {
    obj@parameters$metric <- bin.thresh
  }
  
  if (all(obj@data$Presence %in% c(0, 1))) {
    # Binary data of SDM model
    
    # translate thresholding terms (metric only for backwards compatibility)
    if(!is.null(metric)){
      metric <- switch(metric, Kappa = "maxKappa", CCR = "max.prop.correct",
                       TSS = "max.sensitivity+specificity", SES = "sensitivity=specificity",
                       LW = "min.occurence.prediction", ROC = "min.ROC.plot.distance")
    } else {
      bin.thresh <- switch(bin.thresh, Kappa = "kappa", CCR = "no_omission",
                           TSS = "spec_sens", SES = "equal_sens_spec",
                           EP = "prevalence")
    }


