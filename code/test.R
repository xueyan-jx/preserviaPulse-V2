library(SSDM)
library(raster)
library(here)
library(dplyr)
library(readr)

df <- read.csv(here('data/occurrences/GBIF_BIEN_DP_final.csv'), sep = ",")
colnames(df)

options(digits = 16)

df <- df %>%
  mutate(
    coords = gsub("c\\(|\\)", "", geometry),  # remove c() wrapper
    x = as.double(sapply(strsplit(coords, ","), `[`, 1)),  # first coord
    y = as.double(sapply(strsplit(coords, ","), `[`, 2))   # second coord
  ) 

df <- df[, c("x", "y", "species", "scientificName")]
readr::write_csv(df, 'data/occurrences/plant_for_model_test.csv')


Env <- load_var(path = 'data/env/', format = '.tif', verbose = FALSE)
Env

Occ <- load_occ(path = 'data/occurrences', Env,
                Xcol = 'x', Ycol = 'y',
                file = 'plant_for_model_test.csv', verbose = FALSE)
head(Occ)


SDM <- modelling('GLM', subset(Occ, Occ$species == unique(Occ$species)[1]), 
                 Env, Xcol = 'x', Ycol = 'y', verbose = FALSE)

plot(SDM@projection, main = 'SDM\nwith GLM algorithm')
