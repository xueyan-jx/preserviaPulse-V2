library(SSDM)
library(raster)
library(here)
df <- read.csv(here('data/occurrences/animals/birds_cleaned_0519.csv'), sep=';')
colnames(df)

options(digits = 16)
df <- df %>%
  mutate(
    coords = gsub("c\\(|\\)", "", geometry),  # remove c() wrapper
    x = as.double(sapply(strsplit(coords, ","), `[`, 1)),  # first coord
    y = as.double(sapply(strsplit(coords, ","), `[`, 2))   # second coord
  ) 
write.table(df, 'data/occurrences/animals/birds_cleaned_0519_xy.csv', sep=';')


Env <- load_var(path = 'data/env', format = '.tif', verbose = FALSE)
Env

Occ <- load_occ(path = 'data/occurrences/animals', Env,
                Xcol = 'x', Ycol = 'y',
                file = 'birds_cleaned_0519_xy.csv', sep = ';', verbose = FALSE)
head(Occ)


SDM <- modelling('GLM', subset(Occ, Occ$species == unique(Occ$species)[1]), 
                 Env, Xcol = 'x', Ycol = 'y', verbose = FALSE)

plot(SDM@projection, main = 'SDM\nwith GLM algorithm')