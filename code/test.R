library(SSDM)
library(raster)
df <- read.csv('data/occurrences/animals/birds-cleaned-0515.csv', sep=';')
colnames(df)

options(digits = 16)
df <- df %>%
  mutate(
    coords = gsub("c\\(|\\)", "", geometry),  # remove c() wrapper
    x = as.double(sapply(strsplit(coords, ","), `[`, 1)),  # first coord
    y = as.double(sapply(strsplit(coords, ","), `[`, 2))   # second coord
  ) 
write.table(df, 'data/occurrences/animals/birds_cleaned_0515_xy.csv', sep=';')


Env <- load_var(path = 'data/env', format = '.tif', verbose = FALSE)
Env

Occ <- load_occ(path = 'data/occurrences/animals', Env,
                Xcol = 'x', Ycol = 'y',
                file = 'birds_cleaned_0515_xy.csv', sep = ';', verbose = FALSE)
head(Occ)