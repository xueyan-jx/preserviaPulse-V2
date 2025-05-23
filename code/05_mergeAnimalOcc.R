## Purpose of script: merge multi-source data for animal species
## Authors: GEOG 274
## Date: Spring, 2025

library(here)
library(sf)
library(terra) # read raster
library(dplyr)
library(rgbif)
library(googledrive)
library(googlesheets4)
library(tigris) # to download county boundary
library(stringr)
library(purrr)
library(tidyr)
library(ggplot2)
library(CoordinateCleaner) # add spatial filters
library(readr) # for using write_csv


# Set a directory for data
here() # first check path
anim_dir <- here("data/occurrences/animals") # create a data folder with relative path
file.path(anim_dir) # double check its absolute path
if (!dir.exists(anim_dir)) dir.create(anim_dir) # if it's not there already, create it

# Authorization to Google Drive
drive_auth()
gs4_auth(token = drive_token())

# Get boundaries for areas of interest
ca_counties <- counties(state = "CA", cb = TRUE) # gather CA county boundary file
target_counties <- ca_counties %>%
  filter(NAME %in% c("Santa Barbara", "Ventura", "San Luis Obispo")) # limit to counties of interest
bbox_tgt <- st_bbox(target_counties) # also get the bounding box to add a spatial filter when querying
bbox_wkt <- st_as_text(st_as_sfc(bbox_tgt)) # re-format to use it in rgbif

# ------------------- get a list of all species --------------------
ss <- drive_get("Special Status Species")
all_names <- as.data.frame(matrix(ncol=3, nrow=0))
colnames(all_names) <- c('Name (Latin)', 'Name (common)', 'taxon')
for(taxon in c('Plants', 'Birds', 'Mammals', 'Herps', 'Inverterbrates')){
  dat <- read_sheet(ss, sheet=taxon) %>% dplyr::select(`Name (Latin)`, `Name (common)`)
  dat$taxon <- taxon
  all_names <- rbind(all_names, dat)
}

all_names <- all_names %>% mutate(
  `Name (Latin)` = ifelse(`Name (Latin)`=='Vulpes vulpes ssp.', 'Vulpes vulpes', `Name (Latin)`)
) %>% unique()

# -------------------- download GBIF data from Google Drive ----------------
all_cleaned_files <- drive_find(pattern = "-cleaned.csv", type = "csv")
all_cleaned_files$name
# I then mannually matched with info on gdrive

drive_download(all_cleaned_files[1,]$id, path = here('data/occurrences/animals/birds.csv'))
drive_download(all_cleaned_files[3,]$id, path = here('data/occurrences/animals/herps_invs.csv'))
drive_download(all_cleaned_files[2,]$id, path = here('data/occurrences/animals/mammals.csv'))


## ------------- (1) check problematic birds - skip ---------------
# read all csv starts with bird and rbind
bird_files <- list.files(path='data/occurrences/animals', pattern = "birds*")
gbif_birds <- as.data.frame(matrix(nrow=0, ncol=4))
colnames(gbif_birds) <- c('species', 'decimalLongitude', 'decimalLatitude', 'source')
for(i in 1:length(bird_files)){
  tmp <- read.csv(file.path('data/occurrences/animals/', bird_files[i]), sep=';')
  tmp <- tmp %>% select(species, decimalLongitude, decimalLatitude)
  tmp$source <- paste0('birds', as.character(i))
  gbif_birds <- rbind(gbif_birds, tmp)
}
rm(tmp)
nrow(gbif_birds) # 632019
gbif_birds_unique <- unique(gbif_birds) # 98712

bird_info <- as.data.frame(table(gbif_birds_unique[c('species', 'source')])) %>% 
  filter(Freq>0)
dup_birds <- as.data.frame(table(bird_info$species)) %>% filter(Freq>1)
dup_birds <- as.character(dup_birds$Var1)

# 4 diff!
# see for each bird ...
problem_bird <- gbif_birds_unique %>% filter(species == dup_birds[6])
unq_problem_bird <- unique(problem_bird)
table(unq_problem_bird$source)
nrow(unq_problem_bird)

unq_bird <- unique(unq_problem_bird %>% select(decimalLatitude, decimalLongitude))
nrow(unq_bird)

gbif_birds_unique_final <- unique(gbif_birds[c('species', 'decimalLongitude', 'decimalLatitude')]) # 86369
nrow(gbif_birds_unique_final)

## ------------- (2) further clean birds ---------------
gbif_birds <- read.csv('data/occurrences/animals/birds.csv', sep=';') %>%
  dplyr::select(species, decimalLongitude, decimalLatitude) %>% unique()
birds_name <- unique(gbif_birds$species)
write_csv(gbif_birds, 'data/occurrences/animals/gbif_birds_unique_final.csv')

## ------------- (3) further clean mammals ---------------
gbif_mammals <- read.csv('data/occurrences/animals/mammals.csv', sep=';') %>%
  dplyr::select(species, decimalLongitude, decimalLatitude) %>% unique()
# mammal_final_info <- as.data.frame(table(gbif_mammals$species))
mammal_name <- unique(gbif_mammals$species)
write_csv(gbif_mammals, 'data/occurrences/animals/gbif_mammals_unique_final.csv')

## ------------- (4) further clean herps & invs ---------------
gbif_herps_invs <- read.csv('data/occurrences/animals/herps_invs.csv', sep=';') %>%
  dplyr::select(species, decimalLongitude, decimalLatitude) %>% unique()
# herpamp_final_info <- as.data.frame(table(gbif_herps_invs$species))
herps_invs_name <- unique(gbif_herps_invs$species)
write_csv(gbif_herps_invs, 'data/occurrences/animals/gbif_herps_invs_unique_final.csv')


# ---------------- download shapefile from Wildlife Observations ... ----------------
# this is probably easier to download manually so I did that
wildlife_obs <- st_read('data/occurrences/wildlifeobs/jldp_wildlife_observations_2012_2014.shp') 
unique(wildlife_obs$Grouping)
rpj_wildlife_obs <- st_transform(wildlife_obs, crs=4326) %>% 
  filter(Grouping %in% c('Reptiles', 'Amphibian', 'Mammal', 'Bat', 'Bird')) %>% 
  dplyr::select(Species, Grouping, geometry)

# extract coordinates
wildlife_coords <- st_coordinates(rpj_wildlife_obs)
rpj_wildlife_obs$decimalLongitude <- wildlife_coords[,1]
rpj_wildlife_obs$decimalLatitude <- wildlife_coords[,2]
rpj_wildlife_obs$geometry <- NULL

rpj_wildlife_obs_unique <- unique(rpj_wildlife_obs) %>% 
  mutate(Species = ifelse(Species == 'California red legged frog', 'California red-legged frog', Species))
rpj_wildlife_obs_unique <- rpj_wildlife_obs_unique %>% 
  mutate(Species = ifelse(Species=='Two-striped Gartersnake', 'Two-striped gartersnake', Species))
rpj_wildlife_obs_unique <- rpj_wildlife_obs_unique %>% 
  mutate(Species = ifelse(Species=='American Badger', 'American badger', Species))


rpj_wildlife_obs_unique <- merge(rpj_wildlife_obs_unique, all_names, by.x='Species', by.y='Name (common)') %>% 
  dplyr::select(`Name (Latin)`, decimalLongitude, decimalLatitude, Grouping)
colnames(rpj_wildlife_obs_unique) <- c('species', 'decimalLongitude', 'decimalLatitude', 'taxon')

unique(rpj_wildlife_obs_unique$taxon)
rpj_wildlife_obs_unique <- rpj_wildlife_obs_unique %>% mutate(
  taxon = ifelse(taxon %in% c('Amphibian', 'Reptiles'), 'herps', ifelse(taxon=='Mammal', 'mammals', 'birds'))
)

# --------------- download dgmd portal data -------------
portal_data_real <- read.csv(here('data/occurrences/animals/DP_clean_animals_0519.csv'), sep=';') %>%
  dplyr::select(species, longitude, latitude)
colnames(portal_data_real) <- colnames(gbif_birds)
head(portal_data_real)

# -------------------- merge a couple data sources ------------------
# birds
bird_names <- all_names %>% filter(taxon == 'Birds')
bird_names <- unique(bird_names$`Name (Latin)`)
wildobs_birds <- rpj_wildlife_obs_unique %>% filter(taxon=='birds') %>% 
  dplyr::select(species, decimalLongitude, decimalLatitude)
dp_birds <- portal_data_real %>% filter(species %in% bird_names)
birds_final <- rbind(gbif_birds, wildobs_birds) %>% unique()
birds_final <- rbind(birds_final, dp_birds) %>% unique()
# birds_final_info <- as.data.frame(table(birds_final$species))
write_csv(birds_final, 'data/occurrences/animals/merged_birds_data_0519.csv')

# mammals
mam_names <- all_names %>% filter(taxon == 'Mammals')
mam_names <- unique(mam_names$`Name (Latin)`)
wildobs_mammals <- rpj_wildlife_obs_unique %>% filter(taxon=='mammals') %>% 
  dplyr::select(species, decimalLongitude, decimalLatitude)
dp_mam <- portal_data_real %>% filter(species %in% mam_names)
mammals_final <- rbind(gbif_mammals, wildobs_mammals) %>% unique()
mammals_final <- rbind(mammals_final, dp_mam) %>% unique()
write_csv(mammals_final, 'data/occurrences/animals/merged_mammals_data_0519.csv')

# herps and invs
herps_invs_names <- all_names %>% filter(taxon %in% c('Herps', 'Invertebrates'))
herps_invs_names <- unique(herps_invs_names$`Name (Latin)`)
wildobs_herps <- rpj_wildlife_obs_unique %>% filter(taxon=='herps') %>% 
  dplyr::select(species, decimalLongitude, decimalLatitude)
dp_herps_invs <- portal_data_real %>% filter(species %in% herps_invs_names)
herps_invs_final <- rbind(gbif_herps_invs, wildobs_herps) %>% unique()
herps_invs_final <- rbind(herps_invs_final, dp_herps_invs) %>% unique()
write_csv(herps_invs_final, 'data/occurrences/animals/merged_herps_invs_data_0519.csv')



