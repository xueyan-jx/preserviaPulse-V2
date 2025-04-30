## Purpose of script: grab GBIF data
## Authors: GEOG 274
## Date: Spring, 2025
## Credits to: Wenxin Yang, Jacqueline Vogel, Yanni Zhan, Xue Yan, Dr. Lei Song (lsong@ucsb.edu)

# Materials:
#   https://docs.ropensci.org/rgbif/index.html
#   More details in: https://techdocs.gbif.org/en/openapi/
#   https://docs.gbif.org/course-data-use/en/data-processing-pipeline.html
#  https://data-blog.gbif.org/post/gbif-filtering-guide/
#. https://data-blog.gbif.org/post/downloading-long-species-lists-on-gbif/

# ------------- Setting up --------------
# Install and load packages
install.packages("rgbif")
install.packages("tigris")
library(here)
library(sf)
library(dplyr)
library(rgbif)
library(googledrive)
library(googlesheets4)
library(tigris) # to download county boundary
library(stringr)
library(purrr)
library(tidyr)
library(ggplot2)



# Set name
myname <- 'Wenxin'

# Set a directory for data
here() # first check path
occ_dir <- here("data/occurrences") # create a data folder with relative path
file.path(occ_dir) # double check its absolute path
if (!dir.exists(occ_dir)) dir.create(occ_dir) # if it's not there already, create it

# Authorization to Google Drive
drive_auth()
gs4_auth(token = drive_token())

# Get boundaries for areas of interest
ca_counties <- counties(state = "CA", cb = TRUE) # gather CA county boundary file
target_counties <- ca_counties %>%
  filter(NAME %in% c("Santa Barbara", "Ventura", "San Luis Obispo")) # limit to counties of interest
bbox_tgt <- st_bbox(target_counties) # also get the bounding box to add a spatial filter when querying
bbox_wkt <- st_as_text(st_as_sfc(bbox_tgt)) # re-format to use it in rgbif

# ------------- Get occurrence info for a species list --------------
## ------------ 1. Get species list ------------

# Get data from the spreadsheet
ss <- drive_get("Special Status Species")

## E.g. sheet names
ss_meta <- gs4_get(ss) 
sheet_names <- ss_meta$sheets$name

dat <- read_sheet(ss, sheet = sheet_names[1]) # Select the taxa you're working on
head(dat)

# for birds, we assigned individuals to different species so we subset the data frame
# skip the following line if you're not working on birds...
# dat <- dat %>% filter(Name=='Wenxin') # replace with your name

# Get species' names
scientific_names_df <- dat[, c("Name (common)", "Name (Latin)")]

# Use Latin name for searching
scientific_names <- scientific_names_df$`Name (Latin)`

## ------------ 2. Get information from gbif ------------
# Confirm names used in GBIF database for your species list
gbif_names <- name_backbone_checklist(scientific_names)

# Only keep valid records
gbif_names <- gbif_names[!is.na(gbif_names$usageKey), ]

# ------------- Query and save data --------------
## Add a few internal filters for something that we definitely not want.
## we can do other filters afterwards.
## 1. Zero coordinate : Coordinates are exactly (0,0). null island
## 2. Country coordinate mismatch : The coordinates fall outside of the 
##    given country's polygon.
## 3. Coordinate invalid : GBIF is unable to interpret the coordinates.
## 4. Coordinate out of range : The coordinates are outside of the range 
##    for decimal lat/lon values ((-90,90), (-180,180)).
## 5. Only keep PRESENT records
## 6. Remove FOSSIL_SPECIMEN and LIVING_SPECIMEN

#remove NA values from usage key before query

# Spin up a download request for GBIF occurrence data
# Reference
# https://docs.ropensci.org/rgbif/articles/getting_occurrence_data.html
# https://data-blog.gbif.org/post/downloading-long-species-lists-on-gbif/

# first make sure if your GBIF username and password are stored in REnviron
# if not please add the following lines to the opened file.
# GBIF_USER=your username
# GBIF_PWD=your password
# GBIF_EMAIL=your email
# Reference
# https://docs.ropensci.org/rgbif/articles/gbif_credentials.html
usethis::edit_r_environ()

# apply a spatial filter to limit 

# then download with the species list
occ_meta <- occ_download(
  pred_in("taxonKey", unique(gbif_names$usageKey)),
  pred("occurrenceStatus", "PRESENT"), 
  pred("hasCoordinate", TRUE), pred("hasGeospatialIssue", FALSE),
  pred_not(pred_in("basisOfRecord",
                   c("FOSSIL_SPECIMEN","LIVING_SPECIMEN"))),
  pred("geometry", bbox_wkt),
  format = "SIMPLE_CSV")

# Save out the name catalog using the GBIF query unique ID in name to pair.
write.csv(gbif_names, 
          file.path(occ_dir, sprintf("gbif_names_%s.csv", occ_meta)),
          row.names = FALSE)

# Download the occurrences
status <- occ_download_wait(occ_meta)
if (status$status == "SUCCEEDED"){
  d <- occ_download_get(occ_meta, path = occ_dir, overwrite = TRUE)
} else {
  warning(sprintf("The request %s.", status$status))
  d <- NULL
}

# get zip file name
download_id <- basename(occ_meta)
# if downloaded file is a zipfile, unzip it in the occurrence directory
if(file.exists(file.path(occ_dir, paste0(download_id, '.zip')))){
  # unzip
  unzip(file.path(occ_dir, paste0(download_id, '.zip')), exdir = occ_dir)
}

# ------------- Post-processing --------------
# read in the google spreadsheet to record each step
pt_info <- drive_get("GBIF Occurrence Removal")
missing_species <- read_sheet(pt_info, sheet = "Missing species")
occ_num_info <- read_sheet(pt_info, sheet = "Occurrence number changes")

## ----------- 1. Get species with no records -----------
# read in downloaded data
# note: if you don't have the download id, just copy paste the file name and run the following line
# download_id = '' # your copy pasted file name without the .csv extension
download_gbif <- read.delim(file.path(occ_dir, paste0(download_id, '.csv')))

unique_species <- unique(download_gbif$taxonKey)
length(unique_species)

missing_from_list <- setdiff(gbif_names$usageKey, unique_species)
print(missing_from_list)

# Find missing species
this_info <- gbif_names[gbif_names$usageKey %in% missing_from_list,]
this_info$`Handling person` = myname
print(unique(this_info$canonicalName))

missing_species <- rbind(missing_species, this_info)

# Update the Google spreadsheet
write_sheet(missing_species, pt_info, sheet = "GBIF missing species")

## ----------- 2. Clip data to target counties -----------
# Convert .csv occurrence data to sf object (spatial data)
occ_sf <- st_as_sf(download_gbif, coords = c("decimalLongitude", "decimalLatitude"), crs = 4326)

# Keep same coordinate reference system
occ_sf <- st_transform(occ_sf, st_crs(target_counties))

# Clip to target counties
occ_in_target <- st_intersection(occ_sf, target_counties)

# Check data distribution
plot(st_geometry(target_counties), col = "lightblue", border = "black", main = "GBIF Data Points within Target Counties")
plot(st_geometry(occ_in_target), add = TRUE, col = "red", pch = 20)

coordinates <- st_coordinates(occ_in_target)
occ_in_target_df <- as.data.frame(occ_in_target)
occ_in_target_df$decimalLongitude <- coordinates[, 1] 
occ_in_target_df$decimalLatitude <- coordinates[, 2]  

write.csv(occ_in_target_df, 
          file.path(occ_dir, paste0(download_id, "_3counties_pts.csv")),
          row.names = FALSE)

# Get the target folder first to ensure it exists
speciesObs_folder <- drive_find(pattern = "speciesObs", type = "folder")
if(nrow(speciesObs_folder) > 0) {
  # Upload to Google Drive if folder found
  drive_upload(
    file.path(occ_dir, paste0(download_id, "_3counties_pts.csv")),
    path = as_id(speciesObs_folder$id[1]),
    name = paste0('GBIF_', download_id, "_3counties_pts.csv")
  )
} else {
  warning("Could not find 'specieObs' folder in Google Drive. File saved locally only.")
}

# There is a misalignment of column names in the downloaded file
# Remember to manually change the name for the last two columns 
# (decimalLongitude	& decimalLatitude)

## ----------- 3. Data clean -----------
# We need to exclude records from iNaturalist due to geoprivacy 
# reference: https://help.inaturalist.org/en/support/solutions/articles/151000169938-what-is-geoprivacy-what-does-it-mean-for-an-observation-to-be-obscured-

# Summarize geoprivacy information
## Create bins (cut into intervals of 1000 meters)
occ_gbif_uncertainty_bin <- occ_in_target_df %>%
  filter(!is.na(coordinateUncertaintyInMeters)) %>%
  mutate(uncertainty_bin = cut(
    coordinateUncertaintyInMeters,
    breaks = seq(0, max(coordinateUncertaintyInMeters, na.rm = TRUE) + 1000, by = 1000),
    include.lowest = TRUE,
    right = FALSE
  )) %>%
  group_by(uncertainty_bin) %>%
  summarise(count = n(), .groups = "drop") %>%
  arrange(uncertainty_bin)

print(occ_gbif_uncertainty_bin)

## Remove records with uncertainty >1000 meters
occ_in_target_df <- occ_in_target_df %>%
  filter(coordinateUncertaintyInMeters <=1000)

# Summarize "issue" column (remove any if necessary)
# https://techdocs.gbif.org/en/data-use/occurrence-issues-and-flags

## Separate issues into lists
occ_in_target_df <- occ_in_target_df %>%
  mutate(issue_list = stringr::str_split(issue, ";"))

## Replace any NA in issue_list with an empty character vector
occ_in_target_df$issue_list <- purrr::map(occ_in_target_df$issue_list, 
                                          ~ if (is.null(.x) || anyNA(.x)) character(0) else .x)

## Find all unique issues
all_issues <- occ_in_target_df %>%
  pull(issue_list) %>%
  unlist() %>%
  unique()

all_issues <- all_issues[all_issues != ""] # remove "" (empty character)

## Create a TRUE/FALSE column for each issue
for (iss in all_issues) {
  occ_in_target_df[[iss]] <- purrr::map_lgl(occ_in_target_df$issue_list, ~ iss %in% .x)
}

## Summarize: count how many TRUEs for each issue
issue_counts <- occ_in_target_df %>%
  dplyr::select(all_of(all_issues)) %>%
  summarise(across(everything(), ~sum(.))) %>%
  tidyr::pivot_longer(cols = everything(), names_to = "issue", values_to = "count")

## Plot
ggplot(issue_counts, aes(x = reorder(issue, -count), y = count)) +
  geom_bar(stat = "identity") +
  coord_flip() +  # Flip for easier reading
  labs(
    title = "Count of Issues in GBIF Dataset",
    x = "Issue",
    y = "Number of Records"
  ) +
  theme_minimal(base_size = 14)

## Remove "CONTINENT_COORDINATE_MISMATCH"
occ_in_target_df <- occ_in_target_df %>%
  filter(
    !GEODETIC_DATUM_INVALID
  )

## Convert list columns to character type 
## (for example, by concatenating the elements into a string)
occ_in_target_df$issue_list <- sapply(occ_in_target_df$issue_list, function(x) paste(x, collapse = ","))

write.csv(occ_in_target_df,"../Gbif-plant-clean.csv", row.names=F)

# Upload to google drive
# Get the target folder first to ensure it exists
speciesObs_folder <- drive_find(pattern = "speciesObs", type = "folder")
if(nrow(speciesObs_folder) > 0) {
  # Upload to Google Drive if folder found
  drive_upload(
    file.path(occ_dir, "Gbif-plant-clean.csv"),
    path = as_id(speciesObs_folder$id[1]),
    name = "Gbif-plant-clean.csv"
  )
} else {
  warning("Could not find 'specieObs' folder in Google Drive. File saved locally only.")
}

## ----------- 4. Record numbers of records after each step -----------
# remember to also update numbers in the Google spreadsheet
num_downloaded <- as.data.frame(table(download_gbif$species))
num_3counties <- as.data.frame(table(occ_in_target_df$species))
colnames(num_downloaded) <- c("Species name", "Occ number (downloaded)")
colnames(num_3counties) <- c("Species name", "Occ number (3 counties)")
this_num_info <- merge(num_downloaded, num_3counties, all = TRUE)
# fill na with 0 
this_num_info[is.na(this_num_info)] <- 0
this_num_info$`Handling person` <- myname

occ_num_info <- rbind(occ_num_info, this_num_info)

# update the Google spreadsheet
write_sheet(occ_num_info, pt_info, sheet = "Occurrence number changes")
