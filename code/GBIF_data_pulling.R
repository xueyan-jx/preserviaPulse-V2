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
install.packages("CoordinateCleaner")
install.packages("rnaturalearthdata")
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
library(CoordinateCleaner) # add spatial filters
library(readr) # for using write_csv


# Set name
myname <- 'Your-name'

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

dat <- read_sheet(ss, sheet = sheet_names[2]) # Select the taxa you're working on
head(dat)

# for birds, we assigned individuals to different species so we subset the data frame
# skip the following line if you're not working on birds...
# dat <- dat %>% filter(Name==myname) # replace with your name

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

num_downloaded <- as.data.frame(table(download_gbif$species))
colnames(num_downloaded) <- c("Species name", "Occ number (downloaded)")

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

write.table(occ_in_target_df, 
          file.path(occ_dir, paste0(download_id, "_3counties_pts.csv")),
          row.names = FALSE, 
          sep = ";",
          quote = TRUE) # changed how we save file to properly handle the geometry column

# Xue's comment: In case write.table also doesn't work,
# Try the below code to write the table—it should write the correct file out
# It's weird when I used write.table, the whole table was messed up
write_csv(occ_in_target_df,
          file.path(occ_dir, paste0(download_id, "_3counties_pts.csv")))

# get # of records after this step
num_3counties <- as.data.frame(table(occ_in_target_df$species))
colnames(num_3counties) <- c("Species name", "Occ number (3 counties)")

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

## ----------- 3. Data cleaning steps -----------
# if need to re-read the csv file
# occ_in_target_df <- read.csv(file.path(occ_dir, paste0(download_id, "_3counties_pts.csv")), sep=";")

### ----------- (1) coordinateUncertainty filter -----------
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
# Wenxin comment: we also decided to keep those with NA because they may not be from iNat, right?
occ_in_target_df <- occ_in_target_df %>%
  filter(coordinateUncertaintyInMeters <=1000 | is.na(coordinateUncertaintyInMeters))

num_coordUnc <- as.data.frame(table(occ_in_target_df$species))
colnames(num_coordUnc) <- c("Species name", "Occ number (coord unc)")

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
# disable scientific notation
options(scipen = 100)

ggplot(issue_counts, aes(x = reorder(issue, -count), y = count)) +
  geom_bar(stat = "identity") +
  coord_flip() +  # Flip for easier reading
  labs(
    title = "Count of Issues in GBIF Dataset",
    x = "Issue",
    y = "Number of Records"
  ) +
  theme_minimal(base_size = 14)

### ----------- (2) issues filter -----------
# remove occurrence points with issues we don't want
# Wenxin comment: I though we're keeping "GEODETIC_DATUM_INVALID"?
li_issues_remove <- c("CONTINENT_COORDINATE_MISMATCH") # create a list storing issues to remove, we need to decide on which ones to remove
occ_in_target_df <- occ_in_target_df %>% filter(!grepl(paste(li_issues_remove, collapse = "|"), issue_list))

num_issue <- as.data.frame(table(occ_in_target_df$species))
colnames(num_issue) <- c("Species name", "Occ number (issue)")

## Convert list columns to character type 
## (for example, by concatenating the elements into a string)
occ_in_target_df$issue_list <- sapply(occ_in_target_df$issue_list, function(x) paste(x, collapse = ","))


### ----------- (3) temporal filters -----------
# remove invalid years & those before 1980
occ_in_target_df <- occ_in_target_df %>% filter(!is.na(year) & year>=1980)
unique(occ_in_target_df$year)

num_time <- as.data.frame(table(occ_in_target_df$species))
colnames(num_time) <- c("Species name", "Occ number (time)")

### ----------- (4) other spatial filters -----------
## References: 
## https://onlinelibrary.wiley.com/doi/full/10.1002/ece3.9168
## CoordinateCleaner package
## https://peerj.com/articles/9916/
## 1. Lon and Lat are equal ✗ (already accounted for by the spatial filter)
## 2. Duplicate records (species, lon, lat, year, month, day) ✓
## 3. Country capital ✗ (already accounted for by the spatial filter)
## 4. Country centroids ✗ (already accounted for by the spatial filter)
## 5. GBIF headquarters ✗ (already accounted for by the spatial filter)
## 6. Biodiversity institutions ✓
## 7. Geographic outliers ✓
## 8. Sea coordinates ✗
## 9. Urban areas ✗
## 10. dd.mm to dd.dd conversion errors ✗
## 11. Rasterized collections ✗ (this is a test for individual species)
## 12. Suspicious individual count ✓ (definitely remove 0s)

## Remove duplicate records with the same species, lon, lat, year, month, day
occ_in_target_df <- occ_in_target_df %>%
  group_by(species, decimalLongitude, decimalLatitude, year, month, day) %>%
  filter(row_number() == 1) %>%
  ungroup()

num_dup <- as.data.frame(table(occ_in_target_df$species))
colnames(num_dup) <- c("Species name", "Occ number (dup)")


## CoordinateCleaner Wrapper
#flag problems
flag_cols_keep <- c('gbifID', '.val', '.inst') # for coordinate validity & biodiversity institutions
flags <- clean_coordinates(x = occ_in_target_df, 
                           lon = "decimalLongitude", 
                           lat = "decimalLatitude",
                           countries = "countryCode",
                           species = "species",
                           tests = c("institutions"
                                     # "capitals", "centroids","gbif", 
                                     # "equal", "zeros",
                                     ))[flag_cols_keep]

# remove rows with invalid coordinates
ids_coordinvalid <- flags[flags$`.val`=='FALSE',]$gbifID
occ_in_target_df <- occ_in_target_df %>% filter(!gbifID %in% ids_coordinvalid)
num_inv <- as.data.frame(table(occ_in_target_df$species))
colnames(num_inv) <- c("Species name", "Occ number (invalid coord)")

# remove rows recorded in biodiversity institutions
ids_biodivinst <- flags[flags$`.inst`=='FALSE',]$gbifID
occ_in_target_df <- occ_in_target_df %>% filter(!gbifID %in% ids_biodivinst)
num_biodivinst <- as.data.frame(table(occ_in_target_df$species))
colnames(num_biodivinst) <- c("Species name", "Occ number (biodiv inst)")


## Suspicious individual count
# Wenxin comment: some tutorials also marks records with >=99 individual counts as suspicious, animal group (for now) decides only to remove records with 0 counts
table(occ_in_target_df$individualCount)
occ_in_target_df <- occ_in_target_df %>%
  filter(individualCount > 0 | is.na(individualCount)) 
#%>%
 # filter(individualCount < 99 | is.na(individualCount))
num_suscount <- as.data.frame(table(occ_in_target_df$species))
colnames(num_suscount) <- c("Species name", "Occ number (susp count)")


# save the cleaned dataframe
write.table(occ_in_target_df,
            file.path(occ_dir, paste0('GBIF', download_id, '-cleaned.csv')), row.names=FALSE,
            sep=';', quote=TRUE)

# Run the below code in case write.table doesn't work
write_csv(occ_in_target_df,
          file.path(occ_dir, paste0(download_id, "-cleaned.csv")))

# Upload to google drive
# Get the target folder first to ensure it exists
speciesObs_folder <- drive_find(pattern = "speciesObs", type = "folder")
if(nrow(speciesObs_folder) > 0) {
  # Upload to Google Drive if folder found
  drive_upload(
    file.path(occ_dir, paste0('GBIF', download_id, '-cleaned.csv')),
    path = as_id(speciesObs_folder$id[1]),
    name = paste0('GBIF', download_id, '-cleaned.csv')
  )
} else {
  warning("Could not find 'speciesObs' folder in Google Drive. File saved locally only.")
}

## ----------- 4. Record numbers of records after each step -----------
# remember to also update numbers in the Google spreadsheet
this_num_info <- Reduce(function(x, y) merge(x, y, all = TRUE),
                     list(num_downloaded, num_3counties, num_coordUnc, 
                          num_issue, num_time, num_dup, 
                          num_inv, num_biodivinst, num_suscount))

# Fill NA values with 0
this_num_info[is.na(this_num_info)] <- 0
# Add handling person
this_num_info$`Handling person` <- myname

# Update occ_num_info with this_num_info
# first remove any existing rows with myname
occ_num_info <- occ_num_info %>% filter(`Handling person`!=myname)
occ_num_info <- rbind(occ_num_info, this_num_info)

# Fill any remaining NA values with 0
occ_num_info[is.na(occ_num_info)] <- 0

# move handling_person to the last column
occ_num_info <- occ_num_info %>%
  select(colnames(occ_num_info)[colnames(occ_num_info) != "Handling person"], "Handling person")

# sort handling person column by alphabetical order
occ_num_info <- occ_num_info %>%
  arrange(`Handling person`)

# update the Google spreadsheet
write_sheet(occ_num_info, pt_info, sheet = "Occurrence number changes")
