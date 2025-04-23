
## Code from Lei to connect local R with GBIF API
## Removed personal directory and login information

# Materials:
#   https://docs.ropensci.org/rgbif/index.html
#   More details in: https://techdocs.gbif.org/en/openapi/
#   https://docs.gbif.org/course-data-use/en/data-processing-pipeline.html
#  https://data-blog.gbif.org/post/gbif-filtering-guide/
#. https://data-blog.gbif.org/post/downloading-long-species-lists-on-gbif/

########### Setting ##############
# Install and load packages
install.packages("rgbif")
library(here)
library(sf)
library(dplyr)
library(rgbif)

# Set a directory for data
occ_dir <- here("~/local_directory_path_here") # update with local path
if (!dir.exists(occ_dir)) dir.create(occ_dir)

########### Pre-processing ##############
# Confirm names used in GBIF database for your species list
scientific_names <- c("Deinandra increscens ssp. Villosa",
                      "Eriodictyon capitatum")
gbif_names <- name_backbone_checklist(dat$`Name (Latin)`)

# But if you don't know the exact scientific names
nm <- "Grasshopper sparrow"
suggest_names <- name_suggest(nm)

all_info <- name_lookup(nm)$data

# A key attribute is matchType, which is the best to be EXACT.

########### Query and save data ##############
## Add a few internal filters for something that we definitely not want.
## we can do other filters afterwards.
## 1. Zero coordinate : Coordinates are exactly (0,0). null island
## 2. Country coordinate mismatch : The coordinates fall outside of the 
##    given country’s polygon.
## 3. Coordinate invalid : GBIF is unable to interpret the coordinates.
## 4. Coordinate out of range : The coordinates are outside of the range 
##    for decimal lat/lon values ((-90,90), (-180,180)).
## 5. Only keep PRESENT records
## 6. Remove FOSSIL_SPECIMEN and LIVING_SPECIMEN

#remove NA values from usage key before query
usage_keys <- unique(gbif_names$usageKey)
usage_keys <- usage_keys[!is.na(usage_keys)]

occ_meta <- occ_download(
  pred_in("taxonKey", unique(usage_keys)),
  pred("occurrenceStatus", "PRESENT"), 
  pred("hasCoordinate", TRUE),
  pred("hasGeospatialIssue", FALSE),
  pred_not(pred_in("basisOfRecord",
                   c("FOSSIL_SPECIMEN", "LIVING_SPECIMEN"))),
  format = "SIMPLE_CSV",
  user = "user",      # Replace with your actual username
  pwd = "pwd",       # Replace with your actual password
  email = "email"  # Replace with your registered email
)

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