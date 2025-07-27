## Purpose of script:
## 1. Find overlap between JLDP and CWHR
## 2. Overlap species occurrence and each CWHR and determine species group for each habitat


#1. Overlap JLDP and each CWHR to determine which habitats we want to keep 
# (e.g., if the habitat range covers over 50% of the area of the preserve, then we can keep); 
# 2. Crop each CWHR to the three countries, and then overlap species occurrence and each CWHR and determine species group for each habitat


## Authors: Xue Yan
## Date: Summer, 2025
library(here)
library(dplyr)
library(readxl)
library(terra)
library(sf)
library(purrr)
library(fs) # for path_ext_remove
library(tools)
# ------------- 1. Data preparation -------
# ------------- (1) Env and boundary -------
## Environment data (for getting CRS and resolution)
ref_Env <- rast(here("data", "Stack_Env", "combined_stack_hist.tif"))
ref_crs_res <- sf::st_crs(crs(ref_Env))

## JLDP
JLDP <- st_read(here("data", "JLDP", "JLDP_boundary2024update.shp"))
JLDP_transformed <- st_transform(JLDP, crs = ref_crs_res)

## CWHR
CWHR_dir <- here("data", "CWHR_habitat","CWHR_Habitats_GIS_2014")
CWHR_files <- list.files(
  path = CWHR_dir,
  pattern = "\\.shp$",
  full.names = TRUE,
  ignore.case = TRUE,
  recursive = TRUE
)

CWHR_list <- purrr::map(CWHR_files, ~ {
  st_read(.x, quiet = TRUE) |> 
    st_transform(ref_crs_res)
})

names(CWHR_list) <- path_file(CWHR_files) |> path_ext_remove()

## Test if CRSs are same
crs(JLDP_transformed) == crs (CWHR_list[["ADS"]])
crs (CWHR_list[["ADS"]]) == crs(ref_Env)

# ------------- 2. Check overlap between JLDP and CWHR -------
# --------------(1) Any overlap------------
habitat_names <- tools::file_path_sans_ext(basename(CWHR_files))

shp_list_overlap <- purrr::map2(CWHR_list, habitat_names, ~ {
  st_filter(.x, JLDP_transformed, .predicate = st_intersects) |>
    dplyr::mutate(habitat = .y)
}) #need to check habitat distribution, is there overlap between habitats

all_habitats_overlap <- bind_rows(shp_list_overlap)|>
  select(-any_of(c("id", "ID","Id")))


# --------------2. "Endemism" of the habitats in the Preserve----
# ------------- (2) Species occurrence -------
# All species occurrence (x and y are under EPSG:2229)
Anim_Plant_occ <- read.csv(here("data","occurrences","Anim_Plant_merge_final0711.csv"))
occ_sf <- st_as_sf(Anim_Plant_occ, coords = c("x", "y"), crs = st_crs(JLDP_transformed))

# Filter to species in IRMP
IRMP_species <- read_excel(here("data", "IRMP-species.xlsx")) %>%
  pull(Species)

occ_filtered <- occ_sf %>%
  filter(species %in% IRMP_species)

occ_in_JLDP <- st_filter(occ_filtered , JLDP_transformed, .predicate = st_intersects)
occ_with_habitat <- st_join(occ_in_JLDP, all_habitats_overlap, left = FALSE)


# --------------(2)"Endemism" of the habitats in the Preserve----
species_occ_counts <- occ_with_habitat %>%
  group_by(species) %>%
  summarise(n_occurrences = n())

occ_with_we <- occ_with_habitat %>%
  mutate(n_occurrences = species_occ_counts$n_occurrences[match(species, species_occ_counts$species)],
         WE = 1 / n_occurrences)

habitat_we <- occ_with_we %>%
  group_by(habitat) %>%
  summarise(total_WE = sum(WE, na.rm = TRUE))


