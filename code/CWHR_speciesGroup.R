## Purpose of script:
## 1. Find overlap between JLDP and CWHR
## 2. Overlap species occurrence and each CWHR and determine species group for each habitat


#1. Overlap JLDP and each CWHR to determine which habitats we want to keep 
# (e.g., if the habitat range covers over 50% of the area of the preserve, then we can keep); 
# 2. Crop each CWHR to the three countries, and then overlap species occurrence and each CWHR and determine species group for each habitat


## Authors: Xue Yan
## Date: Summer, 2025
library(here)
library(tigris) 
library(dplyr)
library(readxl)
library(terra)
library(sf)
library(purrr)
library(fs) # for path_ext_remove
library(tools)
library(ggplot2)
# ------------- 1. Data preparation -------
# ------------- (1) Env and boundary -------
## Environment data (for getting CRS and resolution)
ref_Env <- rast(here("data", "Stack_Env", "combined_stack_hist.tif"))
ref_crs_res <- sf::st_crs(crs(ref_Env))

## Three counties
#ca_counties <- tigris::counties(state = "CA", cb = TRUE) # gather CA county boundary file
#target_counties <- ca_counties %>%
#  filter(NAME %in% c("Santa Barbara", "Ventura", "San Luis Obispo")) 
us_counties <- sf::st_read(here("data", "Boundary", "cb_2024_us_county_500k", "cb_2024_us_county_500k.shp"))
ca_counties <- us_counties[us_counties$STATEFP == "06", ]  # 06 is code for California
target_counties <- ca_counties %>%
  filter(NAME %in% c("Santa Barbara", "Ventura", "San Luis Obispo")) %>%
  st_transform(ca_counties, crs = st_crs(ref_crs_res))

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
crs (CWHR_list[["ADS"]]) == crs(target_counties)

# ------------- (2) Species occurrence -------
# All species occurrence (x and y are under EPSG:2229)
Anim_Plant_occ <- read.csv(here("data","occurrences","Anim_Plant_merge_final0711.csv"))
occ_sf <- st_as_sf(Anim_Plant_occ, coords = c("x", "y"), crs = st_crs(JLDP_transformed))

# Filter to species in IRMP
IRMP_species <- read_excel(here("data", "IRMP-species.xlsx")) %>%
  pull(Species)

occ_filtered <- occ_sf %>%
  filter(species %in% IRMP_species) #target species

occ_in_JLDP <- st_filter(occ_filtered , JLDP_transformed, .predicate = st_intersects) #target species appear in JLDP


# ------------- 2. Any overlap -------
# --------------(1) Check overlap between JLDP and CWHR------------
habitat_names <- tools::file_path_sans_ext(basename(CWHR_files))

shp_list_overlap <- purrr::map2(CWHR_list, habitat_names, ~ {
  st_filter(.x, JLDP_transformed, .predicate = st_intersects) |>
    dplyr::mutate(habitat = .y)
}) |>
  purrr::discard(~ nrow(.x) == 0) 


all_habitats_overlap <- bind_rows(shp_list_overlap)|>
  select(-any_of(c("id", "ID","Id")))

ggplot(bind_rows(shp_list_overlap)) +
  geom_sf(fill = "skyblue") +
  facet_wrap(~ habitat) +
  theme_minimal()

############################################
#Check habitat distribution, is there overlap between habitats
# Clip to JLDP range
habitat_names_overlap <- names(shp_list_overlap)

shp_list_clipped <- map2(shp_list_overlap, habitat_names_overlap, ~ {
  st_intersection(.x, st_geometry(JLDP_transformed)) |> 
    mutate(habitat = .y)
})

shp_combined <- bind_rows(shp_list_clipped)
ggplot(shp_combined) +
  geom_sf(fill = "skyblue") +
  facet_wrap(~ habitat) +
  theme_minimal()


# Clip to three counties range
habitat_names_overlap <- names(shp_list_overlap)

shp_list_clipped1 <- map2(shp_list_overlap, habitat_names_overlap, ~ {
  st_intersection(.x, st_geometry(target_counties)) |> 
    mutate(habitat = .y)
})

shp_combined1 <- bind_rows(shp_list_clipped1)
ggplot(shp_combined) +
  geom_sf(fill = "skyblue") +
  facet_wrap(~ habitat) +
  theme_minimal()
  
##########################################
# --------------(2) Find species group------------
# -------------a. For JLDP------------------------
VIN <- shp_combined %>%
  dplyr::filter(habitat == "VIN")

occ_in_VIN <- st_filter(occ_in_JLDP , VIN, .predicate = st_intersects) 
unique_species <- unique(occ_in_VIN$species)

### Take a look at all species
habitat_types <- unique(shp_combined$habitat)
total_species_count_JLDP <- occ_in_JLDP %>%
  st_drop_geometry() %>%
  count(species, name = "total_count")

# Create dataframe to store occ for habitats
JLDPhabitat_species_stats <- data.frame()

for (hab in habitat_types) {
  # Extract habitat polygon
  habitat_poly <- shp_combined %>% filter(habitat == hab)
  
  # occ in the habitat
  occ_in_hab <- st_filter(occ_in_JLDP, habitat_poly, .predicate = st_intersects)
  
  # stat of each species appear in the habitat
  species_count_in_hab <- occ_in_hab %>%
    st_drop_geometry() %>%
    count(species, name = "habitat_count") %>%
    mutate(habitat = hab)
  
  # percentage of the occ in the habitat
  species_count_in_hab <- species_count_in_hab %>%
    left_join(total_species_count, by = "species") %>%
    mutate(proportion = habitat_count / total_count)
  
  JLDPhabitat_species_stats <- bind_rows(JLDPhabitat_species_stats, species_count_in_hab)
}

write.csv(JLDPhabitat_species_stats, here::here("results", "habitat_species_summary_JLDP.csv"), row.names = FALSE)

# -------------a. For three counties------------------------
habitat_types <- unique(shp_combined1$habitat)
total_species_count <- occ_filtered %>%
  st_drop_geometry() %>%
  count(species, name = "total_count")

# Create dataframe to store occ for habitats
habitat_species_stats <- data.frame()

for (hab in habitat_types) {
  # Extract habitat polygon
  habitat_poly <- shp_combined1 %>% filter(habitat == hab)
  
  # occ in the habitat
  occ_in_hab <- st_filter(occ_filtered, habitat_poly, .predicate = st_intersects)
  
  # stat of each species appear in the habitat
  species_count_in_hab <- occ_in_hab %>%
    st_drop_geometry() %>%
    count(species, name = "habitat_count") %>%
    mutate(habitat = hab)
  
  # percentage of the occ in the habitat
  species_count_in_hab <- species_count_in_hab %>%
    left_join(total_species_count, by = "species") %>%
    mutate(proportion = habitat_count / total_count)
  
  habitat_species_stats <- bind_rows(habitat_species_stats, species_count_in_hab)
}

write.csv(habitat_species_stats, here::here("results", "habitat_species_summary_threeCounties.csv"), row.names = FALSE)


# --------------2. "Endemism" of the habitats ----
occ_with_habitat <- st_join(occ_in_JLDP, all_habitats_overlap, left = FALSE)


# --------------(1) Weighted Endemism of the habitats in the Preserve----
species_occ_counts <- occ_with_habitat %>%
  group_by(species) %>%
  summarise(n_occurrences = n())

occ_with_we <- occ_with_habitat %>%
  mutate(n_occurrences = species_occ_counts$n_occurrences[match(species, species_occ_counts$species)],
         WE = 1 / n_occurrences)

habitat_we <- occ_with_we %>%
  group_by(habitat) %>%
  summarise(total_WE = sum(WE, na.rm = TRUE))

# --------------(2)Weighted Endemism of the habitats for three conties----
occAll_with_habitat <- st_join(occ_filtered, all_habitats_overlap, left = FALSE)

species_occAll_counts <- occAll_with_habitat %>%
  group_by(species) %>%
  summarise(n_occurrences = n())

occAll_with_we <- occAll_with_habitat %>%
  mutate(n_occurrences = species_occAll_counts$n_occurrences[match(species, species_occAll_counts$species)],
         WE = 1 / n_occurrences)

habitat_we_All <- occAll_with_we %>%
  group_by(habitat) %>%
  summarise(total_WE = sum(WE, na.rm = TRUE))

# --------------(3) Find species group------------
