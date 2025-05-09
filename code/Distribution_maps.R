## Purpose of script: map occurrence points for a species
## Authors: GEOG 274
## Date: Spring, 2025
## Credits to: Wenxin Yang, Jacqueline Vogel, Yanni Zhan, Xue Yan, Dr. Lei Song (lsong@ucsb.edu)

# ------------- Setting up --------------
# Install and load packages
library(dplyr)
library(tigris) # to download county boundary
library(ggplot2)
library(here)
library(sf)
library(leaflet)
library(htmlwidgets)

# Set a directory for data
here() # first check path
occ_dir <- here("visualization") # create a data folder with relative path
file.path(occ_dir) # double check its absolute path
if (!dir.exists(occ_dir)) dir.create(occ_dir) # if it's not there already, create it


# Get boundaries for areas of interest
ca_counties <- counties(state = "CA", cb = TRUE) # gather CA county boundary file
target_counties <- ca_counties %>%
  filter(NAME %in% c("Santa Barbara", "Ventura", "San Luis Obispo")) # limit to counties of interest


# ------------- Basic map of distribution for a list of species --------------
# search for a csv that ends with -cleaned.csv and read it
file <- list.files(path = "data/occurrences", pattern = "-cleaned.csv$", full.names = TRUE)
df <- read.delim(file, sep = ';')
unique(df$species)

# see number of records for each specie and sort by Freq
count_info <- as.data.frame(table(df$species)) %>% arrange(desc(Freq))
# define a list of species names, in my case I'm selecting the 5 with most counts
# for some folks, it could be selecting species mentioned by the IRMP
# li_interests <- c('Piranga ludoviciana', 'Larus occidentalis', 'Ardea alba')
li_interests <- count_info$Var1[1:5]
# filter records for specie of interest
df1 <- df %>% filter(species %in% li_interests) %>% select(species, decimalLongitude, decimalLatitude)
# convert to sf object, specify the coordinate reference system and the geometry column
geodf1 <- st_as_sf(df1, crs = 4326, coords = c("decimalLongitude", "decimalLatitude"))

# create static maps for target species
p <- ggplot() +
  geom_sf(data = target_counties, fill = "lightblue", color = "black", linewidth = 1) +
  geom_sf(data = geodf1, col = "red", pch = 20) +
  facet_wrap(~species) +
  theme_minimal() +
  theme(legend.position = "none")

# save the static plot
ggsave(filename = "visualization/static_map.png", plot = p, width = 10, height = 10)

# ------------- Create interactive leaflet map --------------
# create color palette (1 color for each specie)
species_colors <- colorFactor(
  palette = "Set1",
  domain = df1$species
)

# reproject target_counties to WGS84
target_counties_wgs84 <- st_transform(target_counties, 4326)

# create the leaflet map
m <- leaflet() %>%
  # add OpenStreetMap base tiles
  addTiles() %>%
  # set the view to the center of the target counties
  setView(lng = mean(st_coordinates(geodf1)[,1]), 
          lat = mean(st_coordinates(geodf1)[,2]), 
          zoom = 8) %>%
  # add county boundaries
  addPolygons(data = target_counties_wgs84,
              fillColor = "lightblue",
              fillOpacity = 0.3,
              weight = 2,
              color = "black",
              label = ~NAME) 

# add each species as a separate layer group
for(sp in unique(df1$species)) {
  sp_points <- geodf1[geodf1$species == sp,]
  m <- m %>%
    addCircleMarkers(data = sp_points,
                    color = ~species_colors(sp),
                    radius = 5,
                    fillOpacity = 0.7,
                    stroke = FALSE,
                    group = sp,
                    label = ~species,
                    popup = ~paste("Species:", species))
}

# add layer control (panel on the right hand side)
m <- m %>%
  addLayersControl(
    overlayGroups = unique(df1$species),
    options = layersControlOptions(collapsed = FALSE)
  )

# save the interactive map
saveWidget(m, file = "visualization/interactive_species_map.html")
