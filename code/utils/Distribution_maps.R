## Purpose of script: map occurrence points for a species
## Authors: GEOG 274
## Date: Spring, 2025

# ------------- Setting up --------------
# Install and load packages
library(dplyr)
library(tigris) # to download county boundary
library(ggplot2)
library(here)
library(sf)
library(leaflet)
library(htmlwidgets)
library(stringr)

# Set a directory for data
here() # first check path
a_dir <- here("visualization/private") # create a data folder with relative path
file.path(a_dir) # double check its absolute path
if (!dir.exists(a_dir)) dir.create(a_dir) # if it's not there already, create it


# Get boundaries for areas of interest
ca_counties <- counties(state = "CA", cb = TRUE) # gather CA county boundary file
target_counties <- ca_counties %>%
  filter(NAME %in% c("Santa Barbara", "Ventura", "San Luis Obispo")) # limit to counties of interest


# ------------------- get a list of all species --------------------
ss <- drive_get("Special Status Species")
all_names <- as.data.frame(matrix(ncol=3, nrow=0))
colnames(all_names) <- c('Name (Latin)', 'Name (common)', 'Where Listed?')
for(taxon in c('Birds')){
  dat <- read_sheet(ss, sheet=taxon) %>% select(`Name (Latin)`, `Name (common)`,
                                                'Where Listed?')
  all_names <- rbind(all_names, dat)
}

all_names <- all_names %>% mutate(
  `Name (Latin)` = ifelse(`Name (Latin)`=='Vulpes vulpes ssp.', 'Vulpes vulpes', `Name (Latin)`)
) %>% unique()

irmp_sp <- all_names %>% filter(`Where Listed?`=='IRMP')

# ------------- Get bird final info --------------
birds_final_info <- read.csv('data/occurrences/animals/birds_final_num.csv')
colnames(birds_final_info) <- c('Name (Latin)', 'Name (common)', 'Number')

irmp_info <- merge(irmp_sp, birds_final_info)
setdiff(irmp_sp$`Name (Latin)`, irmp_info$`Name (Latin)`)


# ------------- Basic map of distribution for a list of species --------------
# search for a csv that ends with -cleaned.csv and read it
#files <- list.files(path = "data/occurrences/animals", pattern = "-cleaned-0515.csv$", full.names = TRUE)
df <- read.delim('data/occurrences/animals/birds-cleaned-0515.csv', sep = ';')
unique(df$species)

# irmp species
df_irmp <- df %>% filter(species %in% unique(irmp_info$`Name (Latin)`))
df_irmp <- merge(df_irmp, all_names, by.x='species', by.y='Name (Latin)')

# Extract coordinates from geometry string
df_irmp <- df_irmp %>%
  mutate(
    x = as.numeric(str_extract(geometry, "(?<=c\\()[0-9.]+")),
    y = as.numeric(str_extract(geometry, "(?<=, )[0-9.]+"))
  )

# Convert to spatial data
gdf_irmp <- st_as_sf(df_irmp, coords = c("x", "y"), crs = 2229 # same crs to fishnet_clipped
)

# create static maps for target species
p <- ggplot() +
  geom_sf(data = target_counties, fill = "lightblue", color = "black", linewidth = 1) +
  geom_sf(data = gdf_irmp, col = "red", pch = 20) +
  facet_wrap(~`Name (common)`) +
  theme_minimal() +
  theme(legend.position = "none")

# save the static plot
ggsave(filename = "visualization/private/static_map.png", plot = p, width = 10, height = 10)

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
saveWidget(m, file = "visualization/private/interactive_species_map.html")
