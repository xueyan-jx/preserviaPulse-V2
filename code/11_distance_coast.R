ecu<-st_read("path to ECU shapefile")

ecu <- st_transform(ecu, crs="EPSG:2229")

variable_name_for_observations<-st_transform(variable_name_for_observations, crs="EPSG:2229") #change "variable_name_for_observations" to the variable name associated with the species observations

dist_to_coast <- sf::st_distance(variable_name_for_observations, ecu) #finds distance to coast (US survey foot)

min_distances <- apply(dist_to_coast, 1, min) #finds minimum distance to coast

variable_name_for_observations$distance_to_coast <- min_distances
