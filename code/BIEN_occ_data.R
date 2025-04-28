install.packages("RBIEN")
library(BIEN)

# `BIEN_occurrence_county` Returns all occurrences records within a given state/province

country_vector<-c("United States","United States","United States")
state_vector<-c("California","California","California")
county_vector<-c("Santa Barbara","Ventura","San Luis Obispo")
BIEN_occ<-BIEN_occurrence_county(country=country_vector, state = state_vector, county = county_vector)

