df <- read.delim('data/occurrences/GBIF0005883-250426092105405-cleaned.csv', sep = ';')
df_3 <- read.delim('data/occurrences/0005883-250426092105405_3counties_pts.csv', sep=';')
colnames(df)

unique(df$species)

df1 <- df %>% filter(species == 'Rynchops niger')
# Icterus bullockiorum - on record missing
# Hydrobates socorroensis
# Urile pelagicus
# Rynchops niger

df2 <- download_gbif %>% filter(species == 'Urile pelagicus')
colnames(download_gbif)

nm_orig <- unique(gbif_names$species)
nm_dwld <- unique(download_gbif$species)

setdiff(nm_orig, nm_dwld) # Icterus bullockiorum removed
setdiff(nm_dwld, nm_orig) # Synthliboramphus hypoleucus removed

"Synthliboramphus hypoleucus" %in% nm_dwld

nm_3counties <- unique(df_3$species)
setdiff(nm_dwld, nm_3counties)

nm_final <- unique(df$species)
setdiff(nm_3counties, nm_final)
