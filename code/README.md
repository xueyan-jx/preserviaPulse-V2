## GBIF_data_pulling.R

The script to query occurrence data from GBIF database.

## Envi_layer_elev_stack.R

The script section 1-5 generates slope, aspect, TRI, and flow accumulation layers from elevation (SRTM 30m) layer
* Load and filter county boundaries
* Mosaic and clip SRTM DEM
* Compute slope, aspect, TRI, flow accumulation, and distance to coast layers

The script section 6 stacks climate and terrain layers
* Load slope, aspect, flow accumulation, distance to coast, solar, and climate layers
* Clips and resamples all layers to match a common raster template provided by Lei to standardizes projection (EPSG:2229), extent, and resolution
* Stacks layers for both current (1980–2010) and future climate scenarios (2040–2070)
* Saves outputs in organized folders under SDM_EnvLayers/Stack_Env/