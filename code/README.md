## Codes

Study area: the Dangermond Preserve in Santa Barbara county located within the tri-county region in southern California (Santa Barbara, Ventura, San Louis Obispo)

### Species occurrence data pulling and cleaning

![test](../visualization/data_pulling_workflow.png "Flowchart")

#### Input data: 
* Integrated Resources Management Plan (IRMP): get special status species
* IUCN red list range maps: get specieal status species
* CalFlora: get specieal status species

#### R codes:
* [01_GBIF_data_pulling.R](../code/01_GBIF_data_pulling.R)
    * Purpose: step 1 -2 grab GBIF data for target species, and step 3 - 4 clean species data.
    * Birds data accessed on 05/18/2025.
    * Mammals data accessed on 05/07/2025.
    * Herps and invertebrates data accessed on 04/30/2025.
    * Plants data accessed on 05/15/2025.
* [02_BIEN_occ_data.R](../code/02_BIEN_occ_data.R)
    * Purpose: to be filled.
* [03_Integrated_occ_dangermond_Portal.R](../code/03_Integrated_occ_dangermond_Portal.R)
    * Purpose: grab data from the Dangermond Data Portal and clean species data.
* [04_GBIF_BIEN_DP_merge.R](../code/04_GBIF_BIEN_DP_Cal_merge.R)
    * Purpose: merge plant species data from multiple sources.
* [05_mergeAnimalOcc.R](../code/05_mergeAnimalOcc.R)
    * Purpose: merge animal species data from multiple sources.
* [06_postProccOcc.R](../code/06_postProccOcc.R)
    * Purpose: final clean up step to get a single model-ready occurrence data file for all species.

#### Output data:
* Final cleaned data for all species, projected to NAD California Zone 5

Example:

| species | x | y |
| -------| --- | --- |
| species latin name | longitude | latitude |

### Envi_layer_elev.R

The script to generate slope, aspect, TRI, and flow accumulation layers from elevation (SRTM 30m) layer
* Load and filter county boundaries
* Mosaic and clip SRTM DEM
* Compute slope, aspect, TRI, and flow accumulation layers
