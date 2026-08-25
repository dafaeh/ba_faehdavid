# Bachelor Thesis: The cumulative water deficit at 30 m across the contiguous United States – Spatial patterns and landscape controls.

**Degree program:** Geography

**Author:** David Fäh

**Supervision:** Benjamin David Stocker, Fabian Bernhard

## Overview
- In this project, a cloud-based Google Earth Engine (GEE) workflow is
  developed to produce a 30 m resolution map of maximum seasonal cumulative
  water deficits (CWDmax) over the course of 3 years across the contiguous United States.
- Based on this product, it is tested how different drivers affect CWDmax. 
  These drivers are: topographic position, geological substrate, land
  cover and irrigation.

## Data sources
- Evapotranspiration: OpenET DisALEXI model providing actual ET at 30 m resolution
  (Melton et al., 2022): https://doi.org/10.1111/1752-1688.12956
- Precipitation: DAYMET gridded daily precipitation at 1 km resolution
  (Thornton et al., 2021): https://doi.org/10.1038/s41597-021-00973-0
- Topography: FABDEM 30 m bare-earth DTM
  (Hawker et al., 2022): https://doi.org/10.1088/1748-9326/ac4d4f
- Geology: USGS State Geologic Map Compilation (SGMC) at 90 m resolution
  (Horton et al,. 2017): https://doi.org/10.3133/ds1052
- Land use: National Land Cover Database (NLCD) at 30 m resolution
  (Dewitz 2021): https://doi.org/10.5066/P9KZCM54, supplemented with the
  Landsat-based Irrigation Dataset (LANID) at 30 m resolution
  (Xie et al. 2021): https://doi.org/10.5194/essd-13-5689-2021 to distinguish 
  irrigated from non-irrigated cropland.
  
## Project structure
- `gee/`: GEE workflow that computes CWDmax (see `gee/README.md` for further 
  information on how the GEE pipeline works)
- `analysis/`: numbered R workflow scripts (run in order)
- `data-raw/`: immutable input data
- `data/`: processed data, including CWDmax exports from GEE
- `fig/`: generated figures
- `R/`: R functions used in multiple scripts

Note that the contents of `data` and `data-raw` are too large for GitHub. The can be obtained
using the 01_download_data.R script or manually downloaded.
  
## How to run this project
The workflow is separated into two sections.

1. The data for the analysis is generated using GEE. See `gee/README.md`
   for details and execute these steps first.
2. Once the GEE outputs are stored in `data/`, run the scripts in
   `analysis/` in their numerical order to produce the plots used in the thesis and 
   assemble the CONUS-wide CWDmax map.
3. The maps that were used in the thesis can be reproduced using the scripts in 
  `analysis/maps/`

