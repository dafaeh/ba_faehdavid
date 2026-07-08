# Bachelor Thesis: High-Resolution Mapping of Cumulative Water Deficits across 
the Contiguous United States

**Degree program:** Geography

**Author:** David Fäh

**Supervision:** Benjamin David Stocker, Fabian Bernhard

## Overview
- In this project, a cloud-based Google Earth Engine (GEE) workflow is
  developed to produce a 30 m resolution map of maximum seasonal cumulative
  water deficits (CWDmax) over the course of 3 years.
- The influence of multiple drivers of CWDmax is analyzed in an explorative
  way. These drivers are: topographic position, geological substrate, land
  cover and irrigation.

## Data sources
- Evapotranspiration: OpenET DisALEXI model providing actual ET at 30 m resolution
  (Melton et al. 2022, https://doi.org/10.1111/1752-1688.12956).
- Precipitation: DAYMET gridded daily precipitation at 1 km resolution
  (Thornton et al. 2021, https://doi.org/10.1038/s41597-021-00973-0).
- Topography: FABDEM 30 m bare-earth DTM
  (Hawker et al. 2022, https://doi.org/10.1088/1748-9326/ac4d4f).
- Geology: USGS State Geologic Map Compilation (SGMC) at 90 m resolution
  (Horton et al. 2017, https://doi.org/10.3133/ds1052).
- Land use: National Land Cover Database (NLCD) at 30 m resolution
  (Dewitz 2023, https://doi.org/10.5066/P9JZ7AO3), supplemented with the
  Landsat-based Irrigation Dataset (LANID) at 30 m resolution
  (Xie et al. 2021, https://doi.org/10.5194/essd-13-5689-2021) to distinguish 
  irrigated from non-irrigated cropland.
  
## Project structure
- `gee/`: GEE workflow that computes CWDmax (see `gee/README.md`)
- `analysis/`: numbered R workflow scripts (run in order)
- `data-raw/`: immutable input data
- `data/`: processed data, including CWDmax exports from GEE
- `fig/`: generated figures
  
## How to run this project
The workflow is separated into two sections.

1. The data for the analysis is generated using GEE. See `gee/README.md`
   for details.
2. Once the GEE outputs are stored in `data/`, run the scripts in
   `analysis/` in the numerical order of their filenames.

