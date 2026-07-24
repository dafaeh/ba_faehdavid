# GEE Workflow 

Computes CWDmax (30 m resolution) across the contiguous United States (CONUS)
using OpenET (for evapotranspiration) and DAYMET (for precipitation) data. Focus 
regions exported in individual files for further analysis in R. The CONUS-wide 
map is computed and exported in tiles and later mosaicked in R. See the top-level 
README for methodology and background.

## Setup

1. Copy all `.js` files into a new GEE repository.
2. Adjust the `require()` paths (`users/dafaeh/ba_cwd:...`) to match your repo name.
3. Check/adjust `config.js` (`SITES`, `ANALYSIS_YEARS`, `REFERENCE_YEARS, `EXPORT_FOLDER`).

### Adding a region to the export
Regions in `SITES` must be axis-aligned rectangles in EPSG:5070 with
`geodesic = false` like in `config.js`. Rectangles drawn in the GEE Code Editor come out in WGS84 and
have to be transformed to EPSG:5070 first.

## Run

### `main.js` for focus regions

Computes CWDmax for the focus regions defined in `cfg.SITES`. Exports to
Google Drive (folder `cfg.EXPORT_FOLDER`); tasks must be started manually
in the Tasks tab.

### `main_conus.js` for the CONUS-wide map

Lays a tile grid over CONUS and computes CWDmax tile by tile, to stay
within GEE's per-task memory/timeout limits.

- Comment/uncomment tile entries in the script to select which tiles get exported.
- Use the Inspector tool in the GEE Code Editor to click on a tile and read
  off its name/ID.
- Some tiles fall entirely outside CONUS (over ocean or bordering
  countries). These are marked in the code under `excluded tiles` and
  should stay commented out.

After all tiles are exported, mosaic them into a single CONUS-wide map using the
`conus_merge.R` script.

## Runtime & budget

Reference: the Appalachia focus region (~197,000 km²) took 76.8 EECU-h
to compute (with 9 ET reference years for gap filling) and 30–40 minutes for the export
to finish. Use this as a rough per-km² baseline (~3.9 × 10⁻⁴ EECU-h/km²) to estimate other runs.

This workflow was run on a GEE contributor account with a quota of
1000 EECU-h/month. This does not appear to be a hard limit. Further usage
beyond the quota was still possible in practice.