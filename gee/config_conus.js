// config_conus.js
// CONUS-wide configuration for the cwd_max workflow.
// Defines a 10x8 tile grid over the contiguous United States in EPSG:5070.
//
// Grid parameters:
//   CONUS_XMIN : -2 357 010 m
//   CONUS_XMAX :  2 258 190 m   (adjusted so X-span is divisible by 10 x 30 m)
//   CONUS_YMIN :    268 980 m
//   CONUS_YMAX :  3 172 020 m
//   tile_w     :    461 520 m  =  15 384 px at 30 m
//   tile_h     :    362 880 m  =  12 096 px at 30 m
//   px / tile  :  ~186 million  ->  no GEE file split
//   tiles      :  80 total  (66 active, 14 excluded)
//
// Tile naming: r<row>_c<col>, row 01 = southernmost, col 01 = westernmost.
//
// Workflow:
//   Before each run, comment out tiles not needed for that run.
//   Commit the state of this file to Git before submitting tasks so
//   the export history is fully reproducible.
//   See CONUS_Tiles.md for the run schedule.

// -- Grid parameters ----------------------------------------------------------
var CONUS_XMIN = -2357010;
var CONUS_YMIN =   268980;
var tile_w     =   461520;  // m  (461 520 / 30 = 15 384 px)
var tile_h     =   362880;  // m  (362 880 / 30 = 12 096 px)

// Computes the EPSG:5070 rectangle for a tile key (e.g. "r03_c07").
// geodesic = false is required for rectangles in a projected CRS.
function build_tile(key) {
  var parts = key.split("_");
  var r     = parseInt(parts[0].slice(1), 10);
  var c     = parseInt(parts[1].slice(1), 10);
  var xmin  = CONUS_XMIN + (c - 1) * tile_w;
  var ymin  = CONUS_YMIN + (r - 1) * tile_h;
  return ee.Geometry.Rectangle(
    [xmin, ymin, xmin + tile_w, ymin + tile_h], "EPSG:5070", false
  );
}

// -- Active tiles (66) --------------------------------------------------------
// Before each run, comment out tiles not needed and commit to Git.
// Tiles are listed north to south, west to east within each row.
var SITES = {};

// Row 8 (northernmost)
/*
SITES['r08_c01'] = build_tile('r08_c01');
SITES['r08_c02'] = build_tile('r08_c02');
SITES['r08_c03'] = build_tile('r08_c03');
SITES['r08_c04'] = build_tile('r08_c04');
SITES['r08_c05'] = build_tile('r08_c05');
SITES['r08_c06'] = build_tile('r08_c06');
*/
SITES['r08_c10'] = build_tile('r08_c10');
/*

// Row 7
SITES['r07_c01'] = build_tile('r07_c01');
SITES['r07_c02'] = build_tile('r07_c02');
SITES['r07_c03'] = build_tile('r07_c03');
SITES['r07_c04'] = build_tile('r07_c04');
SITES['r07_c05'] = build_tile('r07_c05');
SITES['r07_c06'] = build_tile('r07_c06');
SITES['r07_c07'] = build_tile('r07_c07');
SITES['r07_c08'] = build_tile('r07_c08');
*/
SITES['r07_c09'] = build_tile('r07_c09');
SITES['r07_c10'] = build_tile('r07_c10');
/*

// Row 6
SITES['r06_c01'] = build_tile('r06_c01');
SITES['r06_c02'] = build_tile('r06_c02');
SITES['r06_c03'] = build_tile('r06_c03');
SITES['r06_c04'] = build_tile('r06_c04');
SITES['r06_c05'] = build_tile('r06_c05');
SITES['r06_c06'] = build_tile('r06_c06');
SITES['r06_c07'] = build_tile('r06_c07');
SITES['r06_c08'] = build_tile('r06_c08');
*/
SITES['r06_c09'] = build_tile('r06_c09');
SITES['r06_c10'] = build_tile('r06_c10');
/*
// Row 5
SITES['r05_c01'] = build_tile('r05_c01');
SITES['r05_c02'] = build_tile('r05_c02');
SITES['r05_c03'] = build_tile('r05_c03');
SITES['r05_c04'] = build_tile('r05_c04');
SITES['r05_c05'] = build_tile('r05_c05');
SITES['r05_c06'] = build_tile('r05_c06');
SITES['r05_c07'] = build_tile('r05_c07');
SITES['r05_c08'] = build_tile('r05_c08');
*/
SITES['r05_c09'] = build_tile('r05_c09');
SITES['r05_c10'] = build_tile('r05_c10');
/*
// Row 4
SITES['r04_c01'] = build_tile('r04_c01');
SITES['r04_c02'] = build_tile('r04_c02');
SITES['r04_c03'] = build_tile('r04_c03');
SITES['r04_c04'] = build_tile('r04_c04');
SITES['r04_c05'] = build_tile('r04_c05');
SITES['r04_c06'] = build_tile('r04_c06');
SITES['r04_c07'] = build_tile('r04_c07');
SITES['r04_c08'] = build_tile('r04_c08');
*/
SITES['r04_c09'] = build_tile('r04_c09');
SITES['r04_c10'] = build_tile('r04_c10');
/*
// Row 3
SITES['r03_c01'] = build_tile('r03_c01');
SITES['r03_c02'] = build_tile('r03_c02');
SITES['r03_c03'] = build_tile('r03_c03');
SITES['r03_c04'] = build_tile('r03_c04');
SITES['r03_c05'] = build_tile('r03_c05');
SITES['r03_c06'] = build_tile('r03_c06');
SITES['r03_c07'] = build_tile('r03_c07');
SITES['r03_c08'] = build_tile('r03_c08');
*/
SITES['r03_c09'] = build_tile('r03_c09');
/*
// Row 2
SITES['r02_c04'] = build_tile('r02_c04');
SITES['r02_c05'] = build_tile('r02_c05');
SITES['r02_c06'] = build_tile('r02_c06');
SITES['r02_c07'] = build_tile('r02_c07');
SITES['r02_c08'] = build_tile('r02_c08');
*/
SITES['r02_c09'] = build_tile('r02_c09');
/*
// Row 1 (southernmost)
SITES['r01_c05'] = build_tile('r01_c05');
SITES['r01_c06'] = build_tile('r01_c06');
SITES['r01_c08'] = build_tile('r01_c08');
*/
SITES['r01_c09'] = build_tile('r01_c09');


// -- Excluded tiles (14) ------------------------------------------------------
// No CONUS land coverage. Uncomment a line to re-enable the tile if needed.
// SITES['r08_c07'] = build_tile('r08_c07');
// SITES['r08_c08'] = build_tile('r08_c08');
// SITES['r08_c09'] = build_tile('r08_c09');
// SITES['r03_c10'] = build_tile('r03_c10');
// SITES['r02_c01'] = build_tile('r02_c01');
// SITES['r02_c02'] = build_tile('r02_c02');
// SITES['r02_c03'] = build_tile('r02_c03');
// SITES['r02_c10'] = build_tile('r02_c10');
// SITES['r01_c01'] = build_tile('r01_c01');
// SITES['r01_c02'] = build_tile('r01_c02');
// SITES['r01_c03'] = build_tile('r01_c03');
// SITES['r01_c04'] = build_tile('r01_c04');
// SITES['r01_c07'] = build_tile('r01_c07');
// SITES['r01_c10'] = build_tile('r01_c10');

exports.SITES = SITES;

// -- Analysis period ----------------------------------------------------------
exports.SPINUP_YEAR    = 2019;
exports.ANALYSIS_YEARS = [2020, 2021, 2022];

// -- Reference years ----------------------------------------------------------
// Broader baseline than config.js to maximise gap-fill coverage across
// the full CONUS domain where persistent Landsat gaps are more frequent.
exports.REFERENCE_YEARS = [2016, 2017, 2018, 2019, 2020, 2021, 2022, 2023, 2024];

// -- Snow model parameters ----------------------------------------------------
exports.TEMP_THRESHOLD = 1.0;   // degrees C
exports.MAX_MELT_RATE  = 1.0;   // mm per degree C per day

// -- Export settings ----------------------------------------------------------
exports.EXPORT_FOLDER = 'gee_cwd_exports_conus';
exports.EXPORT_CRS    = 'EPSG:5070';
exports.EXPORT_SCALE  = 30;
exports.MAX_PIXELS    = 1e9;

// -- GEE asset paths ----------------------------------------------------------
exports.OPENET_ASSET =
  'projects/openet/assets/disalexi/conus/gridmet/monthly/v2_1';
exports.DAYMET_ASSET = 'NASA/ORNL/DAYMET_V4';

/*
// -- Visualisation (temporary, comment out before committing) ----------------------
Map.setCenter(-96, 38, 4);
Object.keys(SITES).forEach(function(key) {
  Map.addLayer(SITES[key], {color: '0000ff'}, key);
});
*/