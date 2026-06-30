// config.js
// Central configuration for the cwd_max workflow.
// Edit SITES, SPINUP_YEAR, ANALYSIS_YEARS, and REFERENCE_YEARS between runs.


// Study site geometries (WGS84). Add or remove sites here.
exports.SITES = {
  eel: ee.Geometry.Rectangle(
    [-123.74888168761035, 39.29196540299779,
     -121.09019028136035, 41.29344890598067]
  ),
};

// SPINUP_YEAR initialises the snow pool. Its cwd_max output is discarded.
exports.SPINUP_YEAR    = 2019;
exports.ANALYSIS_YEARS = [2020, 2021, 2022];

// Years used to build the per-month ET reference mean for gap-filling.
exports.REFERENCE_YEARS = [2016, 2017, 2018, 2019, 2020, 2021, 2022, 2023, 2024];

// Snow model parameters
exports.TEMP_THRESHOLD = 1.0;   // °C (snow/rain partitioning threshold)
exports.MAX_MELT_RATE  = 1.0;   // mm °C⁻¹ day⁻¹

// Export settings.
exports.EXPORT_FOLDER = 'gee_cwd_exports';
exports.EXPORT_CRS    = 'EPSG:5070';   // CONUS Albers Equal Area
exports.EXPORT_SCALE  = 30;            // metres — native DisALEXI resolution
exports.MAX_PIXELS    = 1e9;

// GEE asset paths.
exports.OPENET_ASSET = 'projects/openet/assets/disalexi/conus/gridmet/monthly/v2_1';
exports.DAYMET_ASSET = 'NASA/ORNL/DAYMET_V4';