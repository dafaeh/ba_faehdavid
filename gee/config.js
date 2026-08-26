// config.js
// Central configuration for the cwd_max workflow.
// Edit SITES, SPINUP_YEAR, ANALYSIS_YEARS, and REFERENCE_YEARS between runs.


// Study site geometries.
//
// Defined in EPSG:5070 (not WGS84) so that the AOI and the export region
// are the same rectangle.

// geodesic = false is required: without it GEE treats the edges as great
// circles instead of straight lines in the projected CRS.

exports.SITES = {
  northern_california: ee.Geometry.Rectangle(
    [-2337810, 2089020, -2060790, 2366130], 'EPSG:5070', false
  ),

  edwards_plateau: ee.Geometry.Rectangle(
    [-418230, 635220, -19830, 1074360], 'EPSG:5070', false
  ),

  high_plains: ee.Geometry.Rectangle(
    [-687480, 1476120, -115890, 1979370], 'EPSG:5070', false
  ),

  appalachia: ee.Geometry.Rectangle(
    [994740, 1635390, 1519050, 2128530], 'EPSG:5070', false
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