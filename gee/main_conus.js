// main_conus.js
// Orchestrates the CONUS-wide cwd_max workflow.
//
// Identical in structure to main.js except:
//   - Loads config_conus (10x8 tile grid, broader REFERENCE_YEARS).
//   - Iterates over all tiles defined in cfg.SITES.
//
// To run a subset of tiles, comment out the unwanted tiles in
// config_conus.js before submitting tasks.

var cfg         = require('users/dafaeh/ba_cwd:config_conus');
var gap_fill    = require('users/dafaeh/ba_cwd:gap_fill');
var cwd_compute = require('users/dafaeh/ba_cwd:cwd_compute');

function export_site(site_name, cwd_max, gap_filled_months, aoi) {
  var combined = cwd_max.addBands(gap_filled_months);
  var export_region = aoi.bounds(
    ee.ErrorMargin(1, 'meters'),
    cfg.EXPORT_CRS
  );
  Export.image.toDrive({
    image:          combined,
    description:    'cwd_' + site_name,
    folder:         cfg.EXPORT_FOLDER,
    fileNamePrefix: 'cwd_' + site_name,
    region:         export_region,
    scale:          cfg.EXPORT_SCALE,
    crs:            cfg.EXPORT_CRS,
    maxPixels:      cfg.MAX_PIXELS
  });
}

var tile_keys = Object.keys(cfg.SITES);

if (tile_keys.length === 0) {
  print('WARNING: SITES is empty -- no tasks will be submitted.');
  print('Check config_conus.js.');
}

tile_keys.forEach(function(site_name) {
  var aoi = cfg.SITES[site_name];
  print('Processing tile:', site_name);

  var ref_mean_et = gap_fill.build_ref_mean_et(aoi, cfg.REFERENCE_YEARS);

  var snow_pool  = ee.Image.constant(0).rename('snow_pool');
  var cwd_images = [];
  var gap_images = [];
  var all_years  = [cfg.SPINUP_YEAR].concat(cfg.ANALYSIS_YEARS);

  all_years.forEach(function(yr) {
    print('  Year:', yr);
    var out = cwd_compute.compute_year(yr, aoi, ref_mean_et, snow_pool);
    snow_pool = out.final_snow_pool;
    if (yr !== cfg.SPINUP_YEAR) {
      cwd_images.push(out.cwd_max);
      gap_images.push(out.gap_filled_months);
    }
  });

  var cwd_global_max  = ee.ImageCollection(cwd_images).max()
                          .toFloat()
                          .rename('cwd_max');
  var gap_total_count = ee.ImageCollection(gap_images).sum()
                          .toFloat()
                          .rename('gap_filled_months');

  export_site(site_name, cwd_global_max, gap_total_count, aoi);
});

print('All', tile_keys.length, 'export tasks submitted. Check the Tasks tab.');