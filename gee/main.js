// main.js  
// Orchestrates the cwd_max workflow across all sites and years.
// Per site: builds the ET reference mean once, runs spinup + analysis years
// sequentially (snow pool carried forward), then exports cwd_max and
// gap_filled_months as a two-band GeoTIFF to Google Drive.
// Export region: aoi.bounds() in EPSG:5070 to avoid rotated NoData footprints.

var cfg         = require('users/dafaeh/ba_cwd:config');
var gap_fill    = require('users/dafaeh/ba_cwd:gap_fill');
var cwd_compute = require('users/dafaeh/ba_cwd:cwd_compute');

// Submits one Google Drive export task for a study site.
// aoi.bounds() returns the axis-aligned EPSG:5070 bounding box of the AOI,
// avoiding diagonal NoData edges from reprojecting WGS84 polygons.
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

Object.keys(cfg.SITES).forEach(function(site_name) {
  var aoi = cfg.SITES[site_name];
  print('Processing site:', site_name);

  // Monthly ET reference mean over REFERENCE_YEARS: built once per site
  // and reused across every analysis year.
  var ref_mean_et = gap_fill.build_ref_mean_et(aoi, cfg.REFERENCE_YEARS);

  // Sequential multi-year loop: spinup year first to initialise the snow
  // pool, then each analysis year with the pool carried forward.
  var snow_pool  = ee.Image.constant(0).rename('snow_pool');
  var cwd_images = [];
  var gap_images = [];

  var all_years = [cfg.SPINUP_YEAR].concat(cfg.ANALYSIS_YEARS);

  all_years.forEach(function(yr) {
    print('  Year:', yr);

    var out = cwd_compute.compute_year(yr, aoi, ref_mean_et, snow_pool);
    snow_pool = out.final_snow_pool;

    // Discard spinup outputs.
    if (yr !== cfg.SPINUP_YEAR) {
      cwd_images.push(out.cwd_max);
      gap_images.push(out.gap_filled_months);
    }
  });

// Both bands cast to Float32 so masked pixels export as NaN.
  // Single two-band export avoids redundant recomputation of the shared
  // upstream graph (reference mean ET, gap-fill logic).
  var cwd_global_max  = ee.ImageCollection(cwd_images).max()  // pixel-wise maximum across analysis years
                          .toFloat()
                          .rename('cwd_max');
  var gap_total_count = ee.ImageCollection(gap_images).sum()  // pixel-wise sum across analysis years
                          .toFloat()
                          .rename('gap_filled_months');

  export_site(site_name, cwd_global_max, gap_total_count, aoi);
});

print('All export tasks submitted. Check the Tasks tab.');