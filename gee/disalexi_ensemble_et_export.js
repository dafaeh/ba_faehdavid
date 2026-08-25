// export_et_comparison.js
// Exports raw mean monthly ET from OpenET Ensemble and DisALEXI for a
// single AOI, to compare DisALEXI against the ensemble mean.
// This is a diagnostic script, not part of the main cwd_max pipeline.

var cfg = require('users/dafaeh/ba_cwd:config');

// Ensemble is only needed here, so it stays out of config.js.
// v2_1 to match cfg.OPENET_ASSET.
var ENSEMBLE_ASSET =
  'projects/openet/assets/ensemble/conus/gridmet/monthly/v2_1';

// --- AOI ---
var LON_TRIM   = 0.0466;
var LON_EXTEND = 0.0233;

var aoi = ee.Geometry.Rectangle(
  [-123.7068926568897  + LON_TRIM,   39.42401674974655,
   -123.55171077212408 + LON_EXTEND, 39.57131505027446]
);
Map.centerObject(aoi, 11);

// --- Period: mean monthly ET over the cwd_max analysis years ---
// End date is excluded.
var startDate = '2020-01-01';
var endDate   = '2023-01-01';

// Mean across all 36 monthly images.
function mean_monthly_et(asset, band, new_name) {
  return ee.ImageCollection(asset)
    .filterDate(startDate, endDate)
    .filterBounds(aoi)
    .select(band)
    .mean()
    .rename(new_name);
}

var disalexi = mean_monthly_et(cfg.OPENET_ASSET, 'et', 'et_disalexi');
var ensemble = mean_monthly_et(ENSEMBLE_ASSET, 'et_ensemble_mad',
                               'et_ensemble');

var combined = disalexi.addBands(ensemble).toFloat();

// --- Percentiles across both products, to set the shared colour scale in R ---
var stats = combined.reduceRegion({
  reducer:   ee.Reducer.percentile([2, 98]),
  geometry:  aoi,
  scale:     cfg.EXPORT_SCALE,
  maxPixels: cfg.MAX_PIXELS
});
print('ET percentile 2/98 (both products):', stats);

// --- Map preview only ---
var etVis = {
  min: 80, max: 120,
  palette: [
    '#f7fbff', '#deebf7', '#c6dbef', '#9ecae1',
    '#6baed6', '#4292c6', '#2171b5', '#084594'
  ]
};
Map.addLayer(disalexi.clip(aoi), etVis, 'DisALEXI ET 2020-2022');
Map.addLayer(ensemble.clip(aoi), etVis, 'Ensemble ET 2020-2022');

// Masked pixels are transparent in the preview and easy to miss, so paint
// them explicitly before submitting the export.
Map.addLayer(disalexi.mask().not().selfMask().clip(aoi),
             {palette: 'red'}, 'DisALEXI NoData');

// --- Export raw ET as a two-band GeoTIFF to Google Drive ---
var export_region = aoi.bounds(
  ee.ErrorMargin(1, 'meters'),
  cfg.EXPORT_CRS
);

Export.image.toDrive({
  image:          combined,
  description:    'et_comparison_2020_2022',
  folder:         cfg.EXPORT_FOLDER,
  fileNamePrefix: 'et_comparison_2020_2022',
  region:         export_region,
  scale:          cfg.EXPORT_SCALE,
  crs:            cfg.EXPORT_CRS,
  maxPixels:      cfg.MAX_PIXELS
});

print('Export task submitted: et_comparison_2020_2022. Check the Tasks tab.');