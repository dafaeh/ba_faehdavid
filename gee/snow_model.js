// snow_model.js 
// Monthly snow accumulation and melt (Orth et al. 2013).
// Partitioning: precip falls as snow if tmean ≤ TEMP_THRESHOLD, else rain.
// Melt: min(snow_pool, MAX_MELT_RATE × days × (T − T0)) when tmean > T0.
// Liquid input to soil (p_in): rain + melt.

var cfg      = require('users/dafaeh/ba_cwd:config');
var gap_fill = require('users/dafaeh/ba_cwd:gap_fill');

// Single monthly update of the snow pool.
// Returns new snow pool and liquid water (rain + melt) reaching the soil.
// Shared by cwd_compute and the spinup loop in main.js.
exports.snow_step = function(snow_pool, img) {
  var tmean         = img.select('tmean');
  var rain          = img.select('rain');
  var snow          = img.select('snow');
  var days_in_month = ee.Number(img.get('days_in_month'));

  var melt_potential = tmean.subtract(cfg.TEMP_THRESHOLD)
                            .multiply(days_in_month)
                            .multiply(cfg.MAX_MELT_RATE)
                            .max(ee.Image.constant(0));
  var melt           = snow_pool.min(melt_potential);
  var new_snow_pool  = snow_pool.add(snow).subtract(melt)
                                .max(ee.Image.constant(0))
                                .rename('snow_pool');
  var p_in           = rain.add(melt).rename('p_in');

  return { snow_pool: new_snow_pool, p_in: p_in };
};

// Assembles a 12-image monthly collection for a given year.
// Each image contains: et, gap_filled, rain, snow, tmean, days_in_month.
exports.build_monthly_col = function(yr, aoi, ref_mean_et) {
  var daymet_col = ee.ImageCollection(cfg.DAYMET_ASSET)
    .filterBounds(aoi)
    .filterDate(
      ee.Date.fromYMD(yr, 1, 1),
      ee.Date.fromYMD(yr + 1, 1, 1)
    )
    .select(['prcp', 'tmax', 'tmin']);

  var months = ee.List.sequence(1, 12);
  return ee.ImageCollection(months.map(function(m) {
    m = ee.Number(m);
    var t0 = ee.Date.fromYMD(yr, m, 1);
    var t1 = t0.advance(1, 'month');
    var days_in_month = t1.difference(t0, 'day');

    // Gap-filled ET for this month (uses the precomputed reference mean).
    var et_bands = gap_fill.get_monthly_et(yr, m, aoi, ref_mean_et);

    var daymet = daymet_col.filterDate(t0, t1);
    var prcp   = daymet.select('prcp').sum().rename('prcp');
    var tmean  = daymet
      .map(function(img) {
        return img.select('tmax').add(img.select('tmin')).divide(2);
      })
      .mean()
      .rename('tmean');

    var is_snow = tmean.lte(cfg.TEMP_THRESHOLD);
    var snow    = prcp.multiply(is_snow).rename('snow');
    var rain    = prcp.subtract(snow).rename('rain');

    return et_bands
      .addBands(rain)
      .addBands(snow)
      .addBands(tmean)
      .set('system:time_start', t0.millis())
      .set('month', m)
      .set('days_in_month', days_in_month);
  }));
};
