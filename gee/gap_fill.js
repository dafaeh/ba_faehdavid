// gap_fill.js 
// Monthly ET gap-filling for DisALEXI/OpenET.
// Missing pixels (cloud cover, SLC-off gaps) are filled with the
// multi-year reference mean for the same calendar month.

var cfg = require('users/dafaeh/ba_cwd:config');

// Returns a 12-image collection (one per calendar month) of mean ET
// across all reference years. Built once per AOI and reused for every
// analysis year.
// Pixels with no data across all reference years remain masked. NA
// propagates through to cwd_max and gap_filled_months downstream.
exports.build_ref_mean_et = function(aoi, years) {
  var collections = years.map(function(yr) {
    return ee.ImageCollection(cfg.OPENET_ASSET)
      .filterBounds(aoi)
      .filterDate(
        ee.Date.fromYMD(yr, 1, 1),
        ee.Date.fromYMD(yr + 1, 1, 1)
      )
      .select('et')
      .map(function(img) {
        return img.set('month', img.date().get('month'));
      });
  });
  var all_scenes = collections.reduce(function(acc, col) {
    return acc.merge(col);
  });

  return ee.ImageCollection(ee.List.sequence(1, 12).map(function(m) {
    return all_scenes
      .filter(ee.Filter.eq('month', m))
      .mean()
      .set('month', m);
  }));
};

// Returns a gap-filled ET image plus a gap-fill flag for one (year, month).
// gap_filled band: 0 = raw valid, 1 = filled from reference mean, NA = no data.
exports.get_monthly_et = function(yr, month, aoi, ref_mean_et) {
  month  = ee.Number(month);
  var t0 = ee.Date.fromYMD(yr, month, 1);
  var t1 = t0.advance(1, 'month');

  var et_raw = ee.ImageCollection(cfg.OPENET_ASSET)
    .filterBounds(aoi)
    .filterDate(t0, t1)
    .select('et')
    .mean();

  var et_ref = ee.Image(ref_mean_et
    .filter(ee.Filter.eq('month', month))
    .first());

  // Explicit validity masks to avoid silent NA→0 from .unmask() propagation.
  var raw_valid = et_raw.mask().unmask(0);
  var ref_valid = et_ref.mask().unmask(0);
  var any_valid = raw_valid.or(ref_valid);

  var et_filled = et_raw
    .unmask(et_ref)
    .updateMask(any_valid)
    .rename('et');

  var gap_filled = raw_valid.not()
    .updateMask(any_valid)
    .toByte()
    .rename('gap_filled');

  return et_filled.addBands(gap_filled);
};