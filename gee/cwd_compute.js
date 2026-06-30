// cwd_compute.js  
// Computes cwd_max for a single year.
// CWD accumulates monthly water deficit: cwd_m = max(cwd_{m-1} + et_m − p_in_m, 0)
// cwd_max is the peak value across all years

var snow_model = require('users/dafaeh/ba_cwd:snow_model');

// Iterates through 12 months, updating the snow pool and accumulating CWD.
// Returns final_snow_pool so main.js can carry it into the next year.
// Returns gap_filled_months as a per-pixel count of gap-filled months (0–12).
exports.compute_year = function(yr, aoi, ref_mean_et, init_snow_pool) {
  var monthly_col  = snow_model.build_monthly_col(yr, aoi, ref_mean_et);
  var monthly_list = monthly_col.toList(12);

  // Accumulator: [snow_pool, [cwd_month_0, cwd_month_1, ...]]
  // Starts with a zero image representing the pre-season state.
  var init = ee.List([
    init_snow_pool,
    ee.List([ee.Image.constant(0).rename('cwd')])
  ]);

  var result = ee.List.sequence(0, 11).iterate(function(i, acc) {
    var state     = ee.List(acc);
    var snow_pool = ee.Image(state.get(0));
    var cwd_list  = ee.List(state.get(1));
    var prev_cwd  = ee.Image(cwd_list.get(-1));

    var img       = ee.Image(monthly_list.get(ee.Number(i)));
    var et        = img.select('et');
    var step      = snow_model.snow_step(snow_pool, img);

    var new_cwd   = prev_cwd.add(et.subtract(step.p_in))
                            .max(ee.Image.constant(0))
                            .rename('cwd');

    return ee.List([step.snow_pool, cwd_list.add(new_cwd)]);
  }, init);

  var final_state     = ee.List(result);
  var final_snow_pool = ee.Image(final_state.get(0));

  // Drop the initial zero image before taking the maximum.
  // No .clip(aoi): clipping to the WGS84 polygon would produce diagonal
  // NoData edges in EPSG:5070. Spatial extent is controlled at export time
  // via aoi.bounds() in main.js.
  var cwd_images = ee.List(final_state.get(1)).slice(1);
  var cwd_max    = ee.ImageCollection(cwd_images).max()
                     .rename('cwd_max');

  // gap_filled_months shares cwd_max's mask: the count is meaningless
  // where cwd_max is NA (a single masked month taints the whole year).
  var gap_filled_months = monthly_col.select('gap_filled')
                                     .sum()
                                     .updateMask(cwd_max.mask())
                                     .rename('gap_filled_months');

  return {
    cwd_max:           cwd_max,
    gap_filled_months: gap_filled_months,
    final_snow_pool:   final_snow_pool
  };
};
