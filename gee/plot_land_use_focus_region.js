// plot_land_use_focus_region.js
// Visualisation of NLCD 2021 land cover over one study region (Great Plains).
// Data: USGS/NLCD_RELEASES/2021_REL/NLCD

// Study area
var geometry = ee.Geometry.Polygon(
    [[[-103.73661176242638, 40.588750552430916],
      [-103.73661176242638, 36.326001518952935],
      [ -97.38102094211388, 36.326001518952935],
      [ -97.38102094211388, 40.588750552430916]]], null, false);

// NLCD collection
var dataset = ee.ImageCollection('USGS/NLCD_RELEASES/2021_REL/NLCD');
print('Products:', dataset.aggregate_array('system:index'));

// Select the 2021 product and its land cover band
var nlcd2021  = dataset.filter(ee.Filter.eq('system:index', '2021')).first();
print('Bands:', nlcd2021.bandNames());
var landcover = nlcd2021.select('landcover');

// Display
Map.setCenter(-95, 38, 5);
Map.addLayer(landcover, null, 'Landcover');

// Area
print('Area (km2):', geometry.area(1).divide(1e6));