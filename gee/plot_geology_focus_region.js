// plot_geology_focus_region.js
// Visualisation of surface lithology and elevation over two study regions.
// Data: CSP/ERGo/1_0/US/lithology, USGS/SRTMGL1_003

// Study areas
var geometry = ee.Geometry.Polygon(
    [[[-100.27776465680033, 32.660183285641914],
      [-100.27776465680033, 28.789320221427214],
      [ -96.21282325055033, 28.789320221427214],
      [ -96.21282325055033, 32.660183285641914]]], null, false);

var geometry2 = ee.Geometry.Polygon(
    [[[-123.74888168761035, 41.29344890598067],
      [-123.74888168761035, 39.29196540299779],
      [-121.09019028136035, 39.29196540299779],
      [-121.09019028136035, 41.29344890598067]]], null, false);

// Lithology
var lithology = ee.Image('CSP/ERGo/1_0/US/lithology').select('b1');

var lithologyVis = {
  min: 0.0,
  max: 20.0,
  palette: [
    '356eff', 'acb6da', 'd6b879', '313131', 'eda800', '616161', 'd6d6d6',
    'd0ddae', 'b8d279', 'd5d378', '141414', '6db155', '9b6d55', 'feeec9',
    'd6b879', '00b7ec', 'ffda90', 'f8b28c'
  ]
};

Map.addLayer(lithology, lithologyVis, 'Lithology');

// Elevation (SRTM)
var dem = ee.Image('USGS/SRTMGL1_003');
Map.addLayer(dem, { min: 0, max: 2000, palette: ['000000', 'ffffff'] }, 'Elevation');

// Areas
print('Area geometry  (km2):', geometry.area(1).divide(1e6));
print('Area geometry2 (km2):', geometry2.area(1).divide(1e6));