// plot_topography_focus_region.js
// Visualisation of elevation over two focus regions (Appalachians and Northern California).
// Data: USGS/SRTMGL1_003

// Focus regions
var geometry = ee.Geometry.Polygon(
    [[[-84.05115497183012, 40.97746916994725],
      [-84.05115497183012, 37.16650782156939],
      [-78.63489520620512, 37.16650782156939],
      [-78.63489520620512, 40.97746916994725]]], null, false);

var geometry2 = ee.Geometry.Polygon(
    [[[-123.74888168761035, 41.29344890598067],
      [-123.74888168761035, 39.29196540299779],
      [-121.09019028136035, 39.29196540299779],
      [-121.09019028136035, 41.29344890598067]]], null, false);

// Elevation (SRTM)
var dem = ee.Image('USGS/SRTMGL1_003');
Map.addLayer(dem, { min: 0, max: 2000, palette: ['000000', 'ffffff'] }, 'Elevation');

// Focus region outlines
Map.addLayer(geometry,  { color: 'bf04c2' }, 'Appalachians');
Map.addLayer(geometry2, { color: '00ff00' }, 'N. California');

// Areas
print('Area geometry  (km2):', geometry.area(1).divide(1e6));
print('Area geometry2 (km2):', geometry2.area(1).divide(1e6));