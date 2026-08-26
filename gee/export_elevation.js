// export_terrain.js
// Exports FABDEM elevation for all sites defined in config.js.

var cfg = require('users/dafaeh/ba_cwd:config');

var fabdem = ee.ImageCollection('projects/sat-io/open-datasets/FABDEM')
  .select('b1')
  .mosaic()
  .rename('elevation');

Object.keys(cfg.SITES).forEach(function(site_name) {
  var aoi = cfg.SITES[site_name];

  var export_region = aoi.bounds(
    ee.ErrorMargin(1, 'meters'),
    cfg.EXPORT_CRS
  );

  Export.image.toDrive({
    image:          fabdem,
    description:    'elevation_' + site_name,
    folder:         cfg.EXPORT_FOLDER,
    fileNamePrefix: 'elevation_' + site_name,
    region:         export_region,
    scale:          cfg.EXPORT_SCALE,
    crs:            cfg.EXPORT_CRS,
    maxPixels:      cfg.MAX_PIXELS
  });

  print('Export submitted: elevation_' + site_name);
});

print('All elevation exports submitted. Check the Tasks tab.');