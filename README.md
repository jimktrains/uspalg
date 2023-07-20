# uspalg

TZ Shapefile from https://github.com/evansiroky/timezone-boundary-builder
Converted to a gmt file via

    ogr2ogr combined-shapefile-with-oceans.gmt combined-shapefile-with-oceans.shp

I would just read the shapefile directly, but there's a dependency conflict
with libgdal-dev on my machine and I got fed up dealing with it.
