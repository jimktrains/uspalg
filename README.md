# uspalg - micro spatial algorithms

This is a collection of spatial, time, and astronomical calculations
intented for use on a microconroller.

Datafiles via

    wget -x https://data.iana.org/time-zones/releases/tzdata2023c.tar.gz
    wget -x https://data.iana.org/time-zones/releases/tzdb-2023c.tar.lz
    wget -x get https://github.com/evansiroky/timezone-boundary-builder/releases/download/2023b/timezones-with-oceans.shapefile.zip

Converted to a gmt file via

    ogr2ogr combined-shapefile-with-oceans.gmt combined-shapefile-with-oceans.shp

I would just read the shapefile directly, but there's a dependency conflict
with libgdal-dev on my machine and I got fed up dealing with it.
