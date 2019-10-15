shapefile = '/Users/nmtarr/Temp/geomtest.shp'
out_txt = '/Users/nmtarr/Temp/wkt.txt'
#def ccw_wkt_from_shp(shapefile, out_txt):
"""
(str, str) = text written to shpefile

Creates wkt with coordinates oriented counter clockwise for a given shapefile.
Shapefiles are oriented clockwise, which is incompatible with spatial queries
in many database management systems.  Use this to generate wkt that you can
copy and paste into queries.

Arguments:
shapefile -- path to the shpefile to read.
out_txt -- path to the text file to write the wkt to.
"""

import fiona
import shapely
from shapely.geometry import shape, Polygon, LinearRing
#from shapely.wkb import dumps, loads

# Read in a shapefile of polygon of interest.  It must be in CRS 4326
# First get a fiona collection
c = fiona.open(shapefile, 'r')

# Next make it a shapely polygon object
poly = shape(c[0]['geometry'])

# Use LinearRing to determine if coordinates are listed clockwise
coords = c[0]["geometry"]["coordinates"][0]
lr = LinearRing(coords)
if lr.is_ccw == False:
    # Reverse coordinates to make them counter clockwise
    print("Points are clockwise")
    #coords.reverse()
    # Make the polygon's outer ring counter clockwise
    poly2 = shapely.geometry.polygon.orient(poly, sign=1.0)
    # Get the well-known text version of the polygon
    wkt = poly2.wkt
else:
    print("Points were already counter clockwise")
    # Get the well-known text version of the polygon
    wkt = poly.wkt

# Write WKT to text file
print(wkt)
with open(out_txt, 'w') as file:
    file.write(wkt)

# close the collections
c.close()
