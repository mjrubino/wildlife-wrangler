#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 15 11:03:56 2019

@author: nmtarr

Description: Use occurrence polygons to evaluate GAP range maps.

TO DO:
1.  max_error_meters -> error_tolerance
2  remove pad?
"""
#############################################################################
#                               Configuration
#############################################################################
sp_id = 'bybcux0'
summary_name = 'cuckoo'
gbif_req_id = 'r001'
gbif_filter_id = 'f001'

workDir = '/Users/nmtarr/Documents/RANGES/'
codeDir = '/Users/nmtarr/Code/Ranger/'
inDir = workDir + 'Inputs/'
outDir = workDir + 'Outputs/'
# Used in file names for output.
SRID_dict = {'WGS84': 4326, 'AlbersNAD83': 102008}


#############################################################################
#                                  Imports
#############################################################################
import pandas as pd
pd.set_option('display.width', 1000)
#%matplotlib inline
import sqlite3
import sciencebasepy
from pygbif import occurrences
import os
os.chdir('/')
os.chdir(codeDir)
import config
import config
import sqlite3
import os


#############################################################################
#                              Species-concept
#############################################################################
os.chdir(codeDir)
# Get species info from requests database
conn2 = sqlite3.connect(inDir + 'requests.sqlite')
cursor2 = conn2.cursor()
sql_tax = """SELECT gbif_id, common_name, scientific_name,
                    error_tolerance, gap_id, pad
             FROM species_concepts
             WHERE species_id = '{0}';""".format(sp_id)
concept = cursor2.execute(sql_tax).fetchall()[0]
gbif_id = concept[0]
common_name = concept[1]
scientific_name = concept[2]
error_toler = concept[3]
gap_id = concept[4]
pad = concept[5]


#############################################################################
#                          Connect to Database
#############################################################################
# Delete the database if it already exists
evdb = outDir + 'range_eval.sqlite'
if os.path.exists(spdb):
    os.remove(spdb)

# Create or connect to the database
conn = sqlite3.connect(evdb)
os.putenv('SPATIALITE_SECURITY', 'relaxed')
conn.enable_load_extension(True)
conn.execute('SELECT load_extension("mod_spatialite")')
cursor = conn.cursor()

# Make db spatial
cursor.execute('SELECT InitSpatialMetadata();')

sql_rngy = """
        /* Make a table for storing range maps for unique species-time period
           combinations, WITH GEOMETRY */
        CREATE TABLE IF NOT EXISTS range_polygons (
                     rng_polygon_id TEXT NOT NULL PRIMARY KEY,
                     alias TEXT UNIQUE,
                     species_id TEXT NOT NULL,
                     months TEXT,
                     years TEXT,
                     method TEXT,
                     max_error_meters INTEGER,
                     pad INTEGER,
                     date_created TEXT
                     );

        SELECT AddGeometryColumn('range_polygons', 'range_4326', 4326,
                                 'MULTIPOLYGON', 'XY');

        SELECT AddGeometryColumn('range_polygons', 'occurrences_4326', 4326,
                                 'MULTIPOLYGON', 'XY');
"""
cursor.executescript(sql_rngy)


#############################################################################
#                          Make Some Range Polygons
#############################################################################
# Function for making range_polygons
def MakeConcaveHull(rng_poly_id, alias, sp_id, months, years, workDir):
    '''
    Function for creating a range polygon entry in range_eval.range_polygons.
    
    Arguments:
    rng_poly_id -- A unique ID to use for the range map record
    alias -- keyword to use for filenames and shorthand reference to polygon
    sp_id -- species id for this project.  Must be in requests.species_concepts.
    months -- tuple of months to include.  For example: (3,4,5,6,7)
    years -- range of years to use.  Format as '1980-2000'
    
    '''
    print('SRID being used is 4326')
    sql = """
    * Attach requests database */
    ATTACH DATABASE '/Users/nmtarr/Documents/RANGES/Inputs/requests.sqlite'
    AS requests;
    
    /* Attach an occurrences database */
    ATTACH DATABASE '/Users/nmtarr/Documents/RANGES/Outputs/sp_id_occurrences.sqlite'
    AS occs;
    
    /* Create range map for the period. */
    INSERT INTO range_polygons (rng_poly_id, alias, species_id, 
                                months, years,
                                method, date_created, range_4326, 
                                occurrences_4326)
                    SELECT '{0}', '{1}', '{2}', '{3}', '{4},
                        concave hull', date('now'),
                        ConcaveHull(CastToMultiPolygon(GUnion(O.circle_wgs84))),
                        CastToMultiPolygon(GUnion(O.circle_wgs84))
                    FROM occs.occurrences AS O
                    WHERE cast(strftime('%m', occurrenceDate) AS INT) IN {3};
    
    /* Update the range tolerance and pad information */
    UPDATE range_polygons
    SET max_error_meters = (SELECT error_tolerance 
                            FROM requests.species_concepts 
                            WHERE species_id = 'sp_id'),
        pad = (SELECT pad 
               FROM requests.species_concepts 
               WHERE species_id = 'sp_id')
    WHERE species_id = 'sp_id';
    
    /* Recover geometry */
    SELECT RecoverGeometryColumn('range_polygons', 'range_4326', 4326, 'MULTIPOLYGON',
                               'XY');
    
    SELECT RecoverGeometryColumn('range_polygons', 'occurrences_4326', 4326, 'MULTIPOLYGON',
                               'XY');
    
    
    /* Pull out the period for mapping */
    CREATE TABLE temp1 AS SELECT * FROM range_polygons
                    WHERE  alias = '{1}';
    
    SELECT RecoverGeometryColumn('temp1', 'range_4326', 4326, 'MULTIPOLYGON',
                                 'XY');
    
    SELECT RecoverGeometryColumn('temp1', 'occurrences_4326', 4326, 'MULTIPOLYGON',
                                 'XY');
    
    /* Export shapefiles */
    SELECT ExportSHP('temp1', 'range_4326',
                     '{5}{1}_range', 'utf-8');
    
    SELECT ExportSHP('temp1', 'occurrences_4326',
                     '{5}{1}_occs', 'utf-8');
    
    DROP TABLE temp1;
    
    DETACH DATABASE requests;
    
    DETACH DATABASE occs;
    """.format(rng_poly_id, alias, sp_id, months, years, outDir)
    
    return 
    

# Make occurrence shapefiles for each month
month_dict = {'january': (1), 'february':(2), 'march':(3), 'april':(4), 
              'may':(5), 'june':(6), 'july':(7), 'august':(8), 
              'september':(9), 'october':(10), 'november':(11), 
              'december':(12)}

MakeConcaveHull(rng_poly_id=sp_id, alias='month', sp_id=sp_id, months=(6), 
                years='1980-2018', workDir=workDir)

for month in month_dict.keys():
    print(month)
    try:
        sql4 = """
        /* Attach an occurrences database

        /* Create tables for each month and export as shapefiles. */
        INSERT INTO range_polygons (species_id, period, range, circles)
                        SELECT species_id, '{0}',
                        ConcaveHull(CastToMultiPolygon(GUnion(circle_albers))),
                        CastToMultiPolygon(GUnion(circle_albers))
                        FROM occs
                        WHERE occurrenceMonth = {1};

        /* Pull out the period for mapping */
        CREATE TABLE temp1 AS SELECT * FROM rangemaps
                        WHERE period='{0}';

        SELECT RecoverGeometryColumn('temp1', 'range', 102008, 'MULTIPOLYGON',
                                     'XY');

        SELECT RecoverGeometryColumn('temp1', 'circles', 102008, 'MULTIPOLYGON',
                                     'XY');

        /* Transform back to WGS84 so that map can be displayed in ipython */
        ALTER TABLE temp1 ADD COLUMN range_wgs84 blob;

        SELECT RecoverGeometryColumn('temp1', 'range_wgs84', 4326,
                                     'MULTIPOLYGON');

        UPDATE temp1 SET range_wgs84 = Transform(range, 4326);

        ALTER TABLE temp1 ADD COLUMN circles_wgs84 blob;

        SELECT RecoverGeometryColumn('temp1', 'circles_wgs84', 4326,
                             'MULTIPOLYGON');

        UPDATE temp1 SET circles_4326 = Transform(circles, 4326);

        /* Export shapefiles */
        SELECT ExportSHP('temp1', 'range_wgs84', '{0}_rng', 'utf-8');

        SELECT ExportSHP('temp1', 'circles_wgs84', '{0}_occs', 'utf-8');

        DROP TABLE temp1;

        """.format(month, month_dict[month])
        cursor.executescript(sql4)
    except:
        print(Exception)

# Make range shapefiles for each season, display them too
period_dict = {"summer": '(5,6,7,8)', "winter": '(11,12,1,2)',
               "spring": '(3,4,5)', "fall": '(8,9,10,11)',
               "yearly": '(1,2,3,4,5,6,7,8,9,10,11,12)'}
for period in period_dict:
    print(period)
    try:
        sql_season = """
            /*  Insert a record for a range map, created by making polygons into
            a multipolygon geometry and then calculating the concave hull.
            Also, insert circles into a column for provenance */
            INSERT INTO rangemaps (species_id, period, range, circles)
                    SELECT species_id, '{0}',
                    ConcaveHull(CastToMultiPolygon(GUnion(circle_albers))),
                    CastToMultiPolygon(GUnion(circle_albers))
                    FROM occs
                    WHERE occurrenceMonth IN {1};

            /* Pull out the period for mapping */
            CREATE TABLE temp2 AS SELECT * FROM rangemaps
                            WHERE period='{0}';

            SELECT RecoverGeometryColumn('temp2', 'range', 102008,
                                         'MULTIPOLYGON',
                                         'XY');

            SELECT RecoverGeometryColumn('temp2', 'circles', 102008,
                                         'MULTIPOLYGON',
                                         'XY');

            /* Transform back to WGS84 so that map can be displayed in ipython */
            ALTER TABLE temp2 ADD COLUMN range_wgs84 blob;

            SELECT RecoverGeometryColumn('temp2', 'range_wgs84', 4326,
                                         'MULTIPOLYGON');

            UPDATE temp2 SET range_wgs84 = Transform(range, 4326);

            ALTER TABLE temp2 ADD COLUMN circles_wgs84 blob;

            SELECT RecoverGeometryColumn('temp2', 'circles_wgs84', 4326,
                                         'MULTIPOLYGON');

            UPDATE temp2 SET circles_wgs84 = Transform(circles, 4326);


            SELECT ExportSHP('temp2', 'range_wgs84', '{0}_rng', 'utf-8');

            SELECT ExportSHP('temp2', 'circles_wgs84', '{0}_occs', 'utf-8');

            DROP TABLE temp2;
        """.format(period, period_dict[period])
        cursor.executescript(sql_season)
    except:
        print(Exception)
conn.commit()


###############################################  DISPLAY MAPS
#############################################################

##########################################################

workDir = '/Users/nmtarr/Documents/RANGES'



season_colors = {'Fall': 'red', 'Winter': 'white', 'Summer': 'magenta',
                    'Spring': 'blue'}
for period in ['Fall', 'Winter', 'Summer', 'Spring']:
     shp1 = {'file': workDir + '/{0}_rng'.format(period),
                    'drawbounds': True, 'linewidth': 1,
                    'linecolor': season_colors[period],
                    'fillcolor': None}
     shp2 = {'file': workDir + '/{0}_occs'.format(period),
                    'drawbounds': True, 'linewidth': .5, 'linecolor': 'k',
                    'fillcolor': None}
     title = "Yellow-billed Cuckoo occurrence polygons - {0}".format(period)
     try:
         config.MapPolygonsFromSHP([shp1, shp2], title)
     except:
         print(period + " FAILED !!!!")

###########################################  GET THE GAP RANGES
#############################################  FROM SCIENCEBASE
