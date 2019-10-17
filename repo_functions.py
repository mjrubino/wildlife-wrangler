# Define a function for displaying the maps that will be created.
def MapShapefilePolygons(map_these, title):
    """
    Displays shapefiles on a simple CONUS basemap.  Maps are plotted in the order
    provided so put the top map last in the listself.  You can specify a column
    to map as well as custom colors for it.  This function may not be very robust
    to other applications.

    NOTE: The shapefiles have to be in WGS84 CRS.

    (dict, str) -> displays maps, returns matplotlib.pyplot figure

    Arguments:
    map_these -- list of dictionaries for shapefiles you want to display in
                CONUS. Each dictionary should have the following format, but
                some are unneccesary if 'column' doesn't = 'None'.  The critical
                ones are file, column, and drawbounds.  Column_colors is needed
                if column isn't 'None'.  Others are needed if it is 'None'.
                    {'file': '/path/to/your/shapfile',
                     'alias': 'my layer'
                     'column': None,
                     'column_colors': {0: 'k', 1: 'r'}
                    'linecolor': 'k',
                    'fillcolor': 'k',
                    'linewidth': 1,
                    'drawbounds': True
                    'marker': 's'}
    title -- title for the map.
    """
    # Packages needed for plotting
    import matplotlib.pyplot as plt
    from mpl_toolkits.basemap import Basemap
    import numpy as np
    from matplotlib.patches import Polygon
    from matplotlib.collections import PatchCollection
    from matplotlib.patches import PathPatch

    # Basemap
    fig = plt.figure(figsize=(15,12))
    ax = plt.subplot(1,1,1)
    map = Basemap(projection='aea', resolution='l', lon_0=-95.5, lat_0=39.0,
                  height=3200000, width=5000000)
    map.drawcoastlines(color='grey')
    map.drawstates(color='grey')
    map.drawcountries(color='grey')
    map.fillcontinents(color='#a2d0a2',lake_color='#a9cfdc')
    map.drawmapboundary(fill_color='#a9cfdc')

    for mapfile in map_these:
        if mapfile['column'] == None:
            # Add shapefiles to the map
            if mapfile['fillcolor'] == None:
                map.readshapefile(mapfile['file'], 'mapfile',
                                  drawbounds=mapfile['drawbounds'],
                                  linewidth=mapfile['linewidth'],
                                  color=mapfile['linecolor'])
                # Empty scatter plot for the legend
                plt.scatter([], [], c='', edgecolor=mapfile['linecolor'],
                            alpha=1, label=mapfile['alias'], s=100,
                            marker=mapfile['marker'])

            else:
                map.readshapefile(mapfile['file'], 'mapfile',
                          drawbounds=mapfile['drawbounds'])
                # Code for extra formatting -- filling in polygons setting border
                # color
                patches = []
                for info, shape in zip(map.mapfile_info, map.mapfile):
                    patches.append(Polygon(np.array(shape), True))
                ax.add_collection(PatchCollection(patches,
                                                  facecolor= mapfile['fillcolor'],
                                                  edgecolor=mapfile['linecolor'],
                                                  linewidths=mapfile['linewidth'],
                                                  zorder=2))
                # Empty scatter plot for the legend
                plt.scatter([], [], c=mapfile['fillcolor'],
                            edgecolors=mapfile['linecolor'],
                            alpha=1, label=mapfile['alias'], s=100,
                            marker=mapfile['marker'])

        else:
            map.readshapefile(mapfile['file'], 'mapfile', drawbounds=mapfile['drawbounds'])
            for info, shape in zip(map.mapfile_info, map.mapfile):
                for thang in mapfile['column_colors'].keys():
                    if info[mapfile['column']] == thang:
                        x, y = zip(*shape)
                        map.plot(x, y, marker=None, color=mapfile['column_colors'][thang])

            # Empty scatter plot for the legend
            for seal in mapfile['column_colors'].keys():
                plt.scatter([], [], c=mapfile['column_colors'][seal],
                            edgecolors=mapfile['column_colors'][seal],
                            alpha=1, label=mapfile['value_alias'][seal],
                            s=100, marker=mapfile['marker'])

    # Legend -- the method that works is ridiculous but necessary; you have
    #           to add empty scatter plots with the symbology you want for
    #           each shapefile legend entry and then call the legend.  See
    #           plt.scatter(...) lines above.
    plt.legend(scatterpoints=1, frameon=True, labelspacing=1, loc='lower left',
               framealpha=1, fontsize='x-large')

    # Title
    plt.title(title, fontsize=20, pad=-40, backgroundcolor='w')
    return


# Define a function for displaying the maps that will be created.  North Carolina version.
def MapShapefilePolygons_NC(map_these, title):
    """
    Displays shapefiles on a simple CONUS basemap.  Maps are plotted in the order
    provided so put the top map last in the listself.  You can specify a column
    to map as well as custom colors for it.  This function may not be very robust
    to other applications.

    NOTE: The shapefiles have to be in WGS84 CRS.

    (dict, str) -> displays maps, returns matplotlib.pyplot figure

    Arguments:
    map_these -- list of dictionaries for shapefiles you want to display in
                CONUS. Each dictionary should have the following format, but
                some are unneccesary if 'column' doesn't = 'None'.  The critical
                ones are file, column, and drawbounds.  Column_colors is needed
                if column isn't 'None'.  Others are needed if it is 'None'.
                    {'file': '/path/to/your/shapfile',
                     'alias': 'my layer'
                     'column': None,
                     'column_colors': {0: 'k', 1: 'r'}
                    'linecolor': 'k',
                    'fillcolor': 'k',
                    'linewidth': 1,
                    'drawbounds': True
                    'marker': 's'}
    title -- title for the map.
    """
    # Packages needed for plotting
    import matplotlib.pyplot as plt
    from mpl_toolkits.basemap import Basemap
    import numpy as np
    from matplotlib.patches import Polygon
    from matplotlib.collections import PatchCollection
    from matplotlib.patches import PathPatch

    # Basemap
    fig = plt.figure(figsize=(15,12))
    ax = plt.subplot(1,1,1)
    map = Basemap(projection='aea', resolution='i', lon_0=-79.8, lat_0=35.5,
                  height=410000, width=900000)
    map.drawcoastlines(color='grey')
    map.drawstates(color='grey')
    map.drawcountries(color='grey')
    map.fillcontinents(color='#a2d0a2',lake_color='#a9cfdc')
    map.drawmapboundary(fill_color='#a9cfdc')

    for mapfile in map_these:
        if mapfile['column'] == None:
            # Add shapefiles to the map
            if mapfile['fillcolor'] == None:
                map.readshapefile(mapfile['file'], 'mapfile',
                                  drawbounds=mapfile['drawbounds'],
                                  linewidth=mapfile['linewidth'],
                                  color=mapfile['linecolor'])
                # Empty scatter plot for the legend
                plt.scatter([], [], c='', edgecolor=mapfile['linecolor'],
                            alpha=1, label=mapfile['alias'], s=100,
                            marker=mapfile['marker'])

            else:
                map.readshapefile(mapfile['file'], 'mapfile',
                          drawbounds=mapfile['drawbounds'])
                # Code for extra formatting -- filling in polygons setting border
                # color
                patches = []
                for info, shape in zip(map.mapfile_info, map.mapfile):
                    patches.append(Polygon(np.array(shape), True))
                ax.add_collection(PatchCollection(patches,
                                                  facecolor= mapfile['fillcolor'],
                                                  edgecolor=mapfile['linecolor'],
                                                  linewidths=mapfile['linewidth'],
                                                  zorder=2))
                # Empty scatter plot for the legend
                plt.scatter([], [], c=mapfile['fillcolor'],
                            edgecolors=mapfile['linecolor'],
                            alpha=1, label=mapfile['alias'], s=100,
                            marker=mapfile['marker'])

        else:
            map.readshapefile(mapfile['file'], 'mapfile', drawbounds=mapfile['drawbounds'])
            for info, shape in zip(map.mapfile_info, map.mapfile):
                for thang in mapfile['column_colors'].keys():
                    if info[mapfile['column']] == thang:
                        x, y = zip(*shape)
                        map.plot(x, y, marker=None, color=mapfile['column_colors'][thang])

            # Empty scatter plot for the legend
            for seal in mapfile['column_colors'].keys():
                plt.scatter([], [], c=mapfile['column_colors'][seal],
                            edgecolors=mapfile['column_colors'][seal],
                            alpha=1, label=mapfile['value_alias'][seal],
                            s=100, marker=mapfile['marker'])

    # Legend -- the method that works is ridiculous but necessary; you have
    #           to add empty scatter plots with the symbology you want for
    #           each shapefile legend entry and then call the legend.  See
    #           plt.scatter(...) lines above.
    plt.legend(scatterpoints=1, frameon=True, labelspacing=1, loc='lower left',
               framealpha=1, fontsize='x-large')

    # Title
    plt.title(title, fontsize=20, pad=-40, backgroundcolor='w')
    return

def download_GAP_range_CONUS2001v1(gap_id, toDir):
    """
    Downloads GAP Range CONUS 2001 v1 file and returns path to the unzipped
    file.  NOTE: doesn't include extension in returned path so that you can
    specify if you want csv or shp or xml when you use the path.
    """
    import sciencebasepy
    import zipfile

    # Connect
    sb = sciencebasepy.SbSession()

    # Search for gap range item in ScienceBase
    gap_id = gap_id[0] + gap_id[1:5].upper() + gap_id[5]
    item_search = '{0}_CONUS_2001v1 Range Map'.format(gap_id)
    items = sb.find_items_by_any_text(item_search)

    # Get a public item.  No need to log in.
    rng =  items['items'][0]['id']
    item_json = sb.get_item(rng)
    get_files = sb.get_item_files(item_json, toDir)

    # Unzip
    rng_zip = toDir + item_json['files'][0]['name']
    zip_ref = zipfile.ZipFile(rng_zip, 'r')
    zip_ref.extractall(toDir)
    zip_ref.close()

    # Return path to range file without extension
    return rng_zip.replace('.zip', '')


def getGBIFcode(name, rank='species'):
    """
    Returns the GBIF species code for a scientific name.

    Example: gbifcode = getGBIFcode(name = "Dendroica ceruleans")
    """
    from pygbif import species
    key = species.name_backbone(name = name, rank='species')['usageKey']
    return key



def retrieve_gbif_occurrences(codeDir, species_id, inDir, spdb, gbif_req_id,
                              gbif_filter_id, default_coordUncertainty, SRID_dict,
                              outDir, summary_name):
    """
    Retrieves GAP range from ScienceBase and occurrence records from APIs. Filters
    occurrence records, stores them in a database, buffers the xy points,
    and filtering occurrence records, saving them in a database.  Finally, exports
    some maps.

    To do:
    1.  Maximize filtering.
    2.  Can we use EPSG:5070?
    3.  Account for possiblity of non-4326 occurrence records in gbif? Solution
        would be to transform before entering into database.
    4. Make platform flexible

    Arguments:
    codeDir -- directory of this code repo.
    species_id -- project id for the species concept.
    inDir -- directory containing key inputs such as downloaded gap ranges.
    spdb -- occurrence record database to be created by this function.
    gbif_req_id -- GBIF request ID for the process.
    gbif_filter_id -- GBIF filter ID for the process.
    default_coordUncertainty -- distance in meters to use if no coordinate
        Uncertainty is specified for a record.
    SRID_dict -- a dictionary of spatial reference name-code pairs.
    outDir -- where to save maps that are exported by this process.
    summary_name -- a short name for some file names.
    """

    import pandas as pd
    pd.set_option('display.width', 1000)
    import sqlite3
    import sciencebasepy
    from pygbif import occurrences
    import os
    os.chdir('/')
    import json
    import platform
    import shapely
    from shapely.wkt import dumps, loads
    from datetime import datetime


    # Environment variables need to be handled
    if platform.system() == 'Windows':
        os.environ['PATH'] = os.environ['PATH'] + ';' + 'C:/Spatialite'
        os.environ['SPATIALITE_SECURITY'] = 'relaxed'# DOES THIS NEED TO BE RUN BEFORE EVERY CONNECTION????? ?NOT WORKING  ???????????

    if platform.system() == 'Darwin':  # DOES THIS NEED TO BE RUN BEFORE EVERY CONNECTION?????????????????
        #os.putenv('SPATIALITE_SECURITY', 'relaxed')
        os.environ['SPATIALITE_SECURITY'] = 'relaxed'

    print("SPATIALITE_SECURITY set to " + os.environ['SPATIALITE_SECURITY'])



    #############################################################################
    #                              Species-concept
    #############################################################################
    os.chdir(codeDir)
    # Get species info from requests database
    conn2 = sqlite3.connect(codeDir + 'parameters.sqlite')
    cursor2 = conn2.cursor()
    sql_tax = """SELECT gbif_id, common_name, scientific_name,
                        detection_distance_meters, gap_id, geometry
                 FROM species_concepts
                 WHERE species_id = '{0}';""".format(species_id)
    concept = cursor2.execute(sql_tax).fetchall()[0]
    gbif_id = concept[0]
    common_name = concept[1]
    scientific_name = concept[2]
    det_dist = concept[3]
    gap_id = concept[4]
    sp_geom =concept[5]


    #############################################################################
    #                           Create Occurrence Database
    #############################################################################
    """
    Description: Create a database for storing occurrence and species-concept
    data.  Needs to have spatial querying functionality.
    """
    makedb1 = datetime.now()
    spdb = spdb
    # Delete the database if it already exists
    if os.path.exists(spdb):
        os.remove(spdb)

    # Create or connect to the database
    conn = sqlite3.connect(spdb)
    conn.enable_load_extension(True)
    conn.execute('SELECT load_extension("mod_spatialite")')
    cursor = conn.cursor()

    # Make database spatial and add the spatial reference system that GAP used
    conn.executescript('''SELECT InitSpatialMetaData(1);

                         INSERT into spatial_ref_sys
                         (srid, auth_name, auth_srid, proj4text, srtext)
                         values (102008, 'ESRI', 102008, '+proj=aea +lat_1=20 +lat_2=60
                         +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m
                         +no_defs ', 'PROJCS["North_America_Albers_Equal_Area_Conic",
                         GEOGCS["GCS_North_American_1983",
                         DATUM["North_American_Datum_1983",
                         SPHEROID["GRS_1980",6378137,298.257222101]],
                         PRIMEM["Greenwich",0],UNIT["Degree",0.017453292519943295]],
                         PROJECTION["Albers_Conic_Equal_Area"],
                         PARAMETER["False_Easting",0],
                         PARAMETER["False_Northing",0],
                         PARAMETER["longitude_of_center",-96],
                         PARAMETER["Standard_Parallel_1",20],
                         PARAMETER["Standard_Parallel_2",60],
                         PARAMETER["latitude_of_center",40],
                         UNIT["Meter",1],AUTHORITY["EPSG","102008"]]');''')
    conn.commit()


    ################################################# Create tables
    ###############################################################
    sql_cdb = """
            /* Create a table for occurrence records, WITH GEOMETRY */
            CREATE TABLE IF NOT EXISTS occurrences (
                    occ_id INTEGER NOT NULL PRIMARY KEY,
                    species_id INTEGER NOT NULL,
                    source TEXT NOT NULL,
                    request_id TEXT NOT NULL,
                    filter_id TEXT NOT NULL,
                    latitude TEXT,
                    longitude TEXT,
                    coordinateUncertaintyInMeters INTEGER,
                    occurrenceDate TEXT,
                    retrievalDate TEXT NOT NULL DEFAULT CURRENT_TIMESTAMP,
                    individualCount INTEGER DEFAULT 1,
                    generalizations TEXT,
                    remarks TEXT,
                    detection_distance INTEGER,
                    radius_meters INTEGER,
                        FOREIGN KEY (species_id) REFERENCES taxa(species_id)
                        ON UPDATE RESTRICT
                        ON DELETE NO ACTION);

            SELECT AddGeometryColumn('occurrences', 'geom_xy4326', 4326, 'POINT',
                                     'XY');
    """
    cursor.executescript(sql_cdb)
    makedb2 = datetime.now()
    print("Created occurrence db: " + str(makedb2 - makedb1))

    #############################################################################
    #                              GBIF Records
    #############################################################################
    """
    Retrieve GBIF records for a species and save appropriate
    attributes in the occurrence db.
    """
    requesttime1 = datetime.now()
    ############################# RETRIEVE REQUEST PARAMETERS
    # Up-front filters are an opportunity to lighten the load from the start.
    sql_twi = """ SELECT lat_range FROM gbif_requests
                  WHERE request_id = '{0}'""".format(gbif_req_id)
    latRange = cursor2.execute(sql_twi).fetchone()[0]

    sql_twi = """ SELECT lon_range FROM gbif_requests
                  WHERE request_id = '{0}'""".format(gbif_req_id)
    lonRange = cursor2.execute(sql_twi).fetchone()[0]

    sql_twi = """ SELECT years_range FROM gbif_requests
                  WHERE request_id = '{0}'""".format(gbif_req_id)
    years = cursor2.execute(sql_twi).fetchone()[0]

    sql_twi = """ SELECT months_range FROM gbif_requests
                  WHERE request_id = '{0}'""".format(gbif_req_id)
    months = cursor2.execute(sql_twi).fetchone()[0]

    sql_twi = """ SELECT geoissue FROM gbif_requests
                  WHERE request_id = '{0}'""".format(gbif_req_id)
    geoIssue = cursor2.execute(sql_twi).fetchone()[0]
    if geoIssue == 'None':
        geoIssue = None

    sql_twi = """ SELECT coordinate FROM gbif_requests
                  WHERE request_id = '{0}'""".format(gbif_req_id)
    coordinate = cursor2.execute(sql_twi).fetchone()[0]

    sql_twi = """ SELECT continent FROM gbif_requests
                  WHERE request_id = '{0}'""".format(gbif_req_id)
    continent = cursor2.execute(sql_twi).fetchone()[0]
    if continent == "None":
        continent = None

    sql_twi = """ SELECT country FROM gbif_requests
                  WHERE request_id = '{0}'""".format(gbif_req_id)
    country = cursor2.execute(sql_twi).fetchone()[0]
    if country == "None":
        country = None

    ########### SORT OUT GEOMETRY FILTERS
    # Get the geometry from the request filter set
    sql_poly = """ SELECT geometry FROM gbif_requests
                  WHERE request_id = '{0}'""".format(gbif_req_id)
    poly0 = cursor2.execute(sql_poly).fetchone()[0]
    # A geometry could also be stated for the species, assess what to do
    if poly0 == None and sp_geom == None:
        poly = None
    elif poly0 != None and sp_geom == None:
        poly = poly0
    elif poly0 == None and sp_geom != None:
        poly = sp_geom
    elif poly0 != None and sp_geom != None:
        # Get/use the intersection of the two polygons
        filter_polygon = shapely.wkt.loads(poly0)
        sp_polygon = shapely.wkt.loads(sp_geom)
        poly_intersection = filter_polygon.intersection(sp_polygon)
        poly = shapely.wkt.dumps(poly_intersection)

    #################### REQUEST RECORDS ACCORDING TO REQUEST PARAMS
    # First, find out how many records there are that meet criteria
    occ_search = occurrences.search(gbif_id,
                                    year=years,
                                    month=months,
                                    decimalLatitude=latRange,
                                    decimalLongitude=lonRange,
                                    hasGeospatialIssue=geoIssue,
                                    hasCoordinate=coordinate,
                                    continent=continent,
                                    country=country,
                                    geometry=poly)
    occ_count=occ_search['count']

    # Get occurrences in batches, saving into master list
    alloccs = []
    batches = range(0, occ_count, 300)
    for i in batches:
        occ_json = occurrences.search(gbif_id,
                                      limit=300,
                                      offset=i,
                                      year=years,
                                      month=months,
                                      decimalLatitude=latRange,
                                      decimalLongitude=lonRange,
                                      hasGeospatialIssue=geoIssue,
                                      hasCoordinate=coordinate,
                                      continent=continent,
                                      country=country,
                                      geometry=poly)
        occs = occ_json['results']
        alloccs = alloccs + occs

    print("Downloaded records: " + str(datetime.now() - requesttime1))
    print('\t{0} records exist with the request parameters'.format(occ_count))

    ######################### CREATE SUMMARY TABLE OF KEYS/FIELDS RETURNED  !!!!!! SLOW PART STARTS
    requestsummarytime1 = datetime.now()
    keys = [list(x.keys()) for x in alloccs]
    keys2 = set([])
    for x in keys:
        keys2 = keys2 | set(x)
    dfK = pd.DataFrame(index=keys2, columns=['included(n)', 'populated(n)'])
    dfK['included(n)'] = 0
    dfK['populated(n)'] = 0
    for t in alloccs:
        for y in t.keys():
            dfK.loc[y, 'included(n)'] += 1
            try:
                int(t[y])
                dfK.loc[y, 'populated(n)'] += 1
            except:
                if t[y] == None:
                    pass
                elif len(t[y]) > 0:
                    dfK.loc[y, 'populated(n)'] += 1
    dfK.sort_index(inplace=True)
    dfK.to_sql(name='gbif_fields_returned', con=conn, if_exists='replace')

    ############################# SAVE SUMMARY OF VALUES RETURNED (REQUEST)
    summary = {'datums': ['WGS84'],
               'issues': set([]),
               'bases': [],
               'institutions': [],
               'collections': [],
               'generalizations': set([]),
               'remarks': set([]),
               'establishment': set([]),
               'IDqualifier': set([]),
               'protocols': set([])}

    value_summaries = {'bases': {},
                      'datums': {'WGS84': 0},
                      'issues': {},
                      'institutions': {},
                      'collections': {},
                      'protocols': {},
                      'samplingProtocols': {}}

    for occdict in alloccs:
        # datums
        if occdict['geodeticDatum'] != 'WGS84':
            summary['datums'] = summary['datums'] + occdict['geodeticDatum']
            if occdict['geodeticDatum'] not in value_summaries['datums'].keys():
                value_summaries['datums'][occdict['geodeticDatum']] = 0
            else:
                value_summaries['datums'][occdict['geodeticDatum']] += 1
        if occdict['geodeticDatum'] == 'WGS84':
            value_summaries['datums']['WGS84'] += 1

        # issues
        summary['issues'] = summary['issues'] | set(occdict['issues'])
        for issue in occdict['issues']:
            if issue not in value_summaries['issues'].keys():
                value_summaries['issues'][issue] = 1
            if issue in value_summaries['issues'].keys():
                value_summaries['issues'][issue] += 1

        # basis or record
        BOR = occdict['basisOfRecord']
        if BOR == "" or BOR == None:
            BOR = 'UNKNOWN'
        summary['bases'] = summary['bases'] + [BOR]

        if BOR in value_summaries['bases'].keys():
            value_summaries['bases'][BOR] += 1
        else:
            value_summaries['bases'][BOR] = 1

        # institution
        try:
            who = occdict['institutionID']
        except:
            try:
                who = occdict['institutionCode']
            except:
                who = 'unknown'

        summary['institutions'] = summary['institutions'] + [who]

        if who in value_summaries['institutions'].keys():
            value_summaries['institutions'][who] += 1
        else:
            value_summaries['institutions'][who] = 1

        # collections
        try:
            co = occdict['collectionCode']
        except:
            co = 'UNKNOWN'

        summary['collections'] = summary['collections'] + [co]

        if co in value_summaries['collections'].keys():
            value_summaries['collections'][co] += 1
        else:
            value_summaries['collections'][co] = 1

        # establishment means
        try:
            est = occdict['establishmentMeans']
            summary['establishment'] = summary['establishment'] | set([est])
        except:
            pass

        # identification qualifier
        try:
            qual = occdict['identificationQualifier']
            summary['IDqualifier'] = summary['IDqualifier'] | set([qual])
        except:
            pass

        # protocols -- NOTE this essentially combines two fields
        try:
            proto = occdict['protocol']
        except:
            proto = 'UNKNOWN'
        summary['protocols'] = summary['protocols'] | set([proto])

        if proto in value_summaries['protocols'].keys():
            value_summaries['protocols'][proto] += 1
        else:
            value_summaries['protocols'][proto] = 1

        try:
            samproto = occdict['samplingProtocol']
        except:
            samproto = 'UKNOWN'
        summary['protocols'] = summary['protocols'] | set([samproto])

        if samproto in value_summaries['samplingProtocols'].keys():
            value_summaries['samplingProtocols'][samproto] += 1
        else:
            value_summaries['samplingProtocols'][samproto] = 1

    # Remove duplicates, make strings for entry into summary table of attributes
    cursor.executescript("""CREATE TABLE record_attributes (step TEXT, field TEXT, vals TEXT);""")
    for x in summary.keys():
        vals = str(list(set(summary[x]))).replace('"', '')
        stmt = """INSERT INTO record_attributes (step, field, vals)
                  VALUES ("request", "{0}", "{1}");""".format(x, vals)
        cursor.execute(stmt)

    # Store the value summary for the selected fields in a table.
    cursor.executescript("""CREATE TABLE post_request_value_counts
                            (attribute TEXT, value TEXT, count INTEGER);""")
    for x in value_summaries.keys():
        attribute = value_summaries[x]
        for y in value_summaries[x].keys():
            z = value_summaries[x][y]
            frog = """INSERT INTO post_request_value_counts (attribute, value, count)
                      VALUES ("{0}", "{1}", "{2}")""".format(x,y,z)
            cursor.execute(frog)
    print("Created summary table of request results: " + str(datetime.now() - requestsummarytime1))
    #                                                                     !!!!!!  SLOW PART HAS ENDED BY HERE

    ##################################################  FILTER MORE
    ###############################################################
    # Pull out relevant attributes from occurrence dictionaries.  Filtering
    # will be performed with info from these keys.
    filtertime1 = datetime.now()
    keykeys = ['basisOfRecord', 'individualCount', 'acceptedTaxonKey',
               'scientificName', 'acceptedScientificName','taxonomicStatus',
               'decimalLongitude', 'decimalLatitude',
               'coordinateUncertaintyInMeters', 'year',
               'month', 'day', 'eventDate', 'issues','geodeticDatum',
               'gbifID', 'type', 'preparations', 'occurrenceStatus',
               'georeferenceProtocol', 'georeferenceVerificationStatus',
               'occurrenceID', 'dataGeneralizations', 'eventRemarks', 'locality',
               'locationRemarks', 'occurrenceRemarks', 'collectionCode',
               'protocol', 'samplingProtocol', 'institutionCode']
    alloccs2 = []
    for x in alloccs:
        alloccs2.append(dict((y,x[y]) for y in x if y in keykeys))

    # Combine remarks FIELDS
    for x in alloccs2:
        remarks = str()
        try:
            put = x['locality']
            remarks = remarks + "; " + put
        except:
            pass

        try:
            these = x['eventRemarks']
            remarks = remarks + "; " + these
        except:
            pass

        try:
            tog = x['locationRemarks']
            remarks = remarks + "; " + tog
        except:
            pass

        try:
            ether = x['occurrenceRemarks']
            remarks = remarks + ether
        except:
            pass

        try:
            x['remarks'] = remarks
        except Exception as e:
            x['remarks'] = ""

    # Identify data generalizations
    for x in alloccs2:
        if 'dataGeneralizations' not in x.keys():
            x['dataGeneralizations'] = ""

    # HAS COORDINATE UNCERTAINTY
    sql_green = """SELECT has_coordinate_uncertainty FROM gbif_filters
                   WHERE filter_id = '{0}';""".format(gbif_filter_id)
    filt_coordUncertainty = cursor2.execute(sql_green).fetchone()[0]

    if filt_coordUncertainty == 1:
        alloccs3 = [x for x in alloccs2 if 'coordinateUncertaintyInMeters'
                    in x.keys()]
    if filt_coordUncertainty == 0:
        alloccs3 = alloccs2
    del alloccs2

    # MAXIMUM COORDINATE UNCERTAINTY
    sql_maxcoord = """SELECT max_coordinate_uncertainty FROM gbif_filters
                   WHERE filter_id = '{0}';""".format(gbif_filter_id)
    filt_maxcoord = cursor2.execute(sql_maxcoord).fetchone()[0]
    alloccs4 = []
    for x in alloccs3:
        if 'coordinateUncertaintyInMeters' not in x.keys():
            alloccs4.append(x)
        elif x['coordinateUncertaintyInMeters'] <= filt_maxcoord:
            alloccs4.append(x)
        else:
            pass
    del alloccs3

    # COLLECTION CODES
    sql_collection = """SELECT collection_codes_omit FROM gbif_filters
                   WHERE filter_id = '{0}';""".format(gbif_filter_id)
    filt_collection = cursor2.execute(sql_collection).fetchone()[0]
    if type(filt_collection) == str:
        filt_collection = list(filt_collection.split(', '))
    else:
        filt_collection = []

    alloccs5 = []
    for x in alloccs4:
        if 'collectionCode' in x.keys():
            if x['collectionCode'] not in list(filt_collection):
                alloccs5.append(x)
            elif 'collectionCode' not in x.keys():
                alloccs5.append(x)
            else:
                pass
        else:
            alloccs5.append(x)
    del alloccs4

    # INSTITUTIONS
    sql_instit = """SELECT institutions_omit FROM gbif_filters
                   WHERE filter_id = '{0}';""".format(gbif_filter_id)
    filt_instit = cursor2.execute(sql_instit).fetchone()[0]
    if type(filt_instit) == str:
        filt_instit = list(filt_instit.split(', '))
    else:
        filt_instit = []

    alloccs6 = []
    for x in alloccs5:
        if 'institutionCode' in x.keys():
            if x['institutionCode'] not in filt_instit:
                alloccs6.append(x)
        else:
            alloccs6.append(x)
    del alloccs5

    # BASES
    sql_bases = """SELECT bases_omit FROM gbif_filters
                   WHERE filter_id = '{0}';""".format(gbif_filter_id)
    filt_bases = cursor2.execute(sql_bases).fetchone()[0]
    if type(filt_bases) == str:
        filt_bases = list(filt_bases.split(', '))
    else:
        filt_bases = []

    alloccs7 = []
    for x in alloccs6:
         if x['basisOfRecord'] not in list(filt_bases):
             alloccs7.append(x)
         elif 'basisOfRecord' not in x.keys():
             alloccs7.append(x)
         else:
             pass
    del alloccs6

    # PROTOCOLS
    sql_protocols = """SELECT protocols_omit FROM gbif_filters
                   WHERE filter_id = '{0}';""".format(gbif_filter_id)
    filt_protocols = cursor2.execute(sql_protocols).fetchone()[0]
    if type(filt_protocols) == str:
        filt_protocols = list(filt_protocols.split(', '))
    else:
        filt_protocols = []

    alloccs8 = []
    for x in alloccs7:
        if x['protocol'] not in filt_protocols:
            alloccs8.append(x)
        elif 'protocol' not in x.keys():
            alloccs8.append(x)
        else:
            pass
    del alloccs7

    # ISSUES
    sql_issues = """SELECT issues_omit FROM gbif_filters
                   WHERE filter_id = '{0}';""".format(gbif_filter_id)
    filt_issues = cursor2.execute(sql_issues).fetchone()[0]

    if type(filt_issues) == str:
        filt_issues = list(filt_issues.split(', '))
    else:
        filt_issues = []

    alloccs9 = []
    for x in alloccs8: # If none of list items are in issues omit list
        if len(set(x['issues']) & set(filt_issues)) == 0:
            alloccs9.append(x)
        elif 'issues' not in x.keys():
            alloccs9.append(x)
        else:
            pass
    del alloccs8

    # SAMPLING PROTOCOL
    sql_sampling = """SELECT sampling_protocols_omit FROM gbif_filters
                   WHERE filter_id = '{0}';""".format(gbif_filter_id)
    filt_sampling = cursor2.execute(sql_sampling).fetchone()[0]
    if type(filt_sampling) == str:
        filt_sampling = list(filt_sampling.split(', '))
    else:
        filt_sampling = []

    alloccsX = []
    for x in alloccs9:
        if 'samplingProtocol' in x.keys():
            if x['samplingProtocol'] not in filt_sampling:
                alloccsX.append(x)
        else:
            alloccsX.append(x)
    del alloccs9
    print("Performed post-request filtering: " + str(datetime.now() - filtertime1))


    ############################# SAVE SUMMARY OF VALUES KEPT (FILTER)
    filtersummarytime1 = datetime.now()
    summary2 = {'datums': ['WGS84'],
               'issues': set([]),
               'bases': [],
               'institutions': [],
               'collections': [],
               'generalizations': set([]),
               'remarks': set([]),
               'establishment': set([]),
               'IDqualifier': set([]),
               'protocols': set([])}

    for occdict in alloccsX:
        # datums
        if occdict['geodeticDatum'] != 'WGS84':
            summary2['datums'] = summary2['datums'] + occdict['geodeticDatum']
        # issues
        summary2['issues'] = summary2['issues'] | set(occdict['issues'])
        # basis of record
        BOR = occdict['basisOfRecord']
        if BOR == "" or BOR == None:
            summary2['bases'] = summary2['bases'] + ["UNKNOWN"]
        else:
            summary2['bases'] = summary2['bases'] + [BOR]
        # institution
        try:
            try:
                who = occdict['institutionID']
                summary2['institutions'] = summary2['institutions'] + [who]
            except:
                who = occdict['institutionCode']
                summary2['institutions'] = summary2['institutions'] + [who]
        except:
            summary2['institutions'] = summary2['institutions'] + ['UNKNOWN']
        # collections
        try:
            co = occdict['collectionCode']
            summary2['collections'] = summary2['collections'] + [co]
        except:
            pass
        # establishment means
        try:
            est = occdict['establishmentMeans']
            summary2['establishment'] = summary2['establishment'] | set([est])
        except:
            pass
        # identification qualifier
        try:
            qual = occdict['identificationQualifier']
            summary2['IDqualifier'] = summary2['IDqualifier'] | set([qual])
        except:
            pass
        # protocols
        try:
            proto = occdict['protocol']
            summary2['protocols'] = summary2['protocols'] | set([proto])
        except:
            pass
        try:
            samproto = occdict['samplingProtocol']
            summary2['protocols'] = summary2['protocols'] | set([samproto])
        except:
            pass

    # Remove duplicates, make strings for entry into table
    for x in summary2.keys():
        stmt = """INSERT INTO record_attributes (step, field, vals)
                  VALUES ("filter", "{0}", "{1}");""".format(x, str(list(set(summary2[x]))).replace('"', ''))
        cursor.execute(stmt)
    print("Summarized results of filtering: " + str(datetime.now() - filtersummarytime1))


    ###############################################  INSERT INTO DB
    ###############################################################
    # Insert the records   !needs to assess if coord uncertainty is present
    # and act accordingly because insert statement depends on if it's present!
    inserttime1 = datetime.now()
    for x in alloccsX:
        try:
            if 'coordinateUncertaintyInMeters' in x.keys() and x['coordinateUncertaintyInMeters'] > 0:
                insert1 = []
                insert1.append((x['gbifID'], species_id, 'gbif',
                                x['decimalLatitude'], x['decimalLongitude'],
                                x['coordinateUncertaintyInMeters'], x['eventDate'],
                                gbif_req_id, gbif_filter_id,
                                x['dataGeneralizations'], x['remarks']))
            else:
                insert1 = []
                insert1.append((x['gbifID'], species_id, 'gbif',
                                x['decimalLatitude'], x['decimalLongitude'],
                                default_coordUncertainty, x['eventDate'],
                                gbif_req_id, gbif_filter_id,
                                x['dataGeneralizations'], x['remarks']))
            insert1 = tuple(insert1)[0]

            sql1 = """INSERT INTO occurrences ('occ_id', 'species_id', 'source',
                                               'latitude', 'longitude',
                                               'coordinateUncertaintyInMeters',
                                               'occurrenceDate', 'request_id',
                                               'filter_id', 'generalizations',
                                               'remarks', 'geom_xy4326')
                        VALUES {0}, GeomFromText('POINT({1} {2})',
                                                    {3}))""".format(str(insert1)[:-1],
                        x['decimalLongitude'], x['decimalLatitude'],
                        SRID_dict[x['geodeticDatum']])
            cursor.executescript(sql1)
        except Exception as e:
            print("\nThere was a problem with the following record:")
            print(e)
            print(x)

    # Update the individual count when it exists
    for e in alloccsX:
        if 'individualCount' in e.keys():
            sql2 = """UPDATE occurrences
                SET individualCount = {0}
                WHERE occ_id = {1};""".format(e['individualCount'], e['gbifID'])
            cursor.execute(sql2)
    conn.commit()
    print("Performed post-request filtering: " + str(datetime.now() - inserttime1))

    ################################################  HANDLE DUPLICATES
    ###################################################################
    OKsql = """SELECT duplicates_OK FROM gbif_filters
                   WHERE filter_id = '{0}';""".format(gbif_filter_id)
    duplicates_OK = cursor2.execute(OKsql).fetchone()[0]
    if duplicates_OK == "False":
        duptime1 = datetime.now()
        conn3 = sqlite3.connect(spdb)
        cursor3 = conn3.cursor()

        # Get a count of duplicates to report
        sql_dupcnt = """SELECT count(occ_id)
                        FROM occurrences
                        WHERE occ_id NOT IN
                            (SELECT occ_id
                             FROM occurrences
                             GROUP BY latitude, longitude, occurrenceDate
                             HAVING max(individualCount));"""
        dupcount = cursor3.execute(sql_dupcnt).fetchone()[0]

        # Delete duplicate records without the highest individualCount among duplicates.
        sql_deldup = """DELETE
                        FROM occurrences
                        WHERE occ_id NOT IN
                            (SELECT occ_id
                             FROM occurrences
                             GROUP BY latitude, longitude, occurrenceDate
                             HAVING max(individualCount));"""
        cursor3.execute(sql_deldup)

        duptime2 = datetime.now()
        print("Removed duplicates: " + str(duptime2 - duptime1))
        print("\t{0} duplicates were deleted".format(dupcount))
    if duplicates_OK == "True":
        print("DUPLICATES ON LATITUDE, LONGITUDE, DATE-TIME INCLUDED")

    ################################################  BUFFER POINTS
    ###############################################################
    # Buffer the x,y locations with the coordinate uncertainty
    # in order to create circles.  Create versions in albers and wgs84.  The
    # wgs84 version will be used in plotting with Basemap.  Buffer radius is
    # the sum of detectiondistance from requests.species_concepts and
    # coordinate uncertainty in meters here.
    buffertime1 = datetime.now()
    requestsDB = inDir + 'requests.sqlite'  #####???????????????????????????????????????
    sql_det = """
            ATTACH DATABASE '{0}' AS requests;

            UPDATE occurrences
            SET detection_distance = {1};

            UPDATE occurrences
            SET radius_meters = detection_distance + coordinateUncertaintyInMeters;

            DETACH DATABASE requests;
    """.format(requestsDB, det_dist)
    cursor.executescript(sql_det)

    sql_buf = """
            /* Transform to albers (102008) and apply buffer */
            ALTER TABLE occurrences ADD COLUMN circle_albers BLOB;

            UPDATE occurrences SET circle_albers = Buffer(Transform(geom_xy4326,
                                                                    102008),
                                                          radius_meters);

            SELECT RecoverGeometryColumn('occurrences', 'circle_albers', 102008,
                                         'POLYGON', 'XY');

            /* Transform back to WGS84 so it can be displayed in iPython */
            ALTER TABLE occurrences ADD COLUMN circle_wgs84 BLOB;

            UPDATE occurrences SET circle_wgs84 = Transform(circle_albers, 4326);

            SELECT RecoverGeometryColumn('occurrences', 'circle_wgs84', 4326,
                                         'POLYGON', 'XY');
    """
    cursor.executescript(sql_buf)
    print("Buffered points: " + str(datetime.now() - buffertime1))


    ##################################################  EXPORT MAPS
    ###############################################################
    exporttime1 = datetime.now()
    # Export occurrence circles as a shapefile (all seasons)
    cursor.execute("""SELECT ExportSHP('occurrences', 'circle_wgs84',
                     '{0}{1}_circles', 'utf-8');""".format(outDir,
                                                           summary_name))

    # Export occurrence 'points' as a shapefile (all seasons)
    cursor.execute("""SELECT ExportSHP('occurrences', 'geom_4326',
                      '{0}{1}_points', 'utf-8');""".format(outDir,
                                                           summary_name))
    conn.commit()
    #conn.close()
    conn2.commit()
    conn2.close()
    print("Exported maps: " + str(datetime.now() - exporttime1))
    print("\nRecords saved in {0}".format(spdb))

def ccw_wkt_from_shp(shapefile, out_txt):
    """
    Creates wkt with coordinates oriented counter clockwise for a given shapefile.
    Shapefiles are oriented clockwise, which is incompatible with spatial queries
    in many database management systems.  Use this to generate wkt that you can
    copy and paste into queries.

    (str, str) = text written to shpefile

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

    if c.crs['init'] == 'epsg:4326':
        # Next make it a shapely polygon object
        poly = shape(c[0]['geometry'])

        # Use LinearRing to determine if coordinates are listed clockwise
        coords = c[0]["geometry"]["coordinates"][0]
        lr = LinearRing(coords)
        if lr.is_ccw == False:
            # Reverse coordinates to make them counter clockwise
            print("Points were clockwise......reversing")
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
        with open(out_txt, 'w+') as file:
            file.write(wkt)
            print("WKT written to {0}".format(out_txt))

        # close the collections
        c.close()
    else:
        print("You need to reproject the shapefile to EPSG:4326")
    return
