# Define a function for displaying the maps that will be created.
def MapShapefilePolygons(map_these, title):
    """
    Displays shapefiles on a simple CONUS basemap.  Maps are plotted in the order
    provided so put the top map last in the list.  You can specify a column
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
                              gbif_filter_id, default_coordUncertainty,
                              outDir, summary_name, username, password, email):
    """
    Retrieves GAP range from ScienceBase and occurrence records from APIs. Filters
    occurrence records, stores them in a database, buffers the xy points,
    and filtering occurrence records, saving them in a database.  Finally, exports
    some maps.

    Arguments:
    codeDir -- directory of this code repo.
    species_id -- project id for the species concept.
    inDir -- directory containing key inputs such as downloaded gap ranges.
    spdb -- occurrence record database to be created by this function.
    gbif_req_id -- GBIF request ID for the process.
    gbif_filter_id -- GBIF filter ID for the process.
    default_coordUncertainty -- distance in meters to use if no coordinate
        Uncertainty is specified for a record.
    outDir -- where to save maps that are exported by this process.
    summary_name -- a short name for some file names.
    sp_geometry -- True or False to use geometry saved with species concept when
        filtering records.  Request geometry is always used if provided.
    """
    sp_geometry = True #  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
    import sys
    import shutil
    from dwca.read import DwCAReader


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
    conn2 = sqlite3.connect(codeDir + 'parameters.sqlite', isolation_level='DEFERRED')
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


    ############################################################################
    #####################    Create Occurrence Database    #####################
    ############################################################################
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
    conn = sqlite3.connect(spdb, isolation_level='DEFERRED')
    conn.enable_load_extension(True)
    conn.execute('SELECT load_extension("mod_spatialite")')
    cursor = conn.cursor()

    # Make database spatial
    conn.executescript('''SELECT InitSpatialMetaData(1);''')
    conn.commit()

    ################################################# Create tables
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
                    footprintWKT TEXT,
                    weight INTEGER DEFAULT 10,
                    weight_notes TEXT,
                    doi_search TEXT,
                        FOREIGN KEY (species_id) REFERENCES taxa(species_id)
                        ON UPDATE RESTRICT
                        ON DELETE NO ACTION);

            SELECT AddGeometryColumn('occurrences', 'geom_xy4326', 4326, 'POINT',
                                     'XY');
    """
    cursor.executescript(sql_cdb)
    makedb2 = datetime.now()
    print("Created occurrence db: " + str(makedb2 - makedb1))

    requesttime1 = datetime.now()

    ############################################################################
    ###########################   Get Filters   ################################
    ############################################################################
    """
    Retrieve filter parameters from the parameters database.
    """
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
    # It could also be that user opted not to use species geometry.
    if sp_geometry == False:
        sp_geom = None

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

    ###################  RETRIEVE POST-REQUEST FILTER PARAMTERS
    sql_green = """SELECT has_coordinate_uncertainty FROM gbif_filters
                   WHERE filter_id = '{0}';""".format(gbif_filter_id)
    filt_coordUncertainty = cursor2.execute(sql_green).fetchone()[0]

    sql_maxcoord = """SELECT max_coordinate_uncertainty FROM gbif_filters
                   WHERE filter_id = '{0}';""".format(gbif_filter_id)
    filt_maxcoord = cursor2.execute(sql_maxcoord).fetchone()[0]

    sql_collection = """SELECT collection_codes_omit FROM gbif_filters
                   WHERE filter_id = '{0}';""".format(gbif_filter_id)
    filt_collection = cursor2.execute(sql_collection).fetchone()[0]
    if type(filt_collection) == str:
        filt_collection = list(filt_collection.split(', '))
    else:
        filt_collection = []

    sql_instit = """SELECT institutions_omit FROM gbif_filters
                   WHERE filter_id = '{0}';""".format(gbif_filter_id)
    filt_instit = cursor2.execute(sql_instit).fetchone()[0]
    if type(filt_instit) == str:
        filt_instit = list(filt_instit.split(', '))
    else:
        filt_instit = []

    sql_bases = """SELECT bases_omit FROM gbif_filters
                   WHERE filter_id = '{0}';""".format(gbif_filter_id)
    filt_bases = cursor2.execute(sql_bases).fetchone()[0]
    if type(filt_bases) == str:
        filt_bases = list(filt_bases.split(', '))
    else:
        filt_bases = []

    sql_issues = """SELECT issues_omit FROM gbif_filters
                   WHERE filter_id = '{0}';""".format(gbif_filter_id)
    filt_issues = cursor2.execute(sql_issues).fetchone()[0]
    if type(filt_issues) == str:
        filt_issues = list(filt_issues.split(', '))
    else:
        filt_issues = []

    sql_sampling = """SELECT sampling_protocols_omit FROM gbif_filters
                   WHERE filter_id = '{0}';""".format(gbif_filter_id)
    filt_sampling = cursor2.execute(sql_sampling).fetchone()[0]
    if type(filt_sampling) == str:
        filt_sampling = list(filt_sampling.split(', '))
    else:
        filt_sampling = []

    print("Got request params and sorted out geometry constraints: " + str(datetime.now() - requesttime1))
    requesttime2 = datetime.now()

    # List of informative df columns/dictionary keys to keep (used later)
    keeper_keys = ['basisOfRecord', 'individualCount', 'scientificName',
                   'decimalLongitude', 'decimalLatitude',
                   'coordinateUncertaintyInMeters',
                   'eventDate', 'issue', 'issues', 'gbifID', 'id',
                   'dataGeneralizations', 'eventRemarks', 'locality',
                   'locationRemarks', 'collectionCode',
                   'samplingProtocol', 'institutionCode', 'institutionID'
                   'establishmentMeans', 'institutionID', 'footprintWKT',
                   'identificationQualifier', 'occurrenceRemarks', 'datasetName']

    ############################################################################
    #######################      Get GBIF Records      #########################
    ############################################################################
    """
    Retrieve GBIF records for a species and save appropriate
    attributes in the occurrence db.
    """
    ####################################### HOW MANY RECORDS EXIST TO PULL FROM?
    ############################################################################
    # First, find out how many records there are that meet criteria
    occ_search = occurrences.search(gbif_id,
                                    year=years,
                                    month=months,
                                    decimalLatitude=latRange,
                                    decimalLongitude=lonRange,
                                    hasGeospatialIssue=geoIssue,
                                    hasCoordinate=coordinate,
                                    country=country,
                                    geometry=poly)
    occ_count=occ_search['count']
    print(str(occ_count) + " records available")


    ############################################################################
    #                         < 100,000 RECORDS (JSON)
    ############################################################################
    if occ_count <100000:
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
                                          country=country,
                                          geometry=poly)
            occs = occ_json['results']
            alloccs = alloccs + occs

        print("Downloaded records: " + str(datetime.now() - requesttime2))
        print('\t{0} records exist with the request parameters'.format(occ_count))

        ##########################  SUMMARY TABLE OF KEYS/FIELDS RETURNED (JSON)
        ########################################################################
        """                                                                         ### TOO SLOW
        keys = [list(x.keys()) for x in alloccs]
        keys2 = set([])
        for x in keys:
            keys2 = keys2 | set(x)
        dfK = pd.DataFrame(index=keys2, columns=['included(n)', 'populated(n)'])
        dfK['included(n)'] = 0
        dfK['populated(n)'] = 0
        """
        requestsummarytime1 = datetime.now()
        """                                    #####  START SLOW
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
        """
        slow1 = datetime.now()
        """                                                       #######
        print("\t Slow part 1 : " + str(slow1 - requestsummarytime1))                 # Timer
        #                                                                     !!!!!!  SLOW PART HAS ENDED BY HERE

        dfK.sort_index(inplace=True)
        dfK.to_sql(name='gbif_fields_returned', con=conn, if_exists='replace')
        """

        ############################# SUMMARY OF VALUES RETURNED (JSON; REQUEST)
        ########################################################################
        summary = {'datums': ['WGS84'],
                   'issues': set([]),
                   'bases': [],
                   'institutions': [],
                   'collections': [],
                   'generalizations': set([]),
                   'remarks': set([]),
                   'establishment': set([]),
                   'IDqualifier': set([]),
                   'samplingProtocol': set([])}

        value_counts = {'bases': {},
                          'datums': {'WGS84': 0},
                          'issues': {},
                          'institutions': {},
                          'collections': {},
                          'samplingProtocols': {}}

        for occdict in alloccs:
            # datums
            if occdict['geodeticDatum'] != 'WGS84':
                summary['datums'] = summary['datums'] + occdict['geodeticDatum']
                if occdict['geodeticDatum'] not in value_counts['datums'].keys():
                    value_counts['datums'][occdict['geodeticDatum']] = 0
                else:
                    value_counts['datums'][occdict['geodeticDatum']] += 1
            if occdict['geodeticDatum'] == 'WGS84':
                value_counts['datums']['WGS84'] += 1

            # issues
            summary['issues'] = summary['issues'] | set(occdict['issues'])
            for issue in occdict['issues']:
                if issue not in value_counts['issues'].keys():
                    value_counts['issues'][issue] = 1
                if issue in value_counts['issues'].keys():
                    value_counts['issues'][issue] += 1

            # basis or record
            BOR = occdict['basisOfRecord']
            if BOR == "" or BOR == None:
                BOR = 'UNKNOWN'
            summary['bases'] = summary['bases'] + [BOR]

            if BOR in value_counts['bases'].keys():
                value_counts['bases'][BOR] += 1
            else:
                value_counts['bases'][BOR] = 1

            # institution
            try:
                who = occdict['institutionID']
            except:
                try:
                    who = occdict['institutionCode']
                except:
                    who = 'unknown'

            summary['institutions'] = summary['institutions'] + [who]

            if who in value_counts['institutions'].keys():
                value_counts['institutions'][who] += 1
            else:
                value_counts['institutions'][who] = 1

            # collections
            try:
                co = occdict['collectionCode']
            except:
                co = 'UNKNOWN'

            summary['collections'] = summary['collections'] + [co]

            if co in value_counts['collections'].keys():
                value_counts['collections'][co] += 1
            else:
                value_counts['collections'][co] = 1

            # dataset
            try:
                dn = occdict['datasetName']
            except:
                dn = 'UNKNOWN'

            summary['datasetName'] = summary['datasetName'] + [dn]

            if dn in value_counts['datasetName'].keys():
                value_counts['datasetName'][dn] += 1
            else:
                value_counts['datasetName'][dn] = 1

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

            # sampling protocol
            try:
                samproto = occdict['samplingProtocol']
            except:
                samproto = 'UKNOWN'
            summary['samplingProtocols'] = set([samproto])

            if samproto in value_counts['samplingProtocols'].keys():
                value_counts['samplingProtocols'][samproto] += 1
            else:
                value_counts['samplingProtocols'][samproto] = 1

        # Remove duplicates, make strings for entry into summary table of attributes
        cursor.executescript("""CREATE TABLE unique_values (step TEXT, field TEXT, vals TEXT);""")
        for x in summary.keys():
            vals = str(list(set(summary[x]))).replace('"', '')
            stmt = """INSERT INTO unique_values (step, field, vals)
                      VALUES ("request", "{0}", "{1}");""".format(x, vals)
            cursor.execute(stmt)

        # Store the value summary for the selected fields in a table.
        cursor.executescript("""CREATE TABLE pre_filter_value_counts
                                (attribute TEXT, value TEXT, count INTEGER);""")
        for x in value_counts.keys():
            attribute = value_counts[x]
            for y in value_counts[x].keys():
                z = value_counts[x][y]
                frog = """INSERT INTO pre_filter_value_counts (attribute, value, count)
                          VALUES ("{0}", "{1}", "{2}")""".format(x,y,z)
                cursor.execute(frog)
        #print("\t Slow part 2 : " + str(datetime.now() - slow1))
        print("Created summary table of request results: " + str(datetime.now() - requestsummarytime1))

        #########################################################  FILTER (JSON)
        ########################################################################
        # Pull out relevant attributes from occurrence dictionaries.  Filtering
        # will be performed with info from these keys.
        filtertime1 = datetime.now()

        alloccs2 = []
        for x in alloccs:
            alloccs2.append(dict((y,x[y]) for y in x if y in keeper_keys))

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
        if filt_coordUncertainty == 1:
            alloccs3 = [x for x in alloccs2 if 'coordinateUncertaintyInMeters'
                        in x.keys()]
        if filt_coordUncertainty == 0:
            alloccs3 = alloccs2
        del alloccs2

        # MAXIMUM COORDINATE UNCERTAINTY
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
        alloccs6 = []
        for x in alloccs5:
            if 'institutionCode' in x.keys():
                if x['institutionCode'] not in filt_instit:
                    alloccs6.append(x)
            else:
                alloccs6.append(x)
        del alloccs5

        # BASES
        alloccs7 = []
        for x in alloccs6:
             if x['basisOfRecord'] not in list(filt_bases):
                 alloccs7.append(x)
             elif 'basisOfRecord' not in x.keys():
                 alloccs7.append(x)
             else:
                 pass
        del alloccs6

        # ISSUES
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
        alloccsX = []
        for x in alloccs9:
            if 'samplingProtocol' in x.keys():
                if x['samplingProtocol'] not in filt_sampling:
                    alloccsX.append(x)
            else:
                alloccsX.append(x)
        del alloccs9
        print("Performed post-request filtering: " + str(datetime.now() - filtertime1))


        ################################## SUMMARY OF VALUES KEPT (FILTER; JSON)
        ########################################################################
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
                   'samplingProtocols': set([])}

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
                samproto = occdict['samplingProtocol']
                summary2['samplingProtocols'] = summary2['samplingProtocols'] | set([samproto])
            except:
                pass

        # Remove duplicates, make strings for entry into table
        for x in summary2.keys():
            stmt = """INSERT INTO unique_values (step, field, vals)
                      VALUES ("filter", "{0}", "{1}");""".format(x, str(list(set(summary2[x]))).replace('"', ''))
            cursor.execute(stmt)
        print("Summarized results of filtering: " + str(datetime.now() - filtersummarytime1))


        #################################################  INSERT INTO DB (JSON)
        ########################################################################
        # Insert the records   !needs to assess if coord uncertainty is present
        # and act accordingly because insert statement depends on if it's present!
        inserttime1 = datetime.now()
        if default_coordUncertainty != False:
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
                                                       'remarks')
                              VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)"""
                    cursor.executemany(sql1, [(insert1)])

                except Exception as e:
                    print("\nThere was a problem with the following record:")
                    print(e)
                    print(x)
            print("Inserted records: " + str(datetime.now() - inserttime1))

        # Update the individual count when it exists
        inserttime3 = datetime.now()
        for e in alloccsX:
            if 'individualCount' in e.keys():
                sql2 = """UPDATE occurrences
                    SET individualCount = {0}
                    WHERE occ_id = {1};""".format(e['individualCount'], e['gbifID'])
                cursor.execute(sql2)
        conn.commit()
        print("Updated individuaCount column: " + str(datetime.now() - inserttime3))

    ############################################################################
    #                         > 100,000 RECORDS (DF)
    ############################################################################
    else:
        ########################################################## DOWNLOAD (DF)
        ########################################################################
        # Make the data request using the download function.  Results are
        # emailed.
        # First, build a query list.  NoneType values cause problems, so only
        # add arguments if their value isn't NoneType.
        download_filters = ['taxonKey = {0}'.format(gbif_id)]

        if coordinate != None:
            download_filters.append('hasCoordinate = {0}'.format(coordinate))

        if country != None:
            download_filters.append('country = {0}'.format(country))

        if years != None:
            download_filters.append('year >= {0}'.format(years.split(",")[0]))
            download_filters.append('year <= {0}'.format(years.split(",")[1]))

        if months != None:
            download_filters.append('month >= {0}'.format(months.split(",")[0]))
            download_filters.append('month <= {0}'.format(months.split(",")[1]))

        if poly != None:
            download_filters.append('geometry within {0}'.format(poly))

        if geoIssue != None:
            download_filters.append('hasGeospatialIssue = {0}'.format(geoIssue))

        if latRange != None:
            download_filters.append('decimalLatitude >= {0}'.format(latRange.split(",")[0]))
            download_filters.append('decimalLatitude <= {0}'.format(latRange.split(",")[1]))
        if lonRange !=None:
            download_filters.append('decimalLongitude >= {0}'.format(lonRange.split(",")[0]))
            download_filters.append('decimalLongitude <= {0}'.format(lonRange.split(",")[1]))

        bigdown1 = datetime.now()
        d = occurrences.download(download_filters,
                                 pred_type='and',
                                 user = username,
                                 pwd = password,
                                 email = email)
        # Get the value of the download key
        dkey = d[0]

        # Now download the actual zip file containing the Darwin Core files
        # NOTE: The download can take a while to generate and is not immediately
        # available once the download_get command has been issued. Use a
        # while and try loop to make sure the download has succeeded.
        # The zipdownload variable will be a dictionary of the path,
        # the file size, and the download key unique code. It can be used
        # to change the file name, unzip the file, etc.
        print("Downloading Darwin Core Archive zip file for this species .....")
        gotit = None
        while gotit is None:
            try:
                zipdownload = occurrences.download_get(key=dkey, path=inDir)
                gotit = 1
            except:
                pass
        print("Download complete: " + str(datetime.now() - bigdown1))

        # Read the "occurrence.txt" file into a Pandas dataframe
        read1 = datetime.now()
        with DwCAReader(inDir + dkey + '.zip') as dwca:
            dfRaw = dwca.pd_read('occurrence.txt', low_memory=False)#, usecols=keeper_keys)

        df0 = dfRaw.filter(items=keeper_keys, axis=1)


        ####################################################  RENAME FIELDS (DF)
        ########################################################################
        df0.rename(mapper={"id": "occ_id",
                           "decimalLatitude": "latitude",
                           "decimalLongitude": "longitude",
                           "eventDate": "occurrenceDate"}, inplace=True, axis='columns')

        df0.to_csv("T:/temp/dfOcc.csv")
        print("Reading and saving downloaded records: " + str(datetime.now() - read1))

        ############################  SUMMARY TABLE OF KEYS/FIELDS RETURNED (DF)
        ########################################################################
        # Count entries per atrribute(column), reformat as new df with appropriate
        # columns.  Finally, insert into db.
        cheese1 = datetime.now()
        df_populated1 = pd.DataFrame(dfRaw.count(axis=0).T.iloc[1:])
        df_populated1['included(n)'] = len(dfRaw)
        df_populated1['populated(n)'] = df_populated1[0]
        df_populated2 = df_populated1.filter(items=['included(n)', 'populated(n)'], axis='columns')
        df_populated2.index.name = 'attribute'
        df_populated2.to_sql(name='gbif_fields_returned', con=conn, if_exists='replace')
        print("Summarized fields returned: " + str(datetime.now() - cheese1))

        ############################### SUMMARY OF VALUES RETURNED (DF; REQUEST)
        ########################################################################
        # Create a table for storing unique attribute values that came back.
        breadtime = datetime.now()
        summary = {'datums': ['WGS84'],
                   'issues': set([]),
                   'bases': [],
                   'institutions': [],
                   'collections': [],
                   'datasets':[],
                   'generalizations': set([]),
                   'remarks': set([]),
                   'establishment': set([]),
                   'IDqualifier': set([]),
                   'samplingProtocols': set([])}

        value_counts = {'bases': {},
                          'datums': {'WGS84': 0},
                          'issues': {},
                          'institutions': {},
                          'collections': {},
                          'datasets': {},
                          'samplingProtocols': {}}

        def get_vals(df, column_name):
            '''
            Return a set of unique values from a column
            '''
            stoat = df[column_name].unique()
            stoat = [str(x).split(";") for x in stoat]
            stoat1 = []
            for x in stoat:
                for y in x:
                    if y == "" or y == None:
                        stoat1.append('UNKNOWN') # ? Keep?
                    else:
                        stoat1.append(y)
            return set(stoat1)

        # datums - ? - couldn't find this info in the table

        # issues
        summary['issues'] = get_vals(df0, 'issue')

        group = df0['occ_id'].groupby(df0['issue'])
        gemstone = group.count()
        for x in gemstone.index:
            value_counts['issues'][x] = gemstone[x]

        # basis or record
        summary['bases'] = get_vals(df0, 'basisOfRecord')

        group = df0['occ_id'].groupby(df0['basisOfRecord'])
        gemstone = group.count()
        for x in gemstone.index:
            value_counts['bases'][x] = gemstone[x]

        # institution
        summary['institutions'] = get_vals(df0, 'institutionCode')

        group = df0['occ_id'].groupby(df0['institutionCode'])
        gemstone = group.count()
        for x in gemstone.index:
            value_counts['institutions'][x] = gemstone[x]

        # collections
        summary['collections'] = get_vals(df0, 'collectionCode')

        group = df0['occ_id'].groupby(df0['collectionCode'])
        gemstone = group.count()
        for x in gemstone.index:
            value_counts['collections'][x] = gemstone[x]

        # datasets
        summary['datasets'] = get_vals(df0, 'datasetName')

        group = df0['occ_id'].groupby(df0['datasetName'])
        gemstone = group.count()
        for x in gemstone.index:
            value_counts['datasets'][x] = gemstone[x]

        # establishment means
        try:
            summary['establishment'] = get_vals(df0, 'establishmentMeans')
        except:
            summary['establishment'] = ""

        # identification qualifier
        summary['IDqualifier'] = get_vals(df0, 'identificationQualifier')

        # protocols
        summary['samplingProtocols'] = get_vals(df0, 'samplingProtocol')

        group = df0['occ_id'].groupby(df0['samplingProtocol'])
        gemstone = group.count()
        for x in gemstone.index:
            value_counts['samplingProtocols'][x] = gemstone[x]


        # Remove duplicates, make strings for entry into summary table of attributes
        cursor.executescript("""CREATE TABLE unique_values (step TEXT, field TEXT, vals TEXT);""")
        for x in summary.keys():
            vals = str(list(set(summary[x]))).replace('"', '')
            stmt = """INSERT INTO unique_values (step, field, vals)
                      VALUES ("request", "{0}", "{1}");""".format(x, vals)
            cursor.execute(stmt)

        # Store the value summary for the selected fields in a table.
        cursor.executescript("""CREATE TABLE pre_filter_value_counts
                                (attribute TEXT, value TEXT, count INTEGER);""")
        for x in value_counts.keys():
            attribute = value_counts[x]
            for y in value_counts[x].keys():
                z = value_counts[x][y]
                frog = """INSERT INTO pre_filter_value_counts (attribute, value, count)
                          VALUES ("{0}", "{1}", "{2}")""".format(x,y,z)
                cursor.execute(frog)
        print("Created summary table of request results: " + str(datetime.now() - breadtime))


        ######################################  SUMMARIZE SOURCES PRE FILTER (DF)
        ########################################################################
        #
        moss = df0.groupby(['institutionCode', 'collectionCode', 'datasetName'])[['occ_id']].size()
        moss.to_sql(name='pre_filter_source_counts', con = conn, if_exists='replace')


        ##########################################  ADD SOME DEFAULT VALUES (DF)
        ########################################################################
        if default_coordUncertainty != False:
            df0.fillna(value={'coordinateUncertaintyInMeters': default_coordUncertainty,
                              'individualCount': int(1)},
                       inplace=True)


        ###########################################################  FILTER (DF)
        ########################################################################
        fiddlertime = datetime.now()
        # HAS COORDINATE UNCERTAINTY
        if filt_coordUncertainty == 1:
            df1 = df0[pd.isnull(df0['coordinateUncertaintyInMeters']) == False]
        if filt_coordUncertainty == 0:
            df1 = df0
        # OTHER FILTERS
        df2 = df1[df1['coordinateUncertaintyInMeters'] <= filt_maxcoord]
        del df1
        df3 = df2[df2['collectionCode'].isin(filt_collection) == False]
        del df2
        df4 = df3[df3['institutionCode'].isin(filt_instit) == False]
        del df3
        df5 = df4[df4['basisOfRecord'].isin(filt_bases) == False]
        del df4
        df7 = df5[df5['samplingProtocol'].isin(filt_sampling) == False]
        del df5
        # ISSUES - this one is more complex because multiple issues can be listed per record
        # Method used is complex, but hopefully faster than simple iteration over all records
        df7.fillna(value={'issue': ""}, inplace=True)
        unique_issue = list(df7['issue'].unique())
        violations = [x for x in unique_issue if len(set(str(x).split(";")) & set(filt_issues)) != 0] # entries that contain violations

        df8 = df7[df7['issue'].isin(violations) == False] # Records without entries that are violations.
        del df7
        print("Performed post-request filtering: " + str(datetime.now() - fiddlertime))

        newstime = datetime.now()
        # Create any new columns needed
        df8["remarks"] = df8['locality'] + ";" + df8['eventRemarks'] + ";" + df8['locationRemarks'] + ";" + df8['occurrenceRemarks']
        df8["species_id"] = species_id
        df8["request_id"] = gbif_req_id
        df8["filter_id"] = gbif_filter_id
        df8["retrievalDate"] = datetime.now()
        df8["detection_distance"] = det_dist
        df8["radius_meters"] = df8["detection_distance"] + df8["coordinateUncertaintyInMeters"]
        df8["source"] = "gbif"
        df8["doi_search"] = ""
        df8["weight"] = 10
        df8["weight_notes"] = ""
        print("Calculated new columns: " + str(datetime.now() - newstime))

        ###################################################  INSERT INTO DB (DF)
        ########################################################################
        biggin = datetime.now()
        '''
        sql1 = """INSERT INTO occurrences ('occ_id', 'species_id', 'source',
                                           'latitude', 'longitude',
                                           'coordinateUncertaintyInMeters',
                                           'occurrenceDate', 'request_id',
                                           'filter_id', 'generalizations',
                                           'remarks')
                  VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?);"""
        for x in df8.index:
            insert2 = [df8.loc[x,"id"], species_id, df8.loc[x,"source"],
                       df8.loc[x,"decimalLatitude"], df8.loc[x,"decimalLongitude"],
                       df8.loc[x,"coordinateUncertaintyInMeters"],
                       df8.loc[x,"eventDate"], request_id, filter_id,
                       df8.loc[x,"dataGeneralizations"], df8.loc[x,"remarks"]]
            cursor.execute(sql1, [(insert2)])
        conn.commit()
        '''
        df8.to_sql(name='occurrences', con = conn, if_exists='replace',
                   chunksize=1000)
        sql_toad = '''SELECT AddGeometryColumn('occurrences', 'geom_xy4326', 4326,
                                        'POINT', 'XY');'''
        cursor.execute(sql_toad)

        print("Inserted records into table: " + str(datetime.now() - biggin))

        ################################## SUMMARY OF VALUES KEPT (FILTER; JSON)
        ########################################################################
        kepttime = datetime.now()
        summary = {'datums': ['WGS84'],
                   'issues': set([]),
                   'bases': [],
                   'institutions': [],
                   'collections': [],
                   'generalizations': set([]),
                   'remarks': set([]),
                   'establishment': set([]),
                   'IDqualifier': set([])}

        # issues
        summary['issues'] = get_vals(df8, 'issue')

        # basis or record
        summary['bases'] = get_vals(df8, 'basisOfRecord')

        # institution
        summary['institutions'] = get_vals(df8, 'institutionID') | get_vals(df8, 'institutionCode')

        # collections
        summary['collections'] = get_vals(df8, 'collectionCode')

        # establishment means
        try:
            summary['establishment'] = get_vals(df8, 'establishmentMeans')
        except:
            summary['establishment'] = ""

        # identification qualifier
        summary['IDqualifier'] = get_vals(df8, 'identificationQualifier')

        # protocols
        summary['samplingProtocols'] = get_vals(df8, 'samplingProtocol')

        # Remove duplicates, make strings for entry into summary table of attributes
        for x in summary.keys():
            vals = str(list(set(summary[x]))).replace('"', '')
            stmt = """INSERT INTO unique_values (step, field, vals)
                      VALUES ("filter", "{0}", "{1}");""".format(x, vals)
            cursor.execute(stmt)
        print("Summarized unique values retained: " + str(datetime.now() - kepttime))

    ################################################  MAKE POINT GEOMETRY COLUMN
    ############################################################################
    inserttime2 = datetime.now()
    try:
        sql2 = """UPDATE occurrences
                  SET geom_xy4326 = GeomFromText('POINT('||"longitude"||' '||"latitude"||')', 4326);"""
        cursor.execute(sql2)
    except Exception as e:
        print(e)

    ###### EVENTUALLY ADD CODE TO OVERIDE POLYGON GEOMETRY WITH FOOTPRINT

    print("Updated occurrences table geometry column: " + str(datetime.now() - inserttime2))

    #########################################################  HANDLE DUPLICATES
    ############################################################################
    OKsql = """SELECT duplicates_OK FROM gbif_filters
                   WHERE filter_id = '{0}';""".format(gbif_filter_id)
    duplicates_OK = cursor2.execute(OKsql).fetchone()[0]

    conn2.commit()
    conn2.close()
    del cursor2

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

        conn3.commit()
        conn3.close()
        del cursor3

    if duplicates_OK == "True":
        print("DUPLICATES ON LATITUDE, LONGITUDE, DATE-TIME INCLUDED")


    #############################################################  BUFFER POINTS
    ############################################################################
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
            /* Transform to albers (5070) and apply buffer */
            ALTER TABLE occurrences ADD COLUMN polygon_5070 BLOB;

            UPDATE occurrences SET polygon_5070 = Buffer(Transform(geom_xy4326,
                                                                    5070),
                                                          radius_meters);

            SELECT RecoverGeometryColumn('occurrences', 'polygon_5070', 5070,
                                         'POLYGON', 'XY');

            /* Transform back to WGS84 so it can be displayed in iPython */
            ALTER TABLE occurrences ADD COLUMN polygon_4326 BLOB;

            UPDATE occurrences SET polygon_4326 = Transform(polygon_5070, 4326);

            SELECT RecoverGeometryColumn('occurrences', 'polygon_4326', 4326,
                                         'POLYGON', 'XY');
    """
    cursor.executescript(sql_buf)
    print("Buffered points: " + str(datetime.now() - buffertime1))


    ###############################################################  EXPORT MAPS
    ############################################################################
    exporttime1 = datetime.now()
    # Export occurrence circles as a shapefile (all seasons)
    cursor.execute("""SELECT ExportSHP('occurrences', 'polygon_4326',
                     '{0}{1}_polygons', 'utf-8');""".format(outDir,
                                                           summary_name))

    # Export occurrence 'points' as a shapefile (all seasons)
    cursor.execute("""SELECT ExportSHP('occurrences', 'geom_4326',
                      '{0}{1}_points', 'utf-8');""".format(outDir,
                                                           summary_name))
    conn.commit()
    #conn.close()

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
