import pandas as pd
file = "T:/temp/dfOcc.csv"
keykeys = ['basisOfRecord', 'individualCount', 'scientificName',
           'decimalLongitude', 'decimalLatitude', 'coordinateUncertaintyInMeters',
           'eventDate', 'issue', 'id', 'occurrenceID', 'dataGeneralizations',
           'eventRemarks', 'locality', 'locationRemarks', 'collectionCode',
           'protocol', 'samplingProtocol', 'institutionCode',
           'establishmentMeans', 'institutionID', 'identificationQualifier']


df0 = pd.read_csv(file, nrows=1000, usecols=keykeys)
print(df0.head())

# Count entries per atrribute(column), reformat as new df with appropriate
# columns.  Finally, insert into db.
df_populated1 = pd.DataFrame(df0.count(axis=0).T.iloc[1:])
df_populated1['included(n)'] = len(df0)
df_populated1['populated(n)'] = df_populated1[0]
df_populated2 = df_populated1.filter(items=['included(n)', 'populated(n)'], axis='columns')
df_populated2.index.name = 'attribute'
print(df_populated2.head(25))
#df_populated2.to_sql(name='gbif_fields_returned', con=conn, if_exists='replace')

# Create a table for storing unique attribute values that came back.
# For each column of interest (listed in summary keys) in the dataframe,
# get a list of unique values and insert table.
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

def get_vals(column_name):
    '''
    Return a set of unique values from a column
    '''
    stoat = df0[column_name].unique()
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
summary['issues'] = get_vals('issue')

# basis or record
summary['bases'] = get_vals('basisOfRecord')

# institution
summary['institutions'] = get_vals('institutionID') | get_vals('institutionCode')

# collections
summary['collections'] = get_vals('collectionCode')

# establishment means
summary['establishment'] = get_vals('establishmentMeans')

# identification qualifier
summary['IDqualifier'] = get_vals('identificationQualifier')

# protocols
summary['protocols'] = get_vals('protocol') | get_vals('samplingProtocol')

print(summary)
#print(get_vals('protocol'))

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
for x in value_counts.keys():
    attribute = value_counts[x]
    for y in value_counts[x].keys():
        z = value_counts[x][y]
        frog = """INSERT INTO post_request_value_counts (attribute, value, count)
                  VALUES ("{0}", "{1}", "{2}")""".format(x,y,z)
        cursor.execute(frog)
#print("\t Slow part 2 : " + str(datetime.now() - slow1))
print("Created summary table of request results: " + str(datetime.now() - requestsummarytime1))

###########
for x in df0.columns:
    if 'abstract' in x:
        print(x)
###########