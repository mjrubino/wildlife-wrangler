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


##################################################  POST REQUEST FILTERING
##########################################################################
# New colummn for compiled remarks and notes
df0["remarks"] = df0['locality'] + ";" + df0['eventRemarks'] + ";" + df0['locationRemarks'] + ";" + df0['occurrenceRemarks']

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
df6 = df5[df5['protocol'].isin(filt_bases) == False]
del df5
df7 = df6[df6['samplingProtocol'].isin(filt_sampling) == False]]
del df6
# ISSUES - this one is more complex because multiple issues can be listed per record
# Method used is complex, but hopefully faster than simple iteration over all records
unique_issue = df7['issue'].unique() # List of unique issue entries
violations = [x for x in unique_issue if len(set(x.split(";")) & set(filt_issues)) == 0] # entries that contain violations
df8 = df7[df7['issue'].isin(violations) == False] # Records without entries that are violations.

print("Performed post-request filtering: " + str(datetime.now() - filtertime1))
