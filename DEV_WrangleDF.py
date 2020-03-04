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

########################################################################
########################################################################
########################################################################





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
df7 = df6[df6['samplingProtocol'].isin(filt_sampling) == False]
del df6
# ISSUES - this one is more complex because multiple issues can be listed per record
# Method used is complex, but hopefully faster than simple iteration over all records
unique_issue = df7['issue'].unique() # List of unique issue entries
violations = [x for x in unique_issue if len(set(x.split(";")) & set(filt_issues)) == 0] # entries that contain violations
df8 = df7[df7['issue'].isin(violations) == False] # Records without entries that are violations.

print("Performed post-request filtering: " + str(datetime.now() - filtertime1))

################################################ SUMMARIZE VALUES PRESENT IN TABLE
################################################################################
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

# issues
summary['issues'] = get_vals(df8, 'issue')

# basis or record
summary['bases'] = get_vals(df8, 'basisOfRecord')

# institution
summary['institutions'] = get_vals(df8, 'institutionID') | get_vals(df8, 'institutionCode')

# collections
summary['collections'] = get_vals(df8, 'collectionCode')

# establishment means
summary['establishment'] = get_vals(df8, 'establishmentMeans')

# identification qualifier
summary['IDqualifier'] = get_vals(df8, 'identificationQualifier')

# protocols
summary['protocols'] = get_vals(df8, 'protocol') | get_vals(df8, 'samplingProtocol')

print(summary)

# Remove duplicates, make strings for entry into summary table of attributes
for x in summary.keys():
    vals = str(list(set(summary[x]))).replace('"', '')
    stmt = """INSERT INTO record_attributes (step, field, vals)
              VALUES ("filter", "{0}", "{1}");""".format(x, vals)
    cursor.execute(stmt)
