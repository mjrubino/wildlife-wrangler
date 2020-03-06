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
