# -*- coding: utf-8 -*-
"""
Created on Thu Mar 19 14:56:25 2020

@author: nmtarr

FInd and delete duplicates.  

Find duplicates at various decimal precisions.

What about when decimal precision isn't equal?

Method: Duplicates are droppped from the first df based on which has greater 
individual count.  Then for each number of decimal places present from lat and
long columns
"""
import pandas as pd
import decimal

pd.set_option('display.max_columns', 20)
pd.set_option('display.max_rows', 100)
pd.set_option('display.width', 600)

df = pd.DataFrame({'latitude': ["10.1234", "10.123", "10.123", "10.12", "10.2345", "10.2345", "10.234", "10.23", "10.3456", "10.3456", "10.3456"],
                   'longitude': ["-12.1234", "-12.123", "-12.123", "-12.12", "-12.2345", "-12.2345", "-12.234", "-12.23", "-12.3456", "-12.3456", "-12.34567"],
                   'occ_id': ["a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "l"],
                   'individualCount': [1,2,3,3,5,6,7,8,9,10,11]})

#df = pd.read_csv("T:/temp/df8.csv", dtype={'latitude': str, 'longitude': str})
df8 = df
print('raw')
print(df8)
# Record df length before removing duplicates
initial_length = len(df8)
print(initial_length)
print('\n')


df9 = df8.copy()
"""
############ RECTIFY UNEQUAL LAT-LONG PRECISION
First, trim decimal length in cases where demical length differs between
lat and long, result is equal lat and long length
"""
df9['dup_latPlaces'] = [len(x.split(".")[1]) for x in df9['latitude']]
df9['dup_lonPlaces'] = [len(x.split(".")[1]) for x in df9['longitude']]
df9['dup_OGprec'] = df9['dup_latPlaces']

df10 = df9[df9['dup_latPlaces'] != df9['dup_lonPlaces']]
for i in df10.index:
    x = df10.loc[i]
    if x['dup_latPlaces'] < x['dup_lonPlaces']:
        trim_len = int(x['dup_latPlaces'])
    else:
        trim_len = int(x['dup_lonPlaces'])
    df9.loc[i, 'latitude'] = x['latitude'][:trim_len + 3]
    df9.loc[i, 'longitude'] = x['longitude'][:trim_len + 4]
    # Record the resulting precision for reference later
    df9.loc[i, 'dup_OGprec'] = trim_len
df9.drop(['dup_latPlaces', 'dup_lonPlaces'], axis=1, inplace=True)
print('trimmed\n')

##### TEMP PRINT
df9.sort_values(by='occ_id', ascending=True, inplace=True, 
                kind='mergesort')
print(df9)
print('\n')

"""
########  INITIAL DROP OF DUPLICATES
Initial drop of duplicates on 'latitude', 'longitude', 'occurrenceDate', 
keeping the first (highest individual count)
Sort so that the highest individual count is first ############ ADD OCCURRENCEDATE BACK IN
"""
df9.sort_values(by=['latitude', 'longitude', 
                    'individualCount'],
                ascending=False, inplace=True, kind='mergesort', 
                na_position='last')

df9.drop_duplicates(subset=['latitude', 'longitude'], 
                    keep='first', inplace=True)
####### TEMP PRINT
print('duplicate drop 1\n')
df9.sort_values(by='occ_id', ascending=True, inplace=True, 
                kind='mergesort')
print(df9)

"""
#########  FIND INPRECISE DUPLICATES


"""
# Duplicates now must be assessed at each decimal precision present in the data
# so get a list of unique precisions.  Order is important: descending.
precisions = list(set(df9['dup_OGprec']))
precisions.sort(reverse=True)
# The highest precisions listed at this point has already been done: drop it.
precisions = precisions[1:]

winners = []
weeds = []

precision = 3
df = df9
#def drop_duplicates(precision, df, winners):
"""
Function to find undesirable duplicates at a particular decimal precision.

The general strategy is to split the input df in two.  
"""
print("\nDecimal precision = " + str(precision))
# Create a df with records from the input df having decimal precision > the
# precision level being assessed.
dfx = df.copy()
dfLonger = dfx[dfx['dup_OGprec'] > precision]
print(dfx)
print(dfLonger)
dfLonger['latitude'] = [x[:precision + 3] for x in dfLonger['latitude']]
dfLonger['longitude'] = [x[:precision + 4] for x in dfLonger['longitude']]
print(dfLonger)


# Create a df with records having the precision being investigated
dfShorter1 = dfx[dfx['dup_OGprec'] == precision]
print('\n')
print(df8)
print(dfLonger)
print(dfShorter1)

# Find records in dfShorter1 with lat, lon ('date') combo existing in dfLonger
# and append to list of weeds


"""

dfx['latitude'] = [x[:precision + 3] for x in dfx['latitude']]
dfx['longitude'] =[x[:precision + 4] for x in dfx['longitude']]

# Drop duplicates, keeping the one with the largest original decimal precision
dfx.sort_values(by=['latitude', 'longitude', 
                    'dup_OGprec'],
                ascending=False, inplace=True, kind='mergesort', 
                na_position='last')

dfx.drop_duplicates(subset=['latitude', 'longitude'], 
                    keep='first', inplace=True)
ids = list(dfx['occ_id'])
for k in ids:
    winners.append(k)
print('\nduplicate drop ' + str(precision) + '\n')
dfx.sort_values(by='occ_id', ascending=True, inplace=True, 
            kind='mergesort')
print(dfx)
#return dfx

for i in precisions:
    drop_duplicates(i, df9, winners)

# Remove duplicate occ_ids from the keeper (non-duplicates) list
winners = list(set(winners))
winners.sort()

# Keep only records from first df that are in the list of winners
df11 = df9[df9['occ_id'].isin(winners) == True]

# Now drop records that have a higher precision equivalent.


# Drop excess columns
df11.drop(['dup_OGprec'], inplace=True, axis=1)

print('\n')
print(df11)
print('\n')
print(winners)
"""