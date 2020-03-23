# -*- coding: utf-8 -*-
"""
Created on Thu Mar 19 14:56:25 2020

@author: nmtarr

Method: The first df is cleaned up by dropping duplicates based on which 
record has greater individual count.  Before doing that, records with unequal
decimal precision in the lat and long fields and those fields are truncated
to the shorter precision present.

An input df likely contains records with equal decimal precision in lat and 
long fields, but that is lower than the rest (i.e. lat and long have 3 places
right of the decimal whereas most records have 4).  Duplication may occur
between lower and higher precision records at the lower precision.  Therefore,
duplication must be assessed at each of the lower precision levels present.  
The strategy for that is to, at each precision level, split the main df in two:
one with records having the precision level of the investigation and another 
with records greater than the precision level. The "greater than" df records'
lat and long values are then truncated to the precision level.  Records are 
identified from the "equals precision" df that have their lat, long, and date
values represented in the "greater than" df, and such records id's are
collected in a list of records to remove from the input/main df.  This process
is iterated over all precision levels present in the data.

# Test df to use for development/troubleshooting.
df = pd.DataFrame({'latitude': ["10.1234", "10.123", "10.123", "10.12", 
                                "10.2345", "10.2345", "10.234", "10.23", 
                                "10.3456", "10.3456", "10.3456", "10.777",
                                "10.9"],
                   'longitude': ["-12.1234", "-12.123", "-12.123", "-12.12", 
                                 "-12.2345", "-12.2345", "-12.234", "-12.23",
                                 "-12.3456", "-12.3456", "-12.34567", 
                                 "-12.777", "-12.9"],
                   'occ_id': ["a", "b", "c", "d", "e", "f", "g", "h", "i",
                              "j", "l", "m", "n"],
                   'occurrenceDate': [1,1,1,1,1,1,1,1,1,1,1,1,1],
                   'individualCount': [1,2,3,3,5,6,7,8,9,10,11,12,13]})
"""
import pandas as pd#########
from datetime import datetime#############
import decimal

startduptime = datetime.now()
pd.set_option('display.max_columns', 20)######################
pd.set_option('display.max_rows', 100)######################
pd.set_option('display.width', 600)######################

df = pd.read_csv("T:/temp/df8.csv", dtype={'latitude': str, 'longitude': str})
df9 = df.copy()

def drop_duplicates_latlongdate(df):
    # Record df length before removing duplicates
    initial_length = len(df)
    
    """
    ############ RECTIFY UNEQUAL LAT-LONG PRECISION
    First, trim decimal length in cases where demical length differs between
    lat and long, result is equal lat and long length.  Record the trimmed
    decimal precision in a temp column for use later as a record to "native"
    precision.
    """
    df['dup_latPlaces'] = [len(x.split(".")[1]) for x in df['latitude']]
    df['dup_lonPlaces'] = [len(x.split(".")[1]) for x in df['longitude']]
    df['dup_OGprec'] = df['dup_latPlaces']
    df22 = df[df['dup_latPlaces'] != df['dup_lonPlaces']]
    for i in df22.index:
        x = df22.loc[i]
        if x['dup_latPlaces'] < x['dup_lonPlaces']:
            trim_len = int(x['dup_latPlaces'])
        else:
            trim_len = int(x['dup_lonPlaces'])
        df.loc[i, 'latitude'] = x['latitude'][:trim_len + 3]
        df.loc[i, 'longitude'] = x['longitude'][:trim_len + 4]
        # Record the resulting precision for reference later
        df.loc[i, 'dup_OGprec'] = trim_len
    df.drop(['dup_latPlaces', 'dup_lonPlaces'], axis=1, inplace=True)
    
    """
    ########  INITIAL DROP OF DUPLICATES
    Initial drop of duplicates on 'latitude', 'longitude', 'occurrenceDate', 
    keeping the first (highest individual count)
    Sort so that the highest individual count is first ############ ADD OCCURRENCEDATE BACK IN
    """
    df.sort_values(by=['latitude', 'longitude', 
                        'individualCount'],
                    ascending=False, inplace=True, kind='mergesort', 
                    na_position='last')
    
    df.drop_duplicates(subset=['latitude', 'longitude'], 
                        keep='first', inplace=True)
    print(str(initial_length - len(df)) + " duplicate records dropped: {0}".format(startduptime))
    """
    #########  FIND IMPRECISE DUPLICATES
    Get a list of "native" precisions that are present in the data to loop through.
    Next, iterate through this list collecting id's of records that need to be 
    removed from the main df.
    """
    # Get list of unique precisions.  Order is important: descending.
    precisions = list(set(df['dup_OGprec']))
    print(precisions)
    precisions.sort(reverse=True)
    # The highest precisions listed at this point has already been done: drop it.
    precisions = precisions[1:]
    
    # List for collecting records that are duplicates
    duplis = []
    
    # The precision-specific duplicate testing happens repeatedly, so make it a 
    # function.  
    def drop_duplicates(precision, df):
        """
        Function to find undesirable duplicates at a particular decimal precision.
        
        The general strategy is to split the input df in two. 
        
        Parameters:
        precision -- The level of precision (places right of decimal) in lat and long
                     for the assessment of duplicates.
        df -- dataframe to assess and drop duplicates from.  This function works
              'inplace'.
        """
        print("\nDecimal precision = " + str(precision))
        # Create a df with records from the input df having decimal precision > the
        # precision level being assessed.
        dfx = df.copy()
        dfLonger = dfx[dfx['dup_OGprec'] > precision].copy()
        # Truncate lat and long values
        dfLonger['latitude'] = [x[:precision + 3] for x in dfLonger['latitude']]
        dfLonger['longitude'] = [x[:precision + 4] for x in dfLonger['longitude']]
        
        # Create a df with records having the precision being 
        # investigated
        dfShorter1 = dfx[dfx['dup_OGprec'] == precision]
    
        # Find records in dfShorter1 with lat, lon ('date') combo 
        # existing in dfLonger and append to list of duplis
        dfduplis = pd.merge(dfShorter1, dfLonger, how='inner', 
                            on=['latitude', 'longitude', 'occurrenceDate'])
        dups_ids = dfduplis['occ_id_x']
        for d in dups_ids:
            duplis.append(d)
            
    # Drop lat long duplicates at lower decimal precisions
    for p in precisions:
        drop_duplicates(p, df)
        
    print('\nduplis')
    print(set(duplis))
    print('\n')
    
    # Drop rows from the current main df that have been identified as duplicates.
    df2 = df[df['occ_id'].isin(duplis) == False].copy()
    
    # Drop excess columns
    df2.drop(['dup_OGprec'], inplace=True, axis=1)
    
    duptime = datetime.now() - startduptime
    print(str(initial_length - len(df2)) + " duplicate records dropped: {0}".format(duptime))
    return df2

drop_duplicates_latlongdate(df9)