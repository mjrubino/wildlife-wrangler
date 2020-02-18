import pandas as pd
file = "T:/temp/dfOcc.csv"

df0 = pd.read_csv(file, nrows=1000)
print(df0)


# Count entries per atrribute(column), reformat as new df with appropriate
# columns.  Finally, insert into db.
df_populated1 = pd.DataFrame(df0.count(axis=0).T.iloc[1:])
df_populated1['included(n)'] = len(df0)
df_populated1['populated(n)'] = df_populated1[0]
df_populated2 = df_populated1.filter(items=['included(n)', 'populated(n)'], axis='columns')
df_populated2.index.name = 'attribute'
print(df_populated2)
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
print(df0.columns)
for x in summary.keys():
    print(x)
    print(set(df0[x]))
