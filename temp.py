import pandas as pd
import numpy as np
df0 = pd.DataFrame({"coordinateUncertaintyInMeters": ['UNKNOWN', '4', '6', '300']})
print(df0)
df0['coordinateUncertaintyInMeters'].replace(to_replace='4', value=np.NaN, inplace=True)
print(df0)
df0['coordinateUncertaintyInMeters'].replace(to_replace='UNKNOWN', value='99', inplace=True)
print(df0)
