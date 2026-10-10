
import pandas as pd
import numpy as np
a={'chiense':30,
   'english':45,
   'math':65,
   'chemistry':54,
   'geo':56}

b = pd.Series(a)
print(b)
print(b.iloc[1:3])
b.loc['chiense':'math']
print(b.loc['chiense':'math'])

print(b.index)
print(b.values)
print(b.size)