import pandas as pd
import numpy as np

def var_to_unique_list(ds, column):
    flagvals = ds[column].values.tolist()[0]
    flagvals = [x for x in flagvals if pd.isnull(x) == False] # Remove nas
    flagvals = list(np.unique(flagvals)) # Get unique values
    flagvals = " ".join(flagvals)
    return flagvals

