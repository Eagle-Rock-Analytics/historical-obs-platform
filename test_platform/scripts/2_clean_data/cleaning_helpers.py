import pandas as pd
import numpy as np

# Given a variable column in an xarray object, this function returns a string with all of the unique variables
# present in that column. This is used to generate lists of qaqc flag values from existing data in the cleaning stage.
# Inputs: ds is xarray dataset, column is the name of the column
# Outputs: flagvals, a string that can be provided as the flag_values attribute for a QA/QC flag.
def var_to_unique_list(ds, column):
    flagvals = ds[column].values.tolist()[0]
    flagvals = [x for x in flagvals if pd.isnull(x) == False] # Remove nas
    flagvals = list(np.unique(flagvals)) # Get unique values
    flagvals = " ".join(flagvals)
    return flagvals

