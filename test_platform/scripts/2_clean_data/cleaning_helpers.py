import pandas as pd
import numpy as np
import boto3

s3 = boto3.resource('s3')

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


# Given a network name, return all relevant AWS filepaths for other functions.
def get_file_paths(network):
    rawdir = "1_raw_wx/{}/".format(network)
    cleandir = "2_clean_wx/{}/".format(network)
    qaqcdir = "3_qaqc_wx/{}/".format(network)
    return rawdir, cleandir, qaqcdir


# function can be renamed to be "delete files from AWS" in future
def rename_MARITIME(bucket_name, cleandir): 
    """
    Ancilliary one-time help to rename MARITIME cleaned files which was previously saved with lowercase station names. 
    Will be fixed in future cleaning update for MARITIME network. 

    *** NOTE ***
    Function moved into cleaning_helpers.py as may be helpful in future for automation/delete empty files from Synoptic update
    in a modified format (eg the delete line, and constructed to look for small file sizes).
    """
    
    # source: https://stackoverflow.com/questions/32501995/boto3-s3-renaming-an-object-using-copy-object

    files = [] # Get files in MARITIME
    for item in s3.Bucket(bucket_name).objects.filter(Prefix = cleandir):
        file = str(item.key)
        files += [file]

    files = list(filter(lambda f: f.endswith(".nc"), files)) # Get list of cleaned file names

    for file in files:

        # example old file name: 2_clean_wx/MARITIME/MARITIME_wptw1.nc
        file_parts = file.split("/")
        stn_name = file_parts[-1]
        stn_name = stn_name.split(".")[0] # remove file extension

        # example new file name: 2_clean_wx/MARITIME/MARITIME_WPTW1.nc
        file_upper = cleandir + stn_name.upper() + ".nc"

        s3.Object(bucket_name, file_upper).copy_from(CopySource=bucket_name+'/'+file) # rename to upper case
        print('{0} renamed to {1}'.format(file, file_upper))

        s3.Object(bucket_name, file).delete() # delete lowercase name files
        print('Old file {} deleted from AWS'.format(file))