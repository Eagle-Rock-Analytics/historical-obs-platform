"""
This is a script where Stage 3: QA/QC related common functions, conversions, and operations is stored for ease of use
for the Historical Observations Platform.
"""

## Import Libraries
import boto3
import geopandas as gp
import numpy as np
import pandas as pd
import requests
import urllib
import datetime
import math
import shapely
import xarray as xr
import matplotlib.pyplot as plt
from io import BytesIO, StringIO
import scipy.stats as stats

## Set AWS credentials
s3 = boto3.resource("s3")
s3_cl = boto3.client('s3') # for lower-level processes

## Set relative paths to other folders and objects in repository.
bucket_name = "wecc-historical-wx"
wecc_terr = "s3://wecc-historical-wx/0_maps/WECC_Informational_MarineCoastal_Boundary_land.shp"
wecc_mar = "s3://wecc-historical-wx/0_maps/WECC_Informational_MarineCoastal_Boundary_marine.shp"


## QA/QC helper functions
#-----------------------------------------------------------------------------
def get_file_paths(network):
    rawdir = "1_raw_wx/{}/".format(network)
    cleandir = "2_clean_wx/{}/".format(network)
    qaqcdir = "3_qaqc_wx/{}/".format(network)
    mergedir = "4_merge_wx/{}/".format(network)
    return rawdir, cleandir, qaqcdir, mergedir

#-----------------------------------------------------------------------------
def get_wecc_poly(terrpath, marpath):
    """
    Identifies a bbox of WECC area to filter stations against
    Input vars: shapefiles for maritime and terrestrial WECC boundaries
    Returns: spatial objects for each shapefile, and bounding box for their union.
    """
    t = gp.read_file(terrpath)  ## Read in terrestrial WECC shapefile.
    m = gp.read_file(marpath)   ## Read in marine WECC shapefile.
    bbox = t.union(m).bounds    ## Combine polygons and get bounding box of union.
    return t,m, bbox

#-----------------------------------------------------------------------------
# Log print auxiliary functions
def printf(*args, verbose=True, log_file=None, **kwargs):
    import datetime
    
    tLog = lambda : datetime.datetime.utcnow().strftime("%m-%d-%Y %H:%M:%S") + " : \t"
    args = [str(a) for a in args]
    
    if verbose:
        if log_file is not None:
            print(" ".join([tLog(), *args]), **kwargs) or \
            print(" ".join([tLog(),*args]), file=log_file, **kwargs)
        else:
            print(" ".join([tLog(), *args]), **kwargs)   
    else:
        if log_file is not None:
            print(" ".join([tLog(), *args]), file=log_file, **kwargs)
        else:
            pass

#-----------------------------------------------------------------------------
def create_bins(data, bin_size=0.25):
    '''Create bins from data covering entire data range'''

    # set up bins
    b_min = np.floor(np.nanmin(data))
    b_max = np.ceil(np.nanmax(data))+bin_size
    bins = np.arange(b_min, b_max, bin_size)

    return bins

#-----------------------------------------------------------------------------
def pdf_bounds(df, mu, sigma, bins):
    '''Calculate pdf distribution, return pdf and threshold bounds'''
    y = stats.norm.pdf(bins, mu, sigma)

    # add vertical lines to indicate thresholds where pdf y=0.1
    pdf_bounds = np.argwhere(y > 0.1).squeeze()
    if len(pdf_bounds) == 0:
        printf('PDF distribution warning: there is a bad value present causing issues with pdf y=0.1 determination')
        return (y, int(0), int(len(y)-1)) # returning furthest edge cases, return to in V2
        # return (y, bnds[0], bnds[-1]) # returning furthest edge cases, return to in V2

    else:
        # find first index
        # left_bnd = round(bins[pdf_bounds[0] - 1])
        # right_bnd = round(bins[pdf_bounds[-1] + 1])
        # rounds +1 and -1 is giving out of bounds error, using ceil, floor, and clip instead
        try:
            bnds = np.clip([np.floor(pdf_bounds[0]), np.ceil(pdf_bounds[-1])], 0, len(bins)-1).astype("int")
            bnds = bins[bnds]
            return (y, bnds[0] -1, bnds[-1] + 1)
        except:
            left_bnd = round(bins[pdf_bounds[0] - 1])
            right_bnd = round(bins[pdf_bounds[-1] + 1])
            return (y, left_bnd, right_bnd)

#-----------------------------------------------------------------------------
def qaqc_dist_whole_stn_bypass_check(df, vars_to_check, min_num_months=5):
    """
    Checks the number of valid observation months in order to proceed through monthly distribution checks. 
    Identifies whether a station record has too few months and produces a fail pass flag. 

    Input:
    ------
        df [pd.DataFrame]: station dataset converted to dataframe through QAQC pipeline
        vars_to_check [list]: list of variables to run whole station bypass check on
        min_num_months [int]: minimum number of months required to pass check, default is 5

    Output:
    -------
        df [pd.DataFrame]: QAQC dataframe with flagged values (see below for flag meaning)
        pass_flag [str]: pass flag indicating whether whole station pass/fail number of minimum months

    Flag meaning:
    -------------
        19,Yellow flag,Warning: Whole station has too few monthly observations to proceed through any monthly distribution check. Observations were not assessed for quality
    """
    # This piece will return a dictionary with the var name as key, and values are pd.Series with the
    # month and the number of years of data
    global stn_length
    stn_length = map(qaqc_var_length_bypass_check, [df]*len(vars_to_check), vars_to_check)
    stn_length = {k:v for k,v in zip(vars_to_check, stn_length)}
    
    nYears = np.array([v.max() for k,v in stn_length.items()])

    for var in vars_to_check:
        if stn_length[var].max()<min_num_months: 
            df.loc[:,var+"_eraqc"] = 19  # see era_qaqc_flag_meanings.csv

    return df, stn_length

#-----------------------------------------------------------------------------
def qaqc_dist_var_bypass_check(df, var, min_num_months=5):     
    """
    Checks the number of valid observation months in order to proceed through monthly distribution checks. 
    Identifies whether a station record has too few months and produces a fail pass flag. 

    Input:
    ------
        df [pd.DataFrame]: station dataset converted to dataframe through QAQC pipeline
        var [string]: variable to run bypass check on
        min_num_months [int]: minimum number of months required to pass check, default is 5

    Output:
    -------
        df [pd.DataFrame]: QAQC dataframe with flagged values (see below for flag meaning)
        pass_flag [str]: pass flag indicating whether whole station pass/fail number of minimum months

    Flag meaning:
    -------------
        20,Yellow flag,Warning: Variable has too few monthly non-nan observations to proceed through any monthly distribution check. Observations were not assessed for quality
    """

    df = df.copy()

    # if all values are null for that month across years
    if df[var].isnull().all() == True:
        df[var+'_eraqc'] = 20 # see era_qaqc_flag_meanings.csv

    # if more than min_num_months have invalid obs
    if df.groupby(by=["year","month"])[var]\
         .aggregate("median").transform(np.isnan).sum()>min_num_months:
        df[var+'_eraqc'] = 20 # see era_qaqc_flag_meanings.csv
        
    return df

#-----------------------------------------------------------------------------
def qaqc_var_length_bypass_check(df, var):
    return df.loc[:,[var, "month","year"]].groupby(by=["month"])['year'].unique().apply(len)

#-----------------------------------------------------------------------------
# Red vs. Yellow flagging
def grab_valid_obs(df, var, var2=None, kind='keep'):
    """
    Observations that have been flagged by QA/QC test should not proceed through any
    other QA/QC test. 
    
    Exception is yellow flag (flags: 19, 20), where the entire obs record is too short
    to evaluate for quality. Yellow flag is based on the qaqc_dist_var_bypass_check and
    qaqc_dist_whole_stn_bypass_check checks (in qaqc_unusual_gaps.py) 

    kind [str]: options are "keep" to include 19 and 20 flags (yellow setting)
                       and "drop" to exclude 19 and 20 flags (red setting)
    """

    # grab obs with no flags
    df_noflags = df.loc[df[var+'_eraqc'].isnull() == True]

    # retains yellow flagged obs for QA/QC checks where distribution not assessed
    if kind == 'keep': 
        df_yellowflag = df.loc[(df[var+'_eraqc'] == 19) | (df[var+'_eraqc'] == 20)] # grab obs with yellow flags
        df_merge = pd.concat([df_noflags, df_yellowflag]) # merge together
        df_valid = df_merge.sort_values(by='time') # sort by time

    # yellow flag set to red flag for distibution checks (gaps + climatological outliers) 
    # these obs do not pass through specific QA/QC functions
    elif kind == "drop": # 
        df_valid = df_noflags

    # only applies to some logic checks
    if var2 != None:
        df_valid = df.loc[(df[var+'_eraqc'].isnull() == True) & (df[var2+'_eraqc'].isnull() == True)]

    return df_valid
