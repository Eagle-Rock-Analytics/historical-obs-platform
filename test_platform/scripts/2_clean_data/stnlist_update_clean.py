"""
This script iterates through a specified network and checks to see what stations have been successfully cleaned,
updating the station list in the 1_raw_wx folder to reflect station availability. Error.csvs in the cleaned bucket are also parsed,
with relevant errors added to the corresponding stations if station files are not cleaned, or if the errors occur during or after the cleaning process.

Note that because errors.csv are parsed, very old errors.csv may want to be removed manually from AWS or thresholded below
(removing those produced during code testing)

As of 01/23, current networks are as follows:
ASOSAWOS, CAHYDRO, CDEC, CIMIS, CNRFC, CRN, CW3E, CWOP, HADS, HNXWFO, HOLFUY, HPWREN, LOXWFO, MAP,
MARITIME, MTRWFO, NCAWOS, NDBC, NOS-NWLON, NOS-PORTS, OtherISD, RAWS, SCAN, SGXWFO, SHASAVAL,
SNOTEL, VCAPCD
"""

import boto3
import pandas as pd
from io import BytesIO, StringIO
import numpy as np
import xarray as xr
import s3fs
from datetime import datetime

# Set environment variables
bucket_name = "wecc-historical-wx"
raw_wx = "1_raw_wx/"
clean_wx = "2_clean_wx/"
s3 = boto3.resource("s3")
s3_cl = boto3.client("s3")


# Function: Given a network name, return a pandas dataframe containing the network's station list from the raw bucket.
def get_station_list(network):
    network_prefix = raw_wx+network+"/"

    # If station list is CIMIS, extension is .xlsx
    if network == "CIMIS":
        station_list = f"stationlist_{network}.xlsx"
        obj = s3_cl.get_object(Bucket= bucket_name, Key= network_prefix+station_list)
        station_list = pd.read_excel(BytesIO(obj['Body'].read()))
    elif network == "ASOSAWOS": # Use merged station list
        network_prefix = clean_wx+network+"/"
        station_list = f"stationlist_{network}_merge.csv"
        obj = s3_cl.get_object(Bucket= bucket_name, Key= network_prefix+station_list)
        station_list =  pd.read_csv(obj['Body'])
    else:
        station_list = f"stationlist_{network}.csv"
        obj = s3_cl.get_object(Bucket= bucket_name, Key= network_prefix+station_list)
        station_list = pd.read_csv(obj['Body'])
    return station_list

# Function: Given a network name, return a pandas dataframe of all cleaned stations in the 2_clean_wx AWS bucket, with the date the file was last modified.
def get_cleaned_stations(network):
    df = {'ID':[], 'Time_Cleaned':[]}
    network_prefix = clean_wx+network+"/"
    for item in s3.Bucket(bucket_name).objects.filter(Prefix = network_prefix+network+"_"):
        if network == 'CW3E': # cleaned CW3E data is stored as one file per year per station, handle separately
            stn = item.key.split("_")[3]
            clean_id = network + "_" + stn
            time_mod = item.last_modified
        else: # all other networks, one file per station
            clean_id = item.key.split("/")[-1].replace(".nc", "") # Get ID from file name
            time_mod = item.last_modified
        df['ID'].append(clean_id)
        df['Time_Cleaned'].append(time_mod)
    return pd.DataFrame(df)

# Function: Given a network name, return a pandas dataframe containing all errors reported for the network in the cleaning stage.
def parse_error_csv(network):
    errordf = []
    errors_prefix = clean_wx+network+"/"+"errors"
    for item in s3.Bucket(bucket_name).objects.filter(Prefix = errors_prefix):
        obj = s3_cl.get_object(Bucket= bucket_name, Key= item.key)
        errors = pd.read_csv(obj['Body'])
        if errors.empty:# If file empty
            continue
        else:
            errors = errors[['File', 'Time', 'Error']]
            errordf.append(errors)
    if not errordf: # If no errors in cleaning
        return pd.DataFrame()
    else:
        errordf = pd.concat(errordf)
        errordf = errordf.drop_duplicates(subset = ['File', 'Error'])
        errordf = errordf[errordf.File != "Whole network"] # Drop any whole network errors

        return errordf

# Function: update station list and save to AWS, adding cleaned status, time of clean and any relevant errors
# Optional: clean_var_add will open every cleaned station file, check which variables are present, and flag in the station list
    # Default is False, time intensive
    # Recommendation: run clean_var_add=True only after a full clean or partial clean update
def clean_qa(network, clean_var_add=False):
    if 'otherisd' in network: # Fixing capitalization issues
        network = "OtherISD"
    else:
        network = network.upper()

    # Call functions
    stations = get_station_list(network)
    cleanids = get_cleaned_stations(network)
    errors = parse_error_csv(network)

    if cleanids.empty:
        print("No cleaned files for this network. Please run the relevant cleaning script and try again.")
        exit()

    # Clean station list
    stations = stations.loc[:, ~stations.columns.str.match('Unnamed')] # Drop double index column

    # Add standardized ID column - ERA-ID
    if 'ASOS' in network or 'OtherISD' in network:
        stations['ERA-ID'] = network+"_"+stations['ISD-ID'].str.replace("-", "")
    elif 'CIMIS' in network:
        stations = stations.dropna(subset = ['Station Number']) # Drop one last row that is not standardized
        stations['ERA-ID'] = network+"_"+stations['Station Number'].astype('int').astype('str')
    elif 'CW3E' in network:
        stations['ERA-ID'] = network+"_"+stations['STID'].str.replace("C3", "")
    elif network in ['CAHYDRO', 'CDEC', 'CNRFC', 'CRN', 'CWOP', 'HADS', 'HNXWFO', 'HOLFUY', 'HPWREN', 'LOXWFO',  # MADIS networks
                    'MAP', 'MTRWFO', 'NCAWOS', 'NOS-NWLON', 'NOS-PORTS', 'RAWS', 'SGXWFO', 'SHASAVAL', 'VCAPCD']:
        stations['ERA-ID'] = network+"_"+stations['STID']
    elif network in ['MARITIME', 'NDBC']:
        stations['ERA-ID'] = network+"_"+stations['STATION_ID']
    elif network in ['SCAN', 'SNOTEL']:
        stations['ERA-ID'] = network+"_"+stations['stationTriplet'].str.split(":").str[0]

    # Make ERA-ID first column
    eraid = stations.pop('ERA-ID')
    eraid = eraid.str.upper()  # Standardize names (one outlier station in CWOP)
    stations.insert(0, 'ERA-ID', eraid)

    # Join cleaned columns to column list
    stations = stations.merge(cleanids, left_on = 'ERA-ID', right_on = 'ID', how = 'outer')
    if 'ID_y' in stations.columns:
        stations['Cleaned'] = np.where(stations.ID_y.isna(), "N", "Y") # Make binary cleaned column
        # Drop ID column
        stations = stations.drop(['ID_x', "ID_y"], axis = 1)
    else:
        stations['Cleaned'] = np.where(stations.ID.isna(), "N", "Y") # Make binary cleaned column
        # Drop ID column
        stations = stations.drop("ID", axis = 1)

    # Move Time_Cleaned to last
    s = stations.pop('Time_Cleaned')
    stations = pd.concat([stations, s], axis = 1)

    # Add errors to column by station - only add error if error occurred at or after file clean, if file cleaned.
    stations['Errors'] = np.nan

    # Remove any NAs from ERA-ID
    stations = stations.loc[stations['ERA-ID'].notnull()]

    # If station not in stations list, add it manually
    # NOT YET TESTED - TEST on CRN and NCAWOS after full clean.
    cleanids_list = [x for x in cleanids.ID if any(y in x for y in stations['ERA-ID'])]
    cleanids_nolist = [x for x in cleanids.ID.tolist() if x not in cleanids_list] # Get any ids that aren't in the stations list already
    cleanids_nolist = cleanids.loc[cleanids.ID.isin(cleanids_nolist)]
    if not cleanids_nolist.empty:
        not_in_list = pd.DataFrame({'ERA-ID': cleanids_nolist.ID, 'Cleaned': 'Y', 'Time_Cleaned' : cleanids_nolist.Time_Cleaned})
        stations = pd.concat([stations, not_in_list])

    # Get list of station IDs
    ids = [id.split("_")[-1] for id in stations['ERA-ID'].tolist()]

    if errors.empty: # If no errors, stop here.
        pass

    else:
        # Add relevant ID to errors csv
        errors['ID'] = np.nan
        errors.reset_index(inplace=True, drop = True)
        errors['Time'] = pd.to_datetime(errors['Time'], format = "%Y%m%d%H%M", utc = True)

        for index, row in errors.iterrows():
            id = []
            id = [x for x in ids if x in row['File']]
            if id:
                errors.loc[index, 'ID'] = network+"_"+id[-1]

        for index, row in stations.iterrows(): # For each station
            error_sta = errors.loc[errors.ID == row['ERA-ID']]
            if error_sta.empty: # if no errors for station
                continue
            else:
                if not pd.isnull(row['Time_Cleaned']): # If file cleaned
                    error_sta = error_sta.loc[(error_sta.Time>= row['Time_Cleaned'])|(error_sta.Time.isna()), :] # Only keep errors from cleaning at or after time of clean

                if len(error_sta)==1:
                    stations.loc[index, 'Errors'] = error_sta['Error'].values[0]
                elif len(error_sta)>1:
                    values = [f'{x.File}: {x.Error}' for index, x in error_sta.iterrows()]
                    value = " ".join(values)
                    stations.loc[index, 'Errors'] = value

    # Print summary
    if network != 'CW3E':
        if 'Y' in stations['Cleaned'].values:
            if 'N' not in stations['Cleaned'].values: # order is important here, if no "N" is present in a cleaned network, it will bark without this
                print("Station list updated for cleaned {} stations. All stations cleaned: {} stations cleaned.".format(network, stations['Cleaned'].value_counts()['Y']))
            else:
                print("Station list updated for cleaned {} stations. {} stations cleaned, {} stations not cleaned.".format(network, stations['Cleaned'].value_counts()['Y'], stations['Cleaned'].value_counts()['N']))
        else:
            print("Station list updated for cleaned {} stations. No stations cleaned successfully. {} stations not yet cleaned.".format(network, stations['Cleaned'].value_counts()['N']))
    else: # network is CW3E
        if 'Y' in stations['Cleaned'].values:
            if 'N' not in stations['Cleaned'].values:
                print("Station list updated for cleaned {} stations. All stations cleaned: {} station-years cleaned.".format(network, int(stations['Cleaned'].value_counts()['Y'])))
            else:
                print("Station list updated for cleaned {} stations. {} station-years cleaned, {} station-years not cleaned.".format(network, int(stations['Cleaned'].value_counts()['Y']), stations['Cleaned'].value_counts()['N']))
        else:
            print("Station list updated for cleaned {} stations. No stations cleaned successfully. {} stations not yet cleaned.".format(network, stations['Cleaned'].value_counts()['N']))

    # clean_var_add
    if clean_var_add == True:
        print('Processing all cleaned files to assess variable coverage -- this may take awhile based on size of network!') # useful warning

        # Set up error handling for identifing cleaned files that can't open
        # Appears that some datetime/index failed to format -- need to reclean
        errors = {'File':[], 'Time':[], 'Error':[]} # Set up error handling.
        end_api = datetime.now().strftime('%Y%m%d%H%M') # Set end time to be current time at beginning of download: for error handling csv.
        timestamp = datetime.utcnow().strftime("%m-%d-%Y, %H:%M:%S") # For attributes of netCDF file.


        # add in default columns of "N" to cleaned station list for all core and associated variables
        # also adds column that counts number of valid/non-nan observations
        core_vars = ['tas', 'tdps', 'tdps_derived',
                     'ps', 'psl', 'ps_altimeter', 'ps_derived',
                     'pr', 'pr_5min', 'pr_1h', 'pr_24h', 'pr_localmid',
                     'hurs',  'sfcWind', 'sfcWind_dir', 'rsds']
        for var in core_vars:
            stations[str(var)] = "N"
            stations[str(var+"_nobs")] = 0      # default of 0 to start

        # add column for total length of each record, valid (non-nan) and nans
        stations['total_nobs'] = 0   # default of 0 to start
                    
        # open cleaned datafile
        network_prefix = clean_wx+network+"/"
        files = []
        for item in s3.Bucket(bucket_name).objects.filter(Prefix = network_prefix):
            file = str(item.key)
            files += [file]
            
        # get list of station filenames successfully cleaned    
        files = list(filter(lambda f: f.endswith(".nc"), files))
        
        for file in files: 
            if file not in files: # dont run qa/qc on a station that isn't cleaned
                continue
            else:
                try:
                    print(file)
                    fs = s3fs.S3FileSystem()
                    aws_url = "s3://wecc-historical-wx/"+file

                    with fs.open(aws_url) as fileObj:
                        ds = xr.open_dataset(fileObj) # setting engine=None (default) uses what is best for system, previously engine='h5netcdf'
                        
                        stations.loc[stations['ERA-ID']==ds.station.values[0], 'total_nobs'] = ds.time.shape[0]

                        # mark each variable as present if in dataset, and count number of valid/non-nan values
                        for var in ds.variables:
                            if var in core_vars:
                                stations.loc[stations['ERA-ID']==ds.station.values[0], str(var)] = 'Y'
                                stations.loc[stations['ERA-ID']==ds.station.values[0], str(var+"_nobs")] = ds[str(var)].count()

                        # close dataset
                        ds.close()
                except Exception as e:
                    print('{} not opening'.format(file))
                    errors['File'].append(file) 
                    errors['Time'].append(end_api)
                    errors['Error'].append('clean_var_add error in opening file: {}'.format(e))
                    continue

        # reset index
        stations = stations.reset_index(drop = True)
        
        
    # Save station file to cleaned bucket
    print(stations) # For testing
    new_buffer = StringIO()
    stations.to_csv(new_buffer, index = False)
    content = new_buffer.getvalue()
    s3_cl.put_object(Bucket=bucket_name, Body=content, Key=clean_wx+network+"/stationlist_{}_cleaned.csv".format(network))

    # Save errors file to cleaned bucket
    errors = pd.DataFrame(errors)
    csv_buffer = StringIO()
    errors.to_csv(csv_buffer)
    content = csv_buffer.getvalue()
    s3_cl.put_object(Bucket=bucket_name, Body=content, Key=clean_wx+network+"/add_clean_var_errors_{}_{}.csv".format(network, end_api))


if __name__ == "__main__":
    clean_qa('SHASAVAL', clean_var_add=True)


    # List of all stations for ease of use here:
    # ASOSAWOS, CAHYDRO, CIMIS, CW3E, CDEC, CNRFC, CRN, CWOP, HADS, HNXWFO, HOLFUY, HPWREN, LOXWFO
    # MAP, MTRWFO, NCAWOS, NOS-NWLON, NOS-PORTS, RAWS, SGXWFO, SHASAVAL, VCAPCD, MARITIME
    # NDBC, SCAN, SNOTEL

    # Note: OtherISD only runs as "otherisd"
    # Note: Make sure there is no space in the name CAHYDRO ("CA HYDRO" will not run)
