"""
This script performs qa/qc protocols for cleaned NDBC/MARITIME data for ingestion into the Historical Observations Platform.
Approach:
(1) Remove duplicate stations
(2) Handle variables that report at different intervals and/or change frequency over time (convert to hourly?)
(3) QA/QC testing, including consistency checks, gaps, checks against climatological distributions, and cross variable checks.
(4) Case study analysis for accuracy -- SHOULD THIS BE A SEPARATE SCRIPT/PROCESS?
Inputs: Cleaned data for an individual network
Outputs: QA/QC-processed data for an individual network, priority variables, all times. Organized by station as .nc file.
"""

# Step 0: Environment set-up
# Import libraries
import os
import datetime
# import numpy as np
import pandas as pd
import xarray as xr
import boto3
import s3fs
from io import BytesIO, StringIO


## Import qaqc stage calc functions
try:
    from calc_qaqc import *
except:
    print("Error importing calc_qaqc.py")


## Set up directory to save files temporarily, if it doesn't already exist.
try:
    os.mkdir('temp') # Make the directory to save data in. Except used to pass through code if folder already exists.
except:
    pass


## Set AWS credentials
s3 = boto3.resource("s3")
s3_cl = boto3.client('s3') # for lower-level processes


## Set relative paths to other folders and objects in repository.
bucket_name = "wecc-historical-wx"


# NDBC and MARITIME only
def spurious_buoy_check(station, df, qc_vars):
    """
    Checks the end date on specific buoys to confirm disestablishment/drifting dates of coverage.
    If station reports data past disestablishment date, data records are flagged as suspect.
    If station reports data during buoy dates, data records are flagged as suspect.
    """
    known_issues = ['NDBC_46023', 'NDBC_46045', 'NDBC_46051', 'NDBC_46044', 'MARITIME_PTAC1', 'MARITIME_PTWW1', 'MARITIME_MTYC1', 'MARITIME_MEYC1',
                    'MARITIME_SMOC1', 'MARITIME_ICAC1']
    potential_issues = ['NDBC_46290', 'NDBC_46404', 'NDBC_46212', 'NDBC_46216', 'NDBC_46220', 'NDBC_46226', 'NDBC_46227', 'NDBC_46228', 
                        'NDBC_46230', 'NDBC_46234', 'NDBC_46245', 'NDBC_46250']
                        
    if station in known_issues:
        print('{0} is flagged as suspect, checking data coverage'.format(station)) # testing

        # buoys with "data" past their disestablishment dates
        if station == 'NDBC_46023': # disestablished 9/8/2010
            for i in range(df.shape[0]):
                for j in qc_vars:
                    if (df.index[i][-1].date() >= datetime.date(2010, 9, 9)) == True:
                        df.loc[df.index[i], j] = 2 # see qaqc_flag_meanings.csv
            
        elif station == "NDBC_46045": # disestablished 11/1997
            for i in range(df.shape[0]):
                for j in qc_vars:
                    if (df.index[i][-1].date() >= datetime.date(1997, 12, 1)) == True:
                        df.loc[df.index[i], j] = 2 # see qaqc_flag_meanings.csv

        elif station == "NDBC_46051": # disestablished 4/1996, and out of range of DEM (past coastal range) but reports nan elevation
            # qaqc_elev_infill sets elevation to be assumed 0.0m, but flagging 
            for i in range(df.shape[0]):
                for j in qc_vars:
                    if (df.index[i][-1].date() >= datetime.date(1996, 5, 1)) == True:
                        df.loc[df.index[i], j] = 2 # see qaqc_flag_meanings.csv 
          
        elif station == "MARITIME_PTAC1": # data currently available 1984-2012, but disestablished 2/9/2022
            # only flag if new data is added after 2022 in a new data pull
            for i in range(df.shape[0]):
                for j in qc_vars:
                    if (df.index[i][-1].date() >= datetime.date(2022, 2, 9)) == True:
                        df.loc[df.index[i], j] = 2 # see qaqc_flag_meanings.csv

        # adrift buoy that reports valid data during adrift period (5/2/2015 1040Z to 5/3/2015 1600Z)
        elif station == "NDBC_46044":
            for i in range(df.shape[0]):
                for j in qc_vars:
                    if (df.index[i][-1] >= datetime.datetime(2015, 5, 2, 10, 40, 0)) and (df.index[i][-1] <= datetime.datetime(2015, 5, 3, 15, 50, 0)):
                        df.loc[df.index[i], j] = 2 # see qaqc_flag_meanings.csv

        # other known issues
        elif station == "MARITIME_PTWW1": # wind data obstructed by ferries docking at pier during day hours
            # only wind vars need flag during "day" hours, currently set for 6am to 8pm every day
            for i in range(df.shape[0]):
                if (df.index[i][-1].time() >= datetime.time(6, 0)) and (df.index[i][-1].time() <= datetime.time(20, 0)):
                    df.loc[df.index[i], "sfcWind_eraqc"] = 1 # see qaqc_flag_meanings.csv
                    df.loc[df.index[i], "sfcWind_dir_eraqc"] = 1 

        # elif station == "MARITIME_MTYC1" or station == "MARITIME_MEYC1": # buoy was renamed, no relocation; MTYC1 2005-2016, MEYC1 2016-2021
        #     # modify attribute/naming with note
        #     # this will get flagged in station proximity tests

        # elif station == "MARITIME_SMOC1" or station == "MARITIME_ICAC1": # buoy was renamed, small relocation (see notes); SMOC1 2005-2010, ICAC1 2010-2021
        #     # modify attribute/naming with note
        #     # this will get flagged in station proximity tests
        return df

    elif station in potential_issues: 
        # other stations have partial coverage of their full data records as well as disestablishment dates
        # if new data is added in the future, needs a manual check and added to known issue list if requires handling
        # most of these should be caught by not having a cleaned data file to begin with, so if this print statement occurs it means new raw data was cleaned and added to 2_clean_wx/
        print("{0} has a reported disestablishment date, requires manual confirmation of dates of coverage".format(station))
        for i in range(df.shape[0]):
            for j in qc_vars:
                df.loc[df.index[i], j] = "1"  ## QC FLAG FOR SUSPECT -- using as a placeholder here
    
        return df

    else: # station is not suspicious, move on
        return df


## Function: Conducts whole station qa/qc checks (lat-lon, within WECC, elevation)
def whole_station_qaqc(network, cleandir, qaqcdir):

    # Set up error handling.
    errors = {'File':[], 'Time':[], 'Error':[]} # Set up error handling
    end_api = datetime.datetime.now().strftime('%Y%m%d%H%M') # Set end time to be current time at beginning of download: for error handling csv
    timestamp = datetime.datetime.utcnow().strftime("%m-%d-%Y, %H:%M:%S")

    try:
        files = [] # Get files
        for item in s3.Bucket(bucket_name).objects.filter(Prefix = cleandir):
            file = str(item.key)
            files += [file]

        # Get cleaned station file and read in metadata
        station_file = [file for file in files if 'stationlist_' in file]
        obj = s3_cl.get_object(Bucket=bucket_name, Key=station_file[0])
        station_file = pd.read_csv(BytesIO(obj['Body'].read()))
        stations = station_file['ERA-ID'].dropna()

        files = list(filter(lambda f: f.endswith(".nc"), files)) # Get list of cleaned file names

    except Exception as e:
        print(e) # testing
        errors['File'].append("Whole network")
        errors['Time'].append(end_api)
        errors['Error'].append("Error in whole network: {}".format(e))

    else: # if files successfully read in
        # for station in stations: # full run
        # for station in stations.sample(5): # TESTING SUBSET
        for station in ["NDBC_46051", "NDBC_46089", "NDBC_46023"]: # known issues, 3 have nan elevations = good for testing
            file_name = cleandir+station+".nc"

            if file_name not in files: # dont run qa/qc on a station that isn't cleaned
                print("{} was not cleaned - skipping qa/qc".format(station))
                errors['File'].append(station)
                errors['Time'].append(end_api)
                errors['Error'].append("No cleaned data for this station, does not proceed to qa/qc: see cleaned station list for reason")
                continue
            else:
                print('Running QA/QC on: ', station) # testing

                try:
                    fs = s3fs.S3FileSystem()
                    aws_url = "s3://wecc-historical-wx/"+file_name

                    with fs.open(aws_url) as fileObj:
                        ds = xr.open_dataset(fileObj, engine='h5netcdf') # CHECK THE ENGINE HERE

                        ## Add qc_flag variable for all variables, including elevation; defaulting to nan for fill value that will be replaced with qc flag
                        exclude_qaqc = ["time", "station", "lat", "lon"] # lat and lon have a different qc check
                        raw_qc_vars = [] # qc_variable for each data variable, will vary station to station
                        era_qc_vars = [] # our qc variable
                        for var in ds.variables:
                             if 'q_code' in var:
                                  raw_qc_vars.append(var) # raw qc variable, need to keep for comparison, then drop

                        for var in ds.variables:
                            if var not in exclude_qaqc and var not in raw_qc_vars:
                                qc_var = var + "_eraqc" # variable/column label
                                era_qc_vars.append(qc_var)
                                ds[qc_var] = xr.full_like(ds[var], np.nan) # adds new variable in shape of original variable with designated nan fill value
                                                
                        stn_to_qaqc = ds.to_dataframe()

                        ## QAQC Functions
                        ## Lat-lon -- does not proceed through qaqc if failure
                        stn_to_qaqc = qaqc_missing_latlon(stn_to_qaqc)
                        if len(stn_to_qaqc.index) == 0:
                            print('{0} has a missing lat-lon, skipping'.format(station)) # testing
                            errors['File'].append(station)
                            errors['Time'].append(end_api)
                            errors['Error'].append('Failure on qaqc_missing_latlon')
                            continue # skipping station
                        print('pass qaqc_missing_latlon') #testing

                        ## Within WECC -- does not proceed through qaqc if failure
                        stn_to_qaqc = qaqc_within_wecc(stn_to_qaqc)
                        if len(stn_to_qaqc.index) == 0:
                            print('{0} lat-lon is out of range for WECC, skipping'.format(station)) # testing
                            errors['File'].append(station)
                            errors['Time'].append(end_api)
                            errors['Error'].append('Failure on qaqc_within_wecc')
                            continue # skipping station
                        print('pass qaqc_within_wecc') #testing

                        ## Buoys with known issues with specific qaqc flags
                        era_qc_vars.remove("elevation_eraqc") # remove elevation_qc var from remainder of analyses so it does not also get flagged -- confirm with final qaqc process
                        try:
                            stn_to_qaqc = spurious_buoy_check(station, stn_to_qaqc, era_qc_vars)

                        except Exception as e:
                            print('Flagging problematic buoy issue for {0}, skipping'.format(station)) # testing
                            errors['File'].append(station)
                            errors['Time'].append(end_api)
                            errors['Error'].append('Failure on spurious_buoy_check: {0}'.format(e))
                            continue # skipping station
                        print('pass spurious_buoy_check') #testing

                        ## Elevation -- if DEM in-filling fails, does not proceed through qaqc
                        stn_to_qaqc = qaqc_elev_infill(stn_to_qaqc) # nan infilling must be before range check
                        if len(stn_to_qaqc.index) == 0:
                            print('DEM in-filling for {} failed, may not mean station does not pass qa/qc -- check'.format(station)) # testing
                            errors['File'].append(station)
                            errors['Time'].append(end_api)
                            errors['Error'].append('DEM in-filling error, may not mean station does not pass qa/qc -- check')
                            continue # skipping station
                        print('Elevation values: {}'.format(stn_to_qaqc['elevation'].unique())) # testing
                        print('Elevation eraqc values: {}'.format(stn_to_qaqc['elevation_eraqc'].unique())) # testing

                        stn_to_qaqc = qaqc_elev_range(stn_to_qaqc)
                        if len(stn_to_qaqc.index) == 0:
                            print('{} elevation out of range for WECC, skipping'.format(station)) # testing
                            errors['File'].append(station)
                            errors['Time'].append(end_api)
                            errors['Error'].append('Failure on qaqc_elev_range')
                            continue
                        print('pass qaqc_elev_range') # testing

                        print(stn_to_qaqc.head(20))


                except Exception as e:
                    # print(e) # testing
                    errors['File'].append(station)
                    errors['Time'].append(end_api)
                    errors['Error'].append("Cannot read files in from AWS: {0}".format(e))
                        


                # last step is to reassign df back to xarray object
                # Sort by time and remove any overlapping timesteps
                # stn_to_qaqc = stn_to_qaqc.sort_values(by='time')
                # stn_to_qaqc = stn_to_qaqc.drop_duplicates()
                # ds = stn_to_qaqc.to_xarray()

                # Update global attributes
                ds = ds.assign_attrs(title = network+" quality controlled")
                ds = ds.assign_attrs(history = 'MARITIME_qaqc.py script run on {} UTC'.format(timestamp))
                ds = ds.assign_attrs(comment = 'Intermediate data product: subject to cleaning but may not be subject to full QA/QC processing.')
                ## need to reassign attributes from cleaning stage here? -- check

                # # Write station file to netcdf format
                # try:
                #     filename = station + ".nc" # Make file name
                #     filepath = qaqcdir + filename # Writes file path

                #     # Write locally
                #     ds.to_netcdf(path = 'temp/temp.nc', engine = 'netcdf4') # Save station file.

                #     # Push file to AWS with correct file name
                #     s3.Bucket(bucket_name).upload_file('temp/temp.nc', filepath)

                #     print('Saving {} with dims {}'.format(filename, ds.dims))
                #     ds.close() # Close dataframe

                # except Exception as e:
                #     print(filename, e)
                #     errors['File'].append(filename)
                #     errors['Time'].append(end_api)
                #     errors['Error'].append('Error saving ds as .nc file to AWS bucket: {}'.format(e))
                #     continue

    # Write errors to csv
    finally:
        errors = pd.DataFrame(errors)
        csv_buffer = StringIO()
        errors.to_csv(csv_buffer)
        content = csv_buffer.getvalue()
        s3_cl.put_object(Bucket=bucket_name, Body=content, Key=qaqcdir+"errors_{}_{}.csv".format(network, end_api)) # Make sure error files save to correct directory


## -------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Run function
if __name__ == "__main__":
    network = "NDBC"
    rawdir, cleandir, qaqcdir, mergedir = get_file_paths(network)
    whole_station_qaqc(network, cleandir, qaqcdir)

# To do:
# reorder variables once entire qaqc is complete before saving
# output csv of flags/consistent flagging
# check the h5netcdf vs. netcdf4 engine
# delete testing notes
