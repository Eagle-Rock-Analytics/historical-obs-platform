"""
This script performs qa/qc protocols for cleaned station data for ingestion into the Historical Observations Platform, and is
independent of network. 
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
import pandas as pd
import xarray as xr
import boto3
import s3fs
from io import BytesIO, StringIO
import argparse

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
        for station in stations.sample(1): # TESTING SUBSET
        # for station in ['ASOSAWOS_72676324198']: # this is the smallest ASOSAWOS file and takes ~10 seconds to run for Victoria
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
                        ds = xr.open_dataset(fileObj) # CHECK THE ENGINE HERE -- setting to default which operates on best with dependencies, previously 'h5netcdf'

                        ## Add qc_flag variable for all variables, including elevation; defaulting to nan for fill value that will be replaced with qc flag
                        exclude_qaqc = ["time", "station", "lat", "lon", "qaqc_process", "sfcWind_method"] # lat and lon have a different qc check
                        raw_qc_vars = [] # qc_variable for each data variable, will vary station to station
                        era_qc_vars = [] # our qc variable
                        for var in ds.variables:
                            if 'q_code' in var:
                                raw_qc_vars.append(var) # raw qc variable, need to keep for comparison, then drop
                            if '_qc' in var:
                                raw_qc_vars.append(var) # raw qc variables, need to keep for comparison, then drop

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


                        ## Elevation -- if DEM in-filling fails, does not proceed through qaqc
                        stn_to_qaqc = qaqc_elev_infill(stn_to_qaqc) # nan infilling must be before range check
                        if len(stn_to_qaqc.index) == 0:
                            print('DEM in-filling for {} failed, may not mean station does not pass qa/qc -- check'.format(station)) # testing
                            errors['File'].append(station)
                            errors['Time'].append(end_api)
                            errors['Error'].append('DEM in-filling error, may not mean station does not pass qa/qc -- check')
                            continue # skipping station

                        stn_to_qaqc = qaqc_elev_range(stn_to_qaqc)
                        if len(stn_to_qaqc.index) == 0:
                            print('{} elevation out of range for WECC, skipping'.format(station)) # testing
                            errors['File'].append(station)
                            errors['Time'].append(end_api)
                            errors['Error'].append('Failure on qaqc_elev_range')
                            continue # skipping station
                        print('pass qaqc_elev_range') # testing
                        
                        
                        ## Missing value checks -- convert any lingering missing value codes to NaNs
                        try:
                            stn_to_qaqc = qaqc_missing_vals(stn_to_qaqc)
                            
                        except:
                            print('Flagging problem on missing values check for {0}, skipping'.format(station)) # testing
                            errors['File'].append(station)
                            errors['Time'].append(end_api)
                            errors['Error'].append('Failure on qaqc_missing_vals')
                            continue # skipping station
                        print('pass qaqc_missing_vals') # testing


                        ## Variable logic checks
                        # precipitation is not negative
                        try:
                            stn_to_qaqc = qaqc_precip_logic_nonegvals(stn_to_qaqc)
                        except Exception as e:
                            print('Flagging problem with negative precipitation values for {0}, skipping'.format(station)) # testing
                            errors['File'].append(station)
                            errors['Time'].append(end_api)
                            errors['Error'].append('Failure on qaqc_precip_logic_nonegvals: {0}'.format(e))
                            continue # skipping station
                        print('pass qaqc_precip_logic_nonegvals') # testing

                        ## precipitation duration logic
                        try:
                            stn_to_qaqc = qaqc_precip_logic_accum_amounts(stn_to_qaqc)
                        except Exception as e:
                            print('Flagging problem with precip duration logic check for {0}, skipping'.format(station)) # testing
                            errors['File'].append(station)
                            errors['Time'].append(end_api)
                            errors['Error'].append('Failure on qaqc_precip_logic_accum_ammounts: {0}'.format(e))
                            continue # skipping station
                        print('pass qaqc_precip_logic_accum_amounts') # testing

                        ## Buoys with known issues with specific qaqc flags
                        if network == 'MARITIME' or network == 'NDBC':
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
                            
                        ## Sensor height: wind
                        try:
                            stn_to_qaqc = qaqc_sensor_height_w(ds, stn_to_qaqc)
                        except Exception as e:
                            print('Flagging problem with anemometer sensor height for {0}, skipping'.format(station)) # testing
                            errors['File'].append(station)
                            errors['Time'].append(end_api)
                            errors['Error'].append('Failure on qaqc_sensor_height_w: {0}'.format(e))
                            continue # skipping station
                        print('complete qaqc_sensor_height_w')
                        
                        ## Sensor height: air temperature
                        try:
                            stn_to_qaqc = qaqc_sensor_height_t(ds, stn_to_qaqc)
                        except Exception as e:
                            print('Flagging problem with thermometer sensor height for {0}, skipping'.format(station)) # testing
                            errors['File'].append(station)
                            errors['Time'].append(end_api)
                            errors['Error'].append('Failure on qaqc_sensor_height_t: {0}'.format(e))
                            continue # skipping station
                        print('complete qaqc_sensor_height_t')
                        
                        ## World record checks: air temperature, dewpoint, wind, pressure
                        try:
                            stn_to_qaqc = qaqc_world_record(stn_to_qaqc)
                        except Exception as e:
                            print('Flagging problem with world record check for {0}, skipping'.format(station)) # testing
                            errors['File'].append(station)
                            errors['Time'].append(end_api)
                            errors['Error'].append('Failure on qaqc_world_record: {0}'.format(e))
                            continue # skipping station
                        print('complete qaqc_world_record')
                        
                        print(stn_to_qaqc.head(10)) # testing


                        ## Variable cross-logic checks
                        # dew point temp cannot exceed air temperature
                        try:
                            stn_to_qaqc = qaqc_crossvar_logic_tdps_to_tas(stn_to_qaqc)
                        except Exception as e:
                            print('Flagging problem with temperature cross-variable logic check for {0}, skipping'.format(station)) # testing
                            errors['File'].append(station)
                            errors['Time'].append(end_api)
                            errors['Error'].append('Failure on qaqc_crossvar_logic_tdps_to_tas: {0}'.format(e))
                            continue # skipping station
                        print('pass qaqc_crossvar_logic_tdps_to_tas') # testing

                        # wind direction should be 0 when wind speed is also 0
                        try:
                            stn_to_qaqc = qaqc_crossvar_logic_calm_wind_dir(stn_to_qaqc)
                        except Exception as e:
                            print('Flagging problem with wind cross-variable logic check for {0}, skipping'.format(station)) # testing
                            errors['File'].append(station)
                            errors['Time'].append(end_api)
                            errors['Error'].append('Failure on qaqc_crossvar_logic_calm_wind_dir: {0}'.format(e))
                            continue # skipping station
                        print('pass qaqc_crossvar_logic_calm_wind_dir') # testing


                except Exception as e:
                    print(e) # testing
                    errors['File'].append(station)
                    errors['Time'].append(end_api)
                    errors['Error'].append("Cannot read files in from AWS: {0}".format(e))
                        


                # last step is to reassign df back to xarray object
                # Sort by time and remove any overlapping timesteps
                # stn_to_qaqc = stn_to_qaqc.sort_values(by='time')
                # stn_to_qaqc = stn_to_qaqc.drop_duplicates()
                # ds = stn_to_qaqc.to_xarray()

                # Update global attributes
                # ds = ds.assign_attrs(title = network+" quality controlled")
                # ds = ds.assign_attrs(history = 'ALLNETWORKS_qaqc.py script run on {} UTC'.format(timestamp))
                # ds = ds.assign_attrs(comment = 'Intermediate data product: subject to cleaning but may not be subject to full QA/QC processing.')
                ## need to reassign attributes from cleaning stage here? -- check

                # # Write station file to netcdf format
                # try:
                #     filename = station + ".nc" # Make file name
                #     filepath = qaqcdir + filename # Writes file path

                #     # Write locally
                #     ds.to_netcdf(path = 'temp/temp.nc', engine = 'netcdf4') # Save station file.

                #     # Push file to AWS with correct file name
                #     s3.Bucket(bucket_name).upload_file('temp/temp.nc', filepath)

                #     print('Saving {0} with dims {1} to {2}'.format(filename, ds.dims, bucket_name+"/"+qaqcdir))
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
        print('errors_{0}_{1}.csv saved to {2}'.format(network, end_api, bucket_name + "/" + qaqcdir))


## -------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Run function
if __name__ == "__main__":
    
    # Create parser
    parser = argparse.ArgumentParser(
        prog="ALLNETWORKS_qaqc",
        description="""This script performs qa/qc protocols for cleaned station data for ingestion into the Historical 
                       observations Platform, and is independent of network.""",
        epilog="""Possible stations:
                  [ASOSAWOS, CAHYDRO, CIMIS, CW3E, CDEC, CNRFC, CRN, CWOP, HADS, HNXWFO, 
                   HOLFUY, HPWREN, LOXWFOMAP, MTRWFO, NCAWOS, NOS-NWLON, NOS-PORTS, OtherISD, 
                   RAWS, SGXWFO, SHASAVAL, VCAPCD, MARITIME, NDBC, SCAN, SNOTEL]
               """
    )
    
    # Define arguments for the program
    parser.add_argument('-n', '--network', default="VCAPCD")
    parser.add_argument('-v', '--verbose', default="True")
    
    # Parse arguments
    args = parser.parse_args()
    network = args.network

    rawdir, cleandir, qaqcdir, mergedir = get_file_paths(network)
    whole_station_qaqc(network, cleandir, qaqcdir)

# Dev to do:
# reorder variables once entire qaqc is complete before saving
# output csv of flags/consistent flagging
# check the h5netcdf vs. netcdf4 engine
# delete testing notes

