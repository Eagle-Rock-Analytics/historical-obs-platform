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
import time
import tempfile

# Import qaqc stage calc functions
try:
    from calc_qaqc_xarray import *
except:
    print("Error importing calc_qaqc.py")

# Set up directory to save files temporarily and save timing, if it doesn't already exist.
dirs = ["./temp/", "./timing/"]
for d in dirs:
    if not os.path.exists(d):
        os.makedirs(d)

#----------------------------------------------------------------------------
## Set AWS credentials
s3 = boto3.resource("s3")
s3_cl = boto3.client('s3') # for lower-level processes

## Set relative paths to other folders and objects in repository.
bucket_name = "wecc-historical-wx"

#============================================================================
# Define global functions and variables

#----------------------------------------------------------------------------
def setup_error_handling():
    """
    """
    errors = {'File':[], 'Time':[], 'Error':[]} # Set up error handling
    end_api = datetime.datetime.now().strftime('%Y%m%d%H%M') # Set end time to be current time at beginning of download: for error handling csv
    timestamp = datetime.datetime.utcnow().strftime("%m-%d-%Y, %H:%M:%S")
    return errors, end_api, timestamp
    
#----------------------------------------------------------------------------
def print_qaqc_failed(errors, station=None, end_api=None,
                      message=None, test=None, verbose=True):
    """
    """
    if verbose:
        print('{0} {1}, skipping'.format(station, message)) # testing
    errors['File'].append(station)
    errors['Time'].append(end_api)
    errors['Error'].append('Failure on {}'.format(test))
    return errors

#----------------------------------------------------------------------------
## Read network nc files
def read_network_files(network, cleandir):
    """
    """
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

    return files, stations

#----------------------------------------------------------------------------
## Assign ds attributes and save
def process_output_ds(ds, network, timestamp, station, qaqcdir,
                      errors, end_api, verbose=True):
    """
    """
    ds = ds.assign_attrs(title = network + " quality controlled")
    ds = ds.assign_attrs(history = 'ALLNETWORKS_qaqc.py script run on {} UTC'.format(timestamp))
    ds = ds.assign_attrs(comment = 'Intermediate data product: subject to cleaning but may not be subject to full QA/QC processing.')
    # need to reassign attributes from cleaning stage here? -- check

    # Write station file to netcdf format
    try:
#     if True:
#         filename = station + ".nc" # Make file name
        filename = station + "_xarray.nc" # Make file name
        filepath = qaqcdir + filename # Writes file path

        tmpFile = tempfile.NamedTemporaryFile(dir = "./temp/", prefix = "_" + station, 
                                              suffix = ".nc", delete = False)
        
        t0 = time.time()
        if verbose:
            print("Saving local file {}".format(tmpFile.name))
        # Write locally
        ds.to_netcdf(tmpFile.name) # Save station file.
        if verbose:
            print("Done saving local file. Ellapsed time: {:.2f} s.".format(time.time()-t0), flush=True)
                  
        # Push file to AWS with correct file name
        t0 = time.time()
        if verbose:
            print("Pushing {} to AWS s3 bucket".format(filepath))
        s3.Bucket(bucket_name).upload_file(tmpFile.name, filepath)
        os.system("mv {} temp/{}_xarray.nc".format(tmpFile.name, station))
        
        if verbose:
            print("Done pushing file to AWS. Ellapsed time: {:.2f} s.".format(time.time()-t0), flush=True)
                  
        print('Saving {0} with dims {1} to {2}'.format(filename, ds.dims, bucket_name+"/"+qaqcdir))

#     else:
    except Exception as e:
        if verbose:
            print("netCDF writing failed for {} with Error: {}".format(filename, e))
        errors = print_qaqc_failed(errors, filename, end_api, 
                                   message='Error saving ds as .nc file to AWS bucket: {}'.format(e), 
                                   test="process_output_ds",
                                   verbose=verbose
                                  )
        return

#----------------------------------------------------------------------------
## Run full QA/QC pipeline
def run_qaqc_pipeline(ds ,network, file_name, 
                      errors, station, end_api, 
                      verbose=True):
    """
    """


    ## Add qc_flag variable for all variables, including elevation; 
    ## defaulting to nan for fill value that will be replaced with qc flag
    exclude_qaqc = ["time", "station", "lat", "lon", 
                    "qaqc_process", "sfcWind_method"] # lat and lon have a different qc check

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
            # adds new variable in shape of original variable with designated nan fill value
            ds = ds.assign({qc_var: xr.ones_like(ds[var])*np.nan})

    ##########################################################
    ## QAQC Functions

    #=========================================================
    ## PART 1

    #---------------------------------------------------------
    ## Lat-lon -- does not proceed through qaqc if failure
    new_ds = qaqc_missing_latlon(ds, verbose=verbose)
    if new_ds is None:
        errors = print_qaqc_failed(errors, station, end_api, 
                                   message="has a missing lat-lon", 
                                   test="qaqc_missing_latlon",
                                   verbose=verbose
                                  )
    else:
        stn_to_qaqc = new_ds
        if verbose:
            print('pass qaqc_missing_latlon') #testing

    #---------------------------------------------------------
    ## Within WECC -- does not proceed through qaqc if failure
    new_ds = qaqc_within_wecc(stn_to_qaqc, verbose=verbose)
    if new_ds is None:
        errors = print_qaqc_failed(errors, station, end_api, 
                                   message="lat-lon is out of range for WECC", 
                                   test="qaqc_within_wecc",
                                   verbose=verbose
                                  )
    else:
        stn_to_qaqc = new_ds
        if verbose:
            print('pass qaqc_within_wecc')

    #---------------------------------------------------------
    ## Elevation -- if DEM in-filling fails, does not proceed through qaqc
    new_ds = qaqc_elev_infill(stn_to_qaqc, verbose=verbose) # nan infilling must be before range check
    if new_ds is None:
        errors = print_qaqc_failed(errors, station, end_api, 
                                   message="DEM in-filling failed for", 
                                   test="DEM in-filling, may not mean station does not pass qa/qc -- check",
                                   verbose=verbose
                                  )
    else:
        stn_to_qaqc = new_ds

    #---------------------------------------------------------
    new_ds = qaqc_elev_range(stn_to_qaqc, verbose=verbose)
    if new_ds is None:
        errors = print_qaqc_failed(errors, station, end_api, 
                                   message="elevation out of range for WECC", 
                                   test="qaqc_elev_range",
                                   verbose=verbose
                                  )
    else:
        stn_to_qaqc = new_ds
        if verbose:
            print('pass qaqc_elev_range')

    #=========================================================
    ## Part 2: Variable logic checks

    #---------------------------------------------------------
    # precipitation is not negative
    new_ds = qaqc_precip_logic_nonegvals(stn_to_qaqc, verbose=verbose)
    if new_ds is None:
        errors = print_qaqc_failed(errors, station, end_api, 
                                   message="Flagging problem with negative precipitation values for", 
                                   test="qaqc_precip_logic_nonegvals",
                                   verbose=verbose
                                  )
    else:
        stn_to_qaqc = new_ds
        if verbose:
            print('pass qaqc_precip_logic_nonegvals')

    #---------------------------------------------------------
    ## precipitation duration logic
    new_ds = qaqc_precip_logic_accum_amounts(stn_to_qaqc, verbose=verbose)
    if new_ds is None:
        errors = print_qaqc_failed(errors, station, end_api, 
                                   message="Flagging problem with precip duration logic check for", 
                                   test="qaqc_precip_logic_accum_amounts",
                                   verbose=verbose
                                  )
    else:
        stn_to_qaqc = new_ds
        if verbose:
            print('pass qaqc_precip_logic_accum_amounts') # testing

    #======================================================================
    ## Part 3 functions (individual variable/timestamp)
    ## NDBC and MARITIME only

    #----------------------------------------------------------------
    ## Buoys with known issues with specific qaqc flags
    if network == 'MARITIME' or network == 'NDBC':
        # remove elevation_qc var from remainder of analyses so it does not also get flagged -- 
        # confirm with final qaqc process
        era_qc_vars.remove("elevation_eraqc") 

        new_ds = spurious_buoy_check(station, stn_to_qaqc, era_qc_vars, verbose=verbose)
        if new_ds is None:
            errors = print_qaqc_failed(errors, station, end_api, 
                                       message="Flagging problematic buoy issue for", 
                                       test="spurious_buoy_check",
                                   verbose=verbose
                                  )
        else:
            stn_to_qaqc = new_ds
            if verbose:
                print('pass spurious_buoy_check') # testing

    #======================================================================
    ## Part 4?

    #-----------------------------------------------------------------
    ## Sensor height: wind
    new_ds = qaqc_sensor_height_w(stn_to_qaqc, verbose=verbose)
    if new_ds is None:
        errors = print_qaqc_failed(errors, station, end_api, 
                                   message="Flagging problem with anemometer sensor height for", 
                                   test="qaqc_sensor_height_w",
                                   verbose=verbose
                                  )
    else:
        stn_to_qaqc = new_ds
        if verbose:
            print('pass qaqc_sensor_height_w')

    #-----------------------------------------------------------------
    ## Sensor height: air temperature
    new_ds = qaqc_sensor_height_t(stn_to_qaqc, verbose=verbose)
    if new_ds is None:
        errors = print_qaqc_failed(errors, station, end_api, 
                                   message="Flagging problem with thermometer sensor height for", 
                                   test="qaqc_sensor_height_t",
                                   verbose=verbose
                                  )
    else:
        stn_to_qaqc = new_ds
        if verbose:
            print('pass qaqc_sensor_height_t')

    #-----------------------------------------------------------------
    ## World record checks: air temperature, dewpoint, wind, pressure
    new_ds = qaqc_world_record(stn_to_qaqc, verbose=verbose)
    if new_ds is None:
        errors = print_qaqc_failed(errors, station, end_api, 
                                   message="Flagging problem with world record check for", 
                                   test="qaqc_world_record",
                                   verbose=verbose
                                  )
    else:
        stn_to_qaqc = new_ds
        if verbose:
            print('pass qaqc_world_record')

    #-----------------------------------------------------------------
    ## Variable cross-logic checks
    # dew point temp cannot exceed air temperature
    new_ds = qaqc_crossvar_logic_tdps_to_tas(stn_to_qaqc, verbose=verbose)
    if new_ds is None:
        errors = print_qaqc_failed(errors, station, end_api, 
                                   message="Flagging problem with temperature cross-variable logic check for", 
                                   test="qaqc_crossvar_logic_tdps_to_tas",
                                   verbose=verbose
                                  )
    else:
        stn_to_qaqc = new_ds
        if verbose:
            print('pass qaqc_crossvar_logic_tdps_to_tas')        

    #-----------------------------------------------------------------
    # wind direction should be 0 when wind speed is also 0
    new_ds = qaqc_crossvar_logic_calm_wind_dir(stn_to_qaqc, verbose=verbose)
    if new_ds is None:
        errors = print_qaqc_failed(errors, station, end_api, 
                                   message="Flagging problem with wind cross-variable logic check for", 
                                   test="qaqc_crossvar_logic_calm_wind_dir",
                                   verbose=verbose
                                  )
    else:
        stn_to_qaqc = new_ds
        if verbose:
            print('pass qaqc_crossvar_logic_calm_wind_dir')        

        #-----------------------------------------------------------------
        # Distribution checks
        # unusual gaps (part 1)
#         new_ds = qaqc_dist_gaps_part1(stn_to_qaqc)
#         if new_ds is None:
#             errors = print_qaqc_failed(errors, station, end_api, 
#                                        message="Flagging problem with unusual gap distribution function for", 
#                                        test="qaqc_dist_gaps_part1",
#                                        verbose=verbose
#                                       )
#         else:
#             stn_to_qaqc = new_ds
#             if verbose:
#                 print('pass qaqc_dist_gaps_part1')
                
        #-----------------------------------------------------------------

    return stn_to_qaqc

#==============================================================================
## Function: Conducts whole station qa/qc checks (lat-lon, within WECC, elevation)
def whole_station_qaqc(network, cleandir, qaqcdir, verbose=True):
    """
    """    
    print()
    # Set up error handling.
    errors, end_api, timestamp = setup_error_handling()

    # Read in network files
    try:
        files, stations = read_network_files(network, cleandir)
    except Exception as e:
        if verbose: 
            print(e) # testing
        errors = print_qaqc_failed(errors, station="Whole network", end_api=end_api, 
                                   message="Error in whole network:", test=e)
    
    # if files not successfully read in:
    else: 
        """
        -----------------------------------
        for station in stations: # full run
        -----------------------------------
        for station in ['ASOSAWOS_72676324198']: 
        this is the smallest ASOSAWOS file
        takes ~10 seconds to run for Victoria
        -----------------------------------
        """
        # # TESTING SUBSET
        # if sample is None:
        #     stations_sample = list(stations)
        # else:
        #     if sample>len(stations):
        #         raise NameError("{} larger than station number in {} ({}).".format(sample, network, len(stations)))
#             # Random select sampled stations
        stations_sample = stations.sample(1)
            # Select stations for timing analysis
            # stations_sample = list(stations.iloc[:sample])
        
        # Loop over stations
        for station in stations_sample:
            
            file_name = cleandir+station+".nc"
            
            if file_name not in files: # dont run qa/qc on a station that isn't cleaned
                if verbose:
                    print("{} was not cleaned - skipping qa/qc".format(station)) # testing
                message = "No cleaned data for this station, does not proceed to qa/qc: see cleaned station list for reason"
                errors = print_qaqc_failed(errors, station="Whole network", end_api=end_api, 
                                           message=message, test="whole_station_qaqc")
                continue
            else:
                if verbose:
                    print('Running QA/QC on: ', station) # testing
                
                fs = s3fs.S3FileSystem()
                aws_url = "s3://wecc-historical-wx/"+file_name

                with fs.open(aws_url) as fileObj:
#                     if True:
                    try:
                        if verbose:
                            print("Reading {}".format(aws_url), flush=True)
                        t0 = time.time()
                        # Drop station dimension
                        ds = xr.open_dataset(fileObj).isel(station=0)
                        # Load data to speed up analysis
                        ds = ds.load()
                        # Drop time duplicates
                        ds = ds.drop_duplicates(dim="time")
                        if verbose:
                            print("Done reading. Ellapsed time: {:.2f} s.".
                                  format(time.time()-t0), flush=True)
                        # CHECK THE ENGINE HERE:
                        # setting to default which operates on best with dependencies, previously 'h5netcdf'
                        
                        # Run full QA/QC pipeline
                        t0 = time.time()
                        if verbose:
                            print("Running QA/QC pipeline on {}".format(aws_url), flush=True)
                        ds = run_qaqc_pipeline(ds, network, file_name, errors, 
                                               station, end_api, verbose=verbose)
                        if verbose:
                            print("Done running QA/QC pipeline. Ellapsed time: {:.2f} s.".
                                  format(time.time()-t0), flush=True) 
                        ## Assign ds attributes and save .nc file
                        if ds is not None:
                            # Re-load station coordinate
                            ds = ds.assign_coords({"station":station}).expand_dims({"station":1})

                            # Process output
                            t0 = time.time()
                            process_output_ds(ds, network, timestamp, station, qaqcdir, 
                                              errors, end_api, verbose=verbose)
                            ds.close()
                        del(ds)
                        
#                     else:
                    except Exception as e:
                        if verbose:
                            print("run_qaqc_pipeline failed with error: {}".format(e))
                            errors = print_qaqc_failed(
                                         errors, station, end_api, 
                                         message="Cannot read files in from AWS: {0}".format(e), 
                                         test="run_qaqc_pipeline",
                                         verbose=verbose)
            
    # Write errors to csv
    finally:
        pass
        errors = pd.DataFrame(errors)
        csv_buffer = StringIO()
        errors.to_csv(csv_buffer)
        content = csv_buffer.getvalue()
        # Make sure error files save to correct directory
        s3_cl.put_object(Bucket=bucket_name, Body=content, Key=qaqcdir+"errors_{}_{}_xarray.csv".format(network, end_api)) 
        if verbose:
            print('errors_{0}_{1}.csv saved to {2}'.format(network, end_api, bucket_name + "/" + qaqcdir))
    return

