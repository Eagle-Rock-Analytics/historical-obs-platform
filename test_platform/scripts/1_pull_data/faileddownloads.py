# For each network, compare station list to existing file names and generate list of failed downloads.
# Then, try to download again.

from MADIS_pull import get_madis_station_csv
from SCAN_pull import get_scan_station_data
import boto3
import pandas as pd
import config

bucket_name = "wecc-historical-wx"
s3 = boto3.resource("s3")  
s3_cl = boto3.client("s3") 

# 10/19: The following networks have 1+ station that failed to download
# CWOP - done
# RAWS
# SNOTEL
# LOXWFO - done
# SCAN - done
# MTRWFO

# For MADIS networks
def madis_retry_downloads(token, bucket_name, network):
    # Get list of files in folder
    prefix = "1_raw_wx/"+network
    files = []
    for item in s3.Bucket(bucket_name).objects.filter(Prefix = prefix): 
        file = str(item.key)
        files += [file]
    files = list(filter(lambda f: f.endswith(".csv"), files)) # Get list of file names

    station_file = [file for file in files if "station" in file] # ID station file
    station_file = str(station_file[0])
    files = [file for file in files if "errors" not in file]
    files = [file for file in files if "station" not in file] # Remove error and station list files

    # Get only station IDs from file names
    stations = [file.split("/")[-1] for file in files]
    stations = [file.replace(".csv", '') for file in stations]

    # Read in station list
    station_list = s3_cl.get_object(Bucket= bucket_name, Key= station_file)
    station_list = pd.read_csv(station_list['Body'])

    # Get list of IDs not in download folder
    missed_stations = [id for id in station_list['STID'] if id not in stations]
    missed_ids = station_list[['STID', 'start']] # Format list in way that MADIS_pull script wants it.
    missed_ids = missed_ids[missed_ids.STID.isin(missed_stations)]

    # Check for duplicate stations in station list
    dup_stations = station_list.loc[station_list.duplicated(subset='STID', keep=False)]
    dup_stations = dup_stations[['STID', 'start']]
    dup_todownload = dup_stations.sort_values(by='start', ascending=True).drop_duplicates( # Keep first start date.
    keep='first', subset=['STID'])

    # If any stations duplicated, take the earlier start data and add to the redownload queue.
    download_ids = pd.merge(missed_ids, dup_todownload, how='outer')
    download_ids = download_ids.sort_values(by='start', ascending=True).drop_duplicates(
    keep='first', subset=['STID']) # If download_ids has duplicate STIDs, take first record.

    # Reorganize start date to meet API specs.
    download_ids['start'] = pd.to_datetime(download_ids['start'], format='%Y-%m-%dT%H:%M:%SZ')
    download_ids['start'] = download_ids['start'].dt.strftime('%Y%m%d%H%M')

    # Print list of stations to download.
    print(download_ids)
    
    # Note here we ignore the end date of files, since we will be trimming the last 2 months of data anyways.
    # This could be changed down the road as these dates diverge.
    errors = get_madis_station_csv(token = token, bucket_name = bucket_name, directory = prefix+"/", ids = download_ids, timeout = False)
    
    # Manually print out errors for immediate verification of success. Will also save to AWS.
    print(errors)

#madis_retry_downloads(token = config.token, bucket_name = bucket_name, network = "LOXWFO")

# For SCAN/SNOTEL
    # Get list of files in folder
def scan_retry_downloads(bucket_name, network, terrpath, marpath):
    prefix = "1_raw_wx/"+network[0]
    files = []
    for item in s3.Bucket(bucket_name).objects.filter(Prefix = prefix): 
        file = str(item.key)
        files += [file]
    files = list(filter(lambda f: f.endswith(".csv"), files)) # Get list of file names

    station_file = [file for file in files if "station" in file] # ID station file
    station_file = str(station_file[0])
    files = [file for file in files if "errors" not in file]
    files = [file for file in files if "station" not in file] # Remove error and station list files

    # Get only station IDs from file names
    stations = [file.split("/")[-1] for file in files]
    stations = [file.replace(".csv", '') for file in stations]
    
    # Read in station list
    station_list = s3_cl.get_object(Bucket= bucket_name, Key= station_file)
    station_list = pd.read_csv(station_list['Body'])
    
    # # Get list of IDs not in download folder
    missed_ids = [id for id in station_list['stationTriplet'] if id not in stations]
    missed_stations = station_list.loc[station_list['stationTriplet'].isin(missed_ids)]
    
    # # Print list of stations to download.
    print(missed_stations)
    
    # # Note here we ignore the end date of files, since we will be trimming the last 2 months of data anyways.
    # # This could be changed down the road as these dates diverge.
    errors = get_scan_station_data(terrpath, marpath, bucket_name, start_date = None, stations = missed_stations, primary = True, networks = network)

    # # Manually print out errors for immediate verification of success. Will also save to AWS.
    print(errors)

# Set paths to WECC shapefiles in AWS bucket.
wecc_terr = "s3://wecc-historical-wx/0_maps/WECC_Informational_MarineCoastal_Boundary_land.shp"
wecc_mar = "s3://wecc-historical-wx/0_maps/WECC_Informational_MarineCoastal_Boundary_marine.shp"

scan_retry_downloads(bucket_name = bucket_name, network = ["SCAN"], terrpath = wecc_terr, marpath = wecc_mar)