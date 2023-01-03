'''
This script contains the functions needed to update all data on a weekly basis.
Data will be downloaded from the last date updated in AWS until 45 days prior to the present day (the window of final data).

'''
from datetime import datetime, timezone, timedelta
import ASOSAWOS_pullftp
from SCANSNOTEL_pull import get_scan_station_data
from OtherISD_pull import get_wecc_stations, get_otherisd_data_ftp
from MADIS_pull import get_madis_metadata, get_madis_station_csv
from CIMIS_pull import get_cimis_update_ftp
from HADS_pull import get_hads_update
from pull_qa import retry_downloads
import boto3
import config

# Set AWS credentials
s3 = boto3.resource('s3')
s3_cl = boto3.client('s3')
bucket_name = 'wecc-historical-wx'

# Set file paths
wecc_terr = "s3://wecc-historical-wx/0_maps/WECC_Informational_MarineCoastal_Boundary_land.shp"
wecc_mar = "s3://wecc-historical-wx/0_maps/WECC_Informational_MarineCoastal_Boundary_marine.shp"


# Set parameters
today = datetime.now(timezone.utc).date() # Get today's day, month and year
download_window = 45 # Preliminary data lag
download_date = (today - timedelta(days = download_window)) # Set to midnight UTC time

# Get last download date of files in AWS folder
def get_last_date(bucket_name, folder, n, file_ext = None):
    # Inputs: 
    # bucket_name: name of AWS bucket
    # folder: file path to relevant folder
    # n: the number of files to compare the most recent download date to (determined as a percentage of no. of stations per network)
    # file_ext: optional, file extension to look for (e.g. .gz, only useful if not .csv or multiple file types)
    files = s3.Bucket(bucket_name).objects.filter(Prefix = folder)
    if file_ext:
        files = [[obj.key, obj.last_modified] for obj in sorted(files, key=lambda x: x.last_modified, 
            reverse=True) if obj.key.endswith(file_ext) and not any(x in obj.key for x in ('errors', 'station'))]
    else:
        files = [[obj.key, obj.last_modified] for obj in sorted(files, key=lambda x: x.last_modified, 
            reverse=True) if not any(x in obj.key for x in ('errors', 'station'))]
    
    # Get most recent and nth most recent time
    # If they differ by more than one day, take the nth most recent time (as an indicator that the most recent file may not represent the rest of the data)
    if (files[0][1] - files[n][1])>timedelta(days=1):
        last_time_modified = files[n][1]
    else:
        last_time_modified = files[0][1]
    
    return last_time_modified.date()

# Update script: SCAN
def update_SCAN():
    network = "SCAN"
    last_time_mod = get_last_date(bucket_name, f'1_raw_wx/{network}/', n = 25)
    if last_time_mod < download_date:
        print(f"Downloading {network} data from {last_time_mod} to {download_date}.")
    #get_scan_station_data(wecc_terr, wecc_mar, bucket_name, start_date = str(last_time_mod), end_date = str(download_date), networks = [network]) 
    # # TO DO: how do we want to deal with file names for update snippets??? give separate name? save in separate folder????

# Update script: SNOTEL
def update_SNOTEL():
    network = "SNOTEL"
    usda_network = "SNTL"
    last_time_mod = get_last_date(bucket_name, f'1_raw_wx/{network}/', n = 25)
    if last_time_mod < download_date:
        print(f"Downloading {network} data from {last_time_mod} to {download_date}.")
    else:
        print(f"{network} station files up to date.")
    #get_scan_station_data(wecc_terr, wecc_mar, bucket_name, start_date = str(last_time_mod), end_date = str(download_date), networks = [usda_network]) 
    # # TO DO: how do we want to deal with file names for update snippets??? give separate name? save in separate folder????

# Update script: OtherISD
# As currently written, this will overwrite all files for the current year.
def update_otherisd():
    network = "OtherISD"
    directory = f'1_raw_wx/{network}/'
    last_time_mod = get_last_date(bucket_name, folder = directory, n = 25, file_ext = '.gz')
    if last_time_mod < download_date:
        print(f"Downloading {network} data from {last_time_mod} to {download_date}.")
        stations = get_wecc_stations(wecc_terr, wecc_mar, bucket_name, directory)
        get_otherisd_data_ftp(stations, bucket_name, directory, start_date = str(last_time_mod), get_all = True)
        retry_downloads(token = config.token, bucket_name= bucket_name, networks = [network])
    else:
        print(f"{network} station files up to date.")

# Update script: ASOS-AWOS
# As currently written, this will overwrite all files for the current year.
def update_asosawos():
    network = "ASOSAWOS"
    directory = f'1_raw_wx/{network}/'
    last_time_mod = get_last_date(bucket_name, folder = directory, n = 25, file_ext = '.gz')
    if last_time_mod < download_date:
        print(f"Downloading {network} data from {last_time_mod} to {download_date}.")
        stations = ASOSAWOS_pullftp.get_wecc_stations(wecc_terr, wecc_mar)
        ASOSAWOS_pullftp.get_asosawos_data_ftp(stations, bucket_name, directory, start_date = str(last_time_mod))
        retry_downloads(token = config.token, bucket_name= bucket_name, networks = [network])
    else:
        print(f"{network} station files up to date.")

# Update script: CIMIS
# No retry download method available.
# This may overwrite the most recent if pull is repeated more frequently than monthly.
def update_cimis():
    network = "CIMIS"
    directory = f'1_raw_wx/{network}/'
    last_time_mod = get_last_date(bucket_name, folder = directory, n = 25, file_ext = '.zip')
    if last_time_mod < download_date:
        print(f"Downloading {network} data from {last_time_mod} to {download_date}.")
        get_cimis_update_ftp(bucket_name, directory, start_date = str(last_time_mod), end_date = str(download_date))
    else:
        print(f"{network} station files up to date.")

# Update script: HADS
def update_hads():
    network = "HADS"
    directory = f'1_raw_wx/{network}/'
    last_time_mod = get_last_date(bucket_name, folder = directory, n = 25, file_ext = '.gz')
    
    # if last_time_mod < download_date:
    #     print(f"Downloading {network} data from {last_time_mod} to {download_date}.")
    get_hads_update(bucket_name, directory, start_date = str(last_time_mod), end_date = str(download_date))

# To do:
# CW3E needs the download date methods to be updated to work correctly here
# HADS names file by day, need to update time filtering method to work like this rather than by year.
# MADIS: add end_time as parameter to input into url.


if __name__ == "__main__":
    #update_SCAN()
    #update_SNOTEL()
    #update_otherisd()
    #update_asosawos()
    #update_cimis()
    update_hads()
