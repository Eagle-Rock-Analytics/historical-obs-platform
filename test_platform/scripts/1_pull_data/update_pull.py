'''
This script contains the functions needed to update all data on a weekly basis.
Data will be downloaded from the last date updated in AWS until 45 days prior to the present day (the window of final data).

'''
from datetime import datetime, timezone, timedelta
import ASOSAWOS_pullftp
from SCANSNOTEL_pull import get_scan_station_data
from OtherISD_pull import get_wecc_stations, get_otherisd_data_ftp
from MADIS_pull import get_madis_metadata, get_madis_station_csv
from CW3E_pull import get_cw3e_metadata, get_cw3e_update
from CIMIS_pull import get_cimis_update_ftp
from HADS_pull import get_hads_update
from MARITIME_pull import get_maritime_update, get_maritime_station_ids
from MADIS_pull import madis_update
from pull_qa import retry_downloads
import boto3
import config
import pandas as pd
from io import BytesIO

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
def get_last_date(bucket_name, folder, file_ext = None):
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
    last_time_modified_list = [x[1] for x in files] # Get most recent time
    last_time_modified = max(last_time_modified_list) # most recent date of data pull from stations that require an updated data
    return last_time_modified.date()
        

# Update script: SCAN
def update_SCAN(last_time_mod = None):
    today = datetime.now(timezone.utc).date() # Get today's day, month and year
    network = "SCAN"
    if last_time_mod is None:
        last_time_mod = get_last_date(bucket_name, f'1_raw_wx/{network}/')
    if last_time_mod < download_date:
        print(f"Downloading {network} data from {last_time_mod} to {download_date}.")
        get_scan_station_data(wecc_terr, wecc_mar, bucket_name, start_date = str(last_time_mod), end_date = str(download_date), networks = [network], fileext = str(today)) 
    else:
        print(f"{network} station files up to date.")

# Update script: SNOTEL
def update_SNOTEL(last_time_mod = None):
    today = datetime.now(timezone.utc).date() # Get today's day, month and year
    network = "SNOTEL"
    usda_network = "SNTL"
    if last_time_mod is None:
        last_time_mod = get_last_date(bucket_name, f'1_raw_wx/{network}/')
    if last_time_mod < download_date:
        print(f"Downloading {network} data from {last_time_mod} to {download_date}.")
        get_scan_station_data(wecc_terr, wecc_mar, bucket_name, start_date = str(last_time_mod), end_date = str(download_date), networks = [usda_network], fileext = str(today))      
    else:
        print(f"{network} station files up to date.")

# Update script: OtherISD
# As currently written, this will overwrite all files for the current year.
# Timeout may occasionaly occur on API end, but resolves upon another attempt after waiting a short amount of time (~hours)
def update_otherisd(last_time_mod = None):
    network = "OtherISD"
    directory = f'1_raw_wx/{network}/'
    if last_time_mod is None:
        last_time_mod = get_last_date(bucket_name, folder = directory, file_ext = '.gz')
    if last_time_mod < download_date:
        print(f"Downloading {network} data from {last_time_mod} to {download_date}.")
        stations = get_wecc_stations(wecc_terr, wecc_mar, bucket_name, directory)
        get_otherisd_data_ftp(stations, bucket_name, directory, start_date = str(last_time_mod), get_all = True)
        retry_downloads(token = config.token, bucket_name= bucket_name, networks = [network])
    else:
        print(f"{network} station files up to date.")

# Update script: ASOS-AWOS
# As currently written, this will overwrite all files for the current year.
def update_asosawos(last_time_mod = None):
    network = "ASOSAWOS"
    directory = f'1_raw_wx/{network}/'
    if last_time_mod is None:
        last_time_mod = get_last_date(bucket_name, folder = directory, file_ext = '.gz')
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
def update_cimis(last_time_mod = None):
    network = "CIMIS"
    directory = f'1_raw_wx/{network}/'
    if last_time_mod is None:
        last_time_mod = get_last_date(bucket_name, folder = directory, file_ext = '.zip')
    if last_time_mod < download_date:
        print(f"Downloading {network} data from {last_time_mod} to {download_date}.")
        get_cimis_update_ftp(bucket_name, directory, start_date = str(last_time_mod), end_date = str(download_date))
    else:
        print(f"{network} station files up to date.")

# Update script: HADS
def update_hads(last_time_mod = None):
    network = "HADS"
    directory = f'1_raw_wx/{network}/'
    if last_time_mod is None:
        last_time_mod = get_last_date(bucket_name, folder = directory, file_ext = '.gz')
    
    if last_time_mod < download_date:
        print(f"Downloading {network} data from {last_time_mod} to {download_date}.")
        get_hads_update(bucket_name, directory, start_date = str(last_time_mod), end_date = str(download_date))

    else:
        print(f"{network} station files up to date.")

# Update script: CW3E
# Script will download multiple byte files for each station for each day selected.
# The LBH station file will always be redownloaded completely. 
# At this point, this station gets dropped during the cleaning phase, so this should not affect anything.
# Massive number of files will make this update script slower than the others.
def update_cw3e(last_time_mod = None):
    network = "CW3E"
    directory = f'1_raw_wx/{network}/'
    if last_time_mod is None:
        last_time_mod = get_last_date(bucket_name, folder = directory, file_ext = 'm')
    
    if last_time_mod < download_date:
        print(f"Downloading {network} data from {last_time_mod} to {download_date}.")
        get_cw3e_metadata(token = config.token, terrpath = wecc_terr, marpath = wecc_mar, bucket_name = bucket_name, directory = directory)
        get_cw3e_update(bucket_name, directory, start_date = str(last_time_mod), end_date = str(download_date))
    else:
        print(f"{network} station files up to date.")

# Update script: MARITIME
# Update delay by (minimum of) 45 days. This means 2022 data isn't available in yearly format until ~ Feb 15 2022, e.g.
# Note: NDBC archive has a different file extenstion for preliminary data that is not updated until the full month has passed the preliminary period.
# Update/automation will therefore need to restrict to only pulling complete months, which may be beyond 45 days
def update_maritime(last_time_mod = None):
    network = "MARITIME"
    directory = f'1_raw_wx/{network}/'
    if last_time_mod is None:
        last_time_mod = get_last_date(bucket_name, folder = directory)
        # Note: File extension is either .gz or .zip, so do not filter by file extension here as errors & station files will be dropped automatically.

    if last_time_mod < download_date:
        print(f"Downloading {network} data from {last_time_mod} to {download_date}.")
        stations = get_maritime_station_ids(wecc_terr, wecc_mar, '1_raw_wx/MARITIME/','1_raw_wx/NDBC/')
        get_maritime_update(stations, bucket_name, network, start_date = str(last_time_mod), end_date = str(download_date))
    else:
        print(f"{network} station files up to date.")

# Update script: MARITIME
# Update delay by (minimum of) 45 days. This means 2022 data isn't available in yearly format until ~ Feb 15 2022, e.g.
# Note: NDBC archive has a different file extenstion for preliminary data that is not updated until the full month has passed the preliminary period.
# Update/automation will therefore need to restrict to only pulling complete months, which may be beyond 45 days
def update_ndbc(last_time_mod = None):
    network = "NDBC"
    directory = f'1_raw_wx/{network}/'
    if last_time_mod is None:
        last_time_mod = get_last_date(bucket_name, folder = directory) 
        # Note: File extension is either .gz or .zip, so do not filter by file extension here as errors & station files will be dropped automatically.
    
    if last_time_mod < download_date:
        print(f"Downloading {network} data from {last_time_mod} to {download_date}.")
        stations = get_maritime_station_ids(wecc_terr, wecc_mar, '1_raw_wx/MARITIME/','1_raw_wx/NDBC/')
        get_maritime_update(stations, bucket_name, network, start_date = str(last_time_mod), end_date = str(download_date))
    else:
        print(f"{network} station files up to date.")

# Update script: MADIS
def update_madis(network, last_time_mod = None):
    directory = f'1_raw_wx/{network}/'
    if network == "CAHYDRO":
        if last_time_mod is None:
            last_time_mod = get_last_date(bucket_name, folder = directory, file_ext = '.csv')
    
        if last_time_mod < download_date:
            print(f"Downloading {network} data from {last_time_mod} to {download_date}.")
            madis_update(token= config.token, networks = ["CA HYDRO"], pause = None, start_date = str(last_time_mod), end_date = str(download_date))
        else:
            print(f"{network} station files up to date.")

    else: # all other MADIS networks
        if last_time_mod is None:
            last_time_mod = get_last_date(bucket_name, folder = directory, n = int(n[network]), file_ext = '.csv')
        
        if last_time_mod < download_date:
            print(f"Downloading {network} data from {last_time_mod} to {download_date}.")
            madis_update(token= config.token, networks = [network], pause = None, start_date = str(last_time_mod), end_date = str(download_date))
        else:
            print(f"{network} station files up to date.")


if __name__ == "__main__":
    last_time_mod = None # option to custom select start date of download (format: datetime(year, month, day).date())
    # update_SCAN(last_time_mod)
    # update_SNOTEL(last_time_mod)
    # update_otherisd(last_time_mod)
    # update_asosawos(last_time_mod)
    update_cimis(last_time_mod)
    # update_hads(last_time_mod)
    # update_cw3e(last_time_mod)
    # update_maritime(last_time_mod)
    # update_ndbc(last_time_mod)
    # update_madis(network = 'CAHYDRO', last_time_mod)
