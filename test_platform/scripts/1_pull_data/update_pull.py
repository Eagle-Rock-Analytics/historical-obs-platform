'''
This script contains the functions needed to update all data on a weekly basis.
Data will be downloaded from the last date updated in AWS until 45 days prior to the present day (the window of final data).

'''
from datetime import datetime, timezone, timedelta
import boto3

# Set AWS credentials
s3 = boto3.resource('s3')
s3_cl = boto3.client('s3')

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
    # file_ext: optional, file extension to look for (e.g. .gz)
    files = s3.Bucket(bucket_name).objects.filter(Prefix = folder)
    if file_ext:
        files = [[obj.key, obj.last_modified] for obj in sorted(files, key=lambda x: x.last_modified, 
            reverse=True) if obj.key.endswith(file_ext) and not any(x in obj.key for x in ('errors', 'station'))]
    else:
        files = [[obj.key, obj.last_modified] for obj in sorted(files, key=lambda x: x.last_modified, 
            reverse=True) if not any(x in obj.key for x in ('errors', 'station'))]
    
    # Get most recent and nh most recent time
    # If they differ by more than one day, take the nth most recent time (as an indicator that the most recent file may not represent the rest of the data)
    if (files[0][1] - files[n][1])>timedelta(days=1):
        last_time_modified = files[n][1]
    else:
        last_time_modified = files[0][1]
    
    return last_time_modified

last_time_mod = get_last_date('wecc-historical-wx', '1_raw_wx/HADS/', n = 25, file_ext = ".gz")
print(last_time_mod)

