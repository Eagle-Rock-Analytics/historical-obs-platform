# %%
import boto3
import pandas as pd
import numpy as np
import matplotlib.pylab as plt
from io import BytesIO

#%%
# Set up variables
bucket_name = "wecc-historical-wx"
s3 = boto3.resource("s3")
s3_cl = boto3.client('s3') # for lower-level processes


# Read in all station lists.
files = []
s3_bucket = s3.Bucket(bucket_name)
station_files = [f.key for f in s3_bucket.objects.all() if 'station' in f.key]

# Remove ISD ASOS/AWOS, madis and isd lists to prevent duplicates
remove_list = ['isd_asosawos_stations.csv', '/isd_stations.csv', 'madis_stations.csv'] # / used to keep otherisd_stations.csv
station_files = [r for r in station_files if not any(z in r for z in remove_list)]

print(station_files)
dffull = pd.DataFrame()

# For each network, coerce data to common format.
for network in station_files:
    # Get data
    df = pd.DataFrame()
    obj = s3.Object(bucket_name, network)
    body = obj.get()['Body'].read()
    temp = pd.read_csv(BytesIO(body), encoding='utf8')  
    networkname = network.split("/")[-2] # Get network name
    

    # Deal with distinct formatting of different station lists.
    # Make column names to lower.
    temp.columns = temp.columns.str.lower()
    
    # Start date: search for start, begin or connect
    # End date: search for end or disconnect
    # Do this by network and add any unique network station lists explicitly here.
    if 'asosawos' in network:
        df['name'] = temp['name']
        df['start-date'] = temp['startdate']
        df['end-date'] = np.nan

    elif 'otherisd' in network:
        df['name'] = temp['station name']
        df['start-date'] = temp['start_time']
        df['end-date'] = temp['end_time']

    elif 'cimis' in network:
        print(temp) # This network is causing issues, troubleshoot.
        df['name'] = temp['name']
        df['start-date'] = temp['connect']
        df['end-date'] = temp['disconnect']

    elif 'hads' in network:
        df['name'] = temp['dcp location name']
        df['start-date'] = np.nan
        df['end-date'] = np.nan

    elif 'maritime' in network:
        pass # Waiting on maritime station list

    elif 'scan' in network:
        df['name'] = temp['name']
        df['start-date'] = temp['begindate']
        df['end-date'] = temp['enddate']

    else: # MADIS networks
        df['name'] = temp['name']
        df['start-date'] = temp['start']
        df['end-date'] = temp['end']

    # Make network column
    df['network'] = networkname
    print(networkname)
    print(df)
    dffull = pd.concat([dffull, df], sort = False)

print(dffull)
    
# Join all together.

# Format start date of each station - round down to month.
dffull['start-date'] = pd.to_datetime(dffull['start-date'])
dffull['end-date'] = pd.to_datetime(dffull['start-date'])
dffull['period'] = [pd.period_range(*v, freq='M') for v in zip(dffull['start-date'], dffull['end-date'])]

print(dffull)

#out = df.explode('period').pivot_table(
#        'Value', 'network', 'period', aggfunc='count', fill_value=0)
#out.columns = out.columns.strftime('%b-%y')

# Group by starting month and network and sum
#sumdf = df.groupby(['start-date','network']).size().reset_index(name="Time")

#### Visualizing station count
# For each month, calculate # of stations with data by network.

#%%
# https://pandas.pydata.org/docs/reference/api/pandas.DataFrame.plot.area.html
# https://www.geeksforgeeks.org/matplotlib-pyplot-stackplot-in-python/


# Eventual extension: make same plot but stations by variable.

# NOTES: To make MASTER STATION LIST (below)
    # Save one unified master station list
#     # Station ID: network ID + unique station ID
    
    
#     # Station name: search for name, station
#     if 'name' in temp.columns:
#         print(temp['name'])

#     break

#     # Latitude: search for latitude or lat

#     # Longitude: search for longitude or lon

#     # Elevation: search for elevation or elev

# # # By network: ASOS/AWOS
# #     if 'asosawos' in network:
# #         df['era-id'] = "ASOSAWOS_"+ # Unique ID - TBD with this network!!!!
# #         df['station-name'] = temp['name']
# #         df['station-alt-name'] = temp['altname']
# #         df['latitude'] = temp['lat']
# #         df['longitude'] = temp['lon']
# #         df['elevation'] = temp['elev']
# #         df['utc'] = temp['utc']
# #         df['start-date'] = temp['startdate']
# #         df['end-date'] = np.nan

#     if 'hads' in network:
#         df['era-id'] = "HADS_"+ # Unique ID - TBD with this network!!!!
#         df['station-name'] = temp['dcp location name']
#         df['station-alt-name'] = np.nan
#         df['latitude'] = temp['latitude']
#         df['longitude'] = temp['longitude']
#         df['elevation'] = np.nan
#         df['utc'] = np.nan
#         df['start-date'] = np.nan
#         df['end-date'] = np.nan
