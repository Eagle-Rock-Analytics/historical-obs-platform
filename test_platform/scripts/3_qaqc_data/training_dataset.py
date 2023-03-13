# -*- coding: utf-8 -*-
"""
Identify representative testing set of stations for Stage 3: QA/QC development

"""

# Import libraries
import boto3
import pandas as pd
import numpy as np
import matplotlib.pylab as plt
import matplotlib
from io import BytesIO, StringIO
from datetime import datetime, timezone, date
import geopandas as gpd
import os
import xarray as xr
import rioxarray

# Function: iterate through raw folder and get all station lists
def get_station_list_paths(bucket_name, directory):
    # Set up variables
    s3 = boto3.resource("s3")
    s3_cl = boto3.client('s3') # for lower-level processes
    get_last_modified = lambda obj: int(obj['LastModified'].strftime('%s')) #  Write method to get last modified file

    # Read in all station lists.
    # Get list of folder prefixes
    response = s3_cl.list_objects_v2(Bucket=bucket_name, Prefix = directory, Delimiter = '/')

    networks = {'Network':[], 'NetworkPath':[], 'StationFile':[]}

    for prefix in response['CommonPrefixes']: # For each folder path
        networkpath = prefix['Prefix'][:-1]
        networkname = networkpath.split("/")[-1]
        station_file = s3.Bucket(bucket_name).objects.filter(Prefix = networkpath+"/"+"stationlist_").all()
        if len(list(station_file))==1: # If one item returned
            for item in station_file: # Get station file
                networks['Network'].append(networkname)
                networks['NetworkPath'].append(networkpath)
                networks['StationFile'].append(item.key)
                break # If more than one file of this format found in folder, just take the most recent.
        elif len(list(station_file))==0: # If no items found using search above
            files = s3.Bucket(bucket_name).objects.filter(Prefix = networkpath+"/") # List all files in folder
            file = [file for file in files if "station" in file.key] # More general search for 'station'
            for item in file: # Keep all files found here. These files may be different (e.g. ISD ASOS/AWOS vs ASOS/AWOS station lists)
                networks['Network'].append(networkname)
                networks['NetworkPath'].append(networkpath)
                networks['StationFile'].append(item.key)
        elif len(list(station_file))>1: # If more than one identically formatted station list returned (shouldn't happen), take most recent
            file = [obj.key for obj in sorted(station_file, key=lambda x: x.last_modified, reverse=True)] # Sort station files by most recent edit
            networks['Network'].append(networkname)
            networks['NetworkPath'].append(networkpath)
            networks['StationFile'].append(file[0]) # Add path to most recently changed file. (Tested and works)
        # Note method currently doesn't have a method for dealing with more than one normally formatted station file 
    return networks        

def get_station_list(bucket_name, directory):
    # Set up variables
    s3 = boto3.resource("s3") # Note these calls are reproduced in get_station_list_paths, but this enables us to run that function separately.
    s3_cl = boto3.client('s3') # for lower-level processes
    dffull = pd.DataFrame()
    
    networks = get_station_list_paths(bucket_name, directory)
    networks = pd.DataFrame(networks)
    
    # Remove ASOS/AWOS, madis and isd lists to prevent duplicates
    # note: for station masterlist, you'll probably want to KEEP the asosawos_stations and remove isd_asosawos or merge both!
    # we do this here because the ISD station list has start/end dates, and asosawos does not.
    # If station has more than one type of station file (e.g. ASOSAWOS), manually remove the one you don't want to use here.
    remove_list = ['/stationlist_ASOSAWOS.csv', '/isd_stations.csv', 'madis_stations.csv'] # / used to keep otherisd_stations.csv
    networks = networks[~networks.StationFile.str.contains('|'.join(remove_list))]

    # Check that no network has >1 station file and break code if it does.
    boolean = networks["Network"].is_unique

    total = 0 

    #print(boolean)
    if boolean is False: # Tested by commenting out two remove lines above, works.
        dupes = networks['Network'].duplicated()
        doubles = list(networks['Network'][dupes])
        print("Error: More than one station file found for the following networks: {}. Inspect folder, add duplicated files to remove_list and re-run function.".format(doubles))   
        return
    
    for index, row in networks.iterrows():
        try:
            #print(row['StationFile'])
            # Get data
            df = pd.DataFrame()
            obj = s3_cl.get_object(Bucket = bucket_name, Key = row['StationFile']) # Get file using StationFile as path.
            body = obj['Body'].read()
            temp = pd.read_csv(BytesIO(body), encoding='utf8')  

        except Exception as e: # If there's an encoding error when reading in file
            #print("Error parsing station file for {}: {}".format(row['Network'], e))
            try:
                if 'xlsx' in row['StationFile']: # If file is .xlsx file (CIMIS)
                    temp = pd.read_excel(BytesIO(body), engine = 'openpyxl') # Use different pandas function to read in.
            except Exception as e: # If there's an encoding error when reading in file
                print("Error parsing station file for {}: {}".format(row['Network'], e))  
        
        try:
            networkname = row['Network']

            print(networkname, len(temp))
            total += len(temp)
            print(total)
            

            #Deal with distinct formatting of different station lists.
            # Make column names to lower.
            temp.columns = temp.columns.str.lower()
            
            # Delete index column
            remove = [col for col in temp.columns if "unnamed" in col]
            if remove is not None:
                temp = temp.drop(remove, axis = 1)
            
            # Each df should have 6 columns.
            # network: network name.
            # name: station name
            # latitude
            # longitude
            # start-date: date station started running
            # end-date: date station stopped running
            # pulled: Y/N/NA, if station assessed as downloaded by pull_qa.py
            # time_checked: time of pull check.
            
            # # name: station name, station_name, name, dcp location name --> use name as filter.
            if any('name' in str for str in temp.columns):
                colname = [col for col in temp.columns if "name" in col]
                if len(colname)>1: # If more than one col returned
                    removelist = set(['countyname', 'name', 'altname'])
                    colname = list(set(colname)-removelist) # Use sets to exclude partial matches (e.g. 'name' in 'countyname')
                    if len(colname)>1:
                        print("Too many options for station name columns. Add manually to removelist: {}".format(colname))
                        break
                
                #print(temp[colname])
                df['name'] = temp[colname].values.reshape(-1, ) # Assign column to df.
                #print(df)

            # latitude: lat or latitude
            if any('lat' in str for str in temp.columns):
                colname = [col for col in temp.columns if "lat" in col]
                if len(colname)>1: # If more than one col returned
                    removelist = set([])
                    colname = list(set(colname)-removelist) # Use sets to exclude partial matches (e.g. 'name' in 'countyname')
                    if len(colname)>1:
                        print("Too many options for latitude columns. Add manually to removelist: {}".format(colname))
                        break
                df['latitude'] = temp[colname].values.reshape(-1, )
            else:
                df['latitude'] = np.nan

            # longitude: lat or latitude
            if any('lon' in str for str in temp.columns):
                colname = [col for col in temp.columns if "lon" in col]
                if len(colname)>1: # If more than one col returned
                    removelist = set([])
                    colname = list(set(colname)-removelist) # Use sets to exclude partial matches (e.g. 'name' in 'countyname')
                    if len(colname)>1:
                        print("Too many options for longitude columns. Add manually to removelist: {}".format(colname))
                        break
                df['longitude'] = temp[colname].values.reshape(-1, )
            else:
                df['longitude'] = np.nan

            # elevation: elev or elevation (TO DO: convert to same unit!!)
            if any('elev' in str for str in temp.columns):
                colname = [col for col in temp.columns if "elev" in col]
                if len(colname)>1: # If more than one col returned
                    removelist = set(['anemometer_elev', 'barometer_elev', 'elev'])
                    colname = list(set(colname)-removelist) # Use sets to exclude partial matches (e.g. 'name' in 'countyname')
                    if len(colname)>1:
                        if 'elev_dem' in colname:
                            colname.remove("elev_dem")
                        if len(colname)>1:
                            print("Too many options for elevation columns. Add manually to removelist: {}".format(colname))
                            break
                df['elevation'] = temp[colname].values.reshape(-1, )
            else:
                df['elevation'] = np.nan
            
            
            # # start-date: search for start, begin or connect
            if any(y in x for x in temp.columns for y in ['begin', 'start', 'connect']):
                colname = [col for col in temp.columns if any(sub in col for sub in ['begin', 'start', 'connect'])]
                if len(colname)>1: # If more than one col returned
                    removelist = set(['start_time', 'begin'])# Add any items to be manually removed here.
                    colname = list(set(colname)-removelist) # Use sets to exclude partial matches (e.g. 'name' in 'countyname')
                    if len(colname)>1:
                        if 'start_time' in colname: # If both start_time (parsed) and begin (not parsed) columns present, remove begin.
                            if 'begin' in colname:
                                colname.remove("begin")
                        if 'disconnect' in colname:
                            colname.remove("disconnect")
                        if len(colname)>1:
                            print("Too many options for start date columns. Add manually to removelist: {}".format(colname))
                            break
                df['start-date'] = temp[colname].values.reshape(-1, )
            else: # If no start date provided
                df['start-date'] = np.nan
            
            # # end-date: search for end or disconnect
            if any(y in x for x in temp.columns for y in ['end', 'disconnect']):
                colname = [col for col in temp.columns if any(sub in col for sub in ['end', 'disconnect'])]
                if len(colname)>1: # If more than one col returned
                    removelist = set([])# Add any items to be manually removed here.
                    colname = list(set(colname)-removelist) # Use sets to exclude partial matches (e.g. 'name' in 'countyname')
                    if len(colname)>1:
                        if 'end_time' in colname: # If both start_time (parsed) and begin (not parsed) columns present, remove begin.
                            if 'end' in colname:
                                colname.remove("end")
                        if len(colname)>1:
                            print("Too many options for end date columns. Add manually to removelist: {}".format(colname))
                            break
                df['end-date'] = temp[colname].values.reshape(-1, )
            else: # If no start date provided
                df['end-date'] = np.nan
                
            # Add cleaned and date checked columns, if they exist
            if any('cleaned' in str for str in temp.columns):
                df['cleaned'] = temp['cleaned'].values.reshape(-1, )
            else:
                df['cleaned'] = np.nan

            if any('time_checked' in str for str in temp.columns):
                df['time_checked'] = temp['time_checked'].values.reshape(-1, )
            else:
                df['time_checked'] = np.nan

            # Make network column
            df['network'] = networkname

            dffull = pd.concat([dffull, df], sort = False)

        except Exception as e:
            print(e)

    # Organize full dataframe.
    # If end date is "active", make this be today's date.
    today = datetime.now()
    
    print(len(dffull))
    dffull['end-date'] = dffull['end-date'].replace('Active',today)

    # Format dates in datetime format.
    dffull['start-date'] = pd.to_datetime(dffull['start-date'], utc = True)
    dffull['end-date'] = pd.to_datetime(dffull['end-date'], utc = True)

    # Drop empty rows - lat/lon as min criteria for inclusion
    dffull.dropna(subset=["latitude"], inplace=True)
    print(len(dffull))
    dffull.dropna(subset=["longitude"], inplace=True)
    print(len(dffull))
    # Remove any duplicates (of network and ID)
    dffull.sort_values(by=['start-date'], ascending = True, inplace = True)# Sort by network and time so oldest network is always first
    dffull.drop_duplicates(subset = ['name', 'latitude', 'longitude', 'network'], inplace = True)
    print(len(dffull))
    # Resort by network
    dffull.sort_values(by=['network'], inplace = True)

    # Reset index.
    dffull = dffull.reset_index(drop = True)

    return dffull

# Create station list data frame from cleaned data
# Run functions - generate station chart
bucket_name = "wecc-historical-wx"
directory = "2_clean_wx/"

dffull = get_station_list(bucket_name, directory)

# Format dates in datetime format (this gets lost in import).
dffull['start-date'] = pd.to_datetime(dffull['start-date'], utc = True)
dffull['end-date'] = pd.to_datetime(dffull['end-date'], utc = True)

# Quality control.
# Fix nas
## Filter out rows w/o start date
subdf = dffull.loc[~dffull['start-date'].isnull()].copy()
# Filter out rows without data between 1980 and now.
subdf = subdf.loc[(subdf['start-date']<=datetime.utcnow().replace(tzinfo=timezone.utc)) & (subdf['end-date']>='1980-01-01')]
# Filter out stations that were not cleaned
subdf = subdf[subdf['cleaned'] == 'Y']
    
# Make a geodataframe.
gdf = gpd.GeoDataFrame(subdf, geometry=gpd.points_from_xy(subdf.longitude, subdf.latitude))
gdf.set_crs(epsg=4326, inplace=True) # Set CRS
    
# Read in BIOCLIM vars as one xarray object
bio = xr.open_mfdataset('map_data/wc2.1_30s_bio_*.tif', combine = 'nested',
                        concat_dim = 'var', engine = 'rasterio')

# Extract BIOCLIM vars at station sites
appended_data = []
for index, row in subdf.iterrows():
    stn_bio = bio.sel(x = row['longitude'], y = row['latitude'], method = 'nearest')
    bio_df = stn_bio.to_dataframe() # why is this so slow??
    bio_df['name'] = row['name'] 
    appended_data.append(bio_df)
    print(row['name'] + ' complete')

subdf.to_csv('cleaned_stations.csv')

# Elevation distribution
num_bins = 20
# the histogram of the data
n, bins, patches = plt.hist(subdf['elevation'], num_bins, facecolor='blue', alpha=0.5)
plt.xlabel('Elevation (m)')
plt.ylabel('Stations (count)')
plt.show()


