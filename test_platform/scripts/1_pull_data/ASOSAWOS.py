### Draft 1: Cleaning script for ASOS/AWOS network.

# ## Load packages
# import requests

# ## Read in ASOS/AWOS data from ISD using API GET call.
### To do this, get all station IDs in WECC, download data for all stations, then filter by source.


## Load packages
from ctypes import cast
from ftplib import FTP
from pydap.client import open_url
import wget
import os
from datetime import datetime, timezone
import xarray as xr
import pandas as pd
import fiona
from shapely.geometry import shape, Point
import numpy as np
import pandas as pd
import geopandas as gp
from geopandas.tools import sjoin
import requests


# Set envr variables
workdir = "/home/ella/Desktop/Eagle-Rock/Historical-Data-Platform/ASOS/"

def get_wecc_stations(): #Could alter script to have shapefile as input also, if there's a use for this.
## Login.

    ## using ftplib, get list of stations as csv
    filename = 'isd-history.csv'
    ftp = FTP('ftp.ncdc.noaa.gov')
    ftp.login() # user anonymous, password anonymous
    ftp.cwd('pub/data/noaa/')  # Change WD.
    file = open(filename, 'wb') # Open destination file.
    ftp.retrbinary('RETR '+ filename, file.write) # Write file
    print('{} saved'.format(filename))
    file.close()                                 

    ## Read in csv and only filter to include US stations.
    stations = pd.read_csv("isd-history.csv")
    weccstations = stations[(stations['CTRY']=="US")]

    # Use spatial geometry to only keep points in wecc marine / terrestrial areas.
    geometry = [Point(xy) for xy in zip(weccstations['LON'], weccstations['LAT'])]
    crs = {'init' :'epsg:4326'}
    weccgeo = gp.GeoDataFrame(weccstations, crs=crs, geometry=geometry)
    
    ## get bbox of WECC to use to filter stations against
    ### FIX to remove hardcoded path var
    wecc_terr = '/home/ella/Desktop/Eagle Rock/Historical Data Platform /historical-obs-platform/test_platform/data/0_maps/WECC_Informational_MarineCoastal_Boundary_land.shp'
    t = gp.read_file(wecc_terr)
    ### FIX to remove hardcoded path var
    wecc_mar = '/home/ella/Desktop/Eagle Rock/Historical Data Platform /historical-obs-platform/test_platform/data/0_maps/WECC_Informational_MarineCoastal_Boundary_marine.shp'
    m = gp.read_file(wecc_mar)

    # Get terrestrial stations.
    weccgeo = weccgeo.to_crs(t.crs)
    terwecc = sjoin(weccgeo, t, how='left')
    terwecc = terwecc.dropna()

    # Get marine stations.
    marwecc = sjoin(weccgeo, m, how='left')
    marwecc = marwecc.dropna()
    
    # Join and remove duplicates.
    
    weccstations = (pd.concat([terwecc.iloc[:,:11], marwecc.iloc[:,:11]], ignore_index=True, sort =False)
            .drop_duplicates(['USAF','WBAN'], keep='first'))

    # Generate ID from USAF/WBAN combo for API call.
    weccstations['ISD-ID'] = weccstations['USAF']+weccstations['WBAN'].astype("str")

    # Reformat time strings for API call.
    weccstations['start_time'] = [datetime.strptime(str(i), '%Y%m%d').strftime('%Y-%m-%d') for i in weccstations['BEGIN']]
    weccstations['end_time'] = [datetime.strptime(str(i), '%Y%m%d').strftime('%Y-%m-%d') for i in weccstations['END']]

    # Filter stations returning data matching time specifications
    ### TO WRITE.

    weccstations.reset_index()
    return weccstations

#stations = get_wecc_stations()

#print(stations)

def get_isd_data(station_list, start_date = None): 
    
    # Remove depracated stations if filtering by time.
    if start_date is not None:
        try:
            station_list = station_list[station_list["end_time"] >= start_date] # Filter to ensure station is not depracated before time period of interest.
            print(station_list)
        except Exception as e:
            print("Oops!", e)
    
## Login.
## using ftplib
ftp = FTP('ftp.ncdc.noaa.gov')
ftp.login() # user anonymous, password anonymous
ftp.cwd('pub/data/noaa')  # Change WD.
pwd = ftp.pwd() # Get base file path.

## ASOS/AWOS FTP only goes until 2021.

# Get list of folders (by year) in main FTP folder.
years = ftp.nlst()
print(years)

    # Prepare variables for FTP call.
    



#get_isd_data(stations, start_date="2020-09-01")
#     for index, station in station_list.iterrows():
#         id = station['ISD-ID']
        
#         if start_time is None:
#             start = station['start_time']
#         else if start_time.:
#             start = start_time

#         if end_time is None:
#             end = station['end_time']
#         else:
#             end = end_time
        
#         print(station)
#         print(start, end, id)
#         #url = "https://www.ncei.noaa.gov/access/services/data/v1?dataset=global-hourly&stations={}&startDate={}&endDate={}&format=csv&includeStationName=true&includeStationLocation=1&units=metric".format(stations, start_time, end_time)

# get_isd_data(stations[1:10], start_time = '1920-01-12', end_time = '9201-01-10')
# #print(url)
# #response = requests.get(url)
# #print(response)
