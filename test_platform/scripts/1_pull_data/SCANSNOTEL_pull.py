"""
This script scrapes SNOTEL and SCAN network data for ingestion into the Historical Observations Platform via SOAP API.
It may in the future be extended to pull other hourly data from networks available through the platform (e.g. BOR, USGS)
Approach:
(1) get_wecc_poly generates a bounding box from WECC shapefiles.
(3) get_SCAN_stations identifies station ids and provides metadata for SCAN stations in WECC.
(3) get SCAN_station_data saves a csv for each station ID, with the option to select a start date for downloads
    (defaults to station start date) and networks of interest, and primary or secondary sensors.
Inputs: API key, paths to WECC shapefiles, networks, start date (optional).
Outputs:
(1) Raw data for the SCAN and/or SNOTEL networks, all variables, all times. Organized by station.
(2) An error csv noting all station IDs where downloads failed, for each network.
"""

# Step 0: Environment set-up
import requests
import pandas as pd
from datetime import datetime
import numpy as np
import re
import boto3
from io import BytesIO, StringIO
from zeep import Client # For calling SOAP APIs
from zeep.helpers import serialize_object
import calc_pull

# Debug logging - for testing only.
#import logging.config
# logging.config.dictConfig({
#     'version': 1,
#     'formatters': {
#         'verbose': {
#             'format': '%(name)s: %(message)s'
#         }
#     },
#     'handlers': {
#         'console': {
#             'level': 'DEBUG',
#             'class': 'logging.StreamHandler',
#             'formatter': 'verbose',
#         },
#     },
#     'loggers': {
#         'zeep.transports': {
#             'level': 'DEBUG',
#             'propagate': True,
#             'handlers': ['console'],
#         },
#     }
# })

## Set AWS credentials
s3 = boto3.resource("s3")
s3_cl = boto3.client('s3') # for lower-level processes
bucket_name = "wecc-historical-wx"

# Connect to SCAN API
url = "https://wcc.sc.egov.usda.gov/awdbWebService/services?WSDL"
client = Client(url)


# Set envr variables
# Set paths to WECC shapefiles in AWS bucket.
wecc_terr = "s3://wecc-historical-wx/0_maps/WECC_Informational_MarineCoastal_Boundary_land.shp"
wecc_mar = "s3://wecc-historical-wx/0_maps/WECC_Informational_MarineCoastal_Boundary_marine.shp"

# Function: get SCAN station list in WECC region from SOAP API and save to AWS.
# Inputs: path to terrestrial and marine WECC shapefiles relative to home directory, AWS bucket name.
# Outputs: station metadata dataframe (also saved to AWS).
def get_SCAN_stations(terrpath, marpath, bucket_name, networks = None):
    try:
        t,m,bbox = calc_pull.get_wecc_poly(terrpath, marpath)
        bbox_api = bbox.loc[0,:].tolist() # [lonmin,latmin,lonmax,latmax]
        
        # Format bounding box for SOAP request
        request_data = {
            'minLatitude': str(bbox_api[1]),
            'maxLatitude': str(bbox_api[3]),
            'minLongitude': str(bbox_api[0]),
            'maxLongitude': str(bbox_api[2]),
            'logicalAnd': 'true' # and add other required parameter

        }

        # get station ids from bounding box
        stations = client.service.getStations(**request_data)

        # use station ids to return metadata
        request_data = {
            'stationTriplets': stations
        }

        # Call SOAP API
        metadata = client.service.getStationMetadataMultiple(**request_data)
        metadata = serialize_object(metadata) # Reformat object
        station_metadata = pd.DataFrame(metadata) # Convert to pandas df

        # If networks is none, set default networks.
        if networks is None:
            networks = ['SNTL', 'SCAN']

        for i in networks:
            subset = station_metadata[station_metadata['stationTriplet'].str.contains(i)]
            if i == 'SNTL': # Fix name differences between SNOTEL/SNTL
                i = "SNOTEL"
            subdir = "1_raw_wx/{}/".format(i)

            # Save to AWS
            csv_buffer = StringIO()
            subset.to_csv(csv_buffer, index = False)
            content = csv_buffer.getvalue()
            s3_cl.put_object(Bucket=bucket_name, Body=content, Key=subdir+"stationlist_{}.csv".format(i))

        return station_metadata

    except Exception as e:
        print("Error: {}".format(e))

# Function: download USDA station data using SOAP API.
# Data is organized by station, by sensor. 
# Inputs:
# (1-2): path to terrestrial and marine WECC shapefiles in AWS
# (3) bucket_name: name of AWS bucket
# (4) start_date: If none, download data starting on 1980-01-01. Otherwise, download data after start date (format "YYYY-MM-DD").
# (5) stations: optional. if not provided get_SCAN_stations is called. 
# # # Required to have columns named "stationTriplet" and "beginDate"
# (6) primary: if primary is False, download secondary sensor data.
# (7) networks: specify which network to download. If blank, get all networks: ["SNTL", "SCAN", "USGS", "BOR"]
# Outputs: CSV for each station saved to savedir, starting at start date and ending at current date.
def get_scan_station_data(terrpath, marpath, bucket_name, start_date = None, end_date = None, stations = None, primary = True, networks = None, fileext = None):

    # Set end time to be current time at beginning of download
    end_api = datetime.now().strftime('%Y%m%d%H%M')

    # Set up request parameters

    # Select network data if not provided
    if networks is None:
        networks = ["SNTL", "SCAN"] # Specify default networks
    
    # Get stations metadata if not provided
    if stations is None:
        stations = get_SCAN_stations(terrpath, marpath, bucket_name, networks)
    stationTriplets = list(stations['stationTriplet'])

    # Get start and end dates
    if start_date is None:
        start_date = "1980-01-01"
        
    if end_date is None:
        end_date = datetime.now().strftime('%Y-%m-%d') #yyyy-mm-dd
    
    # Select primary or secondary sensors
    if primary: # If primary True
        ordinal = "1"
    else:
        ordinal = "2"

    # Select codes for desired variables
    sensorcodes = ['TAVG', # Air temp avg
                    'TOBS', # Air temp observed
                    'PREC', # Precipitation accumulation
                    'PRCP', # Precipitation increment
                    'PRCPSA', # Precipitation increment snow-adjusted
                    'PRES', # Barometric pressure
                    'DPTP', # Dew point temperature
                    'RHUM', # Relative humidity
                    'RHUMV', # Relative humidity (avg),
                    'SRAD', # Solar radiation
                    'SRADV', # Solar radiation (avg)
                    'SRADT', # Solar radiation (total)
                    'PVPV', # Partial vapor pressure
                    'SVPV', # Saturated vapor pressure
                    'WDIRV', # Wind direction (avg)
                    'WDIR', # Wind direction
                    'WSPD', # Wind speed
                    'WSPDV'] # Wind speed (avg)

    
    for i in networks:
        # Set up error handling df.
        errors = {'Station ID':[], 'Time':[], 'Error':[]}

        networkTriplets = [k for k in stationTriplets if i in k]
        if i == "SNTL":
            directory = "1_raw_wx/SNOTEL/" # Maintain consistency across methods
        else:
            directory = "1_raw_wx/"+i+"/" # Define directory path for AWS
        print(directory)

        for j in networkTriplets: 
            try:
                df = pd.DataFrame() # Set up empty pd dataframe.
                # # For each sensor:
                for sensor in sensorcodes:
                    # Instantaneous:
                    # Required parameters:
                    if i in ['SNTL', 'SCAN']:
                        request_data = {
                                'stationTriplets': j, # Station ID(s) # Subset for testing.
                                'elementCd': sensor, # Sensor ID
                                'ordinal': ordinal, # Primary or secondary sensors
                                'beginDate': start_date, # Set earliest date
                                'endDate': end_date, # Set end date as today
                                'filter': 'ALL', # Select time of day to be all
                                'unitSystem': 'ENGLISH' # Convert to consistent unit system
                            }  
                        inst_data = client.service.getInstantaneousData(**request_data)
                    else:
                        request_data = {
                                'stationTriplets': j, # Station ID(s) # Subset for testing.
                                'elementCd': sensor, # Sensor ID
                                'ordinal': ordinal, # Primary or secondary sensors
                                'beginDate': start_date, # Set earliest date
                                'endDate': end_date # Set end date as today
                            }  
                        inst_data = client.service.getHourlyData(**request_data)
                        # Note: hourly data response not yet tested (only empty dfs returned),
                        # but should be used for non snotel/scan networks.

                    data = serialize_object(inst_data[0]) # Reformat object

                    station_data = pd.DataFrame(data['values']) # Convert to pandas df
                    if station_data.empty: # If no data returned for sensor, skip to next sensor
                        #print("no variable data found {}".format(sensor)) # For testing only.
                        continue 
                    else:
                        station_data = station_data.add_prefix(sensor+"_") # Add sensor code to column titles. 
                        if df.empty:
                            df = pd.concat([df, station_data])
                            timecol = sensor+'_time'
                            df['time'] = df[timecol]
                        else:
                            df = df.merge(station_data, left_on = 'time', right_on = sensor+"_time", how = 'outer')
                            timecol = sensor+'_time'
                            df['time'] = np.where(df['time'].isnull(), df[timecol], df['time']) # If any times in timecol are na, update with times from newest data.

                #If df empty, skip to next station
                if df.empty:
                    print("No data found for station {}. Skipping to next station.".format(j))
                    continue

                # Make time column first in order
                time = df.pop('time')
                df.insert(0, 'time', time)

                # Sort by time
                df = df.sort_values(by = 'time')
                
                # Write df to csv
                #print(df)
                csv_buffer = StringIO()
                df.to_csv(csv_buffer, index = False)
                content = csv_buffer.getvalue()
                if fileext is None:
                    s3_cl.put_object(Bucket=bucket_name, Body=content, Key=directory+"{}.csv".format(j))
                else:
                    s3_cl.put_object(Bucket=bucket_name, Body=content, Key=directory+"{}_{}.csv".format(j, fileext))
                print("Saved data for station {} in network {}".format(j, i))
            except Exception as e:
                print(e)
                errors['Station ID'].append(j)
                errors['Time'].append(end_api)
                errors['Error'].append(e)

        # Save errors to network folder
        csv_buffer = StringIO()
        errors = pd.DataFrame(errors)
        errors.to_csv(csv_buffer)
        content = csv_buffer.getvalue()
        if i == "SNTL":
            i == "SNOTEL" # Fix name for errors file
        s3_cl.put_object(Bucket=bucket_name, Body=content,Key=directory+"errors_{}_{}.csv".format(i, end_api))

    return errors
        

if __name__ == "__main__":
    get_scan_station_data(wecc_terr, wecc_mar, bucket_name, networks = ['SNTL'])



# Note: neither BOR nor USGS have hourly data for any of our variables of interest.
# At present, they are removed from our default networks of interest. However, they're left in the code, in the event that we
#  want to explore other networks/frequencies of data in the future.

# Note: code below downloads additional metadata lists (with different station names) for SCAN/SNOTEL stations.
# Refer to if names aren't matching correctly down the line.

# # SCAN
# url = 'https://wcc.sc.egov.usda.gov/nwcc/yearcount?network=scan&counttype=listwithdiscontinued&state='
# html = requests.get(url).content
# df_list = pd.read_html(html)
# scandf = df_list[-1]
# scandf['huccode'] = scandf['huc'].str.split("(").str[-1]
# scandf['huccode'] = scandf['huccode'].str.replace(")", "")
# # Remove all scan stations without huc code
# scandf['huccode'].replace('', np.nan, inplace=True)
# scandf.dropna(subset = ['huccode'], inplace = True) # Dropped one station.

# # SNOTEL
# url = 'https://wcc.sc.egov.usda.gov/nwcc/yearcount?network=sntl&counttype=listwithdiscontinued&state='
# html = requests.get(url).content
# df_list = pd.read_html(html)
# sntldf = df_list[-1]
# sntldf['huccode'] = sntldf['huc'].str.split("(").str[-1]
# sntldf['huccode'] = sntldf['huccode'].str.replace(")", "")
# # Remove all scan stations without huc code
# sntldf['huccode'].replace('', np.nan, inplace=True)
# sntldf.dropna(subset = ['huccode'], inplace = True) # Dropped two stations.
