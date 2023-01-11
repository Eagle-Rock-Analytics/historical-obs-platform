"""
This script downloads MARITIME & NDBC data from NDBC using http.
Approach:
(1) Download data using station list.
Inputs: bucket name in AWS, directory to save file to (folder path), range of years of data to download,
parameter to only download changed files (optional)
Outputs: Raw data for an individual network, all variables, all times. Organized by station, with 1 file per year.
Notes:
1. This function assumes users have configured the AWS CLI such that their access key / secret key pair are stored in ~/.aws/credentials.
See https://docs.aws.amazon.com/cli/latest/userguide/getting-started-install.html for guidance.
"""
## Load packages
import requests
from datetime import datetime, date
import pandas as pd
import boto3
from io import StringIO
import calc_pull
from shapely.geometry import Point
import geopandas as gp
from geopandas.tools import sjoin
import numpy as np


# Set AWS credentials
s3 = boto3.resource("s3")
s3_cl = boto3.client('s3') # for lower-level processes
bucket_name = 'wecc-historical-wx'

# Set paths to directories for each network
directory_mar = '1_raw_wx/MARITIME/'
directory_ndbc = '1_raw_wx/NDBC/'

# Set paths to WECC shapefiles in AWS bucket.
wecc_terr = "s3://wecc-historical-wx/0_maps/WECC_Informational_MarineCoastal_Boundary_land.shp"
wecc_mar = "s3://wecc-historical-wx/0_maps/WECC_Informational_MarineCoastal_Boundary_marine.shp"

# Set envr variables
# years = list(map(str,range(1980,2023))) # Get list of years from 1980 to current year.
years = list(map(str,range(1980,datetime.now().year+1))) # Get list of years from 1980 to current year.
## All years except the current year formatted as: stidhYYYY.txt.gz
## Current year formatted as: stidMYYYY.txt.gz

# ----------------------------------------------------------------------------------------------------------------------------

# Step 0: Get list of IDs to download.
### Get all stations that start with "46" (Pacific Ocean) and all lettered (all CMAN) stations.
### Note: Gliders (non-stationary) location reported here as 30 N -90 W
def get_maritime_station_ids(terrpath, marpath, directory_mar, directory_ndbc):
    """Returns list of station ids as CSV"""
    url = 'https://www.ndbc.noaa.gov/data/stations/station_table.txt'
    r = requests.get(url)
    lines = r.content.split(b'\n') # split by line
    df = []
    for line in lines:
        row = line.split(b'|') # Split on |
        row = [x.decode('utf-8') for x in row] # converts to string
        row = list(map(str.strip, row))
        df.append(row)

    stations = pd.DataFrame(columns=['STATION_ID', 'OWNER', 'TTYPE', 'HULL', 'NAME', 'PAYLOAD', 'LOCATION', 'TIMEZONE', 'FORECAST', 'NOTE'], data=df[1:])
    stations = stations.dropna()
    stations = stations.drop(columns=['TTYPE', 'HULL', 'PAYLOAD', 'FORECAST', 'TIMEZONE'])

    ## Use spatial geometry to only keep points within WECC marine/terrestrial region (should only be marine)
    ## Locations listed as one value: "30.000 N 90.000 W" as an example, split on spaces
    mar_lon = []
    mar_lat = []
    location = stations['LOCATION']
    for line in location:
        line = line.split(' ')
        if "S" in line:
            mar_lat.append(-1 * float(line[0]))
        else:
            mar_lat.append(float(line[0]))

        if "E" in line:
            mar_lon.append(float(line[2]))
        else:
            mar_lon.append(-1 * float(line[2]))

    stations['LATITUDE'] = mar_lat
    stations['LONGITUDE'] = mar_lon

    geometry = [Point(xy) for xy in zip(stations['LONGITUDE'], stations['LATITUDE'])] # Zip lat lon coords
    weccgeo = gp.GeoDataFrame(stations, crs="EPSG:4326", geometry=geometry) # converts to geodataframe

    ## get bbox of WECC to use to filter stations against
    t, m, bbox = calc_pull.get_wecc_poly(terrpath, marpath) # Call get_wecc_poly.

    ## Get terrestrial stations.
    weccgeo = weccgeo.to_crs(t.crs) # Convert to CRS of terrestrial stations.
    terwecc = sjoin(weccgeo.dropna(), t, how='left') # Only keep stations in terrestrial WECC region.
    terwecc = terwecc.dropna() # Drop empty rows.

    ## Get marine stations.
    marwecc = sjoin(weccgeo.dropna(), m, how='left') # Only keep stations in marine WECC region.
    marwecc = marwecc.dropna() # Drop empty rows.

    weccstations = (pd.concat([terwecc, marwecc], ignore_index = True, sort=False)
        .drop_duplicates(['STATION_ID'], keep='first'))

    ## drop columns and rename in_wecc to in_terr
    weccstations.drop(['OBJECTID_1', 'OBJECTID', 'Shape_Leng','geometry', 'FID_WECC_B', 'BUFF_DIST', 'index_right'], axis = 1, inplace = True)
    weccstations.rename(columns = {'in_WECC':'in_terr_wecc', 'in_marine':'in_mar_wecc'}, inplace = True)

    ## Identify which buoys are NDBC and which are CMAN/MARITIME
    network = []
    for item in weccstations['STATION_ID']:
        if item[:2] == "46":
            network.append('NDBC')
        else:
            network.append('MARITIME')

    weccstations['NETWORK'] = network

    ## Moves Note column to last column because of weird lat-lon column overwriting -- metadata issue for cleaning
    weccstations = weccstations.reindex(columns = [col for col in weccstations.columns if col != 'NOTE'] + ['NOTE'])

    ## Splits dataframe into respective networks
    maritime_network = weccstations[weccstations['NETWORK'] == 'MARITIME']
    ndbc_network = weccstations[weccstations['NETWORK'] == 'NDBC']

    ## Write stations to respective AWS bucket - requires separate buffers
    mar_buffer = StringIO()
    maritime_network.to_csv(mar_buffer)
    content = mar_buffer.getvalue()
    s3_cl.put_object(Bucket=bucket_name, Body=content, Key=directory_mar+"stationlist_MARITIME.csv")

    ndbc_buffer = StringIO()
    ndbc_network.to_csv(ndbc_buffer)
    content = ndbc_buffer.getvalue()
    s3_cl.put_object(Bucket=bucket_name, Body=content, Key=directory_ndbc+"stationlist_NDBC.csv")

    ## Purposely returning the full weccstations, instead of two separte dfs for get_maritime function
    return weccstations


## Read in MARITIME data using HTTP access.
## Removing "get_all" functionality for now, but may need to be redesigned in future for current year/preliminary data
def get_maritime(stations, bucket_name, network, years = None):

    # Set up error handling df.
    errors = {'Station ID':[], 'Time':[], 'Error':[]}

    # Set end time to be current time at beginning of download
    end_api = datetime.now().strftime('%Y%m%d%H%M')

    ## HTTP access instead of FTP: https://www.ndbc.noaa.gov/data/historical/stdmet/
    directory = "1_raw_wx/"+network+"/"
    dir_stations = stations.loc[stations['NETWORK']==network]

    ## Default for years
    if years is None:
        years = list(map(str,range(1980,datetime.now().year+1))) # Get list of years from 1980 to current year
    else:
        years = years

    ## Identifies which stations are owned by Canadian Dept of Environment and Climate Change (not qaqc'd by NDBC)
    canadian_owners = list(stations.loc[stations['OWNER'] == "CM"]['STATION_ID']) + list(stations.loc[stations['OWNER'] == "C"]["STATION_ID"]) # Canadian flags

    for year in years:
        if year < str(datetime.now().year):
            for filename in dir_stations['STATION_ID']:
                try:  # Try to get station txt.gz.
                    if filename in canadian_owners: # Environment and Climate Change Canadian Moored buoy
                        url = "https://www.meds-sdmm.dfo-mpo.gc.ca/alphapro/wave/waveshare/csvData/c{}_csv.zip".format(filename) # all years in one file, including current year
                        s3_obj = s3.Object(bucket_name, directory+"{}_csv.zip".format(filename)) # all years file format

                    else: # NOAA NDBC Buoy archive
                        url = "https://www.ndbc.noaa.gov/data/historical/stdmet/{}h{}.txt.gz".format(filename, year)
                        s3_obj = s3.Object(bucket_name, directory+"{}h{}.txt.gz".format(filename, year))

                    with requests.get(url, stream=True) as r:
                        if r.status_code == 404: # Catches any stations that don't have specific years, could be cleaner potentially
                            continue
                        elif r.status_code == 200:
                            s3_obj.put(Body=r.content)
                            ## Note: The Canadian buoy all-years-file still gets downloaded/overwritten for len(years) times, where it could just be downloaded once
                            print("Saving data for station {} for {}".format(filename, year)) 
                        else:
                            errors['Station ID'].append(filename)
                            errors['Time'].append(end_api)
                            errors['Error'].append(r.status_code)
                            print("Error: {}".format(r.status_code))

                except Exception as e:
                    print("Error: {}".format(e))

        else: # only applicable for the NDBC stored data (non-Canadian source)
            for filename in dir_stations['STATION_ID']:
                try:
                    year = int(year)
                    start_month = date(year, 1, 1) # Jan
                    current_month = date.today()
                    i = current_month.month

                    while i >= start_month.month:
                        current_date = date(year, i, 1)
                        x = current_date.strftime("%b") # NDBC folder names is 3-letter month code
                        url = "https://www.ndbc.noaa.gov/data/stdmet/{}/{}{}{}.txt.gz".format(x, filename, i, year)
                        s3_obj = s3.Object(bucket_name, directory+"{}{}{}.txt.gz".format(filename, i, year))

                        with requests.get(url, stream=True) as r:
                            if r.status_code == 404:
                                pass
                            elif r.status_code == 200:
                                s3_obj.put(Body=r.content)
                                print("Saving data for station {} for {} {}".format(filename, x, year)) 
                            else:
                                errors['Station ID'].append(filename)
                                errors['Time'].append(end_api)
                                errors['Error'].append(r.status_code)
                                print("Error: {}".format(r.status_code))

                        i -= 1

                except Exception as e:
                    print("Error: {}".format(e))

    ### Write errors to csv with respective networks
    error_buffer = StringIO()
    errors = pd.DataFrame(errors)
    errors.to_csv(error_buffer)
    content = error_buffer.getvalue()
    s3_cl.put_object(Bucket=bucket_name, Body=content,Key=directory+"errors_{}_{}.csv".format(network, end_api))



## Read in MARITIME data using HTTP access.
# As above, but with start and end time filtering functionality. Note here we can only filter download by start and end month, so we download all
# data for month even if only a few days are requested.
def get_maritime_update(stations, bucket_name, network, start_date = None, end_date = None):

    # Set up error handling df.
    errors = {'Station ID':[], 'Time':[], 'Error':[]}

    # Set end time to be current time at beginning of download
    end_api = datetime.now().strftime('%Y%m%d%H%M')

    ## HTTP access instead of FTP: https://www.ndbc.noaa.gov/data/historical/stdmet/
    directory = "1_raw_wx/"+network+"/"
    dir_stations = stations.loc[stations['NETWORK']==network]

    # Get range of years to download
    if start_date is None:
        start_month = date(year, 1, 1).month # Jan
        if end_date is None:
            years = list(map(str,range(1980, int(datetime.now().year)+1)))
        else:
            years = list(map(str,range(1980, int(end_date[0:4])+1)))
            end_month = datetime.strptime(end_date, '%Y-%m-%d').month
    else:
        start_year = int(start_date[0:4])
        start_month = datetime.strptime(start_date, '%Y-%m-%d').month
        if end_date is None:
            years = list(map(str,range(start_year, int(datetime.now().year)+1)))
        else:
            years = list(map(str,range(start_year, int(end_date[0:4])+1)))
            end_month = datetime.strptime(end_date, '%Y-%m-%d').month


    ### IN PROGRESS
    # Set up cross year flag
    if int(datetime.now().strftime('%j'))<=45: # if download occurring in first 45 days of year
        if (str(datetime.now().year-1)) in years: # and includes previous year's data
            
            month_list = pd.period_range(start=start_date, end=end_date, freq='M')
            month_list = [month.strftime("%b") for month in month_list]
            
            cross_year = 1 # Set flag to be true
    
    ## Identifies which stations are owned by Canadian Dept of Environment and Climate Change (not qaqc'd by NDBC)
    canadian_owners = list(stations.loc[stations['OWNER'] == "CM"]['STATION_ID']) + list(stations.loc[stations['OWNER'] == "C"]["STATION_ID"]) # Canadian flags

    if cross_year == 1:
        if start_date is not None and end_date is not None:
            month_list = pd.period_range(start=start_date, end=end_date, freq='M')
            month_list = [month.strftime("%b") for month in month_list]
            print(month_list)
            exit()

    ### IN PROGRESS ^

    for year in years:
        if year < str(datetime.now().year): # If not in current year
            for filename in dir_stations['STATION_ID']:
                try:  # Try to get station txt.gz.
                    if filename in canadian_owners: # Environment and Climate Change Canadian Moored buoy
                        url = "https://www.meds-sdmm.dfo-mpo.gc.ca/alphapro/wave/waveshare/csvData/c{}_csv.zip".format(filename) # all years in one file, including current year
                        s3_obj = s3.Object(bucket_name, directory+"{}_csv.zip".format(filename)) # all years file format

                    else: # NOAA NDBC Buoy archive
                        url = "https://www.ndbc.noaa.gov/data/historical/stdmet/{}h{}.txt.gz".format(filename, year)
                        s3_obj = s3.Object(bucket_name, directory+"{}h{}.txt.gz".format(filename, year))

                    with requests.get(url, stream=True) as r:
                        if r.status_code == 404: # Catches any stations that don't have specific years, could be cleaner potentially
                            continue
                        elif r.status_code == 200:
                            s3_obj.put(Body=r.content)
                            ## Note: The Canadian buoy all-years-file still gets downloaded/overwritten for len(years) times, where it could just be downloaded once
                            print("Saving data for station {} for {}".format(filename, year)) 
                        else:
                            errors['Station ID'].append(filename)
                            errors['Time'].append(end_api)
                            errors['Error'].append(r.status_code)
                            print("Error: {}".format(r.status_code))

                except Exception as e:
                    print("Error: {}".format(e))


        else: # only applicable for the NDBC stored data (non-Canadian source)
            for filename in dir_stations['STATION_ID']:
                try:
                        
                    # For each month in the current year
                    start_month = date(year, 1, 1) # Jan
                    current_month = date.today()
                    i = current_month.month

                    while i >= start_month.month:
                        current_date = date(year, i, 1)
                        x = current_date.strftime("%b") # NDBC folder names is 3-letter month code
                        url = "https://www.ndbc.noaa.gov/data/stdmet/{}/{}{}{}.txt.gz".format(x, filename, i, year)
                        s3_obj = s3.Object(bucket_name, directory+"{}{}{}.txt.gz".format(filename, i, year))

                        with requests.get(url, stream=True) as r:
                            if r.status_code == 404:
                                pass
                            elif r.status_code == 200:
                                s3_obj.put(Body=r.content)
                                print("Saving data for station {} for {} {}".format(filename, x, year)) 
                            else:
                                errors['Station ID'].append(filename)
                                errors['Time'].append(end_api)
                                errors['Error'].append(r.status_code)
                                print("Error: {}".format(r.status_code))

                        i -= 1

                except Exception as e:
                    print("Error: {}".format(e))


    ### Write errors to csv with respective networks
    error_buffer = StringIO()
    errors = pd.DataFrame(errors)
    errors.to_csv(error_buffer)
    content = error_buffer.getvalue()
    s3_cl.put_object(Bucket=bucket_name, Body=content,Key=directory+"errors_{}_{}.csv".format(network, end_api))

## ----------------------------------------------------------------------------------------------------------------------------
# To download all data, run:
if __name__ == "__main__":
    network_to_run = "NDBC" # "MARITIME" or "NDBC"
    stations = get_maritime_station_ids(wecc_terr, wecc_mar, directory_mar, directory_ndbc)
    get_maritime(stations, bucket_name, network_to_run, years = None)
    
## Full Pull Notes
## 1. Select either "MARITIME" or "NDBC" as network of choice to download for "network_to_run"
