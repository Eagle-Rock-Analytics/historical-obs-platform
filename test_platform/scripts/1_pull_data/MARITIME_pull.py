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
import re
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
years = list(map(str,range(1980,2023))) # Get list of years from 1980 to current year.

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

    ## Splits dataframe into respective networks
    maritime_network = weccstations[weccstations['NETWORK'] == 'MARITIME']
    ndbc_network = weccstations[weccstations['NETWORK'] == 'NDBC']

    ## Write stations to respective AWS bucket - requires separate buffers
    mar_buffer = StringIO()
    maritime_network.to_csv(mar_buffer)
    content = mar_buffer.getvalue()
    s3_cl.put_object(Bucket=bucket_name, Body=content, Key=directory_mar+"MARITIME_stations.csv")

    ndbc_buffer = StringIO()
    ndbc_network.to_csv(ndbc_buffer)
    content = ndbc_buffer.getvalue()
    s3_cl.put_object(Bucket=bucket_name, Body=content, Key=directory_ndbc+"NDBC_stations.csv")

    ## Purposely returning the full weccstations, instead of two separte dfs for get_maritime function
    return weccstations


## Read in MARITIME data using HTTP access.
# network = "NDBC" or "MARITIME"
def get_maritime(stations, bucket_name, network, years, get_all = True):

    # Set up error handling df.
    errors = {'Station ID':[], 'Time':[], 'Error':[]}

    # Set end time to be current time at beginning of download
    end_api = datetime.now().strftime('%Y%m%d%H%M')

    ## HTTP access instead of FTP: https://www.ndbc.noaa.gov/data/historical/stdmet/
    directory = "1_raw_wx/"+network+"/"
    dir_stations = stations.loc[stations['NETWORK']==network]

    ## Identifies which stations are owned by Canadian Dept of Environment and Climate Change (not qaqc'd by NDBC)
    # stations = stations.astype({"OWNER": "string"})
    canadian_owners = list(stations.loc[stations['OWNER'] == "CM"]['STATION_ID']) + list(stations.loc[stations['OWNER'] == "C"]["STATION_ID"]) # Canadian flags
    print('canada', canadian_owners) # testing

    for year in years:
        if year < (datetime.now().year):
            for filename in dir_stations['STATION_ID']:
                print(filename)
                try:  # Try to get station txt.gz.
                    if filename in canadian_owners: # Environment and Climate Change Canadian Moored buoy
                        # url = "https://www.meds-sdmm.dfo-mpo.gc.ca/alphapro/wave/waveshare/fbyears/C{}/c{}_{}.zip".format(str(filename), str(filename), str(year)) # individual year files
                        url = "https://www.meds-sdmm.dfo-mpo.gc.ca/alphapro/wave/waveshare/csvData/c{}_csv.zip".format(str(filename)) # all years in one file, including current year
                        # s3_obj = s3.Object(bucket_name, directory+"c{}_{}.zip") # individual year file format
                        s3_obj = s3.Object(bucket_name, directory+"c{}_csv.zip".format(str(filename))) # all years file format

                    else: # NOAA NDBC Buoy archive
                        url = "https://www.ndbc.noaa.gov/data/historical/stdmet/{}h{}.txt.gz".format(str(filename), str(year))
                        s3_obj = s3.Object(bucket_name, directory+"{}h{}.txt.gz".format(str(filename), str(year)))

                    with requests.get(url, stream=True) as r:
                        if r.status_code == 404: # Catches any stations that don't have specific years, could be cleaner potentially
                            next
                        elif r.status_code == 200: # If API call returns a response
                            if "RESPONSE_MESSAGE" in r.text: # If error response returned. Note that this is formatted differently depending on error type.
                                # Get error message and clean.
                                error_text = str(re.search("(RESPONSE_MESSAGE.*)",r.text).group(0)) # Get response message.
                                error_text = re.sub("RESPONSE_MESSAGE.: ", "", error_text)
                                error_text = re.sub(",.*", "", error_text)
                                error_text = re.sub('"', '', error_text)

                                # Append rows to dictionary
                                errors['Station ID'].append(dir_stations['STATION_ID'])
                                errors['Time'].append(end_api)
                                errors['Error'].append(error_text)
                                next
                            else:
                                s3_obj.put(Body=r.content)
                                print("Saving data for station {} for {}".format(filename, year)) # Nice for testing, remove for full run.
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
                                next
                            elif r.status_code == 200:
                                if "RESPONSE_MESSAGE" in r.text: # If error response returned. Note that this is formatted differently depending on error type.
                                    # Get error message and clean.
                                    error_text = str(re.search("(RESPONSE_MESSAGE.*)",r.text).group(0)) # Get response message.
                                    error_text = re.sub("RESPONSE_MESSAGE.: ", "", error_text)
                                    error_text = re.sub(",.*", "", error_text)
                                    error_text = re.sub('"', '', error_text)

                                    # Append rows to dictionary
                                    errors['Station ID'].append(dir_stations['STATION_ID'])
                                    errors['Time'].append(end_api)
                                    errors['Error'].append(error_text)
                                    next
                                else:
                                    s3_obj.put(Body=r.content)
                                    print("Saving data for station {} for {} {}".format(filename, x, year)) # Nice for testing, remove for full run.
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


## Comparison of which stations downloaded, updating the station_list csv with y/n to download column
## In general, if a station cannot be downloaded (has a N for download) it is an ocean-observing buoy ONLY, or no data is provided (optimization/testing buoy)
## MARITIME: 99 of 123 downloaded (24 are not meteorological)
## NDBC: 136 of 151 downloaded (16 are not meteorological)
def download_comparison(stations, bucket_name, network):

    ## There is a warning that is not relevant to our purposes, this turns it off
    pd.options.mode.chained_assignment = None  # default='warn'

    dir_bucket = s3.Bucket(bucket_name)

    ## Read in station list to compare against
    directory = "1_raw_wx/"+network+"/"
    dir_stations = stations.loc[stations['NETWORK']==network]

    files_downloaded = []
    for object_summary in dir_bucket.objects.filter(Prefix=directory):
        if object_summary.key[-7:] == ".txt.gz":
            files_downloaded.append(object_summary.key)

    downloaded_stns = set()
    for file in files_downloaded:
        dn_file = file.split("/")[2][:5]    # Grabs station_id from each filename
        downloaded_stns.add(dn_file)
    downloaded_stns = list(downloaded_stns)

    ## Adds download column so we can compare post full data pull, will get filled after full pull
    ## Mainly important for the oceanographic buoys that do not contain wx obs but are flagged as a part of WECC
    stn_yes = []

    ## Identifies whether a station in the station_list is a part of the downloaded_stns list
    dir_station_list = dir_stations['STATION_ID'].tolist() # All stations from station_list

    # Add flag
    dir_stations['Download'] = np.where(dir_stations['STATION_ID'].isin(downloaded_stns), "Y", "N")

    ## Reorders the indices in both stationlists
    ## Previously was the full index from station_table, so there was a mismatch in index and actual number of provided stations
    stn_idx = []
    for i in range(len(dir_stations.index)):
        stn_idx.append(i)
    dir_stations.index = stn_idx

    # ## Write stations to respective AWS bucket
    new_buffer = StringIO()
    dir_stations.to_csv(new_buffer)
    content = new_buffer.getvalue()
    s3_cl.put_object(Bucket=bucket_name, Body=content, Key=directory+"{}_stations.csv".format(network))

    # SAVE TO AWS - to do.
    print("{} station_list updated to reflect which weather-observing stations downloaded".format(network)) ## Function may take a while to run, useful to indicate completion

    return dir_stations

## ----------------------------------------------------------------------------------------------------------------------------
# To download all data, run:
network_to_run = "MARITIME" # "MARITIME" or "NDBC"
stations = get_maritime_station_ids(wecc_terr, wecc_mar, directory_mar, directory_ndbc)
get_maritime(stations, bucket_name, network_to_run, years = [2021], get_all = True)
download_comparison(stations, bucket_name, network_to_run)

## Full Pull Notes
## 0. Select either "MARITIME" or "NDBC" as network of choice to download for "network_to_run"
## 1. For first full data pull, set get_all = True
## 2. For all subsequent data pulls/update with newer data, set get_all = False
