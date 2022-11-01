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
from ftplib import FTP
from datetime import datetime, timezone
import pandas as pd
import csv
import boto3
from io import BytesIO, StringIO
import calc_pull
from shapely.geometry import Point
import geopandas as gp
from geopandas.tools import sjoin
from re import search
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

# Set paths to WECC shapefiles in AWS bucket.
wecc_terr = "s3://wecc-historical-wx/0_maps/WECC_Informational_MarineCoastal_Boundary_land.shp"
wecc_mar = "s3://wecc-historical-wx/0_maps/WECC_Informational_MarineCoastal_Boundary_marine.shp"

# Set envr variables
# years = list(map(str,range(1980,datetime.now().year+1))) # Get list of years from 1980 to current year.
years = list(map(str, range(1980,2021)))

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

    ## Reorders the indices in both stationlists
    ## Previously was the full index from station_table, so there was a mismatch in index and actual number of provided stations
    mar_idx = []
    for i in range(len(maritime_network.index)):
        mar_idx.append(i)
    maritime_network.index = mar_idx

    ndbc_idx = []
    for j in range(len(ndbc_network.index)):
        ndbc_idx.append(j)
    ndbc_network.index = ndbc_idx

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

    for year in years:
        for filename in dir_stations['STATION_ID']:
            if dir_stations['OWNER'] == "CM": # Environment and Climate Change Canadian Moored buoy
                url = "https://www.meds-sdmm.dfo-mpo.gc.ca/alphapro/wave/waveshare/fbyears/C{}/c{}_{}.zip".format(str(filename), str(filename), str(year))
            else: # NOAA NDBC Buoy archive
                url = "https://www.ndbc.noaa.gov/data/historical/stdmet/{}h{}.txt.gz".format(str(filename), str(year))

            # Try to get station txt.gz.
            try:
                #request = requests.get(url)
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

    ### Write errors to csv with respective networks
    error_buffer = StringIO()
    errors = pd.DataFrame(errors)
    errors.to_csv(error_buffer)
    content = error_buffer.getvalue()
    s3_cl.put_object(Bucket=bucket_name, Body=content,Key=directory+"errors_{}_{}.csv".format(network, end_api))


## Comparison of which stations downloaded, and adds y/n to stn_download column in station_list
## MARITIME: approx 93 of 123 can be downloaded (are meteorological observers, remaining are ocean only)
## NDBC: approx 122 of 151 can be downloaded
def download_comparison(stations, bucket_name, network):

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
    print(dir_station_list)

    dn_flag = []
    for all_stn in dir_station_list:
        if np.isin(all_stn, downloaded_stns, assume_unique=True) == True:
            dn_flag.append('Y')
            print("{} was downloaded".format(all_stn)) ## Useful for testing, can be deleted

        else:
            dn_flag.append('N')
            print("{} was not downloaded".format(all_stn)) ## Useful for testing, can be deleted

    dir_stations.insert(-1, 'STN_DOWNLOAD', dn_flag)

    print("{} station_list updated to reflect which stations downloaded".format(network)) ## Function may take a while to run, useful to indicate completion

    return dir_stations

## ----------------------------------------------------------------------------------------------------------------------------
# To download all data, run:
stations = get_maritime_station_ids(wecc_terr, wecc_mar, directory_mar, directory_ndbc)
# get_maritime(stations, bucket_name, "MARITIME", years = years, get_all = True)
download_comparison(stations, bucket_name, "MARITIME")

## Full Pull Notes
## 0. Select either "MARITIME" or "NDBC" as network of choice to download
## 1. For first full data pull, set get_all = True
## 2. For all subsequent data pulls/update with newer data, set get_all = False
