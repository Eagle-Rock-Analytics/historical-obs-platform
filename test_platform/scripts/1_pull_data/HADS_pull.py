"""
This script downloads HADS data directly from https://www.ncei.noaa.gov/data/nws-hydrometeorological-automated-data-system/archive/.
Approach:
(1) Get station list (does not need to be re-run constantly)
(2) Download data using station list.
Inputs: bucket name in AWS, directory to save file to (folder path), station list (optional), start date of file pull (optional),
parameter to only download changed files (optional)
Outputs: Raw data for an individual network, all variables, all times. Organized by station, with 1 station per day per year.

Notes:
1. The file for each station-year is updated daily for the current year. 
To pull real-time data, we may want to write just an API call with date ranges and stations and update the most recent year folder only. 
This is a separate function/branch.
2. This function assumes users have configured the AWS CLI such that their access key / secret key pair are stored in ~/.aws/credentials.
See https://docs.aws.amazon.com/cli/latest/userguide/getting-started-install.html for guidance.
"""

# Step 0: Environment set-up
import requests
import pandas as pd
from datetime import datetime, timezone
import re
import boto3
from io import StringIO
import calc_pull
import geopandas as gp
from shapely.geometry import Point
from geopandas.tools import sjoin
import requests 
from bs4 import BeautifulSoup 

## Set AWS credentials
s3 = boto3.resource("s3")
s3_cl = boto3.client('s3') # for lower-level processes
bucket_name = "wecc-historical-wx"
directory = "1_raw_wx/HADS/"

try:
    import config # Import API keys.
except:
    print("Missing config.py file with API token. Make file if necessary.")
    exit()

# Set envr variables
# Set paths to WECC shapefiles in AWS bucket.
wecc_terr = "s3://wecc-historical-wx/0_maps/WECC_Informational_MarineCoastal_Boundary_land.shp"
wecc_mar = "s3://wecc-historical-wx/0_maps/WECC_Informational_MarineCoastal_Boundary_marine.shp"

# Downloading SHEF files directly from NOAA (1997-2022)

# Function to identify HADS stations in WECC and save stations as csv to AWS folder.
# Pulls in WECC shapefile, returns list of all stations and IDs. 
# Inputs: path to terrestrial WECC shapefile, path to marine WECC file. 
# Both paths given relative to home directory for git project.
# This function should be run occasionally.
def get_hads_stations(terrpath, marpath):   
    ## Get list of HADS stations as CSV
    url = 'https://hads.ncep.noaa.gov/compressed_defs/all_dcp_defs.txt'
    r = requests.get(url)
    lines = r.content.split(b'\n') # Split by line
    df = []
    for line in lines:
        row = line.split(b'|') # Split on |
        row = row[0:10] # Only keep first 10 columns
        row = [x.decode('utf-8') for x in row] # Convert to string
        row = list(map(str.strip, row))
        df.append(row)

    # Create pandas dataframe
    stations = pd.DataFrame(columns=['GOES NESDIS ID','NWSLI','DCP Owner','State Location','Hydrologic Service Area','Latitude','Longitude','Initial Daily Transmission Time (UTC)','DCP Transmission Interval (Minutes)','DCP Location Name'], data=df)
    stations = stations.dropna()
    
    # # Use spatial geometry to only keep points in wecc marine / terrestrial areas.
    stations['Longitude'] = [calc_pull._lon_dms_to_dd(i) for i in stations['Longitude']]
    stations['Latitude'] = [calc_pull._lat_dms_to_dd(i) for i in stations['Latitude']]

    geometry = [Point(xy) for xy in zip(stations['Longitude'], stations['Latitude'])] # Zip lat lon coords.
    weccgeo = gp.GeoDataFrame(stations, crs='EPSG:4326', geometry=geometry) # Convert to geodataframe.
    
    ## get bbox of WECC to use to filter stations against
    t, m, bbox = calc_pull.get_wecc_poly(terrpath, marpath) # Call get_wecc_poly.

    # Get terrestrial stations.
    weccgeo = weccgeo.to_crs(t.crs) # Convert to CRS of terrestrial stations.
    terwecc = sjoin(weccgeo.dropna(), t, how='left') # Only keep stations in terrestrial WECC region.
    terwecc = terwecc.dropna() # Drop empty rows.

    # Get marine stations.
    marwecc = sjoin(weccgeo.dropna(), m, how='left') # Only keep stations in marine WECC region.
    marwecc = marwecc.dropna() # Drop empty rows.
    
    
    #Join and remove duplicates using GOES NESDIS ID as combined unique identifier.
    weccstations = (pd.concat([terwecc, marwecc], ignore_index=True, sort =False)
           .drop_duplicates(['GOES NESDIS ID'], keep='first'))

    # Drop columns
    weccstations.drop(['OBJECTID_1', 'OBJECTID', 'Shape_Leng','geometry', 'FID_WECC_B', 'BUFF_DIST', 'index_right'], axis = 1, inplace = True)
        
    # Rename in_wecc to in_terr
    weccstations.rename(columns = {'in_WECC':'in_terr_wecc', 'in_marine':'in_mar_wecc'}, inplace = True)

    # Write stations to AWS bucket
    wecc_buffer = StringIO()
    weccstations.to_csv(wecc_buffer)
    content = wecc_buffer.getvalue()
    s3_cl.put_object(Bucket=bucket_name, Body=content,Key=directory+"HADS_stations.csv")
    return weccstations

test = get_hads_stations(wecc_terr, wecc_mar)
print(test) # For testing.

# Function: Takes URL of given website and returns all file download links contained on it.
# Input: one url
# Output: zipped list of download urls and date last modified.
def get_file_links(url): 
      
    # create response object 
    r = requests.get(url) 
      
    # create beautiful-soup object 
    soup = BeautifulSoup(r.content,'html5lib') 
      
    # find all links on web-page 
    links = soup.findAll('a') 
  
    # filter the link sending with dat.gz
    file_links = [url + link['href'] for link in links if link['href'].endswith('dat.gz')] 

    # get all 'last modified' dates
    dates = soup.find_all(string=re.compile(r'\d{4}-\d{2}-\d{2}')) # Get all dates from website.
    dates = [i.strip() for i in dates]
    links = list(zip(file_links, dates))
    return links 

# Function: takes input list of links generated by get_file_links() and saves files to AWS bucket.
# Inputs: links from get_file_links(), name of AWS bucket, and folder path to save directory.
def link_to_aws(links, bucket_name, directory): 
    s3 = boto3.resource("s3")
    for i, (file_link, dates) in enumerate(links): 
  
        '''iterate through all links in links 
        and download them one by one'''
        # UPDATE TO REFLECT ZIPPED LIST - TO DO.  
        # obtain filename by splitting url and getting 
        # last string 
        file_name = file_link.split('/')[-1] 
        s3_obj = s3.Object(bucket_name, directory+file_name)
    
        # create response object 
        with requests.get(file_link, stream = True) as r:
            if r.status_code == 200: # If API call returns a response
                s3_obj.put(Body=r.content)
                print("Saving file {}".format(file_name))
    return

# Function: query NCEI server for HADS data and download zipped files.
# Run this one time to get all historical data or to update changed files for all years.
# Inputs: 
# bucket_name: name of AWS bucket
# directory: folder path within bucket
# years: format 'YYYY" (optional)
# get_all: True or False. If False, only download files whose last edit date is newer than
#  the most recent files downloaded in the save folder. Only use to update a complete set of files.
def get_hads_dat(bucket_name, directory, start_date = None, get_all = True):
    # Set up error handling
    errors = {'Date':[], 'Time':[], 'Error':[]}
    end_api = datetime.now().strftime('%Y%m%d%H%M') # Set end time to be current time at beginning of download

    # Set up AWS to write to bucket.
    s3 = boto3.client('s3')

    try:
        objects = s3.list_objects(Bucket=bucket_name,Prefix = directory)
        all = objects['Contents']     
        # Get date of last edited file   
        latest = max(all, key=lambda x: x['LastModified'])
        last_edit_time = latest['LastModified']
        # Get list of all file names
        alreadysaved = []
        for item in all:
            files = item['Key']
            alreadysaved.append(files) 
        alreadysaved = [ele.replace(directory, '') for ele in alreadysaved]
    except:
        get_all = True # If folder empty or there's an error with the "last downloaded" metadata, redownload all data.

    # Set up filtering by time.
    if start_date is None:
        years = range(1997, 2023)
    else:
        start_date = int(start_date)
        years = range(start_date, 2023)
    
    for i in years:
        try:
            yearurl = 'https://www.ncei.noaa.gov/data/nws-hydrometeorological-automated-data-system/archive/{}/'.format(str(i))
            links = get_file_links(yearurl)
            if get_all != True:
                try:
                    # Get files in AWS folder.
                    objects = s3.list_objects(Bucket=bucket_name,Prefix = directory)
                    all = objects['Contents']     
                    # Get date of last edited file   
                    latest = max(all, key=lambda x: x['LastModified'])
                    last_edit_time = latest['LastModified']
                    # Get list of all file names
                    alreadysaved = []
                    for item in all:
                        files = item['Key']
                        alreadysaved.append(files) 
                    alreadysaved = [ele.replace(directory, '') for ele in alreadysaved]
                    
                    # Filter links by file names already in folder, and files updated since last download.
                    link_sub = []
                    date_sub = []
                    for k, (link, date) in enumerate(links):
                        if link.split("/")[-1] in alreadysaved:
                            date = datetime.strptime(date, "%Y-%m-%d %H:%M").replace(tzinfo=timezone.utc)
                            if (date>last_edit_time):
                                link_sub.append(link)
                                date_sub.append(date)
                            else:
                                print("{} already saved".format(link.split("/")[-1]))
                                continue
                        else:
                            link_sub.append(link)
                            date_sub.append(date)
                    link_sub = list(zip(link_sub, date_sub))
                    link_to_aws(link_sub, bucket_name, directory)

                except:
                    get_all = True # If folder empty or there's an error with the "last downloaded" metadata, redownload all data.
            if get_all == True:
                link_to_aws(links, bucket_name, directory)
        except Exception as e:
            print("Error in downloading year {}: {}". format(i, e))
            errors['Date'].append(i)
            errors['Time'].append(end_api)
            errors['Error'].append(e)
    
    #Write errors to csv
    csv_buffer = StringIO()
    errors = pd.DataFrame(errors)
    errors.to_csv(csv_buffer)
    content = csv_buffer.getvalue()
    s3.put_object(Bucket=bucket_name, Body=content,Key=directory+"errors_hads_{}.csv".format(end_api))

    return



get_hads_dat(bucket_name, directory, start_date = '2022', get_all = False)
# Tested get_all function: works
# Tested start_date function: works
# Saves to correct directory.
# Error saves correctly.


"""
Alternative method: depracated. Uses Synoptic API and MADIS data, begins in 2003. 
Recommendation: Don't use unless serious problems found when cleaning above data.
"""

# Function: return metadata for stations from synoptic API, filtering stations based on bounding box generated by get_wecc_poly.
# Inputs: Synoptic API public token (imported from config.py), path to terrestrial and marine WECC shapefiles relative to home directory.
# Outputs: Dataframe with list of station IDs and start date of station.
def get_meso_metadata(token, terrpath, marpath):
    try:
        t,m,bbox = calc_pull.get_wecc_poly(terrpath, marpath)
        bbox_api = bbox.loc[0,:].tolist() # [lonmin,latmin,lonmax,latmax]
        bbox_api = ','.join([str(elem) for elem in bbox_api])
        # Access station metadata to get list of IDs in bbox and network
        # Using: https://developers.synopticdata.com/mesonet/v2/stations/timeseries/
        # HADS data network ID in synoptic is 106 (see below)
        url = "https://api.synopticdata.com/v2/stations/metadata?token={}&network=106&bbox={}&recent=20&output=json".format(token, bbox_api)
        request = requests.get(url).json()
        # print(request)
        ids = []
        for each in request['STATION']:
            ids.append([each['STID'], each['PERIOD_OF_RECORD']['start']]) # Keep station ID and start date

        ids = pd.DataFrame(ids, columns = ['STID', 'start']).sort_values('start') # Sort by start date (note some stations return 'None' here)
        # print(ids)
        # Reformat date to match API format
        ids['start'] = pd.to_datetime(ids['start'], format='%Y-%m-%dT%H:%M:%SZ')
        ids['start'] = ids['start'].dt.strftime('%Y%m%d%H%M')
        return ids
    except Exception as e:
        print("Error: {}".format(e))

# Function: download HADS station data from Synoptic API.
# Inputs:
# (1) token: Synoptic API token, stored in config.py file
# (2) ids: Takes dataframe with two columns, station ID and startdate. These are the stations to be downloaded, generated from get_meso_metadata.
# (3) bucket_name: name of AWS bucket
# (4) directory: folder path in AWS bucket
# (5) start_date: If none, download data starting on 1980-01-01. Otherwise, download data after start date ("YYYYMMDDHHMM").
# (6) options: timeout = True will identify and download any station data that timed out the API request.
# Outputs: CSV for each station saved to savedir, starting at start date and ending at current date.
def get_hads_station_csv(token, ids, bucket_name, directory, start_date = None, **options):

    # Set up error handling df.
    errors = {'Station ID':[], 'Time':[], 'Error':[]}

    # Set end time to be current time at beginning of download
    end_api = datetime.now().strftime('%Y%m%d%H%M')

    for index, id in ids.iterrows(): # For each station
        # Set start date
        if start_date is None:
            if id["start"] == "NaN" or pd.isna(id['start']): # Adding error handling for NaN dates.
                start_api = "198001010000" # If any of the stations in the chunk have a start time of NA, pull data from 01-01-1980 OH:OM and on.
                # To check: all stations seen so far with NaN return empty data. If true universally, can skip these instead of saving.
            else:
                start_api = id["start"]
        else:
            start_api = start_date

        # Generate URL
        # Note: decision here to use full flag suite of MesoWest and Synoptic data.
        # See Data Checks section here for more information: https://developers.synopticdata.com/mesonet/v2/stations/timeseries/
        url = "https://api.synopticdata.com/v2/stations/timeseries?token={}&stid={}&start={}&end={}&output=csv&qc=on&qc_remove_data=off&qc_flags=on&qc_checks=synopticlabs,mesowest".format(token, id['STID'], start_api, end_api)
        # print(url) # For testing.

        # Try to get station csv.
        try:
            #request = requests.get(url)
            s3_obj = s3.Object(bucket_name, directory+"{}.csv".format(id["STID"]))

            # If **options timeout = True, save file as STID_2.csv
            if options.get("timeout") == True:
                s3_obj = s3.Object(bucket_name, directory+"{}_2.csv".format(id["STID"]))

            with requests.get(url, stream=True) as r:
                if r.status_code == 200: # If API call returns a response
                    if "RESPONSE_MESSAGE" in r.text: # If error response returned. Note that this is formatted differently depending on error type.
                        # Get error message and clean.
                        error_text = str(re.search("(RESPONSE_MESSAGE.*)",r.text).group(0)) # Get response message.
                        error_text = re.sub("RESPONSE_MESSAGE.: ", "", error_text)
                        error_text = re.sub(",.*", "", error_text)
                        error_text = re.sub('"', '', error_text)

                        # Append rows to dictionary
                        errors['Station ID'].append(id['STID'])
                        errors['Time'].append(end_api)
                        errors['Error'].append(error_text)
                        next
                    else:
                        s3_obj.put(Body=r.content)
                        print("Saving data for station {}".format(id["STID"])) # Nice for testing, remove for full run.

                else:
                    errors['Station ID'].append(id['STID'])
                    errors['Time'].append(end_api)
                    errors['Error'].append(r.status_code)
                    print("Error: {}".format(r.status_code))

        except Exception as e:
            print("Error: {}".format(e))

    # Write errors to csv for AWS
    csv_buffer_err = StringIO()
    errors = pd.DataFrame(errors)
    errors.to_csv(csv_buffer_err)
    content = csv_buffer_err.getvalue()
    s3_cl.put_object(Bucket=bucket_name, Body=content, Key=directory+"errors_hads_{}.csv".format(end_api))

# Quality control: if any files return status 408 error, split request into smaller requests and re-run.
# Note: this approach assumes no file will need more than 2 splits. Test this when fuller data downloaded.
def get_hads_station_timeout_csv(token, bucket_name, directory):
    ids_split = []
    for item in s3.Bucket(bucket_name).objects.all():
        files = item.key
    files = list(filter(lambda f: f.endswith(".csv"), files)) # Get list of file names

    for file in files:
        with open(bucket_name + directory + file,'r') as file:
            data = file.readlines()
            lastRow = data[-1]
            if 'Timeout' in lastRow: # If timeout recorded
                lastrealrow = data[-2] # Get last row of recorded data
                station = lastrealrow.split(",")[0] # Get station ID
                time = lastrealrow.split(",")[1] # Get last completed timestamp
                ids_split.append([station, time]) # Add to dataframe

    ids_split = pd.DataFrame(ids_split, columns = ['STID', 'start'])
    ids_split['start'] = pd.to_datetime(ids_split['start'], format='%Y-%m-%dT%H:%M:%SZ')
    ids_split['start'] = ids_split['start'].dt.strftime('%Y%m%d%H%M')
    if ids_split.empty is False:
        print(ids_split)
        get_hads_station_csv(token, bucket_name, directory, ids = ids_split, timeout = True)

        # Check to see if any of the split files needs to be split again.
        for item in s3.Bucket(bucket_name).objects.all():
            files = item.key
        files = list(filter(lambda f: f.endswith("_2.csv"), files)) # Get list of file names

        for file in files:
            with open(bucket_name + directory + file,'r') as file:
                data = file.readlines()
                lastRow = data[-1]
                if 'Timeout' in lastRow: # If timeout recorded
                    lastrealrow = data[-2] # Get last row of recorded data
                    station = lastrealrow.split(",")[0] # Get station ID
                    time = lastrealrow.split(",")[1] # Get last completed timestamp
                    ids_split.append([station, time]) # Add to dataframe
                    if ids_split.empty is False:
                        print("Attention!: Run this script again on _2.csv files.")

    elif ids_split.empty is True:
        return

# Run script.
#ids = get_meso_metadata(token = config.token, terrpath = wecc_terr, marpath = wecc_mar)
#print(ids)
#get_hads_station_csv(token = config.token, bucket_name = bucket_name, directory = directory, ids = ids.sample(40)) # .Sample() subset is for testing, remove for full run.
#get_hads_station_timeout_csv(token = config.token, bucket_name = bucket_name, directory = directory)
   

###
# Notes / in progress

# ISD 
# Run functions
#stations = get_hads_stations(wecc_terr, wecc_mar)
#print(stations) # For testing.
#get_asosawos_data_ftp(stations, bucket_name, directory, start_date = "2003-01-01", get_all = True)


# Function to get up to date station list of ASOS AWOS stations in WECC.
# Pulls in ISD station list and ASOSAWOS station list (two separate csvs), joins by ICAO and returns list of station IDs.
# Inputs: path to terrestrial WECC shapefile, path to marine WECC file. 
# Both paths given relative to home directory for git project.
#def get_hads_stations(terrpath, marpath): #Could alter script to have shapefile as input also, if there's a use for this.
# #    ## Login.
    
#     ## Get list of HADS stations as CSV
#     url = 'https://hads.ncep.noaa.gov/compressed_defs/all_dcp_defs.txt'
#     r = requests.get(url)
#     lines = r.content.split(b'\n') # Split by line
#     df = []
#     for line in lines:
#         row = line.split(b'|') # Split on |
#         row = row[0:10] # Only keep first 10 columns
#         row = [x.decode('utf-8') for x in row] # Convert to string
#         row = list(map(str.strip, row))
#         df.append(row)

#     # Create pandas dataframe
#     stations = pd.DataFrame(columns=['GOES NESDIS ID','NWSLI','DCP Owner','State Location','Hydrologic Service Area','Latitude','Longitude','Initial Daily Transmission Time (UTC)','DCP Transmission Interval (Minutes)','DCP Location Name'], data=df)
#     stations = stations.dropna()
    
#     # # Use spatial geometry to only keep points in wecc marine / terrestrial areas.
#     stations['Longitude'] = [calc_pull._lon_dms_to_dd(i) for i in stations['Longitude']]
#     stations['Latitude'] = [calc_pull._lat_dms_to_dd(i) for i in stations['Latitude']]

#     geometry = [Point(xy) for xy in zip(stations['Longitude'], stations['Latitude'])] # Zip lat lon coords.
#     weccgeo = gp.GeoDataFrame(stations, crs='EPSG:4326', geometry=geometry) # Convert to geodataframe.
    
#     ## get bbox of WECC to use to filter stations against
#     t, m, bbox = calc_pull.get_wecc_poly(terrpath, marpath) # Call get_wecc_poly.

#     # Get terrestrial stations.
#     weccgeo = weccgeo.to_crs(t.crs) # Convert to CRS of terrestrial stations.
#     terwecc = sjoin(weccgeo.dropna(), t, how='left') # Only keep stations in terrestrial WECC region.
#     terwecc = terwecc.dropna() # Drop empty rows.

#     # Get marine stations.
#     marwecc = sjoin(weccgeo.dropna(), m, how='left') # Only keep stations in marine WECC region.
#     marwecc = marwecc.dropna() # Drop empty rows.
    
    
#     #Join and remove duplicates using GOES NESDIS ID as combined unique identifier.
#     weccstations = (pd.concat([terwecc, marwecc], ignore_index=True, sort =False)
#            .drop_duplicates(['GOES NESDIS ID'], keep='first'))

#     # Join to list of ISD stations.
#     isdstations = pd.read_csv("test_platform/data/1_raw_wx/HADS/isd-history.csv")
#     isdstations = isdstations[(isdstations['CTRY']=="US")]

#     geometry = [Point(xy) for xy in zip(isdstations['LON'], isdstations['LAT'])] # Zip lat lon coords.
#     isdgeo = gp.GeoDataFrame(isdstations, crs='EPSG:4326', geometry=geometry) # Convert to geodataframe.
    
#     terwecc = sjoin(isdgeo.dropna(), t, how='left') # Only keep stations in terrestrial WECC region.
#     terwecc = terwecc.dropna() # Drop empty rows.

#     # Get marine stations.
#     marwecc = sjoin(isdgeo.dropna(), m, how='left') # Only keep stations in marine WECC region.
#     marwecc = marwecc.dropna() # Drop empty rows.
    
#     #Join and remove duplicates using GOES NESDIS ID as combined unique identifier.
#     isdstations = (pd.concat([terwecc, marwecc], ignore_index=True, sort =False)
#            .drop_duplicates(['USAF', 'WBAN'], keep='first')).reset_index()

# One approach: using sjoin (memory killed)
    # weccstations = weccstations.to_crs('EPSG:4269')
    # isdstations = isdstations.to_crs('EPSG:4269')

    # isdbuffer = isdstations.copy()
    # isdbuffer.geometry = isdstations.geometry.buffer(150)

    # weccstations = weccstations.drop(['index_right'], axis=1)
    # isdbuffer = isdbuffer.drop(['index_right'], axis=1)

    # join = gp.sjoin(weccstations, isdbuffer, how='left', op='within')
    # print(join)
    
    # weccstations.reset_index(drop = True, inplace = True)
    
    # # Get closest point
    # known_xy = np.stack([weccstations['Longitude'], weccstations['Latitude']], -1)
    # tree = scipy.spatial.cKDTree(known_xy)
   
    # query_xy = np.stack([isdstations['LON'], isdstations['LAT']], -1)
    # distances, indices = tree.query(query_xy)
    
    # df = pd.DataFrame([distances, indices]).transpose()
    # df.reset_index(inplace = True)
    
    # df.columns = ['ISD-Index', 'Distance', 'WECC-Index']
    # df['WECC-Index'] = df['WECC-Index'].astype(int)

    # df['WECC-Names'] = weccstations['DCP Location Name'][df['WECC-Index']].values
    # df['ISD-Names'] = isdstations.loc[df['ISD-Index'], 'STATION NAME'].values
    # df = df.sort_values(by ='Distance')
    # print(df.head(50))
    # # Not conclusive. Pause here.
