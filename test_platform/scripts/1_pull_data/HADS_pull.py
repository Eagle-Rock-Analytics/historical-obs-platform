"""
This script downloads HADS data directly from https://www.ncei.noaa.gov/data/nws-hydrometeorological-automated-data-system/archive/.
Approach:
(1) Get station list (does not need to be re-run constantly)
(2) Download data using station list.
Inputs: bucket name in AWS, directory to save file to (folder path), station list (optional), start date of file pull (optional),
parameter to only download changed files (optional)
Outputs: Raw data for an individual network, all variables, all times. Organized by time, with 1 file per day per year.

Notes:
1. This function assumes users have configured the AWS CLI such that their access key / secret key pair are stored in ~/.aws/credentials.
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
# Inputs: S3 path to terrestrial WECC shapefile, S3 path to marine WECC file.
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
    s3_cl.put_object(Bucket=bucket_name, Body=content,Key=directory+"stationlist_HADS.csv")
    return weccstations

# Function: Takes URL of given website and returns all file download links contained on it.
# Input: one url
# Output: zipped list of download urls and date last modified.

def get_file_links(url): 
      
    # create response object 
    r = requests.get(url) 
      
    # create beautiful-soup object 
    soup = BeautifulSoup(r.content,'html.parser') 
      
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

# Function: query NCEI server for HADS data and download zipped files.
# Run this to get updated data.
# Inputs:
# bucket_name: name of AWS bucket
# directory: folder path within bucket
# start_date: YYYY-MM-DD (optional)
# end_date: YYYY-MM-DD (optional)
def get_hads_update(bucket_name, directory, start_date = None, end_date = None):
    # Set up error handling
    errors = {'Date':[], 'Time':[], 'Error':[]}
    end_api = datetime.now().strftime('%Y%m%d%H%M') # Set end time to be current time at beginning of download

    # Set up AWS to write to bucket.
    s3 = boto3.client('s3')

    
    # Set up filtering by time.
    if start_date is None:
        if end_date is None:
            years = range(1997, int(datetime.now().year)+1)
        else:
            years = range(1997, int(end_date[0:4])+1)
            end_day = datetime.strptime(end_date, '%Y-%m-%d').strftime('%j')
    else:
        start_year = int(start_date[0:4])
        start_day = datetime.strptime(start_date, '%Y-%m-%d').strftime('%j')
        if end_date is None:
            years = range(start_year, int(datetime.now().year)+1)
        else:
            years = range(start_year, int(end_date[0:4])+1)
            end_day = datetime.strptime(end_date, '%Y-%m-%d').strftime('%j')

    
    for i in years:
        try:
            yearurl = 'https://www.ncei.noaa.gov/data/nws-hydrometeorological-automated-data-system/archive/{}/'.format(str(i))
            links = get_file_links(yearurl)

            days_to_download = list(range(0,367)) # To include leap years
            if start_date is not None and i == start_year:
                # Get rid of links before start date.
                days_to_download = [x for x in days_to_download if x >= int(start_day)]
                
            if end_date is not None and i == int(end_date[0:4]):
                # Get rid of links after end date.    
                days_to_download = [x for x in days_to_download if x <= int(end_day)]
            
            # TO DO: filter file names by day
            #for k, (link, date) in enumerate(links):
            #    if link.split("/")[-1] in alreadysaved:


            # # TO DO: Filter links by year of day

            #         print(links)
            #link_to_aws(links, bucket_name, directory)
        except Exception as e:
            print("Error in downloading year {}: {}". format(i, e))
            errors['Date'].append(i)
            errors['Time'].append(end_api)
            errors['Error'].append(e)

    # #Write errors to csv
    # csv_buffer = StringIO()
    # errors = pd.DataFrame(errors)
    # errors.to_csv(csv_buffer)
    # content = csv_buffer.getvalue()
    # s3.put_object(Bucket=bucket_name, Body=content,Key=directory+"errors_hads_{}.csv".format(end_api))

    return


if __name__ == "__main__":
    get_hads_stations(wecc_terr, wecc_mar)
    get_hads_dat(bucket_name, directory, start_date = None, get_all = True)

# Note, for first full data pull, set get_all = True
# For all subsequent data pulls/update with newer data, set get_all = False
