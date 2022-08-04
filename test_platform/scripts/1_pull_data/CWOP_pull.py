# CWOP pull script
# This script contains 3 functions. 
# get_wecc_poly generates a bounding box from WECC shapefiles.
# get_meso_metadata produces a list of station IDs filtered by bounding box and network (CWOP in this case)
# get cwop_station_csv saves a csv for each station ID provided, with the option to select a start date for downloads (defaults to station start date).
## It also produces an error csv noting all station IDs where downloads failed.

import requests
import geopandas as gp
import pandas as pd
import numpy as np
from datetime import datetime
import os
import re
import csv

# Set envr variables
token = "34e62da9b7c74f15a02ca172e6206bc3" # Synoptic API token - should be moved to AWS secrets manager or .env file asap.
wecc_terr = 'test_platform/data/0_maps/WECC_Informational_MarineCoastal_Boundary_land.shp'
wecc_mar = 'test_platform/data/0_maps/WECC_Informational_MarineCoastal_Boundary_marine.shp' 
save_dir = 'test_platform/data/1_raw_wx/CWOP/'   

# Function to return wecc shapefiles and combined bounding box given path variables.
def get_wecc_poly(terrpath, marpath):
    ## get bbox of WECC to use to filter stations against
    ## Read in terrestrial WECC shapefile.
    t = gp.read_file(terrpath)
    ## Read in marine WECC shapefile.
    m = gp.read_file(marpath)
    ## Combine polygons and get bounding box of union.
    bbox = t.union(m).bounds
    return t,m, bbox

def get_meso_metadata(terrpath, marpath):
    try:
        t,m,bbox = get_wecc_poly(terrpath, marpath)
        bbox_api = bbox.loc[0,:].tolist() # [lonmin,latmin,lonmax,latmax]
        bbox_api = ','.join([str(elem) for elem in bbox_api])
        # Access station metadata to get list of IDs in bbox and network
        # Using: https://developers.synopticdata.com/mesonet/v2/stations/timeseries/
        # CWOP data network ID in synoptic is 65 (see below)
        url = "https://api.synopticdata.com/v2/stations/metadata?token={}&network=65&bbox={}&recent=20&output=json".format(token, bbox_api)
        request = requests.get(url).json()
        #print(request)
        ids = []
        for each in request['STATION']:
            ids.append([each['STID'], each['PERIOD_OF_RECORD']['start']]) # Keep station ID and start date
        ids = pd.DataFrame(ids, columns = ['STID', 'start']).sort_values('start') # Sort by start date (note some stations return 'None' here)
        # Reformat date to match API format
        ids['start'] = pd.to_datetime(ids['start'], format='%Y-%m-%dT%H:%M:%SZ')
        ids['start'] = ids['start'].dt.strftime('%Y%m%d%H%M')
        return ids
    except Exception as e:
        print("Error: {}".format(e))

# Make a csv for each station and save. 
# To do:
# Add error handling.
def get_cwop_station_csv(token, ids, save_dir, start_date = None): 
# Takes dataframe with two columns, station ID and startdate. Optional parameter to specify later start date ("YYYYMMDDHHMM")
    
    try:
        os.mkdir(save_dir) # Make the directory to save data in. Except used to pass through code if folder already exists.
    except:
        pass
    
    # Set up error handling df.
    errors = {'Station ID':[], 'Time':[], 'Error':[]}

    # Set end time to be current time at beginning of download
    end_api = datetime.now().strftime('%Y%m%d%H%M')
        
    for index, id in ids.iterrows(): # For a group of stations
        
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
        url = "https://api.synopticdata.com/v2/stations/timeseries?token={}&stid={}&start={}&end={}&output=csv&qc=on&qc_remove_data=off&qc_flags=on".format(token, id['STID'], start_api, end_api)
        #print(url) # For testing.

        # Try to get station csv.
        try:
            #request = requests.get(url)
            filepath = save_dir+'{}.csv'.format(id["STID"]) # Set path to desired folder. # Write file to name of station ID in synoptic-- Change file name to reflect dates?? 
            #print(filepath)
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
                        with open(filepath, 'wb') as f:
                            print("Saving data for station {}".format(id["STID"])) # Nice for testing, remove for full run.
                            for line in r.iter_lines():
                                f.write(line+'\n'.encode())
                                next
                else:
                    errors['Station ID'].append(id['STID'])
                    errors['Time'].append(end_api)
                    errors['Error'].append(r.status_code)
                    print("Error: {}".format(r.status_code))

                    # If needed, split into smaller chunk? But csv seems to not kick so many errors.
            # if request['SUMMARY']['RESPONSE_CODE'] == -1 & request['SUMMARY']['RESPONSE_MESSAGE'].startswith("Querying too many station hours"):
            # # Response when the API call is too large. Split into smaller chunks.
            #     print("API request too large. Splitting into smaller chunks.")
            #     chunks = np.array_split(chunk, 3)
            #     pass # Write code to rerun w/ smaller chunks here.
            
            # Parse response into dataframe
            #### TO DO: clean and parse at the same time here?

        except Exception as e:
            print("Error: {}".format(e))
    
    # Write errors to csv
    filepath = save_dir+"errors_cwop_{}.csv".format(end_api) # Set path to save error file.
    print(errors)
    with open(filepath, "w") as outfile:
        # pass the csv file to csv.writer function.
        writer = csv.writer(outfile)
    
        # pass the dictionary keys to writerow
        # function to frame the columns of the csv file
        writer.writerow(errors.keys())
    
        # make use of writerows function to append
        # the remaining values to the corresponding
        # columns using zip function.
        writer.writerows(zip(*errors.values()))

# Test!
ids = get_meso_metadata(wecc_terr, wecc_mar)
# Get 3 real rows (or more as desired.)
ids = ids.sample(3)
# And make 3 fake rows that might break our code - wrong station id, wrong time, and NaN in time format.
test = pd.DataFrame({'STID': ["0ier", "E2082", "F2382"],
                     'start': ["200001010000", "2", "NaN"]})
ids = ids.append(test, ignore_index = True)
print(ids)
get_cwop_station_csv(token = token, ids = ids, save_dir = save_dir) # Run this with our test data.


### NOTES/SCRAPS FROM ALONG THE WAY

# with open("test.csv", 'wb') as f, \
#                 requests.get("https://api.synopticdata.com/v2/stations/timeseries?token=34e62da9b7c74f15a02ca172e6206bc3&stid=irgt&start=201905132239&end=202207261229&output=csv&qc=on&qc_remove_data=off&qc_flags=on", stream=True) as r: 
#                 id = "AFtest"
#                 end_api = "200000000000"
#                 if r.status_code == 200: # If API call returns a response
#                     #print("Saving data for station {}".format(id["STID"])) # Nice for testing, remove for full run.
#                     #print(r.text())
#                     if "# RESPONSE_MESSAGE" in r.text:
#                        error_text = str(re.search("(RESPONSE_MESSAGE:.*)",r.text).group(0))
#                        error_text = error_text.replace("RESPONSE_MESSAGE:", "")
#                        error = {'Station ID': id, 'Time': end_api, 'Error': error_text} # Add error to error log. 
#                        print(error)
#                     #for line in r.iter_lines():
                        
#                         #f.write(line+'\n'.encode())
#                         #next




# Option 2: generate JSONS for sets of stations. Still to figure out - how to join etc.
### Not finished!
start_date = None
def get_cwop_stations(token, ids, chunk_size, start_date = None): 
# Takes dataframe with two columns, station ID and startdate. Optional parameter to specify later start date ("YYYYMMDDHHMM")
# no_chunk specifies the number of groups the ID list is split into. Play with this if API starts to break.
    no_chunk = int(len(ids)/chunk_size) # Calculate the number of groups needed to have "chunk_size" stations in each grouo.
    id_chunks = np.array_split(ids, no_chunk)
    for chunk in id_chunks: # For a group of stations
        # Get list of station IDs
        idlist = ",".join(chunk["STID"])

        # Get first start date in chunk
        if start_date is None:
            if chunk["start"].isnull().any()==True: # Adding error handling for NaN dates.
                start_api = "198001010000" # If any of the stations in the chunk have a start time of NA, pull data from 01-01-1980 and on.
            else:
                start_api = chunk["start"].min()
        else:
            start_api = start_date

        end_api = datetime.now().strftime('%Y%m%d%H%M')
        url = "https://api.synopticdata.com/v2/stations/timeseries?token={}&stid={}&start={}&end={}&output=json&qc_remove_data=off&qc_flags=on".format(token, idlist, start_api, end_api)
        print(url) # For testing.
        try:
            request = requests.get(url).json()
            if request.status_code == 200:
                pass
                # WRITE CODE TO SAVE / APPEND.
            if request['SUMMARY']['RESPONSE_CODE'] == -1 & request['SUMMARY']['RESPONSE_MESSAGE'].startswith("Querying too many station hours"):
            # Response when the API call is too large. Split into smaller chunks.
                print("API request too large. Splitting into smaller chunks.")
                chunks = np.array_split(chunk, 3)
                pass # Write code to rerun w/ smaller chunks here.
            
            # Parse response into dataframe
            #### TO DO: clean and parse at the same time here?

        except Exception as e:
            print("Error: {}".format(e))
        # get_cwop_stations(token = token, ids = ids, chunk_size = 5, start_date = "202001010000")


# t,m,bbox = get_wecc_poly(wecc_terr, wecc_mar)
# bbox_api = bbox.loc[0,:].tolist() # [lonmin,latmin,lonmax,latmax]
# bbox_api = ','.join([str(elem) for elem in bbox_api])
        
# url = "https://api.synopticdata.com/v2/stations/metadata?token={}&network=65&bbox={}&stid=F1632&recent=20&output=json".format(token, bbox_api)
# request = requests.get(url).json()
# print(request)

# For testing, grab last 20 minutes of data.
#url = "https://api.synopticdata.com/v2/stations/timeseries?token={}&network=65&bbox={}&recent=20&output=json&qc_remove_data=off&qc_flags=on".format(token, bbox_api)
#request = requests.get(url).json()
#print(request)

#print(bbox_api) For testing

# Get info about CWOP network
#url = "https://api.synopticdata.com/v2/networks?token={}&shortname=APRSWXNET/CWOP".format(token)
#request = requests.get(url).json()
#print(request) # Print information about CWOP network.


# From https://github.com/derksmertzer/madis-dumper/blob/master/madis_class.py
def get_files(self):

        """ CALLED by date_construct: get_files() """

        day_dl = time()
        print('starting download for next day')

        # loop through hours
        for k in range(self.hours):
            dates = next(self.iter)
            http = ARCHIVE_URL + dates[0] + '/' + dates[1] + '/' + dates[2] + URL_MADIS + dates[0] + dates[1] + dates[2] + \
                dates[3] + '.gz'
            gz_file = http.split('/')[-1]
            file = gz_file.split('.')[0]
            hr_dl = time()
            print('downloading file: {}'.format(file))

            # attempt to connect to madis https server
            bck_off = 1
            with open(DL_DIR + gz_file, 'wb') as f:
                for attempt in range(1, 11):
                    try:
                        r = requests.get(http)
                        f.write(r.content)
                        str_error = None
                    except requests.exceptions.ConnectionError:
                        print('...error connecting to host')
                        str_error = True
                    if str_error:
                        sleep(bck_off)
                        bck_off *= 2
                        self.diagnostics(file, attempt, 0)
                    elif attempt == 10:
                        self.diagnostics(file, 1)
                        raise SystemExit('...unable to connect to host, terminating program at file: {}'.format(file))
                    else:
                        break

                print('...connection to host successful')

            try:
                with gzip.open(DL_DIR + gz_file, 'rb') as gz:
                    with open(DL_DIR + file, 'wb') as nc:
                        nc.write(gz.read())
            # catches missing files and continues to next iteration
            except OSError:
                print("...file {0} is not located on https server".format(file))
                self.diagnostics(file)
                os.remove(DL_DIR + gz_file)
                os.remove(DL_DIR + file)
                continue
            # catches corrupt (or possibly no contents to read in?) files and continues to next iteration
            except EOFError:
                print("...file {0} either corrupt or too small".format(file))
                self.diagnostics(file)
                os.remove(DL_DIR + gz_file)
                os.remove(DL_DIR + file)
                continue

            print('...time to download file: {:.2f}'.format(time() - hr_dl))
            os.remove(DL_DIR + gz_file)
