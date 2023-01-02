# From https://www.ncei.noaa.gov/access/homr/api
# Takes any of the following ids:
#qid=[COOP|FAA|GHCND|GHCNMLT|ICAO|IGRA|NCDCSTNID|NWSLI|TRANS|WBAN|WMO]?:[a-z0-9]*

# Returns four tables

import requests
import pandas as pd
from io import StringIO
import boto3

## Set AWS credentials
s3_cl = boto3.client('s3') 
bucket_name = "wecc-historical-wx"


def flatten_data(y):
    out = {}

    def flatten(x, name=''):
        if type(x) is dict:
            for a in x:
                flatten(x[a], name + a + '_')
        elif type(x) is list:
            i = 0
            for a in x:
                flatten(a, name + str(i) + '_')
                i += 1
        else:
            out[name[:-1]] = x

    flatten(y)
    return out


def get_homr_metadata(id):
    testurl = f'https://www.ncei.noaa.gov/access/homr/services/station/{id}'
    request = requests.get(testurl).json()
    for i in request['stationCollection']['stations']:
        #print(i)
        #flat = flatten_data(i) # Ignore dictionaries, keep each row as a station
        pandas = pd.json_normalize(i, max_level = 0)

        # Split into 4 tables:
        # Names
        names = pd.json_normalize(pandas["names"])
        print(names)
        exit()
            
        # Identifiers & Platforms
        # Location
        # Remarks
        # Updates

        print(pandas)
        print(pandas.columns)
        exit()
    #print(homr_metadata)
    #print(homr_metadata.columns)

# get_homr_metadata("coop", "0467")


def get_all_homr_ids():
    # Function to iterate through all WECC states and save header HOMR metadata for all stations.
    # Saves homr_ids.csv to the QAQC folder in AWS.

    dfs = []
    states = ['WA', 'OR', 'CA', 'ID', 'NV', 'AZ', 'UT', 'NM', 'CO', 'TX', 'WY', 'MT', 'SD']
    for state in states:
        testurl = f'https://www.ncei.noaa.gov/access/homr/services/station/search?headersOnly=true&state={state}'
        try:
            request = requests.get(testurl).json()
            df = pd.json_normalize(request['stationCollection']['stations'])
            df.columns = [col.replace("header.", "") for col in df.columns]
            df.columns = [col.replace("por.", "") for col in df.columns]
            dfs.append(df)
            #print(df)
        except Exception as e:
            print(f"Unable to process {state}. {e}")
            continue

    homr_df = pd.concat(dfs)
    # Save station list to AWS
    savedir = "3_qaqc_wx/"
    
    csv_buffer = StringIO()
    homr_df.to_csv(csv_buffer, index = False)
    content = csv_buffer.getvalue()
    
    s3_cl.put_object(Bucket="wecc-historical-wx", Body=content, Key=savedir+"homr_ids.csv")

if __name__ == "__main__":   
    #get_all_homr_ids() # Only run periodically.
    get_homr_metadata('20002078')
    