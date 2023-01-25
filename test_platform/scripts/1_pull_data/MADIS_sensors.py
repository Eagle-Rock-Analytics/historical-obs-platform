"""
This script outputs a sensor list csv for each of the MADIS networks and is saved to the
qaqc_wx folder for each network.

This script only needs to be run occasionally/when a full MADIS pull occurs.
"""

# Step 0: Environment set-up
import requests
import pandas as pd
from datetime import datetime
import re
import boto3
from io import StringIO
import numpy as np
import calc_pull
from MADIS_pull import get_network_metadata

try:
    import config # Import API keys.
except:
    print("Missing config.py file with API token. Make file if necessary.")
    exit()

## Set AWS credentials
s3 = boto3.resource("s3")
s3_cl = boto3.client('s3') # for lower-level processes
bucket_name = "wecc-historical-wx"

# Set envr variables
# Set paths to WECC shapefiles in AWS bucket.
wecc_terr = "s3://wecc-historical-wx/0_maps/WECC_Informational_MarineCoastal_Boundary_land.shp"
wecc_mar = "s3://wecc-historical-wx/0_maps/WECC_Informational_MarineCoastal_Boundary_marine.shp"
directory = "3_qaqc_wx/"

# List of MADIS network names
madis_networks = ['CA HYDRO', 'CDEC', 'CNRFC', 'CRN', 'CWOP', 'HNXWFO', 'HOLFUY', 'HPWREN', 'LOXWFO', 'MAP', 'MTRWFO', 'NCAWOS', 'NOS-NWLON', 'NOS-PORTS',
'RAWS', 'SGXWFO', 'SHASAVAL', 'VCAPCD']

# Get dictionary of network names and short IDs
# Inputs: token, list of network names
def madis_network_name_to_id(token, networks):
    networks = get_network_metadata(token) # Get list of networks and IDs
    shortnames = [name for name in networks['SHORTNAME'] if any(x in name for x in madis_networks)]
    drop = ['MAP-ASN', 'NSRAWS', 'TMPRAWS']
    shortnames = [x for x in shortnames if x not in drop]
    networks = networks[networks['SHORTNAME'].isin(shortnames)]
    networks = networks[['SHORTNAME', 'ID']]
    networks['SHORTNAME'] = networks['SHORTNAME'].replace(to_replace="APRSWXNET/CWOP", value = "CWOP")
    networks['SHORTNAME'] = networks['SHORTNAME'].replace(to_replace="CA HYDRO", value="CAHYDRO") # removes space issue
    return networks


def get_madis_sensor_metadata(token, terrpath, marpath, networks, bucket_name, directory):
    try:
        t,m,bbox = calc_pull.get_wecc_poly(terrpath, marpath)
        bbox_api = bbox.loc[0,:].tolist() # [lonmin,latmin,lonmax,latmax]
        bbox_api = ','.join([str(elem) for elem in bbox_api])

        for index, row in networks.iterrows():
            networkname = row['SHORTNAME']
            networkid = row['ID']

            # Access station metadata to get list of IDs in bbox and network
            # Using: https://developers.synopticdata.com/mesonet/v2/stations/timeseries/
            url = "https://api.synopticdata.com/v2/stations/metadata?token={}&network={}&bbox={}&complete=1&sensorvars=1&recent=20&output=json".format(token, networkid, bbox_api)
            request = requests.get(url).json()

            ids = []
            station_list = pd.DataFrame(request['STATION'])
            station_list = pd.concat([station_list, station_list["PERIOD_OF_RECORD"].apply(pd.Series)], axis=1) # Split Period of Record column
            station_list = pd.concat([station_list, station_list["UNITS"].apply(pd.Series)], axis=1) # Split Period of Record column
            station_list = pd.concat([station_list, station_list["SENSOR_VARIABLES"].apply(pd.Series, dtype='string')], axis=1) # Split Variables column
            station_list = station_list.drop("PERIOD_OF_RECORD", axis =1)
            station_list = station_list.drop("UNITS", axis =1)
            station_list = station_list.drop("SENSOR_VARIABLES", axis =1)

            # Drop all columns with vars that are not of interest.
            coltokeep = ['altimeter', 'air_temp', 'relative_humidity',
                    'wind_speed', 'wind_direction', 'precip_accum_since_local_midnight',
                    'precip_accum_24_hour', 'dew_point_temperature', 'pressure',
                    'precip_accum_one_hour', 'solar_radiation', 'precip_accum', 'precip_accum_five_minute']

            # Get index of 'end', last standard column
            end_ind = station_list.columns.get_loc("end")

            # Drop columns
            station_list = station_list.drop(columns=[col for col in station_list.columns[end_ind+1:] if col not in coltokeep]) # Drop all columns not in coltokeep list.

            # Each sensor has a column of dictionaries. For each sensor, collapse column.
            cols_in_df = [i for i in coltokeep if i in station_list.columns]
            for i in cols_in_df:
                df = pd.json_normalize(station_list[i])
                station_list.pop(i)
                station_list = station_list.join(df)

            # Standardize NAs
            station_list.replace("None", np.nan, inplace = True)
            station_list = station_list.fillna(value=np.nan)

            # Standardize column names
            station_list.rename(columns=lambda s: s.replace("PERIOD_OF_RECORD.", ""), inplace=True)
            station_list.rename(columns=lambda s: s.replace(".", "_"), inplace=True)

            # Save station list to AWS
            savedir = directory+networkname+"/"

            csv_buffer_err = StringIO()
            station_list.to_csv(csv_buffer_err)
            content = csv_buffer_err.getvalue()

            print("Sensor list saved for {}.".format(networkname))
            s3_cl.put_object(Bucket=bucket_name, Body=content, Key=savedir+"sensorlist_{}.csv".format(networkname))

    except Exception as e:
        print("Error: {}".format(e))


if __name__ == "__main__":
    networks = madis_network_name_to_id(token = config.token, networks = madis_networks)
    get_madis_sensor_metadata(token = config.token, terrpath = wecc_terr, marpath= wecc_mar, networks = networks, bucket_name= bucket_name, directory=directory)
