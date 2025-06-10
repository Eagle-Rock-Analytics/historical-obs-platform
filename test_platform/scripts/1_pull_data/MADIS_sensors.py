"""
MADIS_sensors.py

This script outputs a sensor list csv for each of the MADIS networks and is saved to the qaqc_wx folder for each network.

Functions
---------
- madis_network_name_to_id: Get dictionary of network names and short IDs.
- get_madis_sensor_metadata: Retrieves sensor metadata for each network. 

Intended Use
------------
This script only needs to be run occasionally/when a full MADIS pull occurs.
"""

import requests
import pandas as pd
from datetime import datetime
import re
import boto3
from io import StringIO
import numpy as np

try:
    from calc_pull import get_wecc_poly
except:
    print("Error importing calc_pull")

try : 
    from MADIS_pull import get_network_metadata
except: 
    print("Error importing get_network_metadata")

try:
    import config  # Import API keys.
except:
    print("Missing config.py file with API token. Make file if necessary.")
    exit()

s3 = boto3.resource("s3")
s3_cl = boto3.client("s3")  # for lower-level processes
BUCKET_NAME = "wecc-historical-wx"
WECC_TERR = (
    "s3://wecc-historical-wx/0_maps/WECC_Informational_MarineCoastal_Boundary_land.shp"
)
WECC_MAR = "s3://wecc-historical-wx/0_maps/WECC_Informational_MarineCoastal_Boundary_marine.shp"
DIRECTORY = "3_qaqc_wx/"

# List of MADIS network names
MADIS_NETWORKS = [
    "CA HYDRO",
    "CDEC",
    "CNRFC",
    "CRN",
    "CWOP",
    "HNXWFO",
    "HOLFUY",
    "HPWREN",
    "LOXWFO",
    "MAP",
    "MTRWFO",
    "NCAWOS",
    "NOS-NWLON",
    "NOS-PORTS",
    "RAWS",
    "SGXWFO",
    "SHASAVAL",
    "VCAPCD",
]


def madis_network_name_to_id(token: str, networks: list[str]):
    """Get dictionary of network names and short IDs.

    Parameters
    ----------
    token : str
        API token
    networks : list[str]
        network names to retrieve

    Returns
    -------
    networks : dict
        network names to retrieve
    """

    # Get list of networks and IDs
    networks = get_network_metadata(token)  
    shortnames = [
        name for name in networks["SHORTNAME"] if any(x in name for x in madis_networks)
    ]
    drop = ["MAP-ASN", "NSRAWS", "TMPRAWS"]
    shortnames = [x for x in shortnames if x not in drop]
    networks = networks[networks["SHORTNAME"].isin(shortnames)]
    networks = networks[["SHORTNAME", "ID"]]
    networks["SHORTNAME"] = networks["SHORTNAME"].replace(
        to_replace="APRSWXNET/CWOP", value="CWOP"
    )
    
    # space issue between "CA HYDRO"
    networks["SHORTNAME"] = networks["SHORTNAME"].replace(
        to_replace="CA HYDRO", value="CAHYDRO"
    )

    return networks


def get_madis_sensor_metadata(
    token: str, terrpath: str, marpath: str, networks: dict, directory: str
):
    """Retrieves sensor metadata for each network. 
    
    Parameters
    ----------
    token : str
        API token
    terrpath : str

    marpath : str

    networks : dict
        network names to retrieve metadata for
    directory : str
        AWS path to directory
    
    Returns
    -------
    None
    """

    try:
        t, m, bbox = get_wecc_poly(terrpath, marpath)
        bbox_api = bbox.loc[0, :].tolist()  # [lonmin,latmin,lonmax,latmax]
        bbox_api = ",".join([str(elem) for elem in bbox_api])

        for index, row in networks.iterrows():
            networkname = row["SHORTNAME"]
            networkid = row["ID"]

            # Access station metadata to get list of IDs in bbox and network
            # Using: https://developers.synopticdata.com/mesonet/v2/stations/timeseries/
            url = f"https://api.synopticdata.com/v2/stations/metadata?token={token}&network={networkid}&bbox={bbox_api}&complete=1&sensorvars=1&recent=20&output=json"
            request = requests.get(url).json()

            ids = []
            station_list = pd.DataFrame(request["STATION"])
            station_list = pd.concat(
                [station_list, station_list["PERIOD_OF_RECORD"].apply(pd.Series)],
                axis=1,
            ) 

            # Split Period of Record column
            station_list = pd.concat(
                [station_list, station_list["UNITS"].apply(pd.Series)], axis=1
            )  
            station_list = pd.concat(
                [
                    station_list,
                    station_list["SENSOR_VARIABLES"].apply(pd.Series, dtype="string"),
                ],
                axis=1,
            ) 
            station_list = station_list.drop("PERIOD_OF_RECORD", axis=1)
            station_list = station_list.drop("UNITS", axis=1)
            station_list = station_list.drop("SENSOR_VARIABLES", axis=1)

            # Drop all columns with vars that are not of interest
            coltokeep = [
                "altimeter",
                "air_temp",
                "relative_humidity",
                "wind_speed",
                "wind_direction",
                "precip_accum_since_local_midnight",
                "precip_accum_24_hour",
                "dew_point_temperature",
                "pressure",
                "precip_accum_one_hour",
                "solar_radiation",
                "precip_accum",
                "precip_accum_five_minute",
            ]

            # Get index of 'end', last standard column
            end_ind = station_list.columns.get_loc("end")

            # Drop columns
            station_list = station_list.drop(
                columns=[
                    col
                    for col in station_list.columns[end_ind + 1 :]
                    if col not in coltokeep
                ]
            ) 

            # Each sensor has a column of dictionaries. For each sensor, collapse column
            cols_in_df = [i for i in coltokeep if i in station_list.columns]
            for i in cols_in_df:
                df = pd.json_normalize(station_list[i])
                station_list.pop(i)
                station_list = station_list.join(df)

            # Standardize NAs
            station_list.replace("None", np.nan, inplace=True)
            station_list = station_list.fillna(value=np.nan)

            # Standardize column names
            station_list.rename(
                columns=lambda s: s.replace("PERIOD_OF_RECORD.", ""), inplace=True
            )
            station_list.rename(columns=lambda s: s.replace(".", "_"), inplace=True)

            # Save station list to AWS
            savedir = directory + networkname + "/"

            csv_buffer_err = StringIO()
            station_list.to_csv(csv_buffer_err)
            content = csv_buffer_err.getvalue()

            print(f"Sensor list saved for {networkname}.")
            s3_cl.put_object(
                Bucket=BUCKET_NAME,
                Body=content,
                Key=savedir + f"sensorlist_{networkname}.csv",
            )

    except Exception as e:
        print(f"Error: {e}")

    return None


if __name__ == "__main__":
    networks = madis_network_name_to_id(config.token, MADIS_NETWORKS)
    get_madis_sensor_metadata(config.token, WECC_TERR, WECC_MAR, networks, DIRECTORY)
