"""
asos_lookup.py

A user can pass either the ASOS airport code or the “city” and the script/function identifes and returns the name of the HDP station.

This dictionaries developed for this function are specific to a set of stations of interest provided by our partners.

Usage
-----
The user runs the function using an airport code: asos_station_lookup(code='KSAC')
Or using the airport city: asos_station_lookup(city='Sacramento')

The function then identifies this input information as matching with HDP station ID “ASOSAWOS_72483023232”

The user can also input the fill airport name:asos_station_lookup(city='Sacramento Executive Airport')
In which case the function takes the city indicated in that name and returns the matching HDP station ID.

Output
------
The HDP station ID is returned.

"""
import pandas as pd
import boto3

# Set AWS credentials
s3 = boto3.resource("s3")
s3_cl = boto3.client("s3")  # for lower-level processes

# Set relative paths to other folders and objects in repository.
BUCKET_NAME = "wecc-historical-wx"
QAQC_DIR = "3_qaqc_wx"
MERGE_DIR = "4_merge_wx"
stations_csv_path = f"s3://{BUCKET_NAME}/{MERGE_DIR}/all_network_stationlist_merge.csv"


def asos_station_lookup(code: str | None = None, city: str | None = None) -> None:
    """
    Function that returns the HDP station ID for an input of either the
    1. four digit code or (e.g. KSAC)
    2. city (e.g. Sacramento)

    of the airport associated with the station.

    This functon is specific to the 16 stations of interest shared by project partners.

    Paramters
    ---------
    code : str
        ASOS airport code (e.g. KSAC weather station)
    city : str
        ASOS aiport city (e.g. Sacramento Executive Airport weather station)

    Returns
    -------
    None

    """
    # Define dictionaries matching HDP station IDs to airport codes and cities
    ## this was developed beforehand matching HDP ASOSAWOS station locations with those in a list of stations provided by partners

    merge_list = pd.read_csv(
        f"s3://wecc-historical-wx/{MERGE_DIR}/ASOSAWOS/stationlist_ASOSAWOS_merge.csv"
    )
    code_dict = pd.Series(
        merge_list["ERA-ID"].values, index=merge_list["ICAO"]
    ).to_dict()
    city_dict = pd.Series(
        merge_list["ERA-ID"].values, index=merge_list["STATION NAME"]
    ).to_dict()

    # If user inputs an airport code
    if code:
        if code in code_dict:
            hdp_station = code_dict[code]
            print(
                f"The HDP station name for input airport code {code} is {hdp_station}"
            )
        else:
            print(f"Input code '{code}' not in station dictionary.")
    # If user inputs a city
    elif city:
        # this catches cases in which the user inputs the entire airport name
        hdp_station = [value for key, value in city_dict.items() if city.upper() in key]
        if len(hdp_station) == 1:
            # now pull the ID out from the list that is returned above
            hdp_station = hdp_station[0]
            print(
                f"The HDP station name for input airport city '{city}' is {hdp_station}"
            )
        elif len(hdp_station) > 1:
            print(
                f"There are multiple stations associated with '{city}': {hdp_station}"
            )
        else:
            print(f"Input city '{city}' not in station dictionary.")
    else:
        print(
            "Please input either an airport code (e.g. code='KSAC') or city (e.g. city='Sacramento')."
        )

    return None
