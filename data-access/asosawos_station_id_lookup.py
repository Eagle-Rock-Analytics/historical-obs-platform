"""
asosawos_station_id_lookup.py

A user can pass either the ASOS airport code or the “city” and the script/function identifes and returns the name of the HDP station.

Usage
-----
Using the function asos_station_lookup() outside this script:
    The user runs the function using an airport code: asos_station_lookup(code='KSAC')
    Or using the airport city: asos_station_lookup(city='Sacramento')

    The function then identifies this input information as matching with HDP station ID “ASOSAWOS_72483023232”

    The user can also input the fill airport name:asos_station_lookup(city='Sacramento Executive Airport')
    In which case the function takes the city indicated in that name and returns the matching HDP station ID.

The script asosawos_station_id_lookup.py:
    Run "python asosawos_station_id_lookup.py"
    The script asks the user if they want to look up a station based on the airport, city, or airport code.
    The user must input "city", "aiport", or "code"
    The script then asks the user to input the airport, city, or airpot code.
    It returns the ID(s) of the station(s) associated with that input
    Example:
        python asos_lookup.py
        Would you like to input a city, airport name, or code? (type 'city','airport', or 'code'): city
        Please type the city or airport name: sacramento
        There are multiple stations associated with 'sacramento': ['ASOSAWOS_72483023232', 'ASOSAWOS_72483323206']

"""

import pandas as pd

# Set relative paths to other folders and objects in repository.
BUCKET_NAME = "wecc-historical-wx"
MERGE_DIR = "4_merge_wx"
asosawos_station_list_path = (
    f"s3://{BUCKET_NAME}/{MERGE_DIR}/ASOSAWOS/stationlist_ASOSAWOS_merge.csv"
)

def asosawos_station_lookup(code: str | None = None, city: str | None = None) -> None:
    """
    Function that returns the HDP station ID for an input of either the
    1. four digit code or (e.g. KSAC)
    2. city (e.g. Sacramento) or airport name (e.g. Boise Airport)

    of the airport associated with the station.

    Paramters
    ---------
    code : str
        ASOS airport code (e.g. KSAC weather station)
    city : str
        ASOS aiport city (e.g. Sacramento Executive Airport weather station)

    Returns
    -------
    None

    Example:
    -------
    asos_station_lookup(code=None, city='Sacramento')
    """
    # Define dictionaries matching HDP station IDs to airport codes and cities
    merge_list = pd.read_csv(asosawos_station_list_path)
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
    # If user inputs a city or airport
    elif city:
        # this allows for the user to input the city or entire airport name, and is also case-insensitive
        hdp_station = [value for key, value in city_dict.items() if city.upper() in key]
        if (len(hdp_station) == 1):  
            # if there is only one station associated with the input
            # now pull the ID out from the list that is returned above
            hdp_station = hdp_station[0]
            print(
                f"The HDP station name for input airport city '{city}' is {hdp_station}"
            )
        elif (len(hdp_station) > 1):  
            # if there are multuple stations associated with the input
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


if __name__ == "__main__":

    city_or_code = input(
        f"Would you like to input a city, airport name, or code? (type 'city', 'airport', or 'code'): "
    )

    if city_or_code == "code":
        input_code = input(f"Please type the four-letter airport code (ex: KSAC): ")
        asosawos_station_lookup(code=input_code)
    elif city_or_code == "city" or city_or_code == "airport":
        input_city = input(f"Please type the city or airport name: ")
        asosawos_station_lookup(city=input_city)
