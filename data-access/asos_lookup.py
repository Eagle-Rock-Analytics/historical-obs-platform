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

    city_dict = {
        "Arcata": "ASOSAWOS_72594524283",
        "Eureka": "ASOSAWOS_72594524283",
        "Blue Canyon": "ASOSAWOS_72584523225",
        "Burbank": "ASOSAWOS_72288023152",
        "Fresno": "ASOSAWOS_72389093193",
        "Fullerton": "ASOSAWOS_72297603166",
        "Oakland": "ASOSAWOS_72493023230",
        "Palm Springs": "ASOSAWOS_72286893138",
        "Palmdale": "ASOSAWOS_72382023182",
        "Red Bluff": "ASOSAWOS_72591024216",
        "Riverside": "ASOSAWOS_72286903171",
        "Sacramento": "ASOSAWOS_72483023232",
        "San Diego": "ASOSAWOS_72290023188",
        "San Jose": "ASOSAWOS_72494523293",
        "Santa Maria": "ASOSAWOS_72394023273",
        "Santa Rosa": "ASOSAWOS_72495723213",
        "Torrance": "ASOSAWOS_72295503174",
    }

    code_dict = {
        "KACV": "ASOSAWOS_72594524283",
        "KBLU": "ASOSAWOS_72584523225",
        "KBUR": "ASOSAWOS_72288023152",
        "KFAT": "ASOSAWOS_72389093193",
        "KFUL": "ASOSAWOS_72297603166",
        "KOAK": "ASOSAWOS_72493023230",
        "KPSP": "ASOSAWOS_72286893138",
        "KPMD": "ASOSAWOS_72382023182",
        "KRBL": "ASOSAWOS_72591024216",
        "KRAL": "ASOSAWOS_72286903171",
        "KSAC": "ASOSAWOS_72483023232",
        "KSAN": "ASOSAWOS_72290023188",
        "KRHV": "ASOSAWOS_72494523293",
        "KSMX": "ASOSAWOS_72394023273",
        "KSTS": "ASOSAWOS_72495723213",
        "KTOA": "ASOSAWOS_72295503174",
    }

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
        hdp_station = [val for key, val in city_dict.items() if key in city]
        if hdp_station:
            # now pull the ID out from the list that is returned above
            hdp_station = hdp_station[0]
            print(
                f"The HDP station name for input airport city '{city}' is {hdp_station}"
            )
        else:
            print(f"Input city '{city}' not in station dictionary.")
    else:
        print(
            "Please input either an airport code (e.g. code='KSAC') or city (e.g. city='Sacramento')."
        )

    return None
