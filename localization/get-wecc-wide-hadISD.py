""" 
Need two CSVs (included in repo):
WECC-NOAA-stations.csv - has ICAO codes for each WECC weather station as
requested by stakeholders. 
US-NOAA-stations.csv - contains NOAA ISD metadata, including ICAO,
standard station name (important for consistency w/ their metadata),
and numeric station ID (WMO + FAA codes).

This code then matches the WECC stations via ICAO with their numeric
IDs, since HadISD file names for stations contain the IDs. After matching,
it downloads the data from the Met Office after checking that it exists,
then saves it locally.
"""

import pandas as pd
import requests
from urllib.error import HTTPError

# dataframe containing WECC station ICAOs
wecc_stns = pd.read_csv("WECC-NOAA-stations.csv", header=None)
wecc_stns.columns = ["state", "name", "icao"]
wecc_stns = wecc_stns.fillna("not found")
wecc_stns.drop(["state", "name"], axis=1, inplace=True)
wecc_stns = wecc_stns.astype(str)

# NOAA metadata dataframe
isd_stns = pd.read_csv("US-NOAA-stations.csv", header=None)
isd_stns.columns = [
    "wmo",
    "faa",
    "name",
    "country",
    "state",
    "icao",
    "latitude",
    "longitude",
    "elevation",
    "begin",
    "end",
]
isd_stns = isd_stns.fillna("not found")
isd_stns = isd_stns.astype(str)
isd_stns["faa"] = isd_stns["faa"].str.zfill(5)

# match ICAO from wecc_stns to isd_stns to get the ID number
commondf = pd.merge(wecc_stns, isd_stns, on=["icao"])
commondf.sort_values(
    by=["icao", "end"], ascending=[True, False], inplace=True, ignore_index=True
)
commondf["station id"] = commondf[["wmo", "faa"]].agg("-".join, axis=1)

# now check that the file exists, then download and save
# also save one csv that has WECC station metadata only
# this is used in later processing code
save_dir = "./wecc-hadisd3.3.0.202202p/"
url_prefix = "https://www.metoffice.gov.uk/hadobs/hadisd/v330_202202p/data/hadisd.3.3.0.202202p_19310101-20220301_"
file_prefix = "hadisd.3.3.0.202202p_19310101-20220301_"
wecclist = []  # will be populated with WECC metadata
for i in range(len(commondf)):
    to_find = commondf["station id"].values[i]
    icao_match = commondf["icao"].values[i]
    file_url = url_prefix + to_find + ".nc.gz"
    save_name = save_dir + file_prefix + to_find + "_" + icao_match + ".nc.gz"
    r = requests.get(file_url, timeout=1)
    if r.status_code != 200:
        continue
    else:
        wecclist.append([a for a in commondf.iloc[i].values])
        open(save_name, "wb").write(r.content)
    r = None

weccdf.to_csv("wecc-station-data.csv")  # save WECC metadata to file
