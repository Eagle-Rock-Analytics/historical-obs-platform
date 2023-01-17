'''
This function iterates through all networks and checks to see what stations have been successfully cleaned, 
updating the station list in the 1_raw_wx folder to reflect station availability. 
If stations have been cleaned but do not have data available, they are added manually to the station list.

'''
import boto3
import pandas as pd
from io import BytesIO

# Set environment variables
bucket_name = "wecc-historical-wx"
raw_wx = "1_raw_wx/"
clean_wx = "2_clean_wx/"
s3 = boto3.resource("s3")  
s3_cl = boto3.client("s3") 


# Function: Given a network name, return a pandas dataframe containing the network's station list from the raw bucket.
def get_station_list(network):
    network = network.upper()
    network_prefix = raw_wx+network+"/"
    
    # If station list is CIMIS, extension is .xlsx
    if network == "CIMIS":
        station_list = f"stationlist_{network}.xlsx"
        obj = s3_cl.get_object(Bucket= bucket_name, Key= network_prefix+station_list)
        station_list = pd.read_excel(BytesIO(obj['Body'].read()))
    else:
        station_list = f"stationlist_{network}.csv"
        obj = s3_cl.get_object(Bucket= bucket_name, Key= network_prefix+station_list)
        station_list = pd.read_csv(obj['Body'])
    return station_list

# Function: Given a network name, return a pandas dataframe of all cleaned stations in the 2_clean_wx AWS bucket, with the date the file was last modified.
def get_cleaned_stations(network):
    df = {'ID':[], 'Time_Cleaned':[]}
    network = network.upper()
    network_prefix = clean_wx+network+"/"
    for item in s3.Bucket(bucket_name).objects.filter(Prefix = network_prefix+network+"_"): 
        clean_id = item.key.split("/")[-1].split("_")[-1].replace(".nc", "") # Get ID from file name
        time_mod = item.last_modified
        df['ID'].append(clean_id)
        df['Time_Cleaned'].append(time_mod)
    return pd.DataFrame(df)

# Function: Given a network name, return a pandas dataframe containing all errors reported for the network in the cleaning stage.
def parse_error_csv(network):
    errordf = []
    network = network.upper()
    errors_prefix = clean_wx+network+"/"+"errors"
    for item in s3.Bucket(bucket_name).objects.filter(Prefix = errors_prefix): 
        obj = s3_cl.get_object(Bucket= bucket_name, Key= item.key)
        errors = pd.read_csv(obj['Body'])
        if errors.empty:# If file empty
            continue
        else: 
            errors = errors[['File', 'Time', 'Error']]
            errordf.append(errors)
            
    errordf = pd.concat(errordf)
    errordf = errordf.drop_duplicates(subset = ['File', 'Error'])
    errordf = errordf[errordf.File != "Whole network"] # Drop any whole network errors
    return errordf

# Function: 
def clean_qa(network):
    network = network.upper()

    # Call functions
    stations = get_station_list(network)
    cleanids = get_cleaned_stations(network)
    errors = parse_error_csv(network)

    # Clean station list
    if 'Unnamed' in stations.columns:
        stations = stations.drop(["Unnamed: 0"], axis = 1) # Drop double index column

    # If station ID cleaned, add Y to cleaned column and time_cleaned 


    print(stations, cleanids, errors)

clean_qa('ndbc')