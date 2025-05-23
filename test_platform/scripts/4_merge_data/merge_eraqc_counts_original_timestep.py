"""
This function generates and exports a dataframe with counts of unique QAQC flag values per variable.
These tables are used to produce QAQC flag statistics for the QAQC success report.

"""

## Import Libraries
from functools import reduce
import boto3
import pandas as pd
from io import BytesIO, StringIO

# New logger function
from merge_log_config import logger

# Set AWS credentials
s3 = boto3.resource("s3")
s3_cl = boto3.client("s3")  # for lower-level processes

# Set relative paths to other folders and objects in repository.
bucket_name = "wecc-historical-wx"


# -----------------------------------------------------------------------------
def eraqc_counts_original_timestep(
    df: pd.DataFrame, network: str, station: str
) -> None:
    """
    Generates a dataframe of raw qaqc flag value counts for every variable those are generated for.
    Exports the dataframe as a csv to AWS.

    Parameters
    ----------
    df : pd.DataFrame
        station dataset converted to dataframe through QAQC pipeline
    network: str
        network name
    station: str
        station name

    Returns
    -------
    None
    """
    # identify _eraqc variables
    eraqc_vars = [var for var in df.columns if "_eraqc" in var]

    # filter df for only qaqc columns
    # also replace Nan values with 'no_flag' for two reasons:
    #   1. to enable us to count total observations for the success report
    #   2. to clarify what the Nan value indicates
    df = df[eraqc_vars].fillna("no_flag")

    # generate df of counts of each unique flag for each variable
    # fill all Nan values with 0, since Nan = no observations counted
    flag_counts = df.apply(pd.Series.value_counts).fillna(0)

    # rename columns
    flag_counts.columns = flag_counts.columns.str.replace("_eraqc", "", regex=True)

    # rename index (i.e. eraqc values) and then reset index
    flag_counts = flag_counts.rename_axis("eraqc_flag_values").reset_index()

    # send file to AWS
    new_buffer = StringIO()
    flag_counts.to_csv(new_buffer, index=False)
    content = new_buffer.getvalue()
    key = f"4_merge_wx/{network}/eraqc_counts/original_timestep_{station}.csv"

    s3_cl.put_object(
        Bucket=bucket_name,
        Body=content,
        Key=key,
    )

    return None
