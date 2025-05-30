"""
This function generates and exports a dataframe with counts of unique QAQC flag values per variable - in their native timestep, before hourly standardization.
These tables are used to produce QAQC flag statistics for the QAQC success report.

"""

## Import Libraries
import pandas as pd

# Set relative paths to other folders and objects in repository.
bucket_name = "wecc-historical-wx"


# -----------------------------------------------------------------------------
def eraqc_counts_original_timestep(
    df: pd.DataFrame, network: str, station: str
) -> None:
    """
    Generates a dataframe of raw qaqc flag value counts for every variable,
    in their native timestep, before hourly standardization.
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

    # set all counts to integers, for readability
    flag_counts = flag_counts.astype(int)

    # rename columns
    flag_counts.columns = flag_counts.columns.str.replace("_eraqc", "", regex=True)

    # rename index (i.e. eraqc values) and then reset index
    flag_counts = flag_counts.rename_axis("eraqc_flag_values")

    # send file to AWS
    csv_s3_filepath = f"s3://wecc-historical-wx/4_merge_wx/{network}/eraqc_counts/{station}_flag_counts_native_timestep.csv"
    flag_counts.to_csv(csv_s3_filepath, index=True)

    return None

# -----------------------------------------------------------------------------
