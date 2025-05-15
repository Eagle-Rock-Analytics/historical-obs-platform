"""
This is a script where Stage 3: QA/QC related common functions, conversions, and operations is stored for ease of use
for the Historical Observations Platform.
"""

## Import Libraries
import boto3
import numpy as np
import pandas as pd

# New logger function
from merge_log_config import logger

## Set AWS credentials
s3 = boto3.resource("s3")
s3_cl = boto3.client("s3")  # for lower-level processes

## Set relative paths to other folders and objects in repository.
BUCKET_NAME = "wecc-historical-wx"


# ----------------------------------------------------------------------
def delete_vars(df: pd.DataFrame, var_attrs: dict) ->tuple[pd.DataFrame, dict]:




    return df, var_attrs 