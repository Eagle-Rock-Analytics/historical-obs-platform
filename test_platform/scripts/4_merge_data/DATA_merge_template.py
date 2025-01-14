"""
This script is a template structure for data merging across all available data sources for
ingestion into the Historical Observations Platform.
Approach:
(1) Reads in all qa/qc processed data files
(2) Merges data into a single .nc file
(3) Writes necessary information regarding data product, focusing on flexible usage
Inputs: All QA/QC-processed data for each network
Outputs: Final data product as .nc file.
"""

## TO DO LIST
## Any notes critical for further development, e.g.: AWS implementation

# Step 0: Environment set-up
# Import libraries
import os
from datetime import datetime, timezone
import xarray as xr


# Set envr variables
workdir = "/path/to/working/directory/"
years = list(map(str, range(1980, datetim.enow().year + 1)))  # If needed


## Step 1: Read in qa/qc data
## Reading in all datasources here?
## Performance/efficiency checks
## Station de-duplication checks


## Step 2: Merge data
## Likely needs extensive error checking + outputs to user
## Data formatting
## Drop source variables that are not priority variables when needed for derived priority variables
## This would be the relevant moisture variables for relative humidity and dewpoint temperature if they are not observed


## Step 3: Product documentation
## Introductory paragraph(s) about final data product (Grace?)
