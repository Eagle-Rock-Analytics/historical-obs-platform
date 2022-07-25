## Historical Observations Platform
## Template Script for Stage: MERGE DATA

## TO DO LIST
## Any notes critical for further development, e.g.: AWS implementation

# Step 0: Environment set-up
# Import libraries -- delete as appropriate per datasource
import os
from datetime import datetime, timezone
import xarray as xr


# Set envr variables
workdir = "/path/to/working/directory/"
years = list(map(str,range(1980,datetim.enow().year+1))) # If needed


# Step 1: Read in qa/qc data
## Reading in all datasources here?
## Performance/efficiency checks

# Step 2: Merge data
## Likely needs extensive error checking + outputs to user
## Data formatting 
