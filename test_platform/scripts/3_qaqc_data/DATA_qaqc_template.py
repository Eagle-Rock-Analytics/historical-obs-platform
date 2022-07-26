"""
This script is a template structure for data qa/qc protocols for a variety of data sources for
ingestion into the Historical Observations Platform.
Approach:
(1) Remove duplicate stations
(2) Handle variables that report at different intervals and/or change frequency over time (convert to hourly?)
(3) QA/QC testing, including consistency checks, gaps, checks against climatological distributions, and cross variable checks.
(4) Case study analysis for accuracy
Inputs: Cleaned data for an individual network
Outputs: QA/QC-processed data for an individual network, priority variables, all times. Organized by station as .nc file.
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
years = list(map(str,range(1980,datetim.enow().year+1))) # If needed


## Step 1: Remove duplicate stations
# Alternative spellings of station ids (spelling case)
# lat-longitude
# short distance (40-50km as baseline?)


## Step 2: Handle variables that report at different intervals and/or change frequency over time
## Convert to hourly data?


## Step 3: QA/QC (these could potentially be individual steps too)
## Testing with random subsample to ensure acceptably low level of false positive rates for specified thresholds and tests for each variables

# How to handle existing QA/QC flags -- use tracked qa/qc flags from DATA_clean.py

# Internal variable consistency checks
    # logical checks for climate variables
    # Examples:
        # Precipitation should not be negative
        # Relative humidity values between 0 and 100 (or 0 and 1)
        # Dew point temperature does not exceed air temperature
    # Gaps for jumps in a value
        # Use delta distributions to test
    # Climatological checks against distributions
        # WMO: 0.3 - 99.7 as percentiles
    # Repeated values

# Cross-variable checks
    # Logic checks
        # Example: wind speed and direction
        # If speed is 0, direction should be 0

# Drop original data variables to only keep qa/qc processed converted variables

## Step 3: Case study analysis for quality check
## Is this a separate script once a complete product at this stage or after merging?
