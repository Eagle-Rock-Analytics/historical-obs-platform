## ## Historical Observations Platform
## Template Script for Stage: QA/QC DATA

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


# Step 1: Remove duplicate Stations
# Alternative spellings of station ids (spelling case)
# lat-longitude
# short distance (40-50km as baseline?)


# Step 2: QA/QC
# Testing with random subsample to ensure acceptably low level of false positive rates for
# specified thresholds and tests for each variables

# How to handle existing QA/QC flags

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


# Step 3: Case study analysis for quality check
