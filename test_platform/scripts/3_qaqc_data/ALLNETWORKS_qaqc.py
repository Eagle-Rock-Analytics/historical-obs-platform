"""
This script performs qa/qc protocols for cleaned station data for ingestion into the Historical Observations Platform, and is
independent of network. 
Approach:
(1) Remove duplicate stations
(2) Handle variables that report at different intervals and/or change frequency over time (convert to hourly?)
(3) QA/QC testing, including consistency checks, gaps, checks against climatological distributions, and cross variable checks.
(4) Case study analysis for accuracy -- SHOULD THIS BE A SEPARATE SCRIPT/PROCESS?

Inputs: Cleaned data for an individual network
Outputs: QA/QC-processed data for an individual network, priority variables, all times. Organized by station as .nc file.
"""

import os
import tempfile
import argparse

# Import all qaqc script functions
try:
    from qaqc_plot import *
    from qaqc_utils import *
    from qaqc_wholestation import *
    from qaqc_logic_checks import *
    from qaqc_buoy_check import *
    from qaqc_frequent import *
    from qaqc_unusual_gaps import *
    from qaqc_unusual_large_jumps import *
    from qaqc_climatological_outlier import *
    from qaqc_unusual_streaks import *
except Exception as e:
    print("Error importing qaqc script: {}".format(e))

# Import qaqc stage calc functions
try:
    from QAQC_pipeline import *
except:
    print("Error importing QAQC_pipeline.py")
    
# =================================================================================================
# Main Function
if __name__ == "__main__":
     
    # Create parser
    parser = argparse.ArgumentParser(
        prog="ALLNETWORKS_qaqc",
        description="""This script performs qa/qc protocols for cleaned station data for ingestion into the Historical 
                       Observations Platform, and is independent of network.""",
        epilog="""Possible stations:
                  [ASOSAWOS, CAHYDRO, CIMIS, CW3E, CDEC, CNRFC, CRN, CWOP, HADS, HNXWFO, 
                   HOLFUY, HPWREN, LOXWFO, MAP, MTRWFO, NCAWOS, NOS-NWLON, NOS-PORTS, OtherISD, 
                   RAWS, SGXWFO, SHASAVAL, VCAPCD, MARITIME, NDBC, SCAN, SNOTEL]
               """
    )
    
    # Define arguments for the program
    parser.add_argument('-n', '--network', default="TRAINING", help="Network name (default to 'VCAPCD')", type=str)
    parser.add_argument('-l', '--local', default=False, help="Save files and plots locally (default to False)", type=bool)
    parser.add_argument('-r', '--rad_scheme', default="remove_zeros", help="Radiation handling scheme for frequent values check. See qaqc_frequent_values for options (default to 'remove_zeros'", type=str)
    parser.add_argument('-v', '--verbose', default=False, help="Print statements throughout script (default to False)", type=bool)

    # Parse arguments
    args = parser.parse_args()
    network = args.network
    if network.lower()=="training":
        network="TRAINING"
    rad_scheme = args.rad_scheme
    verbose = args.verbose
    local = args.local

    if network=="NETWORK:
        rawdir,cleandir,qaqcdir,mergedir = None,None,None,None
    else:
        rawdir, cleandir, qaqcdir, mergedir = get_file_paths(network)
    whole_station_qaqc(network, cleandir, qaqcdir, rad_scheme, verbose=verbose, local=local)

# Dev to do:
# reorder variables once entire qaqc is complete before saving
# output csv of flags/consistent flagging
# check the h5netcdf vs. netcdf4 engine
# delete testing notes
