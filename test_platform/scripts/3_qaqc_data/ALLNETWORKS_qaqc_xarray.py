import os
import tempfile
import argparse 

# Import qaqc stage calc functions
try:
    from calc_qaqc_xarray import *
except:
    print("Error importing calc_qaqc.py")

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
                       observations Platform, and is independent of network.""",
        epilog="""Possible stations:
                  [ASOSAWOS, CAHYDRO, CIMIS, CW3E, CDEC, CNRFC, CRN, CWOP, HADS, HNXWFO, 
                   HOLFUY, HPWREN, LOXWFOMAP, MTRWFO, NCAWOS, NOS-NWLON, NOS-PORTS, OtherISD, 
                   RAWS, SGXWFO, SHASAVAL, VCAPCD, MARITIME, NDBC, SCAN, SNOTEL]
               """
    )
    
    # Define arguments for the program
    parser.add_argument('-n', '--network', default="VCAPCD", help="Network name", type=str)
    parser.add_argument('-v', '--verbose', default=True, help="printing statemets throughout script", type=bool)
    parser.add_argument('-s', '--sample', default=1, help="number of stations within network to test", type=int)
    parser.add_argument('-f', '--full', default=False, help="if True, test all stations in network", type=bool)
    
    # Parse arguments
    args = parser.parse_args()
    network = args.network
    verbose = args.verbose
    sample = args.sample
    full = args.full
    if full:
        sample = None
        
    rawdir, cleandir, qaqcdir, mergedir = get_file_paths(network)
    
    whole_station_qaqc(network, cleandir, qaqcdir, verbose=verbose, sample=sample)

# Dev to do:
# reorder variables once entire qaqc is complete before saving
# output csv of flags/consistent flagging
# check the h5netcdf vs. netcdf4 engine
# delete testing notes