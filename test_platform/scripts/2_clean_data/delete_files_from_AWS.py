"""This script deletes files in batch from AWS. Options for which files to delete are provided. """

import pandas as pd
import s3fs
import boto3
from cleaning_helpers import get_file_paths

s3 = boto3.resource('s3')


# function to batch delete files from AWS
def delete_files_from_AWS(network, which_to_delete):
    """
    This function should be used CAREFULLY as it will PERMANENTLY batch delete files from AWS!!!

    Options of which files to delete:
    which_to_delete = "empty"
        When grabbing data for an updated pull from Synoptic/MADIS networks for end of 2022 data, 
        our update_pull script produced an update raw data file for every station, regardless if 
        that station had new data to add. Therefore, there are many "empty" raw data files in our 
        /1_raw_wx bucket that will cause issues in cleaning and updating the cleaned station list 
        for cleaned variable coverage. 
        This option *specifically* identifies these empty files and deletes them. 
        Note: we are no longer pulling data from Synoptic in the future.

    which_to_delete = "update_pull"
        When grabbing data for an updated pull from Synoptic/MADIS networks for end of 2022 data, 
        our update_pull script produced an update raw data file for every station. This file has a 
        suffix with the date of the update pull. When the next full pull is performed, the original 
        raw data files will be overwritten, but the old update pull files with the date suffix will
        not, even though that data should be encapsulated in the latest new full pull. 
        This option *specifically* identifies the files with a date suffix and deletes them. 
        Only perform after a new full pull has occured. 

    which_to_delete = "monthly"
        NDBC, MARITIME networks produce a monthly file if the year of data coverage is within the
        current year, rather than a single annual file. Once the annual file is available for the year,
        these monthly files should be deleted.

    which_to_delete = "no_wx_data" ?
        Could be set-up to delete raw files that do not contain any meteorological info, i.e., buoys
        Not continuing for now.

    Note: if which_to_delete is not provided an option, nothing will delete. 

    Modified from: https://stackoverflow.com/questions/32501995/boto3-s3-renaming-an-object-using-copy-object 
    """

    bucket_name  = "wecc-historical-wx"

    # deletes only update pull files with date extensions that are empty
    if which_to_delete == "empty":

        # grab raw data directory for specific network
        rawdir, cleandir, qaqcdir = get_file_paths(network)

        # get files
        all_files = []
        for item in s3.Bucket(bucket_name).objects.filter(Prefix = rawdir):
            file = str(item.key)
            all_files += [file]

        # identify files with date in filename
        files_to_delete = []
        for file in all_files:
            file_parts = file.split("/")
            stn_file = file_parts[-1] # grabs just file name, no path
            stn_file = stn_file.split(".")[0] # removes file extenstion

            if ("_" in stn_file) and "errors" not in stn_file and "stationlist" not in stn_file:
                fn_date = file.split("_")[-1] # grabs part of filename after _ character
                if "-" in fn_date: # must follow YYYY-MM-DD format
                    # print(file)
                    files_to_delete += [file]  

        # open each file, and assess if empty
        empty = []
        for file in files_to_delete: 
            print('Checking {} for empty size'.format(file))
            fs = s3fs.S3FileSystem()
            aws_url = "s3://" + bucket_name + "/" + file

            with fs.open(aws_url) as fileObj:
                df = pd.read_csv(fileObj, header=6, low_memory=False)

                # assess if empty
                if len(df.index) == 1 and df['Station_ID'].isnull().values.any(): # any station that is missing the station_id is empty
                    empty += [file]

        print('Number of empty files to delete in {0}: {1}'.format(rawdir, len(empty)))
        print(empty)

        ## UNCOMMENT WHEN READY TO ACTUALLY RUN -- DO NOT TEST/RUN THIS PART UNTIL FILES ARE CONFIRMED, IT WILL DELETE FILES
        ## delete files in empty -- commenting out for now/safety
        # for file in empty:
        #     s3.Object(bucket_name, file).delete() # delete file
        #     print('File {} deleted from AWS'.format(file))       




    # deletes ALL update pull files with date extensions
    elif which_to_delete == "update_pull":

        # grab raw data directory for specific network
        rawdir, cleandir, qaqcdir = get_file_paths(network)

        # get files
        all_files = []
        for item in s3.Bucket(bucket_name).objects.filter(Prefix = rawdir):
            file = str(item.key)
            all_files += [file]

        # identify files with date in filename
        files_to_delete = []
        for file in all_files:
            file_parts = file.split("/")
            stn_file = file_parts[-1] # grabs just file name, no path
            stn_file = stn_file.split(".")[0] # removes file extenstion

            if ("_" in stn_file) and "errors" not in stn_file and "stationlist" not in stn_file:
                fn_date = file.split("_")[-1] # grabs part of filename after _ character
                if "-" in fn_date: # must follow YYYY-MM-DD format
                    files_to_delete += [file]

        print('Number of files to delete in {0}: {1}'.format(rawdir, len(files_to_delete)))
        print(files_to_delete) # should have date in filename

        ## UNCOMMENT WHEN READY TO ACTUALLY RUN -- DO NOT TEST/RUN THIS PART UNTIL FILES ARE CONFIRMED, IT WILL DELETE FILES
        ## delete files in files_to_delete -- commenting out for now/safety
        # for file in files_to_delete:
        #     s3.Object(bucket_name, file).delete() # delete file
        #     print('File {} deleted from AWS'.format(file))



    # deletes the monthly files from NDBC, MARITIME after year file is in
    # note: CIMIS has monthly files, but should not delete, as these files retain data for all stations
    elif which_to_delete == "monthly":

        # grab raw data directory for specific network
        rawdir, cleandir, qaqcdir = get_file_paths(network)

        if network == "NDBC" or network == "MARITIME":

            # get files
            all_files = []
            for item in s3.Bucket(bucket_name).objects.filter(Prefix = rawdir):
                file = str(item.key)
                all_files += [file]

            all_files = [file for file in all_files if 'stationlist' not in file]
            all_files = [file for file in all_files if 'error' not in file]

            # identify monthly files
            files_to_delete = []
            for file in all_files:
                # print(file)
                file_parts = file.split("/")
                stn_file = file_parts[-1] # grabs just file name, no path
                stn_file = stn_file.split(".")[0] # removes file extenstion (.txt.gz, or .zip)

                if "csv" in stn_file:
                    continue # skip canadian files, they're annual

                else:
                    stn_year = stn_file[-4:] # use to identify if annual file is already downloaded
                    stn_name = stn_file[:5] # NDBC and MARITIME station names are always 5 characters long

                    yr_file_to_check = stn_name + "h" + stn_year
                    if stn_file[-5] != "h":
                        if yr_file_to_check in all_files: 
                            files_to_delete += [file]
                    
                        else: # annual file for that station is not downloaded yet
                            print('Warning! The annual file for this {} is not yet downloaded. Please check before proceeding.'.format(stn_name))

            print('Number of files to delete in {0}: {1}'.format(rawdir, len(files_to_delete)))
            print(files_to_delete) # should have date in filename

            ## UNCOMMENT WHEN READY TO ACTUALLY RUN -- DO NOT TEST/RUN THIS PART UNTIL FILES ARE CONFIRMED, IT WILL DELETE FILES
            ## delete files in files_to_delete -- commenting out for now/safety
            # for file in files_to_delete:
            #     s3.Object(bucket_name, file).delete() # delete file
            #     print('File {} deleted from AWS'.format(file))

                    

    elif which_to_delete == None:
        print('Please provide an option to delete files from AWS. Options are "empty", "monthly", and "update_pull". Check documentation for meanings.')


if __name__ == "__main__":
    network = "HNXWFO"
    delete_files_from_AWS(network, which_to_delete="empty")
    # intentionally going to have user do this network by network, as we should be careful as to what we delete for now


# -----------------------------------------------------------------------------------------------------------------------------

## Old code
## function can be renamed to be "delete files from AWS" in future

# def rename_MARITIME(bucket_name, cleandir): 
#     """
#     Ancilliary one-time help to rename MARITIME cleaned files which was previously saved with lowercase station names. 
#     Will be fixed in future cleaning update for MARITIME network. 

#     *** NOTE ***
#     Function moved into cleaning_helpers.py as may be helpful in future for automation/delete empty files from Synoptic update
#     in a modified format (eg the delete line, and constructed to look for small file sizes).
#     """
    
#     # source: https://stackoverflow.com/questions/32501995/boto3-s3-renaming-an-object-using-copy-object

#     files = [] # Get files in MARITIME
#     for item in s3.Bucket(bucket_name).objects.filter(Prefix = cleandir):
#         file = str(item.key)
#         files += [file]

#     files = list(filter(lambda f: f.endswith(".nc"), files)) # Get list of cleaned file names

#     for file in files:

#         # example old file name: 2_clean_wx/MARITIME/MARITIME_wptw1.nc
#         file_parts = file.split("/")
#         stn_name = file_parts[-1]
#         stn_name = stn_name.split(".")[0] # remove file extension

#         # example new file name: 2_clean_wx/MARITIME/MARITIME_WPTW1.nc
#         file_upper = cleandir + stn_name.upper() + ".nc"

#         s3.Object(bucket_name, file_upper).copy_from(CopySource=bucket_name+'/'+file) # rename to upper case
#         print('{0} renamed to {1}'.format(file, file_upper))

#         s3.Object(bucket_name, file).delete() # delete lowercase name files
#         print('Old file {} deleted from AWS'.format(file))