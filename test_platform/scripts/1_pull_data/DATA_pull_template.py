## Historical Observations Platform
## Template Script for Stage: PULL DATA

## TO DO LIST
## Any notes critical for further development, e.g.: AWS implementation

# Step 0: Environment set-up
# Import libraries -- delete as appropriate per datasource
import os
from datetime import datetime, timezone
import xarray as xr
from ftplib import FTP


# Set envr variables
workdir = "/path/to/working/directory/"
years = list(map(str,range(1980,datetim.enow().year+1))) # If needed

# Step 1: Download data
## This may happen via various methods depending on how the data is stored

# ------------------------------------------------------------------------------------------------
# FTP PROCESS
def get_DATA(workdir, years, OTHER_FLAGS):
    # Login, using ftplib
    ftp = FTP('ftp_link_here')
    ftp.login() # follow user/password for Hist-Observations
    ftp.cwd('working/dir/for/data/goes/here')
    pwd = ftp.pwd() # gets base path

    for i in years:
        ftp.cwd(pwd)
        filenames = ftp.nlst() # Returns list of all files in folder
        filenames = list(filter(lambda f: f.endswith('.nc'), filenames)) # Only keep .nc filenames

        for filename in filenames:
            # DOWNLOAD DATA STEPS

            # Important checks to (optionally) include when testing
            # Check 1: If new data download, download all available data
            #          If updated data download, only download data since last data pull
            local_filename = os.path.join(workdir, filename)
            file = open(local_filename, 'wb') # Open destination file
            ftp.retrbinary('RETR ' + filename, file.write) # Write file -- RENAME CONVENTIONS?

            # Check 2: As data downlaods, print to command line
            print('{} saved'.format(filename))

            # Check 3: Error handling if data missing
            except Exception as e:
                print("Error in downloading date {}/{}: {}".format(FORMAT))

    ftp.quit() # this is the polite way to close a ftp connection

get_DATA(workdir, years=years, OTHER_FLAGS) # To download data, run

# ------------------------------------------------------------------------------------------------
# API PROCESS
# Using ASOS as example here

states = ["STATE", "STATE2", "STATE3"]
def get_DATA(OPTIONAL_FLAGS):
    # DOWNLOAD DATA STEPS
    isd_station_list = pd.DataFrame(columns = ['cols_to_fill_here', 'and_another_one'])
    for i in states:
        try:
            URL = "URL_goes_here".format(i)
            stations = requests.get(url).json()
            print("Stations loaded.")

            for j in stations['stationCollection']['stations']:
                try:
                    StnId = j['LABEL']
                    Name = j['header']['preferredName']
                    Latitude = j['header']['latitude_dec']
                    Longitude = j['header']['longitude_dec']
                    State = i
                    BeginDate = j['header']['por']['beginDate']
                    EndDate = j['header']['por']['endDate']
                    platform = []

            # Check 3: Error handling if data missing
                except Exception as e:
                    print("Error in downloading date {}/{}: {}".format(FORMAT))

get_DATA()

# ------------------------------------------------------------------------------------------------
# GOOGLE DRIVE PROCESS -- using .gz files as example here that were zipped for batch download
def unzip_dir(directory):
    os.chdir(directory)
    extension = ".zip"
    for item in os.listdir(directory):
        if item.endswith(extension):
            zip_name = os.path.abspath(item)
            print("Unzipping year (folder): ", zip_name)
            ZipFile(zip_name).extracall(directory) # Unzips all files in directory ###
            os.remove(zip_file)

unzip_dir(datadir)

def gz_extract(directory):
    extension = ".gz" # may be .gz or other
    os.chdir(directory)
    for item in os.listdir(directory):
        if item.endswith(extension):
            EXT_name = os.path.abspath(item)
            filename = (os.path.basename(EXT_name)).rsplit('.',1)[0]
            with gzip.open(EXT_name, 'rb') as f_in, open(filename, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
            os.remove(EXT_name) # Removes .gz files, cleans up directory for ease

            # Check 2: As data downlaods, print to command line
            print("{} saved".format(filename))

gz_extract(DIRECTORY)
