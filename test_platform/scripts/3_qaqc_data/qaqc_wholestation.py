#======================================================================
## Part 1a functions (whole station/network)
## Note: QA/QC functions in part 1a of whole station checks do not proceed through QA/QC if failure occurs

#----------------------------------------------------------------------
# missing value cehck: double check that all missing value observations are converted to NA before QA/QC
def qaqc_missing_vals(df, verbose=True):
    '''
    Checks data to be qaqc'ed for any errant missing values that made it through cleaning
    Converts those missing values to NAs
    Searches for missing values in 3_qaqc_data/missing_data_flags.csv
    '''

    missing_vals = pd.read_csv('missing_data_flags.csv')

    all_vars = [col for col in df.columns if 'qc' not in col]
    obs_vars = [var for var in all_vars if var not in ['lon','lat','time','elevation','station','anemometer_height_m','thermometer_height_m']]
    
    try:
        for item in obs_vars:
            # pull missing values which are appropriate for the range of real values for each variable 
            missing_codes = missing_vals.loc[missing_vals['variable'].str.contains(item) | missing_vals['variable'].str.contains('all')]

            # values in column that == missing_flag values, replace with NAs
            # note numerical vals converted to strings first to match missing_flag formatting
            df[item] = np.where(df[item].astype(str).isin(missing_codes['missing_flag']), float('NaN'), df[item])

            print('Updating missing values for: {}'.format(item))
    except:
        return None

    return df

#----------------------------------------------------------------------
# missing spatial coords (lat-lon)
def qaqc_missing_latlon(df, verbose=True):
    """
    Checks if latitude and longitude is missing for a station.
    If missing, station is flagged to not proceed through QA/QC.
    """

    # latitude or longitude
    variables = list(df.columns)
    if "lon" not in variables or "lat" not in variables:
        return None

    if df['lat'].isnull().all():
        return None

    if df['lon'].isnull().all():
        return None

    # df['lon'] = df['lon'].fillna(method="pad")
    # df['lat'] = df['lon'].fillna(method="pad")
    
    return df
        
#----------------------------------------------------------------------
# in bounds of WECC
def qaqc_within_wecc(df, verbose=True):
    """
    Checks if station is within terrestrial & marine WECC boundaries.
    If outside of boundaries, station is flagged to not proceed through QA/QC.
    """

    t = gp.read_file(wecc_terr).iloc[0].geometry  ## Read in terrestrial WECC shapefile.
    m = gp.read_file(wecc_mar).iloc[0].geometry   ## Read in marine WECC shapefile.
    pxy = shapely.geometry.Point(df['lon'].mean(), df['lat'].mean())
    if pxy.within(t) or pxy.within(m):
        return df
    else:
        return None

#----------------------------------------------------------------------
# elevation
def _grab_dem_elev_m(lats_to_check, lons_to_check, verbose=True):
    """
    Pulls elevation value from the USGS Elevation Point Query Service, 
    lat lon must be in decimal degrees (which it is after cleaning)
    Modified from: 
    https://gis.stackexchange.com/questions/338392/getting-elevation-for-multiple-lat-long-coordinates-in-python
    """
    url = r'https://epqs.nationalmap.gov/v1/json?'

    dem_elev_short = np.ones_like(lats_to_check)*np.nan
    
    for i,lat,lon in zip(range(len(lats_to_check)), lats_to_check, lons_to_check):
        # define rest query params
        params = {
            'output': 'json',
            'x': lon,
            'y': lat,
            'units': 'Meters'
        }

        # format query string and return value
        result = requests.get((url + urllib.parse.urlencode(params)))
        dem_elev_long = float(result.json()['value'])
        # make sure to round off lat-lon values so they are not improbably precise for our needs
        # dem_elev_short[i] = '{:.2f}'.format(dem_elev_long) 
        dem_elev_short[i] = np.round(dem_elev_long, decimals=2) 

    return dem_elev_short.astype("float")

#----------------------------------------------------------------------
def qaqc_elev_infill(df, verbose=True):
    """
    Checks if elevation is NA/missing. If missing, fill in elevation from either DEM or station.
    Some stations have all nan elevation values (e.g., NDBC, MARITIME)
    Some stations have single/few but not all nan elevation values (e.g., otherisd, asosawos)
    """    
    if verbose:
        print('Elevation values pre-infilling: {}'.format(df['elevation'].unique()))
        print('Elevation eraqc values pre-infilling: {}'.format(df['elevation_eraqc'].unique())) # testing

    # elev values missing
    isNan = df['elevation'].isnull()

    # if all are missing
    if isNan.any():
        if isNan.all():
            dem_elev_values = _grab_dem_elev_m([df['lat'].iloc[0]], 
                                               [df['lon'].iloc[0]],
                                               verbose=verbose)    
            df.loc[:, 'elevation'] = dem_elev_values[0]    
            df.loc[:, 'elevation_eraqc'] = 3    
        else:
            # if some missing
            try:
                if df['lon'].is_unique and df['lat'].is_unique:
                    dem_elev_values = _grab_dem_elev_m([df['lat'].iloc[0]], 
                                                       [df['lon'].iloc[0]],
                                                       verbose=verbose)
                    df.loc[:, 'elevation'] = dem_elev_values[0]
                    df.loc[:, 'elevation_eraqc'] = 3
                else:
                    dem_elev_values = _grab_dem_elev_m([df.loc[isNan, 'lat']], 
                                                       [df.loc[isNan, 'lon']],
                                                        verbose=verbose)
                    df.loc[isNan, 'elevation'] = dem_elev_values
                    df.loc[isNan, 'elevation_eraqc'] = 3
                return df
            
            # elevation cannot be obtained from DEM
            except:
                if verbose:
                    print("Elevation cannot be obtained from DEM")
                return None
    else:
        return df

#----------------------------------------------------------------------
def qaqc_elev_range(df, verbose=True):
    """
    Checks valid values to identify suspicious elevation values that are larger than 10m in difference
    Checks if valid elevation value is outside of range of reasonable values for WECC region.
    If outside range, station is flagged to not proceed through QA/QC.
    """
    # first check for suspicious values
    # elev_vals = df['elevation'].unique() # identify how many different elevation "values" are present

    # elevation values flagged as incorrectly coded
    # uses a threshold of 10m different from the station elevation to identify suspicious elevations
    # control = _grab_dem_elev_m([df.loc['lat'].iloc[0]], 
    #                            [df.loc['lon'].iloc[0]])[0]

    # TO DO:
    # This is the original version, but what about if the first value is bad? (look up to the DEM version)
    control = df['elevation'].iloc[0]
    isOff = np.logical_or(df['elevation'] > control + 10,
                          df['elevation'] < control - 10)
    if isOff.any():
        # in-fill if value is missing
        try:
            if df['lon'].is_unique and df['lat'].is_unique:
                dem_elev_values = _grab_dem_elev_m([df['lat'].iloc[0]], 
                                                   [df['lon'].iloc[0]])
                df.loc[ifOff, 'elevation'] = dem_elev_values[0]
                df.loc[isOff, 'elevation_eraqc'] = 3
            else:
                dem_elev_values = _grab_dem_elev_m([df.loc[ifOff, 'lat']], 
                                                   [df.loc[isOff, 'lon']])
                df.loc[ifOff, 'elevation'] = dem_elev_values
                df.loc[isOff, 'elevation_eraqc'] = 3

        # elevation cannot be obtained from DEM
        except:
            return None
    else:
        return df

    if verbose:
        print('Elevation values post-infilling/correcting: {}'.format(df['elevation'].unique())) # testing
        print('Elevation qaqc values post-infilling/correcting: {}'.format(df['elevation_eraqc'].unique())) # testing
    
    # then check for in range if value is present but outside of reasonable value range
    # death valley is 282 feet (85.9 m) below sea level, denali is ~6190 m
    
    isOut = df['elevation'].iloc[0] < -86.0 or df['elevation'].iloc[0]>6200.0
    if isOut.any():
        return None
    else:
        return df

#======================================================================
## Part 1b functions (whole station/network)
## Note: QA/QC functions in part 1b of whole station checks proceed through QA/QC if failure occurs

#----------------------------------------------------------------------
## flag values outside world records for North America
def qaqc_world_record(df, verbose=True):
    '''
    Checks if temperature, dewpoint, windspeed, or sea level pressure are outside North American world records
    If outside minimum or maximum records, flags values
    '''
    try:
        # world records from HadISD protocol, cross-checked with WMO database
        # https://wmo.asu.edu/content/world-meteorological-organization-global-weather-climate-extremes-archive
        T_X = {"North_America":329.92} #K
        T_N = {"North_America":210.15} #K
        D_X = {"North_America":329.85} #K
        D_N = {"North_America":173.15} #K
        W_X = {"North_America":113.2} #m/s
        W_N = {"North_America":0.} #m/s
        S_X = {"North_America":108330} #Pa
        S_N = {"North_America":87000} #Pa

        maxes = {"tas": T_X, "tdps": D_X, "tdps_derived": D_X, "sfcWind": W_X, "psl": S_X}
        mins = {"tas": T_N, "tdps": D_N, "tdps_derived": D_N, "sfcWind": W_N, "psl": S_N}

        # variable names to check against world record limits
        wr_vars = ['tas', 'tdps_derived', 'tdps', 'sfcWind', 'psl']

        for var in wr_vars:
            if var in list(df.columns):
                isOffRecord = np.logical_or(df[var] < mins[var]['North_America'],
                                            df[var] > maxes[var]['North_America'])
                if isOffRecord.any():
                    df.loc[isOffRecord, var + '_eraqc'] = 11
        return df
    except Exception as e:
        if verbose:
            print("qaqc_world_record failed with Exception: {}".format(e))
        return None

#----------------------------------------------------------------------
## sensor height - air temperature
def qaqc_sensor_height_t(df, verbose=True):
    '''
    Checks if temperature sensor height is within 2 meters above surface +/- 1/3 meter tolerance.
    If missing or outside range, temperature value for station is flagged to not proceed through QA/QC.
    '''
    try:
        # Check if thermometer height is missing   
        isHeightMissing = df['thermometer_height_m'].isnull().any()

        if isHeightMissing:
            # df.loc[:,'tas_eraqc'] = 6 # see era_qaqc_flag_meanings.csv
            df['tas_eraqc'] = 6 # see era_qaqc_flag_meanings.csv
        else:
            isHeightWithin = np.logical_and(df['thermometer_height_m'] >= (2 - 1/3),
                                            df['thermometer_height_m'] <= (2 + 1/3))

            # Thermometer height present but outside 10m +/- tolerance
            if not isHeightWithin:
                # df.loc[:, 'tas_eraqc'] = 7
                df['tas_eraqc'] = 7

        return df
    except Exception as e:
        if verbose:
            print("qaqc_sensor_height_w failed with Exception: {}".format(e))
        return None

#----------------------------------------------------------------------
## sensor height - wind
def qaqc_sensor_height_w(df, verbose=True):
    '''
    Checks if wind sensor height is within 10 meters above surface +/- 1/3 meter tolerance.
    If missing or outside range, wind speed and direction values for station are flagged to not proceed through QA/QC.
    '''
    # try:
    if True:
        # Check if anemometer height is missing
        isHeightMissing = df['anemometer_height_m'].isnull().any()

        if isHeightMissing:
            # df.loc[:,'sfcWind_eraqc'] = 8 # see era_qaqc_flag_meanings.csv
            # df.loc[:,'sfcWind_dir_eraqc'] = 8
            df['sfcWind_eraqc'] = 8 # see era_qaqc_flag_meanings.csv
            df['sfcWind_dir_eraqc'] = 8

        else: # sensor height present
            # Check if anemometer height is within 10 m +/- 1/3 m
            isHeightWithin = df['anemometer_height_m'][0] >= (10 - 1/3) and df['anemometer_height_m'][0] <= (10 + 1/3)
            # Anemometer height present but outside 10m +/- tolerance
            if not isHeightWithin:
                df['sfcWind_eraqc'] = 9
                df['sfcWind_dir_eraqc'] = 9 
        return df

    else:
    # except Exception as e:
        if verbose:
            print("qaqc_sensor_height_w failed with Exception: {}".format(e))
        return None