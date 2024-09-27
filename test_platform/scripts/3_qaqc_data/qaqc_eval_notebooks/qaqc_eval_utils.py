'''Functions used across QA/QC evaluation protocol for Historical Data Platform'''

from pyproj import CRS, Transformer
import geopandas as gpd
from geopandas import GeoDataFrame
from shapely.geometry import Point, Polygon
import s3fs

import xarray as xr
import numpy as np
import pandas as pd
import sys
import os

import matplotlib.pyplot as plt
import cartopy.feature as cf
from matplotlib.ticker import MaxNLocator
import cartopy.crs as ccrs

import datetime

# from qaqc_eval_plot import event_plot, latlon_to_mercator_cartopy

sys.path.append(os.path.expanduser('../'))
from qaqc_plot import flagged_timeseries_plot, _plot_format_helper, id_flag
from QAQC_pipeline import qaqc_ds_to_df


def known_issue_check(network, var, stn):
    '''
    Identifies if station under evaluation has a known network issue.
    At present, only prints out a statement if there is an issue.
    Eventually may want to do <something>

    Note: See "Known Network Issues for QA/QC Validation" planning doc.
    '''
    print('Checking for known station issues...')
    # RAWS
    if network == 'RAWS':
        if var == 'tas':
            print('Known network issue for {} {}: values may be too high (on order of 10°F) if sun is shining strongly and winds are light.'.format(
                network, var))

        elif var == 'pr':
            print('Known network issue for {} {}: stations are not maintained in winter, instrument may freeze. Consider subsetting for May-September.'.format(
                network, var))
            # V2 note: exclude RAWS data during specific notes -- would require new function to flag

    # SNOTEL
    if network == 'SNOTEL':
        if var == 'tas':
            print('Known network issue for {} {}: values may remain at exactly 0.0°C for two or more consecutive days. Should be caught by unusual_streaks.'.format(
                network, var))
            print('Known network issue for {} {}: SNOTEL temperature sensors transition between mid-1990s and mid-2000s to new sensory type produces warm bias at \
            colder temperatures. Min temperature may be too high, max temperature may be too low.'.format(
                network, var))
            # V2 note: trend analysis may identify these issues, nearest neighbor check could identify

    # ASOSAWOS + OtherISD
    if network == 'ASOSAWOS':
        if var == 'tdps':
            print('Known network issue for {} {}: values may be stuck at around 0.0°C, or have excessive mirror contamination. Should be caught by unsusual_streaks.'.format(
                network, var))
    
    if network == 'ASOSAWOS' or network == 'OtherISD':
        if var == 'pr':
            print('Known network issue for {} {}: ASOS network began installation in 1996, with poor instrumentation for measuring snowfall. Precipitation between \
            1980-1996 may be more likely to be flagged.'.format(
                network, var))

    # CIMIS
    if network == 'CIMIS':
        if var == 'pr':
            print('Known network issue for {} {}: stations located in flat agricultural areas, sensor may be detecting sprinkler irrigation events. \
            Network does have stringent QC protocol.'.format(
                network, var))
            # V2 note: nearest neighbor check could confirm

    # NDBC / MARITIME
    if network == 'NDBC' or network == 'MARITIME':
        print('Known network issue for {}: some buoys have data past their known disestablishment dates. Should be caught by spurious_buoy_check.'.format(
            network))

        if stn == 'NDBC_46044':
            print('Known network issue for {} station NDBC_46044: buoy went adrift during reporting period. Confirm if data was flagged by QA/QC.'.format(
                network, stn))
            # V2 note: if not flagged, needs to be -- would require new function
        
        if stn == 'MARITIME_MTYC1' or stn == 'MARITIME_MEYC1' or stn == 'MARITIME_SMOC1' or stn == 'MARITIME_ICAC1':
            print('Known network issue for {} station {}: buoy was renamed and/or relocated. May cause issue for station proximity tests.'.format(
                network, stn))
            # V2 note: noted in qaqc_buoy_check but not handled -- would require new function


def subset_eval_stns(event_to_eval, stn_list, subset=None, return_stn_ids=False):
    '''
    Identifies stations to evaluate for specific V1 QA/QC events.
    Option to subset to a more manageable number of random stations for initial evaluation. 
    '''
    
    # TO DO: validation check on event_to_eval options

    event_flags = []
    event_flags.append('all')
    event_flags.append(event_to_eval) # options: santa_ana_wind, winter_storm, AR, mudslide, aug2020_heatwave, sep2020_heatwave, aug2022_heatwave, offshore_wind

    # grab stations per event
    event_stns = stn_list[stn_list['event_type'].isin(event_flags)]

    # exclude "manual check on end date" for time being -- SNOTEl stations all have 2100 as their end date regardless of when data actually ends
    mask = event_stns['notes'] == 'manual check on end date'
    event_stns = event_stns[~mask]
    # print('{} potential stations available for evaluation for {} event!'.format(len(event_stns), event_to_eval))

    # identify stations in geographic region we are looking for
    census_shp_dir = "s3://wecc-historical-wx/0_maps/ca_counties/" 
    ca_county = gpd.read_file(census_shp_dir)

    # different areas based on events
    ## TODO: Need an option for "WECC wide" (or no spatial subsetting)
    if event_to_eval == 'santa_ana_wind':
        counties_to_grab = ['Los Angeles', 'Orange', 'San Diego', 'San Bernardino', 'Riverside'] ## potentially need to broaden area? 

    elif event_to_eval == 'winter_storm': # focus on Northern/Central/Bay Area to begin with // WECC wide
        counties_to_grab = ['Butte', 'Colusa', 'Del Norte', 'Glenn', 'Humboldt', 'Lake', 'Lassen', 
        'Mendocino', 'Modoc', 'Nevada', 'Plumas', 'Shasta', 'Sierra', 'Siskiyou', 'Tehama','Trinity',
        'Alpine', 'Amador', 'Calaveras', 'El Dorado', 'Fresno', 'Inyo', 'Kings', 'Madera', 'Mariposa',
        'Merced', 'Mono', 'Placer', 'Sacramento', 'San Joaquin', 'Stanislaus', 'Sutter', 'Yuba', 'Tulare',
        'Tuolumne', 'Yolo', 'Alameda', 'Contra Costa','Marin', 'Monterey','Napa', 'San Benito', 'San Francisco',
        'San Mateo', 'Santa Clara', 'Santa Cruz', 'Solano', 'Sonoma']

    elif event_to_eval == 'mudslide':
        counties_to_grab = ['Santa Barbara']

    elif event_to_eval == 'AR':
        counties_to_grab = [] # CA

    elif event_to_eval == 'aug2020_heatwave': # August 2020 "aug2020_heatwave" -- 
        counties_to_grab = [] # CA

    elif event_to_eval == 'sep2020_heatwave': # September 2020 "sep2020_heatwave"
        counties_to_grab = ['San Luis Obispo', 'Kern', 'San Bernadino', 'Santa Barbara', 'Ventura',
        'Los Angeles', 'Orange', 'Riverside', 'San Diego', 'Imperial']

    elif event_to_eval == 'aug2022_heatwave': # August 2022 -- Labor Day Heatwave "aug2022_heatwave"
        counties_to_grab = ['San Luis Obispo', 'Kern', 'San Bernadino', 'Santa Barbara', 'Ventura',
        'Los Angeles', 'Orange', 'Riverside', 'San Diego', 'Imperial']

    elif event_to_eval == 'offshore_wind':
        counties_to_grab = ['San Diego', 'Orange', 'Los Angeles', 'Ventura', 'Santa Barbara',
        'San Luis Obispo', 'Monterey', 'Santa Cruz', 'San Mateo', 'Santa Clara', 'Alameda',
        'San Francisco', 'Contra Costa', 'Solano', 'Marin', 'Sonoma', 'Mendocino', 'Humboldt', 'Del Norte']

    target_counties = ca_county[ca_county['NAME'].isin(counties_to_grab)]
    target_counties = GeoDataFrame(target_counties, geometry=target_counties.geometry)

    geometry = [Point(latlon_to_mercator_cartopy(lat,lon)) for lat,lon in zip (event_stns.latitude, event_stns.longitude)]
    event_stns = GeoDataFrame(event_stns,geometry=geometry).set_crs(crs="EPSG:3857", allow_override=True) # adding geometry column
    event_stns_local = gpd.overlay(event_stns, target_counties, how="intersection") # subsetting for stations within county boundaries
    print('{} potential stations available for evaluation for {} event.'.format(len(event_stns_local), event_to_eval))

    if subset != None:
        if len(event_stns_local) <= subset:
            eval_stns = event_stns_local
        else:
            eval_stns = event_stns_local.sample(subset, replace=False)
            print('{} stations selected for evaluation for {} event!'.format(subset, event_to_eval))
    else:
        eval_stns = event_stns_local

    # lastly, check if there are any known issues --- need to refactor for check
    # check_networks = event_stns.network.unique()
    # known_issue_check(event_stns, var='tas')

    if return_stn_ids:
        print('Stations selected for evaluation:\n', list(eval_stns['era-id']))

    return eval_stns


def id_all_flags(ds):
    '''Prints all unique values of all eraqaqc flags'''
    ds_vars = list(ds.keys())
    qc_vars = [i for i in ds_vars if '_eraqc' in i]
    if len(qc_vars) == 0:
        print('Station has no eraqc variables -- please double check that this station has completed QA/QC!')
    else:
        for var in qc_vars:
            print(var, np.unique(ds[var].data))


def pull_nc_from_aws(fname):
    print('Retrieving data for station...')
    s3 = s3fs.S3FileSystem(anon=False)
    network = fname.split('_')[0]
    s3_url = 's3://wecc-historical-wx/3_qaqc_wx_dev/{}/{}.nc'.format(network, fname)

    try:
        s3_file_obj = s3.open(s3_url, mode='rb')
        ds = xr.open_dataset(s3_file_obj, engine='h5netcdf')
        return ds
        
    except:
        print(f'Station {fname} not found in bucket -- please check if station completed QA/QC.')


def event_info(event, alt_start_date=None, alt_end_date=None):
    start_date = {
        "santa_ana_wind"   : "1988-02-16",
        "winter_storm"     : "1990-12-20", 
        "AR"               : "2017-01-16",
        "mudslide"         : "2018-01-05",
        "aug2020_heatwave" : "2020-08-14",
        "sep2020_heatwave" : "2020-09-05",
        "aug2022_heatwave" : "2022-08-30",
        "offshore_wind"    : "2021-01-15",
        "alternative"      : alt_start_date 
                  }
    end_date = {
        "santa_ana_wind"   : "1988-02-19",
        "winter_storm"     : "1990-12-24",
        "AR"               : "2017-01-20",
        "mudslide"         : "2018-01-09",
        "aug2020_heatwave" : "2020-08-15",
        "sep2020_heatwave" : "2020-09-08",
        "aug2022_heatwave" : "2022-09-09",
        "offshore_wind"    : "2021-01-16",
        "alternative"      : alt_end_date 
    }

    event_start = start_date[event]
    event_end = end_date[event]

    return (event_start, event_end)


def event_subset(df, event, buffer=7, alt_start_date=None, alt_end_date=None):
    """Subsets for the event itself + buffer around to identify event"""
    print('Subsetting station record for event duration with {} day buffer...'.format(str(buffer)))

    df['time'] = pd.to_datetime(df['time']) # set to searchable datetime
    event_start, event_end = event_info(event, alt_start_date, alt_end_date) # grab dates from lookup dictionary
    
    datemask = ((df['time'] >= (pd.Timestamp(event_start) - datetime.timedelta(days=buffer))) & (df['time'] <= (pd.Timestamp(event_end) + datetime.timedelta(days=buffer)))) # subset for event dates + buffer
    event_sub = df.loc[datemask]
    
    return event_sub


def flags_during_event(subset_df, var, event):
    """Provides info on which flags were placed during event for evaluation"""
    print('Flags set on {} during {} event: {}'.format(var, event, subset_df[var+'_eraqc'].unique()))
    all_event_flags = []
    for item in subset_df[var+'_eraqc'].unique():
        all_event_flags.append(item)
    return all_event_flags


def multi_stn_check(list_of_stations, event, buffer=7, alt_start_date=None, alt_end_date=None):
    """this function does all the major identification steps outlined in the notebook"""
    for stn in list_of_stations:
        print('Evaluation on {}...'.format(stn))

        # retrieve data
        try:
            ds = pull_nc_from_aws(stn)
        except:
            continue

        # convert to dataframe
        print('Converting to dataframe...')
        df, MultiIndex, attrs, var_attrs, era_qc_vars = qaqc_ds_to_df(ds)

        # identify vars for evaluation
        vars_to_check = ['tas', 'hurs', 'sfcWind', 'sfcWind_dir']
        vars_to_eval = [var for var in vars_to_check if var in df.columns] # check if variable is not present in the specific station

        for var in vars_to_eval:
            known_issue_check(network=df.station.unique()[0].split('_')[0], 
                            var=var, 
                            stn=df.station.unique()[0]) # check if known issues are present first!
            print('Evaluating: {}'.format(var))
            flagged_timeseries_plot(df, var=var)
        
        # subset for event
        subset_df = event_subset(df, event, buffer, alt_start_date, alt_end_date)

        if len(subset_df) != 0:
            for v in vars_to_eval:
                all_flags = flags_during_event(subset_df, var=v, event=event)
                event_plot(subset_df, var=v, event=event, alt_start_date=alt_start_date, alt_end_date=alt_end_date)
        else:
            ds.close()
        
        # # if all are none or empty, close ds and move on
        # if len(subset_df) == 0 or np.isnan(all_flags[0]):
        #     ds.close()
        #     print('Closing dataset!\n')
        
        # else:
        #     # proceed
        #     print('{} is flagged during {}!'.format(stn, event))

def find_other_events(df, event_start, event_end, buffer=7, subset=None, return_stn_ids=True):
    print('Subsetting station record for event duration with {} day buffer...'.format(str(buffer)))
    
    df['start_date'] = pd.to_datetime(df['start_date'])
    df['end_date'] = pd.to_datetime(df['end_date'])
    event_start = pd.to_datetime(event_start).tz_localize('UTC')
    event_end = pd.to_datetime(event_end).tz_localize('UTC')
    
    event_sub = df.loc[(df['start_date'] <= (event_start - datetime.timedelta(days=buffer))) & (df['end_date'] >= (event_end + datetime.timedelta(days=buffer)))]

    # exclude "manual check on end date" stations since we don't know when they actually end
    event_sub = event_sub.loc[event_sub['notes'] != 'manual check on end date']

    # subset to make more manageable
    if subset != None:
        if len(event_sub) <= subset:
            eval_stns = event_sub
        else:
            eval_stns = event_sub.sample(subset, replace=False)
            print('{} stations selected for evaluation for comparison!'.format(subset))
    else:
        eval_stns = event_sub

    # return station ids for ease
    if return_stn_ids:
        print('Stations selected for evaluation:\n', list(eval_stns['era-id']))

    return eval_stns




# def return_ghcn_vars(ghcn_df, input_var):
#     '''
#     Given an input variable, return GHCNh location variables and all relevant data variables,
#     rather than utilizing the whole 240 cols, or having to know how ghcnh labels the cols.

#     input_var must follow ERA naming scheme (tas, tdps, ps, pr, etc.)
#     '''
#     ghcnh_vars = pd.read_csv('ghcnh_data_headers.csv')

#     # include station-ID, time, loc, elevation (cols 1-10)
#     stn_info_cols = ['Station_ID', 'Station_name',
#                      'Year','Month','Day','Hour','Minute',
#                      'Latitude','Longitude','Elevation']
    
#     var_cols = []
#     if input_var == 'tas':
#         varquery = 'temperature'
        
#     elif input_var == 'tdps' or 'tdps_derived':
#         varquery = 'dew_point_temperature'
        
#     elif input_var == 'ps' or 'psl':
#         varquery = 'station_level_pressure'
        
#     elif input_var == 'sfcWind_dir':
#         varquery = 'wind_direction'
        
#     elif input_var == 'sfcWind':
#         varquery = ['wind_speed', 'wind_gust']

#     elif input_var == 'hurs':
#         varquery = 'relative_humidity'
        
#     elif input_var == 'rsds':
#         print('GHCNh data does not have solar radiation data to evaluate against.')
#         varquery = '' 
        
#     elif input_var == 'pr' or input_var == 'pr_1h' or input_var == 'pr_5min':
#         varquery = 'precipitation'

#     i = ghcn_df.query()
    
#     var_cols = [i for i in ghcnh_vars if varquery in i]
#     cols_to_return = stn_info_cols + var_cols

#     return ghcn_df[[cols_to_return]]

def return_ghcn_vars(ghcn_df, input_var):
    '''
    Given an input variable, return GHCNh location variables and all relevant data variables,
    rather than utilizing the whole 240 cols, or having to know how ghcnh labels the cols.

    input_var must follow ERA naming scheme (tas, tdps, ps, pr, etc.)
    '''
    ghcnh_vars = pd.read_csv('ghcnh_data_headers.csv')

    # include station-ID, time, loc, elevation (cols 1-10)
    stn_info_cols = ['Station_ID', 'Station_name',
                     'Year','Month','Day','Hour','Minute',
                     'Latitude','Longitude','Elevation']
    vars = {
        'tas': 'temperature',
        'tdps': 'dew_point_temperature',
        'tdps_derived': 'dew_point_temperature',
        'ps': 'station_level_pressure',
        'psl': 'station_level_pressure',
        'sfcWind_dir': 'wind_direction',
        'sfcWind': 'wind_speed',
        'tas': 'temperature',
        'hurs': 'relative_humidity',
        'rsds': "N/A",
        'pr': 'precipitation',
        'pr_1h': 'precipitation',
        'pr_5min': 'precipitation',
    }
    if input_var in vars.keys():
        i = ghcn_df.columns.get_loc(vars[input_var])
        j = i+6
        # For wind, include wind gust
        if input_var=="sfcWind":
            j = j+6
        ghcn_df.iloc[:, i:j]
        
        return ghcn_df.iloc[:, i:j]
    else:
        raise Exception(f"Variable {input_var} not in variables' dictionary")



# projection stuffs
census_shp_dir = "s3://wecc-historical-wx/0_maps/ca_counties/" 
ca_county = gpd.read_file(census_shp_dir) # from s3 bucket

def latlon_to_mercator_cartopy(lat, lon):

    proj_latlon = CRS('EPSG:4326')
    proj_mercator = CRS('EPSG:3857')
    
    # Transform the coordinates
    transformer = Transformer.from_crs(proj_latlon, proj_mercator, always_xy=True)
    x,y = transformer.transform(lon, lat)
    
    return x, y

def stn_visualize(stn_id, stn_list, event_to_eval):
    # grab station id info and reproject coords
    stn = stn_list.loc[stn_list['era-id'] == stn_id]
    lon, lat = stn.longitude.values[0], stn.latitude.values[0]
    x,y = latlon_to_mercator_cartopy(lat, lon)

    # figure set-up
    fig, ax = plt.subplots(subplot_kw={'projection':ccrs.epsg(3857)})
    ax.coastlines()
    ax.add_feature(cf.BORDERS)
    ax.add_feature(cf.STATES, lw=0.5)

    ax.set_extent([lon+1, lon-1, lat-1, lat+1])    

    # Obtain the limits of the plot
    x0,x1,y0,y1 = ax.get_extent()

    # Create a polygon with the limits of the plot
    polygon = Polygon(((x0,y0),(x0,y1),(x1,y1),(x1,y0)))

    # Use only the counties that overlap with the actual plot
    counties = ca_county[ca_county.overlaps(polygon)]

    # Plot the counties' geometries
    for geometry in counties.geometry:
        ax.add_geometries(geometry.boundary, crs=ax.projection, 
                          facecolor='none', edgecolor='teal',
                          lw=0.5, 
                         )
    
    ax.xaxis.set_major_locator(MaxNLocator(integer=True))
    ax.yaxis.set_major_locator(MaxNLocator(integer=True))

    ax.plot(lon, lat, 'ok', markersize=8, transform=ccrs.PlateCarree(), mfc='none')
    ax.plot(x, y, '.r', markersize=4)
    ax.annotate('{}'.format(stn_id), xy=(x,y), xytext=(x+10, y+10), fontsize=6) # station name
    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=["bottom", "left"],
                    ls=":", lw=0.5)
    ax.set_title("{} evaluation \nat {}".format(event_to_eval, stn_id))

def event_plot(df, var, event, alt_start_date=None, alt_end_date=None, dpi=None):
    '''Produces timeseries of variables that have flags placed'''
    
    fig, ax = plt.subplots(figsize=(10,3))

    # plot all observations
    df.plot(ax=ax, x='time', y=var, marker=".", ms=4, lw=1, 
        color="k", alpha=0.5, label='Cleaned data')

    # plot event timeline 
    event_start, event_end = event_info(event, alt_start_date, alt_end_date)
    ax.axvspan(event_start, event_end, color='red', alpha=0.1, label='{}'.format(event))

    # ax.axhline(event_start, color='red', lw=2, alpha=0.25)
    # ax.axhline(event_end, color='red', lw=2, alpha=0.25)
    # ax.fill_between(x='time', 0, 1, where=y)

    # plot any flags placed by QA/QC
    if len(df[var+'_eraqc'].dropna().unique()) != 0: 
        
        # identify flagged data, can handle multiple flags
        for flag in df[var+'_eraqc'].dropna().unique():
            flag_name = id_flag(flag)
            flag_label = "{:.3f}% of data flagged by {}".format(
                100*len(df.loc[df[var+'_eraqc'] == flag, var])/len(df), 
                flag_name)

            flagged_data = df[~df[var+'_eraqc'].isna()]
            flagged_data.plot(x="time", y=var, ax=ax, 
                              marker="o", ms=7, lw=0, 
                              mfc="none", color="C3",
                              label=flag_label)

    legend = ax.legend(loc='upper left', prop={'size': 8})    

    # plot aesthetics
    ylab, units, miny, maxy = _plot_format_helper(var)
    plt.ylabel('{} [{}]'.format(ylab, units));
    plt.xlabel('')
    plt.title('QA/QC event evaluation: {}: {}'.format(event, df['station'].unique()[0]), fontsize=10)