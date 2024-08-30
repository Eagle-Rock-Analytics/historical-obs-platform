'''Functions used across QA/QC evaluation protocol for Historical Data Platform'''

from pyproj import CRS, Transformer
import geopandas as gpd
from geopandas import GeoDataFrame
from shapely.geometry import Point



def known_issue_check(network, var, stn):
    '''
    Identifies if station under evaluation has a known network issue.
    At present, only prints out a statement if there is an issue.
    Eventually may want to do <something>

    Note: See "Known Network Issues for QA/QC Validation" planning doc.
    '''

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


def subset_eval_stns(event_to_eval, stn_list, subset=None):
    '''
    Identifies stations to evaluate for specific V1 QA/QC events.
    Option to subset to a more manageable number of random stations for initial evaluation. 
    '''

    # TO DO: validation check on event_to_eval options

    event_flags = []
    event_flags.append('all')
    event_flags.append(event_to_eval) # options: santa_ana_wind, winter_storm, AR, mudslide, heatwave1, heatwave2, heatwave3, offshore_wind

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
    ## Need an option for "WECC wide"
    if event_to_eval == 'santa_ana_wind':
        counties_to_grab = ['Los Angeles', 'Orange']
    elif event_to_eval == 'winter_storm':
        counties_to_grab = []
    elif event_to_eval == 'mudslide':
        counties_to_grab = ['Santa Barbara']
    elif event_to_eval == 'AR':
        counties_to_grab = []
    elif event_to_eval == 'heatwave1': # August 2020
        counties_to_grab = []
    elif event_to_eval == 'heatwave2': # September 2020
        counties_to_grab = ['Los Angeles']
    elif event_to_eval == 'heatwave3': # August 2022
        counties_to_grab = []
    elif event_to_eval == 'offshore_wind':
        counties_to_grab = []

    target_counties = ca_county[ca_county['NAME'].isin(counties_to_grab)]
    target_counties = GeoDataFrame(target_counties, geometry=target_counties.geometry)

    geometry = [Point(latlon_to_mercator_cartopy(lat,lon)) for lat,lon in zip (event_stns.latitude, event_stns.longitude)]
    event_stns = GeoDataFrame(event_stns,geometry=geometry).set_crs(crs="EPSG:3857", allow_override=True) # adding geometry column
    event_stns_local = gpd.overlay(event_stns, target_counties, how="intersection") # subsetting for stations within county boundaries
    print('{} potential stations available for evaluation for {} event!'.format(len(event_stns_local), event_to_eval))

    if subset != None:
        if len(event_stns_local) <= subset:
            eval_stns = event_stns_local
        else:
            eval_stns = event_stns_local.sample(subset, replace=False, random_state=1)
            print('{} stations selected for evaluation for {} event!'.format(subset, event_to_eval))
    else:
        eval_stns = event_stns_local

    # lastly, check if there are any known issues --- need to refactor for check
    # check_networks = event_stns.network.unique()
    # known_issue_check(event_stns, var='tas')

    return eval_stns


# Projection stuffs
def latlon_to_mercator_cartopy(lat, lon):

    proj_latlon = CRS('EPSG:4326')
    proj_mercator = CRS('EPSG:3857')
    
    # Transform the coordinates
    transformer = Transformer.from_crs(proj_latlon, proj_mercator, always_xy=True)
    x,y = transformer.transform(lon, lat)
    
    return x, y