'''Plotting functions used across QA/QC evaluation protocol for Historical Data Platform'''

from pyproj import CRS, Transformer
import geopandas as gpd
from geopandas import GeoDataFrame
from shapely.geometry import Point

import matplotlib.pyplot as plt
import cartopy.feature as cf
from matplotlib.ticker import MaxNLocator
import cartopy.crs as ccrs

from qaqc_eval_utils import event_info

import sys
import os
sys.path.append(os.path.expanduser('../'))
from qaqc_plot import flagged_timeseries_plot, _plot_format_helper

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
    # ax.plot(ca_county)
    ax.set_extent([lon+1, lon-1, lat-1, lat+1])    
    ax.xaxis.set_major_locator(MaxNLocator(integer=True))
    ax.yaxis.set_major_locator(MaxNLocator(integer=True))

    ax.plot(lon, lat, 'ok', markersize=8, transform=ccrs.PlateCarree(), mfc='none')
    ax.plot(x, y, '.r', markersize=4)
    ax.annotate('{}'.format(stn_id), xy=(x,y), xytext=(x+10, y+10), fontsize=6) # station name
    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=["bottom", "left"],
                    ls=":", lw=0.5)
    ax.set_title("{} evaluation \nat {}".format(event_to_eval, stn_id))


def event_plot(df, var, event, dpi=None):
    '''Produces timeseries of variables that have flags placed'''

    fig, ax = plt.subplots(figsize=(10,3))

    # plot all observations
    df.plot(ax=ax, x='time', y=var, marker=".", ms=4, lw=1, 
        color="k", alpha=0.5, label='Cleaned data')

    # plot event timeline 
    event_start, event_end = event_info(event)
    ax.axvspan(event_start, event_end, color='red', alpha=0.1, label='{}'.format(event))

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