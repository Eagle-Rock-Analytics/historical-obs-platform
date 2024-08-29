'''Functions used across QA/QC evaluation protocol for Historical Data Platform'''

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
