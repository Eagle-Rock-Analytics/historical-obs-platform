"""
opens hadISD netcdfs (after expanding via hadisd-decompress-time.py),
then slightly cleans the file.
(1) linearly interpolates 3-hourly and smaller gaps in temp and dpt.
(2) calculates RH.
(3) converts temps from C > F
(4) separates variables and their QC data into individual files -
I find it cleaner to work with these timeseries rather than the large
arrays originally available.
"""

import pandas as pd
import xarray as xr
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

xr.set_options(keep_attrs=True)
import thermo_helpers as th
import datetime
import os

## some definitions
# station identifiers from NOAA
wecc_df = pd.read_csv("wecc-station-data.csv", header=0).drop(["Unnamed: 0"], axis=1)

# order of vars in reporting_stats
reporting_vars = [
    "temperatures",
    "dewpoints",
    "slp",
    "stnlp",
    "windspeeds",
    "winddirs",
    "total_cloud_cover",
    "low_cloud_cover",
    "mid_cloud_cover",
    "high_cloud_cover",
    "precip1_depth",
    "precip2_depth",
    "precip3_depth",
    "precip6_depth",
    "precip9_depth",
    "precip12_depth",
    "precip15_depth",
    "precip18_depth",
    "precip24_depth",
    "relative_humidity",
]

# look-up dictionaries for QC flag columns for different parts of this test
flag_col_dict = {
    "temperatures": [0, 1, 4, 5, 8, 12, 16, 20, 24, 27, 41, 44, 54, 58, 71],
    "original_temperatures": [0, 1, 4, 5, 8, 12, 16, 20, 24, 27, 41, 44, 54, 58],
    "dewpoints": [0, 2, 4, 6, 9, 13, 17, 21, 25, 28, 30, 31, 32, 42, 45, 55, 59, 71],
    "original_dewpoints": [
        0,
        2,
        4,
        6,
        9,
        13,
        17,
        21,
        25,
        28,
        30,
        31,
        32,
        42,
        45,
        55,
        59,
    ],
    "slp": [0, 3, 4, 7, 11, 15, 19, 23, 26, 29, 43, 46, 57, 60],  # 26 should be empty
    "windspeeds": [0, 4, 10, 14, 18, 22, 47, 56, 61, 62, 63, 64, 65],
    "winddirs": [0, 4, 10, 14, 18, 22, 47, 48, 56, 61, 62, 63, 64, 65, 66, 67, 68],
    "clouds": [33, 34, 35, 36, 37, 38, 39, 40],  # for completeness, but unused
    "total_cloud_cover": [33, 37, 40],
    "low_cloud_cover": [34, 38, 40],
    "mid_cloud_cover": [35, 38, 39, 40],
    "high_cloud_cover": [36, 38, 39, 40],
    "stnlp": [69],
    "precip1_depth": [70],
    "precip2_depth": [70],
    "precip3_depth": [70],
    "precip6_depth": [70],
    "precip9_depth": [70],
    "precip12_depth": [70],
    "precip15_depth": [70],
    "precip18_depth": [70],
    "precip24_depth": [70],
    "relative_humidity": [
        0,
        1,
        2,
        4,
        5,
        6,
        8,
        9,
        12,
        13,
        16,
        17,
        20,
        21,
        24,
        25,
        27,
        28,
        30,
        31,
        32,
        41,
        42,
        44,
        45,
        54,
        55,
        58,
        59,
        71,
    ],
}  # t and dpt combined

qc_expanded_strings = [
    "Duplicate value",
    "T - Frequent value",
    "Td - Frequent value",
    "Frequent value",
    "Diurnal cycle",
    "T - Gap",
    "Td - Gap",
    "Gap",
    "T - Record",
    "Td - Record",
    "Record",
    "Record",
    "T - Straight string",
    "Td - Straight string",
    "Straight string",
    "Straight string",
    "T - Hour string",
    "Td - Hour string",
    "Hour string",
    "Hour string",
    "T - Day string",
    "Td - Day string",
    "Day string",
    "Day string",
    "T- Climatological outlier",
    "Td - Climatological outlier",
    "Climatological outlier",
    "T - Spike",
    "Td - Spike",
    "Spike",
    "Supersaturation",
    "Dewpoint depression",
    "Dewpoint cutoff",
    "Unobservable",
    "Unobservable",
    "Unobservable",
    "Unobservable",
    "Small total",
    "Full",
    "Full",
    "Negative value",
    "T - Neighbor outlier",
    "Td - Neighbor outlier",
    "Neighbor outlier",
    "T - Month cleanup",
    "Td - Month cleanup",
    "Month cleanup",
    "Month cleanup",
    "Dir - Month cleanup",
    "Month cleanup",
    "Month cleanup",
    "Month cleanup",
    "Month cleanup",
    "Month cleanup",
    "T - Odd cluster",
    "Td - Odd cluster",
    "Odd cluster",
    "Odd cluster",
    "T - Variance",
    "Td - Variance",
    "Variance",
    "Variance",
    "Speed logic",
    "Direction logic",
    "Wind rose",
    "Spike",
    "Dir - Straight string",
    "Dir - Hour string",
    "Dir - Day string",
    "Station level pressure",
    "Precipitation checks",
    "Post-interpolation supersaturation",
]

######## loop through WECC stations
for s in range(len(wecc_df)):
    stn_id = wecc_df["station id"][s]
    # station name is for metadata
    icao_i = wecc_df["icao"][s]
    stn_in = wecc_df["name"][s]
    # filename has ID and ICAO
    f_id = stn_id + "_" + icao_i
    # save directory defined below
    sdir = "./cleaned-wecc-hadISD/" + icao_i
    # time stamp for metadata later
    ts = datetime.datetime.now().strftime("%Y-%m-%d %H:%M")

    try:
        os.mkdir(sdir)
    except:
        pass

    fdir = "expanded-wecc/"
    bname = "hadisd.3.3.0.202202p_19310101-20220301_" + f_id
    fname = fdir + bname + ".nc"
    ds = xr.open_dataset(fname, decode_times=False, mask_and_scale=True)
    # will refill nans with fill values at the end
    ds = ds.squeeze()
    ds = ds.drop("reporting_v")  # don't want this
    ds = ds.drop("station_id")  # don't want this
    orig_att = ds.attrs  # save those attributes for later

    # remove unused "flagged_value" att

    ############ make individual observation and quality control data
    # files for each variable ######################################

    for var in list(ds.data_vars):
        try:
            ds[var].attrs.pop("flagged_value")  # remove this att
        except:
            continue

    ts = datetime.datetime.now().strftime("%Y-%m-%d %H:%M")

    ##### we need to add and interpolate some variables #####

    # initialize arrays to keep track of t, dpt, and rh
    # measurements which get QC'd out in post due to supersaturation
    arr_sss = np.empty(len(ds["time"]))
    arr_sss[:] = 0.0

    # convert t and dpt to fahrenheit
    # the rh function will convert back to c as necessary

    ds["temperatures"] = th.c_to_f(ds["temperatures"])
    ds["temperatures"].attrs["units"] = "degrees Fahrenheit"
    min_t = ds["temperatures"].attrs["valid_min"]
    ds["temperatures"].attrs["valid_min"] = th.c_to_f(min_t)
    max_t = ds["temperatures"].attrs["valid_max"]
    ds["temperatures"].attrs["valid_max"] = th.c_to_f(max_t)

    ds["dewpoints"] = th.c_to_f(ds["dewpoints"])
    ds["dewpoints"].attrs["units"] = "degrees Fahrenheit"
    min_d = ds["dewpoints"].attrs["valid_min"]
    ds["dewpoints"].attrs["valid_min"] = th.c_to_f(min_d)
    max_d = ds["dewpoints"].attrs["valid_max"]
    ds["dewpoints"].attrs["valid_max"] = th.c_to_f(max_d)

    # interpolate 3-hour gaps in temperature and dewpoint
    # after saving the original variables for the record
    ds["original_temperatures"] = ds["temperatures"]
    ds["original_dewpoints"] = ds["dewpoints"]
    ds["temperatures"] = ds["temperatures"].interpolate_na(
        dim="time", method="linear", limit=3, max_gap=3
    )
    ds["dewpoints"] = ds["dewpoints"].interpolate_na(
        dim="time", method="linear", limit=3, max_gap=3
    )

    # take note on where the above leads to supersaturation
    # make a new QC variable, and put the flagged values
    # into a new variable as well
    supersaturated = np.where(ds["temperatures"].values - ds["dewpoints"].values < 0.0)
    arr_sss[supersaturated] = 1.0
    post_sss = xr.DataArray(data=arr_sss, coords={"time": ds["time"]})
    ds["interp_sss"] = post_sss

    # reuse the array to save some memory
    # save flagged t values
    arr_sss[:] = np.nan
    arr_sss[supersaturated] = ds["temperatures"].values[supersaturated]
    post_sss = xr.DataArray(data=arr_sss, coords={"time": ds["time"]})
    ds["temperatures_sss"] = post_sss
    ds["flagged_obs"].loc[dict(flagged=0, time=ds["time"][supersaturated])] = ds[
        "temperatures_sss"
    ][supersaturated]
    ds = ds.drop("temperatures_sss")

    # save flagged dpt values
    arr_sss[:] = np.nan
    arr_sss[supersaturated] = ds["dewpoints"].values[supersaturated]
    post_sss = xr.DataArray(data=arr_sss, coords={"time": ds["time"]})
    ds["dewpoints_sss"] = post_sss
    ds["flagged_obs"].loc[dict(flagged=1, time=ds["time"][supersaturated])] = ds[
        "dewpoints_sss"
    ][supersaturated]
    ds = ds.drop("dewpoints_sss")

    # now compute rh with the flagged t and dpt values
    # because we want to save those values just in case
    rh = th.relative_humidity_from_dewpoint(
        ds["temperatures"].values, ds["dewpoints"].values
    )
    relative_humidity = xr.DataArray(data=rh, coords={"time": ds["time"]})
    ds["relative_humidity"] = relative_humidity
    ds["relative_humidity"].attrs = {
        "long_name": "Relative humidity at screen height (~2m)",
        "units": "fraction",
        "note": "derived in post-processing",
    }

    # and save the offending rh values in a flagged variable
    arr_sss[:] = np.nan
    arr_sss[supersaturated] = ds["relative_humidity"].values[supersaturated]
    post_sss = xr.DataArray(data=arr_sss, coords={"time": ds["time"]})
    ds["flagged_rh"] = post_sss

    # now make a brand new qc flag variable with the ones we just added
    updated_qc = xr.concat(
        [ds["quality_control_flags"], ds["interp_sss"]], dim="test", coords="all"
    )
    ds = ds.drop("quality_control_flags")
    ds["quality_control_flags"] = updated_qc
    ds = ds.drop("interp_sss")

    # now save the reasonable t, dpt, and rh variables
    ds["relative_humidity"] = ds["relative_humidity"].where(
        ds["temperatures"] - ds["dewpoints"] > 0.0
    )
    ds["temperatures"] = ds["temperatures"].where(
        ds["temperatures"] - ds["dewpoints"] > 0.0
    )
    ds["dewpoints"] = ds["dewpoints"].where(ds["temperatures"] - ds["dewpoints"] > 0.0)

    # get QC data into usable form
    qualList = ["quality_control_flags"]
    qual_ds = ds[qualList].fillna(0.0)
    qual_ds = qual_ds.reset_coords(
        names=["longitude", "latitude", "elevation"], drop=True
    )

    # qual_dsv = qual_ds.isel(test=flag_col_dict[v])
    qual_dff = qual_ds.to_dataframe()
    qual_dff.columns = qual_dff.columns.to_flat_index()

    # a pandas dataframe makes some things easier
    qual_dff.columns = qual_dff.columns.get_level_values(0)
    qual_dff.reset_index(inplace=True)
    qual_dff = qual_dff.pivot(index="time", columns="test")
    qual_dff.columns = qual_dff.columns.to_flat_index()
    colNames = qual_dff.columns

    ####### now go through the variables
    for i in range(len(reporting_vars)):
        v = reporting_vars[i]
        sname = sdir + "/" + bname + "_" + v + ".nc"

        if "cloud" in sname:
            continue
        #     print(sname)
        q_to_keep_list = flag_col_dict.get(v)
        qual_df = qual_dff.iloc[:, q_to_keep_list]
        qual_df = qual_df[qual_df != 0.0]

        ############################## define special cases

        # t and dpt are special cases
        # since I did some interpolating to fill gaps
        if "temperatures" in sname or "dewpoints" in sname:
            extn = True
            rht = False
            pp = False

        # rh is a special case:
        # derived from t and td and thus its QC is dependent on two
        # sets of QC codes. There will be lots of overlap!
        elif "relative_humidity" in sname:
            rht = True
            extn = False
            pp = False

        # precip is another special case because I will put
        # all the accumulation times into one file
        elif "precip" in sname:
            pp = True
            rht = False
            extn = False

        else:
            pp = False
            rht = False
            extn = False

        q_strings = [qc_expanded_strings[c] for c in q_to_keep_list]
        q_ids = np.arange(0, len(q_strings), 1)
        qual_df.columns = q_strings

        result = {}
        resultLen = {}
        total_thrown = 0

        for q in q_strings:
            to_throw = np.isfinite(qual_df[q])
            l = qual_df.index[to_throw]
            to_throw = len(list(l))
            total_thrown = total_thrown + to_throw
            resultLen[str(q)] = to_throw

        # qual_da = qual_ds['quality_control_flags']
        qual_da = qual_ds.isel(test=flag_col_dict[v])
        t_breakdown = [
            k + ": " + str(n) for k, n in list(zip(q_strings, resultLen.values()))
        ]
        t_breakdown = ", \n".join(t_breakdown)
        qual_da.attrs["total_QC_instances"] = str(total_thrown)
        qual_da.attrs["removed_obs_breakdown"] = t_breakdown

        if not pp:
            qual_da.attrs["note"] = "QC reasons may overlap"

        if extn:
            recovered_obs = len(np.where(~np.isnan(ds[v].values))[0]) - len(
                np.where(~np.isnan(ds["original_" + v].values))[0]
            )
            qual_da.attrs["note2"] = "Last QC code does not apply to original values"

        ########################### now get the not reported stats for the variable
        # special cases

        if rht:
            no_report = ds["no_report"][:, 0:2]
            no_report.attrs[
                "long_name"
            ] = "Temperature and/or dewpoint not reported by station"
            no_report.attrs["reporting_v_order"] = "temperatures, dewpoints"
            no_report_t = str(int(sum(ds["no_report"][:, 0].values)))
            no_report.attrs["temperature_obs_not_reported"] = no_report_t
            no_report_d = str(int(sum(ds["no_report"][:, 1].values)))
            no_report.attrs["dewpoint_obs_not_reported"] = no_report_d
            fl_obs = ds["flagged_rh"]
            temp_rs = ds["reporting_stats"][0:2]
            temp_rs.attrs[
                "long_name"
            ] = "Temperature and dewpoint reporting frequency and accuracy for each month"
            temp_rs.attrs["reporting_v_order"] = "temperatures, dewpoints"

        else:
            no_report = ds["no_report"][:, i]
            no_report.attrs["long_name"] = "Not reported by station"
            no_report_total = str(int(sum(ds["no_report"][:, i].values)))
            no_report.attrs["obs_not_reported"] = no_report_total
            fl_obs = ds["flagged_obs"][:, i]
            temp_rs = ds["reporting_stats"][i]
            temp_rs.attrs[
                "long_name"
            ] = "Reporting frequency and accuracy for each month"

        fl_obs.attrs["long_name"] = "Observation values removed by QC flags"
        fl_obs.attrs["units"] = ds[v].attrs["units"]

        if extn:
            merged_ds = xr.merge(
                [ds[v], ds["original_" + v], qual_da, fl_obs, temp_rs, no_report]
            )
            merged_ds[v].attrs[
                "post_processing_applied"
            ] = "3-hour gaps filled via linear interpolation"
            merged_ds["original_" + v].attrs["note"] = "original output by HadISD"
            merged_ds[v].attrs["obs_inserted_by_linear_interpolation"] = recovered_obs
            merged_ds["flagged_obs"].attrs[
                "note"
            ] = "Includes removal after interpolation"

        else:
            merged_ds = xr.merge([ds[v], qual_da, fl_obs, temp_rs, no_report])

        merged_ds["quality_control_flags"].attrs["total_QC_instances"] = str(
            total_thrown
        )
        merged_ds["quality_control_flags"].attrs["removed_obs_breakdown"] = t_breakdown
        merged_ds["quality_control_flags"].attrs["note"] = "QC reasons may overlap"

        merged_atts = {
            "data_reformatting_by": "Beth McClenny, Eagle Rock Analytics",
            "history": ds.attrs["history"] + "Data reformatted " + ts,
            "keywords": "sub-daily, station, extremes, " + v,
            "summary": "Quality-controlled, sub-daily, station dataset containing " + v,
            "station_longname": stn_in,
            "station_ICAO": icao_i,
        }
        to_add = {**orig_att, **merged_atts}
        merged_ds.attrs = to_add
        merged_ds = merged_ds.fillna(-1e30)

        for dv in list(merged_ds.data_vars):
            merged_ds[dv].encoding["_FillValue"] = -1e30

        merged_ds.longitude.encoding["_FillValue"] = None
        merged_ds.latitude.encoding["_FillValue"] = None
        merged_ds.time.encoding["_FillValue"] = None
        merged_ds.elevation.encoding["_FillValue"] = None
        if len(merged_ds["time"]) - int(no_report_total) > 0:
            print("saving " + sname)
            merged_ds.to_netcdf(sname)

    ############ last make files with (1) station ID time series (to see if/
    # when data has come from a merged station) and (2) no report time series
    # (ie, the station was not reporting anything at all. ##################

    the_vars = ["input_station_id", "no_obs"]

    for i in range(len(the_vars)):
        v = the_vars[i]
        sname = sdir + "/" + bname + "_" + v + ".nc"
        merged_ds = ds[v]

        # metadata
        merged_atts = {
            "data_reformatting_by": "Beth McClenny, Eagle Rock Analytics",
            "history": ds.attrs["history"] + "Data reformatted " + ts,
            "keywords": "sub-daily, station, extremes, " + v,
            "summary": "Quality-controlled, sub-daily, station dataset containing " + v,
            "station_longname": stn_in,
            "station_ICAO": icao_i,
        }
        to_add = {**orig_att, **merged_atts}
        merged_ds.attrs = to_add
        merged_ds = merged_ds.fillna(-1e30)

        if "station_id" in sname:
            merged_ds.attrs["description"] = "input station id"
            merged_ds.encoding["dtype"] = "S1"

        else:
            merged_ds.attrs[
                "description"
            ] = "time series showing when station was not reporting any data"
            merged_ds.encoding["_FillValue"] = -1e30

        # some xarray encoding stuff
        merged_ds.longitude.encoding["_FillValue"] = None
        merged_ds.latitude.encoding["_FillValue"] = None
        merged_ds.time.encoding["_FillValue"] = None
        merged_ds.elevation.encoding["_FillValue"] = None
        print("saving " + sname)
        merged_ds.to_netcdf(sname)
