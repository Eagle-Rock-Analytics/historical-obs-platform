"""
HadISD data features a compressed time axis to save server space. 
This time axis excludes times which had no station reporting.
The following code adds these timesteps back in so the time spacing
is consistent.
"""

import pandas as pd
import xarray as xr

xr.set_options(keep_attrs=True)
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

# need WECC station metadata
wecc_df = pd.read_csv("wecc-station-data.csv", header=0).drop(["Unnamed: 0"], axis=1)

# order of vars in reporting_stats
# needed to flag where there are missing datapoints for each var
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
]

# open original hadISD netCDFs, expand time axis,
# count where all data or var-specific is missing
for s in range(len(wecc_df)):

    stn_id = wecc_df["station id"][s]
    icao_i = wecc_df["icao"][s]
    f_id = stn_id + "_" + icao_i

    fname = (
        "wecc-hadisd3.3.0.202202p/hadisd.3.3.0.202202p_19310101-20220301_"
        + f_id
        + ".nc"
    )
    sname = "expanded-wecc/hadisd.3.3.0.202202p_19310101-20220301_" + f_id + ".nc"
    ds = xr.open_dataset(fname, decode_times=True, mask_and_scale=False)
    the_vars = list(ds.data_vars)
    ds_infilled = ds["time"].resample(time="1H").asfreq()  # resamp to hourly
    ds_infilled = ds.reindex({"time": ds_infilled["time"]})  # reindex time

    # initialize array to flag missing obs per variable
    missing_obs = np.empty((len(ds_infilled["time"]), len(reporting_vars)))
    missing_obs[:, :] = 0.0
    # initialize array to flag where all obs are missing
    arr = np.empty(len(ds_infilled["time"]))
    arr[:] = 0.0

    for i in range(len(reporting_vars)):
        var = reporting_vars[i]
        missing_obs[
            np.where(
                ds_infilled[var].values == ds_infilled[var].attrs["missing_value"]
            ),
            i,
        ] = 1.0
        missing_obs[np.isnan(ds_infilled[var].values)] = 1.0
    no_report = xr.DataArray(
        data=missing_obs,
        coords={"time": ds_infilled["time"], "reporting_v": ds_infilled["reporting_v"]},
    )
    # above array is time x var - flags where a variable is missing

    for i in range(len(the_vars)):
        var = the_vars[i]
        try:
            tocheck = ds[var].values[0]

            if isinstance((tocheck), (int, np.integer)):
                ds_infilled[var] = (
                    ds_infilled[var]
                    .where(ds_infilled[var] > int(-887))
                    .fillna(int(-999))
                    .astype(int)
                )
            elif isinstance((tocheck), (float)):
                ds_infilled[var] = (
                    ds_infilled[var]
                    .where(ds_infilled[var] > -1e30)
                    .fillna(-1e30)
                    .astype(float)
                )
            elif isinstance((tocheck), (np.bytes_)):
                ds_infilled[var] = ds_infilled[var].fillna("no report").astype(str)
        except:
            print(var)

    arr[np.where(np.isnan(ds_infilled["quality_control_flags"].values[:, 0]))] = 1.0
    no_obs = xr.DataArray(data=arr, coords={"time": ds_infilled["time"]})
    # above vector is length of time - flags where station did not report at all
    ds_infilled["no_obs"] = no_obs
    ds_infilled["no_report"] = no_report

    # save new expanded netCDFs
    ds_infilled.to_netcdf(sname)
    ds_infilled = None
