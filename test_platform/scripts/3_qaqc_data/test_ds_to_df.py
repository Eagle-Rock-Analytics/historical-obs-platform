import xarray as xr
import numpy as np
import s3fs
from QAQC_pipeline import qaqc_ds_to_df

fname = 'CIMIS_45'
s3 = s3fs.S3FileSystem(anon=False)
network = fname.split('_')[0]
s3_url = 's3://wecc-historical-wx/3_qaqc_wx_dev/{}/{}.nc'.format(network, fname)
s3_file_obj = s3.open(s3_url, mode='rb')
ds = xr.open_dataset(s3_file_obj, engine='h5netcdf')

df, MultiIndex, attrs, var_attrs, era_qc_vars = qaqc_ds_to_df(ds)
print(df.columns)
print(df.anemometer_height_m)
exit()
