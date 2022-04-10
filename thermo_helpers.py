
def preprocessing(ds):
    ds = ds.expand_dims('time')   
    return ds


def c_to_f(temperature):
    return ((temperature * (9/5)) + 32.)


def f_to_c(temperature):
    return ((temperature - 32.) * (5/9))


def dewpoint_from_relative_humidity(temperature, relative_humidity):
    import numpy as np
    r"""Calculate the ambient dewpoint given air temperature and relative humidity.
    Parameters
    ----------
    temperature : `pint.Quantity`
        Air temperature
    relative_humidity : `pint.Quantity`
        Relative humidity expressed as a ratio in the range 0 < relative_humidity <= 1
    Returns
    -------
    `pint.Quantity`
        Dewpoint temperature
    See Also
    --------
    dewpoint, saturation_vapor_pressure
    """
    if np.any(relative_humidity > 1.0):
        warnings.warn('Relative humidity >120%, ensure proper units.')
    return dewpoint(relative_humidity * saturation_vapor_pressure(temperature))


def dewpoint(vapor_pressure):
    import numpy as np
    r"""Calculate the ambient dewpoint given the vapor pressure.
    Parameters
    ----------
    vapor_pressure : `pint.Quantity`
        Water vapor partial pressure
    Returns
    -------
    `pint.Quantity`
        Dewpoint temperature
    See Also
    --------
    dewpoint_from_relative_humidity, saturation_vapor_pressure, vapor_pressure
    Notes
    -----
    This function inverts the [Bolton1980]_ formula for saturation vapor
    pressure to instead calculate the temperature. This yields the following formula for
    dewpoint in degrees Celsius, where :math:`e` is the ambient vapor pressure in millibars:
    .. math:: T = \frac{243.5 \log(e / 6.112)}{17.67 - \log(e / 6.112)}

    """
    val = np.log(vapor_pressure / 6.112)
    return ((243.5 * val) 
                    / (17.67 - val))
    

def saturation_vapor_pressure(temperature):
    import numpy as np
    r"""Calculate the saturation water vapor (partial) pressure.
    Parameters
    ----------
    temperature : `pint.Quantity`
        Air temperature
    Returns
    -------
    `pint.Quantity`
        Saturation water vapor (partial) pressure
    See Also
    --------
    vapor_pressure, dewpoint
    Notes
    -----
    Instead of temperature, dewpoint may be used in order to calculate
    the actual (ambient) water vapor (partial) pressure.
    The formula used is that from [Bolton1980]_ for T in degrees Celsius:
    .. math:: 6.112 e^\frac{17.67T}{T + 243.5}
    """
    return (6.112 * np.exp((17.67 * temperature) 
                                    / (temperature + 243.5)))

def relative_humidity_from_dewpoint(temperature,dewpoint):
    temperature = f_to_c(temperature)
    dewpoint = f_to_c(dewpoint) 
    
    e = saturation_vapor_pressure(dewpoint)
    e_s = saturation_vapor_pressure(temperature)
    return (e / e_s)

# def add_time_dim(xda):
#     xda = xda.expand_dims(time = file_list.index(xda))

# #     import datetime as datetime
# #     xda = xda.expand_dims(time = xda.index)
#     return xda


# def read_netcdfs(files, dim, transform_func=None):
#     # usage 
#     # combined = read_netcdfs('/all/my/files/*.nc', dim='time')
#     from glob import glob
#     def process_one_path(path):
#         # use a context manager, to ensure the file gets closed after use
#         with xr.open_dataset(path) as ds:
#             # transform_func should do some sort of selection or
#             # aggregation
#             if transform_func is not None:
#                 ds = transform_func(ds)
#             # load all data from the transformed dataset, to ensure we can
#             # use it after closing each original file
#             ds.load()
#             return ds

#     paths = sorted(glob(files))
#     datasets = [process_one_path(p) for p in paths]
#     combined = xr.concat(datasets, dim)
#     return combined