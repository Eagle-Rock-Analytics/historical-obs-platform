# Scratch from CWOP_Clean for stage 2

elif np.isnan(ds['pressure_set_1'].values).all(): # If all direct station observations are NA
                        if 'altimeter_set_1' in ds.keys(): 
                            if not np.isnan(ds['altimeter_set_1'].values).all(): # If altimeter readings observed.
                                # Calculate pressure from altimeter and elevation
                                ds['ps'].values = calc_clean._calc_ps_alt(ds['altimeter_set_1'], ds['elevation']) 
                                ds['ps'].attrs['comment'] = "Calculated from altimeter setting, station elevation using calc_clean.py."
                                ds['ps'].attrs['ancillary_variables'] = "altimeter_set_1 elevation" # List other variables associated with variable (QA/QC)
                                #print("Calculated ps from altimeter.")

                            else: # If altimeter readings all NA
                                if not np.isnan(ds['sea_level_pressure_set_1'].values).all():
                                    ds['ps'].values = calc_clean._calc_ps(ds['sea_level_pressure_set_1'], ds['elevation'], ds['tas']) 
                                    ds['ps'].attrs['comment'] = "Calculated from sea level pressure, station elevation, air temperature using calc_clean.py."
                                    ds = ds.rename({'sea_level_pressure_set_1_qc': 'psl'})
                                    ds['ps'].attrs['ancillary_variables'] = "psl elevation tas" # List other variables associated with variable (QA/QC)
                                
                        else: # If no altimeter readings
                            if not np.isnan(ds['sea_level_pressure_set_1'].values).all():
                                ds['ps'].values = calc_clean._calc_ps(ds['sea_level_pressure_set_1'], ds['elevation'], ds['tas']) 
                                ds['ps'].attrs['comment'] = "Calculated from sea level pressure, station elevation, air temperature using calc_clean.py."
                                ds = ds.rename({'sea_level_pressure_set_1_qc': 'psl'})
                                ds['ps'].attrs['ancillary_variables'] = "psl elevation tas" # List other variables associated with variable (QA/QC)
                            
                else: # If no pressure column exists.
                    if 'altimeter_set_1' in ds.keys(): 
                        if not np.isnan(ds['altimeter_set_1'].values).all(): # If altimeter readings observed.
                                
                                # Calculate pressure from altimeter and elevation
                                ds['ps'].values = calc_clean._calc_ps_alt(ds['altimeter_set_1'], ds['elevation']) 
                                ds['ps'].attrs['comment'] = "Calculated from altimeter setting, station elevation using calc_clean.py."
                                ds['ps'].attrs['ancillary_variables'] = "altimeter_set_1 elevation" # List other variables associated with variable (QA/QC)

                        elif np.isnan(ds['altimeter_set_1'].values).all(): # If altimeter readings not observed:
                            if not np.isnan(ds['sea_level_pressure_set_1'].values).all():
                                ds['ps'].values = calc_clean._calc_ps(ds['sea_level_pressure_set_1'], ds['elevation'], ds['tas']) 
                                ds['ps'].attrs['comment'] = "Calculated from sea level pressure, station elevation, air temperature using calc_clean.py."
                                ds = ds.rename({'sea_level_pressure_set_1_qc': 'psl'})
                                ds['ps'].attrs['ancillary_variables'] = "psl elevation tas" # List other variables associated with variable (QA/QC)
                    


# If vapor pressure doesn't exist but tas and hurs do, calculate using opt 1.
                elif 'tas' in ds.keys() and 'relative_humidity_set_1' in ds.keys():
                    # Note we don't check for non-NA values here because there's no alternative calculation method.
                    # If calculation returns all NAs, this is equivalent to manually setting the column to be NA.
                    # Inputs: tas (K) and relative humidity (%)
                    ds['tdps'] = calc_clean._calc_dewpointtemp_opt1(ds['tas'], ds['relative_humidity_set_1'])

                    # Set attributes for calculation.
                    ds['tdps'].attrs['ancillary_variables'] = "tas hurs" # List other variables associated with variable (QA/QC)
                    ds['tdps'].attrs['comment'] = "Calculated from air temperature and relative humidity using calc_clean.py."
                
                # Vapor pressure not collected by CWOP, so can't calculate tdps using opt2.
                # If none of these variables available, make NA column.
                else:
                    tdps = np.full((1, len(ds['time'])), np.nan)
                    ds['tdps'] = (['station', 'time'], tdps)
                
                # Set attributes
                ds['tdps'].attrs['long_name'] = "dew_point_temperature"
                ds['tdps'].attrs['standard_name'] = "dew_point_temperature"
                ds['tdps'].attrs['units'] = "degree_Kelvin"
                