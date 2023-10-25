'''
This file is used to build radiosonde xarray data set using data from
https://ruc.noaa.gov/raobs/

@author: cem
@version: 1.0 (10/10/2023)

'''

import xarray as xr
from datetime import datetime
import numpy as np

def get_radiosonde_dataset(ds, flight_levels):
    '''
    Create xarray dataset from rawinsonde data from https://ruc.noaa.gov/raobs/ in netcdf format. This function differs from
    others in that if it can't find an observation at a required level, it checks nearby levels to see if data there can be used
    as a proxy for the missing data.

    The other two functions either grab only mandatory data levels or significant levels. But this routine tries to fix the missing data.

    :param ds:
    :param flight_levels: An array of integer containing the layers in mb we want. These levels match those in ERA5 reanalysis data
    :return: radiosonde xarray at the same levels of the available model data
    '''
    # handles case where netcdf has seconds since 1970-01-01, as the units but stores time in poorly defined
    # array floating seconds since 1970-01-01
    #
    xr.set_options(keep_attrs=True)  # keep any coordinate information you have when building new xarray dataset
    time_arr = ds['synTime'].values.astype(int)
    print(time_arr[0])
    #print(time_arr[30])
    print(time_arr[-1])

    try:
        # https://numpy.org/doc/stable/reference/arrays.datetime.html
        time_convert = np.array(time_arr, dtype='datetime64[s]')
    except:
        print("Still unable to properly decode time in netcdf, program exits")

    print(time_convert[0])
    print(time_convert[-1])
    now = datetime.now()  # current date and time
    date_time = now.strftime("%m/%d/%Y, %H:%M:%S")

    # we will be using mandatory heights for most of values, except we will be adding the 125mb level if possible
    # by selecting a nearby significant wind level., but we also want

    htMan = ds['htMan'].values
    prMan = ds['prMan'].values
    wdMan = ds['wdMan'].values
    wsMan = ds['wsMan'].values
    sfc = prMan[0][0]
    print(f"Surface pressure now {sfc}")
    prSigW = ds['prSigW'].values
    htSigW = ds['htSigW'].values
    wdSigW = ds['wdSigW'].values
    wsSigW = ds['wsSigW'].values
    levelsSigW = ds['sigWLevel'].values
    # https://glossary.ametsoc.org/wiki/Mandatory_level this dataset includes surface to the 21 mandatory levels
    mandatory_levels = [sfc, 1000, 925, 850, 700, 500, 400, 300, 250, 200, 150, 100, 70, 50, 30, 20, 10, 7, 5, 3, 2, 1]

    # create new arrays but add level 125mb
    rows, cols = prMan.shape

    #flight_levels = [150, 125, 100, 70, 50, 30, 20, 10]
    cols = len(flight_levels)
    wd = np.zeros((rows,cols),dtype=np.float32)
    wd[:] = np.NAN
    ws = np.zeros((rows,cols),dtype=np.float32)
    ws[:] = np.NAN
    pres = np.zeros((rows,cols),dtype=np.float32)
    pres[:] = np.NAN
    height = np.zeros((rows,cols),dtype=np.float32)
    height[:] = np.NAN

    row_length = len(htMan)
    col_length = len(htMan[0])
    #you don't want to go before 150 or after 10mb so switch the start and end
    col_start = 10 # 150mb
    col_end = 17 #10mb

    # First, pick values from the mandatory data and then supplement with significant values.
    # We supplement 125mb because it is not a mandatory pressure; otherwise we supplement missing data; for example,
    # 10 mb may frequently be supplemented. Significant values are randomly spaced, unlike the mandatory. In the lower
    # atmosphere 5mb difference might just be a hundred feet, but in the upper-atmosphere, like around 10mb, a 5mb
    # change is approximately 1650m change in height. In contrast, around 125mb where we always have missing mandatory
    # data, 5mb is around 240m Therefore we use a should use a scaling difference. For now,
    # accept a difference up approximately 1500 feet or 457 meters for each flight level.
    #                   0     1   2    3   4   5   6   7

    mb_tolerance =   [11,  10,   8,  6,  4,  2.2,  2.1, 2.0] #associated tolerance with original order

    for i in range(row_length):
        jlevel = 0
        pres_values = prSigW[i]
        idx = (np.abs(pres_values - 125)).argmin()
        nearest_pres = pres_values[idx]
        if abs(nearest_pres- 125.0) <= mb_tolerance[1]:
            pres[i][1] = nearest_pres
            height[i][1] = htSigW[i][idx]
            wd[i][1] = wdSigW[i][idx]
            ws[i][1] = wsSigW[i][idx]
        else:
            pres[i][1] = np.nan
            height[i][1] = np.nan
            wd[i][1] = np.nan
            ws[i][1] = np.nan
            #pres_corrected_values = np.nan_to_num(prSigW[i], copy=True, nan=-99999.0, posinf=None, neginf=None)
        for j in range(col_start,col_end):
            value = mandatory_levels[j]
            index = flight_levels.index(value)
            windspd = wsMan[i][j]
            winddir =  wdMan[i][j]
            presx = prMan[i][j]
            heightx = htMan[i][j]
            wd[i][index] = winddir
            ws[i][index] = windspd
            pres[i][index] = presx
            height[i][index] = heightx
            if winddir > 999 or windspd > 999:
                idx = (np.abs(pres_values - value)).argmin()
                nearest_pres = pres_values[idx]
                mb_difference = mb_tolerance[jlevel]
                abs_mb_difference = abs(nearest_pres - value)
                if abs_mb_difference <= mb_difference:
                    pres[i][index] = nearest_pres
                    height[i][index] = htSigW[i][idx]
                    wd[i][index] = wdSigW[i][idx]
                    ws[i][index] = wsSigW[i][idx]
                else:
                    pres[i][index] = np.nan
                    height[i][index] = np.nan
                    wd[i][index] = np.nan
                    ws[i][index] = np.nan
            print(f"{i},{j}:man={mandatory_levels[j]},mb={prMan[i][j]}, m={htMan[i][j]},wd={wdMan[i][j]},ws={wsMan[i][j]}")
            jlevel = jlevel+1

    # convert meters per second wind speed to knots
    #ws = ws*1.944
    #https://www.atmos.albany.edu/facstaff/ktyle/atm533/core/week7/01_Xarray_Intro.html some info from here to build datasets
    da_hgt = xr.DataArray(height, coords=[time_convert, flight_levels],
                               dims=['time','level'], attrs={'units':'m', 'long_name':'Geopotential-Height'})
    da_pres = xr.DataArray(pres, coords=[time_convert, flight_levels],
                                dims=['time','level'],attrs={'units':'hPa','long_name':'Pressure'})
    da_wd = xr.DataArray(wd, coords=[time_convert, flight_levels],
                               dims=['time','level'], attrs={'units':'degree', 'long_name':'Wind Direction'})
    da_ws = xr.DataArray(ws, coords=[time_convert, flight_levels],
                                dims=['time','level'],attrs={'units':'m s**-1','long_name':'Wind'})

    ds_new = xr.Dataset(data_vars={'hgt':da_hgt, 'pres': da_pres, 'wd':da_wd, 'ws':da_ws}).astype('float32')
    return ds_new


def create_radiosonde_man(ds):
    '''
    Create xarray dataset from radiosonde data from https://ruc.noaa.gov/raobs/ in netcdf format.Note,
    if you want to see https://glossary.ametsoc.org/wiki/Mandatory_level for discussion of these levels.
    You can what this data looks like in text form or FSL text (see e.g.,
    https://ruc.noaa.gov/raobs/fsl_format-new.html#:~:text=(FSL%20Radiosonde%20Database),additional%20lines%20are%20data%20lines.)

    This routine only builds radiosonde data from the mandatory levels. See get_radiosonde_dataset which gets data from
    the flight levels requested
    :param ds:
    :return ds: A netcdf data containing mandatory levels of wind speed, direction, height and datetime
    '''

    # Note, there are several different types of time units in this data.
    # Here we handle case where netcdf has seconds since 1970-01-01, as the units but stores time in poorly defined
    # array floating seconds since 1970-01-01. There are two time arrays, relTime and synTime.
    # SynTime seems to better match time given when radiosonde output is in FSL format
    xr.set_options(keep_attrs=True)  # keep any coordinate information you have when building new xarray dataset

    #time_arr = ds['relTime'].values.astype(int)
    time_arr = ds['synTime'].values.astype(int)
    print(time_arr[0])
    print(time_arr[-1])

    try:
        # https://numpy.org/doc/stable/reference/arrays.datetime.html
        time_convert = np.array(time_arr, dtype='datetime64[s]')
    except:
        print("Still unable to properly decode time in netcdf, program exits")

    print(time_convert[0])
    print(time_convert[-1])
    now = datetime.now()  # current date and time
    date_time = now.strftime("%m/%d/%Y, %H:%M:%S")

    htMan = ds['htMan'].values
    prMan = ds['prMan'].values
    wdMan = ds['wdMan'].values
    wsMan = ds['wsMan'].values
    sfc = prMan[0][0]
    print(f"Surface pressure now {sfc}")
    # https://glossary.ametsoc.org/wiki/Mandatory_level this dataset includes surface to the 21 mandatory levels
    mandatory_levels = [sfc, 1000, 925, 850, 700, 500, 400, 300, 250, 200, 150, 100, 70, 50, 30, 20, 10, 7, 5, 3, 2, 1]
    #https://www.atmos.albany.edu/facstaff/ktyle/atm533/core/week7/01_Xarray_Intro.html some info from here to build datasets
    da_mhgt = xr.DataArray(htMan, coords=[time_convert, mandatory_levels],
                               dims=['time','level'], attrs={'units':'m', 'long_name':'Geopotential-Mandatory level'})
    da_mpres = xr.DataArray(prMan, coords=[time_convert, mandatory_levels],
                                dims=['time','level'],attrs={'units':'hPa','long_name':'Pressure=Mandatory level'})
    da_mwd = xr.DataArray(wdMan, coords=[time_convert, mandatory_levels],
                               dims=['time','level'], attrs={'units':'degree', 'long_name':'Wind Direction-Mandatory level'})
    da_mws = xr.DataArray(wsMan, coords=[time_convert, mandatory_levels],
                                dims=['time','level'],attrs={'units':'m s**-1','long_name':'Wind Speed-Mandatory level'})

    ds_new = xr.Dataset(data_vars={'hgt':da_mhgt, 'pres': da_mpres, 'wd':da_mwd, 'ws':da_mws}).astype('float32')
    return ds_new


# this function has been superseeded by create_radiosonde_dataset
def create_radiosonde_hires_dataset(ds, month, levels):
    '''
    Create xarray dataset from rawinsonde data from https://ruc.noaa.gov/raobs/ in netcdf format. Note, if you want to see
    what this data looks like in text form or FSL text (see e.g.,
    https://ruc.noaa.gov/raobs/fsl_format-new.html#:~:text=(FSL%20Radiosonde%20Database),additional%20lines%20are%20data%20lines.)

    This program differs from using the mandatory heights, as radiosondes have data inbetween mandatory data. For example,

      4   1500  14260   -685   -785    110     14  150mb
      6   1484  14326  99999  99999    105     17
      5   1480  14340   -694   -784  99999  99999
      6   1338  14935  99999  99999    145     21
      6   1271  15240  99999  99999    165     18
      5   1070  16242   -744   -854  99999  99999
      7   1040  16406   -759   -869    150      9
      4   1000  16640   -754   -864    145     16 100mb
      5    955  16905   -761   -871  99999  99999
      5    835  17685   -754   -884  99999  99999
      6    754  18288  99999  99999     90     318
      5    708  18654   -685   -875  99999  99999
      4    700  18730   -685   -875     80     43

    The mandatory data goes from 150mb to 100mb, but we also really like to have 125mb. If we use this dataset
    we can extrapolate and use 127.1mb as our estimate for 125mb

    :param ds:
    :param month: Rawinsonde are not launched at precise times. In order to make your data start on the month, pass in value
    :param levels: the mb levels we want data from, eg., 150, 125, 100, etc.
    :return:
    '''
    # handles case where netcdf has seconds since 1970-01-01, as the units but stores time in poorly defined
    # array floating seconds since 1970-01-01
    #
    xr.set_options(keep_attrs=True)  # keep any coordinate information you have when building new xarray dataset
    #time_arr = ds['relTime'].values.astype(int)
    time_arr = ds['synTime'].values.astype(int)
    print(time_arr[0])
    print(time_arr[-1])

    #Rawinsonde balloons are not launched at precise times. In order to make your data start on the month you

    try:
        # https://numpy.org/doc/stable/reference/arrays.datetime.html
        time_convert = np.array(time_arr, dtype='datetime64[s]')
    except:
        print("Still unable to properly decode time in netcdf, program exits")

    print(time_convert[0])
    print(time_convert[30])
    print(time_convert[-1])
    now = datetime.now()  # current date and time
    date_time = now.strftime("%m/%d/%Y, %H:%M:%S")

    prSigW = ds['prSigW'].values
    htSigW = ds['htSigW'].values
    wdSigW = ds['wdSigW'].values
    wsSigW = ds['wsSigW'].values
    levelsSigW = ds['sigWLevel'].values
    sfc = prSigW[0][0]
    print(f"Surface pressure now {sfc}")
    # https://glossary.ametsoc.org/wiki/Mandatory_level this dataset includes surface to the 21 mandatory levels

    #https://www.atmos.albany.edu/facstaff/ktyle/atm533/core/week7/01_Xarray_Intro.html some info from here to build datasets
    da_hgt = xr.DataArray(htSigW, coords=[time_convert, levelsSigW],
                               dims=['time','level'], attrs={'units':'m', 'long_name':'Geopotential-Significant level wrt W'})
    da_pres = xr.DataArray(prSigW, coords=[time_convert, levelsSigW],
                                dims=['time','level'],attrs={'units':'hPa','long_name':'Pressure-Significant level wrt W'})
    da_wd = xr.DataArray(wdSigW, coords=[time_convert, levelsSigW],
                               dims=['time','level'], attrs={'units':'degree', 'long_name':'Wind Direction-Significant level wrt W'})
    da_ws = xr.DataArray(wsSigW, coords=[time_convert, levelsSigW],
                                dims=['time','level'],attrs={'units':'m s**-1','long_name':'Wind Speed-Significant level wrt W'})

    ds_new = xr.Dataset(data_vars={'hgt':da_hgt, 'pres': da_pres, 'wd':da_wd, 'ws':da_ws}).astype('float32')

    #remove data that is not in september
    #https://stackoverflow.com/questions/60791186/select-xarray-dataset-based-on-month

    if len(month) > 1:
        ds = ds_new
    else:
        ds = ds_new.where((ds_new['time.month'] == month), drop=True)

    ws1 = ds['ws'].values
    ds['ws'] = ds['ws'].where(lambda x: x < 999999, np.nan)
    ws3 = ds['ws'].values

    ds['wd'] = ds['wd'].where(lambda x: x <= 999999., np.nan)
    ds['hgt'] = ds['hgt'].where(lambda x: x < 999999, np.nan)
    ds['pres'] = ds['pres'].where(lambda x: x < 999999., np.nan)
    pressure = ds['pres'].to_numpy()
    i,j = pressure.shape
    cols = len(levels)
    rows = i
    indexes = [[0 for i in range(cols)] for j in range(rows)]

    for i in range(0,rows):
        row_values = pressure[i]
        row_values = np.nan_to_num(pressure[i], copy=True, nan=-99999.0, posinf=None, neginf=None)
        for j in range(len(levels)):
            level = levels[j]
            idx = (np.abs(row_values -level)).argmin()
            indexes[i][j] = idx


    return ds, indexes



def create_radiosonde_SigW(ds):
    '''
    Create xarray dataset from radiosonde data from https://ruc.noaa.gov/raobs/ in netcdf format. Note,
    See https://ruc.noaa.gov/raobs/netcdf_format.htm for explanation of Significant levels with respect
    to wind.
    :param ds:
    :return ds: A netcdf data containing all significant levels of wind speed, direction, height and datetime
    '''

    xr.set_options(keep_attrs=True)  # keep any coordinate information you have when building new xarray dataset
    time_dimension = False
    try:
        #time_arr = ds['relTime'].values.astype(int)
        time_arr = ds['synTime'].values.astype(int)
        print(time_arr[0])
        print(time_arr[-1])

        try:
            # https://numpy.org/doc/stable/reference/arrays.datetime.html
            time_convert = np.array(time_arr, dtype='datetime64[s]')
            print(time_convert[0])
            print(time_convert[-1])
            now = datetime.now()  # current date and time
            date_time = now.strftime("%m/%d/%Y, %H:%M:%S")
            time_dimension = True
        except:
            print("Still unable to properly decode time in netcdf, program exits")

    except:
        print("No time dimensions, skip")

    height = ds['htSigW'].values
    pressure = ds['prSigW'].values
    wd = ds['wdSigW'].values
    ws = ds['wsSigW'].values
    levelsSigW = ds['sigWLevel'].values
    #sfc = pressure[0][0]

    # sometimes we just want to check data without time units, this allows us to do so
    if time_dimension:
        da_mhgt = xr.DataArray(height, coords=[time_convert, levelsSigW],
                               dims=['time','level'], attrs={'units':'m', 'long_name':'Geopotential-Significant level wrt W'})
        da_mpres = xr.DataArray(pressure, coords=[time_convert, levelsSigW],
                                dims=['time','level'],attrs={'units':'hPa','long_name':'Pressure-Significant level wrt W'})
        da_mwd = xr.DataArray(wd, coords=[time_convert, levelsSigW],
                               dims=['time','level'], attrs={'units':'degree', 'long_name':'Wind Direction-Significant level wrt W'})
        da_mws = xr.DataArray(ws, coords=[time_convert, levelsSigW],
                                dims=['time','level'],attrs={'units':'m s**-1','long_name':'Wind Speed-Significant level wrt W'})
    else:
        da_mhgt = xr.DataArray(height, coords=[levelsSigW],
                               dims=['level'], attrs={'units':'m', 'long_name':'Geopotential-Significant level wrt W'})
        da_mpres = xr.DataArray(pressure, coords=[levelsSigW],
                                dims=['level'],attrs={'units':'hPa','long_name':'Pressure-Significant level wrt W'})
        da_mwd = xr.DataArray(wd, coords=[levelsSigW],
                               dims=['level'], attrs={'units':'degree', 'long_name':'Wind Direction-Significant level wrt W'})
        da_mws = xr.DataArray(ws, coords=[levelsSigW],
                                dims=['level'],attrs={'units':'m s**-1','long_name':'Wind Speed-Significant level wrt W'})

    ds_new = xr.Dataset(data_vars={'hgt':da_mhgt, 'pres': da_mpres, 'wd':da_mwd, 'ws':da_mws}).astype('float32')

    return ds_new