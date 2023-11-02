'''
This file deals with cross-sectional plots to show wind speeds and
directions versus pressure level or time of data for a single station
using either ERA5 reanalysis data or radiosonde data or both

Author: cem
Version: 1.2 21-July-2023
'''
import argparse
import datetime
import time
import matplotlib;
matplotlib.use('Qt5Agg')  # helps pycharm from hanging when doing plt.show()
import math
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator
import metpy.calc
import xarray as xr
import datetime as dt
from datetime import datetime
import numpy as np
import pandas as pd
import sys
import build_radiosonde as br
import logging

logging.basicConfig(level=logging.WARNING)


def get_average_height(ds, key=None):
    '''
    Takes all the geopotential data and converts it to geopotential
    height and then to feet
    :param ds: xarray dataset
    :return: mean height for a given xarray dataset
    '''
    g = 9.80665  # gravitation constant used to convert geopotential height to height
    meter2feet = 3.28084
    k = meter2feet / g
    if key == None:
        da = ds['z'] * k
    elif key == 'man_hgt':
        ds[key] = ds[key].where(lambda x: x < 50000, np.nan)
        test = ds[key].values
        da = ds[key] * meter2feet

    val = da.values
    see_output = np.array2string(val, formatter={'float_kind': lambda x: "%12.2f" % x})
    print(f"Height values in meters: {see_output}")
    da_mean = np.nanmean(val)
    return da_mean


def dayofmonth_vs_height(ds, title, hour, use_day, output_dir, subtitle):
    '''
        Expects your dataset to be a netcdf containing daily ahourly
    data for one month only. If multiple months program will fail. No explicit
    check on this.
    :param ds: netcdf data from ERA5 reanalysis data containing daily values for specific month
            Should be four dimensions, time, level,[ ]latitude, longitude] and the variables should be
                u, v, and geopotential
    :param title: String containing prefix of your title, program will add date to title. Also used
                    as part of filename
    :param hour: integer hour
    :param use_day: can plot either day number (e.g., 1, 2, ... 31) or datetime
    :param output_dir: folder for image to be saved
    :param subtitle: just info about the parameters for the plot
    :return:
    '''

    if isinstance(hour, int):
        logging.info(f"Processing daily data for hour {hour}")
    else:
        hour = int(hour)

    ds['u_wind'] = ds['u'].metpy.convert_units('knots')
    ds['v_wind'] = ds['v'].metpy.convert_units('knots')
    ds['wind_speed'] = metpy.calc.wind_speed(ds['u_wind'], ds['v_wind'])
    ds['wind_dir'] = metpy.calc.wind_direction(ds['u_wind'], ds['v_wind'])

    logging.info(ds)

    try:
        ds_hr = ds.sel(time=dt.time(hour))
    except:
        print("No time attribute for this dataset, this plot requires a time series. Program ends.")
        sys.exit(-1)

    logging.info(ds_hr)
    logging.info(ds['wind_speed'])

    # only need following when we have more than one location
    #ds_pt = ds_hr.sel(longitude=lon, latitude=lat, method="nearest")
    # Now get average height corresponding to pressure levels
    levels = ds_hr['level'].to_numpy()
    y2labels = []

    # iterate through the different pressure levels and calculate their geopotential height
    for level in reversed(levels):
        dslevel = ds_hr.sel(level=level)
        height_ft = get_average_height(dslevel)
        hgt = int(math.ceil(height_ft / 1000.0))
        ft_str = str(hgt)
        y2labels.append(ft_str)


    # Get time information to help label your graph and for x-axis
    xtime = ds_hr['time'].to_numpy()
    if use_day:
        xday = pd.to_datetime(xtime).day
    else:
        xday = xtime
    logging.info(type(xtime[0]))
    year = pd.to_datetime(xtime[0]).year
    m = pd.to_datetime(xtime).month
    #month = util.num_to_month(m[0])
    month = str(m[0])
    if len(month) < 2:
        month = "0"+month
    hour_title = str(hour)
    if len(hour_title) == 1:
        hour_title = "0" + hour_title

    # create a 2-d grid with winds Y-Axis, time X-Axis
    Y, X = np.meshgrid(ds_hr['level'], xday)
    logging.info(f"X shape {X.shape}")
    logging.info(f"Y shape {Y.shape}")
    logging.info(f"X: {X}")
    logging.info(f"Y: {Y}")

    # Plot wind speed using contourf

    Z = ds_hr['wind_speed'].to_numpy()
    logging.info(Z)
    Z = Z.squeeze()    # we square arrays to get ride of elements of size 1. Here we go from 30 x 8 x 1 x 1 to 30 x 8
    logging.info(f"Z after squeeze {Z}")

    zmin = np.min(Z)
    zmax = np.max(Z)
    logging.info(f"zmin={zmin}, zmax={zmax}")

    # Now, we can make the plot.

    plt.figure(facecolor='lightgrey')
    fig = plt.figure(1, figsize=(12., 8.))
    ax = plt.axes()

    plot_contour = True
    if plot_contour:
        handle_contour = ax.contourf(X, Y, Z, levels=np.arange(0, zmax, 5.0), cmap='PuBuGn')
        #                        levels=np.arange(0, 50, 5.0), cmap='YlGnBu') plt.cm.BuPu PuBuGn
        handle_colorbar = fig.colorbar(handle_contour, pad=.07, aspect=40)
        handle_colorbar.set_label('Wind Speed (Kts)')

    U = ds_hr['u_wind'].squeeze()
    V = ds_hr['v_wind'].squeeze()

    # create a Bounding Box and then don't clip barbs
    bb = plt.barbs(X, Y, U, V, pivot="middle", color="brown", length=6, linewidth=.6)
    bb.set_clip_box(None)

    # Adjust the y-axis to be logarithmic
    ax.set_yscale('symlog')
    rlevels = levels[::-1]  # levels goes from 10mb to 150mb, we reverse this array to go from 150 to 10mb
    ax.set_yticklabels(rlevels)
    ax.set_ylim(ds_hr['level'].max(), ds_hr['level'].min())
    ax.set_yticks(rlevels)
    ax.grid(axis="y", color="white", alpha=0.9, linewidth=1.0, linestyle=":")
    ax.grid(axis="x", color="white", alpha=0.9, linewidth=1.0, linestyle=":")

    # adjust x -axis
    ax.tick_params(axis='x', pad=6)
    ax.tick_params(axis='y', pad=6)

    # width and length in pixels
    ax.tick_params(axis='both', which='major', labelsize=12, width=1.5, length=5)
    ax.xaxis.set_minor_locator(AutoMinorLocator())

    # Adding Twin Axes to plot using dataset_2
    ax2 = ax.twinx()
    ax2.set_yscale('symlog')

    t = ax2.text(1.08, 0.5, 'Mean Height (kft)', rotation=90,
                 verticalalignment='center', horizontalalignment='right',
                 transform=ax.transAxes)
    ax2.set_yticklabels(y2labels)
    ax2.set_yticks(rlevels)
    ax2.set_ylim(ds_hr['level'].max(), ds_hr['level'].min())
    ax2.tick_params(axis='both', which='major', labelsize=12, width=1.5, length=6)

    pre_title = title
    pos_title = ", " + month + '/' + str(year) + ' at ' + hour_title + 'Z'
    ax.set_title(pre_title + pos_title, y=1.02)
    ax.set_ylabel('Pressure (hPa)')
    xlabel_str = 'Day of the Month'
    ax.set_xlabel(xlabel_str)

    #
    logging.info("Finished plotting Pressure level versus day for single hour")

    # Now let's add your additional information
    if use_day:
        subtitle = subtitle+"-use_day=T"
    else:
        subtitle = subtitle+"-use_day=F"
    # add a little side note to the plot, so we know what input parameters created the plot
    ax.annotate("(" + subtitle + ")",
                xy=(0.5, 0.01), xytext=(0, 10),
                xycoords=('axes fraction', 'figure fraction'),
                textcoords='offset points',
                size=8, ha='center', va='bottom')
    file_explanation = title.replace(' ','_')
    filestr = str(year) + "_" + month + "_" + hour_title + "-" + file_explanation + ".png"
    filename = output_dir + filestr

    fig.show()
    input(' After viewing plot do you want to continue?')
    plt.savefig(filename)
    plt.close()
    return "success"




def create_radiosonde_dataset(ds, month, hour):
    '''
    Create xarray dataset from radiosonde data from https://ruc.noaa.gov/raobs/ in netcdf format
    :param ds:
    :param month: Rawinsonde are not launched at precise times. In order to make your data start on the month, pass in value
    :return:
    '''
    # handles case where netcdf has seconds since 1970-01-01, as the units but stores time in poorly defined
    # array floating seconds since 1970-01-01
    #
    xr.set_options(keep_attrs=True)  # keep any coordinate information you have when building new xarray dataset
    time_arr = ds['synTime'].values.astype(int)
    print(time_arr[0])
    print(time_arr[-1])

    # Rawinsonde balloons are not launched at precise times. In order to make your data start on the month you

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
    # https://www.atmos.albany.edu/facstaff/ktyle/atm533/core/week7/01_Xarray_Intro.html some info from here to build datasets
    da_mhgt = xr.DataArray(htMan, coords=[time_convert, mandatory_levels],
                           dims=['time', 'level'], attrs={'units': 'm', 'long_name': 'Geopotential-Mandatory level'})
    da_mpres = xr.DataArray(prMan, coords=[time_convert, mandatory_levels],
                            dims=['time', 'level'], attrs={'units': 'hPa', 'long_name': 'Pressure=Mandatory level'})
    da_mwd = xr.DataArray(wdMan, coords=[time_convert, mandatory_levels],
                          dims=['time', 'level'],
                          attrs={'units': 'degree', 'long_name': 'Wind Direction-Mandatory level'})
    da_mws = xr.DataArray(wsMan, coords=[time_convert, mandatory_levels],
                          dims=['time', 'level'], attrs={'units': 'm s**-1', 'long_name': 'Wind Speed-Mandatory level'})

    ds_new = xr.Dataset(
        data_vars={'man_hgt': da_mhgt, 'man_press': da_mpres, 'man_wd': da_mwd, 'man_ws': da_mws}).astype('float32')

    # remove data that is not in september
    # https://stackoverflow.com/questions/60791186/select-xarray-dataset-based-on-month
    if len(month) > 1:
        # ds_month = ds_new.where((ds_new['time.month'] == month), drop=True)
        # return(ds_month)
        # https://stackoverflow.com/questions/48214357/select-xarray-pandas-index-based-on-a-list-of-months
        # ds_month = ds.sel(time=np.in1d(ds['time.month'], month))
        ds_month = ds_new
    else:
        ds_month = ds_new.where((ds_new['time.month'] == month), drop=True)
        if hour != None:
            ds_hour = ds_month.where((ds_month['time.hour'] == hour), drop=True)
            return ds_hour
        else:
            return ds_month




##############################

def process_file(ncfile, option, output_dir, title, ncfile2):
    '''

    :param ncfile:
    :param option:
    :param output_dir:
    :param title:
    :param ncfile2:
    :return:
    '''

    subtitle = 'cross_section.py with option:' + option + ' ' + ncfile

    filetype1 = "nc"
    try:
        ds = xr.open_dataset(ncfile, engine="netcdf4")
    except:
        try:
            print("Radiosonde data has different time encryption, so modify open statement.")
            ds = xr.open_dataset(ncfile, engine="netcdf4", decode_times=False)
            filetype1 = "cdf"
        except:
            print("unable to decode time correctly, program aborts")
            sys.exit(-1)

    filetype2 = "nc"
    if ncfile2 != None:
        try:
            ds2 = xr.open_dataset(ncfile2, engine="netcdf4")
        except:
            try:
                print("Radiosonde data has different time encryption, so modify open statement.")
                ds2 = xr.open_dataset(ncfile2, engine="netcdf4", decode_times=False)
                filetype2 = "cdf"
            except:
                print("unable to decode time correctly, program aborts")
                sys.exit(-1)

    if option == 'n':
        ds_subset = ds.sel(level=slice(10, 150))
        use_day = False
        hour = 0
        results = dayofmonth_vs_height(ds_subset, title, hour, use_day, output_dir, subtitle)
        hour = 6
        results = dayofmonth_vs_height(ds_subset, title, hour, use_day, output_dir, subtitle)
        hour = 12
        results = dayofmonth_vs_height(ds_subset, title, hour, use_day, output_dir, subtitle)
        hour = 18
        results = dayofmonth_vs_height(ds_subset, title, hour, use_day, output_dir, subtitle)
    else:
        print("Only option n | netcdf is available in this version, program aborts")
        sys.exit(-2)


def main(argv):
    '''
    Either use defaults or command line arguments to run program
    @param argv:
    @return:
    '''

    start_time = time.time()
    print("--- %s seconds ---" % (start_time))

    if len(argv) == 0:
        print("No arguments were passed will use default values!")

    parser = argparse.ArgumentParser()
    switch_help = ("switch defines what we want to plot: "
                  "switch=mandatory | m plot the mandatory levels of radiosonde data which almost matchs the model levels"
                  "switch=significant | s plot the significant wind changes of radiosonde data"
                  "switch=netcdf | n plot model levels from NetCDF data"
                  "switch=comparison | c compare radiosonde data to model data, first file should be radiosonde")
    infile_help = ("input file, if ERA5 filetype must be '.nc'"
                   "if netcdf radiosonde data, must be '.cdf'")
    infile_help2 = ("input file 2, if ERA5 filetype must be '.nc'"
                   "if netcdf radiosonde data, must be '.cdf'")
    parser.add_argument("-i", "--infile", help=infile_help)
    parser.add_argument("-i2", "--infile2", help=infile_help2)
    parser.add_argument("-r", "--radiosonde", help="netcdf radiosonde input file, must end in '.cdf'")
    parser.add_argument("-s","--switch", help=switch_help)
    parser.add_argument("-o","--output_dir", help="where to store the image file produced")
    parser.add_argument("-t", "--title", help="prefix title of plot, will add date, also used for part of the saved image name")
    args = parser.parse_args()
    logging.info(args)
    infile = args.infile
    infile2 = args.infile2
    switch = args.switch
    output_dir = args.output_dir
    title = args.title

    # this program doesn't have defaults, so need to use your input arguments

    if infile == None:
        if switch == 'netcdf' or switch == 'n':
            logging.error("Need an input file, program aborts")
            print("Need an input file, program aborts, like -i c:\data\hilo\hilo.nc")
            sys.exit(-1)
    if infile2 == None:
        print("Only one file set, cannot run the comparison method...")

    if output_dir == None:
        logging.error("Program needs an output_dir like '-o c:\data\hilo\'")
        print("Program needs an output_dir like '-o c:\data\hilo\'")
        sys.exit(-1)
        # add a test to see if folder exists

    if title == None:
        title = "Wind Cross Section"

    if switch == None:
        logging.error("Program needs a switch option like '-s netcdf', so aborting")
        print("Program needs a switch option like '-s netcdf', so aborting")
        sys.exit(-1)
    else:
        switch = switch[0].lower()

    process_file(infile, switch, output_dir, title, infile2)
    print("--- runtime is %s seconds ---" % (time.time() - start_time))


if __name__ == "__main__":
    main(sys.argv[1:])
