'''
    Program used to study the relationship between ERA5 reanalysis winds
    and radiosonde winds.

    This program requires that you download at least one day of NetCDF data from NOAA's
    radiosonde archive site at https://ruc.noaa.gov/raobs/ (that is from hour 00z to 23Z)
    and one day of ERA5 NetCDF reanalysis data from:
    https://cds.climate.copernicus.eu/cdsapp#!/dataset/reanalysis-era5-pressure-levels?tab=form

    When downloading ERA5 reanalysis data you want to use all available pressure levels
    from 150mb to 10mb and download the following variables: U wind component, V wind component,
    and geopotential.

    Because the above two dependencies can be confusing, refer to the accompanying Readme.md file.


    Program can be run in either command mode or if no arguments are present
    using pre-defined  model and radiosonde datasets. These datasets must be present
    in the default directory if default is used (i.e., C:\Data\wind_rose\) and
    include hilo_02sep2024.nc and hilo_02sep2023.cdf. Filetype cdf is radiosonde data,
    and the other model data.

    The program is currently written to process one year of data, but it will pickup in your
    files the starting dates. Both the model and radiosonde data should be available
    because this program will make 3 sets of wind roses:
    1. model data at 00z and 12z for all available heights
    2. radiosonde data at 00z and 12z for the mandatory pressure levels which corresponde
    to all but one model pressure level
    3. All available pressure levels at "significant wind heights" (i.e., when winds change by
    either 10 degrees or 10 knots

    The end results of these wind rose comparisons is that we can a subjective view of how
    model data compares to radiosonde data in terms of prevailing directions and what model
    data misses by excluding pressure levels between 150mb and 10mb.

    One final note, model data is available hourly, while radiosonde data is available
    usually at 00z and 12z; however, there are times when say a hurricane is passing
    near a radiosonde site there will be additional launches and different sounding times.

    This comparison routine requires both 00z and 12z be available and will ignore other times.

    We require in this version that the NetCDF datasets here require one month of data
    at 00Z and 12Z. In the case of radiosonde data, e.g., during extreme weather events such
    as a hurricane, radiosonde observations may be made at different times or be missing.
    Only 00Z and 12Z datasets are used, and we require both 00Z and 12Z data or day is skipped.

    Default sample data is from Hilo, Hawaii located at 19.72N, 155.07W for 02 September 2023.
    Normally the datasets would consist of thirty days of data for monthly analysis but for demonstration
    purposes data are not analysed just displayed as an image.


Author: cem
Version: 1.0 28-January 2022
        2.0 12-August-2023 (added ability to write files to csv) with summary stats
        2.1 09-September-2023 (removed csv file writes for simplicity)
        2.1.1 09-November-2023 (added some documentation highlight binning)

see also:

'''

import sys
import pandas as pd
import xarray as xr
import numpy as np
from numpy import pi
import argparse
import math
from math import radians

import matplotlib;
matplotlib.use('Qt5Agg')  # stops figure from hanging when plotted from pycharm
import matplotlib.pyplot as plt
from matplotlib import colors as mcolors
from matplotlib.gridspec import GridSpec
#import pathlib
from datetime import datetime
import re
import time
import metpy.calc

#local libraries
import build_radiosonde as br

import logging
logging.basicConfig(level=logging.WARNING)


def firstNan(lists):
    '''
    Find index for first NaN
    @param lists:
    @return:
    '''
    index = 0
    for item in lists:
        if math.isnan(item) == True:
            return index
        else:
            index = index + 1



def cycle_model(ds, title, output_dir):
    '''
    Using all reanalysis data between 150mb to 10mb, create wind roses at 00z and 12z.
    Program is designed to iterate through all available data for one year. If more than one year
    of data, will ignore data after the first year.

    For each year and each month, plot daily windroses for the two time periods.

    @param ds: Contains u, v (m/s) and geopotential (not heights just geopotential) xarray data-array
    @param title: Overall directory
    @param output_dir: Output folder
    @return:
    '''

    # need a global variables to store wind speeds between calls
    global ns_winds
    ns_winds = np.zeros((16, 4), dtype=np.float32)

    now = datetime.now()
    date_time = now.strftime("%Y-%m-%d %H:%M:%S")

    # We only care about from 150 to 10mb, note this may be the only data available anyway
    ds = ds.sel(level=slice(10, 150))

    mb = ds['level'].values

    # Get start and end days of each month for this year
    year = ds.time.dt.year # https://docs.xarray.dev/en/latest/user-guide/time-series.html (see section on datetime components)
    yr = year[0].values
    months, start_days, end_days, num_days = get_dates(yr)

    iter = 0  # total iterations
    title_prefix = title
    for month in range(1, 13, 1):
        start = start_days[month - 1]
        end = end_days[month - 1]
        num = num_days[month - 1] + 1

        ds_month = ds.sel(time=slice(start, end))
        shape_of_time_coord = tuple(ds_month.coords[d] for d in ['time', 'level'])
        shape_of_time_dims = tuple(ds_month.dims[d] for d in ['time', 'level'])
        observations = shape_of_time_dims[0]
        if observations == 0:
            print(f"For month {month} no data was found, skipping to next month.")
            continue

        for day in range(1, num, 1):
            # title = "Ishigakijima Island, Japan for "+months[month+1]+"  2021, Day "+str(day)
            title = title_prefix + " for " + months[month - 1] + " " + str(yr) + ", Day " + str(day)
            title = title.replace("_", ", ")
            today = str(yr)+'-'+str(month)+'-'+str(day)
            #https://stackoverflow.com/questions/51976126/subset-an-xarray-dataset-or-data-array-in-the-time-dimension
            ds_dy = ds_month.sel(time=slice(today,today))
            #ds_dy = ds_month.where((ds_month['time.day'] == day), drop=True)    # get ride of data except for current day
            size = ds_dy.time.sizes['time']
            bad_values = [0, 0]
            if size > 2:  # sometimes more than 00z and 12z data, then remove
                print("Got more than 00z and 12z, removing other times...")
                times0 = ds_dy['time'].values
                wd0 = ds_dy['wd'].values
                ds_dy = reduce_time(ds_dy, size)
                wd1 = ds_dy['wd'].values
                times1 = ds_dy['time'].values
                # print(f"before: {times0}, {wd0}")
                # print(f"after: {times1}, {wd1}")
                size = ds_dy.time.sizes['time']
            elif size == 2:  # normally every 12 hours or twice per day
                # for k in range(size):
                #    print(ds_dy['time'][k].values)
                good_data = True
            elif size == 1:  # need 00z and 12z data or skip
                print("Only got one hour of data, skipping ...")
                continue
            else:
                continue

            plt_0012Z_windrose(ds_dy, title, output_dir)

            iter = iter + 1


def cycle_raob(ds, switch, title, output_dir):
    '''
    Function processes specifically for radiosonde data

    @param ds: xarray dataset containing radiosonde data. Here we assume that the dataset only contains
               one year of data at most. But we don't check. Would be easy to just let program go
               on and plot any number of days of data.
    @param year: string of year (yyyy)
    @param month: string of month (mm)
    @param title: string title of plot
    @param output_dir: string, output directory
    @return:
    '''

    # need a global variables to store wind speeds between calls
    global ns_winds
    ns_winds = np.zeros((16, 4), dtype=np.float32)
    year = ds.time.dt.year # https://docs.xarray.dev/en/latest/user-guide/time-series.html (see section on datetime components)
    yr = year[0].values
    months, start_days, end_days, num_days = get_dates(yr)

    iter = 1  # total iterations
    title_prefix = title
    for month in range(1, 13, 1):
        start = start_days[month - 1]
        end = end_days[month - 1]
        num = num_days[month - 1] + 1
        ds_month = ds.sel(time=slice(start, end))
        #shape_of_time_coord = tuple(ds_month.coords[d] for d in ['time', 'level'])
        shape_of_time_dims = tuple(ds_month.dims[d] for d in ['time', 'level'])
        observations = shape_of_time_dims[0]
        if observations == 0:
            print(f"For month {month} no data was found, skipping to next month.")
            continue
        good_data = 0

        for day in range(1, num, 1):
            title = title_prefix + ' on '+months[month - 1] + " " + str(yr) + ", Day " + str(day)
            title = title.replace("_", ", ")
            ds_dy = ds_month.where((ds_month['time.day'] == day), drop=True)
            size = ds_dy.time.sizes['time']
            # bad_values = [0, 0] (when we use this we want to make sure we have over 1 observation in time period
            if size > 2:  # sometimes more than two radiosondes per 24 hours
                times0 = ds_dy['time'].values
                ds_dy = reduce_time(ds_dy, size)
            elif size == 2:
                date_time = np.empty(size, dtype=object)
                # normally every 12 hours or twice per day
                for k in range(size):
                    date_time[k] = ds_dy['time'][k].values
                # Here we are testing to get either 00Z or 12Z and not some odd time like 18Z
                hr_test = ["00", "12"]
                hr_index = []
                for k in range(size):
                    dtString = str(date_time[k]).split('T')[1]
                    hr = dtString[0:2]
                    hr_index.append(hr.find(hr_test[k]))

                if hr_index[0] == -1 or hr_index[1] == -1:
                    bad_hour = True
                for k in range(size):
                    print(ds_dy['time'][k].values)
            elif size == 1:  # once in a while just one launch per day
                print("Only one of the two 12 hourly soundings were available, skip")
                continue
            elif size == 0:
                print("No valid data found")
                continue

            #if bad_values[0] > 0 and bad_values[1] > 0:  # if both 12 hour periods of data are bad
            #    bad_data = bad_data + 1
            #else:
            good_data = good_data + 1
            plt_raob(ds_dy, switch, title, output_dir)
            iter = iter + 1


def dual_wind_rose(fig, df1, df2, kft):
    '''

    :param df:
    :param subtitle:
    :return:
    '''

    # Creates DataFrame.
    # df = pd.DataFrame(data)
    # df1 = df.iloc[0]
    # df2 = df.iloc[1]
    # print(df1)
    # print(df2)

    # get max wind for wind rose radii
    ws_max1 = max(df1.ws)
    ws_max = max(df2.ws)
    if ws_max1 > ws_max:
        ws_max = ws_max1

    bar_color = ['yellow', 'orange', 'red', 'lime', 'green', 'lightblue', 'blue', 'purple']
    N = len(df2.index)
    nbins = 16
    deg2rad = pi / 180.0
    # Fixing random state for reproducibility
    bins = np.arange(- (180 / nbins), 360 + (180 / nbins), 360 / nbins)
    bins_in_radians = bins * deg2rad
    increment_per_bin = (360 / nbins) * 0.5
    bin_center_pt = bins + increment_per_bin
    bin_center_pt_radians = bin_center_pt * deg2rad
    theta_angles = bin_center_pt_radians[:-1]
    quadrant0_threshold = 360 - (180 / nbins)  # anything greater needs to shifted to fit into bins
    angle_size = (2 * np.pi / nbins) - 0.02
    # width = angle_size

    # bin = np.zeros(nbins)
    # http://www.chiark.greenend.org.uk/~peterb/python/polar/
    data = pd.DataFrame({'compass': ['N', 'NNE', 'NE', 'ENE', 'E', 'ESE', 'SE', 'SSE', 'S', 'SSW', 'SW', 'WSW', 'W',
                                     'WNW', 'NW', 'NNW'],
                         'degrees': [0, 22.5, 45, 67.5, 90, 112.5, 135, 157.5, 180, 202.5, 225, 247.5, 270, 292.5, 315,
                                     337.5]})

    df = df1.copy()
    calm = []
    width = []
    slice = pi / 9
    for i in range(len(df['wd'])):
        wspd = df.loc[i, ('ws')]
        if quadrant0_threshold <= df.loc[i, ('wd')] <= 360:
            df.loc[i, ('wd')] = df.loc[i, ('wd')] - 360
        # print(f"{i}: {df.loc[i, ('wd')]}, {df.loc[i, ('ws')]}")
        n, _ = np.histogram(df.loc[i, ('wd')], bins=bins)
        bin_index = n.argmax()
        df.loc[i, ('bin')] = bin_index
        df.loc[i, ('theta')] = theta_angles[bin_index]

        if wspd < 5.0:
            calm.append(1)
            width.append(2 * pi)
            df.loc[i, ('ws')] = 5
        else:
            calm.append(0)
            width.append(slice)

    # df = df.assign(wind_dir = lambda x: (x['wd']) * deg2rad)
    radii = df['ws']
    theta = df['theta']

    # fig = plt.figure(figsize=(8, 3))
    # fig.patch.set_facecolor('lightgoldenrodyellow')
    gs = GridSpec(nrows=1, ncols=2, width_ratios=[1, 1])

    ax1 = fig.add_subplot(gs[0, 0], projection='polar')
    ax1.patch.set_facecolor('white')
    ax1.set_title('00Z run', fontsize=12)
    xticks = 16
    ax1.set_xticks([radians(x) for x in np.arange(0, 360, 360 / xticks)])
    ax1.set_xticklabels(data.degrees)
    ax1.set_theta_zero_location('N')
    ax1.set_theta_direction(-1)
    ax1.set_rlabel_position(325)
    round = int((ws_max + 4) / 5) * 5

    bars1 = ax1.bar(x=theta, height=radii, width=width, edgecolor='k', linewidth=1.0)
    # fig.show()
    ax1.set_rlim(0, round)
    # fig.show()
    df = df2.copy()
    calm = []
    width = []
    for i in range(len(df['wd'])):
        wspd = df.loc[i, ('ws')]
        if quadrant0_threshold <= df.loc[i, ('wd')] <= 360:
            df.loc[i, ('wd')] = df.loc[i, ('wd')] - 360
        # print(f"{i}: {df.loc[i, ('wd')]}, {df.loc[i, ('ws')]}")
        n, _ = np.histogram(df.loc[i, ('wd')], bins=bins)
        bin_index = n.argmax()
        df.loc[i, ('bin')] = bin_index
        df.loc[i, ('theta')] = theta_angles[bin_index]
        if wspd < 5.0:
            calm.append(1)
            width.append(2 * pi)
            df.loc[i, ('ws')] = 5
        else:
            calm.append(0)
            width.append(slice)

    # df = df.assign(wind_dir = lambda x: (x['wd']) * deg2rad)
    radii = df['ws']
    theta = df['theta']

    ax2 = fig.add_subplot(gs[0, 1], projection='polar')
    ax2.patch.set_facecolor('white')
    ax2.set_title('12Z run', fontsize=12)
    xticks = 16
    ax2.set_xticks([radians(x) for x in np.arange(0, 360, 360 / xticks)])
    ax2.set_xticklabels(data.compass)
    ax2.set_theta_zero_location('N')
    ax2.set_theta_direction(-1)
    ax2.set_rlabel_position(325)
    round = int((ws_max + 4) / 5) * 5
    bars2 = ax2.bar(x=theta, height=radii, width=width, edgecolor='k', linewidth=1.0)
    # fig.show()
    ax2.set_rlim(0, round)
    # fig.show()

    # only plot if you have calm data

    labels = []
    for i in range(len(df['level'])):
        level = df.loc[i, ('level')]
        legend_label = str(level) + 'mb (' + kft[i] + 'kft)'
        # x = ax.annotate(legend_label, xy=(0.5, 0.5 - i*.1),
        #                ha='center', va='bottom', xycoords='figure fraction',
        #                bbox=dict(boxstyle='round', fc=bar_color[i], alpha=0.5))

        labels.append(legend_label)
    # labels = ['10','20','30','50','70','100','125','150']
    i = 0
    for r, bar, ll in zip(theta, bars1, labels):
        # bar.set_facecolor(plt.cm.jet(r / np.pi / 2))
        bar.set_facecolor(bar_color[i])
        bar.set_alpha(0.45)
        bar.set_label(ll)
        i += 1

    # fig.show()
    i = 0
    for r, bar, ll in zip(theta, bars2, labels):
        # bar.set_facecolor(plt.cm.jet(r / np.pi / 2))
        bar.set_facecolor(bar_color[i])
        bar.set_alpha(0.45)
        bar.set_label(ll)
        i += 1

    ax2.legend(title="Winds by Level", fancybox=True, bbox_to_anchor=[0.98, 0.9], loc='center', framealpha=0.5,
               prop=dict(size=8))
    # fig.suptitle(subtitle, fontsize=12)
    # fig.show()
    return ax1, ax2

def dual_wind_rose_sigW(fig, df1, df2):
    '''
    Plots a wind rose for 00z and 12z for a given day.
    This plot is designed to plot as many as n Levels.
    It was designed to plot Significant radiosonde levels.
    But also used to plot mandatory levels; however, when plotting
    mandatory levels we place the legend in an easier to read location

    :param df1: 00z data
    :param df2: 12Z data
    :param fig: matplotlib fig object
    :return:
    '''

    # get max wind for wind rose radii
    ws_max1 = max(df1.ws)
    ws_max = max(df2.ws)
    if ws_max1 > ws_max:
        ws_max = ws_max1

    N = len(df2.index)
    nbins = 16
    deg2rad = pi / 180.0
    # Fixing random state for reproducibility
    bins = np.arange(- (180 / nbins), 360 + (180 / nbins), 360 / nbins)
    bins_in_radians = bins * deg2rad
    increment_per_bin = (360 / nbins) * 0.5
    bin_center_pt = bins + increment_per_bin
    bin_center_pt_radians = bin_center_pt * deg2rad
    theta_angles = bin_center_pt_radians[:-1]
    quadrant0_threshold = 360 - (180 / nbins)  # anything greater needs to shifted to fit into bins
    angle_size = (2 * np.pi / nbins) - 0.02
    # width = angle_size

    # http://www.chiark.greenend.org.uk/~peterb/python/polar/
    data = pd.DataFrame({'compass': ['N', 'NNE', 'NE', 'ENE', 'E', 'ESE', 'SE', 'SSE', 'S', 'SSW', 'SW', 'WSW', 'W',
                                     'WNW', 'NW', 'NNW'],
                         'degrees': [0, 22.5, 45, 67.5, 90, 112.5, 135, 157.5, 180, 202.5, 225, 247.5, 270, 292.5, 315,
                                     337.5]})

    df = df1.copy()
    calm = []
    width = []
    slice = pi / 9
    for i in range(len(df['wd'])):
        wspd = df.loc[i, ('ws')]
        if quadrant0_threshold <= df.loc[i, ('wd')] <= 360:
            df.loc[i, ('wd')] = df.loc[i, ('wd')] - 360
        n, _ = np.histogram(df.loc[i, ('wd')], bins=bins)   # use numpy's histogram to bin wind observations
        bin_index = n.argmax()
        df.loc[i, ('bin')] = bin_index
        df.loc[i, ('theta')] = theta_angles[bin_index]

        if wspd < 5.0:
            calm.append(1)
            width.append(2 * pi)
            df.loc[i, ('ws')] = 5
        else:
            calm.append(0)
            width.append(slice)

    # df, calm, width = get_theta(df1)

    radii1 = df['ws']
    theta1 = df['theta']

    gs = GridSpec(nrows=1, ncols=2, width_ratios=[1, 1])

    ax1 = fig.add_subplot(gs[0, 0], projection='polar')
    ax1.patch.set_facecolor('white')
    ax1.set_title('00Z run', fontsize=12)
    xticks = 16
    ax1.set_xticks([radians(x) for x in np.arange(0, 360, 360 / xticks)])
    ax1.set_xticklabels(data.degrees)
    ax1.set_theta_zero_location('N')
    ax1.set_theta_direction(-1)
    ax1.set_rlabel_position(325)
    round = int((ws_max + 4) / 5) * 5

    bars1 = ax1.bar(x=theta1, height=radii1, width=width, edgecolor='k', linewidth=1.0)
    # fig.show()
    ax1.set_rlim(0, round)

    df = df2.copy()
    calm = []
    width = []
    for i in range(len(df['wd'])):
        wspd = df.loc[i, ('ws')]
        if quadrant0_threshold <= df.loc[i, ('wd')] <= 360:
            df.loc[i, ('wd')] = df.loc[i, ('wd')] - 360
        # print(f"{i}: {df.loc[i, ('wd')]}, {df.loc[i, ('ws')]}")
        n, _ = np.histogram(df.loc[i, ('wd')], bins=bins)
        bin_index = n.argmax()
        df.loc[i, ('bin')] = bin_index
        df.loc[i, ('theta')] = theta_angles[bin_index]
        if wspd < 5.0:
            calm.append(1)
            width.append(2 * pi)
            df.loc[i, ('ws')] = 5
        else:
            calm.append(0)
            width.append(slice)

    # df = df.assign(wind_dir = lambda x: (x['wd']) * deg2rad)
    radii2 = df['ws']
    theta2 = df['theta']

    ax2 = fig.add_subplot(gs[0, 1], projection='polar')
    ax2.patch.set_facecolor('white')
    ax2.set_title('12Z run', fontsize=12)
    xticks = 16
    ax2.set_xticks([radians(x) for x in np.arange(0, 360, 360 / xticks)])
    ax2.set_xticklabels(data.compass)
    ax2.set_theta_zero_location('N')
    ax2.set_theta_direction(-1)
    ax2.set_rlabel_position(325)
    round = int((ws_max + 4) / 5) * 5
    bars2 = ax2.bar(x=theta2, height=radii2, width=width, edgecolor='k', linewidth=1.0)
    fig.show()
    ax2.set_rlim(0, round)

    color_names = get_colors()
    labels = []

    # If we only have a few layers, place legend off to the top right
    # otherwise we need to center the legend cause we may have over 60 layers
    legend_location = 'center'

    N = len(df1.index)
    if N < 9:
        legend_pos = [1.05, 1.05]
    else:
        legend_pos = [0.98, 0.5]

    for i in range(len(df1['height_ft'])):
        level = df1.loc[i, ('press')]
        hgt = df1.loc[i, ('height_ft')]
        kft = str(int(hgt / 1000))
        legend_label = str(level) + 'mb (' + kft + 'kft)'
        labels.append(legend_label)

    i = 0
    for r, bar, ll in zip(theta1, bars1, labels):
        if i == 37:
            i = 0
        bar.set_facecolor(color_names[i])
        bar.set_alpha(0.45)
        bar.set_label(ll)
        i += 1

    ax1.legend(title="Winds by Level", fancybox=True, bbox_to_anchor=legend_pos, loc=legend_location, framealpha=0.5,
               prop=dict(size=8))
    fig.show()

    N = len(df2.index)
    if N < 9:
        legend_pos = [1.05, 1.05]
    else:
        legend_pos = [0.98, 0.5]

    labels = []
    for i in range(len(df2['height_ft'])):
        level = df2.loc[i, ('press')]
        hgt = df2.loc[i, ('height_ft')]
        kft = str(int(hgt / 1000))
        legend_label = str(level) + 'mb (' + kft + 'kft)'
        labels.append(legend_label)

    i = 0
    for r, bar, ll in zip(theta2, bars2, labels):
        if i == 37:
            i = 0  # start over
        bar.set_facecolor(color_names[i])
        bar.set_alpha(0.45)
        bar.set_label(ll)
        i += 1

    ax2.legend(title="Winds by Level", fancybox=True, bbox_to_anchor=legend_pos, loc=legend_location, framealpha=0.5,
               prop=dict(size=8))

    return ax1, ax2



def get_average_height(ds):
    '''
    Takes all the geopotential data and converts it to geopotential
    height and then to feet
    :param ds: xarray dataset
    :return: mean height for a given xarray dataset
    '''
    g = 9.80665  # gravitation constant used to convert geopotential height to height
    meter2feet = 3.28084
    k = meter2feet / g
    da = ds['z'] * k
    da_mean = da.mean()
    return da_mean.data


def get_average_snd_height(ds):
    '''
    Convert geopotential based heights to a mean
    height in feet

    @param ds: contains the geopotential height in feet
    @return: mean height in feet
    '''
    meter2feet = 3.28084
    da = ds['hgt'] * meter2feet
    da_mean = da.mean()
    logging.info(f"mean height: {da_mean.values}")
    return da_mean.data


def get_colors():
    '''
    This function gives the names of common colors. It allows us to
    get a list of names and then build colormaps without knowning how
    many different colors we will need. And if we change
    the boolean 'show_colors' in this function it will actually plot
    the colors
    @return: A dictionary of color names
    '''
    colors = dict(mcolors.BASE_COLORS, **mcolors.CSS4_COLORS)

    # Sort colors by hue, saturation, value and name.
    by_hsv = sorted((tuple(mcolors.rgb_to_hsv(mcolors.to_rgba(color)[:3])), name)
                    for name, color in colors.items())
    sorted_names = [name for hsv, name in by_hsv]
    color_names = []

    n = len(sorted_names)
    ncols = 4
    nrows = n // ncols + 1

    show_colors = False
    if show_colors:
        fig, ax = plt.subplots(figsize=(8, 5))

        # Get height and width
        X, Y = fig.get_dpi() * fig.get_size_inches()
        h = Y / (nrows + 1)
        w = X / ncols

    for i, name in enumerate(sorted_names):
        col = i % ncols
        if col != 0:
            continue

        row = i // ncols
        if row < 2:
            continue

        if show_colors:
            y = Y - (row * h) - h

            xi_line = w * (col + 0.05)
            xf_line = w * (col + 0.25)
            xi_text = w * (col + 0.3)

            ax.text(xi_text, y, name, fontsize=(h * 0.8),
                    horizontalalignment='left',
                    verticalalignment='center')

            ax.hlines(y + h * 0.1, xi_line, xf_line,
                      color=colors[name], linewidth=(h * 0.6))

            ax.set_xlim(0, X)
            ax.set_ylim(0, Y)
            ax.set_axis_off()

            fig.subplots_adjust(left=0, right=1,
                                top=1, bottom=0,
                                hspace=0, wspace=0)
            plt.show()
        else:
            color_names.append(name)

    return color_names


def get_dates(year):
    '''
    Function returns the number of days for a given month and year.
    This function is just a short-cut and nothing changes unless we are working with a leap year
    @param year: string for year (yyy)
    @return: array of 3 letter months names, array of start date string, end date string, number of days per month
    '''

    months = ["Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"]
    y = str(year)
    start_days = []
    end_days = []
    num_days = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
    yr = int(year)
    leap_year = False
    if (yr % 400 == 0) and (year % 100 == 0):
        leap_year = True
    elif (yr % 4 == 0) and (yr % 100 != 0):
        leap_year = True


    if leap_year:
        num_days[1] = 29

    for i in range(1, 13):
        mon = str(i)
        if len(mon) < 2:
            mon = "0" + mon
        start_days.append(y + "-" + mon + "-01")
        end_days.append(y + "-" + mon + "-" + str(num_days[i - 1]))

    return months, start_days, end_days, num_days


# Note, the following four functions return either arrays or dataframes. We could standardize
# but here we are playing with combinations for learning which is most convenient
def get_df_wind_at(ds):
    '''
    converts u and v wind components to wind speed and direction
    from model data. Then in place, converts a xarray dataset to
    a dataframe and returns with addition parameters using as
    u and v wind components in knots (u_wind, v_wind), plus
    wind_speed and direction.

    :param ds: xarray dataset containing u and v wind components
    :return: panda dataframe with wind speed and direction added


    '''
    ds['u_wind'] = ds['u'].metpy.convert_units('knots')
    ds['v_wind'] = ds['v'].metpy.convert_units('knots')
    ds['wind_speed'] = metpy.calc.wind_speed(ds['u_wind'], ds['v_wind'])
    ds['wind_dir'] = metpy.calc.wind_direction(ds['u_wind'], ds['v_wind'])
    df = ds.to_dataframe()
    df.reset_index(inplace=True)

    return df


def get_data_at(ds):
    '''
    Converts xarray data-arrays from a dataset to numpy arrays.
    Along the way, converts wind speed in meters per second to knots,
    and sets missing data to numpy NaN.

    :param ds: xarray dataset containing wind speed and direction, and pressure (mb) and heights
    :return: returns data arrays containing wind_dir and wind_spd in kts, pressure (mb), heights
    '''

    ws = ds['ws'].to_numpy() * 1.94384
    wd = ds['wd'].to_numpy()
    pres = ds['pres'].to_numpy()
    hgt = ds['hgt'].to_numpy()

    ws = ws.squeeze()
    wd = wd.squeeze()
    pres = pres.squeeze()
    hgt = hgt.squeeze()

    # standardize missing number to NaN
    missing = 99999999
    ws[ws > missing] = np.nan
    hgt[hgt > missing] = np.nan
    wd[wd > missing] = np.nan
    pres[pres > missing] = np.nan
    return ws, wd, pres, hgt




def get_wind_at(ds, radiosonde):
    '''
    Converts xarray dataset containing wind speed and direction
    to numpy arrays. Along the way, converts meters per second to knots.

    :param radiosonde: boolean, true if radiosonde data
    :param ds: xarray dataset containing wind speed and direction
    :return: returns data arrays containing wind_dir and wind_spd in kts
    '''
    if radiosonde:
        ws = ds['ws'].to_numpy() * 1.94384
        wd = ds['wd'].to_numpy()
    else:
        ds['u_wind'] = ds['u'].metpy.convert_units('knots')
        ds['v_wind'] = ds['v'].metpy.convert_units('knots')
        ds['wind_speed'] = metpy.calc.wind_speed(ds['u_wind'], ds['v_wind'])
        ds['wind_dir'] = metpy.calc.wind_direction(ds['u_wind'], ds['v_wind'])
        ws = ds['wind_speed'].to_numpy()
        wd = ds['wind_dir'].to_numpy()

    ws = ws.squeeze()
    wd = wd.squeeze()
    return ws, wd



def get_raob_df(plot_num, ws, wd, pres, hgt):
    '''
    We want to plot wind roses at 00z and 12z for all available levels. Because here we deal with winds and speeds
    at "significant" wind levels we are dealing with a variable amount of data. Significant wind levels occur when there
    is a significant change in wind in speed (10kts or more) or direction (10 degrees or more).

    We are only interested in winds between approximately 150mb to 10mb (roughly 47,000 ft to 102,000ft).
    These are close to tropics and lower-subtropics height levels for era5 data over the Pacific.
    These winds are in the stratosphere of this region of the globe.

    Here we create a panda dataframe to hold the plotting data

    @param plot_num:
    @param ws:
    @param wd:
    @param pres:
    @param hgt:
    @param title:
    @param month:
    @param year:
    @return: df
    '''

    length_ws = firstNan(ws) - 1
    length_wd = firstNan(wd) - 1
    length_pres = firstNan(pres) - 1
    length_hgt = firstNan(hgt) - 1
    length = min(length_wd, length_ws, length_pres, length_hgt)

    nbins = 16
    quadrant0_threshold = 360 - (180 / nbins)

    if length == 0:
        row0Bad = True

    n_calm = np.zeros(length + 1, dtype=int)
    df = pd.DataFrame(columns=['tau', 'height_m', 'height_ft', 'press', 'ws', 'wd'])
    frames = []

    if plot_num == 0:
        tau = '00Z'
    else:
        tau = '12Z'

    ft_str1_mb = "NA"
    ft_str0_mb = "NA"
    end_level = None
    beg_level = None
    height_ft = None
    good_data = 0
    for col in range(0, length + 1):
        z = hgt[col]
        height_ft = hgt[col] * 3.28084
        if height_ft < 45000.0:  # skip heights to we get to 45,000 ft
            continue
        elif height_ft > 102000:
            height_ft = hgt[col - 1] * 3.28084
            end_level = pres[col - 1]
            break

        level = pres[col]
        if level > 150:
            continue

        good_data += 1
        height = int(round(height_ft / 1000.0))
        hgt_str = str(height) + ' kft'

        if ft_str0_mb == "NA":
            ft_str0_mb = hgt_str
            beg_level = level

        spd = ws[col]
        if quadrant0_threshold <= wd[col] <= 360:
            dir = wd[col] - 360.0  # here we remap to different wind direction coordinates
        else:
            dir = wd[col]

        if spd < 5.0:
            dir = np.NaN
            n_calm[col] = 1

        dic = {'tau': tau, 'height_m': hgt[col], 'height_ft': height_ft, 'press': pres[col], 'ws': spd, 'wd': dir}
        frames.append(dic)

    if good_data < 6:
        print("No suitable data found today. Usually because of missing upper-levels at 10 and 20mb")
        return df, "no_data", n_calm
    else:
        h = int(round(height_ft / 1000))
        ft_str1_mb = str(h) + 'kft'
        if end_level == None:
            end_level = pres[col]
        if beg_level == None:
            beg_level = pres[col]

        df = pd.concat([df, pd.DataFrame(frames)], axis=0, ignore_index=True)

        subplot_title = 'Wind rose from ' + str(beg_level) + 'mb (' + ft_str0_mb + ') to ' + str(
            end_level) + 'mb (' + ft_str1_mb + ')'

        if end_level > 150.0:
            print("No suitable data found today")

    return df, subplot_title, n_calm


def process_file(infile, option, output_dir, title):
    '''
    Function reads in data, looking for data at two times, 00z and 12z.
    It will process the data to make wind roses at either 00z or 12z from
    either ERA5 reanalysis or radiosonde data from the NOAA archive

    @param infile: string with model full path filename
    @param option: string with value 'model', or 'significant' or 'mandatory'
    @param outputdir: string with full path for output directory
    @param title: string with title
    @return:
    '''

    xr.set_options(keep_attrs=True) # xarray method to keep original attributes when reading or modifying a dataset

    switch = option.lower()

    if switch == 's' or switch == 'm':
        print("Selected option to compare radiosonde data to ERA5 model data")
        if infile.find(".cdf") != -1:
            print(f"Processing radiosonde data: {infile}")
            ds = xr.open_dataset(infile, engine="netcdf4", decode_times=False)

            try:
                lat = ds['staLat'][0].item(0)
                lon = ds['staLon'][0].item(0)  # best way to convert to a scalar
            except:
                lat = ds['staLat'].item(0)
                lon = ds['staLon'].item(0)

            if switch == 's':
                ds = br.create_radiosonde_SigW(ds)
            else:
                ds = br.create_radiosonde_man(ds)

            mytime = ds.coords['time'].values
            dt = mytime[0]
            ts = pd.to_datetime(str(dt))
            dateit = ts.strftime('%Y-%m-%d')
            raob_year = re.search(r"(\d{4})", dateit).group(1)
            raob_month = re.search(r"-(\d{2})", dateit).group(1)
            print(f"The datasets starting year {raob_year} and starting month {raob_month}")

            cycle_raob(ds, switch, title, output_dir)

        else:
            logging.error("Raob dataset does not have the right file type, will exit")
            print("Raob dataset does not have the right file type, will exit")
            sys.exit(-1)
    else:
        print(f"Processing era data: {infile} with option {option}")
        ds2 = xr.open_dataset(infile, engine="netcdf4")

        lat = ds2['latitude'][0].values
        lon = ds2['longitude'][0].values

        mytime = ds2.coords['time'].values
        dt = mytime[0]
        ts = pd.to_datetime(str(dt))

        dateit = ts.strftime('%Y-%m-%d')
        #model_year = re.search(r"(\d{4})", dateit).group(1)
        #model_month = re.search(r"(-\d{2})", dateit).group(1)

        lat_str = str(int(lat))
        lon_str = str(int(lon))
        print(f"Latitude={lat_str}, longitude={lon_str} from data at {dateit}")

        cycle_model(ds2, title, output_dir)

    print(f"Finished runing wind_rose with option {option} using file {infile}")




def plt_0012Z_windrose(ds, title, output_dir):
    '''
    Here we take model data from just one day. Your model at this point (ds) should
    only contain a day of data. Actually, can contain more I guess but purpose is to
    show plot for single data at 00z and 12z and number of calm observations.

    :param ds:  xarray data for one day of model data
    :param title: title for your plot
    :param ouput_dir:
    :return:

    '''

    ws, wd = get_wind_at(ds, False)
    # print(f"wind direction: {wd}")
    # print(f"wind speed: {ws}")

    nbins = 16
    levels = ds['level'].values
    nlevels = len(levels)

    w_dir = np.zeros((2, nlevels), dtype=np.float32)
    w_spd = np.zeros((2, nlevels), dtype=np.float32)
    n_calm = np.zeros(nlevels, dtype=int)
    n_calm1 = np.zeros(nlevels, dtype=int)
    n_calm2 = np.zeros(nlevels, dtype=int)

    df1 = pd.DataFrame(columns=['tau', 'level', 'ws', 'wd'])
    df2 = pd.DataFrame(columns=['tau', 'level', 'ws', 'wd'])
    frames1 = []
    frames2 = []

    row = 0  # tau 00Z
    quadrant0_threshold = 360 - (180 / nbins)
    for col in range(0, nlevels):
        spd = ws[row][col]  # should probably delete this don't see why I have this
        if quadrant0_threshold <= wd[row][col] <= 360:
            dir = wd[row][col] - 360.0  # here we remap to different wind direction coordinates
        else:
            dir = wd[row][col]
        if spd < 5.0:
            dir = np.NaN
        dic = {'tau': '00Z', 'level': levels[col], 'ws': spd, 'wd': dir}
        frames1.append(dic)
    df1 = pd.concat([df1, pd.DataFrame(frames1)], axis=0, ignore_index=True)

    row = 1  # tau 12Z
    for col in range(0, nlevels):
        spd = ws[row][col]  # should probably delete this don't see why I have this
        if quadrant0_threshold <= wd[row][col] <= 360:
            dir = wd[row][col] - 360.0  # here we remap to different wind direction coordinates
        else:
            dir = wd[row][col]
        if spd < 5.0:
            dir = np.NaN
        dic = {'tau': '12Z', 'level': levels[col], 'ws': spd, 'wd': dir}
        frames2.append(dic)

    # https://stackoverflow.com/questions/72657415/fix-futurewarning-related-to-the-pandas-append-function
    df2 = pd.concat([df2, pd.DataFrame(frames2)], axis=0, ignore_index=True)

    # Here we use numpy's histogram program to bin wind observations into 16 quadrants
    bins = np.arange(- (180 / nbins), 360 + (180 / nbins), 360 / nbins)

    n1, bins = np.histogram(df1.wd, bins=bins)
    n2, bins = np.histogram(df2.wd, bins=bins)

    # now check for calm events
    for cols in range(0, nlevels):
        dir1 = df1.loc[cols, ('wd')]
        dir2 = df2.loc[cols, ('wd')]

        spd1 = df1.loc[cols, ('ws')]
        spd2 = df2.loc[cols, ('ws')]

        if spd1 <= 5.0 and spd2 <= 5.0:
            n_calm[cols] = 1
            n_calm1[cols] = 1
            n_calm2[cols] = 1
        elif spd1 <= 5.0 and spd2 > 5.0:
            n_calm1[cols] = 1
        elif spd2 <= 5.0 and spd1 > 5.0:
            n_calm2[cols] = 1

    # loop to see what the average height (ft) is at each mb level
    levels_kft = []
    for level in levels:
        dslevel = ds.sel(level=level)
        df = get_df_wind_at(dslevel)
        height_ft = get_average_height(dslevel)
        # hgt = int(math.ceil(height_ft / 1000.0))
        hgt = int(round(height_ft / 1000.0))
        hgt_str = str(hgt)
        levels_kft.append(hgt_str)

    fig = plt.figure(figsize=(8, 6))
    fig.patch.set_facecolor('lightgoldenrodyellow')

    ax1, ax2 = dual_wind_rose(fig, df1, df2, levels_kft)

    logging.info("Finished plotting wind roses at each pressure level")

    calm1 = np.sum(n_calm1)
    calm2 = np.sum(n_calm2)

    plt.suptitle(title + '\n\n', fontsize=13, fontweight='bold', y=0.88)
    # Now let's add your additional information, to remember what program options created your plot
    subtitle = "(Plotted with switch=netcdf)"
    ax1.annotate("(" + subtitle + ")",
                     xy=(0.75, 0.005), xytext=(0, 10),
                     xycoords=('axes fraction', 'figure fraction'),
                     textcoords='offset points',
                     size=8, ha='center', va='bottom')
    show_calm(ax1, ax2, calm1, calm2)

    filestr = title.replace(",", "")
    filestr = filestr.replace(" ", "_") + '.png'
    filename = output_dir + filestr
    fig.show()
    plt.savefig(filename)
    print(f"Saving an image file with name: {filename}")
    plt.close()
    return


def plt_raob(ds, switch, title, output):
    # https://github.com/marcia-marques/wind-rose
    '''
    We want to plot wind roses at 00z and 12z for all available levels. Because here we deal with winds and speeds
    at "significant" wind levels we are dealing with a variable amount of data. Significant wind levels occur when there
    is a significant change in wind (I don't know what the criteria is for significant change in wind with increasing
    altitude).
    :param ds:  xarray data for one day of model data
    :param title: title for your plot
    :param month:
    :param year:
    :return:

    '''

    ws, wd, pres, hgt = get_data_at(ds)

    # The two times (00Z and 12Z) will usually have different number of observations

    plt_no = 0
    titles =['00Z','12Z']
    df1, subtitle1, n_calm1 = get_raob_df(plt_no, ws[plt_no,:], wd[plt_no,:], pres[plt_no,:], hgt[plt_no,:])
    if len(df1) == 0:
        return

    rows = len(df1.axes[0])
    if rows == 0:
        return

    plt_no = 1
    df2, subtitle2, n_calm2 = get_raob_df(plt_no, ws[plt_no,:], wd[plt_no,:], pres[plt_no,:], hgt[plt_no,:])
    if len(df2) == 0:
        return

    rows = len(df2.axes[0])
    if rows == 0:
        return

    nbins = 16
    quadrant0_threshold = 360 - (180 / nbins)

    bins = np.arange(- (180 / nbins), 360 + (180 / nbins), 360 / nbins)

    n1, bins = np.histogram(df1.wd, bins=bins)
    n2, bins = np.histogram(df2.wd, bins=bins)

    # ns_total 2d array of 16 different wind quadrants and 5 wind speed ranges
    #ns_total = np.zeros((16, 5), dtype=np.int64)

    calm1 = np.sum(n_calm1)
    calm2 = np.sum(n_calm2)

    fig = plt.figure(figsize=(8, 6))
    fig.patch.set_facecolor('lightgoldenrodyellow')
    # https://stackoverflow.com/questions/16476924/how-to-iterate-over-rows-in-a-dataframe-in-pandas
    ax1, ax2 = dual_wind_rose_sigW(fig, df1, df2)

    logging.info("Finished plotting wind roses at each pressure level")
    plt.suptitle(title + '\n\n', fontsize=13, fontweight='bold', x=0.3, y=0.92)
    show_calm(ax1, ax2, calm1, calm2)
    fig.show()
    subtitle = "Plotted with switch="+switch
    # Now let's add your additional information, to remember what program created your plot
    ax1.annotate("(" + subtitle + ")",
        xy=(0.75, 0.002), xytext=(0, 10),
        xycoords=('axes fraction', 'figure fraction'),
        textcoords='offset points',
        size=8, ha='center', va='bottom')

    filestr = title.replace(",", "")
    filestr = filestr.replace(" ", "_") + '.png'

    # This just saves files in subdirectory img in your current working directory, it you want to use, uncomment
    #current = pathlib.Path.cwd().parent
    #current_str = current.__str__()
    #filename = current_str + "\\img\\" + filestr
    filename = output+filestr
    plt.savefig(filename)
    print(f"Saving an image file with name: {filename}")
    plt.close()



def reduce_time(ds, size):
    '''
    Reduce or get rid of the zulu times we don't want. For example, radiosondes can be launched
    at 01Z or anything other random times. Also, when storms around, often launched at 06Z and 18Z
    but we only want 00Z and 12Z radiosonde data
    :param ds:
    :param size:
    :return: ds: with only times 00Z and 12Z returned
    '''
    date_time = np.empty(size, dtype=object)
    # store the strings separate from xarray otherwise get indexing problem. Since
    # once we remove a dimension variable with our drop, our size of time decreases
    for k in range(size):
        date_time[k] = ds['time'][k].values

    hr_test = ["00", "12"]
    for k in range(size):
        dtString = str(date_time[k]).split('T')[1]
        hr = dtString[0:2]
        hr_index1 = hr.find(hr_test[0])
        hr_index2 = hr.find(hr_test[1])
        if hr_index1 == 0 or hr_index2 == 0:
            bad_hour = False
        else:
            bad_hour = True

        if bad_hour:
            ds = ds.drop_sel(time=date_time[k])

    return ds


def show_calm(ax1, ax2, n_calm1, n_calm2):
    '''
    Plot number of calm observations for each time period
    @param ax1: plot axis 1
    @param ax2: plot axis 2
    @param n_calm1: number of calm observations during 00z
    @param n_calm2: number of calm observations during 12z

    @return:
    '''


    # Stats for 00Z
    calm1 = np.sum(n_calm1)

    # Stats for 12Z
    calm2 = np.sum(n_calm2)
    calm = calm1 + calm2

    calm_label = f"Total Calm Layers = {calm}"
    calm_label1 = f"Calm Layers at 00Z = {calm1}"
    calm_label2 = f"Calm Layers at 12Z= {calm2}"

    xoff = 0.50
    yoff = 0.185
    z1 = ax2.annotate(calm_label, xy=(xoff, yoff - 0.12),
                      ha='center', va='bottom', xycoords='figure fraction', size=9)
    xoff = 0.28
    yoff = 0.24
    z2 = ax2.annotate(calm_label1, xy=(xoff, yoff - 0.12),
                      ha='center', va='bottom', xycoords='figure fraction', size=9)
    xoff = 0.7
    z3 = ax2.annotate(calm_label2, xy=(xoff, yoff - 0.12),
                      ha='center', va='bottom', xycoords='figure fraction', size=9)
    return

def main(argv):
    '''
    Either use defaults or command line arguments to run program
    @param argv:
    @return:
    '''
    global subtitle

    start_time = time.time()
    print("--- %s seconds ---" % (start_time))

    if len(argv) == 0:
        print("No arguments were passed will use default values!")

    parser = argparse.ArgumentParser()
    switch_help = ("switch defines what we want to plot: "
                  "switch=mandatory | m plot the mandatory levels of radiosonde data which almost matchs the model levels"
                  "switch=significant | s plot the significant wind changes of radiosonde data"
                  "switch=netcdf | n plot model levels from NetCDF data")
    infile_help = ("input file, if ERA5 filetype must be '.nc'"
                   "if netcdf radiosonde data, must be '.cdf'")
    parser.add_argument("-i", "--infile", help=infile_help)
    parser.add_argument("-r", "--radiosonde", help="netcdf radiosonde input file, must end in '.cdf'")
    parser.add_argument("-s","--switch", help=switch_help)
    parser.add_argument("-o","--output_dir", help="where to store the image file produced")
    parser.add_argument("-t", "--title", help="title")
    args = parser.parse_args()
    logging.info(args)
    infile = args.infile
    switch = args.switch
    output_dir = args.output_dir
    title = args.title

    # Below are default settings but the require you put our two sample data sets in a windows
    # folder called "c:\data\Hilo\"
    # If you use command line arguments you won't need the defaults.
    if switch == None:
        switch = "netcdf"
    if output_dir == None:
        output_dir = "c:\\data\\Hilo\\"
    else:
        slash = "\\"
        if output_dir.endswith(slash):
            directory = True
        else:
            output_dir = output_dir + slash
    if title == None:
        title = "Winds for Hilo"
    switch = switch[0].lower()

    if infile == None:
        if switch == 'netcdf' or switch == 'n':
            infile = "c:\\data\\Hilo\\hilo_sep2023.nc"
        else:
            infile = "c:\\data\\Hilo\\raob_hilo_sep2023.cdf"

    process_file(infile, switch, output_dir, title)
    print("--- runtime is %s seconds ---" % (time.time() - start_time))


if __name__ == "__main__":
    main(sys.argv[1:])
