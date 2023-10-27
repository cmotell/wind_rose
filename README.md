# Introduction

In this repository we show start building a gallary of Python programs
that look at modeled stratospheric winds and how well modeled winds
compare to radiosonde observations. 

Note, we start with NOAA archived radiosonde observations and ERA5 Reanalysis winds.
In the Appendix of this readme is a discussion on how to obtain these data.

Also, attached to this repository are two sample netcdf files containing
radiosonde data (file type of cdf), and ERA5 reanalysis data (file type of nc).

All datasources discussed here will be netcdf (or shown much latter, how)
to convert Grib files to Netcdf), and they will contain the U and V wind components 
as well as either Geopotential (from which Geopotential Height is derived), or
actual heights of the U and V wind components.

# Acknowledgements

Often in my code, you will see a URL link where I got the hint to solve
what was for me, a more difficult coding solution. 

Here, we are going to be looking at a variety of python
program techniques to display stratospheric winds. In the first example, 
I highlight one of my wind rose solutions. The code for wind roses came from
a combination of what I saw others do, modified for my needs. Unfortunately,
I don't remember where I got some of my wind rose ideas. But I do know I liked
the work of https://github.com/marcia-marques/wind-rose; however, I needed to
modify her work to suite my needs.

## Example 1: demonstrations the use of wind roses 

Run the program *wind_rose.py* to show wind roses. Let's say we have two following
datasets: *hilo_sep2023.nc*, *raob_hilo_sep2023.cdf*.

The wind rose program has the following command switches:

<ul>
    <li>-i or --infile which is the full path to your data file</li>
    <li>-o or --output_dir which is full path to your output directory</li>
    <li>-s or --switch which is your options: 
        <ul><li>n or netcdf for ERA5 reanalysis model data</li>
            <li>m or mandatory to display radiosonde mandatory pressure levels</li>
            <li>s or signficant to display significant wind level changes</li>
        </ul></li>
    <li>-t or --title which is the title to your program</li>
</ul>


Note, we also can do default runs assuming that you place our two sample data
files in the directory:

c:\data\hilo

### 1a. Model Data

In the first example lets look at ERA5 modeled winds, in the form of a wind rose,
between 150mb and 10mb. Run *wind_rose.py* with following command parameters:

<img src="./img/netcdf_conf.png"/>

This will produce the following image:

<img src="./img/Hilo_ERA5_Wind_for_Sep_2023_Day_1.png"/>

### 1b. Radiosonde at the Mandatory Pressure Levels

Next, run *wind_rose.py* with the following command parameters to see the mandatory
pressure level radiosonde observations.

<img src="./img/mandatory_conf.png"/>

resulting in

<img src="./img/Hilo_Radiosonde_Wind_on_Sep_2023_Day_1.png"/>

We see only 6 model pressure levels between approximately 46,000 to 87,000 ft.
Also, we didn't get any observations above the 20mb level which occurs frequently.

Lastly, modify the configuration as follows:

<img src="./img/sigw_conf.png"/>

Resulting in the following image:

<img src="./img/Hilo_Radiosonde_sigWinds_on_Sep_2023_Day_1.png"/>

## Example 2: demonstration of a time series of upper wind data for one month

(Working on ... coming soon) Here in this example site we will to a cross-sectional time series for specified 
location for a given month. Here we make the assumption that all data consist of
a single month of data (could be easily modified using techique in Example 1
to display as many images as there our months of data)

# Appendix


## A.1 Installing the Software Packages

Here we show the steps I used to install the Python packages. 
First create a conda environment:

<img src="./img/conda_create.png"/>

Next, install your base of xarray and NetCDF:

<img src="./img/conda_install_xarray.png"/>

## A.1 Accessing data from copernicus climate store:

First select the fields and heights as shown below (for both regional and site studies):

<img src="./img/copernicus_base_pull.png"/>

Next, select the times, here we selected only 00z and 12z to correspond to normal
radiosonde data (note, for some weather events such as hurricanes near a site, more
frequest radiosonde is often available)

<img src="./img/copernicus_base_pull.png"/>

Finally, select a latitude, and longitude to match your radiosonde data. Here, we
selected the latitude and longitude of the Hilo, Hawaii radiosonde site

<img src="./img/copernicus_hilo_pull.png"/>

Note, when selecting a point, you must have a very slight difference and your north-south, 
east-west locations. I believe the model uses nearest-neighbor to select the data from a
grid-point (should confirm).

## Accessing radiosonde data

The data for this comes from https://ruc.noaa.gov/raobs/

To grab data in NetCDF format select as follows:

<img src="./img/sound_pull.png"/>

Note, don't worry about the units are anything here. Netcdf data has a standard format and 
doesn't follow what this first pages say.

Select continue

<img src="./img/sound_pull_netcdf.png"/>

Now search for the site you want, but selecting continue on this page and
selecting your site, like Hilo as shown below:

<img src="./img/sound_hilo.png"/>

The filetype produced is "cdf" but actually "nc", and this file will contain
the data we want to analyze. Note, if you want to see some of the contents in
radiosonde data, select a text format and visually inspect the data.**





