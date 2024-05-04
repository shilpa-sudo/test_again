#!/usr/bin/env python
# coding: utf-8


#standard imports
#get_ipython().run_line_magic('matplotlib', 'notebook')
import netCDF4 as nc
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from netCDF4 import Dataset
from matplotlib.dates import num2date
import glob
from matplotlib import dates as mdates
import datetime
from matplotlib import animation
from datetime import date
import marineHeatWaves as mhw
from cartopy import config
import cartopy.crs as ccrs
import xarray as xr
import pandas as pd
import geopandas as gpd
import matplotlib.colors as mcolors
import os
import natsort

dtime = pd.date_range('01/25/1998','01/27/1998') #test with 1998 jan
dx = xr.open_dataset('/home/shilpa/glory_mat_analysis/glorys_quaterdeg/2018_glorys_quarterdegree.nc')
dx
lati = dx.nav_lat.isel(x=10).values
longi = dx.nav_lon.isel(y=10).values
dtimepd = pd.to_datetime(dtime)

v_extent = np.zeros((len(dtimepd),len(lati),len(longi)))

#v_extent = np.zeros((1,len(lati),len(longi)))

for lat_ind in range(len(lati[:])):
    for lon_ind in range (len(longi[:])):
        print(f'lat = {lati[lat_ind]}, lon = {longi[lon_ind]}')

        fn = (f'/home/shilpa/glory_mat_analysis/south_pacific/active_dates_at_surface_vertical_extent{lati[lat_ind]}_{longi[lon_ind]}.csv')

        if os.path.exists(fn):

            print(fn)

            df0 = pd.read_csv(fn)

            df0['surface_active_dates'] = pd.to_datetime(df0['surface_active_dates'])


            #surface_active_dates_jan_feb_march_2016 = df0[(df0.surface_active_dates.dt.strftime('%Y-%m') == '2016-01') | (df0.surface_active_dates.dt.strftime('%Y-%m') == '2016-02') | (df0.surface_active_dates.dt.strftime('%Y-%m') == '2016-03')]
            #print(surface_active_dates_jan_feb_march_2016)
            #surface_active_dates_jan_feb_march_2016['surface_active_dates'] = pd.to_datetime(surface_active_dates_jan_feb_march_2016['surface_active_dates'])
            for i in range (len(dtimepd)):

                
                surface_dates_filter = df0[df0['surface_active_dates'] == dtimepd[i] ]
                print(surface_dates_filter)
                print(len(surface_dates_filter))

                if len(surface_dates_filter) != 0:

                    v_ext = np.array(surface_dates_filter.vertical_extent)
                    print(v_ext)

                    v_extent[i,lat_ind,lon_ind] = v_ext[:]
                    #v_extent[:,lat_ind,lon_ind] = v_ext[:]



                else:
                    print('date in dt does not match date in surface active dates at this location')
        else:
            print('path does not exist')


# In[56]:


xrds = xr.Dataset(coords = dict(time = dtime,lat = lati,lon = longi),data_vars = dict(v_extent = (['time','lat','lon'],v_extent)))

xrds.to_netcdf(f'vertical_extent_south_pacific.nc')

