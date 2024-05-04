#!/usr/bin/env python
# coding: utf-8


#standard imports
#get_ipython().run_line_magic('matplotlib', 'notebook')
import netCDF4 as nc
import numpy as np
#import matplotlib as mpl

from netCDF4 import Dataset

import glob

import xarray as xr
import pandas as pd
#import geopandas as gpd
#import matplotlib.colors as mcolors
import os
#import natsort

dtime = pd.date_range('01/01/1993','12/31/2019')
dx = xr.open_dataset('/home2/datawork/slal/glorys_quaterdeg/1993_glorys_quarterdegree.nc')
dx
lati = dx.nav_lat.isel(x=10).values
longi = dx.nav_lon.isel(y=10).values
dtimepd = pd.to_datetime(dtime)

#v_extent = np.zeros((len(dtimepd),len(lati),len(longi)))

for tx in range (2000,5000):

    v_extent = np.zeros((1,len(lati),len(longi)))

    for lat_ind in range(len(lati[:])):
        for lon_ind in range (len(longi[:])):
            print(f'lat = {lati[lat_ind]}, lon = {longi[lon_ind]}')

            fn = (f'/home2/datawork/slal/south_pacific/active_dates_at_surface_vertical_extent{lati[lat_ind]}_{longi[lon_ind]}.csv')

            if os.path.exists(fn):

                print(fn)

                df0 = pd.read_csv(fn)

                df0['surface_active_dates'] = pd.to_datetime(df0['surface_active_dates'])
                surface_dates_filter = df0[df0['surface_active_dates'] == dtimepd[tx] ]
                print(surface_dates_filter)
                print(len(surface_dates_filter))

                if len(surface_dates_filter) != 0:

                    v_ext = np.array(surface_dates_filter.vertical_extent)
                    print(v_ext)

                    v_extent[:,lat_ind,lon_ind] = v_ext[:]
                    #v_extent[:,lat_ind,lon_ind] = v_ext[:]



                else:
                    print('date in dt does not match date in surface active dates at this location')
            else:
                print('path does not exist')


# In[56]:


    xrds = xr.Dataset(coords = dict(time = [dtime[tx]],lat = lati,lon = longi),data_vars = dict(v_extent = (['time','lat','lon'],v_extent)))

    xrds.to_netcdf(f'/home2/datawork/slal/daily_v_ext_glorys_quarter/daily_vertical_extent_south_pacific_{tx}.nc')
