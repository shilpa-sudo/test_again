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






dtime = pd.date_range('01/01/1993','12/31/2019')
# In[10]:

def fsorted2 (lat,lon):

    fsorted2 = [(f'/home/shilpa/glory_mat_analysis/glorys_quarter_depth0/all_mhws_at_{lat}_{lon}.csv'),
 (f'/home/shilpa/glory_mat_analysis/glorys_quarter_depth1m/all_mhws_at_{lat}_{lon}.csv'),
 (f'/home/shilpa/glory_mat_analysis/glorys_quarter_depth2m/all_mhws_at_{lat}_{lon}.csv'),
 (f'/home/shilpa/glory_mat_analysis/glorys_quarter_depth3m/all_mhws_at_{lat}_{lon}.csv'),
 (f'/home/shilpa/glory_mat_analysis/glorys_quarter_depth5m/all_mhws_at_{lat}_{lon}.csv'),
 (f'/home/shilpa/glory_mat_analysis/glorys_quarter_depth6m/all_mhws_at_{lat}_{lon}.csv'),
 (f'/home/shilpa/glory_mat_analysis/glorys_quarter_depth7m/all_mhws_at_{lat}_{lon}.csv'),
 (f'/home/shilpa/glory_mat_analysis/glorys_quarter_depth10m/all_mhws_at_{lat}_{lon}.csv'),
 (f'/home/shilpa/glory_mat_analysis/glorys_quarter_depth11m/all_mhws_at_{lat}_{lon}.csv'),
 (f'/home/shilpa/glory_mat_analysis/glorys_quarter_depth13m/all_mhws_at_{lat}_{lon}.csv'),
 (f'/home/shilpa/glory_mat_analysis/glorys_quarter_depth15m/all_mhws_at_{lat}_{lon}.csv'),
 (f'/home/shilpa/glory_mat_analysis/glorys_quarter_depth18m/all_mhws_at_{lat}_{lon}.csv'),
 (f'/home/shilpa/glory_mat_analysis/glorys_quarter_depth21m/all_mhws_at_{lat}_{lon}.csv'),
 (f'/home/shilpa/glory_mat_analysis/glorys_quarter_depth25m/all_mhws_at_{lat}_{lon}.csv'),
 (f'/home/shilpa/glory_mat_analysis/glorys_quarter_depth29m/all_mhws_at_{lat}_{lon}.csv'),
 (f'/home/shilpa/glory_mat_analysis/glorys_quarter_depth34m/all_mhws_at_{lat}_{lon}.csv'),
 (f'/home/shilpa/glory_mat_analysis/glorys_quarter_depth40m/all_mhws_at_{lat}_{lon}.csv'),
 (f'/home/shilpa/glory_mat_analysis/glorys_quarter_depth47m/all_mhws_at_{lat}_{lon}.csv'),
 (f'/home/shilpa/glory_mat_analysis/glorys_quarter_depth55m/all_mhws_at_{lat}_{lon}.csv'),
 (f'/home/shilpa/glory_mat_analysis/glorys_quarter_depth65m/all_mhws_at_{lat}_{lon}.csv'),
 (f'/home/shilpa/glory_mat_analysis/glorys_quarter_depth77m/all_mhws_at_{lat}_{lon}.csv'),
 (f'/home/shilpa/glory_mat_analysis/glorys_quarter_depth92m/all_mhws_at_{lat}_{lon}.csv'),
 (f'/home/shilpa/glory_mat_analysis/glorys_quarter_depth100m/all_mhws_at_{lat}_{lon}.csv'),
 (f'/home/shilpa/glory_mat_analysis/glorys_quarter_depth130m/all_mhws_at_{lat}_{lon}.csv'),
 (f'/home/shilpa/glory_mat_analysis/glorys_quarter_depth155m/all_mhws_at_{lat}_{lon}.csv'),
 (f'/home/shilpa/glory_mat_analysis/glorys_quarter_depth186m/all_mhws_at_{lat}_{lon}.csv'),
 (f'/home/shilpa/glory_mat_analysis/glorys_quarter_depth222m/all_mhws_at_{lat}_{lon}.csv'),
 (f'/home/shilpa/glory_mat_analysis/glorys_quarter_depth266m/all_mhws_at_{lat}_{lon}.csv'),
 (f'/home/shilpa/glory_mat_analysis/glorys_quarter_depth318m/all_mhws_at_{lat}_{lon}.csv'),
 (f'/home/shilpa/glory_mat_analysis/glorys_quarter_depth380m/all_mhws_at_{lat}_{lon}.csv'),
 (f'/home/shilpa/glory_mat_analysis/glorys_quarter_depth453m/all_mhws_at_{lat}_{lon}.csv'),
 (f'/home/shilpa/glory_mat_analysis/glorys_quarter_depth500m/all_mhws_at_{lat}_{lon}.csv'),
 (f'/home/shilpa/glory_mat_analysis/glorys_quarter_depth600m/all_mhws_at_{lat}_{lon}.csv'),
 (f'/home/shilpa/glory_mat_analysis/glorys_quarter_depth700m/all_mhws_at_{lat}_{lon}.csv'),
 (f'/home/shilpa/glory_mat_analysis/glorys_quarter_depth900m/all_mhws_at_{lat}_{lon}.csv'),
 (f'/home/shilpa/glory_mat_analysis/glorys_quarter_depth1000m/all_mhws_at_{lat}_{lon}.csv'),
 (f'/home/shilpa/glory_mat_analysis/glorys_quarter_depth1200m/all_mhws_at_{lat}_{lon}.csv'),
 (f'/home/shilpa/glory_mat_analysis/glorys_quarter_depth1500m/all_mhws_at_{lat}_{lon}.csv')]
    return fsorted2
# In[ ]:


dx = xr.open_dataset('/home/shilpa/glory_mat_analysis/glorys_quaterdeg/2018_glorys_quarterdegree.nc')
dx

lati = dx.nav_lat.isel(x=10).values
longi = dx.nav_lon.isel(y=10).values
#print(len(lati))
#print(len(longi))


x = (dx.deptht.values)
xlist = list(x)
depth = xlist

rdlist = []
for i in range(len(depth)):

    rd = round(depth[i],2)

    rdlist.append(rd)

print(rdlist)


# In[ ]:


def get_active_dates_at_surface(df0):

    active_date_range_depth0 = []

    for i in range(len(df0)):
        trange = pd.date_range(start=(df0['date_start']).iloc[i], end=(df0['date_end']).iloc[i])
        active_date_range_depth0.append(trange)

    active_date_range_depth0f = [element for sublist in active_date_range_depth0 for element in sublist]

    return active_date_range_depth0f

# In[ ]:


def get_active_dates_at_eachdepth_intersects_activedates_at_surface(df):

    active_date_range_depth = []

    for i in range(len(df)):
        trange = pd.date_range(start=(df['date_start']).iloc[i], end=(df['date_end']).iloc[i])
        active_date_range_depth.append(trange)

    active_date_range_depthf = [element for sublist in active_date_range_depth for element in sublist]
    xx = set(active_date_range_depth0f).intersection(set(active_date_range_depthf))#check if dates at surface are present at this depth
    samedatelist = list(xx)

    return samedatelist


# In[ ]:


def depth_extent_for_eachday1 (active_date_range_depth0f,surface_depth_extent_list):
    depth_extend_for_eachday = []

    for d in range(len(active_date_range_depth0f)): #for each day at the surface

        depthextend = []
        for i in range(38):
                        #get depth extent of surface events, break otherwise
            if surface_depth_extent_list[i] == [] : #none of the dates are common to surface dates
                break
            elif active_date_range_depth0f[d] in list(surface_depth_extent_list[i]):
                depthextend.append(rdlist[i])
            else:
                break
        depth_extend_for_eachday.append(depthextend)

    return depth_extend_for_eachday



for lat in lati[:]:
    for lon in longi[:]:

        fn = (f'/home/shilpa/glory_mat_analysis/glorys_quarter_depth0/all_mhws_at_{lat}_{lon}.csv')

        if os.path.exists(fn):

            print(fn)

            df0 =  pd.read_csv(fn)
            #df0['date_peak'] = pd.to_datetime(ncd_data['date_peak'])

            active_date_range_depth0f  = get_active_dates_at_surface(df0)
            #dtf_active_surface_dates = pd.DataFrame(active_date_range_depth0f, columns=['surface_active_dates'])
            #df.to_csv(f'active_dates_at_surface_{lat}_{lon}.csv', index=False)

            #print(len(active_date_range_depth0f))
            #lensurfcaeactivedates = len(active_date_range_depth0f)
            #sum_surface_activedateslist.append(lensurfcaeactivedates)

            fsortedlist = fsorted2 (lat,lon)

            surface_depth_extent_list = []

            for fnx in fsortedlist:
                if os.path.exists(fnx):
                    df = pd.read_csv(fnx)
                    samedates_to_surface = get_active_dates_at_eachdepth_intersects_activedates_at_surface(df)
                    surface_depth_extent_list.append(samedates_to_surface)
                else:
                    surface_depth_extent_list.append([])

            print(len(active_date_range_depth0f))
            print(len(surface_depth_extent_list))

            texet = depth_extent_for_eachday1(active_date_range_depth0f,surface_depth_extent_list)

	        #print(len(texet))

            maxdepthtimeseries = []
            for i in range(len(active_date_range_depth0f)):
                #print(len(active_date_range_depth0f)
                #print(len(depth_extent_for_eachday))  

                print(texet[i])

                maxdepthextent_eachday = np.nanmax(texet[i])

                print(maxdepthextent_eachday)

                maxdepthtimeseries.append(maxdepthextent_eachday)

            dtf_active_surface_dates_maxv_extent = pd.DataFrame({'surface_active_dates': active_date_range_depth0f,'vertical_extent': maxdepthtimeseries})
            dtf_active_surface_dates_maxv_extent.to_csv(f'/home/shilpa/glory_mat_analysis/south_pacific/active_dates_at_surface_vertical_extent{lat}_{lon}.csv', index=False)



        else:
            print(f'surface file name does not exist at lat {lat} and lon {lon} ')

            #sum_surface_activedateslist.append(np.nan)

