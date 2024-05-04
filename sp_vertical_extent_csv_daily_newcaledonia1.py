#!/usr/bin/env python
# coding: utf-8


#standard imports
#get_ipython().run_line_magic('matplotlib', 'notebook')
import netCDF4 as nc
#import numpy as np
#import matplotlib as mpl

from netCDF4 import Dataset

#import glob

#import xarray as xr
#import pandas as pd
#import geopandas as gpd
#import matplotlib.colors as mcolors
import os
#import natsort

import numpy as np
from cartopy import config
import cartopy.crs as ccrs
import xarray as xr
import pandas as pd
import geopandas as gpd
import glob as glob
from shapely.geometry import Polygon




dtime = pd.date_range('01/01/1993','12/31/2019')
#dx = xr.open_dataset('/home2/datawork/slal/glorys_quaterdeg/1993_glorys_quarterdegree.nc')
#dx
latis = dx.nav_lat.isel(x=10).values
longis = dx.nav_lon.isel(y=10).values
lati = np.array(latis[30:85])

longi = np.array(longis[59:129])


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



    eez = gpd.read_file("/home2/datahome/slal/World_EEZ_v11_20191118/eez_v11.shp") # get the EEZs


    ncd = eez[eez["TERRITORY1"] == "New Caledonia"]


    #df = pd.read_csv('/home2/datahome/slal/lonlist_clip_eez.csv')
    #dfarrno = df.to_numpy()
#dfarrno.shape
    #nlon = dfarrno[:,1]

    def country_mean_vertical_extent(country):
        edgar_gdf_vertical_extent.crs = country.crs
        country_data = gpd.clip(edgar_gdf_vertical_extent, country)
        country_data_filtered = country_data[(country_data['v_extent'] != 0)]
        mean_vertical_extent_country = (np.nanmean(country_data_filtered['v_extent']))
        return mean_v_extent_country

    def country_mean_plus2sd_v_extent(country):
        edgar_gdf_vertical_extent.crs = country.crs
        country_data = gpd.clip(edgar_gdf_vertical_extent, country)
        country_data_filtered = country_data[(country_data['v_extent'] != 0)]
        mean_vertical_extent_country = (np.nanmean(country_data_filtered['v_extent']))
        std = (np.nanstd(country_data_filtered['v_extent']))
        upper_bound = mean_vertical_extent_country + 2 * std
        return upper_bound


# In[ ]:


    def country_mean_minus2sd_v_extent(country):
        edgar_gdf_vertical_extent.crs = country.crs
        country_data = gpd.clip(edgar_gdf_vertical_extent, country)
        country_data_filtered = country_data[(country_data['v_extent'] != 0)]
        mean_vertical_extent_country = (np.nanmean(country_data_filtered['v_extent']))
        std = (np.nanstd(country_data_filtered['v_extent']))
        lower_bound = mean_vertical_extent_country - 2 * std
        return lower_bound




    edgar_v_extent = xrds.v_extent[:,:,:].to_dataframe()
    #edgar_spatial['newlons'] = nlon
    #edgar_spatial = edgar_spatial.reset_index()
    edgar_gdf_v_extent = gpd.GeoDataFrame(edgar_v_extent, geometry=gpd.points_from_xy(edgar_v_extent.lon, edgar_v_extent.lat))

    ncd_mean_v_extent = country_mean_vertical_extent(ncd)
    ncd_mean_plus2sd_v_extent = country_mean_plus2sd_v_extent(ncd)
    ncd_mean_minus2sd_v_extent = country_mean_minus2sd_v_extent(ncd)
    
    d = {'time_index': tx, 'time' : dtime[tx],
        'New_Caledonia_mean_vertical_extent':ncd_mean_v_extent,'New_Caledonia_mean_plus2sd_v_extent':ncd_mean_plus2sd_v_extent,
        'New_Caledonia_mean_minus2sd_intensity': ncd_mean_minus2sd_v_extent}

    
    dfnew = pd.DataFrame(data=d,index=[0])
    dfnew.to_csv(f"/home2/datawork/slal/new_caledonia_vertical_extent_daily/new_caledonia_vertical_extent_{tx}.csv")






