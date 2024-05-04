#!/usr/bin/env python
# coding: utf-8

# In[1]:


# standard imports
#%matplotlib notebook
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


# In[2]:


fnames = glob.glob('/home/shilpa/glory_mat_analysis/glorys_quaterdeg/*glorys*')  #""/home/shilpa/glory_mat_analysis/glorys_quaterdeg

fnames.sort()


# In[3]:


filelist = []
for i in fnames:
    filelist.append(i)


# In[4]:


from functools import partial
def _preprocess(x):
    return x.isel(deptht=[0])  #0.49m

partial_func = partial(_preprocess)


# In[5]:


ds = xr.open_mfdataset(filelist[:], combine='nested',concat_dim="time_counter", preprocess= _preprocess)


# In[6]:


time = ds.time_counter
strtime = time.data[:].astype(str)
pdtime = pd.to_datetime(strtime)
t = np.array([date.toordinal(x) for x in pdtime])


# In[7]:


#get our working array
temparr = np.array(ds.votemper[:,0,:,:]) #0.49m


# In[12]:


print('temparr created')
print('temparr shape',temparr.shape)


# In[9]:


#just get the lats and lons
longi = ds.nav_lon.isel(y=10).values
lati = ds.nav_lat.isel(x=10).values


# In[10]:


def f(x):
    # return math.sqrt(x)
    mhws, clim = mhw.detect(t, x)
    
    spatiallist = []
    
    events = mhws['n_events']
    if events is None or events == [] or events == () or events == 0:
        y = np.zeros(9861)
            
        spatiallist.append(y)
    
    else:
             
        index_start = mhws['index_start'] #= []
             
        index_end = mhws['index_end']# = []
             
        
        mhw_state_list = []
        z = np.zeros(9861)
                
        for i in range(events):    
            z[index_start[i]:index_end[i]]=1
            mhw_state_list.append(z)
    
        spatiallist.append(mhw_state_list[0])
    
    return spatiallist


# In[11]:


def along_axis(testarr):
    
    return np.apply_along_axis(f, 0, testarr)  #number specifies the axis


# In[ ]:


mhwspatial = along_axis(temparr[:,:,:])
print('mhwspatial created')
print('mhwspatial shape', mhwspatial.shape)


# In[ ]:


xrds = xr.Dataset(
coords = dict(time = pd.date_range("1993-01-01", periods=9861),lats = lati,longs = longi),
data_vars = dict(spextent = (['time','lats','longs'],mhwspatial)))
    
xrds.to_netcdf(f'Om_glorys_quarter_spatial_extent_pdtime.nc')

print('the deed is done for depth 0m')

