#!/usr/bin/env python
# coding: utf-8

# In[ ]:


# get spatial extent and intensity


# In[1]:


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
import mpld3
from matplotlib.ticker import FormatStrFormatter
from cartopy import config
import cartopy.crs as ccrs
from matplotlib.backends.backend_pdf import PdfPages
import xarray as xr
import pandas as pd
import matplotlib.animation as animation
import geopandas as gpd
import scipy.stats as ss
import seaborn as sb


# In[2]:


flist = ['/home/shilpa/glory_mat_analysis/noaa_oi_sst/1993.nc',
         '/home/shilpa/glory_mat_analysis/noaa_oi_sst/1994.nc',
         '/home/shilpa/glory_mat_analysis/noaa_oi_sst/1995.nc',
         '/home/shilpa/glory_mat_analysis/noaa_oi_sst/1996.nc',
         '/home/shilpa/glory_mat_analysis/noaa_oi_sst/1997.nc',
         '/home/shilpa/glory_mat_analysis/noaa_oi_sst/1998.nc',
         '/home/shilpa/glory_mat_analysis/noaa_oi_sst/1999.nc',
         '/home/shilpa/glory_mat_analysis/noaa_oi_sst/2000.nc',
         '/home/shilpa/glory_mat_analysis/noaa_oi_sst/2001.nc',
         '/home/shilpa/glory_mat_analysis/noaa_oi_sst/2002.nc',
         '/home/shilpa/glory_mat_analysis/noaa_oi_sst/2003.nc',
         '/home/shilpa/glory_mat_analysis/noaa_oi_sst/2004.nc',
         '/home/shilpa/glory_mat_analysis/noaa_oi_sst/2005.nc',
         '/home/shilpa/glory_mat_analysis/noaa_oi_sst/2006.nc',
         '/home/shilpa/glory_mat_analysis/noaa_oi_sst/2007.nc',
         '/home/shilpa/glory_mat_analysis/noaa_oi_sst/2008.nc',
         '/home/shilpa/glory_mat_analysis/noaa_oi_sst/2009.nc',
         '/home/shilpa/glory_mat_analysis/noaa_oi_sst/2010.nc',
         '/home/shilpa/glory_mat_analysis/noaa_oi_sst/2011.nc',
         '/home/shilpa/glory_mat_analysis/noaa_oi_sst/2012.nc',
         '/home/shilpa/glory_mat_analysis/noaa_oi_sst/2013.nc',
         '/home/shilpa/glory_mat_analysis/noaa_oi_sst/2014.nc',
         '/home/shilpa/glory_mat_analysis/noaa_oi_sst/2015.nc',
         '/home/shilpa/glory_mat_analysis/noaa_oi_sst/2016.nc',
         '/home/shilpa/glory_mat_analysis/noaa_oi_sst/2017.nc',
         '/home/shilpa/glory_mat_analysis/noaa_oi_sst/2018.nc',
         '/home/shilpa/glory_mat_analysis/noaa_oi_sst/2019.nc']


# In[3]:


#ds = xr.open_mfdataset(flist)

from functools import partial

def _preprocess(x, lat_bnds):
    return x.sel(lat=slice(*lat_bnds))
lat_bnds = (-39.875, 39.875)

partial_func = partial(_preprocess, lat_bnds=lat_bnds)

ds = xr.open_mfdataset(flist[:],combine='nested', concat_dim="time",preprocess=partial_func)



# In[4]:


ds


# In[5]:


t = [pd.to_datetime(x).toordinal() for x in ds.time.values]
#print(t)


# In[6]:


#print(len(t))


# In[7]:


#print(ds.sst.shape)

sstar = np.array(ds.sst[:,:,:])
# In[ ]:


intensitylist = []
mhw_state_spatial_list = []
for j in range(320):
    print(j)
    for k in range(1440):
        print(k)
        #sstraw = ds.sst.isel(lat=j,lon = k)
        sst = sstar[:,j,k]
        mhws, clim = mhw.detect(t, sst)
        events = mhws['n_events']
        print(events)
        if events is None or events == [] or events == () or events == 0:
             y = np.zeros(9861)
             mhw_state_spatial_list.append(y)
             intensitylist.append(y)
        else:
             print(events)
             index_start = mhws['index_start'] #= []
             print(index_start)
             index_end = mhws['index_end']# = []
             print(index_end)
        
             mhw_state_list = []
             mhw_intensity_list = []
             x = np.zeros(9861)
             xx = np.zeros(9861)
                
             for i in range(events):    
                x[index_start[i]:index_end[i]]=1
                mhw_state_list.append(x)
                xx[index_start[i]:index_end[i]]= ((sst[index_start[i]:index_end[i]]) -  (clim['seas'][index_start[i]:index_end[i]]))
              
                mhw_intensity_list.append(xx)
               
               
             mhw_state_spatial_list.append(mhw_state_list[0])  
            
             intensitylist.append(mhw_intensity_list[0])


# In[ ]:


spatialextent2D = np.array(maxverticalextentlist).reshape(720,1440)
intensity2D = np.array(intensitylist).reshape(720,1440)


# In[ ]:


xrds = xr.Dataset(
       coords = dict(latitude = ds.lat,longitude = ds.lon),
       data_vars = dict(spatial2D = (['latitude','longitude'],spatialextent2D),
                       intensity2D = (['latitude','longitude'],intensity2D)))


# In[ ]:


xrds.to_netcdf('intensity_spatialextent2D_noaa.nc')

