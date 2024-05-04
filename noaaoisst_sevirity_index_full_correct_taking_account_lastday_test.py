#!/usr/bin/env python
# coding: utf-8

# In[1]:


# standard imports

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
import folium
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
from matplotlib.patches import Rectangle


# In[2]:


dtf = xr.open_dataset('/home/shilpa/glory_mat_analysis/noaaoisst93_19_sst_AOI_quarterdegree.nc',decode_times=False) 
lati= dtf.lati #degrees north
longi = dtf.longi #degrees east
time = dtf.time.values
temp = dtf.temp


# In[3]:


temparr = np.array(temp)


# In[4]:


tt = np.array(time)
t = tt.astype(int)
datessfull = [date.fromordinal(ttt.astype(int)) for ttt in t]


# In[5]:


datessfull


# In[6]:


dtf


# In[7]:


mhw_severity_index = np.zeros((9861,130,260))


# In[8]:


mhw_severity_index.shape


# In[17]:


for j in range(130):
    for k in range(260):
        print(f'j={j},k={k}')
        
        sst = temparr[:,j,k]
        mhws, clim = mhw.detect(t, sst)
        
        events = mhws['n_events']
        print(f'events = {events}')
        if events is None or events == [] or events == () or events == 0:
             print('no events at this location')
            
        else:
             print(events)
             index_start = mhws['index_start'] #= []
             print(index_start)
             index_end = mhws['index_end']# = []
             print(index_end)
             
             for i in range(events):  

                if index_end[i] != 9861: #for all events not ending on the last day that is 31 dec 2019

                    mhw_severity_index[index_start[i]:(index_end[i])+1,j,k] = ((sst[index_start[i]:(index_end[i])+1]) -  (clim['seas'][index_start[i]:(index_end[i])+1]))/ ((clim['thresh'][index_start[i]:(index_end[i])+1]) -  (clim['seas'][index_start[i]:(index_end[i])+1])) # this will ensure that the last day in mhw is included

                    
                #elif index_end[i] == 9861:                       #for all events ending on the last day that is 31 dec 2019
                   
                    #mhw_intensity_spatial[j,k,index_start[i]:index_end[i]] = ((sst[index_start[i]:index_end[i]]) -  (clim['seas'][index_start[i]:index_end[i]]))
                  
                    


# In[18]:


date_range = pd.date_range(start='1993-01-01', end='2019-12-31')


# In[19]:


len(date_range)


# In[20]:


longi = np.array(dtf.longi.values)
lati = np.array(dtf.lati.values)


# In[21]:


longi.shape


# In[25]:


xrds = xr.Dataset(
       coords = dict(time = date_range,lons = longi,lats = lati),
       data_vars = dict(
                   spatial_intensity = (['time','lats','lons'],mhw_severity_index)))
                
                   


# In[27]:


xrds.to_netcdf('full_daily_severity_index_noaaoisst_aoi_lastdayincluded.nc')


# In[28]:


s2D2 = xr.open_dataset('full_daily_severity_index_noaaoisst_aoi_lastdayincluded.nc')
s2D2


# In[ ]:




