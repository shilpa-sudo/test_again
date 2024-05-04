#!/usr/bin/env python
# coding: utf-8

# In[1]:


#get_ipython().run_line_magic('matplotlib', 'notebook')
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from netCDF4 import Dataset
from matplotlib.dates import num2date
import glob
from matplotlib import dates as mdates
import datetime

from datetime import date

import marineHeatWaves as mhw

#from cartopy import config
#import cartopy.crs as ccrs

import xarray as xr
import pandas as pd

#import geopandas as gpd


# In[ ]:





# In[2]:


fnames = glob.glob('/home/shilpa/glory_mat_analysis/glorys_quaterdeg/*glorys*')  #""/home/shilpa/glory_mat_analysis/glorys_quaterdeg

fnames.sort()


# In[ ]:





# In[3]:


filelist = []
for i in fnames:
    filelist.append(i)


# In[4]:


filelist


# In[5]:


#filelist[2]


# In[6]:


#filelist[-152]


# In[7]:


#print(len(filelist[2:-152]))


# In[8]:


from functools import partial
def _preprocess(x):
    return x.isel(deptht=[30])  #453m

partial_func = partial(_preprocess)


# In[9]:


ds = xr.open_mfdataset(filelist[:], combine='nested',concat_dim="time_counter", preprocess= _preprocess)


# In[10]:


ds


# In[11]:


ds.deptht.values


# In[12]:


time = ds.time_counter


# In[13]:


time


# In[14]:


time.data[0]


# In[15]:


strtime = time.data[:].astype(str)


# In[16]:


pdtime = pd.to_datetime(strtime)


# In[17]:


pdtime


# In[18]:


t = np.array([date.toordinal(x) for x in pdtime])


# In[19]:


t


# In[20]:


temparr = np.array(ds.votemper[:,0,:,:]) 


# In[21]:


temparr.shape


# In[22]:


tempmean = np.nanmean(temparr)
print(tempmean)


# In[23]:


longi = ds.nav_lon.isel(y=10).values
lati = ds.nav_lat.isel(x=10).values


# In[24]:


longi.shape


# In[25]:


lati.shape


# In[ ]:


for j in range(141):
    for k in range(604):

        mhws, clim = mhw.detect(t, temparr[:,j,k])
        print("I am location ", "j =",j,",","k =",k)
        
        events = mhws['n_events']
        print(events)
        date_start = mhws['date_start'] 
        date_end = mhws['date_end'] 
        date_peak = mhws['date_peak'] 
        duration = mhws['duration'] 
        intensity_max = mhws['intensity_max']
        intensity_mean = mhws['intensity_mean'] 
        intensity_cumulative = mhws['intensity_cumulative'] 
        category = mhws['category'] 
        rateonset = mhws['rate_onset']
        ratedecline = mhws['rate_decline'] 
        lat = lati[j]
        lon = longi[k]
        
        if events == 0:
            print('no events at this location')
            
        else:
        
            d = {'latitude':lat, 'longitude':lon,'date_start': date_start, 'date_end': date_end, 'date_peak': date_peak, 
                 'duration' : duration, 'intensity_max': intensity_max, 'intensity_mean' : intensity_mean,
                 'intensity_cumulative':intensity_cumulative, 'category':category, 'rateonset': rateonset, 
                 'ratedecline':ratedecline}

            df = pd.DataFrame(data=d)
            df.to_csv('/home/shilpa/glory_mat_analysis/glorys_quarter_depth400m/all_mhws_at_%s_%s.csv'%(lat,lon))


# In[ ]:


fnames = glob.glob('/home/shilpa/glory_mat_analysis/glorys_quarter_depth400m/all_mhws_at*') 


# In[ ]:


df_concat2 = pd.concat([pd.read_csv(f) for f in fnames], ignore_index = True)
df_concat2


# In[ ]:


df_concat2.to_csv('/home/shilpa/glory_mat_analysis/glorys_quarter_depth400m/mhws_fullset_depth400m.csv')


# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:




