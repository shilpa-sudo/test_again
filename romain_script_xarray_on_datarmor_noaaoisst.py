#!/usr/bin/env python
# coding: utf-8

# In[2]:


import numpy as np
#import matplotlib as mpl
#import matplotlib.pyplot as plt
#from matplotlib.dates import num2date
#from matplotlib import dates as mdates
import xarray as xr
import pandas as pd
import marineHeatWaves as mhw

import glob #is this already installed in virtual envt?

import datetime #is this already installed in virtual envt?

from datetime import date

from functools import partial #is this already installed in virtual envt?


# In[2]:


fnames = glob.glob('/home/datawork-lead/datarmor-only/shilpa/noaaoisst/*.nc')  
fnames.sort()


# In[3]:


filelist = []
for i in fnames:
    filelist.append(i)


# In[3]:


# get files from 1993 to 2019


# In[11]:


filelist[11:]
#filelist[:2]


# In[24]:


def _preprocess(x, lon_bnds, lat_bnds):

    return x.sel(lon=slice(*lon_bnds), lat=slice(*lat_bnds))


lon_bnds, lat_bnds = (145.125, 291.875), (-37.875,-2.125)

partial_func = partial(_preprocess, lon_bnds=lon_bnds, lat_bnds=lat_bnds)

ds = xr.open_mfdataset(filelist[:2], combine='nested',concat_dim="time",preprocess=partial_func)  


# In[25]:


ds


# In[ ]:


time = ds.time
strtime = time.data[:].astype(str)
pdtime = pd.to_datetime(strtime)


# In[ ]:


t = np.array([date.toordinal(x) for x in pdtime])


# In[ ]:


temparr = np.array(ds.sst[:,:,:]) 


# In[ ]:


spatial_extent = np.zeros((9861, 144, 588))
date_start_arr = np.zeros((9861, 144, 588))
date_end_arr = np.zeros((9861, 144, 588))
date_peak_arr = np.zeros((9861, 144, 588))
index_start_arr = np.zeros((9861, 144, 588))
index_end_arr = np.zeros((9861, 144, 588))
index_peak_arr = np.zeros((9861, 144, 588))
duration_arr = np.zeros((9861, 144, 588))
intensity_max_arr = np.zeros((9861, 144, 588))
intensity_mean_arr = np.zeros((9861, 144, 588)) 
intensity_cumulative_arr = np.zeros((9861, 144, 588))
category_arr = np.zeros((9861, 144, 588))
rateonset_arr = np.zeros((9861, 144, 588))
ratedecline_arr = np.zeros((9861, 144, 588))

intensity_t_minus_clim_arr = np.zeros((9861, 144, 588)) 
intensity_t_minus_thresh_arr = np.zeros((9861, 144, 588))

reaction_window_days_arr = np.zeros((9861, 144, 588))
coping_window_days_arr = np.zeros((9861, 144, 588))
recovery_window_days_arr = np.zeros((9861, 144, 588))

severity_index = np.zeros((9861, 144, 588))


# In[ ]:


for j in range(144):
    for k in range(588):

        mhws, clim = mhw.detect(t, temparr[:,j,k])
        print("I am location ", "j =",j,",","k =",k)
        
        events = mhws['n_events']
        print(events)
        date_start = mhws['date_start'] 
        date_end = mhws['date_end'] 
        date_peak = mhws['date_peak']
        index_start = mhws['index_start']
        index_end = mhws['index_end']
        index_peak = mhws['index_peak']
        duration = mhws['duration'] 
        intensity_max = mhws['intensity_max']
        intensity_mean = mhws['intensity_mean'] 
        intensity_cumulative = mhws['intensity_cumulative'] 
        category = mhws['category'] 
        rateonset = mhws['rate_onset']
        ratedecline = mhws['rate_decline']  
        
        intensity_t_minus_clim = temparr[:,j,k] - clim['seas'][:,j,k]  
        intensity_t_minus_thresh = temparr[:,j,k] - clim['thresh'][:,j,k] 
        
        intensity_t_minus_clim_arr[:,j,k] = intensity_t_minus_clim       # array is filled here
        intensity_t_minus_thresh_arr[:,j,k] = intensity_t_minus_thresh   # array is filled here
        
        severity_index[:,j,k] = intensity_t_minus_clim / ((clim['thresh'][:,j,k] - clim['seas'][:,j,k]))

        
        
        #fill array below
        for i in range(events):
            if index_end[i] != 9860: #if events are not ending on the last day
                spatial_extent[index_start[i]:(index_end[i])+1,j,k] = 1
                date_start_arr[index_start[i]:(index_end[i])+1,j,k] =  date_start[i]
                date_end_arr[index_start[i]:(index_end[i])+1,j,k] =  date_end[i]
                date_peak_arr[index_start[i]:(index_end[i])+1,j,k] =  date_peak[i] 
                index_start_arr[index_start[i]:(index_end[i])+1,j,k] =  index_start[i]
                index_end_arr[index_start[i]:(index_end[i])+1,j,k] =  index_end[i]
                index_peak_arr[index_start[i]:(index_end[i])+1,j,k] =  index_peak[i]
                duration_arr[index_start[i]:(index_end[i])+1,j,k] =  duration[i]
                intensity_max_arr[index_start[i]:(index_end[i])+1,j,k] =  intensity_max[i]
                intensity_mean_arr[index_start[i]:(index_end[i])+1,j,k] =  intensity_mean[i]
                intensity_cumulative_arr[index_start[i]:(index_end[i])+1,j,k] =  intensity_cumulative[i]
                category_arr[index_start[i]:(index_end[i])+1,j,k] =  category[i]
                rateonset_arr[index_start[i]:(index_end[i])+1,j,k] =  rateonset[i]
                ratedecline_arr[index_start[i]:(index_end[i])+1,j,k] =  ratedecline[i]
                
                reaction_window_days_arr[index_start[i]:(index_end[i])+1,j,k] = ((index_peak[i]) - (index_start[i]))
                coping_window_days_arr[index_start[i]:(index_end[i])+1,j,k]   = ((index_end[i]) - (index_peak[i]))

            #else: #if events are ending on the last day
                #spatial_extent[index_start[i]:index_end[i],j,k] = 1
                #date_start_arr[index_start[i]:index_end[i],j,k] =  date_start[i]
                #date_end_arr[index_start[i]:index_end[i],j,k] =  date_end[i]
                #date_peak_arr[index_start[i]:index_end[i],j,k] =  date_peak[i] 
                #index_start_arr[index_start[i]:index_end[i],j,k] =  index_start[i]
                #index_end_arr[index_start[i]:index_end[i],j,k] =  index_end[i]
                #index_peak_arr[index_start[i]:index_end[i],j,k] =  index_peak[i]
                #duration_arr[index_start[i]:index_end[i],j,k] =  duration[i]
                #intensity_max_arr[index_start[i]:index_end[i],j,k] =  intensity_max[i]
                #intensity_mean_arr[index_start[i]:index_end[i],j,k] =  intensity_mean[i]
                #intensity_cumulative_arr[index_start[i]:index_end[i],j,k] =  intensity_cumulative[i]
                #category_arr[index_start[i]:index_end[i],j,k] =  category[i]
                #rateonset_arr[index_start[i]:index_end[i],j,k] =  rateonset[i]
                #ratedecline_arr[index_start[i]:index_end[i],j,k] =  ratedecline[i]
                
                #reaction_window_days_arr[index_start[i]:index_end[i],j,k] = ((index_peak[i]) - (index_start[i]))
                #coping_window_days_arr[index_start[i]:index_end[i],j,k]   = ((index_end[i]) - (index_peak[i]))

        # special case for recovery window days
        for i in range(events-1):
                if index_end[i] != 9860: #if events are not ending on the last day
                    recovery_window_days_arr[index_start[i]:(index_end[i])+1,j,k] = ((index_start[i+1])-(index_end[i]))
                    
                #else:#if events are ending on the last day
                    #recovery_window_days_arr[index_start[i]:index_end[i],j,k] = ((index_start[i+1])-(index_end[i]))
                    
                    
                    


# In[28]:


longi = np.array(ds.lon.values)
lati = np.array(ds.lat.values)
time = time


# In[ ]:


xrds = xr.Dataset(
       coords = dict(time = time,lon = longi,lat = lati),
       data_vars = dict(
                   spatial_extent = (['time','lat','lon',], spatial_extent),
                   date_start_arr = (['time','lat','lon',], date_start_arr),    
                   date_end_arr = (['time','lat','lon',], date_end_arr),  
                   date_peak_arr = (['time','lat','lon',], date_peak_arr),  
                   
                   index_start_arr = (['time','lat','lon',], index_start_arr),  
                    
                   index_end_arr = (['time','lat','lon',], index_end_arr),  
                   index_peak_arr = (['time','lat','lon',], index_peak_arr),  
                   duration_arr = (['time','lat','lon',], duration_arr),  
                   intensity_max_arr = (['time','lat','lon',], intensity_max_arr),  
                   intensity_mean_arr = (['time','lat','lon',], intensity_mean_arr),  
                   intensity_cumulative_arr = (['time','lat','lon',], intensity_cumulative_arr),  
                   category_arr = (['time','lat','lon',], category_arr),  
                   rateonset_arr = (['time','lat','lon',],  rateonset_arr),  
                   ratedecline_arr = (['time','lat','lon',], ratedecline_arr),  
                   intensity_t_minus_clim_arr = (['time','lat','lon',], intensity_t_minus_clim_arr),  
                   intensity_t_minus_thresh_arr = (['time','lat','lon',], intensity_t_minus_thresh_arr),  

                   reaction_window_days_arr = (['time','lat','lon',], reaction_window_days_arr),  
                   coping_window_days_arr = (['time','lat','lon',], coping_window_days_arr),  
                   recovery_window_days_arr = (['time','lat','lon',], recovery_window_days_arr),
               
                   severity_index = (['time','lat','lon',], severity_index)
                       ))
             



# In[ ]:




xrds.to_netcdf('/home2/datawork/slal/noaaoisst_mhw_output_test.nc')

