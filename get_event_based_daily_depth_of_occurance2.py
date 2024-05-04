#!/usr/bin/env python
# coding: utf-8

# In[ ]:


# standard imports

import netCDF4 as nc
import numpy as np
import xarray as xr
import pandas as pd
import os
import natsort


# In[ ]:





# In[ ]:


# get full range of dates
datessfull = pd.date_range(start=('01/01/1993'), end=('12/31/2019'))


# In[ ]:


# all events in same file

fname = f'/home2/datawork/slal/all_subsurface_mhws.csv'
df1 = pd.read_csv(fname)  


# In[ ]:


# for each row (which represents each event) get the trange and append the trange to active dates list
for i in range(len(df1)):
    active_depth_list = []
    active_date_list = []
    lat_list = []
    lon_list = []
    
    trange = pd.date_range(start=(df1['date_start']).iloc[i], end=(df1['date_end']).iloc[i])
    active_date_list.append(trange)
    depth = [(df1['depth']).iloc[i]] * (len(trange))
    active_depth_list.append(depth)
    lat = [(df1['latitude']).iloc[i]] * (len(trange))
    lat_list.append(lat)    
    lon = [(df1['longitude']).iloc[i]] * (len(trange))
    lon_list.append(lon)                

    each_day_each_event = pd.DataFrame({'time': active_date_list,'depth': active_depth_list, 'latitude':lat_list, 'longitude':lon_list})
    each_day_each_event.to_csv(f'/home2/datawork/slal/all_subsurface_events_glorys/{i}.csv', index=False)
   

