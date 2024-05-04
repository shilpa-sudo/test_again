#!/usr/bin/env python
# coding: utf-8

# In[2]:


# standard imports

import netCDF4 as nc
import numpy as np
import xarray as xr
import pandas as pd
import os
import natsort


# In[3]:


# get depth, lon and lat values
dx = xr.open_dataset('/home/shilpa/glory_mat_analysis/glorys_quaterdeg/2018_glorys_quarterdegree.nc')
dx

lati = dx.nav_lat.isel(x=10).values
longi = dx.nav_lon.isel(y=10).values
#print(lati)
#print(longi)
lats = np.array(lati)
lons = np.array(longi)

x = (dx.deptht.values)
xlist = list(x)
depth = xlist

rdlist = []
for i in range(len(depth)):
    
    rd = round(depth[i],2)
    
    rdlist.append(rd)
    
print(rdlist)


# In[4]:


# get full range of dates
datessfull = pd.date_range(start=('01/01/1993'), end=('12/31/2019'))


# In[5]:


# function to read all files for a particular lat and lon

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


# In[6]:


# for each day create a .nc file , 1 file per day with 3Dimensions, depth, lat and lon. 
# then do surface plots of each depth level
# so in one frame you can see all the depths where MHWs were detected


# In[ ]:


test_array = np.zeros((1,len(rdlist),len(lats),len(lons))) # make a test array to fill

for dt in range (len(datessfull)):
    print(f'processing {datessfull[dt]}')
    
    for latind in range(len(lats[:])):
        print(f'at latind = {latind}')
        
        for lonind in range(len(lons[:])):

            # for a particular location, read csv file
            lat = lats[latind]
            #print(f'at latind = {latind}')
            lon = lons[lonind]
            #print(f'at lonind = {lonind}')
            
            fsortedlist = fsorted2 (lat,lon) # get files for all depths for this lat, lon

            for depthind in range(len(rdlist)):
                #print(f'processing {rdlist[depthind]}')
                fn = (fsortedlist[depthind]) 

                if os.path.exists(fn):

                    df0 = pd.read_csv(fn) # for each depth read file

                    active_dates_list = []
                    # for each row (which represents each event) get the trange and append the trange to active dates list
                    for i in range(len(df0)):
                        trange = pd.date_range(start=(df0['date_start']).iloc[i], end=(df0['date_end']).iloc[i])
                        active_dates_list.append(trange)

                    active_date_range_depth = [element for sublist in active_dates_list for element in sublist]   
                    #make one list out of the several lists, at this particular depth , what are the active dates

                    if datessfull[dt] in list(active_date_range_depth): # check if each day falls in list
                        test_array[:,depthind,latind,lonind] = rdlist[depthind]

    xrds = xr.Dataset(coords = dict(time = datessfull[dt],depth = rdlist,latitude = lats,longitude = lons),data_vars = dict(verticalevents = (['time','depth','latitude','longitude'],test_array)))
    xrds.to_netcdf(f'/home/shilpa/glory_mat_analysis/subsurface_glorys/vertical_events_glorys_quarter_{dt:04d}.nc')


# In[11]:


#lats


# In[12]:


#lons


# In[13]:


#len(lons)


# In[14]:


#len(lats)


# In[ ]:




