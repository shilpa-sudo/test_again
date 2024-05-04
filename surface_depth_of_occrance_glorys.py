#!/usr/bin/env python
# coding: utf-8

# In[2]:


import netCDF4 as nc
import numpy as np
import xarray as xr
import pandas as pd
import os
import natsort


# In[3]:


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


datessfull = pd.date_range(start=('01/01/1993'), end=('12/31/2019'))


# In[ ]:


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


for depthind in range(len(rdlist)):
    
    test_array = np.zeros((len(datessfull),len(lats),len(lons)))

    for latind in range(len(lats[:])):
        for lonind in range(len(lons[:])):

            # for a particular location, read csv file
            lat = lats[latind]
            lon = lons[lonind]
            fsortedlist = fsorted2 (lat,lon)
            fn = (fsortedlist[depthind]) 
            
            if os.path.exists(fn):

                #print(fn)

                df0 = pd.read_csv(fn)

                active_dates_list = []
                # for each row (which represents each event) get the trange and append the trange to active dates list
                for i in range(len(df0)):
                    trange = pd.date_range(start=(df0['date_start']).iloc[i], end=(df0['date_end']).iloc[i])
                    #print(trange)
                    active_dates_list.append(trange)
                    #print(active_dates_list)

                active_date_range_depth = [element for sublist in active_dates_list for element in sublist]   
                #print(active_date_range_depth)

                for dt in range (len(datessfull)):
                    if datessfull[dt] in list(active_date_range_depth):
                        test_array[dt,latind,lonind] = rdlist[depthind]
    xrds = xr.Dataset(coords = dict(time = datessfull,latitude = lats,longitude = lons),data_vars = dict(verticalevents = (['time','latitude','longitude'],test_array)))
    
    xrds.to_netcdf(f'/home/shilpa/glory_mat_analysis/subsurface_glorys/vertical_events_glorys_quarter_{rdlist[depthind]}.nc')


# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:




