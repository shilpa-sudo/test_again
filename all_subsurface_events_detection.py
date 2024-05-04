#!/usr/bin/env python
# coding: utf-8

# In[1]:


# standard imports

import netCDF4 as nc
import numpy as np
import xarray as xr
import pandas as pd
import os
import natsort


# In[ ]:


# get depth, lon and lat values
dx = xr.open_dataset('/home/datawork-lead/datarmor-only/shilpa/glorys_quarter_datarmor/2018_glorys_quarterdegree_computed_on_datarmor.nc')


lati = dx.nav_lat.isel(x=10).values
longi = dx.nav_lon.isel(y=10).values

lats = np.array(lati)
lons = np.array(longi)

x = (dx.deptht.values)
xlist = list(x)
depth = xlist

rdlist = []
for i in range(len(depth)):
    
    rd = round(depth[i],2)
    
    rdlist.append(rd)
    


# In[ ]:


# get full range of dates
datessfull = pd.date_range(start=('01/01/1993'), end=('12/31/2019'))


# In[ ]:


dep = [0,1,2,3,5,6,7,10,11,13,15,18,21,25,29,34,40,47,55,65,77,92,100,130,155,186,222,266,318,380,453,500,600,700,900,1000,1200,1500]


# In[ ]:


# list of 38 files
fsortedlist = []
for depp in dep:
    fname = f'/home2/datawork/slal/mhws_fullset_depth{depp}m_withdepth.csv'
    fsortedlist.append(fname)


# In[ ]:


test_array = np.zeros((1,len(rdlist),len(lats),len(lons))) # make a test array to fill

for dt in range (len(datessfull[400])):
    print(f'processing {datessfull[dt]}')
    
    for latind in range(len(lats[:])):
        print(f'at latind = {latind}')
        
        for lonind in range(len(lons[:])):

            # for a particular location, read csv file
            lat = lats[latind]
            #print(f'at latind = {latind}')
            lon = lons[lonind]
            #print(f'at lonind = {lonind}')
            
            #fsortedlist = fsorted2 (lat,lon) # get files for all depths for this lat, lon

            for depthind in range(len(rdlist)):
                #print(f'processing {rdlist[depthind]}')
                fn = (fsortedlist[depthind]) 

                if os.path.exists(fn):

                    df0 = pd.read_csv(fn) 
                    
                    df1 = df0[(df0['latitude'] == lat) & (df0['longitude'] == lon)]

                    active_dates_list = []
                    # for each row (which represents each event) get the trange and append the trange to active dates list
                    for i in range(len(df1)):
                        trange = pd.date_range(start=(df1['date_start']).iloc[i], end=(df1['date_end']).iloc[i])
                        active_dates_list.append(trange)

                    active_date_range_depth = [element for sublist in active_dates_list for element in sublist]   
                    #make one list out of the several lists, at this particular depth , what are the active dates

                    if datessfull[dt] in list(active_date_range_depth): # check if each day falls in list
                        test_array[:,depthind,latind,lonind] = rdlist[depthind]

    xrds = xr.Dataset(coords = dict(time = datessfull[dt],depth = rdlist,latitude = lats,longitude = lons),data_vars = dict(verticalevents = (['time','depth','latitude','longitude'],test_array)))
    xrds.to_netcdf(f'/home2/datawork/slal/all_subsurface_events_glorys/vertical_events_glorys_quarter_{dt:04d}.nc')
 


# In[ ]:





# In[ ]:





# In[ ]:




