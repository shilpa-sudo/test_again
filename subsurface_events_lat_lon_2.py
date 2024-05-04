#!/usr/bin/env python
# coding: utf-8

# In[2]:


import netCDF4 as nc
import numpy as np
import xarray as xr
import pandas as pd
import os


# In[3]:


dtime = pd.date_range('01/01/1993','12/31/2019')


# In[4]:


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


# In[5]:


dx = xr.open_dataset('/home/shilpa/glory_mat_analysis/glorys_quaterdeg/2018_glorys_quarterdegree.nc')
dx

lati = dx.nav_lat.isel(x=10).values
longi = dx.nav_lon.isel(y=10).values
#print(len(lati))
#print(len(longi))


x = (dx.deptht.values)
xlist = list(x)
depth = xlist

rdlist = []
for i in range(len(depth)):

    rd = round(depth[i],2)

    rdlist.append(rd)

print(rdlist)


# In[ ]:


for lat in lati[70:]:
    for lon in longi[300:]:
        
        fsortedlist = fsorted2 (lat,lon)
        time = []
        depths = []
        for fnx in range(len(rdlist)):
                if os.path.exists(fsortedlist[fnx]):
                    df1 = pd.read_csv(fsortedlist[fnx])
                    
                    for i in range(len(df1)):
    
                        trange = pd.date_range(start=(df1['date_start']).iloc[i], end=(df1['date_end']).iloc[i])
                        time.extend(trange)
                        depth = [rdlist[fnx]] * (len(trange))
                        depths.extend(depth)
        dtf = pd.DataFrame({'time': time,'depths': depths})
        dtf.to_csv(f'/home/shilpa/glory_mat_analysis/south_pacific_subsurface/subsurface_{lat}_{lon}.csv', index=False)


# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:




