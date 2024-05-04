#!/usr/bin/env python
# coding: utf-8

# In[1]:


#in this notebook we will regrid glorys using noaa quarter degree


# In[2]:


import sys
import os
import numpy as np
import pandas as pd
import xarray as xr
import netCDF4
import matplotlib.pyplot as plt


# In[3]:


path = "/media/shilpa/LaCie/GLORYS12V1/"
#files = path+ "GLORYS12V1_1dAV_[1-2][0-9][0-9][0-9]*gridT*.nc"
files = path+ "GLORYS12V1_1dAV_2019*gridT*.nc"


# In[4]:


ens=xr.open_mfdataset(paths=files,combine='nested',concat_dim="time_counter",data_vars='all',
                       decode_times=True)


# In[5]:


ens


# In[6]:


ens["nav_lon"].values


# In[7]:


ens["nav_lonc"] = ens["nav_lon"].where(ens["nav_lon"] > 0,ens["nav_lon"]+360) 
#select for values greater than 0, if less then 0 add 360


# In[8]:


ens["nav_lonc"].values


# In[9]:


longitude = ens.nav_lonc.isel(y=10).values


# In[10]:


longitude


# In[11]:


longitude = ens.nav_lonc.isel(y=10).values
latitude = ens.nav_lat.isel(x=10).values
ens = ens.assign_coords({"x": longitude, "y": latitude})
sst = ens.votemper.isel(deptht=0)

noaa = xr.open_dataset("/media/shilpa/LaCie/noaaoisst/noaaoisst_seasonal_mean.nc")
noaa_zoom= noaa.sel(lat = slice(-35, 0.25), lon = slice(140,291), drop = True)

lons = noaa_zoom.lon.values
lats = noaa_zoom.lat.values


# In[12]:


lons


# In[13]:



dds2 = ens.interp(x=lons, y=lats, method="linear")

#sst = dds2.votemper.isel(deptht=0)


# In[14]:


dds2


# In[15]:


#dds2.to_netcdf('1993_glorys_quarterdegree.nc', 'w')


# In[16]:


from datetime import datetime
from pathlib import Path

from dask.diagnostics import ProgressBar


# In[ ]:


data = dds2
out_file = '2019_glorys_quarterdegree.nc'
write_job = data.to_netcdf(out_file, compute=False)
with ProgressBar():
    print(f"Writing to {out_file}")
    write_job.compute()


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





# In[ ]:




