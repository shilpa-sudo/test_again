#!/usr/bin/env python
# coding: utf-8

# In[1]:


import xarray as xr
import numpy as np
import pandas as pd


# In[ ]:


ds = xr.open_dataset('/home/datawork-lead/datarmor-only/shilpa/spatial_extent_noaa_aoi_masking.nc')


# In[ ]:


sp_mask = np.array(ds.spatial_extent_pacific)


# In[ ]:


mask = sp_mask == 1.0


# In[ ]:


ds_intensity = xr.open_dataset('/home/datawork-lead/datarmor-only/shilpa/intensity_t_minus_clim_noaa_aoi_masking.nc')


# In[ ]:


variable = np.array(ds_intensity.intensity_t_minus_clim_pacific)


# In[ ]:


variable_mask = variable * mask


# In[ ]:


lati = ds.lat.values
longi = ds.lon.values


# In[ ]:


date_range = pd.date_range(start='1981-09-01', end='2023-06-26')


# In[ ]:


xrds = xr.Dataset(
       coords = dict(time = date_range,lat = lati,lon = longi),
       data_vars = dict(
                   mhw_intensity_t_minus_clim = (['time','lat','lon'],variable_mask)))

xrds.to_netcdf('/home/datawork-lead/datarmor-only/shilpa/aoi_mhw_intensity_noaaoisst.nc')


# In[ ]:





# In[ ]:




