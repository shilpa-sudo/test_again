#!/usr/bin/env python
# coding: utf-8

# In[1]:


# check if there is a rising trend in sst in the AOI


# In[2]:


get_ipython().run_line_magic('matplotlib', 'notebook')
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
import math


# In[3]:


data_crs = ccrs.PlateCarree(central_longitude=0)

# But you want the map to include the dateline.
proj_crs = ccrs.PlateCarree(central_longitude=180)


# In[4]:


eez = gpd.read_file("/home/shilpa/Downloads/World_EEZ_v11_20191118/eez_v11.shp") # /home/shilpa/Downloads/World_EEZ_v11_20191118/eez_v11.shp, /home/shilpa/Downloads/World_EEZ_v11_20191118/eez_boundaries_v11.shp
#print(eez)


# In[5]:


ds = xr.open_dataset('mhw_intensity_noaa.nc')
ds


# In[ ]:


dtf_p90 = Dataset('mhw_state_AOI_quarterdegree.nc') 
mhwstate_p90 = dtf_p90.variables['mhwstate']
lat = dtf_p90.variables['lat']
lon = dtf_p90.variables['lon']
time = dtf_p90.variables['time']
tt = np.array(time)
t = tt.astype(int)
datessfull = [date.fromordinal(tt.astype(int)) for tt in t]

