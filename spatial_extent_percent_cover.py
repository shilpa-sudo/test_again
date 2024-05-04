#!/usr/bin/env python
# coding: utf-8

# In[1]:


# standard imports

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


# In[2]:


#get dataset
dtf = Dataset('noaaoisst93_19_sst_AOI_quarterdegree.nc') 
lati= dtf.variables['lati'] #degrees north
longi = dtf.variables['longi'] #degrees east
time = dtf.variables['time']
temp = dtf.variables['temp']


# In[3]:


tt = np.array(time)
t = tt.astype(int)
datessfull = [date.fromordinal(ttt.astype(int)) for ttt in t]
#print(datessfull)


# In[4]:


dtf_p50 = Dataset('mhw_state_AOI_quarterdegree_50p.nc') 
mhwstate_p50 = dtf_p50.variables['mhwstate']

dtf_p85 = Dataset('mhw_state_AOI_quarterdegree_85p.nc') 
mhwstate_p85 = dtf_p85.variables['mhwstate']

dtf_p90 = Dataset('mhw_state_AOI_quarterdegree.nc') 
mhwstate_p90 = dtf_p90.variables['mhwstate']

dtf_p99 = Dataset('mhw_state_AOI_quarterdegree_99p.nc') 
mhwstate_p99 = dtf_p99.variables['mhwstate']


# In[ ]:


percent_mhwstatelist_p50 = []
percent_mhwstatelist_p85 = []
percent_mhwstatelist_p99 = []


for i in range(9861):
    percent_mhwstate_p50 = ((np.count_nonzero(mhwstate_p50[:,:,i]==1)/31810))*100   #31810 = number of cells with data
    pmhws_p50 = round(percent_mhwstate_p50,2)
    percent_mhwstatelist.append(pmhws_p50)
    
    percent_mhwstate_p85 = ((np.count_nonzero(mhwstate_p85[:,:,i]==1)/31810))*100   #31810 = number of cells with data
    pmhws_p85 = round(percent_mhwstate_p85,2)
    percent_mhwstatelist.append(pmhws_p85)
    
    percent_mhwstate_p99 = ((np.count_nonzero(mhwstate_p99[:,:,i]==1)/31810))*100   #31810 = number of cells with data
    pmhws_p99 = round(percent_mhwstate_p99,2)
    percent_mhwstatelist.append(pmhws_p99)


# In[ ]:


pmhsarr_p50 = np.array(percent_mhwstatelist_p50)
dataframe_p50 = pd.DataFrame(pmhsarr_p50) 
dataframe_p50.to_csv("percent_mhw_state_noaa_quarterdegree_p50.csv")

pmhsarr_p85 = np.array(percent_mhwstatelist_p85)
dataframe_p85 = pd.DataFrame(pmhsarr_p85) 
dataframe_p85.to_csv("percent_mhw_state_noaa_quarterdegree_p85.csv")

pmhsarr_p99 = np.array(percent_mhwstatelist_p99)
dataframe_p99 = pd.DataFrame(pmhsarr_p99) 
dataframe_p99.to_csv("percent_mhw_state_noaa_quarterdegree_p99.csv")


# In[ ]:


df_p90 = pd.read_csv('percent_mhw_state_noaa_quarterdegree.csv')
mhws_percentage_p90 = df_p90.to_numpy()
#print(mhws_percentage_p90.shape)
#print(mhws_percentage_p90[0,1])


# In[ ]:


fig,ax = plt.subplots(figsize=(12,5))
ax.plot(datessfull,mhws_percentage_p99[:,1],'k')
ax.plot(datessfull,mhws_percentage_p90[:,1],'b')
ax.plot(datessfull,mhws_percentage_p85[:,1],'r')
ax.plot(datessfull,mhws_percentage_p50[:,1],'g')
ax.set_ylabel('Percentage area coverage MHW ')


# In[ ]:


import cartopy
import cartopy.mpl.geoaxes
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

fig,ax = plt.subplots(figsize=(12,8))

ax.plot(datessfull,mhws_percentage_p99[:,1],'k')
ax.plot(datessfull,mhws_percentage_p90[:,1],'b')
ax.plot(datessfull,mhws_percentage_p85[:,1],'r')
ax.plot(datessfull,mhws_percentage_p50[:,1],'g')


ax.set_ylabel('Percentage area coverage MHW ')

axins = inset_axes(ax, width="30%", height="30%", loc="upper left", 
                   axes_class=cartopy.mpl.geoaxes.GeoAxes, 
                   axes_kwargs=dict(map_projection=cartopy.crs.PlateCarree(central_longitude=180)))
axins.add_feature(cartopy.feature.COASTLINE)
axins.set_xlim(-35,29)
axins.set_ylim(-34.875,-2.375)

plt.show()

