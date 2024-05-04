#!/usr/bin/env python
# coding: utf-8

# In[1]:


## saving 0.25 degree MHW sensitive areas


# In[2]:


# standard imports

#%matplotlib notebook
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


# In[3]:


import scipy as sp
from scipy import linalg
from scipy import stats
import scipy.ndimage as ndimage
from datetime import date


# In[6]:


#!ncdump -h mhw_state_AOI_quarterdegree_50p.nc


# In[7]:


dtf_p50 = Dataset('mhw_state_AOI_quarterdegree_50p.nc') 
mhwstate_p50 = dtf_p50.variables['mhwstate']
lati= dtf_p50.variables['lat'] #degrees north
longi = dtf_p50.variables['lon'] #degrees east
time = dtf_p50.variables['time']
tt = np.array(time)
t = tt.astype(int)
datessfull = [date.fromordinal(ttt.astype(int)) for ttt in t]

dtf_p75 = Dataset('mhw_state_AOI_quarterdegree_75p.nc') 
mhwstate_p75 = dtf_p75.variables['mhwstate']

dtf_p85 = Dataset('mhw_state_AOI_quarterdegree_85p.nc') 
mhwstate_p85 = dtf_p85.variables['mhwstate']

dtf_p90 = Dataset('mhw_state_AOI_quarterdegree.nc') 
mhwstate_p90 = dtf_p90.variables['mhwstate']

dtf_p99 = Dataset('mhw_state_AOI_quarterdegree_99p.nc') 
mhwstate_p99 = dtf_p99.variables['mhwstate']


# In[8]:


df_p90 = pd.read_csv('percent_mhw_state_noaa_quarterdegree.csv')
mhws_percentage_p90 = df_p90.to_numpy()
#print(mhws_percentage_p90.shape)
#print(mhws_percentage_p90[0,1])
df_p99 = pd.read_csv('percent_mhw_state_noaa_quarterdegree_p99.csv')
mhws_percentage_p99 = df_p99.to_numpy()

df_p85 = pd.read_csv('percent_mhw_state_noaa_quarterdegree_p85.csv')
mhws_percentage_p85 = df_p85.to_numpy()

df_p75 = pd.read_csv('percent_mhw_state_noaa_quarterdegree_p75.csv')
mhws_percentage_p75 = df_p75.to_numpy()

df_p50 = pd.read_csv('percent_mhw_state_noaa_quarterdegree_p50.csv')
mhws_percentage_p50 = df_p50.to_numpy()


# In[9]:


data_crs = ccrs.PlateCarree(central_longitude=0)

# But you want the map to include the dateline.
proj_crs = ccrs.PlateCarree(central_longitude=180)


# In[10]:


eez = gpd.read_file("/home/shilpa/Downloads/World_EEZ_v11_20191118/eez_v11.shp")


# In[12]:


xx,yy = np.meshgrid(longi,lati)


# In[26]:


#print('%s pm is the time right now. In one hour, it will be %s pm.'%(3,4))


# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:


#do this slow and steady....it might be painfully slow now but it will be worth in the end when u see the gif

for i in range(353,700,1):
    
    fig = plt.figure(figsize=(13,5))

    ax = plt.axes(projection=proj_crs)

    clevs = np.arange(0.5, 1.5, 0.5) 
    
    im = ax.contourf(xx[:,:], yy[:,:], mhwstate_p50[:,:,i], 
                   levels=clevs,
                   cmap= "Greens",#"RdBu_r",#"plasma",
                     alpha=0.5,
                   transform=data_crs,
                   transform_first=True)  # This kwarg makes it much faster.

    im = ax.contourf(xx[:,:], yy[:,:], mhwstate_p75[:,:,i], 
                   levels=clevs,
                   cmap= "Oranges",#"RdBu_r",#"plasma",
                     alpha=0.5,
                   transform=data_crs,
                   transform_first=True)  # This kwarg makes it much faster.


    im = ax.contourf(xx[:,:], yy[:,:], mhwstate_p85[:,:,i], 
                   levels=clevs,
                   cmap= "PuRd",#"RdBu_r",#"plasma",
                     alpha=0.9,
                   transform=data_crs,
                   transform_first=True)  # This kwarg makes it much faster.


    im = ax.contourf(xx[:,:], yy[:,:], mhwstate_p90[:,:,i], 
                   levels=clevs,
                   cmap= "plasma",#"RdBu_r",#"plasma",
                   transform=data_crs,
                   transform_first=True)  
    
    im = ax.contourf(xx[:,:], yy[:,:], mhwstate_p99[:,:,i], 
               levels=clevs,
               cmap= "Reds",#"RdBu_r",#"plasma",
               transform=data_crs,
               transform_first=True)  
    # This kwarg makes it much faster.
#cbar = plt.colorbar(im, ax=ax, orientation="horizontal", shrink=0.5)
#cbar.set_label("No.of MHW events")
    #ax.coastlines()
    
    eez[0:1].plot(ax=ax,color='none',edgecolor='black',transform=data_crs)
    eez[2:3].plot(ax=ax,color='none',edgecolor='black',transform=data_crs)
    eez[4:5].plot(ax=ax,color='none',edgecolor='black',transform=data_crs)
    eez[7:9].plot(ax=ax,color='none',edgecolor='black',transform=data_crs)
    eez[11:12].plot(ax=ax,color='none',edgecolor='black',transform=data_crs)
    eez[15:16].plot(ax=ax,color='none',edgecolor='black',transform=data_crs)
    eez[20:21].plot(ax=ax,color='none',edgecolor='black',transform=data_crs)
    eez[33:35].plot(ax=ax,color='none',edgecolor='black',transform=data_crs)
    eez[37:38].plot(ax=ax,color='none',edgecolor='black',transform=data_crs)
    eez[47:48].plot(ax=ax,color='none',edgecolor='black',transform=data_crs)
    eez[51:52].plot(ax=ax,color='none',edgecolor='black',transform=data_crs)
    eez[53:54].plot(ax=ax,color='none',edgecolor='black',transform=data_crs)
    eez[248:250].plot(ax=ax,color='none',edgecolor='black',transform=data_crs)
    eez[243:245].plot(ax=ax,color='none',edgecolor='black',transform=data_crs)
    eez[246:248].plot(ax=ax,color='none',edgecolor='black',transform=data_crs)
    eez[251:253].plot(ax=ax,color='none',edgecolor='black',transform=data_crs)
    eez[259:260].plot(ax=ax,color='none',edgecolor='black',transform=data_crs)
    eez[261:262].plot(ax=ax,color='none',edgecolor='black',transform=data_crs)
    
    gl = ax.gridlines(draw_labels=True)
    gl.xlabels_top = False
    gl.ylabels_right = False
    fig.suptitle(datessfull[i].strftime("%Y-%m-%d"), fontsize=14)
    #ax.set_title('Percent area in MHW state:%s'%percent_mhwstatelist[i])
    #ax.set_title('Area in MHW state: %s Percent'% mhws_percentage[i,1])
    ax.set_title('Percent Area in MHW state: %s (P99), %s (P90), %s (P85), %s (P75), %s (P50)'%(mhws_percentage_p99[i,1], mhws_percentage_p90[i,1] ,mhws_percentage_p85[i,1] ,mhws_percentage_p75[i,1] ,mhws_percentage_p50[i,1]))
    ax.set_xlim(-35,29)
    ax.set_ylim(-34.875,-2.375)
    plt.tight_layout()

    plt.savefig('noaa_0.25degree_mhw_sensitive_areas_mhw_areas/mhw_state_%s.png'%i)

