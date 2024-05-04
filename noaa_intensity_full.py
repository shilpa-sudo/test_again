#!/usr/bin/env python
# coding: utf-8

# In[1]:


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


# In[4]:


ds = xr.open_dataset('mhw_intensity_noaa.nc')
ds


# In[5]:


lat = np.array(ds.lat)
lon = np.array(ds.lon)
intensity3D = np.array(ds.mhw_intensity)


# In[6]:


xx,yy = np.meshgrid(lon,lat)
print(xx.shape,yy.shape)


# In[7]:


t = np.array(ds.time)

datessfull = [date.fromordinal(tt.astype(int)) for tt in t]
print(datessfull)


# In[8]:


intensity3D.shape


# In[ ]:


data_crs = ccrs.PlateCarree(central_longitude=0)

# But you want the map to include the dateline.
proj_crs = ccrs.PlateCarree(central_longitude=180)


eez = gpd.read_file("/home/shilpa/Downloads/World_EEZ_v11_20191118/eez_v11.shp")

for i in range(9861):
    
    fig = plt.figure(figsize=(11,4))

    ax = plt.axes(projection=proj_crs)

    clevs = np.arange(0.00001, 8, 0.5) 
    im=ax.contourf(xx[:,:], yy[:,:], intensity3D[:,:,i], 
                   levels=clevs,
                   cmap= "turbo",#"RdBu_r",#"RdBu_r",#"plasma",
                   transform=data_crs,
                   transform_first=True)  # This kwarg makes it much faster.
    cbar = plt.colorbar(im, ax=ax, orientation="horizontal", shrink=0.5)
    cbar.set_label("Intensity [°C]")
    #ax.coastlines()
    #eez[:].plot(ax=ax,color='none',edgecolor='black',transform=data_crs)

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
    gl.top_labels = False
    gl.right_labels = False
    #fig.suptitle('MHW events 1993-2019 (NOAAOISST)', fontsize=14)
    fig.suptitle(datessfull[i].strftime("%Y-%m-%d"), fontsize=14)
    ax.set_xlim(-35,29)
    ax.set_ylim(-34.875,-2.375)
    plt.tight_layout()
    plt.savefig(f'/home/shilpa/glory_mat_analysis/noaa_aoi_intensity_full_0_8_deg/mhw_intensity_aoi_{i}.png') 
    plt.close(fig)
    print(f'just finished {i}')


# In[ ]:





# In[ ]:




