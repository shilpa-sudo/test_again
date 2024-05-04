#!/usr/bin/env python
# coding: utf-8

# In[1]:


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
from matplotlib.patches import Rectangle
import matplotlib.path as mpath
from shapely.geometry import Polygon
from sklearn.linear_model import LinearRegression
import matplotlib.colors as mcolors



from matplotlib.patches import Polygon as mplPolygon
from shapely.geometry import Polygon as ShapelyPolygon
from shapely.geometry import Point as ShapelyPoint
from shapely.ops import cascaded_union
from shapely.geometry import MultiPolygon

from matplotlib.animation import FuncAnimation
import os


# In[2]:


# get spatial extent noaa
sp_ext_noaa = xr.open_dataset('/media/shilpa/Expansion/spatial_extent_noaa_aoi_masking.nc')


# In[3]:


sp_ext_noaa


# In[4]:


# get spatial extent glorys
sp_ext_glorys = xr.open_dataset('/media/shilpa/Expansion/spatial_extent_glorys_aoi_masking.nc')


# In[5]:


sp_ext_glorys


# In[6]:


data_crs = ccrs.PlateCarree(central_longitude=0)

# But you want the map to include the dateline.
proj_crs = ccrs.PlateCarree(central_longitude=180)


# In[7]:


eez = gpd.read_file("/home/shilpa/Downloads/World_EEZ_v11_20191118/eez_v11.shp")


# In[8]:


fiji = eez[eez["TERRITORY1"] == "Fiji"] 
ncd = eez[eez["TERRITORY1"] == "New Caledonia"]
vanuatu = eez[eez["TERRITORY1"] == "Vanuatu"]
wf = eez[eez["TERRITORY1"] == "Wallis and Futuna"]
solo = eez[eez["TERRITORY1"] == "Solomon Islands"]
tonga = eez[eez["TERRITORY1"] == "Tonga"]
samoa = eez[eez["TERRITORY1"] == "Samoa"]
tuvalu = eez[eez["TERRITORY1"] == "Tuvalu"]
niue = eez[eez["TERRITORY1"] == "Niue"]
cooks = eez[eez["TERRITORY1"] == "Cook Islands"]
tokelau = eez[eez["TERRITORY1"] == "Tokelau"]
amsam  = eez[eez["TERRITORY1"] == "American Samoa"]


# In[9]:


#noaa_sp_ext = np.array(sp_ext_noaa.spatial_extent_pacific)
glorys_sp_ext = np.array(sp_ext_glorys.spatial_extent_pacific)


# In[10]:


start_time = '1993-01-01'  # Replace with your actual start time
end_time = '2019-12-31' 
noaa_1993_2019 =  sp_ext_noaa.sel(time=slice(start_time, end_time))


# In[11]:


noaa_sp_ext = np.array(noaa_1993_2019.spatial_extent_pacific)


# In[12]:


longi = sp_ext_glorys.lon
lati = sp_ext_glorys.lat
xx,yy = np.meshgrid(longi,lati)


# In[13]:


datessfull = pd.date_range(start='1993-01-01', end='2019-12-31', freq='D')


# In[ ]:


def update(frame, ax):
    ax.clear()
    ax.coastlines()
    clevs = np.arange(0.5, 1.5, 0.5) 


    im = ax.contourf(xx[:,:], yy[:,:], noaa_sp_ext[frame,:,:], 
                   levels=clevs,
                   cmap= "Blues",#"RdBu_r",#"plasma",
                   alpha = 0.5,
                   transform=data_crs,
                   transform_first=True) 

    im = ax.contourf(xx[:,:], yy[:,:], glorys_sp_ext[frame,:,:], 
                   levels=clevs,
                   cmap= "Reds",#"RdBu_r",#"plasma",
                   alpha = 0.5,
                   transform=data_crs,
                   transform_first=True)  


    #cbar = plt.colorbar(im, ax=ax, orientation="horizontal", shrink=0.5)
    #cbar.set_label("No.of MHW events")
    #ax.coastlines()
    #eez[:].plot(ax=ax,color='none',edgecolor='black',transform=data_crs)



    fiji.boundary.plot(ax=ax, color="black",alpha = 0.4,transform=data_crs)
    ncd.boundary.plot(ax=ax, color="black",alpha = 0.4,transform=data_crs)
    vanuatu.boundary.plot(ax=ax, color="black",alpha = 0.4,transform=data_crs)
    wf.boundary.plot(ax=ax, color="black",alpha = 0.4,transform=data_crs)
    solo.boundary.plot(ax=ax, color="black",alpha = 0.4,transform=data_crs)
    tonga.boundary.plot(ax=ax, color="black",alpha = 0.4,transform=data_crs)
    samoa.boundary.plot(ax=ax, color="black",alpha = 0.4,transform=data_crs)
    tuvalu.boundary.plot(ax=ax, color="black",alpha = 0.4,transform=data_crs)
    niue.boundary.plot(ax=ax, color="black",alpha = 0.4,transform=data_crs)
    cooks.boundary.plot(ax=ax, color="black",alpha = 0.4,transform=data_crs)
    tokelau.boundary.plot(ax=ax, color="black",alpha = 0.4,transform=data_crs)
    amsam.boundary.plot(ax=ax, color="black",alpha = 0.4,transform=data_crs)

    

    gl = ax.gridlines(draw_labels=True)
    gl.top_labels = False
    gl.right_labels = False

    #fig.suptitle('MHW events 1993-2019 (NOAAOISST)', fontsize=14)
    fig.suptitle(datessfull[frame].strftime("%Y-%m-%d"),fontsize=14)
    ax.set_title('MHW Spatial extent NOAAOISST (Blue), GLORYS (Red)')
    ax.set_xlim(-35,29)
    ax.set_ylim(-34.875,-2.375)
    #plt.tight_layout()
    #plt.savefig('mhwstate_aoi_231197.png')

# Create a figure
fig, ax = plt.subplots(subplot_kw={'projection': proj_crs}, figsize=(13, 5))
ax.coastlines()
clevs = np.arange(0.5, 1.5, 0.5) 


im = ax.contourf(xx[:,:], yy[:,:], noaa_sp_ext[0,:,:], 
               levels=clevs,
               cmap= "Blues",#"RdBu_r",#"plasma",
               alpha = 0.5,
               transform=data_crs,
               transform_first=True) 

im = ax.contourf(xx[:,:], yy[:,:], glorys_sp_ext[0,:,:], 
               levels=clevs,
               cmap= "Reds",#"RdBu_r",#"plasma",
               alpha = 0.5,
               transform=data_crs,
               transform_first=True)  


#cbar = plt.colorbar(im, ax=ax, orientation="horizontal", shrink=0.5)
#cbar.set_label("No.of MHW events")
#ax.coastlines()
#eez[:].plot(ax=ax,color='none',edgecolor='black',transform=data_crs)


fiji.boundary.plot(ax=ax, color="black",alpha = 0.4,transform=data_crs)
ncd.boundary.plot(ax=ax, color="black",alpha = 0.4,transform=data_crs)
vanuatu.boundary.plot(ax=ax, color="black",alpha = 0.4,transform=data_crs)
wf.boundary.plot(ax=ax, color="black",alpha = 0.4,transform=data_crs)
solo.boundary.plot(ax=ax, color="black",alpha = 0.4,transform=data_crs)
tonga.boundary.plot(ax=ax, color="black",alpha = 0.4,transform=data_crs)
samoa.boundary.plot(ax=ax, color="black",alpha = 0.4,transform=data_crs)
tuvalu.boundary.plot(ax=ax, color="black",alpha = 0.4,transform=data_crs)
niue.boundary.plot(ax=ax, color="black",alpha = 0.4,transform=data_crs)
cooks.boundary.plot(ax=ax, color="black",alpha = 0.4,transform=data_crs)
tokelau.boundary.plot(ax=ax, color="black",alpha = 0.4,transform=data_crs)
amsam.boundary.plot(ax=ax, color="black",alpha = 0.4,transform=data_crs)


gl = ax.gridlines(draw_labels=True)
gl.top_labels = False
gl.right_labels = False

#fig.suptitle('MHW events 1993-2019 (NOAAOISST)', fontsize=14)
fig.suptitle(datessfull[0].strftime("%Y-%m-%d"),fontsize=14)
ax.set_title('MHW Spatial extent NOAAOISST (Blue), GLORYS (Red)')
ax.set_xlim(-35,29)
ax.set_ylim(-34.875,-2.375)
#plt.tight_layout()

# Create a FuncAnimation object
animation = FuncAnimation(fig, update, frames=9861, interval=100, fargs=(ax,))  # 100 frames, 100ms per frame

# Create a directory to save the frames as PNG images
output_dir = 'new_testframe'
#os.makedirs(output_dir, exist_ok=True)

# Save the frames as individual PNG images
for frame in range(3000):
    update(frame, ax)
    plt.savefig(f'{output_dir}/frame_{frame:04d}.png', dpi=100)

# Use FFmpeg to create an MP4 video from the PNG frames
#os.system('ffmpeg -framerate 10 -i new_testframe/frame_%04d.png -c:v libx264 -r 30 -pix_fmt yuv420p mhw_glorys_noaa_sp_ext_93_19.mp4')

# Clean up the temporary frame files
#for frame in range(9861):
#    os.remove(f'{output_dir}/frame_{frame:03d}.png')

# Optionally, remove the frames directory
#os.rmdir(output_dir)

# Show the video in Matplotlib (optional)
# You can also save the animation as a GIF or other formats if needed.
#plt.show()


# In[ ]:




