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
from matplotlib.patches import Rectangle
import matplotlib.path as mpath
from shapely.geometry import Polygon
from sklearn.linear_model import LinearRegression
import matplotlib.colors as mcolors
from natsort import natsorted


from matplotlib.patches import Polygon as mplPolygon
from shapely.geometry import Polygon as ShapelyPolygon
from shapely.geometry import Point as ShapelyPoint
from shapely.ops import cascaded_union
from shapely.geometry import MultiPolygon

from matplotlib.animation import FuncAnimation
import os

from matplotlib.colors import Normalize


# In[2]:


date_range = pd.date_range(start='1993-01-01', end='2019-12-31', freq='D')


# In[3]:


eez = gpd.read_file("/home/shilpa/Downloads/World_EEZ_v11_20191118/eez_v11.shp")


# In[4]:


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


# In[5]:


data_crs = ccrs.PlateCarree(central_longitude=0)

# But you want the map to include the dateline.
proj_crs = ccrs.PlateCarree(central_longitude=180)


# In[6]:


ds = xr.open_mfdataset('aoi_vertical_extent/*.nc')


# In[7]:


ds


# In[8]:


longi = ds.lon
longi = xr.where(longi >  0, longi, (longi+360))
lati = ds.lat
xx,yy = np.meshgrid(longi,lati)


# In[9]:


#lngi = np.array(longi)
#print(lngi)


# In[10]:


#longilabels = xr.where(lngi < 180, lngi, (lngi-360))


# In[11]:


#longilabels


# In[10]:


#uneven_levels = [0.0000001,10,20,30,40,50,60,70,80,90,100,200,300,400,500,550,600,650,700,750,800,850,900,1000,1500]


# In[11]:


glorys_v_extent = ds.v_extent


# In[12]:


ve_array = np.array(glorys_v_extent)


# In[13]:


ve_array[ve_array == 0] = np.nan


# In[48]:


norm = Normalize(vmin=0, vmax=1500) 


# 

# In[42]:


def easy_plot_single_test(ax,frame, date_range):
    
    surf = ax.plot_surface(xx, yy, ve_array[frame,:,:],cmap = 'turbo',norm=norm, alpha = 0.8)
    ax.set_zlim(0, 1600)

    ax.invert_zaxis()
    # Add labels and a color bar
    ax.set_xlabel('longitude')
    ax.set_ylabel('latitude')
    ax.set_zlabel('depth')
    #cbar = fig.colorbar(surf, shrink=0.5, orientation='horizontal', label='Depth (meters)')
    datet = date_range[frame].strftime('%Y-%m-%d')
    ax.set_title(f'GLORYS vertical extent {datet}', fontsize=14)


# In[50]:


def update(frame, ax):
    ax.clear()
    
    easy_plot_single_test(ax,frame=frame, date_range=date_range)
  
# Create a figure

fig = plt.figure(figsize=(19, 19))
ax = fig.add_subplot(111, projection='3d',)
surf = ax.plot_surface(xx, yy, ve_array[0,:,:],cmap = 'turbo',norm=norm, alpha = 0.8)
ax.set_zlim(0, 1600)
ax.set_xlim(-35,29)
ax.set_ylim(-30.875,-2.375)
ax.invert_zaxis()
# Add labels and a color bar
ax.set_xlabel('longitude')
ax.set_ylabel('latitude')
ax.set_zlabel('depth')
cbar = fig.colorbar(surf, shrink=0.5, orientation='horizontal', label='Depth (meters)')



# Create a FuncAnimation object
animation = FuncAnimation(fig, update, frames=9861, interval=100, fargs=(ax,))  # 100 frames, 100ms per frame

# Create a directory to save the frames as PNG images
output_dir = 'glorys_vertical_extent_full'
#os.makedirs(output_dir, exist_ok=True)

# Save the frames as individual PNG images
for frame in range(9861):
    update(frame, ax)
    plt.savefig(f'{output_dir}/frame_{frame:04d}.png', dpi=100)

# Use FFmpeg to create an MP4 video from the PNG frames
os.system('ffmpeg -framerate 10 -i glorys_vertical_extent_full/frame_%04d.png -c:v libx264 -r 30 -pix_fmt yuv420p 3D_mhw_glorys_maxarea_with_vertical_extent_full.mp4')


# In[ ]:





# In[ ]:




