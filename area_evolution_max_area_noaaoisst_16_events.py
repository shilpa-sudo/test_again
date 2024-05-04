#!/usr/bin/env python
# coding: utf-8

# In[ ]:


# make time evolution of 16 max area events in noaaoisst max area event


# In[1]:


#get_ipython().run_line_magic('matplotlib', 'notebook')
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


#get noaaoisst


# In[3]:


noaa_polygon_ds_full = xr.open_dataset('noaaoisst_full_mhw_areas_attempt.nc')
start_time = '1993-01-01'  # Replace with your actual start time
end_time = '2019-12-31' 
noaa_polygon_ds =  noaa_polygon_ds_full.sel(time=slice(start_time, end_time))


# In[4]:


noaaarea = np.array(noaa_polygon_ds.mhw_spatial_extent_area)


# In[5]:


noaaarea.shape


# In[6]:


longi_n = noaa_polygon_ds.lon
lati_n = noaa_polygon_ds.lat
xx,yy = np.meshgrid(longi_n,lati_n)
lat = np.array(lati_n)
lon = np.array(longi_n)


# In[7]:


indexlist  = [1527,1944,3285,4839,5376,5847,6223,6473,6748,7656,7727,8098,8432,8678,8860,9071]


# In[ ]:





# In[ ]:





# In[ ]:





# In[8]:


date_range = pd.date_range(start='1993-01-01', end='2019-12-31', freq='D')


# In[ ]:





# In[ ]:





# In[ ]:


# create function to run the script


# In[ ]:





# In[9]:


def easy_plot_single_test(ax,lon,lat,uneven_levels,data,colors_sch, frame, date_range, index):
    uneven_levels = uneven_levels
    cmap_rb = plt.get_cmap(colors_sch)#'turbo')#'RdBu_r')
    colors = cmap_rb(np.linspace(0, 1, (len(uneven_levels) - 1)))
    cmap, norm = mcolors.from_levels_and_colors(uneven_levels, colors)
    
    im = ax.contourf(lon,lat,data[frame,:,:], 
               levels=uneven_levels ,           
               cmap= colors_sch,#"RdBu_r",#"plasma",
               transform=data_crs,
               transform_first=True)
    
    ax.set_xlim(-35,29)
    ax.set_ylim(-30.875,-2.375)
    datet = date_range[index+frame-150].strftime('%Y-%m-%d')
    ax.set_title(f'max area occurence NOAAOISST {datet}', fontsize=14)


# In[10]:


data_crs = ccrs.PlateCarree(central_longitude=0)

# But you want the map to include the dateline.
proj_crs = ccrs.PlateCarree(central_longitude=180)


# In[11]:


eez = gpd.read_file("/home/shilpa/Downloads/World_EEZ_v11_20191118/eez_v11.shp")


# In[12]:


def update(frame, ax):
    ax.clear()
    #eez[:].plot(ax=ax,color='none',alpha=0.3, edgecolor='black',transform=data_crs)
    ax.coastlines()
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
    #ax.set_xlim(-35,29)
    #ax.set_ylim(-30.875,-2.375)
    easy_plot_single_test(ax = ax, lon=xx, lat=yy, data=test_max_area_evolution, frame=frame, uneven_levels=uneven_levels, colors_sch='turbo', date_range=date_range, index=index)
  


# In[13]:


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


# In[14]:


uneven_levels = [0.0000001,10,20,30,40,50,60,70,80,90,100,200,300,400,500,550,600,650,700,750,800,850]


# In[ ]:


#make the figures


# In[15]:


for index in indexlist[2:3]:
    ds = xr.open_dataset(f'noaa_mask_area_evolution_{index}.nc')
    test_max_area_evolution = np.array(ds.maxarea_evolution)
 
    # Create a figure
    fig, ax = plt.subplots(subplot_kw={'projection': proj_crs}, figsize=(20, 7))
    colors_sch='turbo'
    uneven_levels = uneven_levels
    cmap_rb = plt.get_cmap(colors_sch)#'turbo')#'RdBu_r')
    colors = cmap_rb(np.linspace(0, 1, (len(uneven_levels) - 1)))
    cmap, norm = mcolors.from_levels_and_colors(uneven_levels, colors)

    im = ax.contourf(xx,yy,test_max_area_evolution[0,:,:], 
                   levels=uneven_levels ,           
                   cmap= colors_sch,#"RdBu_r",#"plasma",
                   transform=data_crs,
                   transform_first=True)

    cbar = fig.colorbar(im, ax=ax,orientation="horizontal", shrink=0.5)
    cbar.set_label('[square degrees]')

    #eez[:].plot(ax=ax,color='none',alpha=0.3, edgecolor='black',transform=data_crs)
    ax.coastlines()
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

    ax.set_xlim(-35,29)
    ax.set_ylim(-30.875,-2.375)

    # Create a FuncAnimation object
    animation = FuncAnimation(fig, update, frames=301, interval=100, fargs=(ax,))  # 100 frames, 100ms per frame

    # Create a directory to save the frames as PNG images
    output_dir = f'test_frames_maxarea_{index}'
    os.makedirs(output_dir, exist_ok=True)

    # Save the frames as individual PNG images
    for frame in range(301):
        update(frame, ax)
        plt.savefig(f'{output_dir}/frame_{frame:03d}.png', dpi=100)

    # Use FFmpeg to create an MP4 video from the PNG frames
    os.system(f'ffmpeg -framerate 10 -i test_frames_maxarea_{index}/frame_%03d.png -c:v libx264 -r 30 -pix_fmt yuv420p mhw_noaaoisst_maxarea_{index}.mp4')

    # Clean up the temporary frame files
    #for frame in range(301):
        #os.remove(f'{output_dir}/frame_{frame:03d}.png')

    # Optionally, remove the frames directory
    #os.rmdir(output_dir)

    # Show the video in Matplotlib (optional)
    # You can also save the animation as a GIF or other formats if needed.
    plt.close(fig)


# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:




