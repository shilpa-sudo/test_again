#!/usr/bin/env python
# coding: utf-8

# In[1]:


# in this notebook we plot glorys, noaa, polygon areas overlayed with sun et al. output from zizhie


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


# In[3]:


# open polygon areas for noaa, glorys and sun et al method


# In[4]:


sun_mhw_areas = xr.open_dataset('/media/shilpa/Expansion/mmhwfiji.nc')


# In[5]:


sun_mhw_areas


# In[7]:


datessfull = pd.date_range(start='1982-01-01', end='2022-12-31', freq='D')
len(datessfull) # sun_mhw_areas are from 1982 tO 2022 but we need it from 1993 to 2019 to compare with GLORYS


# In[12]:


date_to_find = '1993-01-01'
s_index_of_date = datessfull.get_loc(date_to_find)
print(s_index_of_date)


# In[13]:


date_to_find = '2020-01-01'
e_index_of_date = datessfull.get_loc(date_to_find)
print(e_index_of_date)


# In[14]:


#start_time = '1993-01-01'  
#end_time = '2019-12-31' 
sun_1993_2019_areas =  np.array(sun_mhw_areas.mmhw_ts.sel(time=slice(s_index_of_date, e_index_of_date)))


# In[15]:


noaa_areas =  xr.open_dataset('noaaoisst_full_mhw_areas_attempt.nc')


# In[16]:


noaa_areas


# In[17]:


start_time = '1993-01-01'  # Replace with your actual start time
end_time = '2019-12-31' 
noaa_1993_2019_areas =  np.array(noaa_areas.mhw_spatial_extent_area.sel(time=slice(start_time, end_time)))


# In[18]:


glorys_areas =  xr.open_dataset('glorys_full_mhw_areas_attempt.nc')


# In[19]:


glorys_areas


# In[20]:


g_areas = np.array(glorys_areas.mhw_spatial_extent_area)


# In[ ]:


# for each day get polygon areas, different hatching style for glorys and noaa, single color in light green
# for sun et al


# In[21]:


data_crs = ccrs.PlateCarree(central_longitude=0)

# But you want the map to include the dateline.
proj_crs = ccrs.PlateCarree(central_longitude=180)


# In[22]:


eez = gpd.read_file("/home/shilpa/Downloads/World_EEZ_v11_20191118/eez_v11.shp")


# In[23]:


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


# In[24]:


longi = glorys_areas.lon
lati = glorys_areas.lat
xx,yy = np.meshgrid(longi,lati)


# In[25]:


datessfull = pd.date_range(start='1993-01-01', end='2019-12-31', freq='D')


# In[26]:


frame = 8000


# In[28]:


uneven_levels = [0.0000001,10,20,30,40,50,60,70,80,90,100,200,300,400,500,550,600,650,700,750,800,850]


# In[29]:


masked_sun = xr.where(sun_1993_2019_areas > 0, 1, np.nan)


# In[38]:


output_dir = 'test_frames_mhwpolyareas'
#os.makedirs(output_dir)


# In[ ]:


for frame in range(9860,9861):
    
    fig, ax = plt.subplots(subplot_kw={'projection': proj_crs}, figsize=(13, 5))

    colors_sch='turbo'
    uneven_levels = uneven_levels
    cmap_rb = plt.get_cmap(colors_sch)#'turbo')#'RdBu_r')
    colors = cmap_rb(np.linspace(0, 1, (len(uneven_levels) - 1)))
    cmap, norm = mcolors.from_levels_and_colors(uneven_levels, colors)

    im = ax.contourf(xx[1:-1,1:-1],yy[1:-1,1:-1],g_areas[frame,1:-1,1:-1],
                   levels=uneven_levels ,     
                   cmap= colors_sch,alpha = 0.5,
                   transform=data_crs,
                   transform_first=True)

    im_= ax.contourf(xx[1:-1,1:-1],yy[1:-1,1:-1],g_areas[frame,1:-1,1:-1],
                   levels=uneven_levels , hatches=['/'],     
                   cmap= colors_sch,alpha = 0.5,
                   transform=data_crs,
                   transform_first=True)


    im_= ax.contourf(xx[1:-1,1:-1],yy[1:-1,1:-1],noaa_1993_2019_areas[frame,1:-1,1:-1],
                   levels=uneven_levels , hatches=['-'],     
                   cmap= colors_sch,alpha = 0.5,
                   transform=data_crs,
                   transform_first=True)

    im2= ax.contourf(xx[1:-1,1:-1],yy[1:-1,1:-1],masked_sun[frame,:,:],     
                   cmap= 'Greens_r', alpha = 0.5,
                   transform=data_crs,
                   transform_first=True)


    cbar = fig.colorbar(im, ax=ax,orientation="horizontal", shrink=0.5)
    cbar.set_label('[square degrees]')


    eez[0:1].plot(ax=ax,color='none',alpha=0.3, edgecolor='black',transform=data_crs)
    eez[2:3].plot(ax=ax,color='none',alpha=0.3, edgecolor='black',transform=data_crs)
    eez[4:5].plot(ax=ax,color='none',alpha=0.3, edgecolor='black',transform=data_crs)
    eez[7:9].plot(ax=ax,color='none',alpha=0.3, edgecolor='black',transform=data_crs)
    eez[11:12].plot(ax=ax,color='none',alpha=0.3, edgecolor='black',transform=data_crs)
    eez[15:16].plot(ax=ax,color='none',alpha=0.3, edgecolor='black',transform=data_crs)
    eez[20:21].plot(ax=ax,color='none',alpha=0.3, edgecolor='black',transform=data_crs)
    eez[33:35].plot(ax=ax,color='none', alpha=0.3, edgecolor='black',transform=data_crs)
    eez[37:38].plot(ax=ax,color='none',alpha=0.3, edgecolor='black',transform=data_crs)
    eez[47:48].plot(ax=ax,color='none',alpha=0.3, edgecolor='black',transform=data_crs)
    eez[51:52].plot(ax=ax,color='none',alpha=0.3, edgecolor='black',transform=data_crs)
    eez[53:54].plot(ax=ax,color='none',alpha=0.3, edgecolor='black',transform=data_crs)
    eez[248:250].plot(ax=ax,color='none',alpha=0.3, edgecolor='black',transform=data_crs)
    eez[243:245].plot(ax=ax,color='none',alpha=0.3, edgecolor='black',transform=data_crs)
    eez[246:248].plot(ax=ax,color='none',alpha=0.3, edgecolor='black',transform=data_crs)
    eez[251:253].plot(ax=ax,color='none',alpha=0.3, edgecolor='black',transform=data_crs)
    eez[259:260].plot(ax=ax,color='none',alpha=0.3, edgecolor='black',transform=data_crs)
    eez[261:262].plot(ax=ax,color='none',alpha=0.3, edgecolor='black',transform=data_crs)

    fiji.boundary.plot(ax=ax, color="black",alpha = 0.3,transform=data_crs)
    ncd.boundary.plot(ax=ax, color="black",alpha = 0.3,transform=data_crs)
    vanuatu.boundary.plot(ax=ax, color="black",alpha = 0.3,transform=data_crs)
    wf.boundary.plot(ax=ax, color="black",alpha = 0.3,transform=data_crs)
    solo.boundary.plot(ax=ax, color="black",alpha = 0.3,transform=data_crs)
    tonga.boundary.plot(ax=ax, color="black",alpha = 0.3,transform=data_crs)
    samoa.boundary.plot(ax=ax, color="black",alpha = 0.3,transform=data_crs)
    tuvalu.boundary.plot(ax=ax, color="black",alpha = 0.3,transform=data_crs)
    niue.boundary.plot(ax=ax, color="black",alpha = 0.3,transform=data_crs)
    cooks.boundary.plot(ax=ax, color="black",alpha = 0.3,transform=data_crs)
    tokelau.boundary.plot(ax=ax, color="black",alpha = 0.3,transform=data_crs)
    amsam.boundary.plot(ax=ax, color="black",alpha = 0.3,transform=data_crs)


    gl = ax.gridlines(draw_labels=True)
    gl.top_labels = False
    gl.right_labels = False

    #fig.suptitle('MHW events 1993-2019 (NOAAOISST)', fontsize=14)
    fig.suptitle(datessfull[frame].strftime("%Y-%m-%d"),fontsize=14)
    ax.set_title('MHW polygon araes NOAAOISST (-), GLORYS(/), Sun et al.(green)')
    ax.set_xlim(-35,29)
    ax.set_ylim(-34.875,-2.375)
    plt.tight_layout()
    #plt.savefig('mhwstate_aoi_231197.png')
    plt.savefig(f'{output_dir}/frame_{frame:04d}.png')#, dpi=100)
    plt.close(fig)

# In[ ]:


os.system('ffmpeg -framerate 10 -i test_frames_mhwpolyareas/frame_%04d.png -c:v libx264 -r 30 -pix_fmt yuv420p mhw_glorys_noaa_sun_comparison_93_19.mp4')


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





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:




