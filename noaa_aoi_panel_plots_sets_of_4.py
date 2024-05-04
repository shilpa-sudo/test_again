#!/usr/bin/env python
# coding: utf-8

# In[1]:


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
from shapely.geometry import Polygon as pPolygon
import matplotlib.colors as mcolors


# In[ ]:


maxinten = xr.open_dataset('mhw_3Doutput_maxintensity_noaa.nc')
maxinten_p90 = maxinten.mhw_max_intensity


# In[ ]:


duration = xr.open_dataset('mhw_3Doutput_duration_noaa.nc')
duration_p90 = duration.mhw_duration


# In[ ]:


cumulative_intensity = xr.open_dataset('mhw_3Doutputs_noaa_cumulative_intensity.nc')
cumulative_intensity_p90 = cumulative_intensity.mhw_cum_intensity


# In[ ]:


ds = xr.open_dataset('mhw_intensity_noaa.nc')
intensity3D = ds.mhw_intensity


# In[ ]:


# get our projection sorted

data_crs = ccrs.PlateCarree(central_longitude=0)

# But you want the map to include the dateline.
proj_crs = ccrs.PlateCarree(central_longitude=180)


# In[ ]:


eez = gpd.read_file("/home/shilpa/Downloads/World_EEZ_v11_20191118/eez_v11.shp")


# In[ ]:


longi = ds.lon
lati = ds.lat


# In[ ]:


xx,yy = np.meshgrid(longi,lati)


# In[ ]:


datessfull = pd.date_range(start='01/01/1993', end='31/12/2019', periods=9861)
print(len(datessfull))


# In[ ]:


def easy_plot(tx,lon,lat,data,clevs,i,j,name,units):
    im = axs[i,j].contourf(lon,lat,data, 
               levels=clevs,
               cmap= "turbo",#"RdBu_r",#"plasma",
               transform=data_crs,
               transform_first=True)
    #axs[i,j].coastlines()
    #eez[:].plot(ax=axs[i,j],color='none',edgecolor='black',transform=data_crs)
    
    
    
    eez[0:1].plot(ax=axs[i,j],color='none',edgecolor='black',transform=data_crs)
    eez[2:3].plot(ax=axs[i,j],color='none',edgecolor='black',transform=data_crs)
    eez[4:5].plot(ax=axs[i,j],color='none',edgecolor='black',transform=data_crs)
    eez[7:9].plot(ax=axs[i,j],color='none',edgecolor='black',transform=data_crs)
    eez[11:12].plot(ax=axs[i,j],color='none',edgecolor='black',transform=data_crs)
    eez[15:16].plot(ax=axs[i,j],color='none',edgecolor='black',transform=data_crs)
    eez[20:21].plot(ax=axs[i,j],color='none',edgecolor='black',transform=data_crs)
    eez[33:35].plot(ax=axs[i,j],color='none',edgecolor='black',transform=data_crs)
    eez[37:38].plot(ax=axs[i,j],color='none',edgecolor='black',transform=data_crs)
    eez[47:48].plot(ax=axs[i,j],color='none',edgecolor='black',transform=data_crs)
    eez[51:52].plot(ax=axs[i,j],color='none',edgecolor='black',transform=data_crs)
    eez[53:54].plot(ax=axs[i,j],color='none',edgecolor='black',transform=data_crs)
    eez[248:250].plot(ax=axs[i,j],color='none',edgecolor='black',transform=data_crs)
    eez[243:245].plot(ax=axs[i,j],color='none',edgecolor='black',transform=data_crs)
    eez[246:248].plot(ax=axs[i,j],color='none',edgecolor='black',transform=data_crs)
    eez[251:253].plot(ax=axs[i,j],color='none',edgecolor='black',transform=data_crs)
    eez[259:260].plot(ax=axs[i,j],color='none',edgecolor='black',transform=data_crs)
    eez[261:262].plot(ax=axs[i,j],color='none',edgecolor='black',transform=data_crs)
    
    cbar = plt.colorbar(im, ax=axs[i,j],orientation="horizontal", shrink=0.5)
    cbar.set_label(units)
    gl = axs[i,j].gridlines(draw_labels=True)
    gl.top_labels = False
    gl.right_labels = False
    axs[i,j].set_xlim(-35,29) #-35
    axs[i,j].set_ylim(-34.875,-2.375)
    axs[i,j].set_title(name, fontsize=9) 
    fig.suptitle(datessfull[tx].strftime("%Y-%m-%d"), fontsize=14)
    


# In[ ]:


def easy_plot_unevencolorbars(tx,lon,lat,data,i,j,name,units,uneven_levels):
    uneven_levels = uneven_levels
    cmap_rb = plt.get_cmap('turbo')#'RdBu_r')
    colors = cmap_rb(np.linspace(0, 1, (len(uneven_levels) - 1)))
    cmap, norm = mcolors.from_levels_and_colors(uneven_levels, colors)
    
    im = axs[i,j].contourf(lon,lat,data, 
               levels=uneven_levels,
               cmap=cmap, norm=norm,
               transform=data_crs,
               transform_first=True)
    #axs[i,j].coastlines()
    #eez[:].plot(ax=axs[i,j],color='none',edgecolor='black',transform=data_crs)
    
    
    
    eez[0:1].plot(ax=axs[i,j],color='none',edgecolor='black',transform=data_crs)
    eez[2:3].plot(ax=axs[i,j],color='none',edgecolor='black',transform=data_crs)
    eez[4:5].plot(ax=axs[i,j],color='none',edgecolor='black',transform=data_crs)
    eez[7:9].plot(ax=axs[i,j],color='none',edgecolor='black',transform=data_crs)
    eez[11:12].plot(ax=axs[i,j],color='none',edgecolor='black',transform=data_crs)
    eez[15:16].plot(ax=axs[i,j],color='none',edgecolor='black',transform=data_crs)
    eez[20:21].plot(ax=axs[i,j],color='none',edgecolor='black',transform=data_crs)
    eez[33:35].plot(ax=axs[i,j],color='none',edgecolor='black',transform=data_crs)
    eez[37:38].plot(ax=axs[i,j],color='none',edgecolor='black',transform=data_crs)
    eez[47:48].plot(ax=axs[i,j],color='none',edgecolor='black',transform=data_crs)
    eez[51:52].plot(ax=axs[i,j],color='none',edgecolor='black',transform=data_crs)
    eez[53:54].plot(ax=axs[i,j],color='none',edgecolor='black',transform=data_crs)
    eez[248:250].plot(ax=axs[i,j],color='none',edgecolor='black',transform=data_crs)
    eez[243:245].plot(ax=axs[i,j],color='none',edgecolor='black',transform=data_crs)
    eez[246:248].plot(ax=axs[i,j],color='none',edgecolor='black',transform=data_crs)
    eez[251:253].plot(ax=axs[i,j],color='none',edgecolor='black',transform=data_crs)
    eez[259:260].plot(ax=axs[i,j],color='none',edgecolor='black',transform=data_crs)
    eez[261:262].plot(ax=axs[i,j],color='none',edgecolor='black',transform=data_crs)
    
    cbar = plt.colorbar(im, ax=axs[i,j],orientation="horizontal", shrink=0.5)
    cbar.set_label(units)
    gl = axs[i,j].gridlines(draw_labels=True)
    gl.top_labels = False
    gl.right_labels = False
    axs[i,j].set_xlim(-35,29) #-35
    axs[i,j].set_ylim(-34.875,-2.375)
    axs[i,j].set_title(name, fontsize=9) 
    fig.suptitle(datessfull[tx].strftime("%Y-%m-%d"), fontsize=14)
    


# In[ ]:


for tx in range(9861):

    nrows=2
    ncols=2

    fig, axs = plt.subplots(nrows=nrows,ncols=ncols,subplot_kw={'projection': proj_crs},figsize=(22,10))


    clevs = np.arange(0.00000001, 8, 0.5)
    easy_plot(lon=xx,lat=yy, tx = tx, data=intensity3D[:,:,tx],
                  clevs=clevs,name='Daily MHW Intensity (1993-2019) NOAAOISST',units='[°C]',i=0,j=0)

    uneven_levels = [0.0000001,1,2,3,4,5,6,7,8,9,10,11,12,13,14, 15,16,17,18,19, 20,21,22,23,24, 25,26,27,28,29, 30,32,34,36,38,40,42,44,46,48, 50, 55, 60,65,70,75,80,85,90,95,100,150,200,300,400,500,600,700]
    easy_plot_unevencolorbars(lon = xx,lat=yy,tx = tx, data=cumulative_intensity_p90[:,:,tx],
                  uneven_levels=uneven_levels,name='MHW Cumulative Intensity (1993-2019) NOAAOISST',units='[°C Days]',i=0,j=1)

    uneven_levels = [0.000001,1,2,3,4,5,6,7,8,9, 10,11,12,13,14, 15,16,17,18,19, 20, 22,24,26,28, 30,32,34,36,38,40, 45, 50, 55, 60,65,70,75,80,85,90,95,100,150,200,300,350]
    easy_plot_unevencolorbars(lon=xx,lat=yy,tx = tx, data=duration_p90[:,:,tx],
                  uneven_levels=uneven_levels,name='MHW Duration (1993-2019) NOAAOISST',units='[Days]',i=1,j=0)

    clevs = np.arange(0.00001, 8, 0.5) 
    easy_plot(lon = xx,lat=yy,tx = tx, data=maxinten_p90[:,:,tx],
                  clevs=clevs,name='MHW Max Intensity (1993-2019) NOAAOISST',units='[°C]',i=1,j=1)

    plt.tight_layout()
    
    plt.savefig(f'/home/shilpa/glory_mat_analysis/noaa_aoi_panel_plots/mhw_panelplots_aoi_{tx}.png') 
    plt.close(fig)
    print(f'just finished {tx}')


# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:




