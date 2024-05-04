#!/usr/bin/env python
# coding: utf-8

# In[1]:


#standard imports
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
import marineHeatWaves as mhw
from cartopy import config
import cartopy.crs as ccrs
import xarray as xr
import pandas as pd
import geopandas as gpd
import matplotlib.colors as mcolors


# In[2]:


#get csv file with information
df1 = pd.read_csv('/home/shilpa/glory_mat_analysis/glorys_quarter_depth1500m/mhws_fullset_depth1500m_withdepth.csv')


depth = 1500

# In[3]:


#get eez
eez = gpd.read_file("/home/shilpa/Downloads/World_EEZ_v11_20191118/eez_v11.shp")


# In[4]:


#get projection sorted
data_crs = ccrs.PlateCarree(central_longitude=0)

#get map to include the dateline.
proj_crs = ccrs.PlateCarree(central_longitude=180)


# In[5]:


gdf = gpd.GeoDataFrame(df1, geometry=gpd.points_from_xy(df1.longitude, df1.latitude))


# In[6]:


gdf['date_peak'] = pd.to_datetime(gdf['date_peak'])
gdf['date_start'] = pd.to_datetime(gdf['date_start'])
gdf['date_end'] = pd.to_datetime(gdf['date_end'])


# In[ ]:


#datelist = []
#for i in range(1,10):
#    y = '2019-0%s'%i
#    datelist.append(y)

#for i in range(10,13):
#    x = '2019-%s'%i
#    datelist.append(x)


# In[27]:


datelist = []
for k in range(1993,2020):
    
    for i in range(1,10):
        y = f'{k}-0{i}'
        datelist.append(y)

    for i in range(10,13):
        x = f'{k}-{i}'
        datelist.append(x)


# In[25]:


#datelist = []
#for k in range(1993,2020):
    
#    for i in range(1,10):
#        y = f'{k}-0%s'%i
#        datelist.append(y)

#    for i in range(10,13):
#        x = f'{k}-%s'%i
#        datelist.append(x)


# In[28]:


datelist


# In[12]:


#gdf_2019_dec_peak = gdf[gdf.date_peak.dt.strftime('%Y-%m') == '2019-01'] 


# In[13]:


#gdf_2019_dec_peak


# In[ ]:


for i in datelist[:]:
    gdf_test = gdf[gdf.date_peak.dt.strftime('%Y-%m') == str(i)] 
    
    fig = plt.figure(figsize=(15,5))

    ax = plt.axes(projection=proj_crs)

    uneven_levels = [0, 5, 10, 15, 20, 50, 3000] #(# here you define the levels you want)
    cmap_rb = plt.get_cmap('OrRd')#'RdBu_r')
    colors = cmap_rb(np.linspace(0, 1, (len(uneven_levels) - 1)))
    cmap, norm = mcolors.from_levels_and_colors(uneven_levels, colors)

    cbar = plt.cm.ScalarMappable(norm=norm, cmap=cmap)

    gdf_test.plot(column='intensity_cumulative',ax=ax,transform=data_crs,alpha = 0.3,cmap=cmap, norm = norm,
            legend=True,legend_kwds={'label': "Cumulative Intensity [Degrees C * Days]", 'orientation': "horizontal",'shrink': 0.3})

    eez[:].plot(ax=ax,color='none',alpha=0.5, edgecolor='black',transform=data_crs)

#ncd.boundary.plot(ax=ax, color="blue",alpha = 0.5,transform=data_crs,linewidth = 2)

    gl = ax.gridlines(draw_labels=True)
    gl.top_labels = False
    gl.right_labels = False   
    fig.suptitle(f'Cumulative Intensity of MHWs that peaked in {i}  (GLORYS regridded to quarter degree) {depth}m', fontsize=14)  
    ax.set_xlim(-38,110) 
    ax.set_ylim(-35,0)
    plt.tight_layout()
    plt.savefig(f'/home/shilpa/glory_mat_analysis/vertical_anal_gl_quarter/{i}_{depth}m.png')
    plt.close(fig)
    print(f'just finished {i}')


# In[ ]:




