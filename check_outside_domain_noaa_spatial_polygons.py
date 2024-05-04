#!/usr/bin/env python
# coding: utf-8

# In[1]:


# quick script to see when and where mhws are going outside of the domain


# In[2]:


#matplotlib notebook
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


# In[3]:


dtf_p90 = xr.open_dataset('mhw_state_AOI_quarterdegree.nc')
mhwstate_p90 = dtf_p90.mhwstate


# In[4]:


data_crs = ccrs.PlateCarree(central_longitude=0)

# But you want the map to include the dateline.
proj_crs = ccrs.PlateCarree(central_longitude=180)


# In[5]:


eez = gpd.read_file("/home/shilpa/Downloads/World_EEZ_v11_20191118/eez_v11.shp")


# In[6]:


longi = dtf_p90.lon
lati = dtf_p90.lat


# In[7]:


xx,yy = np.meshgrid(longi,lati)


# In[8]:


datessfull = pd.date_range(start='01/01/1993', end='31/12/2019', periods=9861)
print(len(datessfull))


# In[9]:


#for each day get stats on patches:

total_area_list =[]
total_patches_list = []
max_areas_list = []
min_area_list = []

allareaslist = []
allstartlocationslist = []
allpolygonslist = []
allvtcslist = []

for i in range(9861): 
    
    fig = plt.figure(figsize=(13,5))

    ax = plt.axes(projection=proj_crs)

    clevs = np.arange(0.5, 1.5, 0.5) 


    im = ax.contourf(xx[:,:], yy[:,:], mhwstate_p90[:,:,i], 
                       levels=clevs,
                       cmap= "plasma",#"RdBu_r",#"plasma",#plasma
                       transform=data_crs,
                       transform_first=True)

    contours = ax.contour(xx[:,:], yy[:,:], mhwstate_p90[:,:,i], 
                       levels=clevs,
                       cmap= "viridis",#"RdBu_r",#"plasma",#Oranges_r
                       transform=data_crs,
                       transform_first=True) 
    
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
    fig.suptitle(datessfull[i].strftime("%Y-%m-%d"), fontsize=14)
    ax.set_xlim(-35,29)
    ax.set_ylim(-34.875,-2.375)
    plt.tight_layout()

   
    areas = []
    startlocations = []
    polygonslist = []
    vtcslist = []
        
    for contour in contours.collections:
        paths = contour.get_paths()
        for path in paths:
            if len(path.vertices)>=3:
                vtcs = path.vertices
                vtcslist.append(vtcs)
                polygon = Polygon(path.vertices)
                polygonslist.append(polygon)
                area = polygon.area
                areas.append(area)
                x, y = path.vertices[0]
                startlocations.append(path.vertices[0])
                plt.text(x, y, f"{area:.3f}",color = 'blue', fontsize=8, ha="center", va="center")

        
    if  areas == [] or areas == None :  
        total_area_list.append(np.nan)
        total_patches_list.append(np.nan)
        max_areas_list.append(np.nan)
        min_area_list.append(np.nan)

        allareaslist.append(np.nan)
        allstartlocationslist.append(np.nan)
        allpolygonslist.append(np.nan)
        allvtcslist.append(np.nan)
        ax.set_title(f'MHW spatial extent P90, Total area covered [square degrees] = 0,Total no. of patches = 0')  
        
        
    else :
        
        total_area = sum(areas)
        total_patches = len(areas)
        max_area = np.max(areas)
        min_area = np.min(areas)

        total_area_list.append(total_area)
        total_patches_list.append(total_patches)
        max_areas_list.append(max_area)
        min_area_list.append(min_area)

        allareaslist.append(areas)
        allstartlocationslist.append(startlocations)
        allpolygonslist.append(polygonslist)
        allvtcslist.append(vtcslist)

        ax.set_title(f'MHW spatial extent P90, Total area covered [square degrees] = {total_area},Total no. of patches = {total_patches}')
        #ax.set_title('Area in MHW state: %s Percent'% mhws_percentage[i,1])
        #ax.set_title('Percent Area in MHW state: %s (P99), %s (P90), %s (P85), %s (P75), %s (P50)'%(mhws_percentage_p99[i,1], mhws_percentage_p90[i,1] ,mhws_percentage_p85[i,1] ,mhws_percentage_p75[i,1] ,mhws_percentage_p50[i,1]))
    
    #plt.savefig(f'/home/shilpa/glory_mat_analysis/noaa_polygons_spatial/{i}.png')
    #plt.close(fig)
    #print(f'just finished {i}')


# In[10]:


from matplotlib.patches import Polygon
from shapely.geometry import Polygon as ShapelyPolygon
from shapely.ops import cascaded_union


# In[15]:


for j in range (9861):

    fig, ax = plt.subplots()
# Add the polygons and the intersection to the axes

    for i in range(len(allvtcslist[j])):
        ax.add_patch(Polygon(allvtcslist[j][i], facecolor='blue', alpha=0.5))

#ax.add_patch(Polygon(list(union_poly.exterior.coords), facecolor='gray', alpha=0.5))

# Set the limits of the axes
    ax.set_xlim(-40, 40)
    ax.set_ylim(-38, 0)
    fig.suptitle(datessfull[i].strftime("%Y-%m-%d"), fontsize=14)

# Show the plot
    plt.show()

    plt.savefig(f'/home/shilpa/glory_mat_analysis/noaa_polygons_spatial/check_outside_domain_{j}.png')
    plt.close(fig)
    print(f'just finished {j}')


# In[ ]:




