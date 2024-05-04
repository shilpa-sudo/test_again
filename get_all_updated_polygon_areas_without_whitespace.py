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
from shapely.geometry import Polygon as pPolygon

from matplotlib.patches import Polygon as mplPolygon
from shapely.geometry import Polygon as ShapelyPolygon
from shapely.ops import cascaded_union
from shapely.geometry import MultiPolygon


# In[2]:


dtf_p90 = xr.open_dataset('mhw_state_AOI_quarterdegree.nc')
mhwstate_p90_ = dtf_p90.mhwstate


# In[3]:


mhwstate_p90 = np.zeros((132, 262,9861))


# In[4]:


mhwstate_p90[1:131,1:261,:] = mhwstate_p90_


# In[5]:


data_crs = ccrs.PlateCarree(central_longitude=0)

# But you want the map to include the dateline.
proj_crs = ccrs.PlateCarree(central_longitude=180)


# In[6]:


eez = gpd.read_file("/home/shilpa/Downloads/World_EEZ_v11_20191118/eez_v11.shp")


# In[7]:


longi = dtf_p90.lon
lati = dtf_p90.lat


# In[8]:


newlongi = np.zeros((262))
newlati = np.zeros((132))


# In[9]:


newlati[0] = lati[0] - 0.25
newlati[1:131] = lati[:]
newlati[-1] = lati[-1] + 0.25


# In[10]:


newlongi[0] = longi[0] - 0.25
newlongi[1:261] = longi[:]
newlongi[-1] = longi[-1] + 0.25


# In[11]:


xx,yy = np.meshgrid(newlongi,newlati)


# In[12]:


datessfull = pd.date_range(start='01/01/1993', end='31/12/2019', periods=9861)
print(len(datessfull))


# In[13]:


#val = pd.Timestamp('1998-07-01') # the value you want to find the index for

#dindex = datessfull.get_loc(val)

#print(dindex)


# In[14]:


#tx = dindex


# In[15]:


allpolygons_eachday = []


for tx in range(len(datessfull)):

    


    total_area_list =[]
    total_patches_list = []
    max_areas_list = []
    min_area_list = []

    allareaslist = []
    allstartlocationslist = []
    allpolygonslist = []
    allvtcslist = []








    fig = plt.figure(figsize=(13,5))

    ax = plt.axes(projection=proj_crs)

    clevs = np.arange(0.5, 1.5, 0.5) 


    im = ax.contourf(xx[:,:], yy[:,:], mhwstate_p90[:,:,tx], 
                          levels=clevs,
                           cmap= "plasma",#"RdBu_r",#"plasma",#plasma
                          transform=data_crs,
                           transform_first=True)

    contours = ax.contour(xx[:,:], yy[:,:], mhwstate_p90[:,:,tx], 
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
    fig.suptitle(datessfull[tx].strftime("%Y-%m-%d"), fontsize=14)
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
                polygon = pPolygon(path.vertices)
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
        plt.close(fig)
        #print(f'just finished {i}')
        
    

    polylist = []

    for j in range(len(allvtcslist[0])):

        poly = ShapelyPolygon(allvtcslist[0][j])
        polylist.append(poly)

     # each polygon intersect with all other polygon but its self

    whitepolygonindex = []
    whitepolygonarea = []


    area_poly_without_whitespace_index = []
    area_poly_without_whitespace_list = []

    for i in range(len(polylist)):
        intersectionarea_list = []
        for j in range(len(polylist)):
            if i ==j:
                print('self')

            else:
                intersectarea = polylist[i].intersection(polylist[j]).area
                intersectionarea_list.append(intersectarea)

        if all(elem == 0 for elem in intersectionarea_list):
            print("All elements in the list are zeros")
        else:
            count = 0
            for element in intersectionarea_list:
                if element != 0:
                    count += 1
            print('polygon index = ',i ,'has intersection with',count,'others')

            area_poly_without_whitespace = (polylist[i].area) - (sum(intersectionarea_list))

            print('polygon_area=',polylist[i].area, 'sum_intesection=',sum(intersectionarea_list),
                  'area_poly_without_whitespace=',area_poly_without_whitespace)

            if area_poly_without_whitespace == 0:
                whitepolygonindex.append(i)
                whitepolygonarea.append(0)

            else:
                area_poly_without_whitespace_index.append(i)
                area_poly_without_whitespace_list.append(area_poly_without_whitespace)

    print('no. of polygons=',len(polylist))
    
    print('whitepolygonindex =',whitepolygonindex)
    print('whitepolygonarea=',whitepolygonarea)

    print('area_poly_without_whitespace_index=',area_poly_without_whitespace_index)
    print('area_poly_without_whitespace_list=',area_poly_without_whitespace_list)
    
    realpolygonareas_original = []
    for i in range(len(polylist)):
        area = polylist[i].area
        realpolygonareas_original.append(area)

    print('no. of polygons=',len(realpolygonareas_original))
    
    
    for i in whitepolygonindex:
        index = whitepolygonindex.index(i)
        realpolygonareas_original[i] = whitepolygonarea[index]
    
    print(len(realpolygonareas_original))
    
    for i in area_poly_without_whitespace_index:
        index = area_poly_without_whitespace_index.index(i)
        realpolygonareas_original[i] = area_poly_without_whitespace_list[index]

    print(len(realpolygonareas_original))    
    
    allpolygons_eachday.append(realpolygonareas_original)
    
    count = 0
    for element in realpolygonareas_original:
        if element == 0:
            count += 1
    print('This day has',count,'white spaces inside MHW areas.')  
    
    count = 0
    for element in realpolygonareas_original:
        if element != 0:
            count += 1
    print('This day has',count,' MHW polygons.')
    
    maxrealpolygonareas = max(realpolygonareas_original)
    totalrealpolygonareas = sum(realpolygonareas_original)
    
    totalrealpolygonareasr = round(totalrealpolygonareas, 2)
    maxrealpolygonareasr = round(maxrealpolygonareas, 2)
    
    print('maximum_area_polygon =',maxrealpolygonareasr, 'total MHW area =',totalrealpolygonareasr)
    
    fig = plt.figure(figsize=(13,5))

    ax = plt.axes(projection=proj_crs)

    clevs = np.arange(0.5, 1.5, 0.5) 


    im = ax.contourf(xx[:,:], yy[:,:], mhwstate_p90[:,:,tx], 
                          levels=clevs,
                           cmap= "plasma",#"RdBu_r",#"plasma",#plasma
                          transform=data_crs,
                           transform_first=True)

    contours = ax.contour(xx[:,:], yy[:,:], mhwstate_p90[:,:,tx], 
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
    fig.suptitle(datessfull[tx].strftime("%Y-%m-%d"), fontsize=14)
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
                polygon = pPolygon(path.vertices)
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
        #ax.set_title(f'MHW spatial extent P90, Total area covered [square degrees] = 0,Total no. of patches = 0')  
        ax.set_title(f'MHW spatial extent P90, Total area covered by MHW [square degrees] = 0, No.of MHW patches = 0, Area of largest MHW patch [square degrees] = 0')

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

        ax.set_title(f'Total area covered by MHW(P90) [square degrees] = {totalrealpolygonareasr}, No.of MHW patches = {count}, Area of largest MHW patch [square degrees] = {maxrealpolygonareasr}')
            #ax.set_title('Area in MHW state: %s Percent'% mhws_percentage[i,1])
            #ax.set_title('Percent Area in MHW state: %s (P99), %s (P90), %s (P85), %s (P75), %s (P50)'%(mhws_percentage_p99[i,1], mhws_percentage_p90[i,1] ,mhws_percentage_p85[i,1] ,mhws_percentage_p75[i,1] ,mhws_percentage_p50[i,1]))

        plt.savefig(f'/home/shilpa/glory_mat_analysis/noaa_polygons_spatial_updated/{tx}.png')
        plt.close(fig)
        print(f'just finished {tx}')


# In[ ]:


d = {'time':datessfull, 'polygon_areas':allpolygons_eachday}


# In[17]:


df = pd.DataFrame(data=d)


# In[18]:


df.to_csv('/home/shilpa/glory_mat_analysis/updated_spatial_polygon_areas_noaaoisst')


# In[19]:





# In[20]:





# In[ ]:





# In[26]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[35]:





# In[36]:





# In[ ]:





# In[ ]:




