#!/usr/bin/env python
# coding: utf-8

# In[1]:



import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.dates import num2date
from matplotlib import dates as mdates
mpl.use('Agg') #nonguifor use on datarmor
import cartopy.crs as ccrs
import xarray as xr
import pandas as pd
import matplotlib.animation as animation
import geopandas as gpd
from matplotlib.patches import Polygon as mplPolygon
from shapely.geometry import Polygon as ShapelyPolygon
from shapely.geometry import Point as ShapelyPoint
from shapely.ops import cascaded_union
from shapely.geometry import MultiPolygon


# In[2]:


dtf_p90 = xr.open_dataset('/home/shilpa/glory_mat_analysis/noaaoisst_lastday_included/full_spatial_noaaoistt_aoi_lastdayincluded.nc')
mhwstate_p90_ = dtf_p90.spatial_extent
mhwstate_p90 = np.zeros((132, 262,9861)) #trying to close all polygons so adding zeros on the sides
mhwstate_p90[1:131,1:261,:] = mhwstate_p90_


# In[3]:


longi = dtf_p90.lons
lati = dtf_p90.lats


# In[4]:


newlongi = np.zeros((262))
newlati = np.zeros((132))
newlati[0] = lati[0] - 0.25  # making new lat by adding quarter degree
newlati[1:131] = lati[:]
newlati[-1] = lati[-1] + 0.25
newlongi[0] = longi[0] - 0.25 # making new lon by adding quarter degree
newlongi[1:261] = longi[:]
newlongi[-1] = longi[-1] + 0.25
xx,yy = np.meshgrid(newlongi,newlati)


# In[5]:


datessfull = pd.date_range(start='01/01/1993', end='12/31/2019', periods=9861)
#print(len(datessfull))


# In[64]:


for tx in range (len(datessfull[:2])):

    spatial_with_areas = np.zeros((1,132, 262))

    fig = plt.figure(figsize=(10,5))
    contours = plt.contour(xx[:,:], yy[:,:], mhwstate_p90[:,:,tx])

    areas = []
    polygonslist = []
    vtcslist = []

    paths = contours.collections[0].get_paths() #get the first contour polygon

    for path in paths:
        if len(path.vertices)>=3: 
            vtcs = path.vertices
            vtcslist.append(vtcs)
            polygon = ShapelyPolygon(path.vertices)
            polygonslist.append(polygon)
            print(polygon)
            area = polygon.area
            print(area)
            areas.append(area)

    plt.close(fig)

    if  areas == [] or areas == None : #means there are no mhw polygons on this day
        print(f'no mhw events on this day {datessfull[tx]}')

        xrds = xr.Dataset(coords = dict(time = [datessfull[tx]],lat = newlati,lon = newlongi),data_vars = dict(mhw_spatial_extent_area = (['time','lat','lon'],spatial_with_areas)))

        #xrds.to_netcdf(f'/home2/datahome/slal/mhw_polygon_areas_{tx}.nc')
        xrds.to_netcdf(f'noaaoisst_mhw_polygon_areas_{tx}.nc')

    else: #find real area of polygon

        whitepolygonindex = []
        whitepolygonarea = []

        area_poly_without_whitespace_index = []
        area_poly_without_whitespace_list = []

        for i in range(len(polygonslist)):
            intersectionarea_list = []
            for j in range(len(polygonslist)):
                if i ==j:
                    print('self')
                    # each polygon intersect with all other polygon but its self

                else:
                    intersectarea = polygonslist[i].intersection(polygonslist[j]).area
                    intersectionarea_list.append(intersectarea)

            if all(elem == 0 for elem in intersectionarea_list):
                print("All elements in the list are zeros")  #there is no intersection, there r no white spaces inside mhw areas
            else:
                count = 0
                for element in intersectionarea_list:
                    if element != 0:
                        count += 1
                print('polygon index = ',i ,'has intersection with',count,'others')

                area_poly_without_whitespace = (polygonslist[i].area) - (sum(intersectionarea_list))  #remove white space from polygon area

                print('polygon_area=',polygonslist[i].area, 'sum_intesection=',sum(intersectionarea_list),
                      'area_poly_without_whitespace=',area_poly_without_whitespace)

                if area_poly_without_whitespace == 0: #this applies only to white polygons inside mhw polygons because when they insect mhw polygons, the area of intersection is same as the size of the white polygon, hence the area difference comes to zero.
                    whitepolygonindex.append(i)
                    whitepolygonarea.append(0)

                else:
                    area_poly_without_whitespace_index.append(i)
                    area_poly_without_whitespace_list.append(area_poly_without_whitespace)




        print('whitepolygonindex =',whitepolygonindex)
        print('whitepolygonarea=',whitepolygonarea)

        print('area_poly_without_whitespace_index=',area_poly_without_whitespace_index)
        print('area_poly_without_whitespace_list=',area_poly_without_whitespace_list)

        realpolygonareas_original = []

        for i in range(len(polygonslist)):
            area = polygonslist[i].area
            realpolygonareas_original.append(area)

        print('no. of polygons=',len(realpolygonareas_original))


        for i in whitepolygonindex:
            index = whitepolygonindex.index(i)
            realpolygonareas_original[i] = whitepolygonarea[index] #put zero where white polygons are

        print(len(realpolygonareas_original))

        for i in area_poly_without_whitespace_index:
            index = area_poly_without_whitespace_index.index(i)
            realpolygonareas_original[i] = area_poly_without_whitespace_list[index]#put real area for mhw polygons

        print(len(realpolygonareas_original))    

        #allpolygons_eachday.append(realpolygonareas_original)

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

        for l in range (132):
                    for m in range (262):
                        if mhwstate_p90[l,m,tx] != 0.0:
                            point = ShapelyPoint(newlongi[m],newlati[l]) #turn each location in mhw state i.e. 1 to a point
                            transformed_point = point
                            for pp in range (len(polygonslist)):
                                is_inside = polygonslist[pp].contains(transformed_point)
                                print(is_inside)
                                if  is_inside == True:
                                    print(f'this point lat ={l}, lon = {m}, insects polygon no. {pp}, with real area = {realpolygonareas_original[pp]}')
                                    spatial_with_areas[:,l,m] = realpolygonareas_original[pp] 

        xrds = xr.Dataset(coords = dict(time = [datessfull[tx]],lat = newlati,lon = newlongi),data_vars = dict(
                       mhw_spatial_extent_area = (['time','lat','lon'],spatial_with_areas)))

        #xrds.to_netcdf(f'/home2/datahome/slal/mhw_polygon_areas_{tx}.nc')

        xrds.to_netcdf(f'noaaoisst_mhw_polygon_areas_{tx}.nc')






# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:




