#!/usr/bin/env python
# coding: utf-8

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


# In[2]:


dtf = Dataset('noaaoisst93_19_sst_AOI_quarterdegree.nc') 
lati= dtf.variables['lati'] #degrees north
longi = dtf.variables['longi'] #degrees east
time = dtf.variables['time']
temp = dtf.variables['temp']


# In[3]:


temparr = np.array(temp)


# In[ ]:


tt = np.array(time)
t = tt.astype(int)


# In[ ]:


#intensitylist = []
#durationlist = []
#maxintensitylist = []
cumulativeintensitylist = []
#rateonsetlist = []
#ratedeclinelist = []


for j in range(130):
    for k in range(260):
        
        sst = temparr[:,j,k]
        
        mhws, clim = mhw.detect(t, sst)
        events = mhws['n_events']
        print(f'location {j},{k}, events = {events}')
        
        if events is None or events == [] or events == () or events == 0:
             y = np.zeros(9861)
             #intensitylist.append(y)
             #durationlist.append(y)
             #maxintensitylist.append(y)
             cumulativeintensitylist.append(y)
             #rateonsetlist.append(y)
             #ratedeclinelist.append(y)
        else:
             
             index_start = mhws['index_start'] #= []
             print(index_start)
             index_end = mhws['index_end']# = []
             print(index_end)
        
             #duration = mhws['duration']
             #intensity_max = mhws['intensity_max']
             intensity_cumulative = mhws['intensity_cumulative']
             #rateonset = mhws['rate_onset']
             #ratedecline = mhws['rate_decline']




             #mhw_intensity_list = []
             #mhw_duration_list = []
             #mhw_maxintensity_list = []
             mhw_cumulativeintensity_list = []
             #mhw_rateonsetlist = []
             #mhw_ratedeclinelist = []
             
             #x = np.zeros(9861)
             #x1 = np.zeros(9861)
             #x2 = np.zeros(9861)
             x3 = np.zeros(9861)
             #x4 = np.zeros(9861)
             #x5 = np.zeros(9861)
            
             
             for i in range(events):    
                
               # x[index_start[i]:index_end[i]]= ((sst[index_start[i]:index_end[i]]) -  (clim['seas'][index_start[i]:index_end[i]]))
                #mhw_intensity_list.append(x)
                
                #x1[index_start[i]:index_end[i]]= duration[i]
                #mhw_duration_list.append(x1)
        
                #x2[index_start[i]:index_end[i]]= intensity_max[i]
                #mhw_maxintensity_list.append(x2)
                
                x3[index_start[i]:index_end[i]]= intensity_cumulative[i]
                mhw_cumulativeintensity_list.append(x3)
                
                #x4[index_start[i]:index_end[i]]= rateonset[i]
                #mhw_rateonsetlist.append(x4)
                
                #x5[index_start[i]:index_end[i]]= ratedecline[i]
                #mhw_ratedeclinelist.append(x5)
                
               
                
             #intensitylist.append(mhw_intensity_list[0])
             #durationlist.append(mhw_duration_list[0])
             #maxintensitylist.append(mhw_maxintensity_list[0])
             cumulativeintensitylist.append(mhw_cumulativeintensity_list[0])
             #rateonsetlist.append(mhw_rateonsetlist[0])
             #ratedeclinelist.append(mhw_ratedeclinelist[0])


# In[ ]:


#iarr = np.array(intensitylist)

#intensity3D = iarr.reshape(130,260,9861)

#darr = (np.array(durationlist)).reshape(130,260,9861)

#maxintenarr =  (np.array(maxintensitylist)).reshape(130,260,9861)
 
cumintenarr =  (np.array(cumulativeintensitylist)).reshape(130,260,9861)    
    
#rateonsetarr =  (np.array(rateonsetlist)).reshape(130,260,9861)      
 
#ratedecarr = (np.array(ratedeclinelist)).reshape(130,260,9861) 
            


# In[ ]:


datessfull = [date.fromordinal(tt.astype(int)) for tt in t]
print(datessfull)


# In[ ]:





# In[ ]:


#save output as netcdf file
# save the climatology ouput as a netcdf file


try: ncfile.close()  # just to be safe, make sure dataset is not already open.
except: pass
ncfile = Dataset('mhw_3Doutputs_noaa_cumulative_intensity.nc',mode='w',format='NETCDF4_CLASSIC') 
print(ncfile)
lat_dim = ncfile.createDimension('lat', 130)     # latitude axis
lon_dim = ncfile.createDimension('lon', 260)    # longitude axis
time_dim = ncfile.createDimension('time', 9861) # time axis 

ncfile.title='nooaoisst mhw cumulative intensity 1993-2019'
ncfile.subtitle="cumulative intensity"

lat = ncfile.createVariable('lat', np.float32, ('lat',))
lat.units = 'degrees_north'
lat.long_name = 'latitude'
lon = ncfile.createVariable('lon', np.float32, ('lon',))
lon.units = 'degrees_east'
lon.long_name = 'longitude'
time = ncfile.createVariable('time', np.float64, ('time',))
time.units = 'days_from_ordinal'
time.long_name = 'time'
# Define a 3D variable to hold the data

#temp = ncfile.createVariable('mhw_intensity',np.float64,('lat','lon','time')) # note: unlimited dimension is leftmost
#temp.units = 'degrees celsius'
#temp.standard_name = 'mhw intensity' # this is a CF standard name

#temp1 = ncfile.createVariable('mhw_duration',np.float64,('lat','lon','time')) # note: unlimited dimension is leftmost
#temp1.units = 'days_in_mhw'
#temp1.standard_name = 'mhw duration' # this is a CF standard name

#temp2 = ncfile.createVariable('mhw_max_intensity',np.float64,('lat','lon','time')) # note: unlimited dimension is leftmost
#temp2.units = 'degrees celsius'
#temp2.standard_name = 'mhw max intensity' # this is a CF standard name

temp3 = ncfile.createVariable('mhw_cum_intensity',np.float64,('lat','lon','time')) # note: unlimited dimension is leftmost
temp3.units = 'degrees celsius days'
temp3.standard_name = 'mhw cumulative intensity' # this is a CF standard name

#temp4 = ncfile.createVariable('mhw_rateonset',np.float64,('lat','lon','time')) # note: unlimited dimension is leftmost
#temp4.units = 'degrees celsius per day'
#temp4.standard_name = 'mhw rate onset' # this is a CF standard name

#temp5 = ncfile.createVariable('mhw_ratedecline',np.float64,('lat','lon','time')) # note: unlimited dimension is leftmost
#temp5.units = 'degrees celsius per day'
#temp5.standard_name = 'mhw rate decline' # this is a CF standard name


lat[:] = lati[:]
lon[:] = longi[:]

#temp[:,:,:] =  intensity3D # Appends data along unlimited dimension
#temp1[:,:,:] =  darr
#temp2[:,:,:] =  maxintenarr # Appends data along unlimited dimension
temp3[:,:,:] =  cumintenarr
#temp4[:,:,:] =  rateonsetarr # Appends data along unlimited dimension
#temp5[:,:,:] =  ratedecarr

time[:] = t
# first print the Dataset object to see what we've got
print(ncfile)
# close the Dataset.
ncfile.close(); print('Dataset is closed!')


# In[ ]:





# In[ ]:





# In[ ]:




