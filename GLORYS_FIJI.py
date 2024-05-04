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


# In[2]:


fnames = glob.glob('/media/shilpa/LaCie/GLORYS12V1/*gridT_*')
fnames.sort()


# In[3]:


filelist = []
for i in fnames:
    filelist.append(i)


# In[4]:


only_ncfiles = []
for each_file in filelist[2:-152]:
    if each_file.endswith('.nc'):  
        print (each_file)
        only_ncfiles.append(each_file)


# In[5]:


#lets get time
tlist = []

for eachfile in only_ncfiles:
    print(eachfile)
    d=Dataset(eachfile)

    tlist.append(d.variables['time_counter'][:])
    d.close()

time=np.hstack(tlist) 


# In[6]:


t_days= time.astype(np.float64)/24
print(t_days)

t_mpl=t_days+datetime.date(1950,1,1).toordinal()
print(t_mpl)

datess = [date.fromordinal(tt.astype(int)) for tt in t_mpl]
print(datess)


# In[7]:


#get dataset
dt = Dataset('/media/shilpa/LaCie/GLORYS12V1/GLORYS12V1_1dAV_19930814_19930815_gridT_R19930818.nc')
lat= dt.variables['nav_lat'] #degrees north  :500:5,500:1200:5]
lon = dt.variables['nav_lon']#[::5]#degrees east   :500:5,500:1200:5]
#time = dt.variables['time']
temp_test = dt.variables['votemper'][0,0,:,:]


# In[8]:


data_crs = ccrs.PlateCarree(central_longitude=0)

# But you want the map to include the dateline.
proj_crs = ccrs.PlateCarree(central_longitude=180)


# In[9]:


eez = gpd.read_file("/home/shilpa/Downloads/World_EEZ_v11_20191118/eez_v11.shp")


# In[10]:


fig = plt.figure(figsize=(10,4))

ax = plt.axes(projection=proj_crs)

clevs = np.arange(0, 35, 1) 
im=ax.contourf(lon[110:335,560:710], lat[110:335,560:710], temp_test[110:335,560:710], 
               levels=clevs,
               cmap= "coolwarm",#"RdBu_r",#"plasma",
               transform=data_crs,
               transform_first=True)  # This kwarg makes it much faster.
cbar = plt.colorbar(im, ax=ax, orientation="horizontal", shrink=0.5)
cbar.set_label(r'Temperature [$^\circ$C]')
eez[:].plot(ax=ax,color='none',edgecolor='black',transform=data_crs)
ax.set_xlim(-9,6) #-35
ax.set_ylim(-29,-8)
#ax.coastlines()
gl = ax.gridlines(draw_labels=True)
gl.xlabels_top = False
gl.ylabels_right = False
#fig.suptitle('MHW events 1993-2019 (NOAAOISST)', fontsize=14)
plt.tight_layout()


# In[11]:


#get the lats and longs for fiji, and then run mhw block analysis on it


# In[12]:


templist = []

for eachfile in only_ncfiles:
    print(eachfile)
    d=Dataset(eachfile)

    templist.append(d.variables['votemper'][0,0,110:335,560:710])  #only for AOI at depth 0
    d.close()
    


# In[13]:


temp=np.ma.concatenate(templist)


# In[14]:


temparr = np.array(temp)
tt = np.array(t_mpl)
t = tt.astype(int)


# In[15]:


print(lon[110:335,560:710].shape, lat[110:335,560:710].shape)


# In[ ]:





# In[17]:


temparr_rs = temparr.reshape(9861,225,150)


# In[ ]:


yearscenterlist = []
countlist = []
intensitymaxlist = []
durationlist = []
intensity_cumulativelist = []
rateonsetlist = []
ratedeclinelist = []
totaldayslist =[]
total_cum_int_list =[]

for j in range(225):
    for k in range(150):

        mhws, clim = mhw.detect(t, temparr_rs[:,j,k])
        mhwBlock = mhw.blockAverage(t, mhws)
        
        yearscentre = mhwBlock['years_centre']
        count = mhwBlock['count']
        mean_max_intensity = mhwBlock['intensity_max']
        duration = mhwBlock['duration']
        cum_intensity = mhwBlock['intensity_cumulative']
        rate_onset = mhwBlock['rate_onset']
        rate_decline = mhwBlock['rate_decline']
        total_days = mhwBlock['total_days']
        total_cum_intensity = mhwBlock['total_icum']
        
        yearscenterlist.append(yearscentre)
        countlist.append(count)
        intensitymaxlist.append(mean_max_intensity)
        durationlist.append(duration)
        intensity_cumulativelist.append(cum_intensity)
        rateonsetlist.append(rate_onset)
        ratedeclinelist.append(rate_decline)
        totaldayslist.append(total_days)
        total_cum_int_list.append(total_cum_intensity)
        
        
yearscenterseries = np.array(yearscenterlist)
countlistseries = np.array(countlist)
intensitymaxseries = np.array(intensitymaxlist)
durationseries = np.array(durationlist)
intensity_cumulative_series = np.array(intensity_cumulativelist)
rateonset_series = np.array(rateonsetlist)
ratedecline_series = np.array(ratedeclinelist)
totaldays_series = np.array(totaldayslist)
total_cum_int_series = np.array(total_cum_int_list)


# In[18]:


years = yearscenterseries.reshape((225,150,27))
countevent = countlistseries.reshape((225,150,27))
inte = intensitymaxseries.reshape((225,150,27))
durat = durationseries.reshape((225,150,27)) 
intensi_cum = intensity_cumulative_series.reshape((225,150,27)) 
ronset = rateonset_series.reshape((225,150,27)) 
rdec = ratedecline_series.reshape((225,150,27)) 
tdays = totaldays_series.reshape((225,150,27)) 
totalcumint = total_cum_int_series.reshape((225,150,27)) 


# In[ ]:


yrs = np.array([1993., 1994., 1995., 1996., 1997., 1998., 1999., 2000., 2001.,
       2002., 2003., 2004., 2005., 2006., 2007., 2008., 2009., 2010.,
       2011., 2012., 2013., 2014., 2015., 2016., 2017., 2018., 2019.])


# In[ ]:


# save ouput as a netcdf file.....this one contains 14 variables

try: ncfile.close()  # just to be safe, make sure dataset is not already open.
except: pass
ncfile = Dataset('glorys_block_properties_fiji.nc',mode='w',format='NETCDF4_CLASSIC') 
print(ncfile)
y_dim = ncfile.createDimension('y', 225)     # latitude axis
x_dim = ncfile.createDimension('x', 150)    # longitude axis
time_dim =ncfile.createDimension('time_center', 27)

ncfile.title='glorys AOI 1993-2019'
ncfile.subtitle="MHW analysis for potential temperature depth 0m"

lati = ncfile.createVariable('lati', np.float32, ('y','x'))
lati.units = 'degrees_north'
lati.long_name = 'latitude'

longi = ncfile.createVariable('longi', np.float32, ('y','x'))
longi.units = 'degrees_east'
longi.long_name = 'longitude'

time_center = ncfile.createVariable('time_center', np.float32, ('time_center'))
time_center.units = 'years_'
time_center.long_name = 'years_center'

var22 = ncfile.createVariable('countevent',np.float64,('y','x','time_center')) # note: unlimited dimension is leftmost
var22.units = 'no.of events' # degrees celsius
var22.standard_name = 'No.of events' # this is a CF standard name

var23 = ncfile.createVariable('intensity',np.float64,('y','x','time_center')) # note: unlimited dimension is leftmost
var23.units = 'Deg.C' # degrees celsius
var23.standard_name = 'mean_of_max_intensity' # this is a CF standard name

var24 = ncfile.createVariable('duration_by_year',np.float64,('y','x','time_center')) # note: unlimited dimension is leftmost
var24.units = 'Days' # degrees celsius
var24.standard_name = 'Average_MHW_duration_by_year' # this is a CF standard name

var25 = ncfile.createVariable('cum_intensity_by_year',np.float64,('y','x','time_center')) # note: unlimited dimension is leftmost
var25.units = 'Days_deg.C' # degrees celsius
var25.standard_name = 'Average_MHW_cum_intensity_by_year' # this is a CF standard name

var26 = ncfile.createVariable('rate_onset_by_year',np.float64,('y','x','time_center')) # note: unlimited dimension is leftmost
var26.units = 'deg.C/day' # degrees celsius
var26.standard_name = 'Average_MHW_onset_rate_by_year' # this is a CF standard name

var27 = ncfile.createVariable('declinerate',np.float64,('y','x','time_center')) # note: unlimited dimension is leftmost
var27.units = 'deg.C/day' # degrees celsius
var27.standard_name = 'Average_MHW_decline_rate_by_year' # this is a CF standard name

var28 = ncfile.createVariable('total_mhw_days',np.float64,('y','x','time_center')) # note: unlimited dimension is leftmost
var28.units = 'total_days' # degrees celsius
var28.standard_name = 'Total_number_of_mhw_days_by_year' # this is a CF standard name


var29 = ncfile.createVariable('total_cumulative_intensity',np.float64,('y','x','time_center')) # note: unlimited dimension is leftmost
var29.units = 'degc*days' # degrees celsius
var29.standard_name = 'Total_cumulative_intensity_by_year' # this is a CF standard name

lati[:,:] = lat[110:335,560:710]
longi[:,:] = lon[110:335,560:710]
time_center[:] = yrs[:]


var22[:,:,:] =  countevent
var23[:,:,:] = inte 
var24[:,:,:] = durat 
var25[:,:,:] = intensi_cum 
var26[:,:,:] = ronset
var27[:,:,:] = rdec 
var28[:,:,:] = tdays 
var29[:,:,:] = totalcumint  

# first print the Dataset object to see what we've got
print(ncfile)
# close the Dataset.
ncfile.close(); print('Dataset is closed!')

