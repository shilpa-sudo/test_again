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


# In[2]:


ds = xr.open_dataset('mhw_intensity_noaa.nc')
intensity3D = ds.mhw_intensity


# In[3]:


# get our projection sorted

data_crs = ccrs.PlateCarree(central_longitude=0)

# But you want the map to include the dateline.
proj_crs = ccrs.PlateCarree(central_longitude=180)


# In[5]:


eez = gpd.read_file("/home/shilpa/Downloads/World_EEZ_v11_20191118/eez_v11.shp")


# In[6]:


lat = np.array(ds.lat)
lon = np.array(ds.lon)
intensity3D = np.array(ds.mhw_intensity)


xx,yy = np.meshgrid(lon,lat)
print(xx.shape,yy.shape)


tm = np.array(ds.time)
t = tm.astype(int)
datessfull = [date.fromordinal(tt.astype(int)) for tt in tm]
print(datessfull)


# In[7]:


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


# In[8]:


df = pd.read_csv('lonlist_clip_eez.csv') 
dfarrno = df.to_numpy()
dfarrno.shape


# In[10]:


nlon = dfarrno[:,1]


# In[14]:


edgar = ds.mhw_intensity[:,:,2007].to_dataframe()
edgar['newlons'] = nlon
edgar = edgar.reset_index()
edgar_gdf = gpd.GeoDataFrame(edgar, geometry=gpd.points_from_xy(edgar.newlons, edgar.lat))
edgar_gdf.crs = ncd.crs
ncd_data = gpd.clip(edgar_gdf, ncd)
max_temp = ncd_data ['mhw_intensity'].max()
mean_temp = ncd_data ['mhw_intensity'].mean()
#ncdarr = np.array(ncd_data)
#nanmean_ncd = (np.nanmean(ncdarr[:,3]))
#print(nanmean_ncd)
print(max_temp)
print(mean_temp)


# In[22]:


(((8*9861)/60)/60)/10


# In[25]:


def country_mean_intensity(country):    
    edgar_gdf.crs = country.crs
    country_data = gpd.clip(edgar_gdf, country)
    #max_temp = ncd_data ['mhw_intensity'].max()
    mean_temp = ncd_data ['mhw_intensity'].mean()
    return mean_temp


# In[26]:


def country_max_intensity(country):    
    edgar_gdf.crs = country.crs
    country_data = gpd.clip(edgar_gdf, country)
    max_temp = ncd_data ['mhw_intensity'].max()
    #mean_temp = ncd_data ['mhw_intensity'].mean()
    return max_temp


# In[ ]:


mean_intensity_all_countries = [] 
max_intensity_all_countries = [] 
for i in range (0,9861,10):
    print('i am at', i)
    edgar = ds.mhw_intensity[:,:,i].to_dataframe()
    edgar['newlons'] = nlon
    edgar = edgar.reset_index()
    edgar_gdf = gpd.GeoDataFrame(edgar, geometry=gpd.points_from_xy(edgar.newlons, edgar.lat))
    
    fiji_mean = country_mean_intensity(country = fiji)
    mean_intensity_all_countries.append(fiji_mean) 
    fiji_max = country_max_intensity(country = fiji)
    max_intensity_all_countries.append(fiji_max)
    
    ncd_mean = country_mean_intensity(country = ncd)
    mean_intensity_all_countries.append(ncd_mean) 
    ncd_max = country_max_intensity(country = ncd)
    max_intensity_all_countries.append(ncd_max)
    
    vanuatu_mean = country_mean_intensity(country = vanuatu)
    mean_intensity_all_countries.append(vanuatu_mean) 
    vanuatu_max = country_max_intensity(country = vanuatu)
    max_intensity_all_countries.append(vanuatu_max)
    
    wf_mean = country_mean_intensity(country = wf)
    mean_intensity_all_countries.append(wf_mean) 
    wf_max = country_max_intensity(country = wf)
    max_intensity_all_countries.append(wf_max)
    
    solo_mean = country_mean_intensity(country = solo)
    mean_intensity_all_countries.append(solo_mean) 
    solo_max = country_max_intensity(country = solo)
    max_intensity_all_countries.append(solo_max)
    
    tonga_mean = country_mean_intensity(country = tonga)
    mean_intensity_all_countries.append(tonga_mean) 
    tonga_max = country_max_intensity(country = tonga)
    max_intensity_all_countries.append(tonga_max)
    
    samoa_mean = country_mean_intensity(country = samoa)
    mean_intensity_all_countries.append(samoa_mean) 
    samoa_max = country_max_intensity(country = samoa)
    max_intensity_all_countries.append(samoa_max)
    
    tuvalu_mean = country_mean_intensity(country = tuvalu)
    mean_intensity_all_countries.append(tuvalu_mean) 
    tuvalu_max = country_max_intensity(country = tuvalu)
    max_intensity_all_countries.append(tuvalu_max)
    
    niue_mean = country_mean_intensity(country = niue)
    mean_intensity_all_countries.append(niue_mean) 
    niue_max = country_max_intensity(country = niue)
    max_intensity_all_countries.append(niue_max)
    
    cooks_mean = country_mean_intensity(country = cooks)
    mean_intensity_all_countries.append(cooks_mean) 
    cooks_max = country_max_intensity(country = cooks)
    max_intensity_all_countries.append(cooks_max)
    
    tokelau_mean = country_mean_intensity(country = tokelau)
    mean_intensity_all_countries.append(tokelau_mean) 
    tokelau_max = country_max_intensity(country = tokelau)
    max_intensity_all_countries.append(tokelau_max)
    
    amsam_mean = country_mean_intensity(country = amsam)
    mean_intensity_all_countries.append(amsam_mean) 
    amsam_max = country_max_intensity(country = amsam)
    max_intensity_all_countries.append(amsam_max)
    
 


# In[27]:


#len(range (0,9861,10))


# In[ ]:


max_intensity_all_countriesarr = np.array(max_intensity_all_countries)
print(max_intensity_all_countriesarr.shape)
max_intensity_countriesarr_rs = max_intensity_all_countriesarr.reshape((987,12))
dataframe = pd.DataFrame(max_intensity_countriesarr_rs) 
dataframe.to_csv("daily_max_intensity_countries_noaaoisst.csv")


# In[ ]:


mean_intensity_all_countriesarr = np.array(mean_intensity_all_countries)
print(mean_intensity_all_countriesarr.shape)
mean_intensity_countriesarr_rs = mean_intensity_all_countriesarr.reshape((987,12))
dataframe = pd.DataFrame(mean_intensity_countriesarr_rs) 
dataframe.to_csv("daily_mean_intensity_countries_noaaoisst.csv")

