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
import scipy.stats as ss
import seaborn as sb
import geopandas as gpd
from shapely.geometry import Polygon
import scipy as sp
from scipy import linalg
from scipy import stats
import scipy.ndimage as ndimage
from datetime import date


# In[2]:


eez = gpd.read_file("/home/shilpa/Downloads/World_EEZ_v11_20191118/eez_v11.shp") # /home/shilpa/Downloads/World_EEZ_v11_20191118/eez_v11.shp, /home/shilpa/Downloads/World_EEZ_v11_20191118/eez_boundaries_v11.shp
#print(eez)


# In[3]:


# get eezs
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


# In[4]:


df = pd.read_csv('/home/shilpa/glory_mat_analysis/lonlist_clip_eez.csv') 
dfarrno = df.to_numpy()
dfarrno.shape


# In[5]:


nlon = dfarrno[:,1]


# In[6]:


spatial = "full_spatial_noaaoistt_aoi_lastdayincluded.nc"


# In[7]:


intensity = "full_intensity_spatial_noaaoisst_aoi_lastdayincluded.nc"


# In[8]:


ds_spatial = xr.open_dataset(spatial)
ds_intensity = xr.open_dataset(intensity)


# In[9]:


ds_intensity


# In[10]:


def country_percent_area(country):    
    edgar_gdf_spatial.crs = country.crs
    country_data = gpd.clip(edgar_gdf_spatial, country)
    percent_area_cov_country = (np.nanmean(country_data.spatial_extent))*100
    return percent_area_cov_country


# In[11]:


def country_max_intensity(country):    
    edgar_gdf_intensity.crs = country.crs
    country_data = gpd.clip(edgar_gdf_intensity, country)
    max_intensity_country = (np.nanmax(country_data.spatial_intensity))
    return max_intensity_country


# In[12]:


fiji_spatial = []
ncd_spatial = []
vanuatu_spatial = []
wf_spatial = []
solo_spatial = []
tonga_spatial = []
samoa_spatial =  []
tuvalu_spatial = []
niue_spatial = []
cooks_spatial = []
tokelau_spatial = []
amsam_spatial = []


fiji_intensity = []
ncd_intensity = []
vanuatu_intensity = []
wf_intensity = []
solo_intensity = []
tonga_intensity = []
samoa_intensity =  []
tuvalu_intensity = []
niue_intensity = []
cooks_intensity = []
tokelau_intensity = []
amsam_intensity = []


for i in range (9861):
    print('i =',i)
    edgar_spatial = ds_spatial.spatial_extent[:,:,i].to_dataframe()
    edgar_spatial['newlons'] = nlon
    edgar_spatial = edgar_spatial.reset_index()
    edgar_gdf_spatial = gpd.GeoDataFrame(edgar_spatial, geometry=gpd.points_from_xy(edgar_spatial.newlons, edgar_spatial.lats))
    
    edgar_intensity = ds_intensity.spatial_intensity[:,:,i].to_dataframe()
    edgar_intensity['newlons'] = nlon
    edgar_intensity = edgar_intensity.reset_index()
    edgar_gdf_intensity = gpd.GeoDataFrame(edgar_intensity, geometry=gpd.points_from_xy(edgar_intensity.newlons, edgar_intensity.lats))
    
    fiji_area = country_percent_area(fiji)
    fiji_spatial.append(fiji_area)
    fiji_max_intensity = country_max_intensity(fiji)
    fiji_intensity.append(fiji_max_intensity)
    
    ncd_area = country_percent_area(ncd)
    ncd_spatial.append(ncd_area)
    ncd_max_intensity = country_max_intensity(ncd)
    ncd_intensity.append(ncd_max_intensity)
    
    vanuatu_area = country_percent_area(vanuatu)
    vanuatu_spatial.append(vanuatu_area)
    vanuatu_max_intensity = country_max_intensity(vanuatu)
    vanuatu_intensity.append(vanuatu_max_intensity)
    
    wf_area = country_percent_area(wf)
    wf_spatial.append(wf_area)
    wf_max_intensity = country_max_intensity(wf)
    wf_intensity.append(wf_max_intensity)
    
    solo_area = country_percent_area(solo)
    solo_spatial.append(solo_area)
    solo_max_intensity = country_max_intensity(solo)
    solo_intensity.append(solo_max_intensity)
    
    tonga_area = country_percent_area(tonga)
    tonga_spatial.append(tonga_area)
    tonga_max_intensity = country_max_intensity(tonga)
    tonga_intensity.append(tonga_max_intensity)

    
    samoa_area = country_percent_area(samoa)
    samoa_spatial.append(samoa_area)
    samoa_max_intensity = country_max_intensity(samoa)
    samoa_intensity.append(samoa_max_intensity)

    tuvalu_area = country_percent_area(tuvalu)
    tuvalu_spatial.append(tuvalu_area)
    tuvalu_max_intensity = country_max_intensity(tuvalu)
    tuvalu_intensity.append(tuvalu_max_intensity)

    niue_area = country_percent_area(niue)
    niue_spatial.append(niue_area)
    niue_max_intensity = country_max_intensity(niue)
    niue_intensity.append(niue_max_intensity)

    cooks_area = country_percent_area(cooks)
    cooks_spatial.append(cooks_area)
    cooks_max_intensity = country_max_intensity(cooks)
    cooks_intensity.append(cooks_max_intensity)
 
    tokelau_area = country_percent_area(tokelau)
    tokelau_spatial.append(tokelau_area)
    tokelau_max_intensity = country_max_intensity(tokelau)
    tokelau_intensity.append(tokelau_max_intensity)
    
    amsam_area = country_percent_area(amsam)
    amsam_spatial.append(amsam_area)
    amsam_max_intensity = country_max_intensity(amsam)
    amsam_intensity.append(amsam_max_intensity)
    

    


# In[ ]:


t = pd.date_range(start='1/1/1993', end='31/12/2019')


# In[ ]:


d = {'time': t, 
    'Fiji_spatial_extent':fiji_spatial,'Fiji_intensity': fiji_intensity,
    'New_Caledonia_spatial_extent':ncd_spatial,'New_Caledonia_intensity': ncd_intensity,
    'Vanuatu_spatial_extent':vanuatu_spatial,'Vanuatu_intensity': vanuatu_intensity,
    'Wallis_and_Futuna_spatial_extent':wf_spatial,'Wallis_and_Futuna_intensity': wf_intensity,
    'Solomon_Islands_spatial_extent':solo_spatial,'Solomon_Islands_intensity': solo_intensity,
    'Tonga_spatial_extent':tonga_spatial,'Tonga_intensity': tonga_intensity,
    'Samoa_spatial_extent':samoa_spatial,'Samoa_intensity': samoa_intensity,
    'Tuvalu_spatial_extent':tuvalu_spatial,'Tuvalu_intensity': tuvalu_intensity,
    'Niue_spatial_extent':niue_spatial,'Niue_intensity': niue_intensity,
    'Cook_Islands_spatial_extent':cooks_spatial,'Cook_Islands_intensity': cooks_intensity,
    'Tokelau_spatial_extent':tokelau_spatial,'Tokelau_intensity': tokelau_intensity,
    'American_Samoa_spatial_extent':amsam_spatial,'American_Samoa_intensity': amsam_intensity}

dfnew = pd.DataFrame(data=d)

dfnew


# In[ ]:


dfnew.to_csv("percent_mhwstate_intensities_countries_noaaoisst.csv")


# In[ ]:




