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


# In[ ]:


df = pd.read_csv('percent_mhwstate_countries_noaaoisst.csv') 
percent_cover_all_countriesarr_rs = df.to_numpy()
percent_cover_all_countriesarr_rs.shape


# In[ ]:


ds = xr.open_dataset('mhw_intensity_noaa.nc')
ds


# In[ ]:


lat = np.array(ds.lat)
lon = np.array(ds.lon)
intensity3D = np.array(ds.mhw_intensity)


xx,yy = np.meshgrid(lon,lat)
print(xx.shape,yy.shape)


tm = np.array(ds.time)
t = tm.astype(int)
datessfull = [date.fromordinal(tt.astype(int)) for tt in tm]
print(datessfull)


# In[ ]:


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


# In[ ]:


df = pd.read_csv('lonlist_clip_eez.csv') 
dfarrno = df.to_numpy()
dfarrno.shape


# In[ ]:


nlon = dfarrno[:,1]


# In[ ]:


# to make script shorter
def country_percent_intensity_above_certain_value(country,threshold):    
    edgar_gdf.crs = country.crs
    country_data = gpd.clip(edgar_gdf, country)
    above_threshold = country_data[country_data['mhw_intensity'] >= threshold]
    percent = ((len(above_threshold))/(len(country_data)))*100
    return percent


# In[ ]:



#arr = percent_cover_all_countriesarr_rs[:,11]  #index 1 is for fiji
# Get the index/indices for a particular value
#value = 90
#indices = np.where(arr >= value)[0]

#ilist = []
#datesflist = []
#nanmaxlist = []
#percentabove1list = []
#percentabove2list = []
#percentabove3list = []
#percentabove4list = []

#for i in indices:
#    edgar = ds.mhw_intensity[:,:,i].to_dataframe()
#    edgar['newlons'] = nlon
#    edgar = edgar.reset_index()
#    edgar_gdf = gpd.GeoDataFrame(edgar, geometry=gpd.points_from_xy(edgar.newlons, edgar.lat))
#    edgar_gdf.crs = tokelau.crs
#    tokelau_data = gpd.clip(edgar_gdf, tokelau)
#    tokelauarr = np.array(tokelau_data)
#    nanmax_tokelau = (np.nanmax(tokelauarr[:,3]))
    
#    above_threshold_1 = tokelau_data[tokelau_data['mhw_intensity'] >= 1]
#    percent_above1 = ((len(above_threshold_1))/(len(tokelau_data)))*100
    
#    above_threshold_2 = tokelau_data[tokelau_data['mhw_intensity'] >= 2]
#    percent_above2 = ((len(above_threshold_2))/(len(tokelau_data)))*100
    
#    above_threshold_3 = tokelau_data[tokelau_data['mhw_intensity'] >= 3]
#    percent_above3 = ((len(above_threshold_3))/(len(tokelau_data)))*100
    
#    above_threshold_4 = tokelau_data[tokelau_data['mhw_intensity'] >= 4]
#    percent_above4 = ((len(above_threshold_4))/(len(tokelau_data)))*100
    
#    ilist.append(i)
#    datesflist.append(datessfull[i])
#    nanmaxlist.append(nanmax_tokelau)
#    percentabove1list.append(percent_above1)
#    percentabove2list.append(percent_above2)
#    percentabove3list.append(percent_above3)
#    percentabove4list.append(percent_above4)
    

#    print(i,datessfull[i],'max_intensity =',nanmax_tokelau,'percent_above1 =',percent_above1,
#         'percent_above2 =',percent_above2, 'percent_above3 =',percent_above3, 'percent_above4 =',percent_above4)


# In[ ]:


morethan_1deg_intensity_all_countries = []  
morethan_2deg_intensity_all_countries = [] 
morethan_3deg_intensity_all_countries = [] 
morethan_4deg_intensity_all_countries = [] 
for i in range (9861):
    edgar = ds.mhw_intensity[:,:,i].to_dataframe()
    edgar['newlons'] = nlon
    edgar = edgar.reset_index()
    edgar_gdf = gpd.GeoDataFrame(edgar, geometry=gpd.points_from_xy(edgar.newlons, edgar.lat))
    
    fiji_temp_1 = country_percent_intensity_above_certain_value(country = fiji,threshold = 1)
    morethan_1deg_intensity_all_countries.append(fiji_temp_1)
    fiji_temp_2 = country_percent_intensity_above_certain_value(country = fiji,threshold = 2)
    morethan_2deg_intensity_all_countries.append(fiji_temp_2)
    fiji_temp_3 = country_percent_intensity_above_certain_value(country = fiji,threshold = 3)
    morethan_3deg_intensity_all_countries.append(fiji_temp_3)
    fiji_temp_4 = country_percent_intensity_above_certain_value(country = fiji,threshold = 4)
    morethan_4deg_intensity_all_countries.append(fiji_temp_4)
    
    
    ncd_temp_1 = country_percent_intensity_above_certain_value(country = ncd,threshold = 1)
    morethan_1deg_intensity_all_countries.append(ncd_temp_1)
    ncd_temp_2 = country_percent_intensity_above_certain_value(country = ncd,threshold = 2)
    morethan_2deg_intensity_all_countries.append(ncd_temp_2)
    ncd_temp_3 = country_percent_intensity_above_certain_value(country = ncd,threshold = 3)
    morethan_3deg_intensity_all_countries.append(ncd_temp_3)
    ncd_temp_4 = country_percent_intensity_above_certain_value(country = ncd,threshold = 4)
    morethan_4deg_intensity_all_countries.append(ncd_temp_4)
    
    vanuatu_temp_1 = country_percent_intensity_above_certain_value(country = vanuatu,threshold = 1)
    morethan_1deg_intensity_all_countries.append(vanuatu_temp_1)
    vanuatu_temp_2 = country_percent_intensity_above_certain_value(country = vanuatu,threshold = 2)
    morethan_2deg_intensity_all_countries.append(vanuatu_temp_2)
    vanuatu_temp_3 = country_percent_intensity_above_certain_value(country = vanuatu,threshold = 3)
    morethan_3deg_intensity_all_countries.append(vanuatu_temp_3)
    vanuatu_temp_4 = country_percent_intensity_above_certain_value(country = vanuatu,threshold = 4)
    morethan_4deg_intensity_all_countries.append(vanuatu_temp_4)
    
    wf_temp_1 = country_percent_intensity_above_certain_value(country = wf,threshold = 1)
    morethan_1deg_intensity_all_countries.append(wf_temp_1)
    wf_temp_2 = country_percent_intensity_above_certain_value(country = wf,threshold = 2)
    morethan_2deg_intensity_all_countries.append(wf_temp_2)
    wf_temp_3 = country_percent_intensity_above_certain_value(country = wf,threshold = 3)
    morethan_3deg_intensity_all_countries.append(wf_temp_3)
    wf_temp_4 = country_percent_intensity_above_certain_value(country = wf,threshold = 4)
    morethan_4deg_intensity_all_countries.append(wf_temp_4)
    
    solo_temp_1 = country_percent_intensity_above_certain_value(country = solo,threshold = 1)
    morethan_1deg_intensity_all_countries.append(solo_temp_1)
    solo_temp_2 = country_percent_intensity_above_certain_value(country = solo,threshold = 2)
    morethan_2deg_intensity_all_countries.append(solo_temp_2)
    solo_temp_3 = country_percent_intensity_above_certain_value(country = solo,threshold = 3)
    morethan_3deg_intensity_all_countries.append(solo_temp_3)
    solo_temp_4 = country_percent_intensity_above_certain_value(country = solo,threshold = 4)
    morethan_4deg_intensity_all_countries.append(solo_temp_4)
    
    tonga_temp_1 = country_percent_intensity_above_certain_value(country = tonga,threshold = 1)
    morethan_1deg_intensity_all_countries.append(tonga_temp_1)
    tonga_temp_2 = country_percent_intensity_above_certain_value(country = tonga,threshold = 2)
    morethan_2deg_intensity_all_countries.append(tonga_temp_2)
    tonga_temp_3 = country_percent_intensity_above_certain_value(country = tonga,threshold = 3)
    morethan_3deg_intensity_all_countries.append(tonga_temp_3)
    tonga_temp_4 = country_percent_intensity_above_certain_value(country = tonga,threshold = 4)
    morethan_4deg_intensity_all_countries.append(tonga_temp_4)
    
    samoa_temp_1 = country_percent_intensity_above_certain_value(country = samoa,threshold = 1)
    morethan_1deg_intensity_all_countries.append(samoa_temp_1)
    samoa_temp_2 = country_percent_intensity_above_certain_value(country = samoa,threshold = 2)
    morethan_2deg_intensity_all_countries.append(samoa_temp_2)
    samoa_temp_3 = country_percent_intensity_above_certain_value(country = samoa,threshold = 3)
    morethan_3deg_intensity_all_countries.append(samoa_temp_3)
    samoa_temp_4 = country_percent_intensity_above_certain_value(country = samoa,threshold = 4)
    morethan_4deg_intensity_all_countries.append(samoa_temp_4)
    
    tuvalu_temp_1 = country_percent_intensity_above_certain_value(country = tuvalu,threshold = 1)
    morethan_1deg_intensity_all_countries.append(tuvalu_temp_1)
    tuvalu_temp_2 = country_percent_intensity_above_certain_value(country = tuvalu,threshold = 2)
    morethan_2deg_intensity_all_countries.append(tuvalu_temp_2)
    tuvalu_temp_3 = country_percent_intensity_above_certain_value(country = tuvalu,threshold = 3)
    morethan_3deg_intensity_all_countries.append(tuvalu_temp_3)
    tuvalu_temp_4 = country_percent_intensity_above_certain_value(country = tuvalu,threshold = 4)
    morethan_4deg_intensity_all_countries.append(tuvalu_temp_4)
    
    niue_temp_1 = country_percent_intensity_above_certain_value(country = niue,threshold = 1)
    morethan_1deg_intensity_all_countries.append(niue_temp_1)
    niue_temp_2 = country_percent_intensity_above_certain_value(country = niue,threshold = 2)
    morethan_2deg_intensity_all_countries.append(niue_temp_2)
    niue_temp_3 = country_percent_intensity_above_certain_value(country = niue,threshold = 3)
    morethan_3deg_intensity_all_countries.append(niue_temp_3)
    niue_temp_4 = country_percent_intensity_above_certain_value(country = niue,threshold = 4)
    morethan_4deg_intensity_all_countries.append(niue_temp_4)
    
    cooks_temp_1 = country_percent_intensity_above_certain_value(country = cooks,threshold = 1)
    morethan_1deg_intensity_all_countries.append(cooks_temp_1)
    cooks_temp_2 = country_percent_intensity_above_certain_value(country = cooks,threshold = 2)
    morethan_2deg_intensity_all_countries.append(cooks_temp_2)
    cooks_temp_3 = country_percent_intensity_above_certain_value(country = cooks,threshold = 3)
    morethan_3deg_intensity_all_countries.append(cooks_temp_3)
    cooks_temp_4 = country_percent_intensity_above_certain_value(country = cooks,threshold = 4)
    morethan_4deg_intensity_all_countries.append(cooks_temp_4)
    
    tokelau_temp_1 = country_percent_intensity_above_certain_value(country = tokelau,threshold = 1)
    morethan_1deg_intensity_all_countries.append(tokelau_temp_1)
    tokelau_temp_2 = country_percent_intensity_above_certain_value(country = tokelau,threshold = 2)
    morethan_2deg_intensity_all_countries.append(tokelau_temp_2)
    tokelau_temp_3 = country_percent_intensity_above_certain_value(country = tokelau,threshold = 3)
    morethan_3deg_intensity_all_countries.append(tokelau_temp_3)
    tokelau_temp_4 = country_percent_intensity_above_certain_value(country = tokelau,threshold = 4)
    morethan_4deg_intensity_all_countries.append(tokelau_temp_4)
    
    amsam_temp_1 = country_percent_intensity_above_certain_value(country = amsam,threshold = 1)
    morethan_1deg_intensity_all_countries.append(amsam_temp_1)
    amsam_temp_2 = country_percent_intensity_above_certain_value(country = amsam,threshold = 2)
    morethan_2deg_intensity_all_countries.append(amsam_temp_2)
    amsam_temp_3 = country_percent_intensity_above_certain_value(country = amsam,threshold = 3)
    morethan_3deg_intensity_all_countries.append(amsam_temp_3)
    amsam_temp_4 = country_percent_intensity_above_certain_value(country = amsam,threshold = 4)
    morethan_4deg_intensity_all_countries.append(amsam_temp_4)


# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:




