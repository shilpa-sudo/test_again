#!/usr/bin/env python
# coding: utf-8

# In[1]:


# standard imports
#%matplotlib notebook
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from cartopy import config
import cartopy.crs as ccrs
import xarray as xr
import pandas as pd
import geopandas as gpd
import seaborn as sns


# In[2]:


df = pd.read_csv('/home/shilpa/glory_mat_analysis/glorys_quarter_depth0/mhws_fullset_depth0m_withdepth.csv')


# In[3]:


# clip for new caledonia eez
eez = gpd.read_file("/home/shilpa/Downloads/World_EEZ_v11_20191118/eez_v11.shp") # get the EEZs


# In[ ]:


edgar_gdf_mhw_events = gpd.GeoDataFrame(df, geometry=gpd.points_from_xy(df.longitude, df.latitude))


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


country = ncd
country_name = "New Caledonia"
edgar_gdf_mhw_events.crs = country.crs
country_data = gpd.clip(edgar_gdf_mhw_events, country)
country_data['country'] = country_name
country_data.to_csv(f"{country_name}_mhw_events_0m_glorys.csv")

