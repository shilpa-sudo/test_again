#!/usr/bin/env python
# coding: utf-8

# In[1]:


#%matplotlib notebook
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from netCDF4 import Dataset
from matplotlib.dates import num2date
import glob
from matplotlib import dates as mdates
import datetime

from datetime import date

import marineHeatWaves as mhw

from cartopy import config
import cartopy.crs as ccrs

import xarray as xr
import pandas as pd

import geopandas as gpd


# In[2]:


df = pd.read_csv('/home/shilpa/glory_mat_analysis/glorys_quarter_depth700m/mhws_fullset_depth700m.csv')


# In[3]:


df


# In[4]:


df['depth'] = 763.33


# In[5]:


df


# In[6]:


#df1 = pd.read_csv('/home/shilpa/glory_mat_analysis/glorys_quarter_depth0/mhws_fullset_depth0m.csv')


# In[7]:


#df1


# In[11]:


df.to_csv('/home/shilpa/glory_mat_analysis/glorys_quarter_depth700m/mhws_fullset_depth700m_withdepth.csv')








# In[12]:


#df2 = pd.read_csv('/home/shilpa/glory_mat_analysis/glorys_quarter_depth0/mhws_fullset_depth0m_withdepth.csv')


# In[13]:


#df2


# In[ ]:




