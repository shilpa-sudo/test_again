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
import marineHeatWaves as mhw
from cartopy import config
import cartopy.crs as ccrs
import xarray as xr
import pandas as pd
import geopandas as gpd
import matplotlib.colors as mcolors

fnames = ['/home/shilpa/glory_mat_analysis/glorys_quarter_depth0/mhws_fullset_depth0m_withdepth.csv',
 '/home/shilpa/glory_mat_analysis/glorys_quarter_depth1000m/mhws_fullset_depth1000m_withdepth.csv',
 '/home/shilpa/glory_mat_analysis/glorys_quarter_depth100m/mhws_fullset_depth100m_withdepth.csv',
 '/home/shilpa/glory_mat_analysis/glorys_quarter_depth10m/mhws_fullset_depth10m_withdepth.csv',
 '/home/shilpa/glory_mat_analysis/glorys_quarter_depth11m/mhws_fullset_depth11m_withdepth.csv',
 '/home/shilpa/glory_mat_analysis/glorys_quarter_depth1200m/mhws_fullset_depth1200m_withdepth.csv',
 '/home/shilpa/glory_mat_analysis/glorys_quarter_depth130m/mhws_fullset_depth130m_withdepth.csv',
 '/home/shilpa/glory_mat_analysis/glorys_quarter_depth13m/mhws_fullset_depth13m_withdepth.csv',
 '/home/shilpa/glory_mat_analysis/glorys_quarter_depth1500m/mhws_fullset_depth1500m_withdepth.csv',
 '/home/shilpa/glory_mat_analysis/glorys_quarter_depth155m/mhws_fullset_depth155m_withdepth.csv',
 '/home/shilpa/glory_mat_analysis/glorys_quarter_depth15m/mhws_fullset_depth15m_withdepth.csv',
 '/home/shilpa/glory_mat_analysis/glorys_quarter_depth186m/mhws_fullset_depth186m_withdepth.csv',
 '/home/shilpa/glory_mat_analysis/glorys_quarter_depth18m/mhws_fullset_depth18m_withdepth.csv',
 '/home/shilpa/glory_mat_analysis/glorys_quarter_depth1m/mhws_fullset_depth1m_withdepth.csv',
 '/home/shilpa/glory_mat_analysis/glorys_quarter_depth21m/mhws_fullset_depth21m_withdepth.csv',
 '/home/shilpa/glory_mat_analysis/glorys_quarter_depth222m/mhws_fullset_depth222m_withdepth.csv',
 '/home/shilpa/glory_mat_analysis/glorys_quarter_depth25m/mhws_fullset_depth25m_withdepth.csv',
 '/home/shilpa/glory_mat_analysis/glorys_quarter_depth266m/mhws_fullset_depth266m_withdepth.csv',
 '/home/shilpa/glory_mat_analysis/glorys_quarter_depth29m/mhws_fullset_depth29m_withdepth.csv',
 '/home/shilpa/glory_mat_analysis/glorys_quarter_depth2m/mhws_fullset_depth2m_withdepth.csv',
 '/home/shilpa/glory_mat_analysis/glorys_quarter_depth318m/mhws_fullset_depth318m_withdepth.csv',
 '/home/shilpa/glory_mat_analysis/glorys_quarter_depth34m/mhws_fullset_depth34m_withdepth.csv',
 '/home/shilpa/glory_mat_analysis/glorys_quarter_depth380m/mhws_fullset_depth380m_withdepth.csv',
 '/home/shilpa/glory_mat_analysis/glorys_quarter_depth3m/mhws_fullset_depth3m_withdepth.csv',
 '/home/shilpa/glory_mat_analysis/glorys_quarter_depth40m/mhws_fullset_depth40m_withdepth.csv',
 '/home/shilpa/glory_mat_analysis/glorys_quarter_depth453m/mhws_fullset_depth453m_withdepth.csv',
 '/home/shilpa/glory_mat_analysis/glorys_quarter_depth47m/mhws_fullset_depth47m_withdepth.csv',
 '/home/shilpa/glory_mat_analysis/glorys_quarter_depth500m/mhws_fullset_depth500m_withdepth.csv',
 '/home/shilpa/glory_mat_analysis/glorys_quarter_depth55m/mhws_fullset_depth55m_withdepth.csv',
 '/home/shilpa/glory_mat_analysis/glorys_quarter_depth5m/mhws_fullset_depth5m_withdepth.csv',
 '/home/shilpa/glory_mat_analysis/glorys_quarter_depth600m/mhws_fullset_depth600m_withdepth.csv',
 '/home/shilpa/glory_mat_analysis/glorys_quarter_depth65m/mhws_fullset_depth65m_withdepth.csv',
 '/home/shilpa/glory_mat_analysis/glorys_quarter_depth6m/mhws_fullset_depth6m_withdepth.csv',
 '/home/shilpa/glory_mat_analysis/glorys_quarter_depth700m/mhws_fullset_depth700m_withdepth.csv',
 '/home/shilpa/glory_mat_analysis/glorys_quarter_depth77m/mhws_fullset_depth77m_withdepth.csv',
 '/home/shilpa/glory_mat_analysis/glorys_quarter_depth7m/mhws_fullset_depth7m_withdepth.csv',
 '/home/shilpa/glory_mat_analysis/glorys_quarter_depth900m/mhws_fullset_depth900m_withdepth.csv',
 '/home/shilpa/glory_mat_analysis/glorys_quarter_depth92m/mhws_fullset_depth92m_withdepth.csv']


df_concat2 = pd.concat([pd.read_csv(f) for f in fnames], ignore_index = True)

df_concat2.to_csv('/home/shilpa/glory_mat_analysis/vertical_anal_gl_quarter/mhws_fullset_vertical_gl_quarterdeg.csv')



