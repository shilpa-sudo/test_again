#standard imports
from functools import partial
import xarray as xr
import pandas as pd
import numpy as np
from datetime import date
import statistics
import marineHeatWaves as mhw
import scipy as sp
from scipy import linalg
from scipy import stats
import scipy.ndimage as ndimage
from datetime import date
from cartopy import config
import cartopy.crs as ccrs
import geopandas as gpd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

var = xr.open_dataset('hot_season_no_of_events_polyareas_specific_sizes_range_noaa_93_22.nc')

#No_mhw_events

var_0_1 = np.array(var.hot_season_no_mhw_events_pacific_0_1)
var_1_10 = np.array(var.hot_season_no_mhw_events_pacific_1_10)
var_10_25 = np.array(var.hot_season_no_mhw_events_pacific_10_25)
var_25_50 = np.array(var.hot_season_no_mhw_events_pacific_25_50)
var_50_100 = np.array(var.hot_season_no_mhw_events_pacific_50_100)
var_100_200 = np.array(var.hot_season_no_mhw_events_pacific_100_200)
var_200_500 = np.array(var.hot_season_no_mhw_events_pacific_200_500)
var_500_700 = np.array(var.hot_season_no_mhw_events_pacific_500_700)



print (np.nanmean(var_0_1),np.nanmin(var_0_1),np.nanmax(var_0_1))
print (np.nanmean(var_1_10),np.nanmin(var_1_10),np.nanmax(var_1_10))
print (np.nanmean(var_10_25),np.nanmin(var_10_25),np.nanmax(var_10_25))
print (np.nanmean(var_25_50),np.nanmin(var_25_50),np.nanmax(var_25_50))
print (np.nanmean(var_50_100),np.nanmin(var_50_100),np.nanmax(var_50_100))
print (np.nanmean(var_100_200),np.nanmin(var_100_200),np.nanmax(var_100_200))
print (np.nanmean(var_200_500),np.nanmin(var_200_500),np.nanmax(var_200_500))
print (np.nanmean(var_500_700),np.nanmin(var_500_700),np.nanmax(var_500_700))


min_events=min([np.nanmin(var_0_1),np.nanmin(var_1_10),np.nanmin(var_10_25),np.nanmin(var_25_50),np.nanmin(var_50_100),np.nanmin(var_100_200),np.nanmin(var_200_500),np.nanmin(var_500_700)])
max_events=max([np.nanmax(var_0_1),np.nanmax(var_1_10),np.nanmax(var_10_25),np.nanmax(var_25_50),np.nanmax(var_50_100),np.nanmax(var_100_200),np.nanmax(var_200_500),np.nanmax(var_500_700)])

print(min_events,max_events)


#projection
data_crs = ccrs.PlateCarree(central_longitude=0)

# But you want the map to include the dateline.
proj_crs = ccrs.PlateCarree(central_longitude=180)

# eez and country boundaries
eez = gpd.read_file("/home/shilpa/Downloads/World_EEZ_v11_20191118/eez_v11.shp")

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


# lat ,lon
lat = var.lat
lon = var.lon

xx,yy = np.meshgrid(lon,lat)

# function to make plots
def easy_plot_unevencolorbars(lon,lat,data,i,j,name,units,uneven_levels):
    uneven_levels = uneven_levels
    cmap_rb = plt.get_cmap('turbo')#'RdBu_r')
    colors = cmap_rb(np.linspace(0, 1, (len(uneven_levels) - 1)))
    cmap, norm = mcolors.from_levels_and_colors(uneven_levels, colors)
    
    im = axs[i,j].contourf(lon,lat,data, 
               levels=uneven_levels,
               cmap=cmap, norm=norm,
               transform=data_crs,
               transform_first=True)
    #axs[i,j].coastlines()
    #eez[:].plot(ax=axs[i,j],color='none',edgecolor='black',transform=data_crs)
    axs[i,j].coastlines()
    fiji.boundary.plot(ax=axs[i,j], color="black",alpha = 0.4,transform=data_crs)
    ncd.boundary.plot(ax=axs[i,j], color="black",alpha = 0.4,transform=data_crs)
    vanuatu.boundary.plot(ax=axs[i,j], color="black",alpha = 0.4,transform=data_crs)
    wf.boundary.plot(ax=axs[i,j], color="black",alpha = 0.4,transform=data_crs)
    solo.boundary.plot(ax=axs[i,j], color="black",alpha = 0.4,transform=data_crs)
    tonga.boundary.plot(ax=axs[i,j], color="black",alpha = 0.4,transform=data_crs)
    samoa.boundary.plot(ax=axs[i,j], color="black",alpha = 0.4,transform=data_crs)
    tuvalu.boundary.plot(ax=axs[i,j], color="black",alpha = 0.4,transform=data_crs)
    niue.boundary.plot(ax=axs[i,j], color="black",alpha = 0.4,transform=data_crs)
    cooks.boundary.plot(ax=axs[i,j], color="black",alpha = 0.4,transform=data_crs)
    tokelau.boundary.plot(ax=axs[i,j], color="black",alpha = 0.4,transform=data_crs)
    amsam.boundary.plot(ax=axs[i,j], color="black",alpha = 0.4,transform=data_crs)

    
    cbar = plt.colorbar(im, ax=axs[i,j],orientation="horizontal", shrink=0.75)
    cbar.set_label(units,fontsize=14)
    gl = axs[i,j].gridlines(draw_labels=True)
    gl.top_labels = False
    gl.right_labels = False
    axs[i,j].set_xlim(-35,29) #-35
    axs[i,j].set_ylim(-34.875,-2.375)
    axs[i,j].set_title(name, fontsize=18) 



nrows = 4
ncols = 2

varname = 'Hot season No. of MHW events'
units_ = '[events]'

fig, axs = plt.subplots(nrows=nrows,ncols=ncols,subplot_kw={'projection': proj_crs},figsize=(20,30))
#uneven_levels = np.arange(min_events,max_events,1)
uneven_levels = [1,2,3,4,5,6,7,8,9,10,12,14,16,17,18,19,20,21,22,24,26,28,30,32,34,36,38,40,45,50,60,70,89]

easy_plot_unevencolorbars(lon=xx,lat=yy, data=var_0_1,
                  uneven_levels=uneven_levels,name=f'{varname} 0-1 sqr degs',units=f'units_',i=0,j=0)

easy_plot_unevencolorbars(lon=xx,lat=yy, data=var_1_10,
                  uneven_levels=uneven_levels,name=f'{varname} 1-10 sqr degs',units=f'units_',i=0,j=1)


easy_plot_unevencolorbars(lon=xx,lat=yy, data=var_10_25,
                  uneven_levels=uneven_levels,name=f'{varname} 10-25 sqr degs',units=f'units_',i=1,j=0)

easy_plot_unevencolorbars(lon=xx,lat=yy, data=var_25_50,
                  uneven_levels=uneven_levels,name=f'{varname} 25-50 sqr degs',units=f'units_',i=1,j=1)


easy_plot_unevencolorbars(lon=xx,lat=yy, data=var_50_100,
                  uneven_levels=uneven_levels,name=f'{varname} 50-100 sqr degs',units=f'units_',i=2,j=0)

easy_plot_unevencolorbars(lon=xx,lat=yy, data=var_100_200,
                  uneven_levels=uneven_levels,name=f'{varname} 100-200 sqr degs',units=f'units_',i=2,j=1)


easy_plot_unevencolorbars(lon=xx,lat=yy, data=var_200_500,
                  uneven_levels=uneven_levels,name=f'{varname} 200-500 sqr degs',units=f'units_',i=3,j=0)

easy_plot_unevencolorbars(lon=xx,lat=yy, data=var_500_700,
                  uneven_levels=uneven_levels,name=f'{varname} 500-700 sqr degs',units=f'units_',i=3,j=1)




axs[0, 0].text(-0.2, 1.05, 'a', transform=axs[0, 0].transAxes, fontsize=20, fontweight='bold')
axs[0, 1].text(-0.2, 1.05, 'b', transform=axs[0, 1].transAxes, fontsize=20, fontweight='bold')
axs[1, 0].text(-0.2, 1.05, 'C', transform=axs[1, 0].transAxes, fontsize=20, fontweight='bold')
axs[1, 1].text(-0.2, 1.05, 'd', transform=axs[1, 1].transAxes, fontsize=20, fontweight='bold')
axs[2, 0].text(-0.2, 1.05, 'e', transform=axs[2, 0].transAxes, fontsize=20, fontweight='bold')
axs[2, 1].text(-0.2, 1.05, 'f', transform=axs[2, 1].transAxes, fontsize=20, fontweight='bold')
axs[3, 0].text(-0.2, 1.05, 'g', transform=axs[3, 0].transAxes, fontsize=20, fontweight='bold')
axs[3, 1].text(-0.2, 1.05, 'h', transform=axs[3, 1].transAxes, fontsize=20, fontweight='bold')


plt.tight_layout()
plt.savefig('hot_total_num_events_noaa_polygon_size_range.png',bbox_inches='tight')
