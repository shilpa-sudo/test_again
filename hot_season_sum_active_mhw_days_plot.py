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

var = xr.open_dataset('hot_season_polyareas_specific_sizes_range_noaa_93_22_sum_active_mhw_days.nc')

#sum active mhw days plot

sum_var_0_1 = np.array(var.hot_sum_mhw_days_pacific_0_1)
sum_var_1_10 = np.array(var.hot_sum_mhw_days_pacific_1_10)
sum_var_10_25 = np.array(var.hot_sum_mhw_days_pacific_10_25)
sum_var_25_50 = np.array(var.hot_sum_mhw_days_pacific_25_50)
sum_var_50_100 = np.array(var.hot_sum_mhw_days_pacific_50_100)
sum_var_100_200 = np.array(var.hot_sum_mhw_days_pacific_100_200)
sum_var_200_500 = np.array(var.hot_sum_mhw_days_pacific_200_500)
sum_var_500_700 = np.array(var.hot_sum_mhw_days_pacific_500_700)



print (np.nanmean(sum_var_0_1),np.nanmin(sum_var_0_1),np.nanmax(sum_var_0_1))
print (np.nanmean(sum_var_1_10),np.nanmin(sum_var_1_10),np.nanmax(sum_var_1_10))
print (np.nanmean(sum_var_10_25),np.nanmin(sum_var_10_25),np.nanmax(sum_var_10_25))
print (np.nanmean(sum_var_25_50),np.nanmin(sum_var_25_50),np.nanmax(sum_var_25_50))
print (np.nanmean(sum_var_50_100),np.nanmin(sum_var_50_100),np.nanmax(sum_var_50_100))
print (np.nanmean(sum_var_100_200),np.nanmin(sum_var_100_200),np.nanmax(sum_var_100_200))
print (np.nanmean(sum_var_200_500),np.nanmin(sum_var_200_500),np.nanmax(sum_var_200_500))
print (np.nanmean(sum_var_500_700),np.nanmin(sum_var_500_700),np.nanmax(sum_var_500_700))


min_mean=min([np.nanmin(sum_var_0_1),np.nanmin(sum_var_1_10),np.nanmin(sum_var_10_25),np.nanmin(sum_var_25_50),np.nanmin(sum_var_50_100),np.nanmin(sum_var_100_200),np.nanmin(sum_var_200_500),np.nanmin(sum_var_500_700)])
max_mean=max([np.nanmax(sum_var_0_1),np.nanmax(sum_var_1_10),np.nanmax(sum_var_10_25),np.nanmax(sum_var_25_50),np.nanmax(sum_var_50_100),np.nanmax(sum_var_100_200),np.nanmax(sum_var_200_500),np.nanmax(sum_var_500_700)])




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

fig, axs = plt.subplots(nrows=nrows,ncols=ncols,subplot_kw={'projection': proj_crs},figsize=(20,30))

uneven_levels = [0,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100,105,110,115,120,125,150,175,200,225,250,275,300,350,400,450,500,550,600,650,700]

easy_plot_unevencolorbars(lon=xx,lat=yy, data=sum_var_0_1,uneven_levels=uneven_levels,name='Hot season No. of MHW days 1-10 sqr degs',units=' [Days]',i=0,j=0)


easy_plot_unevencolorbars(lon=xx,lat=yy, data=sum_var_1_10,
                  uneven_levels=uneven_levels,name='Hot season No. of MHW days 1-10 sqr degs',units=' [Days]',i=0,j=1)


easy_plot_unevencolorbars(lon=xx,lat=yy, data=sum_var_10_25,
                  uneven_levels=uneven_levels,name='Hot season No. of MHW days 10-25 sqr degs',units=' [Days]',i=1,j=0)

easy_plot_unevencolorbars(lon=xx,lat=yy, data=sum_var_25_50,
                  uneven_levels=uneven_levels,name='Hot season No. of MHW days 25-50 sqr degs',units=' [Days]',i=1,j=1)


easy_plot_unevencolorbars(lon=xx,lat=yy, data=sum_var_50_100,
                  uneven_levels=uneven_levels,name='Hot season No. of MHW days 50-100 sqr degs',units=' [Days]',i=2,j=0)

easy_plot_unevencolorbars(lon=xx,lat=yy, data=sum_var_100_200,
                  uneven_levels=uneven_levels,name='Hot season No. of MHW days 100-200 sqr degs',units=' [Days]',i=2,j=1)


easy_plot_unevencolorbars(lon=xx,lat=yy, data=sum_var_200_500,
                  uneven_levels=uneven_levels,name='Hot season No. of MHW days 200-500 sqr degs',units=' [Days]',i=3,j=0)

easy_plot_unevencolorbars(lon=xx,lat=yy, data=sum_var_500_700,
                  uneven_levels=uneven_levels,name='Hot season No. of MHW days 500-700 sqr degs',units=' [Days]',i=3,j=1)




axs[0, 0].text(-0.2, 1.05, 'a', transform=axs[0, 0].transAxes, fontsize=20, fontweight='bold')
axs[0, 1].text(-0.2, 1.05, 'b', transform=axs[0, 1].transAxes, fontsize=20, fontweight='bold')
axs[1, 0].text(-0.2, 1.05, 'C', transform=axs[1, 0].transAxes, fontsize=20, fontweight='bold')
axs[1, 1].text(-0.2, 1.05, 'd', transform=axs[1, 1].transAxes, fontsize=20, fontweight='bold')
axs[2, 0].text(-0.2, 1.05, 'e', transform=axs[2, 0].transAxes, fontsize=20, fontweight='bold')
axs[2, 1].text(-0.2, 1.05, 'f', transform=axs[2, 1].transAxes, fontsize=20, fontweight='bold')
axs[3, 0].text(-0.2, 1.05, 'g', transform=axs[3, 0].transAxes, fontsize=20, fontweight='bold')
axs[3, 1].text(-0.2, 1.05, 'h', transform=axs[3, 1].transAxes, fontsize=20, fontweight='bold')


plt.tight_layout()
plt.savefig('hot_sum_mhw_active_cells_noaa_polygon_size_range.png',bbox_inches='tight')
