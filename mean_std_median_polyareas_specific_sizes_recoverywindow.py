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

var = xr.open_dataset('polyareas_specific_sizes_range_noaa_93_22_mean_std_median_recoverywindowdayspacific.nc')

mean_var_0_1 = np.array(var.mean_recovery_window_days_pacific_0_1)
mean_var_1_10 = np.array(var.mean_recovery_window_days_pacific_1_10)
mean_var_10_25 = np.array(var.mean_recovery_window_days_pacific_10_25)
mean_var_25_50 = np.array(var.mean_recovery_window_days_pacific_25_50)
mean_var_50_100 = np.array(var.mean_recovery_window_days_pacific_50_100)
mean_var_100_200 = np.array(var.mean_recovery_window_days_pacific_100_200)
mean_var_200_500 = np.array(var.mean_recovery_window_days_pacific_200_500)
mean_var_500_700 = np.array(var.mean_recovery_window_days_pacific_500_700)

std_var_0_1 = np.array(var.std_recovery_window_days_pacific_0_1)
std_var_1_10 = np.array(var.std_recovery_window_days_pacific_1_10)
std_var_10_25 = np.array(var.std_recovery_window_days_pacific_10_25)
std_var_25_50 = np.array(var.std_recovery_window_days_pacific_25_50)
std_var_50_100 = np.array(var.std_recovery_window_days_pacific_50_100)
std_var_100_200 = np.array(var.std_recovery_window_days_pacific_100_200)
std_var_200_500 = np.array(var.std_recovery_window_days_pacific_200_500)
std_var_500_700 = np.array(var.std_recovery_window_days_pacific_500_700)

median_var_0_1 = np.array(var.median_recovery_window_days_pacific_0_1)
median_var_1_10 = np.array(var.median_recovery_window_days_pacific_1_10)
median_var_10_25 = np.array(var.median_recovery_window_days_pacific_10_25)
median_var_25_50 = np.array(var.median_recovery_window_days_pacific_25_50)
median_var_50_100 = np.array(var.median_recovery_window_days_pacific_50_100)
median_var_100_200 = np.array(var.median_recovery_window_days_pacific_100_200)
median_var_200_500 = np.array(var.median_recovery_window_days_pacific_200_500)
median_var_500_700 = np.array(var.median_recovery_window_days_pacific_500_700)


print (np.nanmean(mean_var_0_1),np.nanmin(mean_var_0_1),np.nanmax(mean_var_0_1))
print (np.nanmean(mean_var_1_10),np.nanmin(mean_var_1_10),np.nanmax(mean_var_1_10))
print (np.nanmean(mean_var_10_25),np.nanmin(mean_var_10_25),np.nanmax(mean_var_10_25))
print (np.nanmean(mean_var_25_50),np.nanmin(mean_var_25_50),np.nanmax(mean_var_25_50))
print (np.nanmean(mean_var_50_100),np.nanmin(mean_var_50_100),np.nanmax(mean_var_50_100))
print (np.nanmean(mean_var_100_200),np.nanmin(mean_var_100_200),np.nanmax(mean_var_100_200))
print (np.nanmean(mean_var_200_500),np.nanmin(mean_var_200_500),np.nanmax(mean_var_200_500))
print (np.nanmean(mean_var_500_700),np.nanmin(mean_var_500_700),np.nanmax(mean_var_500_700))

print (np.nanmean(std_var_0_1),np.nanmin(std_var_0_1),np.nanmax(std_var_0_1))
print (np.nanmean(std_var_1_10),np.nanmin(std_var_1_10),np.nanmax(std_var_1_10))
print (np.nanmean(std_var_10_25),np.nanmin(std_var_10_25),np.nanmax(std_var_10_25))
print (np.nanmean(std_var_25_50),np.nanmin(std_var_25_50),np.nanmax(std_var_25_50))
print (np.nanmean(std_var_50_100),np.nanmin(std_var_50_100),np.nanmax(std_var_50_100))
print (np.nanmean(std_var_100_200),np.nanmin(std_var_100_200),np.nanmax(std_var_100_200))
print (np.nanmean(std_var_200_500),np.nanmin(std_var_200_500),np.nanmax(std_var_200_500))
print (np.nanmean(std_var_500_700),np.nanmin(std_var_500_700),np.nanmax(std_var_500_700))

print (np.nanmean(median_var_0_1),np.nanmin(median_var_0_1),np.nanmax(median_var_0_1))
print (np.nanmean(median_var_1_10),np.nanmin(median_var_1_10),np.nanmax(median_var_1_10))
print (np.nanmean(median_var_10_25),np.nanmin(median_var_10_25),np.nanmax(median_var_10_25))
print (np.nanmean(median_var_25_50),np.nanmin(median_var_25_50),np.nanmax(median_var_25_50))
print (np.nanmean(median_var_50_100),np.nanmin(median_var_50_100),np.nanmax(median_var_50_100))
print (np.nanmean(median_var_100_200),np.nanmin(median_var_100_200),np.nanmax(median_var_100_200))
print (np.nanmean(median_var_200_500),np.nanmin(median_var_200_500),np.nanmax(median_var_200_500))
print (np.nanmean(median_var_500_700),np.nanmin(median_var_500_700),np.nanmax(median_var_500_700))


min_mean=min([np.nanmin(mean_var_0_1),np.nanmin(mean_var_1_10),np.nanmin(mean_var_10_25),np.nanmin(mean_var_25_50),np.nanmin(mean_var_50_100),np.nanmin(mean_var_100_200),np.nanmin(mean_var_200_500),np.nanmin(mean_var_500_700)])
max_mean=max([np.nanmax(mean_var_0_1),np.nanmax(mean_var_1_10),np.nanmax(mean_var_10_25),np.nanmax(mean_var_25_50),np.nanmax(mean_var_50_100),np.nanmax(mean_var_100_200),np.nanmax(mean_var_100_200),np.nanmax(mean_var_200_500),np.nanmax(mean_var_500_700)])

min_std=min([np.nanmin(std_var_0_1),np.nanmin(std_var_1_10),np.nanmin(std_var_10_25),np.nanmin(std_var_25_50),np.nanmin(std_var_50_100),np.nanmin(std_var_100_200),np.nanmin(std_var_200_500),np.nanmin(std_var_500_700)])
max_std=max([np.nanmax(std_var_0_1),np.nanmax(std_var_1_10),np.nanmax(std_var_10_25),np.nanmax(std_var_25_50),np.nanmax(std_var_50_100),np.nanmax(std_var_100_200),np.nanmax(std_var_100_200),np.nanmax(std_var_200_500),np.nanmax(std_var_500_700)])

min_median=min([np.nanmin(median_var_0_1),np.nanmin(median_var_1_10),np.nanmin(median_var_10_25),np.nanmin(median_var_25_50),np.nanmin(median_var_50_100),np.nanmin(median_var_100_200),np.nanmin(median_var_200_500),np.nanmin(median_var_500_700)])
max_median=max([np.nanmax(median_var_0_1),np.nanmax(median_var_1_10),np.nanmax(median_var_10_25),np.nanmax(median_var_25_50),np.nanmax(median_var_50_100),np.nanmax(median_var_100_200),np.nanmax(median_var_100_200),np.nanmax(median_var_200_500),np.nanmax(median_var_500_700)])



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



nrows = 8
ncols = 3
variable = 'Recovery window'
units = '[Days]'

fig, axs = plt.subplots(nrows=nrows,ncols=ncols,subplot_kw={'projection': proj_crs},figsize=(20,30))

uneven_levels = [4,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100,125,150,175,200,210,215,220,225,230,235,240,245,250,255,260,265,270,280,290,300,310,320,340,350,360,370,380,390,400,450,550,650,750,850,950,1050,1100,1200,1300,1400,1500,1600,1700]
uneven_levels2 = [min_std,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100,120,140,160,180,200,300,400,500,600,700,800,900,1000,1100,1200,1300,1400,1500,1600,max_std]

size = '0-1 sqr.degs.'
easy_plot_unevencolorbars(lon=xx,lat=yy, data=mean_var_0_1,
                  uneven_levels=uneven_levels,name= f'Mean {variable} {size} ',units=f'{units}',i=0,j=0)

#uneven_levels2 = np.arange(min_std,max_std,0.1)
size = '0-1 sqr. degs.'
easy_plot_unevencolorbars(lon=xx,lat=yy, data=std_var_0_1,
                  uneven_levels=uneven_levels2,name= f'Std {size} ',units='Standard Deviation',i=0,j=1)

#uneven_levels = np.arange(min_median, max_median, 1)
size = '0-1 sqr. degs.'
easy_plot_unevencolorbars(lon=xx,lat=yy, data=median_var_0_1,
                  uneven_levels=uneven_levels,name=f'Median {size} ',units=f'{units}',i=0,j=2)

#uneven_levels = np.arange(min_mean, max_mean, 1)
size = '1-10 sqr. degs.'
easy_plot_unevencolorbars(lon=xx,lat=yy, data=mean_var_1_10,
                  uneven_levels=uneven_levels,name= f'Mean {variable} {size} ',units=f'{units}',i=1,j=0)

#uneven_levels = np.arange(min_std,max_std,0.1)
size = '1-10 sqr. degs.'
easy_plot_unevencolorbars(lon=xx,lat=yy, data=std_var_1_10,
                  uneven_levels=uneven_levels2,name= f'Std {size} ',units='Standard Deviation',i=1,j=1)

#uneven_levels = np.arange(min_median, max_median, 1)
size = '1-10 sqr. degs.'
easy_plot_unevencolorbars(lon=xx,lat=yy, data=median_var_1_10,
                  uneven_levels=uneven_levels,name=f'Median {size} ',units=f'{units}',i=1,j=2)

#uneven_levels = np.arange(min_mean, max_mean, 1)
size = '10-25 sqr. degs.'
easy_plot_unevencolorbars(lon=xx,lat=yy, data=mean_var_10_25,
                  uneven_levels=uneven_levels,name= f'Mean {variable} {size} ',units=f'{units}',i=2,j=0)

#uneven_levels = np.arange(min_std,max_std,0.1)
size = '10-25 sqr. degs.'
easy_plot_unevencolorbars(lon=xx,lat=yy, data=std_var_10_25,
                  uneven_levels=uneven_levels2,name= f'Std {size} ',units='Standard Deviation',i=2,j=1)

#uneven_levels = np.arange(min_median, max_median, 1)
size = '10-25 sqr. degs.'
easy_plot_unevencolorbars(lon=xx,lat=yy, data=median_var_10_25,
                  uneven_levels=uneven_levels,name=f'Median {size} ',units=f'{units}',i=2,j=2)

#uneven_levels = np.arange(min_mean, max_mean, 1)
size = '25-50 sqr. degs.'
easy_plot_unevencolorbars(lon=xx,lat=yy, data=mean_var_25_50,
                  uneven_levels=uneven_levels,name= f'Mean {variable} {size} ',units=f'{units}',i=3,j=0)

#uneven_levels = np.arange(min_std,max_std,0.1)
size = '25-50 sqr. degs.'
easy_plot_unevencolorbars(lon=xx,lat=yy, data=std_var_25_50,
                  uneven_levels=uneven_levels2,name= f'Std {size} ',units='Standard Deviation',i=3,j=1)

#uneven_levels = np.arange(min_median, max_median, 1)
size = '25-50 sqr. degs.'
easy_plot_unevencolorbars(lon=xx,lat=yy, data=median_var_25_50,
                  uneven_levels=uneven_levels,name=f'Median {size}',units=f'{units}',i=3,j=2)

#uneven_levels = np.arange(min_mean, max_mean, 1)
size = '50-100 sqr. degs.'
easy_plot_unevencolorbars(lon=xx,lat=yy, data=mean_var_50_100,
                  uneven_levels=uneven_levels,name= f'Mean {variable} {size} ',units=f'{units}',i=4,j=0)

#uneven_levels = np.arange(min_std,max_std,0.1)
size = '50-100 sqr. degs.'
easy_plot_unevencolorbars(lon=xx,lat=yy, data=std_var_50_100,
                  uneven_levels=uneven_levels2,name= f'Std {size} ',units='Standard Deviation',i=4,j=1)

#uneven_levels = np.arange(min_median, max_median, 1)
size = '50-100 sqr. degs.'
easy_plot_unevencolorbars(lon=xx,lat=yy, data=median_var_50_100,
                  uneven_levels=uneven_levels,name=f'Median {size} ',units=f'{units}',i=4,j=2)

#uneven_levels = np.arange(min_mean, max_mean, 1)
size = '100-200 sqr. degs.'
easy_plot_unevencolorbars(lon=xx,lat=yy, data=mean_var_100_200,
                  uneven_levels=uneven_levels,name= f'Mean {variable} {size} ',units=f'{units}',i=5,j=0)

#uneven_levels = np.arange(min_std,max_std,0.1)
size = '100-200 sqr. degs.'
easy_plot_unevencolorbars(lon=xx,lat=yy, data=std_var_100_200,
                  uneven_levels=uneven_levels2,name= f'Std {size} ',units='Standard Deviation',i=5,j=1)

#uneven_levels = np.arange(min_median, max_median, 1)
size = '100-200 sqr. degs.'
easy_plot_unevencolorbars(lon=xx,lat=yy, data=median_var_100_200,
                  uneven_levels=uneven_levels,name=f'Median {size} ',units=f'{units}',i=5,j=2)

#uneven_levels = np.arange(min_mean, max_mean, 1)
size = '200-500 sqr. degs.'
easy_plot_unevencolorbars(lon=xx,lat=yy, data=mean_var_200_500,
                  uneven_levels=uneven_levels,name= f'Mean {variable} {size} ',units=f'{units}',i=6,j=0)

#uneven_levels = np.arange(min_std,max_std,0.1)
size = '200-500 sqr. degs.'
easy_plot_unevencolorbars(lon=xx,lat=yy, data=std_var_200_500,
                  uneven_levels=uneven_levels2,name= f'Std {size} ',units='Standard Deviation',i=6,j=1)

#uneven_levels = np.arange(min_median, max_median, 1)
size = '200-500 sqr. degs.'
easy_plot_unevencolorbars(lon=xx,lat=yy, data=median_var_200_500,
                  uneven_levels=uneven_levels,name=f'Median {size} ',units=f'{units}',i=6,j=2)

#uneven_levels = np.arange(min_mean, max_mean, 1)
size = '500-700 sqr. degs.'
easy_plot_unevencolorbars(lon=xx,lat=yy, data=mean_var_500_700,
                  uneven_levels=uneven_levels,name= f'Mean {variable} {size} ',units=f'{units}',i=7,j=0)

#uneven_levels = np.arange(min_std,max_std,0.1)
size = '500-700 sqr. degs.'
easy_plot_unevencolorbars(lon=xx,lat=yy, data=std_var_500_700,
                  uneven_levels=uneven_levels2,name= f'Std {size} ',units='Standard Deviation',i=7,j=1)

#uneven_levels = np.arange(min_median, max_median, 1)
size = '500-700 sqr. degs.'
easy_plot_unevencolorbars(lon=xx,lat=yy, data=median_var_500_700,
                  uneven_levels=uneven_levels,name=f'Median {size} ',units=f'{units}',i=7,j=2)




axs[0, 0].text(-0.2, 1.05, 'a', transform=axs[0, 0].transAxes, fontsize=20, fontweight='bold')
axs[0, 1].text(-0.2, 1.05, 'b', transform=axs[0, 1].transAxes, fontsize=20, fontweight='bold')
axs[0, 2].text(-0.2, 1.05, 'C', transform=axs[0, 2].transAxes, fontsize=20, fontweight='bold')

axs[1, 0].text(-0.2, 1.05, 'd', transform=axs[1, 0].transAxes, fontsize=20, fontweight='bold')
axs[1, 1].text(-0.2, 1.05, 'e', transform=axs[1, 1].transAxes, fontsize=20, fontweight='bold')
axs[1, 2].text(-0.2, 1.05, 'f', transform=axs[1, 2].transAxes, fontsize=20, fontweight='bold')

axs[2, 0].text(-0.2, 1.05, 'g', transform=axs[2, 0].transAxes, fontsize=20, fontweight='bold')
axs[2, 1].text(-0.2, 1.05, 'h', transform=axs[2, 1].transAxes, fontsize=20, fontweight='bold')
axs[2, 2].text(-0.2, 1.05, 'i', transform=axs[2, 2].transAxes, fontsize=20, fontweight='bold')

axs[3, 0].text(-0.2, 1.05, 'j', transform=axs[3, 0].transAxes, fontsize=20, fontweight='bold')
axs[3, 1].text(-0.2, 1.05, 'k', transform=axs[3, 1].transAxes, fontsize=20, fontweight='bold')
axs[3, 2].text(-0.2, 1.05, 'l', transform=axs[3, 2].transAxes, fontsize=20, fontweight='bold')


axs[4, 0].text(-0.2, 1.05, 'm', transform=axs[4, 0].transAxes, fontsize=20, fontweight='bold')
axs[4, 1].text(-0.2, 1.05, 'n', transform=axs[4, 1].transAxes, fontsize=20, fontweight='bold')
axs[4, 2].text(-0.2, 1.05, 'o', transform=axs[4, 2].transAxes, fontsize=20, fontweight='bold')


axs[5, 0].text(-0.2, 1.05, 'p', transform=axs[5, 0].transAxes, fontsize=20, fontweight='bold')
axs[5, 1].text(-0.2, 1.05, 'q', transform=axs[5, 1].transAxes, fontsize=20, fontweight='bold')
axs[5, 2].text(-0.2, 1.05, 'r', transform=axs[5, 2].transAxes, fontsize=20, fontweight='bold')

axs[6, 0].text(-0.2, 1.05, 's', transform=axs[6, 0].transAxes, fontsize=20, fontweight='bold')
axs[6, 1].text(-0.2, 1.05, 't', transform=axs[6, 1].transAxes, fontsize=20, fontweight='bold')
axs[6, 2].text(-0.2, 1.05, 'u', transform=axs[6, 2].transAxes, fontsize=20, fontweight='bold')


axs[7, 0].text(-0.2, 1.05, 'v', transform=axs[7, 0].transAxes, fontsize=20, fontweight='bold')
axs[7, 1].text(-0.2, 1.05, 'w', transform=axs[7, 1].transAxes, fontsize=20, fontweight='bold')
axs[7, 2].text(-0.2, 1.05, 'x', transform=axs[7, 2].transAxes, fontsize=20, fontweight='bold')


plt.tight_layout()
plt.savefig('mean_median_std_recoverywindow_poly_sizerange_noaa_1.png',bbox_inches='tight')
