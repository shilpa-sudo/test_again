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

mean_var = xr.open_dataset('polyareas_removed_noaa_93_22_mean_std_median_ratedecline.nc')

#recovery_window_days_pacific
mean_var_all = np.array(mean_var.mean_ratedecline)
mean_var_0_1_rem = np.array(mean_var.mean_ratedecline_0_1_rem)
mean_var_0_10_rem = np.array(mean_var.mean_ratedecline_0_10_rem)
mean_var_0_25_rem = np.array(mean_var.mean_ratedecline_0_25_rem)
mean_var_0_50_rem = np.array(mean_var.mean_ratedecline_0_50_rem)
mean_var_0_100_rem = np.array(mean_var.mean_ratedecline_0_100_rem)

std_var_all = np.array(mean_var.std_ratedecline)
std_var_0_1_rem = np.array(mean_var.std_ratedecline_0_1_rem)
std_var_0_10_rem = np.array(mean_var.std_ratedecline_0_10_rem)
std_var_0_25_rem = np.array(mean_var.std_ratedecline_0_25_rem)
std_var_0_50_rem = np.array(mean_var.std_ratedecline_0_50_rem)
std_var_0_100_rem = np.array(mean_var.std_ratedecline_0_100_rem)

median_var_all = np.array(mean_var.median_ratedecline)
median_var_0_1_rem = np.array(mean_var.median_ratedecline_0_1_rem)
median_var_0_10_rem = np.array(mean_var.median_ratedecline_0_10_rem)
median_var_0_25_rem = np.array(mean_var.median_ratedecline_0_25_rem)
median_var_0_50_rem = np.array(mean_var.median_ratedecline_0_50_rem)
median_var_0_100_rem = np.array(mean_var.median_ratedecline_0_100_rem)


print (np.nanmean(mean_var_all),np.nanmin(mean_var_all),np.nanmax(mean_var_all))
print (np.nanmean(mean_var_0_1_rem),np.nanmin(mean_var_0_1_rem),np.nanmax(mean_var_0_1_rem))
print (np.nanmean(mean_var_0_10_rem),np.nanmin(mean_var_0_10_rem),np.nanmax(mean_var_0_10_rem))
print (np.nanmean(mean_var_0_25_rem),np.nanmin(mean_var_0_25_rem),np.nanmax(mean_var_0_25_rem))
print (np.nanmean(mean_var_0_50_rem),np.nanmin(mean_var_0_50_rem),np.nanmax(mean_var_0_50_rem))
print (np.nanmean(mean_var_0_100_rem),np.nanmin(mean_var_0_100_rem),np.nanmax(mean_var_0_100_rem))

print (np.nanmean(std_var_all),np.nanmin(std_var_all),np.nanmax(std_var_all))
print (np.nanmean(std_var_0_1_rem),np.nanmin(std_var_0_1_rem),np.nanmax(std_var_0_1_rem))
print (np.nanmean(std_var_0_10_rem),np.nanmin(std_var_0_10_rem),np.nanmax(std_var_0_10_rem))
print (np.nanmean(std_var_0_25_rem),np.nanmin(std_var_0_25_rem),np.nanmax(std_var_0_25_rem))
print (np.nanmean(std_var_0_50_rem),np.nanmin(std_var_0_50_rem),np.nanmax(std_var_0_50_rem))
print (np.nanmean(std_var_0_100_rem),np.nanmin(std_var_0_100_rem),np.nanmax(std_var_0_100_rem))

print (np.nanmean(median_var_all),np.nanmin(median_var_all),np.nanmax(median_var_all))
print (np.nanmean(median_var_0_1_rem),np.nanmin(median_var_0_1_rem),np.nanmax(median_var_0_1_rem))
print (np.nanmean(median_var_0_10_rem),np.nanmin(median_var_0_10_rem),np.nanmax(median_var_0_10_rem))
print (np.nanmean(median_var_0_25_rem),np.nanmin(median_var_0_25_rem),np.nanmax(median_var_0_25_rem))
print (np.nanmean(median_var_0_50_rem),np.nanmin(median_var_0_50_rem),np.nanmax(median_var_0_50_rem))
print (np.nanmean(median_var_0_100_rem),np.nanmin(median_var_0_100_rem),np.nanmax(median_var_0_100_rem))


min_mean=min([np.nanmin(mean_var_all),np.nanmin(mean_var_0_1_rem),np.nanmin(mean_var_0_10_rem),np.nanmin(mean_var_0_25_rem),np.nanmin(mean_var_0_50_rem),np.nanmin(mean_var_0_100_rem)])
max_mean=max([np.nanmax(mean_var_all),np.nanmax(mean_var_0_1_rem),np.nanmax(mean_var_0_10_rem),np.nanmax(mean_var_0_25_rem),np.nanmax(mean_var_0_50_rem),np.nanmax(mean_var_0_100_rem)])
min_std=min([np.nanmin(std_var_all),np.nanmin(std_var_0_1_rem),np.nanmin(std_var_0_10_rem),np.nanmin(std_var_0_25_rem),np.nanmin(std_var_0_50_rem),np.nanmin(std_var_0_100_rem)])
max_std=max([np.nanmax(std_var_all),np.nanmax(std_var_0_1_rem),np.nanmax(std_var_0_10_rem),np.nanmax(std_var_0_25_rem),np.nanmax(std_var_0_50_rem),np.nanmax(std_var_0_100_rem)])
min_median=min([np.nanmin(median_var_all),np.nanmin(median_var_0_1_rem),np.nanmin(median_var_0_10_rem),np.nanmin(median_var_0_25_rem),np.nanmin(median_var_0_50_rem),np.nanmin(median_var_0_100_rem)])
max_median=max([np.nanmax(median_var_all),np.nanmax(median_var_0_1_rem),np.nanmax(median_var_0_10_rem),np.nanmax(median_var_0_25_rem),np.nanmax(median_var_0_50_rem),np.nanmax(median_var_0_100_rem)])



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
lat = mean_var.lat
lon = mean_var.lon

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



nrows = 6
ncols = 3

fig, axs = plt.subplots(nrows=nrows,ncols=ncols,subplot_kw={'projection': proj_crs},figsize=(20,30))

uneven_levels = np.arange(min_mean, max_mean, 1)
easy_plot_unevencolorbars(lon=xx,lat=yy, data=mean_var_all,
                  uneven_levels=uneven_levels,name='Mean MHW Rate decline (All sizes)',units=' [°C/Day]',i=0,j=0)

uneven_levels2 = np.arange(min_std,max_std,0.1)
easy_plot_unevencolorbars(lon=xx,lat=yy, data=std_var_all,
                  uneven_levels=uneven_levels2,name='Std (All sizes)',units='Standard Deviation',i=0,j=1)

uneven_levels3 = np.arange(min_median, max_median, 1)
easy_plot_unevencolorbars(lon=xx,lat=yy, data=median_var_all,
                  uneven_levels=uneven_levels3,name='Median (All sizes)',units=' [°C/Day]',i=0,j=2)



uneven_levels = np.arange(min_mean, max_mean, 1)

easy_plot_unevencolorbars(lon=xx,lat=yy, data=mean_var_0_1_rem,
                  uneven_levels=uneven_levels,name='Mean MHW Rate decline (0-1 square deg. removed)',units=' [°C/Day]',i=1,j=0)

uneven_levels = np.arange(min_std,max_std,0.1)

easy_plot_unevencolorbars(lon=xx,lat=yy, data=std_var_0_1_rem,
                  uneven_levels=uneven_levels,name='Std (0-1 square deg. removed)',units='Standard Deviation',i=1,j=1)

uneven_levels = np.arange(min_median, max_median, 1)
easy_plot_unevencolorbars(lon=xx,lat=yy, data=median_var_0_1_rem,
                  uneven_levels=uneven_levels,name='Median ',units=' [°C/Day]',i=1,j=2)


uneven_levels = np.arange(min_mean, max_mean, 1)
easy_plot_unevencolorbars(lon=xx,lat=yy, data=mean_var_0_10_rem,
                  uneven_levels=uneven_levels,name='Mean MHW Rate decline (0-10 square deg. removed)',units=' [°C/Day]',i=2,j=0)

uneven_levels = np.arange(min_std,max_std,0.1)

easy_plot_unevencolorbars(lon=xx,lat=yy, data=std_var_0_10_rem,
                  uneven_levels=uneven_levels,name='Std (0-10 square deg. removed)',units='Standard Deviation',i=2,j=1)

uneven_levels = np.arange(min_median, max_median, 1)
easy_plot_unevencolorbars(lon=xx,lat=yy, data=median_var_0_10_rem,
                  uneven_levels=uneven_levels,name='Median ',units=' [°C/Day]',i=2,j=2)


uneven_levels = np.arange(min_mean, max_mean, 1)
easy_plot_unevencolorbars(lon=xx,lat=yy, data=mean_var_0_25_rem,
                  uneven_levels=uneven_levels,name='Mean MHW Rate decline (0-25 square deg. removed)',units=' [°C/Day]',i=3,j=0)

uneven_levels = np.arange(min_std,max_std,0.1)

easy_plot_unevencolorbars(lon=xx,lat=yy, data=std_var_0_25_rem,
                  uneven_levels=uneven_levels,name='Std (0-25 square deg. removed)',units='Standard Deviation',i=3,j=1)

uneven_levels = np.arange(min_median, max_median, 1)
easy_plot_unevencolorbars(lon=xx,lat=yy, data=median_var_0_25_rem,
                  uneven_levels=uneven_levels,name='Median ',units=' [°C/Day]',i=3,j=2)


uneven_levels = np.arange(min_mean, max_mean, 1)

easy_plot_unevencolorbars(lon=xx,lat=yy, data=mean_var_0_50_rem,
                  uneven_levels=uneven_levels,name='Mean MHW Rate decline (0-50 square deg. removed)',units=' [°C/Day]',i=4,j=0)

uneven_levels = np.arange(min_std,max_std,0.1)

easy_plot_unevencolorbars(lon=xx,lat=yy, data=std_var_0_50_rem,
                  uneven_levels=uneven_levels,name='Std (0-50 square deg. removed)',units='Standard Deviation',i=4,j=1)

uneven_levels = np.arange(min_median, max_median, 1)
easy_plot_unevencolorbars(lon=xx,lat=yy, data=median_var_0_50_rem,
                  uneven_levels=uneven_levels,name='Median ',units=' [°C/Day]',i=4,j=2)



uneven_levels = np.arange(min_mean, max_mean, 1)

easy_plot_unevencolorbars(lon=xx,lat=yy, data=mean_var_0_100_rem,
                  uneven_levels=uneven_levels,name='Mean MHW Rate decline (0-100 square deg. removed)',units=' [°C/Day]',i=5,j=0)

uneven_levels = np.arange(min_std,max_std,0.1)

easy_plot_unevencolorbars(lon=xx,lat=yy, data=std_var_0_100_rem,
                  uneven_levels=uneven_levels,name='Std (0-100 square deg. removed)',units='Standard Deviation',i=5,j=1)

uneven_levels = np.arange(min_median, max_median, 1)
easy_plot_unevencolorbars(lon=xx,lat=yy, data=median_var_0_100_rem,
                  uneven_levels=uneven_levels,name='Median ',units=' [°C/Day]',i=5,j=2)

axs[0, 0].text(-0.2, 1.05, 'a', transform=axs[0, 0].transAxes, fontsize=20, fontweight='bold')
axs[0, 1].text(-0.2, 1.05, 'b', transform=axs[0, 1].transAxes, fontsize=20, fontweight='bold')
axs[0, 2].text(-0.2, 1.05, 'C', transform=axs[0, 2].transAxes, fontsize=20, fontweight='bold')

axs[1, 0].text(-0.2, 1.05, 'd', transform=axs[1, 0].transAxes, fontsize=20, fontweight='bold')
axs[1, 1].text(-0.2, 1.05, 'e', transform=axs[1, 1].transAxes, fontsize=20, fontweight='bold')
axs[1, 2].text(-0.2, 1.05, 'f', transform=axs[1, 2].transAxes, fontsize=20, fontweight='bold')

axs[2, 0].text(-0.2, 1.05, 'g', transform=axs[2, 0].transAxes, fontsize=20, fontweight='bold')
axs[2, 1].text(-0.2, 1.05, 'h', transform=axs[2, 1].transAxes, fontsize=20, fontweight='bold')
axs[2, 2].text(-0.2, 1.05, 'i', transform=axs[2, 2].transAxes, fontsize=20, fontweight='bold')

axs[3, 0].text(-0.2, 1.05, 'g', transform=axs[3, 0].transAxes, fontsize=20, fontweight='bold')
axs[3, 1].text(-0.2, 1.05, 'h', transform=axs[3, 1].transAxes, fontsize=20, fontweight='bold')
axs[3, 2].text(-0.2, 1.05, 'h', transform=axs[3, 2].transAxes, fontsize=20, fontweight='bold')


axs[4, 0].text(-0.2, 1.05, 'i', transform=axs[4, 0].transAxes, fontsize=20, fontweight='bold')
axs[4, 1].text(-0.2, 1.05, 'j', transform=axs[4, 1].transAxes, fontsize=20, fontweight='bold')
axs[4, 2].text(-0.2, 1.05, 'h', transform=axs[4, 2].transAxes, fontsize=20, fontweight='bold')


axs[5, 0].text(-0.2, 1.05, 'k', transform=axs[5, 0].transAxes, fontsize=20, fontweight='bold')
axs[5, 1].text(-0.2, 1.05, 'l', transform=axs[5, 1].transAxes, fontsize=20, fontweight='bold')
axs[5, 2].text(-0.2, 1.05, 'h', transform=axs[5, 2].transAxes, fontsize=20, fontweight='bold')


plt.tight_layout()
plt.savefig('mean_median_ratedecline_polygonareasremoved_noaa.png',bbox_inches='tight')
