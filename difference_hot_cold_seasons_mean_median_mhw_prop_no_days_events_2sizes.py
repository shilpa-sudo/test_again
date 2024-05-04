# standard imports
import numpy as np
import xarray as xr
import numpy as np
from scipy.spatial import cKDTree
import matplotlib.pyplot as plt
from cartopy import config
import cartopy.crs as ccrs
from scipy.ndimage import convolve
import geopandas as gpd
import matplotlib.colors as mcolors
import seaborn as sns
import sys
import os
import netCDF4

from datetime import datetime
from pathlib import Path

from dask.diagnostics import ProgressBar

# get files


cold_duration = xr.open_dataset('/home/shilpa/glory_mat_analysis/50_700_cold_season_polyareas_specific_sizes_range_noaa_93_22_mean_std_median_duration.nc')
cold_events = xr.open_dataset('/home/shilpa/glory_mat_analysis/50_less_more_cold_season_no_of_events_polyareas_specific_sizes_range_noaa_93_22.nc')
cold_days = xr.open_dataset('/home/shilpa/glory_mat_analysis/50_less_more_cold_season_polyareas_specific_sizes_range_noaa_93_22_sum_active_mhw_days.nc')
cold_maxintensity= xr.open_dataset('/home/shilpa/glory_mat_analysis/50_less_more_cold_season_polyareas_specific_sizes_range_noaa_93_22_mean_std_median_intensitymax.nc')

hot_duration = xr.open_dataset('/home/shilpa/glory_mat_analysis/50_less_more_hot_season_polyareas_specific_sizes_range_noaa_93_22_mean_std_median_duration.nc')
hot_maxintensity = xr.open_dataset('/home/shilpa/glory_mat_analysis/50_less_more_hot_season_polyareas_specific_sizes_range_noaa_93_22_mean_std_median_intensity_max.nc')
hot_days = xr.open_dataset('/home/shilpa/glory_mat_analysis/50_less_more_hot_season_polyareas_specific_sizes_range_noaa_93_22_sum_active_mhw_days.nc')
hot_events = xr.open_dataset('/home/shilpa/glory_mat_analysis/50_more_less_hot_season_no_of_events_polyareas_specific_sizes_range_noaa_93_22.nc')


duration = xr.open_dataset('/home/shilpa/glory_mat_analysis/50_less_more_polyareas_specific_sizes_range_noaa_93_22_mean_std_median_duration.nc')
maxintensity = xr.open_dataset('/home/shilpa/glory_mat_analysis/50_less_more_polyareas_specific_sizes_range_noaa_93_22_mean_std_median_intensity_max_pacific.nc')
events = xr.open_dataset('/home/shilpa/glory_mat_analysis/50_mess_more_no_of_events_polyareas_specific_sizes_range_noaa_93_22.nc')
days = xr.open_dataset('/home/shilpa/glory_mat_analysis/50_700_polyareas_specific_sizes_range_noaa_93_22_sum_active_mhw_days.nc')

cold_events_all = xr.open_dataset('/home/shilpa/glory_mat_analysis/0_700_cold_season_no_of_events_polyareas_specific_sizes_range_noaa_93_22.nc')
cold_duration_all = xr.open_dataset('/home/shilpa/glory_mat_analysis/0_700_cold_season_polyareas_specific_sizes_range_noaa_93_22_mean_std_median_duration.nc')
cold_maxintensity_all = xr.open_dataset('/home/shilpa/glory_mat_analysis/0_700_cold_season_polyareas_specific_sizes_range_noaa_93_22_mean_std_median_intensitymax.nc')
cold_days_all = xr.open_dataset('/home/shilpa/glory_mat_analysis/0_700_cold_season_polyareas_specific_sizes_range_noaa_93_22_sum_active_mhw_days.nc')
hot_events_all = xr.open_dataset('/home/shilpa/glory_mat_analysis/0_700_hot_season_no_of_events_polyareas_specific_sizes_range_noaa_93_22.nc')
hot_duration_all = xr.open_dataset('/home/shilpa/glory_mat_analysis/0_700_hot_season_polyareas_specific_sizes_range_noaa_93_22_mean_std_median_duration.nc')
hot_maxintensity_all = xr.open_dataset('/home/shilpa/glory_mat_analysis/0_700_hot_season_polyareas_specific_sizes_range_noaa_93_22_mean_std_median_intensity_max.nc')
hot_days_all = xr.open_dataset('/home/shilpa/glory_mat_analysis/0_700_hot_season_polyareas_specific_sizes_range_noaa_93_22_sum_active_mhw_days.nc')
events_all = xr.open_dataset('/home/shilpa/glory_mat_analysis/0_700_no_of_events_polyareas_specific_sizes_range_noaa_93_22.nc')
maxintensity_all = xr.open_dataset('/home/shilpa/glory_mat_analysis/0_700_polyareas_specific_sizes_range_noaa_93_22_mean_std_median_intensity_max_pacific.nc')
duration_all = xr.open_dataset('/home/shilpa/glory_mat_analysis/0_700_polyareas_specific_sizes_range_noaa_93_22_mean_std_median_duration.nc')
days_all = xr.open_dataset('/home/shilpa/glory_mat_analysis/0_700_polyareas_specific_sizes_range_noaa_93_22_sum_active_mhw_days.nc')


                          
longi = cold_duration.lon
lati = cold_duration.lat


xx,yy = np.meshgrid(longi,lati)


# get EEZ
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

# get projections
data_crs = ccrs.PlateCarree(central_longitude=0)

# But you want the map to include the dateline.
proj_crs = ccrs.PlateCarree(central_longitude=180)



# get script
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



    eez[0:1].plot(ax=axs[i,j],color='none',edgecolor='black',transform=data_crs)
    eez[2:3].plot(ax=axs[i,j],color='none',edgecolor='black',transform=data_crs)
    eez[4:5].plot(ax=axs[i,j],color='none',edgecolor='black',transform=data_crs)
    eez[7:9].plot(ax=axs[i,j],color='none',edgecolor='black',transform=data_crs)
    eez[11:12].plot(ax=axs[i,j],color='none',edgecolor='black',transform=data_crs)
    eez[15:16].plot(ax=axs[i,j],color='none',edgecolor='black',transform=data_crs)
    eez[20:21].plot(ax=axs[i,j],color='none',edgecolor='black',transform=data_crs)
    eez[33:35].plot(ax=axs[i,j],color='none',edgecolor='black',transform=data_crs)
    eez[37:38].plot(ax=axs[i,j],color='none',edgecolor='black',transform=data_crs)
    eez[47:48].plot(ax=axs[i,j],color='none',edgecolor='black',transform=data_crs)
    eez[51:52].plot(ax=axs[i,j],color='none',edgecolor='black',transform=data_crs)
    eez[53:54].plot(ax=axs[i,j],color='none',edgecolor='black',transform=data_crs)
    eez[248:250].plot(ax=axs[i,j],color='none',edgecolor='black',transform=data_crs)
    eez[243:245].plot(ax=axs[i,j],color='none',edgecolor='black',transform=data_crs)
    eez[246:248].plot(ax=axs[i,j],color='none',edgecolor='black',transform=data_crs)
    eez[251:253].plot(ax=axs[i,j],color='none',edgecolor='black',transform=data_crs)
    eez[259:260].plot(ax=axs[i,j],color='none',edgecolor='black',transform=data_crs)
    eez[261:262].plot(ax=axs[i,j],color='none',edgecolor='black',transform=data_crs)

    cbar = plt.colorbar(im, ax=axs[i,j],orientation="horizontal", shrink=0.5)
    cbar.set_label(units)
    gl = axs[i,j].gridlines(draw_labels=True)
    gl.top_labels = False
    gl.right_labels = False
    axs[i,j].set_xlim(-35,29) #-35
    axs[i,j].set_ylim(-34.875,-2.375)
    axs[i,j].set_title(name, fontsize=9)


def easy_plot_unevencolorbars2(lon,lat,data,i,j,name,units,uneven_levels):
    uneven_levels = uneven_levels
    cmap_rb = plt.get_cmap('RdBu_r')#'RdBu_r')
    colors = cmap_rb(np.linspace(0, 1, (len(uneven_levels) - 1)))
    cmap, norm = mcolors.from_levels_and_colors(uneven_levels, colors)
    
    im = axs[i,j].contourf(lon,lat,data, 
               levels=uneven_levels,
               cmap=cmap, norm=norm,
               transform=data_crs,
               transform_first=True)
    #axs[i,j].coastlines()
    #eez[:].plot(ax=axs[i,j],color='none',edgecolor='black',transform=data_crs)
    
    
    
    eez[0:1].plot(ax=axs[i,j],color='none',edgecolor='black',transform=data_crs)
    eez[2:3].plot(ax=axs[i,j],color='none',edgecolor='black',transform=data_crs)
    eez[4:5].plot(ax=axs[i,j],color='none',edgecolor='black',transform=data_crs)
    eez[7:9].plot(ax=axs[i,j],color='none',edgecolor='black',transform=data_crs)
    eez[11:12].plot(ax=axs[i,j],color='none',edgecolor='black',transform=data_crs)
    eez[15:16].plot(ax=axs[i,j],color='none',edgecolor='black',transform=data_crs)
    eez[20:21].plot(ax=axs[i,j],color='none',edgecolor='black',transform=data_crs)
    eez[33:35].plot(ax=axs[i,j],color='none',edgecolor='black',transform=data_crs)
    eez[37:38].plot(ax=axs[i,j],color='none',edgecolor='black',transform=data_crs)
    eez[47:48].plot(ax=axs[i,j],color='none',edgecolor='black',transform=data_crs)
    eez[51:52].plot(ax=axs[i,j],color='none',edgecolor='black',transform=data_crs)
    eez[53:54].plot(ax=axs[i,j],color='none',edgecolor='black',transform=data_crs)
    eez[248:250].plot(ax=axs[i,j],color='none',edgecolor='black',transform=data_crs)
    eez[243:245].plot(ax=axs[i,j],color='none',edgecolor='black',transform=data_crs)
    eez[246:248].plot(ax=axs[i,j],color='none',edgecolor='black',transform=data_crs)
    eez[251:253].plot(ax=axs[i,j],color='none',edgecolor='black',transform=data_crs)
    eez[259:260].plot(ax=axs[i,j],color='none',edgecolor='black',transform=data_crs)
    eez[261:262].plot(ax=axs[i,j],color='none',edgecolor='black',transform=data_crs)
    
    cbar = plt.colorbar(im, ax=axs[i,j],orientation="horizontal", shrink=0.5)
    cbar.set_label(units)
    gl = axs[i,j].gridlines(draw_labels=True)
    gl.top_labels = False
    gl.right_labels = False
    axs[i,j].set_xlim(-35,29) #-35
    axs[i,j].set_ylim(-34.875,-2.375)
    axs[i,j].set_title(name, fontsize=9) 



mean_cold_dur_0_700 = np.array(cold_duration_all.mean_cold_season_duration_pacific_0_700)
mean_cold_dur_0_50 = np.array(cold_duration.mean_cold_season_duration_pacific_0_50)
mean_cold_dur_50_700 = np.array(cold_duration.mean_cold_season_duration_pacific_50_700)

median_cold_dur_0_700 = np.array(cold_duration_all.median_cold_season_duration_pacific_0_700)
median_cold_dur_0_50 = np.array(cold_duration.median_cold_season_duration_pacific_0_50)
median_cold_dur_50_700 = np.array(cold_duration.median_cold_season_duration_pacific_50_700)


mean_hot_dur_0_700 = np.array(hot_duration_all.mean_hot_season_duration_pacific_0_700)
mean_hot_dur_0_50 = np.array(hot_duration.mean_hot_season_duration_pacific_0_50)
mean_hot_dur_50_700 = np.array(hot_duration.mean_hot_season_duration_pacific_50_700)


median_hot_dur_0_700 = np.array(hot_duration_all.median_hot_season_duration_pacific_0_700)
median_hot_dur_0_50 = np.array(hot_duration.median_hot_season_duration_pacific_0_50)
median_hot_dur_50_700 = np.array(hot_duration.median_hot_season_duration_pacific_50_700)


mean_cold_max_int_0_700 = np.array(cold_maxintensity_all.mean_cold_season_intensity_max_pacific_0_700)
mean_cold_max_int_0_50 = np.array(cold_maxintensity.mean_cold_season_intensity_max_pacific_0_50)
mean_cold_max_int_50_700 = np.array(cold_maxintensity.mean_cold_season_intensity_max_pacific_50_700)

median_cold_max_int_0_700 = np.array(cold_maxintensity_all.median_cold_season_intensity_max_pacific_0_700)
median_cold_max_int_0_50 = np.array(cold_maxintensity.median_cold_season_intensity_max_pacific_0_50)
median_cold_max_int_50_700 = np.array(cold_maxintensity.median_cold_season_intensity_max_pacific_50_700)


mean_hot_max_int_0_700 = np.array(hot_maxintensity_all.mean_hot_season_intensity_max_pacific_0_700)
mean_hot_max_int_0_50 = np.array(hot_maxintensity.mean_hot_season_intensity_max_pacific_0_50)
mean_hot_max_int_50_700 = np.array(hot_maxintensity.mean_hot_season_intensity_max_pacific_50_700)


median_hot_max_int_0_700 = np.array(hot_maxintensity_all.median_hot_season_intensity_max_pacific_0_700)
median_hot_max_int_0_50 = np.array(hot_maxintensity.median_hot_season_intensity_max_pacific_0_50)
median_hot_max_int_50_700 = np.array(hot_maxintensity.median_hot_season_intensity_max_pacific_50_700)


hot_days_0_700 = np.array(hot_days_all.hot_sum_mhw_days_pacific_0_700)
hot_days_0_50 = np.array(hot_days.hot_sum_mhw_days_pacific_0_50)
hot_days_50_700 = np.array(hot_days.hot_sum_mhw_days_pacific_50_700)


cold_days_0_700 = np.array(cold_days_all.cold_sum_mhw_days_pacific_0_700)
cold_days_0_50 = np.array(cold_days.cold_sum_mhw_days_pacific_0_50)
cold_days_50_700 = np.array(cold_days.cold_sum_mhw_days_pacific_50_700)


hot_events_0_700 = np.array(hot_events_all.hot_season_no_mhw_events_pacific_0_700)
hot_events_0_50 = np.array(hot_events.hot_season_no_mhw_events_pacific_0_50)
hot_events_50_700 = np.array(hot_events.hot_season_no_mhw_events_pacific_50_700)

cold_events_0_700 = np.array(cold_events_all.cold_season_no_mhw_events_pacific_0_700)
cold_events_0_50 = np.array(cold_events.cold_season_no_mhw_events_pacific_0_50)
cold_events_50_700 = np.array(cold_events.cold_season_no_mhw_events_pacific_50_700)

mean_dur_0_700 = np.array(duration_all.mean_duration_pacific_0_700) 
mean_dur_0_50 = np.array(duration.mean_duration_pacific_0_50) 
mean_dur_50_700 = np.array(duration.mean_duration_pacific_50_700)

median_dur_0_700 = np.array(duration_all.median_duration_pacific_0_700)
median_dur_0_50 = np.array(duration.median_duration_pacific_0_50)
median_dur_50_700 = np.array(duration.median_duration_pacific_50_700)

mean_maxintensity_0_700 = np.array(maxintensity_all.mean_maxintensity_pacific_0_700)
mean_maxintensity_0_50 = np.array(maxintensity.mean_maxintensity_pacific_0_50)
mean_maxintensity_50_700 = np.array(maxintensity.mean_maxintensit_pacific_50_700)

median_maxintensity_0_700 = np.array(maxintensity_all.median_maxintensity_pacific_0_700)
median_maxintensity_0_50 = np.array(maxintensity.median_maxintensity_pacific_0_50)
median_maxintensity_50_700 = np.array(maxintensity.median_maxintensit_pacific_50_700)


events_0_700 = np.array(events_all.No_mhw_events_pacific_0_700)
events_0_50 = np.array(events.No_mhw_events_pacific_0_50)
events_50_700 = np.array(events.No_mhw_events_pacific_50_700)

days_0_700 = np.array(days_all.sum_mhw_days_pacific_0_700)
days_0_50 = np.array(days.sum_mhw_days_pacific_0_50)
days_50_700 = np.array(days.sum_mhw_days_pacific_50_700)





# plot
nrows=3
ncols=4

cbarunits = 'Days'
variablename = 'Mean Duration'

fig, axs = plt.subplots(nrows=nrows,ncols=ncols,subplot_kw={'projection': proj_crs},figsize=(22,30))

uneven_levels1 = [5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,40,45,50,55,60,65,70,80,90,100,150,200,250,300,350]

uneven_levels2 = [-350,-100,-50,-10,-5,0,5,10,50,100,350]


size = '0-700 sqr.degs.'

easy_plot_unevencolorbars(lon=xx,lat=yy, data=mean_dur_0_700,
                  uneven_levels=uneven_levels1,name=f' {variablename} {size}',units=f'{cbarunits}',i=0,j=0)

easy_plot_unevencolorbars(lon=xx,lat=yy, data=mean_hot_dur_0_700,
                  uneven_levels=uneven_levels1,name=f' Hot Season {variablename} {size}',units=f'{cbarunits}',i=0,j=1)

easy_plot_unevencolorbars(lon=xx,lat=yy, data=mean_cold_dur_0_700,
                  uneven_levels=uneven_levels1,name= f' Cold Season {variablename} {size}',units=f'{cbarunits}',i=0,j=2)

easy_plot_unevencolorbars2(lon=xx,lat=yy, data=mean_hot_dur_0_700-mean_cold_dur_0_700,
                  uneven_levels=uneven_levels2,name=f'Difference {size}',units=f'{cbarunits}',i=0,j=3)

size = '0-50 sqr.degs.'

easy_plot_unevencolorbars(lon=xx,lat=yy, data=mean_dur_0_50,
                  uneven_levels=uneven_levels1,name=f' {variablename} {size}',units=f'{cbarunits}',i=1,j=0)

easy_plot_unevencolorbars(lon=xx,lat=yy, data=mean_hot_dur_0_50,
                  uneven_levels=uneven_levels1,name=f' Hot Season {variablename} {size}',units=f'{cbarunits}',i=1,j=1)

easy_plot_unevencolorbars(lon=xx,lat=yy, data=mean_cold_dur_0_50,
                  uneven_levels=uneven_levels1,name= f' Cold Season {variablename} {size}',units=f'{cbarunits}',i=1,j=2)

easy_plot_unevencolorbars2(lon=xx,lat=yy, data=mean_hot_dur_0_50-mean_cold_dur_0_50,
                  uneven_levels=uneven_levels2,name=f'Difference {size}',units=f'{cbarunits}',i=1,j=3)

size = '50-700 sqr.degs.'

easy_plot_unevencolorbars(lon=xx,lat=yy, data=mean_dur_50_700,
                  uneven_levels=uneven_levels1,name=f' {variablename} {size}',units=f'{cbarunits}',i=2,j=0)

easy_plot_unevencolorbars(lon=xx,lat=yy, data=mean_hot_dur_50_700,
                  uneven_levels=uneven_levels1,name=f' Hot Season {variablename} {size}',units=f'{cbarunits}',i=2,j=1)

easy_plot_unevencolorbars(lon=xx,lat=yy, data=mean_cold_dur_50_700,
                  uneven_levels=uneven_levels1,name= f' Cold Season {variablename} {size}',units=f'{cbarunits}',i=2,j=2)

easy_plot_unevencolorbars2(lon=xx,lat=yy, data=mean_hot_dur_50_700-mean_cold_dur_50_700,
                  uneven_levels=uneven_levels2,name=f'Difference {size}',units=f'{cbarunits}',i=2,j=3)

plt.tight_layout()
plt.savefig(f'/home/shilpa/glory_mat_analysis/compare_mean_duration_hot_cold_seasons2.png')


# plot
nrows=3
ncols=4

cbarunits = 'Days'
variablename = 'Median Duration'

fig, axs = plt.subplots(nrows=nrows,ncols=ncols,subplot_kw={'projection': proj_crs},figsize=(22,30))

uneven_levels1 = [5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,40,45,50,55,60,65,70,80,90,100,150,200,250,300,350]

uneven_levels2 = [-350,-100,-50,-10,-5,0,5,10,50,100,350]

size = '0-1 sqr.degs.'

easy_plot_unevencolorbars(lon=xx,lat=yy, data=median_hot_dur_0_1,
                  uneven_levels=uneven_levels1,name=f' Hot Season {variablename} {size}',units=f'{cbarunits}',i=0,j=0)

easy_plot_unevencolorbars(lon=xx,lat=yy, data=median_cold_dur_0_1,
                  uneven_levels=uneven_levels1,name= f' Cold Season {variablename} {size}',units=f'{cbarunits}',i=0,j=1)

easy_plot_unevencolorbars2(lon=xx,lat=yy, data=median_hot_dur_0_1-mean_cold_dur_0_1,
                  uneven_levels=uneven_levels2,name=f'Difference {size}',units=f'{cbarunits}',i=0,j=2)

size = '1-10 sqr.degs.'

easy_plot_unevencolorbars(lon=xx,lat=yy, data=median_hot_dur_1_10,
                  uneven_levels=uneven_levels1,name= f' Hot Season {variablename} {size}',units=f'{cbarunits}',i=1,j=0)

easy_plot_unevencolorbars(lon=xx,lat=yy, data=median_cold_dur_1_10,
                  uneven_levels=uneven_levels1,name=f' Cold Season {variablename} {size}',units=f'{cbarunits}',i=1,j=1)

easy_plot_unevencolorbars2(lon=xx,lat=yy, data=median_hot_dur_1_10-mean_cold_dur_1_10,
                  uneven_levels=uneven_levels2,name=f'Difference {size}',units = f'{cbarunits}',i=1,j=2)

size = '10-25 sqr.degs.'

easy_plot_unevencolorbars(lon=xx,lat=yy, data=median_hot_dur_10_25,
                  uneven_levels=uneven_levels1,name=f' Hot Season {variablename} {size}',units= f'{cbarunits}',i=2,j=0)

easy_plot_unevencolorbars(lon=xx,lat=yy, data=median_cold_dur_10_25,
                  uneven_levels=uneven_levels1,name=f' Cold Season {variablename} {size}',units= f'{cbarunits}',i=2,j=1)

easy_plot_unevencolorbars2(lon=xx,lat=yy, data=median_hot_d
plt.tight_layout()
plt.savefig(f'/home/shilpa/glory_mat_analysis/compare_median_duration_hot_cold_seasons2.png')


# plot
nrows=8
ncols=3

cbarunits = '°C'
variablename = 'Mean Max Intensity'

fig, axs = plt.subplots(nrows=nrows,ncols=ncols,subplot_kw={'projection': proj_crs},figsize=(22,30))

uneven_levels1 = [0.5,0.6,0.7,0.8,0.9,1,1.2,1.4,1.6,1.8,2,2.2,2.4,2.6,2.8,3,3.5]

uneven_levels2 = [-3.5,-2,-1.25,-1,-0.8,-0.6,-0.4,-0.2,0,0.2,0.4,0.6,0.8,1,1.25,2,3.5]

size = '0-1 sqr.degs.'

easy_plot_unevencolorbars(lon=xx,lat=yy, data=mean_hot_max_int_0_1,
                  uneven_levels=uneven_levels1,name=f' Hot Season {variablename} {size}',units=f'{cbarunits}',i=0,j=0)

easy_plot_unevencolorbars(lon=xx,lat=yy, data=mean_cold_max_int_0_1,
                  uneven_levels=uneven_levels1,name= f' Cold Season {variablename} {size}',units=f'{cbarunits}',i=0,j=1)

easy_plot_unevencolorbars2(lon=xx,lat=yy, data=mean_hot_max_int_0_1-mean_cold_max_int_0_1,
                  uneven_levels=uneven_levels2,name=f'Difference {size}',units=f'{cbarunits}',i=0,j=2)

size = '1-10 sqr.degs.'

easy_plot_unevencolorbars(lon=xx,lat=yy, data=mean_hot_max_int_1_10,
                  uneven_levels=uneven_levels1,name= f' Hot Season {variablename} {size}',units=f'{cbarunits}',i=1,j=0)

easy_plot_unevencolorbars(lon=xx,lat=yy, data=mean_cold_max_int_1_10,
                  uneven_levels=uneven_levels1,name=f' Cold Season {variablename} {size}',units=f'{cbarunits}',i=1,j=1)

easy_plot_unevencolorbars2(lon=xx,lat=yy, data=mean_hot_max_int_1_10-mean_cold_max_int_1_10,
                  uneven_levels=uneven_levels2,name=f'Difference {size}',units = f'{cbarunits}',i=1,j=2)

size = '10-25 sqr.degs.'

easy_plot_unevencolorbars(lon=xx,lat=yy, data=mean_hot_max_int_10_25,
                  uneven_levels=uneven_levels1,name=f' Hot Season {variablename} {size}',units= f'{cbarunits}',i=2,j=0)

easy_plot_unevencolorbars(lon=xx,lat=yy, data=mean_cold_max_int_10_25,
                  uneven_levels=uneven_levels1,name=f' Cold Season {variablename} {size}',units= f'{cbarunits}',i=2,j=1)

easy_plot_unevencolorbars2(lon=xx,lat=yy, data=mean_hot_max_int_10_25-mean_cold_max_int_10_25,
                  uneven_levels=uneven_levels2,name=f'Difference {size}',units= f'{cbarunits}',i=2,j=2)

size = '25-50 sqr.degs.'

easy_plot_unevencolorbars(lon=xx,lat=yy, data=mean_hot_max_int_25_50,
                  uneven_levels=uneven_levels1,name=f' Hot Season {variablename} {size}',units= f'{cbarunits}',i=3,j=0)

easy_plot_unevencolorbars(lon=xx,lat=yy, data=mean_cold_max_int_25_50,
                  uneven_levels=uneven_levels1,name=f' Cold Season {variablename} {size}',units= f'{cbarunits}',i=3,j=1)

easy_plot_unevencolorbars2(lon=xx,lat=yy, data=mean_hot_max_int_25_50-mean_cold_max_int_25_50,
                  uneven_levels=uneven_levels2,name=f'Difference {size}', units = f'{cbarunits}',i=3,j=2)

size = '50-100 sqr.degs.'

easy_plot_unevencolorbars(lon=xx,lat=yy, data=mean_hot_max_int_50_100,
                  uneven_levels=uneven_levels1,name=f' Hot Season {variablename} {size}',units= f'{cbarunits}',i=4,j=0)

easy_plot_unevencolorbars(lon=xx,lat=yy, data=mean_cold_max_int_50_100,
                  uneven_levels=uneven_levels1,name=f' Cold Season {variablename} {size}',units= f'{cbarunits}',i=4,j=1)

easy_plot_unevencolorbars2(lon=xx,lat=yy, data=mean_hot_max_int_50_100-mean_cold_max_int_50_100,
                  uneven_levels=uneven_levels2,name=f'Difference {size}',units = f'{cbarunits}',i=4,j=2)

size = '100-200 sqr.degs.'

easy_plot_unevencolorbars(lon=xx,lat=yy, data=mean_hot_max_int_100_200,
                  uneven_levels=uneven_levels1,name=f' Hot Season {variablename} {size}',units= f'{cbarunits}',i=5,j=0)

easy_plot_unevencolorbars(lon=xx,lat=yy, data=mean_cold_max_int_100_200,
                  uneven_levels=uneven_levels1,name=f' Cold Season {variablename} {size}',units= f'{cbarunits}',i=5,j=1)

easy_plot_unevencolorbars2(lon=xx,lat=yy, data=mean_hot_max_int_100_200-mean_cold_max_int_100_200,
                  uneven_levels=uneven_levels2,name=f'Difference {size}',units= f'{cbarunits}',i=5,j=2)

size = '200-500 sqr.degs.'

easy_plot_unevencolorbars(lon=xx,lat=yy, data=mean_hot_max_int_200_500,
                  uneven_levels=uneven_levels1,name=f' Hot Season {variablename} {size}',units= f'{cbarunits}',i=6,j=0)

easy_plot_unevencolorbars(lon=xx,lat=yy, data=mean_cold_max_int_200_500,
                  uneven_levels=uneven_levels1,name=f' Cold Season {variablename} {size}',units= f'{cbarunits}',i=6,j=1)

easy_plot_unevencolorbars2(lon=xx,lat=yy, data=mean_hot_max_int_200_500-mean_cold_max_int_200_500,
                  uneven_levels=uneven_levels2,name=f'Difference {size}',units=f'{cbarunits}',i=6,j=2)

size = '500-700 sqr.degs.'

easy_plot_unevencolorbars(lon=xx,lat=yy, data=mean_hot_max_int_500_700,
                  uneven_levels=uneven_levels1,name= f' Hot Season {variablename} {size}',units= f'{cbarunits}',i=7,j=0)

easy_plot_unevencolorbars(lon=xx,lat=yy, data=mean_cold_max_int_500_700,
                  uneven_levels=uneven_levels1,name= f'Cold Season {variablename} {size}',units= f'{cbarunits}',i=7,j=1)

easy_plot_unevencolorbars2(lon=xx,lat=yy, data=mean_hot_max_int_500_700-mean_cold_max_int_500_700,
                  uneven_levels=uneven_levels2,name=f'Difference {size}',units= f'{cbarunits}',i=7,j=2)


plt.tight_layout()
plt.savefig(f'/home/shilpa/glory_mat_analysis/compare_mean_max_intensity_hot_cold_seasons1.png')



# plot
nrows=8
ncols=3

cbarunits = '°C'
variablename = 'Median Max Intensity'

fig, axs = plt.subplots(nrows=nrows,ncols=ncols,subplot_kw={'projection': proj_crs},figsize=(22,30))

uneven_levels1 = [0.5,0.6,0.7,0.8,0.9,1,1.2,1.4,1.6,1.8,2,2.2,2.4,2.6,2.8,3,3.5]

uneven_levels2 = [-3.5,-2,-1.25,-1,-0.8,-0.6,-0.4,-0.2,0,0.2,0.4,0.6,0.8,1,1.25,2,3.5]



size = '0-1 sqr.degs.'

easy_plot_unevencolorbars(lon=xx,lat=yy, data=median_hot_max_int_0_1,
                  uneven_levels=uneven_levels1,name=f' Hot Season {variablename} {size}',units=f'{cbarunits}',i=0,j=0)

easy_plot_unevencolorbars(lon=xx,lat=yy, data=median_cold_max_int_0_1,
                  uneven_levels=uneven_levels1,name= f' Cold Season {variablename} {size}',units=f'{cbarunits}',i=0,j=1)

easy_plot_unevencolorbars2(lon=xx,lat=yy, data=median_hot_max_int_0_1-median_cold_max_int_0_1,
                  uneven_levels=uneven_levels2,name=f'Difference {size}',units=f'{cbarunits}',i=0,j=2)

size = '1-10 sqr.degs.'

easy_plot_unevencolorbars(lon=xx,lat=yy, data=median_hot_max_int_1_10,
                  uneven_levels=uneven_levels1,name= f' Hot Season {variablename} {size}',units=f'{cbarunits}',i=1,j=0)

easy_plot_unevencolorbars(lon=xx,lat=yy, data=median_cold_max_int_1_10,
                  uneven_levels=uneven_levels1,name=f' Cold Season {variablename} {size}',units=f'{cbarunits}',i=1,j=1)

easy_plot_unevencolorbars2(lon=xx,lat=yy, data=median_hot_max_int_1_10-median_cold_max_int_1_10,
                  uneven_levels=uneven_levels2,name=f'Difference {size}',units = f'{cbarunits}',i=1,j=2)

size = '10-25 sqr.degs.'

easy_plot_unevencolorbars(lon=xx,lat=yy, data=median_hot_max_int_10_25,
                  uneven_levels=uneven_levels1,name=f' Hot Season {variablename} {size}',units= f'{cbarunits}',i=2,j=0)

easy_plot_unevencolorbars(lon=xx,lat=yy, data=median_cold_max_int_10_25,
                  uneven_levels=uneven_levels1,name=f' Cold Season {variablename} {size}',units= f'{cbarunits}',i=2,j=1)

easy_plot_unevencolorbars2(lon=xx,lat=yy, data=median_hot_max_int_10_25-median_cold_max_int_10_25,
                  uneven_levels=uneven_levels2,name=f'Difference {size}',units= f'{cbarunits}',i=2,j=2)

size = '25-50 sqr.degs.'

easy_plot_unevencolorbars(lon=xx,lat=yy, data=median_hot_max_int_25_50,
                  uneven_levels=uneven_levels1,name=f' Hot Season {variablename} {size}',units= f'{cbarunits}',i=3,j=0)

easy_plot_unevencolorbars(lon=xx,lat=yy, data=median_cold_max_int_25_50,
                  uneven_levels=uneven_levels1,name=f' Cold Season {variablename} {size}',units= f'{cbarunits}',i=3,j=1)

easy_plot_unevencolorbars2(lon=xx,lat=yy, data=median_hot_max_int_25_50-median_cold_max_int_25_50,
                  uneven_levels=uneven_levels2,name=f'Difference {size}', units = f'{cbarunits}',i=3,j=2)

size = '50-100 sqr.degs.'

easy_plot_unevencolorbars(lon=xx,lat=yy, data=median_hot_max_int_50_100,
                  uneven_levels=uneven_levels1,name=f' Hot Season {variablename} {size}',units= f'{cbarunits}',i=4,j=0)

easy_plot_unevencolorbars(lon=xx,lat=yy, data=median_cold_max_int_50_100,
                  uneven_levels=uneven_levels1,name=f' Cold Season {variablename} {size}',units= f'{cbarunits}',i=4,j=1)

easy_plot_unevencolorbars2(lon=xx,lat=yy, data=median_hot_max_int_50_100-median_cold_max_int_50_100,
                  uneven_levels=uneven_levels2,name=f'Difference {size}',units = f'{cbarunits}',i=4,j=2)

size = '100-200 sqr.degs.'

easy_plot_unevencolorbars(lon=xx,lat=yy, data=median_hot_max_int_100_200,
                  uneven_levels=uneven_levels1,name=f' Hot Season {variablename} {size}',units= f'{cbarunits}',i=5,j=0)

easy_plot_unevencolorbars(lon=xx,lat=yy, data=median_cold_max_int_100_200,
                  uneven_levels=uneven_levels1,name=f' Cold Season {variablename} {size}',units= f'{cbarunits}',i=5,j=1)

easy_plot_unevencolorbars2(lon=xx,lat=yy, data=median_hot_max_int_100_200-median_cold_max_int_100_200,
                  uneven_levels=uneven_levels2,name=f'Difference {size}',units= f'{cbarunits}',i=5,j=2)

size = '200-500 sqr.degs.'

easy_plot_unevencolorbars(lon=xx,lat=yy, data=median_hot_max_int_200_500,
                  uneven_levels=uneven_levels1,name=f' Hot Season {variablename} {size}',units= f'{cbarunits}',i=6,j=0)

easy_plot_unevencolorbars(lon=xx,lat=yy, data=median_cold_max_int_200_500,
                  uneven_levels=uneven_levels1,name=f' Cold Season {variablename} {size}',units= f'{cbarunits}',i=6,j=1)

easy_plot_unevencolorbars2(lon=xx,lat=yy, data=median_hot_max_int_200_500-median_cold_max_int_200_500,
                  uneven_levels=uneven_levels2,name=f'Difference {size}',units=f'{cbarunits}',i=6,j=2)

size = '500-700 sqr.degs.'

easy_plot_unevencolorbars(lon=xx,lat=yy, data=median_hot_max_int_500_700,
                  uneven_levels=uneven_levels1,name= f' Hot Season {variablename} {size}',units= f'{cbarunits}',i=7,j=0)

easy_plot_unevencolorbars(lon=xx,lat=yy, data=median_cold_max_int_500_700,
                  uneven_levels=uneven_levels1,name= f'Cold Season {variablename} {size}',units= f'{cbarunits}',i=7,j=1)

easy_plot_unevencolorbars2(lon=xx,lat=yy, data=median_hot_max_int_500_700-median_cold_max_int_500_700,
                  uneven_levels=uneven_levels2,name=f'Difference {size}',units= f'{cbarunits}',i=7,j=2)


plt.tight_layout()
plt.savefig(f'/home/shilpa/glory_mat_analysis/compare_median_max_intensity_hot_cold_seasons1.png')



# no. of mhw days 

# plot
nrows=8
ncols=3

cbarunits = 'Days'
variablename = 'No. of MHW Days'

fig, axs = plt.subplots(nrows=nrows,ncols=ncols,subplot_kw={'projection': proj_crs},figsize=(22,30))

uneven_levels1 = [0,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100,105,110,115,120,125,150,175,200,225,250,275,300,350,400,450,500,550,600,650,700]

uneven_levels2 = [-700,-500,-200,-100,-90,-80,-70,-60,-50,-40,-30,-20,-10,0,10,20,30,40,50,60,70,80,90,100,200,500,700]



size = '0-1 sqr.degs.'

easy_plot_unevencolorbars(lon=xx,lat=yy, data=hot_days_0_1,
                  uneven_levels=uneven_levels1,name=f' Hot Season {variablename} {size}',units=f'{cbarunits}',i=0,j=0)

easy_plot_unevencolorbars(lon=xx,lat=yy, data=cold_days_0_1,
                  uneven_levels=uneven_levels1,name= f' Cold Season {variablename} {size}',units=f'{cbarunits}',i=0,j=1)

easy_plot_unevencolorbars2(lon=xx,lat=yy, data=hot_days_0_1-cold_days_0_1,
                  uneven_levels=uneven_levels2,name=f'Difference {size}',units=f'{cbarunits}',i=0,j=2)

size = '1-10 sqr.degs.'

easy_plot_unevencolorbars(lon=xx,lat=yy, data=hot_days_1_10,
                  uneven_levels=uneven_levels1,name= f' Hot Season {variablename} {size}',units=f'{cbarunits}',i=1,j=0)

easy_plot_unevencolorbars(lon=xx,lat=yy, data=cold_days_1_10,
                  uneven_levels=uneven_levels1,name=f' Cold Season {variablename} {size}',units=f'{cbarunits}',i=1,j=1)

easy_plot_unevencolorbars2(lon=xx,lat=yy, data=hot_days_1_10-cold_days_1_10,
                  uneven_levels=uneven_levels2,name=f'Difference {size}',units = f'{cbarunits}',i=1,j=2)

size = '10-25 sqr.degs.'

easy_plot_unevencolorbars(lon=xx,lat=yy, data=hot_days_10_25,
                  uneven_levels=uneven_levels1,name=f' Hot Season {variablename} {size}',units= f'{cbarunits}',i=2,j=0)

easy_plot_unevencolorbars(lon=xx,lat=yy, data=cold_days_10_25,
                  uneven_levels=uneven_levels1,name=f' Cold Season {variablename} {size}',units= f'{cbarunits}',i=2,j=1)

easy_plot_unevencolorbars2(lon=xx,lat=yy, data=hot_days_10_25-cold_days_10_25,
                  uneven_levels=uneven_levels2,name=f'Difference {size}',units= f'{cbarunits}',i=2,j=2)

size = '25-50 sqr.degs.'

easy_plot_unevencolorbars(lon=xx,lat=yy, data=hot_days_25_50,
                  uneven_levels=uneven_levels1,name=f' Hot Season {variablename} {size}',units= f'{cbarunits}',i=3,j=0)

easy_plot_unevencolorbars(lon=xx,lat=yy, data=cold_days_25_50,
                  uneven_levels=uneven_levels1,name=f' Cold Season {variablename} {size}',units= f'{cbarunits}',i=3,j=1)

easy_plot_unevencolorbars2(lon=xx,lat=yy, data=hot_days_25_50-cold_days_25_50,
                  uneven_levels=uneven_levels2,name=f'Difference {size}', units = f'{cbarunits}',i=3,j=2)

size = '50-100 sqr.degs.'

easy_plot_unevencolorbars(lon=xx,lat=yy, data=hot_days_50_100,
                  uneven_levels=uneven_levels1,name=f' Hot Season {variablename} {size}',units= f'{cbarunits}',i=4,j=0)

easy_plot_unevencolorbars(lon=xx,lat=yy, data=cold_days_50_100,
                  uneven_levels=uneven_levels1,name=f' Cold Season {variablename} {size}',units= f'{cbarunits}',i=4,j=1)

easy_plot_unevencolorbars2(lon=xx,lat=yy, data=hot_days_50_100-cold_days_50_100,
                  uneven_levels=uneven_levels2,name=f'Difference {size}',units = f'{cbarunits}',i=4,j=2)

size = '100-200 sqr.degs.'
easy_plot_unevencolorbars(lon=xx,lat=yy, data=hot_days_100_200,
                  uneven_levels=uneven_levels1,name=f' Hot Season {variablename} {size}',units= f'{cbarunits}',i=5,j=0)

easy_plot_unevencolorbars(lon=xx,lat=yy, data=cold_days_100_200,
                  uneven_levels=uneven_levels1,name=f' Cold Season {variablename} {size}',units= f'{cbarunits}',i=5,j=1)

easy_plot_unevencolorbars2(lon=xx,lat=yy, data=hot_days_100_200-cold_days_100_200,
                  uneven_levels=uneven_levels2,name=f'Difference {size}',units= f'{cbarunits}',i=5,j=2)

size = '200-500 sqr.degs.'

easy_plot_unevencolorbars(lon=xx,lat=yy, data=hot_days_200_500,
                  uneven_levels=uneven_levels1,name=f' Hot Season {variablename} {size}',units= f'{cbarunits}',i=6,j=0)

easy_plot_unevencolorbars(lon=xx,lat=yy, data=cold_days_200_500,
                  uneven_levels=uneven_levels1,name=f' Cold Season {variablename} {size}',units= f'{cbarunits}',i=6,j=1)

easy_plot_unevencolorbars2(lon=xx,lat=yy, data=hot_days_200_500-cold_days_200_500,
                  uneven_levels=uneven_levels2,name=f'Difference {size}',units=f'{cbarunits}',i=6,j=2)

size = '500-700 sqr.degs.'

easy_plot_unevencolorbars(lon=xx,lat=yy, data=hot_days_500_700,
                  uneven_levels=uneven_levels1,name= f' Hot Season {variablename} {size}',units= f'{cbarunits}',i=7,j=0)

easy_plot_unevencolorbars(lon=xx,lat=yy, data=cold_days_500_700,
                  uneven_levels=uneven_levels1,name= f'Cold Season {variablename} {size}',units= f'{cbarunits}',i=7,j=1)

easy_plot_unevencolorbars2(lon=xx,lat=yy, data=hot_days_500_700-cold_days_500_700,
                  uneven_levels=uneven_levels2,name=f'Difference {size}',units= f'{cbarunits}',i=7,j=2)


plt.tight_layout()
plt.savefig(f'/home/shilpa/glory_mat_analysis/compare_sum_mhw_days_hot_cold_seasons1.png')


# no. of mhw events 

# plot
nrows=8
ncols=3

cbarunits = 'Events'
variablename = 'No. of MHW Events'

fig, axs = plt.subplots(nrows=nrows,ncols=ncols,subplot_kw={'projection': proj_crs},figsize=(22,30))

#uneven_levels = [1,2,3,4,5,6,7,8,9,10,12,14,16,17,18,19,20,21,22,24,26,28,30,32,34,36,38,40,45,50,60,70,89]
#uneven_levels = [1,2,3,4,5,6,7,8,9,10,12,14,16,17,18,19,20,21,22,24,26,28,30,32,34,36,38,40,45,50,60,70,89]
uneven_levels1 = [1,2,3,4,5,6,7,8,9,10,12,14,16,17,18,19,20,21,22,24,26,28,30,32,34,36,38,40,45,50,60,70,89]


uneven_levels2 = [-89,-60,-50,-40,-30,-28,-26,-24,-22,-20,-18,-16,-14,-12,-10,-8,-6,-4,-2,0,2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,40,50,60,89]



size = '0-1 sqr.degs.'

easy_plot_unevencolorbars(lon=xx,lat=yy, data=hot_events_0_1,
                  uneven_levels=uneven_levels1,name=f' Hot Season {variablename} {size}',units=f'{cbarunits}',i=0,j=0)

easy_plot_unevencolorbars(lon=xx,lat=yy, data=cold_events_0_1,
                  uneven_levels=uneven_levels1,name= f' Cold Season {variablename} {size}',units=f'{cbarunits}',i=0,j=1)

easy_plot_unevencolorbars2(lon=xx,lat=yy, data=hot_events_0_1-cold_events_0_1,
                  uneven_levels=uneven_levels2,name=f'Difference {size}',units=f'{cbarunits}',i=0,j=2)

size = '1-10 sqr.degs.'

easy_plot_unevencolorbars(lon=xx,lat=yy, data=hot_events_1_10,
                  uneven_levels=uneven_levels1,name= f' Hot Season {variablename} {size}',units=f'{cbarunits}',i=1,j=0)

easy_plot_unevencolorbars(lon=xx,lat=yy, data=cold_events_1_10,
                  uneven_levels=uneven_levels1,name=f' Cold Season {variablename} {size}',units=f'{cbarunits}',i=1,j=1)

easy_plot_unevencolorbars2(lon=xx,lat=yy, data=hot_events_1_10-cold_events_1_10,
                  uneven_levels=uneven_levels2,name=f'Difference {size}',units = f'{cbarunits}',i=1,j=2)

size = '10-25 sqr.degs.'

easy_plot_unevencolorbars(lon=xx,lat=yy, data=hot_events_10_25,
                  uneven_levels=uneven_levels1,name=f' Hot Season {variablename} {size}',units= f'{cbarunits}',i=2,j=0)

easy_plot_unevencolorbars(lon=xx,lat=yy, data=cold_events_10_25,
                  uneven_levels=uneven_levels1,name=f' Cold Season {variablename} {size}',units= f'{cbarunits}',i=2,j=1)

easy_plot_unevencolorbars2(lon=xx,lat=yy, data=hot_events_10_25-cold_events_10_25,
                  uneven_levels=uneven_levels2,name=f'Difference {size}',units= f'{cbarunits}',i=2,j=2)

size = '25-50 sqr.degs.'

easy_plot_unevencolorbars(lon=xx,lat=yy, data=hot_events_25_50,
                  uneven_levels=uneven_levels1,name=f' Hot Season {variablename} {size}',units= f'{cbarunits}',i=3,j=0)

easy_plot_unevencolorbars(lon=xx,lat=yy, data=cold_events_25_50,
                  uneven_levels=uneven_levels1,name=f' Cold Season {variablename} {size}',units= f'{cbarunits}',i=3,j=1)

easy_plot_unevencolorbars2(lon=xx,lat=yy, data=hot_events_25_50-cold_events_25_50,
                  uneven_levels=uneven_levels2,name=f'Difference {size}', units = f'{cbarunits}',i=3,j=2)

size = '50-100 sqr.degs.'

easy_plot_unevencolorbars(lon=xx,lat=yy, data=hot_events_50_100,
                  uneven_levels=uneven_levels1,name=f' Hot Season {variablename} {size}',units= f'{cbarunits}',i=4,j=0)

easy_plot_unevencolorbars(lon=xx,lat=yy, data=cold_events_50_100,
                  uneven_levels=uneven_levels1,name=f' Cold Season {variablename} {size}',units= f'{cbarunits}',i=4,j=1)

easy_plot_unevencolorbars2(lon=xx,lat=yy, data=hot_events_50_100-cold_events_50_100,
                  uneven_levels=uneven_levels2,name=f'Difference {size}',units = f'{cbarunits}',i=4,j=2)

size = '100-200 sqr.degs.'

easy_plot_unevencolorbars(lon=xx,lat=yy, data=hot_events_100_200,
                  uneven_levels=uneven_levels1,name=f' Hot Season {variablename} {size}',units= f'{cbarunits}',i=5,j=0)

easy_plot_unevencolorbars(lon=xx,lat=yy, data=cold_events_100_200,
                  uneven_levels=uneven_levels1,name=f' Cold Season {variablename} {size}',units= f'{cbarunits}',i=5,j=1)

easy_plot_unevencolorbars2(lon=xx,lat=yy, data=hot_events_100_200-cold_events_100_200,
                  uneven_levels=uneven_levels2,name=f'Difference {size}',units= f'{cbarunits}',i=5,j=2)

size = '200-500 sqr.degs.'

easy_plot_unevencolorbars(lon=xx,lat=yy, data=hot_events_200_500,
                  uneven_levels=uneven_levels1,name=f' Hot Season {variablename} {size}',units= f'{cbarunits}',i=6,j=0)

easy_plot_unevencolorbars(lon=xx,lat=yy, data=cold_events_200_500,
                  uneven_levels=uneven_levels1,name=f' Cold Season {variablename} {size}',units= f'{cbarunits}',i=6,j=1)

easy_plot_unevencolorbars2(lon=xx,lat=yy, data=hot_events_200_500-cold_events_200_500,
                  uneven_levels=uneven_levels2,name=f'Difference {size}',units=f'{cbarunits}',i=6,j=2)

size = '500-700 sqr.degs.'

easy_plot_unevencolorbars(lon=xx,lat=yy, data=hot_events_500_700,
                  uneven_levels=uneven_levels1,name= f' Hot Season {variablename} {size}',units= f'{cbarunits}',i=7,j=0)

easy_plot_unevencolorbars(lon=xx,lat=yy, data=cold_events_500_700,
                  uneven_levels=uneven_levels1,name= f'Cold Season {variablename} {size}',units= f'{cbarunits}',i=7,j=1)

easy_plot_unevencolorbars2(lon=xx,lat=yy, data=hot_events_500_700-cold_events_500_700,
                  uneven_levels=uneven_levels2,name=f'Difference {size}',units= f'{cbarunits}',i=7,j=2)


plt.tight_layout()
plt.savefig(f'/home/shilpa/glory_mat_analysis/compare_no_of_events_hot_cold_seasons1.png')


