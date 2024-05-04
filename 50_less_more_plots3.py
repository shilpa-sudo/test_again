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


duration = xr.open_dataset('/home/shilpa/glory_mat_analysis/means_std_median_polyareas_noaaoisst/noaaoisst_mean_std_mhw_prop/noaaoisst_1981_2023_mean_std_median_duration.nc')
mean_intensity = xr.open_dataset('/home/shilpa/glory_mat_analysis/means_std_median_polyareas_noaaoisst/noaaoisst_mean_std_mhw_prop/noaaoisst_1981_2023_mean_std_median_intensity_mean_pacific.nc')
max_intensity = xr.open_dataset('/home/shilpa/glory_mat_analysis/means_std_median_polyareas_noaaoisst/noaaoisst_mean_std_mhw_prop/noaaoisst_1981_2023_mean_std_median_intensity_max_pacific.nc')
cumulative_intensity = xr.open_dataset('/home/shilpa/glory_mat_analysis/means_std_median_polyareas_noaaoisst/noaaoisst_mean_std_mhw_prop/noaaoisst_1981_2023_mean_std_median_intensity_cumulative_pacific.nc')
rateonset = xr.open_dataset('/home/shilpa/glory_mat_analysis/means_std_median_polyareas_noaaoisst/noaaoisst_mean_std_mhw_prop/noaaoisst_1981_2023_mean_std_median_rateonset_pacific.nc')
ratedecline = xr.open_dataset('/home/shilpa/glory_mat_analysis/means_std_median_polyareas_noaaoisst/noaaoisst_mean_std_mhw_prop/noaaoisst_1981_2023_mean_std_median_ratedecline_arr.nc')
reactionwindow = xr.open_dataset('/home/shilpa/glory_mat_analysis/means_std_median_polyareas_noaaoisst/noaaoisst_mean_std_mhw_prop/noaaoisst_1981_2023_mean_std_median_reaction_window_days_pacific.nc')
copingwindow = xr.open_dataset('/home/shilpa/glory_mat_analysis/means_std_median_polyareas_noaaoisst/noaaoisst_mean_std_mhw_prop/noaaoisst_1981_2023_mean_std_median_coping_window_days_pacific.nc')
recoverywindow = xr.open_dataset('/home/shilpa/glory_mat_analysis/means_std_median_polyareas_noaaoisst/noaaoisst_mean_std_mhw_prop/noaaoisst_1981_2023_mean_std_median_recovery_window_days_pacific.nc')
days = xr.open_dataset('/home/shilpa/glory_mat_analysis/means_std_median_polyareas_noaaoisst/noaaoisst_mean_std_mhw_prop/noaaoisst_1981_2023_total_no_mhwdays.nc')
events = xr.open_dataset('/home/shilpa/glory_mat_analysis/means_std_median_polyareas_noaaoisst/noaaoisst_mean_std_mhw_prop/noaaoisst_1981_2023_total_no_mhwevents.nc')
prob = xr.open_dataset('/home/shilpa/glory_mat_analysis/means_std_median_polyareas_noaaoisst/noaaoisst_mean_std_mhw_prop/noaaoisst_1981_2023_mhw_probability_of_occurrance_sizes_seasons.nc')

longi = duration.lon
lati = duration.lat


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

    cbar = plt.colorbar(im, ax=axs[i,j],orientation="horizontal", shrink=0.75)
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

    cbar = plt.colorbar(im, ax=axs[i,j],orientation="horizontal", shrink=0.75)
    cbar.set_label(units)
    gl = axs[i,j].gridlines(draw_labels=True)
    gl.top_labels = False
    gl.right_labels = False
    axs[i,j].set_xlim(-35,29) #-35
    axs[i,j].set_ylim(-34.875,-2.375)
    axs[i,j].set_title(name, fontsize=9)



mean_cold_dur_0_700 = np.array(duration.mean_duration_pacific_0_700_cold)
mean_cold_dur_0_25 = np.array(duration.mean_duration_pacific_0_25_cold)
mean_cold_dur_25_700 = np.array(duration.mean_duration_pacific_25_700_cold)
median_cold_dur_0_700 = np.array(duration.median_duration_pacific_0_700_cold)
median_cold_dur_0_25 = np.array(duration.median_duration_pacific_0_25_cold)
median_cold_dur_25_700 = np.array(duration.median_duration_pacific_25_700_cold)


mean_hot_dur_0_700 = np.array(duration.mean_duration_pacific_0_700_hot)
mean_hot_dur_0_25 = np.array(duration.mean_duration_pacific_0_25_hot)
mean_hot_dur_25_700 = np.array(duration.mean_duration_pacific_25_700_hot)
median_hot_dur_0_700 = np.array(duration.median_duration_pacific_0_700_hot)
median_hot_dur_0_25 = np.array(duration.median_duration_pacific_0_25_hot)
median_hot_dur_25_700 = np.array(duration.median_duration_pacific_25_700_hot)

mean_cold_maxintensity_0_700 = np.array(max_intensity.mean_intensity_max_pacific_0_700_cold)
mean_cold_maxintensity_0_25 = np.array(max_intensity.mean_intensity_max_pacific_0_25_cold)
mean_cold_maxintensity_25_700 = np.array(max_intensity.mean_intensity_max_pacific_25_700_cold)
median_cold_maxintensity_0_700 = np.array(max_intensity.median_intensity_max_pacific_0_700_cold)
median_cold_maxintensity_0_25 = np.array(max_intensity.median_intensity_max_pacific_0_25_cold)
median_cold_maxintensity_25_700 = np.array(max_intensity.median_intensity_max_pacific_25_700_cold)


mean_hot_maxintensity_0_700 = np.array(max_intensity.mean_intensity_max_pacific_0_700_hot)
mean_hot_maxintensity_0_25 = np.array(max_intensity.mean_intensity_max_pacific_0_25_hot)
mean_hot_maxintensity_25_700 = np.array(max_intensity.mean_intensity_max_pacific_25_700_hot)
median_hot_maxintensity_0_700 = np.array(max_intensity.median_intensity_max_pacific_0_700_hot)
median_hot_maxintensity_0_25 = np.array(max_intensity.median_intensity_max_pacific_0_25_hot)
median_hot_maxintensity_25_700 = np.array(max_intensity.median_intensity_max_pacific_25_700_hot)

mean_cold_meanintensity_0_700 = np.array(mean_intensity.mean_intensity_mean_pacific_0_700_cold)
mean_cold_meanintensity_0_25 = np.array(mean_intensity.mean_intensity_mean_pacific_0_25_cold)
mean_cold_meanintensity_25_700 = np.array(mean_intensity.mean_intensity_mean_pacific_25_700_cold)
median_cold_meanintensity_0_700 = np.array(mean_intensity.median_intensity_mean_pacific_0_700_cold)
median_cold_meanintensity_0_25 = np.array(mean_intensity.median_intensity_mean_pacific_0_25_cold)
median_cold_meanintensity_25_700 = np.array(mean_intensity.median_intensity_mean_pacific_25_700_cold)


mean_hot_meanintensity_0_700 = np.array(mean_intensity.mean_intensity_mean_pacific_0_700_hot)
mean_hot_meanintensity_0_25 = np.array(mean_intensity.mean_intensity_mean_pacific_0_25_hot)
mean_hot_meanintensity_25_700 = np.array(mean_intensity.mean_intensity_mean_pacific_25_700_hot)
median_hot_meanintensity_0_700 = np.array(mean_intensity.median_intensity_mean_pacific_0_700_hot)
median_hot_meanintensity_0_25 = np.array(mean_intensity.median_intensity_mean_pacific_0_25_hot)
median_hot_meanintensity_25_700 = np.array(mean_intensity.median_intensity_mean_pacific_25_700_hot)

mean_cold_cumintensity_0_700 = np.array(cumulative_intensity.mean_intensity_cumulative_pacific_0_700_cold)
mean_cold_cumintensity_0_25 = np.array(cumulative_intensity.mean_intensity_cumulative_pacific_0_25_cold)
mean_cold_cumintensity_25_700 = np.array(cumulative_intensity.mean_intensity_cumulative_pacific_25_700_cold)
median_cold_cumintensity_0_700 = np.array(cumulative_intensity.median_intensity_cumulative_pacific_0_700_cold)
median_cold_cumintensity_0_25 = np.array(cumulative_intensity.median_intensity_cumulative_pacific_0_25_cold)
median_cold_cumintensity_25_700 = np.array(cumulative_intensity.median_intensity_cumulative_pacific_25_700_cold)


mean_hot_cumintensity_0_700 = np.array(cumulative_intensity.mean_intensity_cumulative_pacific_0_700_hot)
mean_hot_cumintensity_0_25 = np.array(cumulative_intensity.mean_intensity_cumulative_pacific_0_25_hot)
mean_hot_cumintensity_25_700 = np.array(cumulative_intensity.mean_intensity_cumulative_pacific_25_700_hot)
median_hot_cumintensity_0_700 = np.array(cumulative_intensity.median_intensity_cumulative_pacific_0_700_hot)
median_hot_cumintensity_0_25 = np.array(cumulative_intensity.median_intensity_cumulative_pacific_0_25_hot)
median_hot_cumintensity_25_700 = np.array(cumulative_intensity.median_intensity_cumulative_pacific_25_700_hot)

mean_cold_rateonset_0_700 = np.array(rateonset.mean_rateonset_pacific_0_700_cold)
mean_cold_rateonset_0_25 = np.array(rateonset.mean_rateonset_pacific_0_25_cold)
mean_cold_rateonset_25_700 = np.array(rateonset.mean_rateonset_pacific_25_700_cold)
median_cold_rateonset_0_700 = np.array(rateonset.median_rateonset_pacific_0_700_cold)
median_cold_rateonset_0_25 = np.array(rateonset.median_rateonset_pacific_0_25_cold)
median_cold_rateonset_25_700 = np.array(rateonset.median_rateonset_pacific_25_700_cold)


mean_hot_rateonset_0_700 = np.array(rateonset.mean_rateonset_pacific_0_700_hot)
mean_hot_rateonset_0_25 = np.array(rateonset.mean_rateonset_pacific_0_25_hot)
mean_hot_rateonset_25_700 = np.array(rateonset.mean_rateonset_pacific_25_700_hot)
median_hot_rateonset_0_700 = np.array(rateonset.median_rateonset_pacific_0_700_hot)
median_hot_rateonset_0_25 = np.array(rateonset.median_rateonset_pacific_0_25_hot)
median_hot_rateonset_25_700 = np.array(rateonset.median_rateonset_pacific_25_700_hot)

mean_cold_ratedecline_0_700 = np.array(ratedecline.mean_ratedecline_arr_0_700_cold)
mean_cold_ratedecline_0_25 = np.array(ratedecline.mean_ratedecline_arr_0_25_cold)
mean_cold_ratedecline_25_700 = np.array(ratedecline.mean_ratedecline_arr_25_700_cold)
median_cold_ratedecline_0_700 = np.array(ratedecline.median_ratedecline_arr_0_700_cold)
median_cold_ratedecline_0_25 = np.array(ratedecline.median_ratedecline_arr_0_25_cold)
median_cold_ratedecline_25_700 = np.array(ratedecline.median_ratedecline_arr_25_700_cold)


mean_hot_ratedecline_0_700 = np.array(ratedecline.mean_ratedecline_arr_0_700_hot)
mean_hot_ratedecline_0_25 = np.array(ratedecline.mean_ratedecline_arr_0_25_hot)
mean_hot_ratedecline_25_700 = np.array(ratedecline.mean_ratedecline_arr_25_700_hot)
median_hot_ratedecline_0_700 = np.array(ratedecline.median_ratedecline_arr_0_700_hot)
median_hot_ratedecline_0_25 = np.array(ratedecline.median_ratedecline_arr_0_25_hot)
median_hot_ratedecline_25_700 = np.array(ratedecline.median_ratedecline_arr_25_700_hot)


mean_cold_reaction_window_0_700 = np.array(reactionwindow.mean_reaction_window_days_pacific_0_700_cold)
mean_cold_reaction_window_0_25 = np.array(reactionwindow.mean_reaction_window_days_pacific_0_25_cold)
mean_cold_reaction_window_25_700 = np.array(reactionwindow.mean_reaction_window_days_pacific_25_700_cold)
median_cold_reaction_window_0_700 = np.array(reactionwindow.median_reaction_window_days_pacific_0_700_cold)
median_cold_reaction_window_0_25 = np.array(reactionwindow.median_reaction_window_days_pacific_0_25_cold)
median_cold_reaction_window_25_700 = np.array(reactionwindow.median_reaction_window_days_pacific_25_700_cold)


mean_hot_reaction_window_0_700 = np.array(reactionwindow.mean_reaction_window_days_pacific_0_700_hot)
mean_hot_reaction_window_0_25 = np.array(reactionwindow.mean_reaction_window_days_pacific_0_25_hot)
mean_hot_reaction_window_25_700 = np.array(reactionwindow.mean_reaction_window_days_pacific_25_700_hot)
median_hot_reaction_window_0_700 = np.array(reactionwindow.median_reaction_window_days_pacific_0_700_hot)
median_hot_reaction_window_0_25 = np.array(reactionwindow.median_reaction_window_days_pacific_0_25_hot)
median_hot_reaction_window_25_700 = np.array(reactionwindow.median_reaction_window_days_pacific_25_700_hot)

mean_cold_recovery_window_0_700 = np.array(recoverywindow.mean_recovery_window_days_pacific_0_700_cold)
mean_cold_recovery_window_0_25 = np.array(recoverywindow.mean_recovery_window_days_pacific_0_25_cold)
mean_cold_recovery_window_25_700 = np.array(recoverywindow.mean_recovery_window_days_pacific_25_700_cold)
median_cold_recovery_window_0_700 = np.array(recoverywindow.median_recovery_window_days_pacific_0_700_cold)
median_cold_recovery_window_0_25 = np.array(recoverywindow.median_recovery_window_days_pacific_0_25_cold)
median_cold_recovery_window_25_700 = np.array(recoverywindow.median_recovery_window_days_pacific_25_700_cold)


mean_hot_recovery_window_0_700 = np.array(recoverywindow.mean_recovery_window_days_pacific_0_700_hot)
mean_hot_recovery_window_0_25 = np.array(recoverywindow.mean_recovery_window_days_pacific_0_25_hot)
mean_hot_recovery_window_25_700 = np.array(recoverywindow.mean_recovery_window_days_pacific_25_700_hot)
median_hot_recovery_window_0_700 = np.array(recoverywindow.median_recovery_window_days_pacific_0_700_hot)
median_hot_recovery_window_0_25 = np.array(recoverywindow.median_recovery_window_days_pacific_0_25_hot)
median_hot_recovery_window_25_700 = np.array(recoverywindow.median_recovery_window_days_pacific_25_700_hot)

mean_cold_coping_window_0_700 = np.array(copingwindow.mean_coping_window_days_pacific_0_700_cold)
mean_cold_coping_window_0_25 = np.array(copingwindow.mean_coping_window_days_pacific_0_25_cold)
mean_cold_coping_window_25_700 = np.array(copingwindow.mean_coping_window_days_pacific_25_700_cold)
median_cold_coping_window_0_700 = np.array(copingwindow.median_coping_window_days_pacific_0_700_cold)
median_cold_coping_window_0_25 = np.array(copingwindow.median_coping_window_days_pacific_0_25_cold)
median_cold_coping_window_25_700 = np.array(copingwindow.median_coping_window_days_pacific_25_700_cold)


mean_hot_coping_window_0_700 = np.array(copingwindow.mean_coping_window_days_pacific_0_700_hot)
mean_hot_coping_window_0_25 = np.array(copingwindow.mean_coping_window_days_pacific_0_25_hot)
mean_hot_coping_window_25_700 = np.array(copingwindow.mean_coping_window_days_pacific_25_700_hot)
median_hot_coping_window_0_700 = np.array(copingwindow.median_coping_window_days_pacific_0_700_hot)
median_hot_coping_window_0_25 = np.array(copingwindow.median_coping_window_days_pacific_0_25_hot)
median_hot_coping_window_25_700 = np.array(copingwindow.median_coping_window_days_pacific_25_700_hot)



hot_days_0_700 = np.array(days.total_days_mhw_0_700_hot)
hot_days_0_25 = np.array(days.total_days_mhw_0_25_hot)
hot_days_25_700 = np.array(days.total_days_mhw_25_700_hot)


cold_days_0_700 = np.array(days.total_days_mhw_0_700_cold)
cold_days_0_25 = np.array(days.total_days_mhw_0_25_cold)
cold_days_25_700 = np.array(days.total_days_mhw_25_700_cold)


hot_events_0_700 = np.array(events.total_events_mhw_0_700_hot)
hot_events_0_25 = np.array(events.total_events_mhw_0_25_hot)
hot_events_25_700 = np.array(events.total_events_mhw_25_700_hot)

cold_events_0_700 = np.array(events.total_events_mhw_0_700_cold)
cold_events_0_25 = np.array(events.total_events_mhw_0_25_cold)
cold_events_25_700 = np.array(events.total_events_mhw_25_700_cold)



hot_prob_0_700_alldays = np.array(prob.hot_all_time_probability_of_occurance_0_700)
hot_prob_0_25_alldays = np.array(prob.hot_all_time_probability_of_occurance_0_25)
hot_prob_25_700_alldays = np.array(prob.hot_all_time_probability_of_occurance_25_700)

cold_prob_0_700_alldays = np.array(prob.cold_all_time_probability_of_occurance_0_700)
cold_prob_0_25_alldays = np.array(prob.cold_all_time_probability_of_occurance_0_25)
cold_prob_25_700_alldays = np.array(prob.cold_all_time_probability_of_occurance_25_700)

hot_prob_0_700_mhwdays = np.array(prob.hot_mhw_only_probability_of_occurance_0_700)
hot_prob_0_25_mhwdays = np.array(prob.hot_mhw_only_probability_of_occurance_0_25)
hot_prob_25_700_mhwdays = np.array(prob.hot_mhw_only_probability_of_occurance_25_700)

cold_prob_0_700_mhwdays = np.array(prob.cold_mhw_only_probability_of_occurance_0_700)
cold_prob_0_25_mhwdays = np.array(prob.cold_mhw_only_probability_of_occurance_0_25)
cold_prob_25_700_mhwdays = np.array(prob.cold_mhw_only_probability_of_occurance_25_700)

mean_dur_0_700 = np.array(duration.mean_duration_pacific_0_700)
mean_dur_0_25 = np.array(duration.mean_duration_pacific_0_25)
mean_dur_25_700 = np.array(duration.mean_duration_pacific_25_700)

median_dur_0_700 = np.array(duration.median_duration_pacific_0_700)
median_dur_0_25 = np.array(duration.median_duration_pacific_0_25)
median_dur_25_700 = np.array(duration.median_duration_pacific_25_700)

mean_maxintensity_0_700 = np.array(max_intensity.mean_intensity_max_pacific_0_700)
mean_maxintensity_0_25 = np.array(max_intensity.mean_intensity_max_pacific_0_25)
mean_maxintensity_25_700 = np.array(max_intensity.mean_intensity_max_pacific_25_700)

median_maxintensity_0_700 = np.array(max_intensity.median_intensity_max_pacific_0_700)
median_maxintensity_0_25 = np.array(max_intensity.median_intensity_max_pacific_0_25)
median_maxintensity_25_700 = np.array(max_intensity.median_intensity_max_pacific_25_700)

mean_meanintensity_0_700 = np.array(mean_intensity.mean_intensity_mean_pacific_0_700)
mean_meanintensity_0_25 = np.array(mean_intensity.mean_intensity_mean_pacific_0_25)
mean_meanintensity_25_700 = np.array(mean_intensity.mean_intensity_mean_pacific_25_700)

median_meanintensity_0_700 = np.array(mean_intensity.median_intensity_mean_pacific_0_700)
median_meanintensity_0_25 = np.array(mean_intensity.median_intensity_mean_pacific_0_25)
median_meanintensity_25_700 = np.array(mean_intensity.median_intensity_mean_pacific_25_700)

mean_cumintensity_0_700 = np.array(cumulative_intensity.mean_intensity_cumulative_pacific_0_700)
mean_cumintensity_0_25 = np.array(cumulative_intensity.mean_intensity_cumulative_pacific_0_25)
mean_cumintensity_25_700 = np.array(cumulative_intensity.mean_intensity_cumulative_pacific_25_700)

median_cumintensity_0_700 = np.array(cumulative_intensity.median_intensity_cumulative_pacific_0_700)
median_cumintensity_0_25 = np.array(cumulative_intensity.median_intensity_cumulative_pacific_0_25)
median_cumintensity_25_700 = np.array(cumulative_intensity.median_intensity_cumulative_pacific_25_700)


mean_rateonset_0_700 = np.array(rateonset.mean_rateonset_pacific_0_700)
mean_rateonset_0_25 = np.array(rateonset.mean_rateonset_pacific_0_25)
mean_rateonset_25_700 = np.array(rateonset.mean_rateonset_pacific_25_700)

median_rateonset_0_700 = np.array(rateonset.median_rateonset_pacific_0_700)
median_rateonset_0_25 = np.array(rateonset.median_rateonset_pacific_0_25)
median_rateonset_25_700 = np.array(rateonset.median_rateonset_pacific_25_700)

mean_ratedecline_0_700 = np.array(ratedecline.mean_ratedecline_arr_0_700)
mean_ratedecline_0_25 = np.array(ratedecline.mean_ratedecline_arr_0_25)
mean_ratedecline_25_700 = np.array(ratedecline.mean_ratedecline_arr_25_700)

median_ratedecline_0_700 = np.array(ratedecline.median_ratedecline_arr_0_700)
median_ratedecline_0_25 = np.array(ratedecline.median_ratedecline_arr_0_25)
median_ratedecline_25_700 = np.array(ratedecline.median_ratedecline_arr_25_700)


mean_reaction_window_0_700 = np.array(reactionwindow.mean_reaction_window_days_pacific_0_700)
mean_reaction_window_0_25 = np.array(reactionwindow.mean_reaction_window_days_pacific_0_25)
mean_reaction_window_25_700 = np.array(reactionwindow.mean_reaction_window_days_pacific_25_700)

median_reaction_window_0_700 = np.array(reactionwindow.median_reaction_window_days_pacific_0_700)
median_reaction_window_0_25 = np.array(reactionwindow.median_reaction_window_days_pacific_0_25)
median_reaction_window_25_700 = np.array(reactionwindow.median_reaction_window_days_pacific_25_700)

mean_recovery_window_0_700 = np.array(recoverywindow.mean_recovery_window_days_pacific_0_700)
mean_recovery_window_0_25 = np.array(recoverywindow.mean_recovery_window_days_pacific_0_25)
mean_recovery_window_25_700 = np.array(recoverywindow.mean_recovery_window_days_pacific_25_700)

median_recovery_window_0_700 = np.array(recoverywindow.median_recovery_window_days_pacific_0_700)
median_recovery_window_0_25 = np.array(recoverywindow.median_recovery_window_days_pacific_0_25)
median_recovery_window_25_700 = np.array(recoverywindow.median_recovery_window_days_pacific_25_700)

mean_coping_window_0_700 = np.array(copingwindow.mean_coping_window_days_pacific_0_700)
mean_coping_window_0_25 = np.array(copingwindow.mean_coping_window_days_pacific_0_25)
mean_coping_window_25_700 = np.array(copingwindow.mean_coping_window_days_pacific_25_700)

median_coping_window_0_700 = np.array(copingwindow.median_coping_window_days_pacific_0_700)
median_coping_window_0_25 = np.array(copingwindow.median_coping_window_days_pacific_0_25)
median_coping_window_25_700 = np.array(copingwindow.median_coping_window_days_pacific_25_700)


events_0_700 = np.array(events.total_events_mhw_0_700)
events_0_25 = np.array(events.total_events_mhw_0_25)
events_25_700 = np.array(events.total_events_mhw_25_700)

days_0_700 = np.array(days.total_days_mhw_0_700)
days_0_25 = np.array(days.total_days_mhw_0_25)
days_25_700 = np.array(days.total_days_mhw_25_700)

prob_0_700_alldays = np.array(prob.all_time_probability_of_occurance_0_700)
prob_0_25_alldays = np.array(prob.all_time_probability_of_occurance_0_25)
prob_25_700_alldays = np.array(prob.all_time_probability_of_occurance_25_700)

prob_0_700_mhwdays = np.array(prob.mhw_only_probability_of_occurance_0_700)
prob_0_25_mhwdays = np.array(prob.mhw_only_probability_of_occurance_0_25)
prob_25_700_mhwdays = np.array(prob.mhw_only_probability_of_occurance_25_700)




# plot
nrows=3
ncols=4

cbarunits = 'Days'
variablename = 'Mean Duration'

fig, axs = plt.subplots(nrows=nrows,ncols=ncols,subplot_kw={'projection': proj_crs},figsize=(22,15))

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

size = '0-25 sqr.degs.'

easy_plot_unevencolorbars(lon=xx,lat=yy, data=mean_dur_0_25,
                  uneven_levels=uneven_levels1,name=f' {variablename} {size}',units=f'{cbarunits}',i=1,j=0)

easy_plot_unevencolorbars(lon=xx,lat=yy, data=mean_hot_dur_0_25,
                  uneven_levels=uneven_levels1,name=f' Hot Season {variablename} {size}',units=f'{cbarunits}',i=1,j=1)

easy_plot_unevencolorbars(lon=xx,lat=yy, data=mean_cold_dur_0_25,
                  uneven_levels=uneven_levels1,name= f' Cold Season {variablename} {size}',units=f'{cbarunits}',i=1,j=2)

easy_plot_unevencolorbars2(lon=xx,lat=yy, data=mean_hot_dur_0_25-mean_cold_dur_0_25,
                  uneven_levels=uneven_levels2,name=f'Difference {size}',units=f'{cbarunits}',i=1,j=3)

size = '25-700 sqr.degs.'

easy_plot_unevencolorbars(lon=xx,lat=yy, data=mean_dur_25_700,
                  uneven_levels=uneven_levels1,name=f' {variablename} {size}',units=f'{cbarunits}',i=2,j=0)

easy_plot_unevencolorbars(lon=xx,lat=yy, data=mean_hot_dur_25_700,
                  uneven_levels=uneven_levels1,name=f' Hot Season {variablename} {size}',units=f'{cbarunits}',i=2,j=1)

easy_plot_unevencolorbars(lon=xx,lat=yy, data=mean_cold_dur_25_700,
                  uneven_levels=uneven_levels1,name= f' Cold Season {variablename} {size}',units=f'{cbarunits}',i=2,j=2)

easy_plot_unevencolorbars2(lon=xx,lat=yy, data=mean_hot_dur_25_700-mean_cold_dur_25_700,
                  uneven_levels=uneven_levels2,name=f'Difference {size}',units=f'{cbarunits}',i=2,j=3)

plt.tight_layout()
plt.savefig(f'/home/shilpa/glory_mat_analysis/compare_mean_duration_hot_cold_seasons2.png')


# plot
nrows=3
ncols=4

cbarunits = 'Days'
variablename = 'Median Duration'

fig, axs = plt.subplots(nrows=nrows,ncols=ncols,subplot_kw={'projection': proj_crs},figsize=(22,15))

uneven_levels1 = [5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,40,45,50,55,60,65,70,80,90,100,150,200,250,300,350]

uneven_levels2 = [-350,-100,-50,-10,-5,0,5,10,50,100,350]

size = '0-700 sqr.degs.'

easy_plot_unevencolorbars(lon=xx,lat=yy, data=median_dur_0_700,
                  uneven_levels=uneven_levels1,name=f' {variablename} {size}',units=f'{cbarunits}',i=0,j=0)

easy_plot_unevencolorbars(lon=xx,lat=yy, data=median_hot_dur_0_700,
                  uneven_levels=uneven_levels1,name=f' Hot Season {variablename} {size}',units=f'{cbarunits}',i=0,j=1)

easy_plot_unevencolorbars(lon=xx,lat=yy, data=median_cold_dur_0_700,
                  uneven_levels=uneven_levels1,name= f' Cold Season {variablename} {size}',units=f'{cbarunits}',i=0,j=2)

easy_plot_unevencolorbars2(lon=xx,lat=yy, data=median_hot_dur_0_700-median_cold_dur_0_700,
                  uneven_levels=uneven_levels2,name=f'Difference {size}',units=f'{cbarunits}',i=0,j=3)

size = '0-25 sqr.degs.'

easy_plot_unevencolorbars(lon=xx,lat=yy, data=median_dur_0_25,
                  uneven_levels=uneven_levels1,name=f' {variablename} {size}',units=f'{cbarunits}',i=1,j=0)

easy_plot_unevencolorbars(lon=xx,lat=yy, data=median_hot_dur_0_25,
                  uneven_levels=uneven_levels1,name=f' Hot Season {variablename} {size}',units=f'{cbarunits}',i=1,j=1)

easy_plot_unevencolorbars(lon=xx,lat=yy, data=median_cold_dur_0_25,
                  uneven_levels=uneven_levels1,name= f' Cold Season {variablename} {size}',units=f'{cbarunits}',i=1,j=2)

easy_plot_unevencolorbars2(lon=xx,lat=yy, data=median_hot_dur_0_25-median_cold_dur_0_25,
                  uneven_levels=uneven_levels2,name=f'Difference {size}',units=f'{cbarunits}',i=1,j=3)

size = '25-700 sqr.degs.'

easy_plot_unevencolorbars(lon=xx,lat=yy, data=median_dur_25_700,
                  uneven_levels=uneven_levels1,name=f' {variablename} {size}',units=f'{cbarunits}',i=2,j=0)

easy_plot_unevencolorbars(lon=xx,lat=yy, data=median_hot_dur_25_700,
                  uneven_levels=uneven_levels1,name=f' Hot Season {variablename} {size}',units=f'{cbarunits}',i=2,j=1)

easy_plot_unevencolorbars(lon=xx,lat=yy, data=median_cold_dur_25_700,
                  uneven_levels=uneven_levels1,name= f' Cold Season {variablename} {size}',units=f'{cbarunits}',i=2,j=2)

easy_plot_unevencolorbars2(lon=xx,lat=yy, data=median_hot_dur_25_700-median_cold_dur_25_700,
                  uneven_levels=uneven_levels2,name=f'Difference {size}',units=f'{cbarunits}',i=2,j=3)

plt.tight_layout()
plt.savefig(f'/home/shilpa/glory_mat_analysis/compare_median_duration_hot_cold_seasons2.png')

# plot
nrows=3
ncols=4

cbarunits = '°C'
variablename = 'Mean Max Intensity'

fig, axs = plt.subplots(nrows=nrows,ncols=ncols,subplot_kw={'projection': proj_crs},figsize=(22,15))

uneven_levels1 = [0.5,0.6,0.7,0.8,0.9,1,1.2,1.3,1.4,1.5,1.6,1.8,2,2.2,2.4,2.6,2.8,3,3.5]

uneven_levels2 = [-3.5,-2,-1.25,-1,-0.8,-0.6,-0.4,-0.2,0,0.2,0.4,0.6,0.8,1,1.25,2,3.5]

size = '0-700 sqr.degs.'

easy_plot_unevencolorbars(lon=xx,lat=yy, data=mean_maxintensity_0_700,
                  uneven_levels=uneven_levels1,name=f' {variablename} {size}',units=f'{cbarunits}',i=0,j=0)

easy_plot_unevencolorbars(lon=xx,lat=yy, data=mean_hot_maxintensity_0_700,
                  uneven_levels=uneven_levels1,name=f' Hot Season {variablename} {size}',units=f'{cbarunits}',i=0,j=1)

easy_plot_unevencolorbars(lon=xx,lat=yy, data=mean_cold_maxintensity_0_700,
                  uneven_levels=uneven_levels1,name= f' Cold Season {variablename} {size}',units=f'{cbarunits}',i=0,j=2)

easy_plot_unevencolorbars2(lon=xx,lat=yy, data=mean_hot_maxintensity_0_700-mean_cold_maxintensity_0_700,
                  uneven_levels=uneven_levels2,name=f'Difference {size}',units=f'{cbarunits}',i=0,j=3)

size = '0-25 sqr.degs.'

easy_plot_unevencolorbars(lon=xx,lat=yy, data=mean_maxintensity_0_25,
                  uneven_levels=uneven_levels1,name=f' {variablename} {size}',units=f'{cbarunits}',i=1,j=0)

easy_plot_unevencolorbars(lon=xx,lat=yy, data=mean_hot_maxintensity_0_25,
                  uneven_levels=uneven_levels1,name=f' Hot Season {variablename} {size}',units=f'{cbarunits}',i=1,j=1)

easy_plot_unevencolorbars(lon=xx,lat=yy, data=mean_cold_maxintensity_0_25,
                  uneven_levels=uneven_levels1,name= f' Cold Season {variablename} {size}',units=f'{cbarunits}',i=1,j=2)

easy_plot_unevencolorbars2(lon=xx,lat=yy, data=mean_hot_maxintensity_0_25-mean_cold_maxintensity_0_25,
                  uneven_levels=uneven_levels2,name=f'Difference {size}',units=f'{cbarunits}',i=1,j=3)

size = '25-700 sqr.degs.'

easy_plot_unevencolorbars(lon=xx,lat=yy, data=mean_maxintensity_25_700,
                  uneven_levels=uneven_levels1,name=f' {variablename} {size}',units=f'{cbarunits}',i=2,j=0)

easy_plot_unevencolorbars(lon=xx,lat=yy, data=mean_hot_maxintensity_25_700,
                  uneven_levels=uneven_levels1,name=f' Hot Season {variablename} {size}',units=f'{cbarunits}',i=2,j=1)

easy_plot_unevencolorbars(lon=xx,lat=yy, data=mean_cold_maxintensity_25_700,
                  uneven_levels=uneven_levels1,name= f' Cold Season {variablename} {size}',units=f'{cbarunits}',i=2,j=2)

easy_plot_unevencolorbars2(lon=xx,lat=yy, data=mean_hot_maxintensity_25_700-mean_cold_maxintensity_25_700,
                  uneven_levels=uneven_levels2,name=f'Difference {size}',units=f'{cbarunits}',i=2,j=3)

plt.tight_layout()
plt.savefig(f'/home/shilpa/glory_mat_analysis/compare_mean_maxintensity_hot_cold_seasons2.png')


# plot
nrows=3
ncols=4

cbarunits = '°C'
variablename = 'Median Max Intensity'

fig, axs = plt.subplots(nrows=nrows,ncols=ncols,subplot_kw={'projection': proj_crs},figsize=(22,15))

uneven_levels1 = [0.5,0.6,0.7,0.8,0.9,1,1.1,1.2,1.3,1.4,1.5,1.6,1.8,2,2.2,2.4,2.6,2.8,3,3.5]

uneven_levels2 = [-3.5,-2,-1.25,-1,-0.8,-0.6,-0.4,-0.2,0,0.2,0.4,0.6,0.8,1,1.25,2,3.5]

size = '0-700 sqr.degs.'

easy_plot_unevencolorbars(lon=xx,lat=yy, data=median_maxintensity_0_700,
                  uneven_levels=uneven_levels1,name=f' {variablename} {size}',units=f'{cbarunits}',i=0,j=0)

easy_plot_unevencolorbars(lon=xx,lat=yy, data=median_hot_maxintensity_0_700,
                  uneven_levels=uneven_levels1,name=f' Hot Season {variablename} {size}',units=f'{cbarunits}',i=0,j=1)

easy_plot_unevencolorbars(lon=xx,lat=yy, data=median_cold_maxintensity_0_700,
                  uneven_levels=uneven_levels1,name= f' Cold Season {variablename} {size}',units=f'{cbarunits}',i=0,j=2)

easy_plot_unevencolorbars2(lon=xx,lat=yy, data=median_hot_maxintensity_0_700-median_cold_maxintensity_0_700,
                  uneven_levels=uneven_levels2,name=f'Difference {size}',units=f'{cbarunits}',i=0,j=3)

size = '0-25 sqr.degs.'

easy_plot_unevencolorbars(lon=xx,lat=yy, data=median_maxintensity_0_25,
                  uneven_levels=uneven_levels1,name=f' {variablename} {size}',units=f'{cbarunits}',i=1,j=0)

easy_plot_unevencolorbars(lon=xx,lat=yy, data=median_hot_maxintensity_0_25,
                  uneven_levels=uneven_levels1,name=f' Hot Season {variablename} {size}',units=f'{cbarunits}',i=1,j=1)

easy_plot_unevencolorbars(lon=xx,lat=yy, data=median_cold_maxintensity_0_25,
                  uneven_levels=uneven_levels1,name= f' Cold Season {variablename} {size}',units=f'{cbarunits}',i=1,j=2)

easy_plot_unevencolorbars2(lon=xx,lat=yy, data=median_hot_maxintensity_0_25-median_cold_maxintensity_0_25,
                  uneven_levels=uneven_levels2,name=f'Difference {size}',units=f'{cbarunits}',i=1,j=3)

size = '25-700 sqr.degs.'

easy_plot_unevencolorbars(lon=xx,lat=yy, data=median_maxintensity_25_700,
                  uneven_levels=uneven_levels1,name=f' {variablename} {size}',units=f'{cbarunits}',i=2,j=0)

easy_plot_unevencolorbars(lon=xx,lat=yy, data=median_hot_maxintensity_25_700,
                  uneven_levels=uneven_levels1,name=f' Hot Season {variablename} {size}',units=f'{cbarunits}',i=2,j=1)

easy_plot_unevencolorbars(lon=xx,lat=yy, data=median_cold_maxintensity_25_700,
                  uneven_levels=uneven_levels1,name= f' Cold Season {variablename} {size}',units=f'{cbarunits}',i=2,j=2)

easy_plot_unevencolorbars2(lon=xx,lat=yy, data=median_hot_maxintensity_25_700-median_cold_maxintensity_25_700,
                  uneven_levels=uneven_levels2,name=f'Difference {size}',units=f'{cbarunits}',i=2,j=3)

plt.tight_layout()
plt.savefig(f'/home/shilpa/glory_mat_analysis/compare_median_maxintensity_hot_cold_seasons2.png')

# no. of mhw days 
print('max days =',np.nanmax(days_0_700))
# plot
nrows=3
ncols=4

cbarunits = 'Days'
variablename = 'No. of MHW Days'

fig, axs = plt.subplots(nrows=nrows,ncols=ncols,subplot_kw={'projection': proj_crs},figsize=(22,15))

uneven_levels1 = [1,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100,105,110,115,120,125,150,175,200,225,250,275,300,350,400,450,500,550,600,650,700,900,1000,1200,1400,1600,1800,2000,2500,3000]

uneven_levels2 = [-1000,-700,-500,-200,-100,-90,-80,-70,-60,-50,-40,-30,-20,-10,0,10,20,30,40,50,60,70,80,90,100,200,500,700,1000]
size = '0-700 sqr.degs.'

easy_plot_unevencolorbars(lon=xx,lat=yy, data=days_0_700,
                  uneven_levels=uneven_levels1,name=f' {variablename} {size}',units=f'{cbarunits}',i=0,j=0)

easy_plot_unevencolorbars(lon=xx,lat=yy, data=hot_days_0_700,
                  uneven_levels=uneven_levels1,name=f' Hot Season {variablename} {size}',units=f'{cbarunits}',i=0,j=1)

easy_plot_unevencolorbars(lon=xx,lat=yy, data=cold_days_0_700,
                  uneven_levels=uneven_levels1,name= f' Cold Season {variablename} {size}',units=f'{cbarunits}',i=0,j=2)

easy_plot_unevencolorbars2(lon=xx,lat=yy, data=hot_days_0_700-cold_days_0_700,
                  uneven_levels=uneven_levels2,name=f'Difference {size}',units=f'{cbarunits}',i=0,j=3)

size = '0-25 sqr.degs.'

easy_plot_unevencolorbars(lon=xx,lat=yy, data=days_0_25,
                  uneven_levels=uneven_levels1,name=f' {variablename} {size}',units=f'{cbarunits}',i=1,j=0)

easy_plot_unevencolorbars(lon=xx,lat=yy, data=hot_days_0_25,
                  uneven_levels=uneven_levels1,name=f' Hot Season {variablename} {size}',units=f'{cbarunits}',i=1,j=1)

easy_plot_unevencolorbars(lon=xx,lat=yy, data=cold_days_0_25,
                  uneven_levels=uneven_levels1,name= f' Cold Season {variablename} {size}',units=f'{cbarunits}',i=1,j=2)

easy_plot_unevencolorbars2(lon=xx,lat=yy, data=hot_days_0_25-cold_days_0_25,
                  uneven_levels=uneven_levels2,name=f'Difference {size}',units=f'{cbarunits}',i=1,j=3)

size = '25-700 sqr.degs.'

easy_plot_unevencolorbars(lon=xx,lat=yy, data=days_25_700,
                  uneven_levels=uneven_levels1,name=f' {variablename} {size}',units=f'{cbarunits}',i=2,j=0)

easy_plot_unevencolorbars(lon=xx,lat=yy, data=hot_days_25_700,
                  uneven_levels=uneven_levels1,name=f' Hot Season {variablename} {size}',units=f'{cbarunits}',i=2,j=1)

easy_plot_unevencolorbars(lon=xx,lat=yy, data=cold_days_25_700,
                  uneven_levels=uneven_levels1,name= f' Cold Season {variablename} {size}',units=f'{cbarunits}',i=2,j=2)

easy_plot_unevencolorbars2(lon=xx,lat=yy, data=hot_days_25_700-cold_days_25_700,
                  uneven_levels=uneven_levels2,name=f'Difference {size}',units=f'{cbarunits}',i=2,j=3)

plt.tight_layout()
plt.savefig(f'/home/shilpa/glory_mat_analysis/compare_sum_mhw_days_hot_cold_seasons2.png')


# no. of mhw events 

# plot
nrows=3
ncols=4

cbarunits = 'Events'
variablename = 'No. of MHW Events'

fig, axs = plt.subplots(nrows=nrows,ncols=ncols,subplot_kw={'projection': proj_crs},figsize=(22,15))

#uneven_levels = [1,2,3,4,5,6,7,8,9,10,12,14,16,17,18,19,20,21,22,24,26,28,30,32,34,36,38,40,45,50,60,70,89]
#uneven_levels = [1,2,3,4,5,6,7,8,9,10,12,14,16,17,18,19,20,21,22,24,26,28,30,32,34,36,38,40,45,50,60,70,89]
uneven_levels1 = [1,2,3,4,5,6,7,8,9,10,12,14,16,17,18,19,20,21,22,24,26,28,30,32,34,36,38,40,45,50,60,70,80,90,100,110,120]


uneven_levels2 = [-89,-60,-50,-40,-30,-28,-26,-24,-22,-20,-18,-16,-14,-12,-10,-8,-6,-4,-2,0,2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,40,50,60,89]

size = '0-700 sqr.degs.'

easy_plot_unevencolorbars(lon=xx,lat=yy, data=events_0_700,
                  uneven_levels=uneven_levels1,name=f' {variablename} {size}',units=f'{cbarunits}',i=0,j=0)

easy_plot_unevencolorbars(lon=xx,lat=yy, data=hot_events_0_700,
                  uneven_levels=uneven_levels1,name=f' Hot Season {variablename} {size}',units=f'{cbarunits}',i=0,j=1)

easy_plot_unevencolorbars(lon=xx,lat=yy, data=cold_events_0_700,
                  uneven_levels=uneven_levels1,name= f' Cold Season {variablename} {size}',units=f'{cbarunits}',i=0,j=2)

easy_plot_unevencolorbars2(lon=xx,lat=yy, data=hot_events_0_700-cold_events_0_700,
                  uneven_levels=uneven_levels2,name=f'Difference {size}',units=f'{cbarunits}',i=0,j=3)

size = '0-25 sqr.degs.'

easy_plot_unevencolorbars(lon=xx,lat=yy, data=events_0_25,
                  uneven_levels=uneven_levels1,name=f' {variablename} {size}',units=f'{cbarunits}',i=1,j=0)

easy_plot_unevencolorbars(lon=xx,lat=yy, data=hot_events_0_25,
                  uneven_levels=uneven_levels1,name=f' Hot Season {variablename} {size}',units=f'{cbarunits}',i=1,j=1)

easy_plot_unevencolorbars(lon=xx,lat=yy, data=cold_events_0_25,
                  uneven_levels=uneven_levels1,name= f' Cold Season {variablename} {size}',units=f'{cbarunits}',i=1,j=2)

easy_plot_unevencolorbars2(lon=xx,lat=yy, data=hot_events_0_25-cold_events_0_25,
                  uneven_levels=uneven_levels2,name=f'Difference {size}',units=f'{cbarunits}',i=1,j=3)

size = '25-700 sqr.degs.'

easy_plot_unevencolorbars(lon=xx,lat=yy, data=events_25_700,
                  uneven_levels=uneven_levels1,name=f' {variablename} {size}',units=f'{cbarunits}',i=2,j=0)

easy_plot_unevencolorbars(lon=xx,lat=yy, data=hot_events_25_700,
                  uneven_levels=uneven_levels1,name=f' Hot Season {variablename} {size}',units=f'{cbarunits}',i=2,j=1)

easy_plot_unevencolorbars(lon=xx,lat=yy, data=cold_events_25_700,
                  uneven_levels=uneven_levels1,name= f' Cold Season {variablename} {size}',units=f'{cbarunits}',i=2,j=2)

easy_plot_unevencolorbars2(lon=xx,lat=yy, data=hot_events_25_700-cold_events_25_700,
                  uneven_levels=uneven_levels2,name=f'Difference {size}',units=f'{cbarunits}',i=2,j=3)

plt.tight_layout()
plt.savefig(f'/home/shilpa/glory_mat_analysis/compare_no_of_events_hot_cold_seasons2.png')




# plot
nrows=3
ncols=4

cbarunits = '°C'
variablename = 'Mean Mean Intensity'

fig, axs = plt.subplots(nrows=nrows,ncols=ncols,subplot_kw={'projection': proj_crs},figsize=(22,15))

uneven_levels1 = [0.5,0.6,0.7,0.8,0.9,1,1.1,1.2,1.3,1.4,1.5,1.6,1.8,2,2.2,2.4,2.6,2.8,3,3.5]

uneven_levels2 = [-3.5,-2,-1.25,-1,-0.8,-0.6,-0.4,-0.2,0,0.2,0.4,0.6,0.8,1,1.25,2,3.5]


size = '0-700 sqr.degs.'

easy_plot_unevencolorbars(lon=xx,lat=yy, data=mean_meanintensity_0_700,
                  uneven_levels=uneven_levels1,name=f' {variablename} {size}',units=f'{cbarunits}',i=0,j=0)

easy_plot_unevencolorbars(lon=xx,lat=yy, data=mean_hot_meanintensity_0_700,
                  uneven_levels=uneven_levels1,name=f' Hot Season {variablename} {size}',units=f'{cbarunits}',i=0,j=1)

easy_plot_unevencolorbars(lon=xx,lat=yy, data=mean_cold_meanintensity_0_700,
                  uneven_levels=uneven_levels1,name= f' Cold Season {variablename} {size}',units=f'{cbarunits}',i=0,j=2)

easy_plot_unevencolorbars2(lon=xx,lat=yy, data=mean_hot_meanintensity_0_700-mean_cold_meanintensity_0_700,
                  uneven_levels=uneven_levels2,name=f'Difference {size}',units=f'{cbarunits}',i=0,j=3)

size = '0-25 sqr.degs.'

easy_plot_unevencolorbars(lon=xx,lat=yy, data=mean_meanintensity_0_25,
                  uneven_levels=uneven_levels1,name=f' {variablename} {size}',units=f'{cbarunits}',i=1,j=0)

easy_plot_unevencolorbars(lon=xx,lat=yy, data=mean_hot_meanintensity_0_25,
                  uneven_levels=uneven_levels1,name=f' Hot Season {variablename} {size}',units=f'{cbarunits}',i=1,j=1)

easy_plot_unevencolorbars(lon=xx,lat=yy, data=mean_cold_meanintensity_0_25,
                  uneven_levels=uneven_levels1,name= f' Cold Season {variablename} {size}',units=f'{cbarunits}',i=1,j=2)

easy_plot_unevencolorbars2(lon=xx,lat=yy, data=mean_hot_meanintensity_0_25-mean_cold_meanintensity_0_25,
                  uneven_levels=uneven_levels2,name=f'Difference {size}',units=f'{cbarunits}',i=1,j=3)

size = '25-700 sqr.degs.'

easy_plot_unevencolorbars(lon=xx,lat=yy, data=mean_meanintensity_25_700,
                  uneven_levels=uneven_levels1,name=f' {variablename} {size}',units=f'{cbarunits}',i=2,j=0)

easy_plot_unevencolorbars(lon=xx,lat=yy, data=mean_hot_meanintensity_25_700,
                  uneven_levels=uneven_levels1,name=f' Hot Season {variablename} {size}',units=f'{cbarunits}',i=2,j=1)

easy_plot_unevencolorbars(lon=xx,lat=yy, data=mean_cold_meanintensity_25_700,
                  uneven_levels=uneven_levels1,name= f' Cold Season {variablename} {size}',units=f'{cbarunits}',i=2,j=2)

easy_plot_unevencolorbars2(lon=xx,lat=yy, data=mean_hot_meanintensity_25_700-mean_cold_meanintensity_25_700,
                  uneven_levels=uneven_levels2,name=f'Difference {size}',units=f'{cbarunits}',i=2,j=3)

plt.tight_layout()
plt.savefig(f'/home/shilpa/glory_mat_analysis/compare_mean_meanintensity_hot_cold_seasons2.png')


# plot
nrows=3
ncols=4

cbarunits = '°C'
variablename = 'Median Mean Intensity'

fig, axs = plt.subplots(nrows=nrows,ncols=ncols,subplot_kw={'projection': proj_crs},figsize=(22,15))

uneven_levels1 = [0.5,0.6,0.7,0.8,0.9,1,1.1,1.2,1.3,1.4,1.5,1.6,1.8,2,2.2,2.4,2.6,2.8,3,3.5]

uneven_levels2 = [-3.5,-2,-1.25,-1,-0.8,-0.6,-0.4,-0.2,0,0.2,0.4,0.6,0.8,1,1.25,2,3.5]

size = '0-700 sqr.degs.'

easy_plot_unevencolorbars(lon=xx,lat=yy, data=median_meanintensity_0_700,
                  uneven_levels=uneven_levels1,name=f' {variablename} {size}',units=f'{cbarunits}',i=0,j=0)

easy_plot_unevencolorbars(lon=xx,lat=yy, data=median_hot_meanintensity_0_700,
                  uneven_levels=uneven_levels1,name=f' Hot Season {variablename} {size}',units=f'{cbarunits}',i=0,j=1)

easy_plot_unevencolorbars(lon=xx,lat=yy, data=median_cold_meanintensity_0_700,
                  uneven_levels=uneven_levels1,name= f' Cold Season {variablename} {size}',units=f'{cbarunits}',i=0,j=2)

easy_plot_unevencolorbars2(lon=xx,lat=yy, data=median_hot_meanintensity_0_700-median_cold_meanintensity_0_700,
                  uneven_levels=uneven_levels2,name=f'Difference {size}',units=f'{cbarunits}',i=0,j=3)

size = '0-25 sqr.degs.'

easy_plot_unevencolorbars(lon=xx,lat=yy, data=median_meanintensity_0_25,
                  uneven_levels=uneven_levels1,name=f' {variablename} {size}',units=f'{cbarunits}',i=1,j=0)

easy_plot_unevencolorbars(lon=xx,lat=yy, data=median_hot_meanintensity_0_25,
                  uneven_levels=uneven_levels1,name=f' Hot Season {variablename} {size}',units=f'{cbarunits}',i=1,j=1)

easy_plot_unevencolorbars(lon=xx,lat=yy, data=median_cold_meanintensity_0_25,
                  uneven_levels=uneven_levels1,name= f' Cold Season {variablename} {size}',units=f'{cbarunits}',i=1,j=2)

easy_plot_unevencolorbars2(lon=xx,lat=yy, data=median_hot_meanintensity_0_25-median_cold_meanintensity_0_25,
                  uneven_levels=uneven_levels2,name=f'Difference {size}',units=f'{cbarunits}',i=1,j=3)

size = '25-700 sqr.degs.'
easy_plot_unevencolorbars(lon=xx,lat=yy, data=median_meanintensity_25_700,
                  uneven_levels=uneven_levels1,name=f' {variablename} {size}',units=f'{cbarunits}',i=2,j=0)

easy_plot_unevencolorbars(lon=xx,lat=yy, data=median_hot_meanintensity_25_700,
                  uneven_levels=uneven_levels1,name=f' Hot Season {variablename} {size}',units=f'{cbarunits}',i=2,j=1)

easy_plot_unevencolorbars(lon=xx,lat=yy, data=median_cold_meanintensity_25_700,
                  uneven_levels=uneven_levels1,name= f' Cold Season {variablename} {size}',units=f'{cbarunits}',i=2,j=2)

easy_plot_unevencolorbars2(lon=xx,lat=yy, data=median_hot_meanintensity_25_700-median_cold_meanintensity_25_700,
                  uneven_levels=uneven_levels2,name=f'Difference {size}',units=f'{cbarunits}',i=2,j=3)

plt.tight_layout()
plt.savefig(f'/home/shilpa/glory_mat_analysis/compare_median_meanintensity_hot_cold_seasons2.png')



# plot
nrows=3
ncols=4

cbarunits = '°C * Days'
variablename = 'Mean Cumulative Intensity'

fig, axs = plt.subplots(nrows=nrows,ncols=ncols,subplot_kw={'projection': proj_crs},figsize=(22,15))

uneven_levels1 = [5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,40,45,50,55,60,65,70,80,90,100,150,200,250,300,350]

uneven_levels2 = [-350,-100,-50,-10,-5,0,5,10,50,100,350]

size = '0-700 sqr.degs.'

easy_plot_unevencolorbars(lon=xx,lat=yy, data=mean_cumintensity_0_700,
                  uneven_levels=uneven_levels1,name=f' {variablename} {size}',units=f'{cbarunits}',i=0,j=0)

easy_plot_unevencolorbars(lon=xx,lat=yy, data=mean_hot_cumintensity_0_700,
                  uneven_levels=uneven_levels1,name=f' Hot Season {variablename} {size}',units=f'{cbarunits}',i=0,j=1)

easy_plot_unevencolorbars(lon=xx,lat=yy, data=mean_cold_cumintensity_0_700,
                  uneven_levels=uneven_levels1,name= f' Cold Season {variablename} {size}',units=f'{cbarunits}',i=0,j=2)

easy_plot_unevencolorbars2(lon=xx,lat=yy, data=mean_hot_cumintensity_0_700-mean_cold_cumintensity_0_700,
                  uneven_levels=uneven_levels2,name=f'Difference {size}',units=f'{cbarunits}',i=0,j=3)

size = '0-25 sqr.degs.'

easy_plot_unevencolorbars(lon=xx,lat=yy, data=mean_cumintensity_0_25,
                  uneven_levels=uneven_levels1,name=f' {variablename} {size}',units=f'{cbarunits}',i=1,j=0)

easy_plot_unevencolorbars(lon=xx,lat=yy, data=mean_hot_cumintensity_0_25,
                  uneven_levels=uneven_levels1,name=f' Hot Season {variablename} {size}',units=f'{cbarunits}',i=1,j=1)

easy_plot_unevencolorbars(lon=xx,lat=yy, data=mean_cold_cumintensity_0_25,
                  uneven_levels=uneven_levels1,name= f' Cold Season {variablename} {size}',units=f'{cbarunits}',i=1,j=2)

easy_plot_unevencolorbars2(lon=xx,lat=yy, data=mean_hot_cumintensity_0_25-mean_cold_cumintensity_0_25,
                  uneven_levels=uneven_levels2,name=f'Difference {size}',units=f'{cbarunits}',i=1,j=3)

size = '25-700 sqr.degs.'

easy_plot_unevencolorbars(lon=xx,lat=yy, data=mean_cumintensity_25_700,
                  uneven_levels=uneven_levels1,name=f' {variablename} {size}',units=f'{cbarunits}',i=2,j=0)

easy_plot_unevencolorbars(lon=xx,lat=yy, data=mean_hot_cumintensity_25_700,
                  uneven_levels=uneven_levels1,name=f' Hot Season {variablename} {size}',units=f'{cbarunits}',i=2,j=1)

easy_plot_unevencolorbars(lon=xx,lat=yy, data=mean_cold_cumintensity_25_700,
                  uneven_levels=uneven_levels1,name= f' Cold Season {variablename} {size}',units=f'{cbarunits}',i=2,j=2)

easy_plot_unevencolorbars2(lon=xx,lat=yy, data=mean_hot_cumintensity_25_700-mean_cold_cumintensity_25_700,
                  uneven_levels=uneven_levels2,name=f'Difference {size}',units=f'{cbarunits}',i=2,j=3)

plt.tight_layout()
plt.savefig(f'/home/shilpa/glory_mat_analysis/compare_mean_cumintensity_hot_cold_seasons2.png')


# plot
nrows=3
ncols=4

cbarunits = '°C * Days'
variablename = 'Median Cumulative Intensity'

fig, axs = plt.subplots(nrows=nrows,ncols=ncols,subplot_kw={'projection': proj_crs},figsize=(22,15))

uneven_levels1 = [5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,40,45,50,55,60,65,70,80,90,100,150,200,250,300,350]

uneven_levels2 = [-350,-100,-50,-10,-5,0,5,10,50,100,350]

size = '0-700 sqr.degs.'

easy_plot_unevencolorbars(lon=xx,lat=yy, data=median_cumintensity_0_700,
                  uneven_levels=uneven_levels1,name=f' {variablename} {size}',units=f'{cbarunits}',i=0,j=0)

easy_plot_unevencolorbars(lon=xx,lat=yy, data=median_hot_cumintensity_0_700,
                  uneven_levels=uneven_levels1,name=f' Hot Season {variablename} {size}',units=f'{cbarunits}',i=0,j=1)

easy_plot_unevencolorbars(lon=xx,lat=yy, data=median_cold_cumintensity_0_700,
                  uneven_levels=uneven_levels1,name= f' Cold Season {variablename} {size}',units=f'{cbarunits}',i=0,j=2)

easy_plot_unevencolorbars2(lon=xx,lat=yy, data=median_hot_cumintensity_0_700-median_cold_cumintensity_0_700,
                  uneven_levels=uneven_levels2,name=f'Difference {size}',units=f'{cbarunits}',i=0,j=3)

size = '0-25 sqr.degs.'

easy_plot_unevencolorbars(lon=xx,lat=yy, data=median_cumintensity_0_25,
                  uneven_levels=uneven_levels1,name=f' {variablename} {size}',units=f'{cbarunits}',i=1,j=0)

easy_plot_unevencolorbars(lon=xx,lat=yy, data=median_hot_cumintensity_0_25,
                  uneven_levels=uneven_levels1,name=f' Hot Season {variablename} {size}',units=f'{cbarunits}',i=1,j=1)

easy_plot_unevencolorbars(lon=xx,lat=yy, data=median_cold_cumintensity_0_25,
                  uneven_levels=uneven_levels1,name= f' Cold Season {variablename} {size}',units=f'{cbarunits}',i=1,j=2)

easy_plot_unevencolorbars2(lon=xx,lat=yy, data=median_hot_cumintensity_0_25-median_cold_cumintensity_0_25,
                  uneven_levels=uneven_levels2,name=f'Difference {size}',units=f'{cbarunits}',i=1,j=3)

size = '25-700 sqr.degs.'

easy_plot_unevencolorbars(lon=xx,lat=yy, data=median_cumintensity_25_700,
                  uneven_levels=uneven_levels1,name=f' {variablename} {size}',units=f'{cbarunits}',i=2,j=0)

easy_plot_unevencolorbars(lon=xx,lat=yy, data=median_hot_cumintensity_25_700,
                  uneven_levels=uneven_levels1,name=f' Hot Season {variablename} {size}',units=f'{cbarunits}',i=2,j=1)

easy_plot_unevencolorbars(lon=xx,lat=yy, data=median_cold_cumintensity_25_700,
                  uneven_levels=uneven_levels1,name= f' Cold Season {variablename} {size}',units=f'{cbarunits}',i=2,j=2)

easy_plot_unevencolorbars2(lon=xx,lat=yy, data=median_hot_cumintensity_25_700-median_cold_cumintensity_25_700,
                  uneven_levels=uneven_levels2,name=f'Difference {size}',units=f'{cbarunits}',i=2,j=3)

plt.tight_layout()
plt.savefig(f'/home/shilpa/glory_mat_analysis/compare_median_cumintensity_hot_cold_seasons2.png')




# plot
nrows=3
ncols=4

cbarunits = '°C/Day'
variablename = 'Mean Rate onset'

fig, axs = plt.subplots(nrows=nrows,ncols=ncols,subplot_kw={'projection': proj_crs},figsize=(22,15))


uneven_levels1 = [0,0.05,0.06,0.07,0.08,0.09,0.1,0.12,0.14,0.16,0.18,0.2,0.22,0.24,0.26,0.28,0.3,0.35,0.45, 0.7]


uneven_levels2 = [-0.7,-0.5,-0.1,-0.05,0,0.05,0.1,0.5,0.7]

size = '0-700 sqr.degs.'

easy_plot_unevencolorbars(lon=xx,lat=yy, data=mean_rateonset_0_700,
                  uneven_levels=uneven_levels1,name=f' {variablename} {size}',units=f'{cbarunits}',i=0,j=0)

easy_plot_unevencolorbars(lon=xx,lat=yy, data=mean_hot_rateonset_0_700,
                  uneven_levels=uneven_levels1,name=f' Hot Season {variablename} {size}',units=f'{cbarunits}',i=0,j=1)

easy_plot_unevencolorbars(lon=xx,lat=yy, data=mean_cold_rateonset_0_700,
                  uneven_levels=uneven_levels1,name= f' Cold Season {variablename} {size}',units=f'{cbarunits}',i=0,j=2)

easy_plot_unevencolorbars2(lon=xx,lat=yy, data=mean_hot_rateonset_0_700-mean_cold_rateonset_0_700,
                  uneven_levels=uneven_levels2,name=f'Difference {size}',units=f'{cbarunits}',i=0,j=3)

size = '0-25 sqr.degs.'

easy_plot_unevencolorbars(lon=xx,lat=yy, data=mean_rateonset_0_25,
                  uneven_levels=uneven_levels1,name=f' {variablename} {size}',units=f'{cbarunits}',i=1,j=0)

easy_plot_unevencolorbars(lon=xx,lat=yy, data=mean_hot_rateonset_0_25,
                  uneven_levels=uneven_levels1,name=f' Hot Season {variablename} {size}',units=f'{cbarunits}',i=1,j=1)

easy_plot_unevencolorbars(lon=xx,lat=yy, data=mean_cold_rateonset_0_25,
                  uneven_levels=uneven_levels1,name= f' Cold Season {variablename} {size}',units=f'{cbarunits}',i=1,j=2)

easy_plot_unevencolorbars2(lon=xx,lat=yy, data=mean_hot_rateonset_0_25-mean_cold_rateonset_0_25,
                  uneven_levels=uneven_levels2,name=f'Difference {size}',units=f'{cbarunits}',i=1,j=3)

size = '25-700 sqr.degs.'

easy_plot_unevencolorbars(lon=xx,lat=yy, data=mean_rateonset_25_700,
                  uneven_levels=uneven_levels1,name=f' {variablename} {size}',units=f'{cbarunits}',i=2,j=0)

easy_plot_unevencolorbars(lon=xx,lat=yy, data=mean_hot_rateonset_25_700,
                  uneven_levels=uneven_levels1,name=f' Hot Season {variablename} {size}',units=f'{cbarunits}',i=2,j=1)

easy_plot_unevencolorbars(lon=xx,lat=yy, data=mean_cold_rateonset_25_700,
                  uneven_levels=uneven_levels1,name= f' Cold Season {variablename} {size}',units=f'{cbarunits}',i=2,j=2)

easy_plot_unevencolorbars2(lon=xx,lat=yy, data=mean_hot_rateonset_25_700-mean_cold_rateonset_25_700,
                  uneven_levels=uneven_levels2,name=f'Difference {size}',units=f'{cbarunits}',i=2,j=3)

plt.tight_layout()
plt.savefig(f'/home/shilpa/glory_mat_analysis/compare_mean_rateonset_hot_cold_seasons2.png')


# plot
nrows=3
ncols=4

cbarunits = '°C/Day'
variablename = 'Median Rate onset'

fig, axs = plt.subplots(nrows=nrows,ncols=ncols,subplot_kw={'projection': proj_crs},figsize=(22,15))


uneven_levels1 = [0,0.05,0.06,0.07,0.08,0.09,0.1,0.12,0.14,0.16,0.18,0.2,0.22,0.24,0.26,0.28,0.3,0.35,0.45, 0.7]


uneven_levels2 = [-0.7,-0.5,-0.1,-0.05,0,0.05,0.1,0.5,0.7]

size = '0-700 sqr.degs.'

easy_plot_unevencolorbars(lon=xx,lat=yy, data=median_rateonset_0_700,
                  uneven_levels=uneven_levels1,name=f' {variablename} {size}',units=f'{cbarunits}',i=0,j=0)

easy_plot_unevencolorbars(lon=xx,lat=yy, data=median_hot_rateonset_0_700,
                  uneven_levels=uneven_levels1,name=f' Hot Season {variablename} {size}',units=f'{cbarunits}',i=0,j=1)

easy_plot_unevencolorbars(lon=xx,lat=yy, data=median_cold_rateonset_0_700,
                  uneven_levels=uneven_levels1,name= f' Cold Season {variablename} {size}',units=f'{cbarunits}',i=0,j=2)

easy_plot_unevencolorbars2(lon=xx,lat=yy, data=median_hot_rateonset_0_700-median_cold_rateonset_0_700,
                  uneven_levels=uneven_levels2,name=f'Difference {size}',units=f'{cbarunits}',i=0,j=3)

size = '0-25 sqr.degs.'

easy_plot_unevencolorbars(lon=xx,lat=yy, data=median_rateonset_0_25,
                  uneven_levels=uneven_levels1,name=f' {variablename} {size}',units=f'{cbarunits}',i=1,j=0)

easy_plot_unevencolorbars(lon=xx,lat=yy, data=median_hot_rateonset_0_25,
                  uneven_levels=uneven_levels1,name=f' Hot Season {variablename} {size}',units=f'{cbarunits}',i=1,j=1)

easy_plot_unevencolorbars(lon=xx,lat=yy, data=median_cold_rateonset_0_25,
                  uneven_levels=uneven_levels1,name= f' Cold Season {variablename} {size}',units=f'{cbarunits}',i=1,j=2)

easy_plot_unevencolorbars2(lon=xx,lat=yy, data=median_hot_rateonset_0_25-median_cold_rateonset_0_25,
                  uneven_levels=uneven_levels2,name=f'Difference {size}',units=f'{cbarunits}',i=1,j=3)

size = '25-700 sqr.degs.'

easy_plot_unevencolorbars(lon=xx,lat=yy, data=median_rateonset_25_700,
                  uneven_levels=uneven_levels1,name=f' {variablename} {size}',units=f'{cbarunits}',i=2,j=0)

easy_plot_unevencolorbars(lon=xx,lat=yy, data=median_hot_rateonset_25_700,
                  uneven_levels=uneven_levels1,name=f' Hot Season {variablename} {size}',units=f'{cbarunits}',i=2,j=1)

easy_plot_unevencolorbars(lon=xx,lat=yy, data=median_cold_rateonset_25_700,
                  uneven_levels=uneven_levels1,name= f' Cold Season {variablename} {size}',units=f'{cbarunits}',i=2,j=2)

easy_plot_unevencolorbars2(lon=xx,lat=yy, data=median_hot_rateonset_25_700-median_cold_rateonset_25_700,
                  uneven_levels=uneven_levels2,name=f'Difference {size}',units=f'{cbarunits}',i=2,j=3)

plt.tight_layout()
plt.savefig(f'/home/shilpa/glory_mat_analysis/compare_median_rateonset_hot_cold_seasons2.png')



# plot
nrows=3
ncols=4

cbarunits = '°C/Day'
variablename = 'Mean Rate Decline'

fig, axs = plt.subplots(nrows=nrows,ncols=ncols,subplot_kw={'projection': proj_crs},figsize=(22,15))

#uneven_levels1 = [0,0.05,0.06,0.07,0.08,0.09,0.1,0.11,0.12,0.13,0.14,0.15,0.16,0.17,0.18,0.19,0.2,0.21,0.22,0.23,0.24, 0.3]

#uneven_levels2 = [-0.3,-0.2,-0.15,-0.1,-0.01,0,0.01,0.1,0.15,0.2,0.3]

uneven_levels1 = [0,0.05,0.06,0.07,0.08,0.09,0.1,0.12,0.14,0.16,0.18,0.2,0.22,0.24,0.26,0.28,0.3,0.35,0.45, 0.7]


uneven_levels2 = [-0.7,-0.5,-0.1,-0.05,0,0.05,0.1,0.5,0.7]


size = '0-700 sqr.degs.'

easy_plot_unevencolorbars(lon=xx,lat=yy, data=mean_ratedecline_0_700,
                  uneven_levels=uneven_levels1,name=f' {variablename} {size}',units=f'{cbarunits}',i=0,j=0)

easy_plot_unevencolorbars(lon=xx,lat=yy, data=mean_hot_ratedecline_0_700,
                  uneven_levels=uneven_levels1,name=f' Hot Season {variablename} {size}',units=f'{cbarunits}',i=0,j=1)

easy_plot_unevencolorbars(lon=xx,lat=yy, data=mean_cold_ratedecline_0_700,
                  uneven_levels=uneven_levels1,name= f' Cold Season {variablename} {size}',units=f'{cbarunits}',i=0,j=2)

easy_plot_unevencolorbars2(lon=xx,lat=yy, data=mean_hot_ratedecline_0_700-mean_cold_ratedecline_0_700,
                  uneven_levels=uneven_levels2,name=f'Difference {size}',units=f'{cbarunits}',i=0,j=3)

size = '0-25 sqr.degs.'

easy_plot_unevencolorbars(lon=xx,lat=yy, data=mean_ratedecline_0_25,
                  uneven_levels=uneven_levels1,name=f' {variablename} {size}',units=f'{cbarunits}',i=1,j=0)

easy_plot_unevencolorbars(lon=xx,lat=yy, data=mean_hot_ratedecline_0_25,
                  uneven_levels=uneven_levels1,name=f' Hot Season {variablename} {size}',units=f'{cbarunits}',i=1,j=1)

easy_plot_unevencolorbars(lon=xx,lat=yy, data=mean_cold_ratedecline_0_25,
                  uneven_levels=uneven_levels1,name= f' Cold Season {variablename} {size}',units=f'{cbarunits}',i=1,j=2)

easy_plot_unevencolorbars2(lon=xx,lat=yy, data=mean_hot_ratedecline_0_25-mean_cold_ratedecline_0_25,
                  uneven_levels=uneven_levels2,name=f'Difference {size}',units=f'{cbarunits}',i=1,j=3)

size = '25-700 sqr.degs.'

easy_plot_unevencolorbars(lon=xx,lat=yy, data=mean_ratedecline_25_700,
                  uneven_levels=uneven_levels1,name=f' {variablename} {size}',units=f'{cbarunits}',i=2,j=0)

easy_plot_unevencolorbars(lon=xx,lat=yy, data=mean_hot_ratedecline_25_700,
                  uneven_levels=uneven_levels1,name=f' Hot Season {variablename} {size}',units=f'{cbarunits}',i=2,j=1)

easy_plot_unevencolorbars(lon=xx,lat=yy, data=mean_cold_ratedecline_25_700,
                  uneven_levels=uneven_levels1,name= f' Cold Season {variablename} {size}',units=f'{cbarunits}',i=2,j=2)

easy_plot_unevencolorbars2(lon=xx,lat=yy, data=mean_hot_ratedecline_25_700-mean_cold_ratedecline_25_700,
                  uneven_levels=uneven_levels2,name=f'Difference {size}',units=f'{cbarunits}',i=2,j=3)

plt.tight_layout()
plt.savefig(f'/home/shilpa/glory_mat_analysis/compare_mean_ratedecline_hot_cold_seasons2.png')



# plot
nrows=3
ncols=4

cbarunits = '°C/Day'
variablename = 'Median Rate Decline'

fig, axs = plt.subplots(nrows=nrows,ncols=ncols,subplot_kw={'projection': proj_crs},figsize=(22,15))

#uneven_levels1 = [0,0.05,0.06,0.07,0.08,0.09,0.1,0.11,0.12,0.13,0.14,0.15,0.16,0.17,0.18,0.19,0.2,0.21,0.22,0.23,0.24, 0.3]

#uneven_levels2 = [-0.3,-0.2,-0.15,-0.1,-0.01,0,0.01,0.1,0.15,0.2,0.3]
uneven_levels1 = [0,0.05,0.06,0.07,0.08,0.09,0.1,0.12,0.14,0.16,0.18,0.2,0.22,0.24,0.26,0.28,0.3,0.35,0.45, 0.7]


uneven_levels2 = [-0.7,-0.5,-0.1,-0.05,0,0.05,0.1,0.5,0.7]


size = '0-700 sqr.degs.'

easy_plot_unevencolorbars(lon=xx,lat=yy, data=median_ratedecline_0_700,
                  uneven_levels=uneven_levels1,name=f' {variablename} {size}',units=f'{cbarunits}',i=0,j=0)

easy_plot_unevencolorbars(lon=xx,lat=yy, data=median_hot_ratedecline_0_700,
                  uneven_levels=uneven_levels1,name=f' Hot Season {variablename} {size}',units=f'{cbarunits}',i=0,j=1)

easy_plot_unevencolorbars(lon=xx,lat=yy, data=median_cold_ratedecline_0_700,
                  uneven_levels=uneven_levels1,name= f' Cold Season {variablename} {size}',units=f'{cbarunits}',i=0,j=2)

easy_plot_unevencolorbars2(lon=xx,lat=yy, data=median_hot_ratedecline_0_700-median_cold_ratedecline_0_700,
                  uneven_levels=uneven_levels2,name=f'Difference {size}',units=f'{cbarunits}',i=0,j=3)

size = '0-25 sqr.degs.'

easy_plot_unevencolorbars(lon=xx,lat=yy, data=median_ratedecline_0_25,
                  uneven_levels=uneven_levels1,name=f' {variablename} {size}',units=f'{cbarunits}',i=1,j=0)

easy_plot_unevencolorbars(lon=xx,lat=yy, data=median_hot_ratedecline_0_25,
                  uneven_levels=uneven_levels1,name=f' Hot Season {variablename} {size}',units=f'{cbarunits}',i=1,j=1)

easy_plot_unevencolorbars(lon=xx,lat=yy, data=median_cold_ratedecline_0_25,
                  uneven_levels=uneven_levels1,name= f' Cold Season {variablename} {size}',units=f'{cbarunits}',i=1,j=2)

easy_plot_unevencolorbars2(lon=xx,lat=yy, data=median_hot_ratedecline_0_25-median_cold_ratedecline_0_25,
                  uneven_levels=uneven_levels2,name=f'Difference {size}',units=f'{cbarunits}',i=1,j=3)

size = '25-700 sqr.degs.'

easy_plot_unevencolorbars(lon=xx,lat=yy, data=median_ratedecline_25_700,
                  uneven_levels=uneven_levels1,name=f' {variablename} {size}',units=f'{cbarunits}',i=2,j=0)

easy_plot_unevencolorbars(lon=xx,lat=yy, data=median_hot_ratedecline_25_700,
                  uneven_levels=uneven_levels1,name=f' Hot Season {variablename} {size}',units=f'{cbarunits}',i=2,j=1)

easy_plot_unevencolorbars(lon=xx,lat=yy, data=median_cold_ratedecline_25_700,
                  uneven_levels=uneven_levels1,name= f' Cold Season {variablename} {size}',units=f'{cbarunits}',i=2,j=2)

easy_plot_unevencolorbars2(lon=xx,lat=yy, data=median_hot_ratedecline_25_700-median_cold_ratedecline_25_700,
                  uneven_levels=uneven_levels2,name=f'Difference {size}',units=f'{cbarunits}',i=2,j=3)

plt.tight_layout()
plt.savefig(f'/home/shilpa/glory_mat_analysis/compare_median_ratedecline_hot_cold_seasons2.png')



# plot
nrows=3
ncols=4

cbarunits = 'Days'
variablename = 'Mean Reaction Window'

fig, axs = plt.subplots(nrows=nrows,ncols=ncols,subplot_kw={'projection': proj_crs},figsize=(22,15))

uneven_levels1 = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,40,45,50,55,60,65,70,80,90,100,150,200,250,300,350,600]

uneven_levels2 = [-600,-350,-100,-50,-10,-5,0,5,10,50,100,350,600]

size = '0-700 sqr.degs.'

easy_plot_unevencolorbars(lon=xx,lat=yy, data=mean_reaction_window_0_700,
                  uneven_levels=uneven_levels1,name=f' {variablename} {size}',units=f'{cbarunits}',i=0,j=0)

easy_plot_unevencolorbars(lon=xx,lat=yy, data=mean_hot_reaction_window_0_700,
                  uneven_levels=uneven_levels1,name=f' Hot Season {variablename} {size}',units=f'{cbarunits}',i=0,j=1)

easy_plot_unevencolorbars(lon=xx,lat=yy, data=mean_cold_reaction_window_0_700,
                  uneven_levels=uneven_levels1,name= f' Cold Season {variablename} {size}',units=f'{cbarunits}',i=0,j=2)

easy_plot_unevencolorbars2(lon=xx,lat=yy, data=mean_hot_reaction_window_0_700-mean_cold_reaction_window_0_700,
                  uneven_levels=uneven_levels2,name=f'Difference {size}',units=f'{cbarunits}',i=0,j=3)

size = '0-25 sqr.degs.'

easy_plot_unevencolorbars(lon=xx,lat=yy, data=mean_reaction_window_0_25,
                  uneven_levels=uneven_levels1,name=f' {variablename} {size}',units=f'{cbarunits}',i=1,j=0)

easy_plot_unevencolorbars(lon=xx,lat=yy, data=mean_hot_reaction_window_0_25,
                  uneven_levels=uneven_levels1,name=f' Hot Season {variablename} {size}',units=f'{cbarunits}',i=1,j=1)

easy_plot_unevencolorbars(lon=xx,lat=yy, data=mean_cold_reaction_window_0_25,
                  uneven_levels=uneven_levels1,name= f' Cold Season {variablename} {size}',units=f'{cbarunits}',i=1,j=2)

easy_plot_unevencolorbars2(lon=xx,lat=yy, data=mean_hot_reaction_window_0_25-mean_cold_reaction_window_0_25,
                  uneven_levels=uneven_levels2,name=f'Difference {size}',units=f'{cbarunits}',i=1,j=3)

size = '25-700 sqr.degs.'

easy_plot_unevencolorbars(lon=xx,lat=yy, data=mean_reaction_window_25_700,
                  uneven_levels=uneven_levels1,name=f' {variablename} {size}',units=f'{cbarunits}',i=2,j=0)

easy_plot_unevencolorbars(lon=xx,lat=yy, data=mean_hot_reaction_window_25_700,
                  uneven_levels=uneven_levels1,name=f' Hot Season {variablename} {size}',units=f'{cbarunits}',i=2,j=1)

easy_plot_unevencolorbars(lon=xx,lat=yy, data=mean_cold_reaction_window_25_700,
                  uneven_levels=uneven_levels1,name= f' Cold Season {variablename} {size}',units=f'{cbarunits}',i=2,j=2)

easy_plot_unevencolorbars2(lon=xx,lat=yy, data=mean_hot_reaction_window_25_700-mean_cold_reaction_window_25_700,
                  uneven_levels=uneven_levels2,name=f'Difference {size}',units=f'{cbarunits}',i=2,j=3)

plt.tight_layout()
plt.savefig(f'/home/shilpa/glory_mat_analysis/compare_mean_reactionwindow_hot_cold_seasons2.png')



# plot
nrows=3
ncols=4

cbarunits = 'Days'
variablename = 'Median Reaction Window'

fig, axs = plt.subplots(nrows=nrows,ncols=ncols,subplot_kw={'projection': proj_crs},figsize=(22,15))

uneven_levels1 = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,40,45,50,55,60,65,70,80,90,100,150,200,250,300,350,600]

uneven_levels2 = [-600,-350,-100,-50,-10,-5,0,5,10,50,100,350,600]

size = '0-700 sqr.degs.'

easy_plot_unevencolorbars(lon=xx,lat=yy, data=median_reaction_window_0_700,
                  uneven_levels=uneven_levels1,name=f' {variablename} {size}',units=f'{cbarunits}',i=0,j=0)

easy_plot_unevencolorbars(lon=xx,lat=yy, data=median_hot_reaction_window_0_700,
                  uneven_levels=uneven_levels1,name=f' Hot Season {variablename} {size}',units=f'{cbarunits}',i=0,j=1)

easy_plot_unevencolorbars(lon=xx,lat=yy, data=median_cold_reaction_window_0_700,
                  uneven_levels=uneven_levels1,name= f' Cold Season {variablename} {size}',units=f'{cbarunits}',i=0,j=2)

easy_plot_unevencolorbars2(lon=xx,lat=yy, data=median_hot_reaction_window_0_700-median_cold_reaction_window_0_700,
                  uneven_levels=uneven_levels2,name=f'Difference {size}',units=f'{cbarunits}',i=0,j=3)

size = '0-25 sqr.degs.'

easy_plot_unevencolorbars(lon=xx,lat=yy, data=median_reaction_window_0_25,
                  uneven_levels=uneven_levels1,name=f' {variablename} {size}',units=f'{cbarunits}',i=1,j=0)

easy_plot_unevencolorbars(lon=xx,lat=yy, data=median_hot_reaction_window_0_25,
                  uneven_levels=uneven_levels1,name=f' Hot Season {variablename} {size}',units=f'{cbarunits}',i=1,j=1)

easy_plot_unevencolorbars(lon=xx,lat=yy, data=median_cold_reaction_window_0_25,
                  uneven_levels=uneven_levels1,name= f' Cold Season {variablename} {size}',units=f'{cbarunits}',i=1,j=2)

easy_plot_unevencolorbars2(lon=xx,lat=yy, data=median_hot_reaction_window_0_25-median_cold_reaction_window_0_25,
                  uneven_levels=uneven_levels2,name=f'Difference {size}',units=f'{cbarunits}',i=1,j=3)

size = '25-700 sqr.degs.'


easy_plot_unevencolorbars(lon=xx,lat=yy, data=median_reaction_window_25_700,
                  uneven_levels=uneven_levels1,name=f' {variablename} {size}',units=f'{cbarunits}',i=2,j=0)

easy_plot_unevencolorbars(lon=xx,lat=yy, data=median_hot_reaction_window_25_700,
                  uneven_levels=uneven_levels1,name=f' Hot Season {variablename} {size}',units=f'{cbarunits}',i=2,j=1)

easy_plot_unevencolorbars(lon=xx,lat=yy, data=median_cold_reaction_window_25_700,
                  uneven_levels=uneven_levels1,name= f' Cold Season {variablename} {size}',units=f'{cbarunits}',i=2,j=2)

easy_plot_unevencolorbars2(lon=xx,lat=yy, data=median_hot_reaction_window_25_700-median_cold_reaction_window_25_700,
                  uneven_levels=uneven_levels2,name=f'Difference {size}',units=f'{cbarunits}',i=2,j=3)

plt.tight_layout()
plt.savefig(f'/home/shilpa/glory_mat_analysis/compare_median_reactionwindow_hot_cold_seasons2.png')

# plot
nrows=3
ncols=4

cbarunits = 'Days'
variablename = 'Mean Coping Window'

fig, axs = plt.subplots(nrows=nrows,ncols=ncols,subplot_kw={'projection': proj_crs},figsize=(22,15))

uneven_levels1 = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,40,45,50,55,60,65,70,80,90,100,150,200,250,300,350,600]

uneven_levels2 = [-600,-350,-100,-50,-10,-5,0,5,10,50,100,350,600]

size = '0-700 sqr.degs.'

easy_plot_unevencolorbars(lon=xx,lat=yy, data=mean_coping_window_0_700,
                  uneven_levels=uneven_levels1,name=f' {variablename} {size}',units=f'{cbarunits}',i=0,j=0)

easy_plot_unevencolorbars(lon=xx,lat=yy, data=mean_hot_coping_window_0_700,
                  uneven_levels=uneven_levels1,name=f' Hot Season {variablename} {size}',units=f'{cbarunits}',i=0,j=1)

easy_plot_unevencolorbars(lon=xx,lat=yy, data=mean_cold_coping_window_0_700,
                  uneven_levels=uneven_levels1,name= f' Cold Season {variablename} {size}',units=f'{cbarunits}',i=0,j=2)

easy_plot_unevencolorbars2(lon=xx,lat=yy, data=mean_hot_coping_window_0_700-mean_cold_coping_window_0_700,
                  uneven_levels=uneven_levels2,name=f'Difference {size}',units=f'{cbarunits}',i=0,j=3)

size = '0-25 sqr.degs.'

easy_plot_unevencolorbars(lon=xx,lat=yy, data=mean_coping_window_0_25,
                  uneven_levels=uneven_levels1,name=f' {variablename} {size}',units=f'{cbarunits}',i=1,j=0)

easy_plot_unevencolorbars(lon=xx,lat=yy, data=mean_hot_coping_window_0_25,
                  uneven_levels=uneven_levels1,name=f' Hot Season {variablename} {size}',units=f'{cbarunits}',i=1,j=1)

easy_plot_unevencolorbars(lon=xx,lat=yy, data=mean_cold_coping_window_0_25,
                  uneven_levels=uneven_levels1,name= f' Cold Season {variablename} {size}',units=f'{cbarunits}',i=1,j=2)

easy_plot_unevencolorbars2(lon=xx,lat=yy, data=mean_hot_coping_window_0_25-mean_cold_coping_window_0_25,
                  uneven_levels=uneven_levels2,name=f'Difference {size}',units=f'{cbarunits}',i=1,j=3)

size = '25-700 sqr.degs.'

easy_plot_unevencolorbars(lon=xx,lat=yy, data=mean_coping_window_25_700,
                  uneven_levels=uneven_levels1,name=f' {variablename} {size}',units=f'{cbarunits}',i=2,j=0)

easy_plot_unevencolorbars(lon=xx,lat=yy, data=mean_hot_coping_window_25_700,
                  uneven_levels=uneven_levels1,name=f' Hot Season {variablename} {size}',units=f'{cbarunits}',i=2,j=1)

easy_plot_unevencolorbars(lon=xx,lat=yy, data=mean_cold_coping_window_25_700,
                  uneven_levels=uneven_levels1,name= f' Cold Season {variablename} {size}',units=f'{cbarunits}',i=2,j=2)

easy_plot_unevencolorbars2(lon=xx,lat=yy, data=mean_hot_coping_window_25_700-mean_cold_coping_window_25_700,
                  uneven_levels=uneven_levels2,name=f'Difference {size}',units=f'{cbarunits}',i=2,j=3)

plt.tight_layout()
plt.savefig(f'/home/shilpa/glory_mat_analysis/compare_mean_copingwindow_hot_cold_seasons2.png')

# plot
nrows=3
ncols=4

cbarunits = 'Days'
variablename = 'Median Coping Window'

fig, axs = plt.subplots(nrows=nrows,ncols=ncols,subplot_kw={'projection': proj_crs},figsize=(22,15))

uneven_levels1 = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,40,45,50,55,60,65,70,80,90,100,150,200,250,300,350,600]

uneven_levels2 = [-600,-350,-100,-50,-10,-5,0,5,10,50,100,350,600]

size = '0-700 sqr.degs.'

easy_plot_unevencolorbars(lon=xx,lat=yy, data=median_coping_window_0_700,
                  uneven_levels=uneven_levels1,name=f' {variablename} {size}',units=f'{cbarunits}',i=0,j=0)

easy_plot_unevencolorbars(lon=xx,lat=yy, data=median_hot_coping_window_0_700,
                  uneven_levels=uneven_levels1,name=f' Hot Season {variablename} {size}',units=f'{cbarunits}',i=0,j=1)

easy_plot_unevencolorbars(lon=xx,lat=yy, data=median_cold_coping_window_0_700,
                  uneven_levels=uneven_levels1,name= f' Cold Season {variablename} {size}',units=f'{cbarunits}',i=0,j=2)

easy_plot_unevencolorbars2(lon=xx,lat=yy, data=median_hot_coping_window_0_700-median_cold_coping_window_0_700,
                  uneven_levels=uneven_levels2,name=f'Difference {size}',units=f'{cbarunits}',i=0,j=3)

size = '0-25 sqr.degs.'

easy_plot_unevencolorbars(lon=xx,lat=yy, data=median_coping_window_0_25,
                  uneven_levels=uneven_levels1,name=f' {variablename} {size}',units=f'{cbarunits}',i=1,j=0)

easy_plot_unevencolorbars(lon=xx,lat=yy, data=median_hot_coping_window_0_25,
                  uneven_levels=uneven_levels1,name=f' Hot Season {variablename} {size}',units=f'{cbarunits}',i=1,j=1)

easy_plot_unevencolorbars(lon=xx,lat=yy, data=median_cold_coping_window_0_25,
                  uneven_levels=uneven_levels1,name= f' Cold Season {variablename} {size}',units=f'{cbarunits}',i=1,j=2)

easy_plot_unevencolorbars2(lon=xx,lat=yy, data=median_hot_coping_window_0_25-median_cold_coping_window_0_25,
                  uneven_levels=uneven_levels2,name=f'Difference {size}',units=f'{cbarunits}',i=1,j=3)

size = '25-700 sqr.degs.'

easy_plot_unevencolorbars(lon=xx,lat=yy, data=median_coping_window_25_700,
                  uneven_levels=uneven_levels1,name=f' {variablename} {size}',units=f'{cbarunits}',i=2,j=0)

easy_plot_unevencolorbars(lon=xx,lat=yy, data=median_hot_coping_window_25_700,
                  uneven_levels=uneven_levels1,name=f' Hot Season {variablename} {size}',units=f'{cbarunits}',i=2,j=1)

easy_plot_unevencolorbars(lon=xx,lat=yy, data=median_cold_coping_window_25_700,
                  uneven_levels=uneven_levels1,name= f' Cold Season {variablename} {size}',units=f'{cbarunits}',i=2,j=2)

easy_plot_unevencolorbars2(lon=xx,lat=yy, data=median_hot_coping_window_25_700-median_cold_coping_window_25_700,
                  uneven_levels=uneven_levels2,name=f'Difference {size}',units=f'{cbarunits}',i=2,j=3)

plt.tight_layout()
plt.savefig(f'/home/shilpa/glory_mat_analysis/compare_median_copingwindow_hot_cold_seasons2.png')



# plot
nrows=3
ncols=4

cbarunits = 'Days'
variablename = 'Mean Recovery Window'

fig, axs = plt.subplots(nrows=nrows,ncols=ncols,subplot_kw={'projection': proj_crs},figsize=(22,15))

uneven_levels1 = [1,14,30,60,100,140,180,220,270,300,320,360,400,600]

uneven_levels2 = [-600,-400,-350,-200,-100,-50,-10,-5,0,5,10,50,100,200,350,400,600]

size = '0-700 sqr.degs.'

easy_plot_unevencolorbars(lon=xx,lat=yy, data=mean_recovery_window_0_700,
                  uneven_levels=uneven_levels1,name=f' {variablename} {size}',units=f'{cbarunits}',i=0,j=0)

easy_plot_unevencolorbars(lon=xx,lat=yy, data=mean_hot_recovery_window_0_700,
                  uneven_levels=uneven_levels1,name=f' Hot Season {variablename} {size}',units=f'{cbarunits}',i=0,j=1)

easy_plot_unevencolorbars(lon=xx,lat=yy, data=mean_cold_recovery_window_0_700,
                  uneven_levels=uneven_levels1,name= f' Cold Season {variablename} {size}',units=f'{cbarunits}',i=0,j=2)

easy_plot_unevencolorbars2(lon=xx,lat=yy, data=mean_hot_recovery_window_0_700-mean_cold_recovery_window_0_700,
                  uneven_levels=uneven_levels2,name=f'Difference {size}',units=f'{cbarunits}',i=0,j=3)

size = '0-25 sqr.degs.'

easy_plot_unevencolorbars(lon=xx,lat=yy, data=mean_recovery_window_0_25,
                  uneven_levels=uneven_levels1,name=f' {variablename} {size}',units=f'{cbarunits}',i=1,j=0)

easy_plot_unevencolorbars(lon=xx,lat=yy, data=mean_hot_recovery_window_0_25,
                  uneven_levels=uneven_levels1,name=f' Hot Season {variablename} {size}',units=f'{cbarunits}',i=1,j=1)

easy_plot_unevencolorbars(lon=xx,lat=yy, data=mean_cold_recovery_window_0_25,
                  uneven_levels=uneven_levels1,name= f' Cold Season {variablename} {size}',units=f'{cbarunits}',i=1,j=2)

easy_plot_unevencolorbars2(lon=xx,lat=yy, data=mean_hot_recovery_window_0_25-mean_cold_recovery_window_0_25,
                  uneven_levels=uneven_levels2,name=f'Difference {size}',units=f'{cbarunits}',i=1,j=3)

size = '25-700 sqr.degs.'

easy_plot_unevencolorbars(lon=xx,lat=yy, data=mean_recovery_window_25_700,
                  uneven_levels=uneven_levels1,name=f' {variablename} {size}',units=f'{cbarunits}',i=2,j=0)

easy_plot_unevencolorbars(lon=xx,lat=yy, data=mean_hot_recovery_window_25_700,
                  uneven_levels=uneven_levels1,name=f' Hot Season {variablename} {size}',units=f'{cbarunits}',i=2,j=1)

easy_plot_unevencolorbars(lon=xx,lat=yy, data=mean_cold_recovery_window_25_700,
                  uneven_levels=uneven_levels1,name= f' Cold Season {variablename} {size}',units=f'{cbarunits}',i=2,j=2)

easy_plot_unevencolorbars2(lon=xx,lat=yy, data=mean_hot_recovery_window_25_700-mean_cold_recovery_window_25_700,
                  uneven_levels=uneven_levels2,name=f'Difference {size}',units=f'{cbarunits}',i=2,j=3)

plt.tight_layout()
plt.savefig(f'/home/shilpa/glory_mat_analysis/compare_mean_recoverywindow_hot_cold_seasons2.png')



# plot
nrows=3
ncols=4

cbarunits = 'Days'
variablename = 'Median Recovery Window'

fig, axs = plt.subplots(nrows=nrows,ncols=ncols,subplot_kw={'projection': proj_crs},figsize=(22,15))

uneven_levels1 = [1,14,30,60,100,140,180,220,270,300,320,360,400,600]

uneven_levels2 = [-600,-400,-350,-200,-100,-50,-10,-5,0,5,10,50,100,200,350,400,600]



size = '0-700 sqr.degs.'

easy_plot_unevencolorbars(lon=xx,lat=yy, data=median_recovery_window_0_700,
                  uneven_levels=uneven_levels1,name=f' {variablename} {size}',units=f'{cbarunits}',i=0,j=0)

easy_plot_unevencolorbars(lon=xx,lat=yy, data=median_hot_recovery_window_0_700,
                  uneven_levels=uneven_levels1,name=f' Hot Season {variablename} {size}',units=f'{cbarunits}',i=0,j=1)

easy_plot_unevencolorbars(lon=xx,lat=yy, data=median_cold_recovery_window_0_700,
                  uneven_levels=uneven_levels1,name= f' Cold Season {variablename} {size}',units=f'{cbarunits}',i=0,j=2)

easy_plot_unevencolorbars2(lon=xx,lat=yy, data=median_hot_recovery_window_0_700-median_cold_recovery_window_0_700,
                  uneven_levels=uneven_levels2,name=f'Difference {size}',units=f'{cbarunits}',i=0,j=3)

size = '0-25 sqr.degs.'

easy_plot_unevencolorbars(lon=xx,lat=yy, data=median_recovery_window_0_25,
                  uneven_levels=uneven_levels1,name=f' {variablename} {size}',units=f'{cbarunits}',i=1,j=0)

easy_plot_unevencolorbars(lon=xx,lat=yy, data=median_hot_recovery_window_0_25,
                  uneven_levels=uneven_levels1,name=f' Hot Season {variablename} {size}',units=f'{cbarunits}',i=1,j=1)

easy_plot_unevencolorbars(lon=xx,lat=yy, data=median_cold_recovery_window_0_25,
                  uneven_levels=uneven_levels1,name= f' Cold Season {variablename} {size}',units=f'{cbarunits}',i=1,j=2)

easy_plot_unevencolorbars2(lon=xx,lat=yy, data=median_hot_recovery_window_0_25-median_cold_recovery_window_0_25,
                  uneven_levels=uneven_levels2,name=f'Difference {size}',units=f'{cbarunits}',i=1,j=3)

size = '25-700 sqr.degs.'

easy_plot_unevencolorbars(lon=xx,lat=yy, data=median_recovery_window_25_700,
                  uneven_levels=uneven_levels1,name=f' {variablename} {size}',units=f'{cbarunits}',i=2,j=0)

easy_plot_unevencolorbars(lon=xx,lat=yy, data=median_hot_recovery_window_25_700,
                  uneven_levels=uneven_levels1,name=f' Hot Season {variablename} {size}',units=f'{cbarunits}',i=2,j=1)

easy_plot_unevencolorbars(lon=xx,lat=yy, data=median_cold_recovery_window_25_700,
                  uneven_levels=uneven_levels1,name= f' Cold Season {variablename} {size}',units=f'{cbarunits}',i=2,j=2)

easy_plot_unevencolorbars2(lon=xx,lat=yy, data=median_hot_recovery_window_25_700-median_cold_recovery_window_25_700,
                  uneven_levels=uneven_levels2,name=f'Difference {size}',units=f'{cbarunits}',i=2,j=3)

plt.tight_layout()
plt.savefig(f'/home/shilpa/glory_mat_analysis/compare_median_recoverywindow_hot_cold_seasons2.png')




# plot
nrows=3
ncols=4

cbarunits = 'Probability of occurrence'
variablename = 'All MHW days'

fig, axs = plt.subplots(nrows=nrows,ncols=ncols,subplot_kw={'projection': proj_crs},figsize=(22,15))

uneven_levels1 = [0,0.005,0.0125,0.025,0.05,0.075,0.1,0.25,0.5,0.75,1]

uneven_levels2 = [-0.1,-0.075,-0.05,-0.025,-0.0125,-0.005,0,0.005,0.0125,0.025,0.05,0.075,0.1]



size = '0-700 sqr.degs.'

easy_plot_unevencolorbars(lon=xx,lat=yy, data=prob_0_700_alldays,
                  uneven_levels=uneven_levels1,name=f' {variablename} {size}',units=f'{cbarunits}',i=0,j=0)

easy_plot_unevencolorbars(lon=xx,lat=yy, data=hot_prob_0_700_alldays,
                  uneven_levels=uneven_levels1,name=f' Hot Season {variablename} {size}',units=f'{cbarunits}',i=0,j=1)

easy_plot_unevencolorbars(lon=xx,lat=yy, data=cold_prob_0_700_alldays,
                  uneven_levels=uneven_levels1,name= f' Cold Season {variablename} {size}',units=f'{cbarunits}',i=0,j=2)

easy_plot_unevencolorbars2(lon=xx,lat=yy, data=hot_prob_0_700_alldays-cold_prob_0_700_alldays,
                  uneven_levels=uneven_levels2,name=f'Difference {size}',units=f'{cbarunits}',i=0,j=3)

size = '0-25 sqr.degs.'

easy_plot_unevencolorbars(lon=xx,lat=yy, data=prob_0_25_alldays,
                  uneven_levels=uneven_levels1,name=f' {variablename} {size}',units=f'{cbarunits}',i=1,j=0)

easy_plot_unevencolorbars(lon=xx,lat=yy, data=hot_prob_0_25_alldays,
                  uneven_levels=uneven_levels1,name=f' Hot Season {variablename} {size}',units=f'{cbarunits}',i=1,j=1)

easy_plot_unevencolorbars(lon=xx,lat=yy, data=cold_prob_0_25_alldays,
                  uneven_levels=uneven_levels1,name= f' Cold Season {variablename} {size}',units=f'{cbarunits}',i=1,j=2)

easy_plot_unevencolorbars2(lon=xx,lat=yy, data=hot_prob_0_25_alldays-cold_prob_0_25_alldays,
                  uneven_levels=uneven_levels2,name=f'Difference {size}',units=f'{cbarunits}',i=1,j=3)

size = '25-700 sqr.degs.'

easy_plot_unevencolorbars(lon=xx,lat=yy, data=prob_25_700_alldays,
                  uneven_levels=uneven_levels1,name=f' {variablename} {size}',units=f'{cbarunits}',i=2,j=0)

easy_plot_unevencolorbars(lon=xx,lat=yy, data=hot_prob_25_700_alldays,
                  uneven_levels=uneven_levels1,name=f' Hot Season {variablename} {size}',units=f'{cbarunits}',i=2,j=1)

easy_plot_unevencolorbars(lon=xx,lat=yy, data=cold_prob_25_700_alldays,
                  uneven_levels=uneven_levels1,name= f' Cold Season {variablename} {size}',units=f'{cbarunits}',i=2,j=2)

easy_plot_unevencolorbars2(lon=xx,lat=yy, data=hot_prob_25_700_alldays-cold_prob_25_700_alldays,
                  uneven_levels=uneven_levels2,name=f'Difference {size}',units=f'{cbarunits}',i=2,j=3)

plt.tight_layout()
plt.savefig(f'/home/shilpa/glory_mat_analysis/prob_of_occurence_all_days_hot_cold.png')


# plot
nrows=3
ncols=4

cbarunits = 'Probability of occurence'
variablename = 'Only MHW Days'

fig, axs = plt.subplots(nrows=nrows,ncols=ncols,subplot_kw={'projection': proj_crs},figsize=(22,15))

uneven_levels1 = [0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1]

uneven_levels2 = [-1,-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1]



size = '0-700 sqr.degs.'

easy_plot_unevencolorbars(lon=xx,lat=yy, data=prob_0_700_mhwdays,
                  uneven_levels=uneven_levels1,name=f' {variablename} {size}',units=f'{cbarunits}',i=0,j=0)

easy_plot_unevencolorbars(lon=xx,lat=yy, data=hot_prob_0_700_mhwdays,
                  uneven_levels=uneven_levels1,name=f' Hot Season {variablename} {size}',units=f'{cbarunits}',i=0,j=1)

easy_plot_unevencolorbars(lon=xx,lat=yy, data=cold_prob_0_700_mhwdays,
                  uneven_levels=uneven_levels1,name= f' Cold Season {variablename} {size}',units=f'{cbarunits}',i=0,j=2)

easy_plot_unevencolorbars2(lon=xx,lat=yy, data=hot_prob_0_700_mhwdays-cold_prob_0_700_mhwdays,
                  uneven_levels=uneven_levels2,name=f'Difference {size}',units=f'{cbarunits}',i=0,j=3)

size = '0-25 sqr.degs.'

easy_plot_unevencolorbars(lon=xx,lat=yy, data=prob_0_25_mhwdays,
                  uneven_levels=uneven_levels1,name=f' {variablename} {size}',units=f'{cbarunits}',i=1,j=0)

easy_plot_unevencolorbars(lon=xx,lat=yy, data=hot_prob_0_25_mhwdays,
                  uneven_levels=uneven_levels1,name=f' Hot Season {variablename} {size}',units=f'{cbarunits}',i=1,j=1)

easy_plot_unevencolorbars(lon=xx,lat=yy, data=cold_prob_0_25_mhwdays,
                  uneven_levels=uneven_levels1,name= f' Cold Season {variablename} {size}',units=f'{cbarunits}',i=1,j=2)

easy_plot_unevencolorbars2(lon=xx,lat=yy, data=hot_prob_0_25_mhwdays-cold_prob_0_25_mhwdays,
                  uneven_levels=uneven_levels2,name=f'Difference {size}',units=f'{cbarunits}',i=1,j=3)

size = '25-700 sqr.degs.'

easy_plot_unevencolorbars(lon=xx,lat=yy, data=prob_25_700_mhwdays,
                  uneven_levels=uneven_levels1,name=f' {variablename} {size}',units=f'{cbarunits}',i=2,j=0)

easy_plot_unevencolorbars(lon=xx,lat=yy, data=hot_prob_25_700_mhwdays,
                  uneven_levels=uneven_levels1,name=f' Hot Season {variablename} {size}',units=f'{cbarunits}',i=2,j=1)

easy_plot_unevencolorbars(lon=xx,lat=yy, data=cold_prob_25_700_mhwdays,
                  uneven_levels=uneven_levels1,name= f' Cold Season {variablename} {size}',units=f'{cbarunits}',i=2,j=2)

easy_plot_unevencolorbars2(lon=xx,lat=yy, data=hot_prob_25_700_mhwdays-cold_prob_25_700_mhwdays,
                  uneven_levels=uneven_levels2,name=f'Difference {size}',units=f'{cbarunits}',i=2,j=3)



plt.tight_layout()
plt.savefig(f'/home/shilpa/glory_mat_analysis/compare_prob_of_occurence_all_mhw_days_hot_cold.png')


