# standard imports

from cartopy import config
import cartopy.crs as ccrs
import xarray as xr
import geopandas as gpd
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np
import pandas as pd
import os

# glorys linear original
glorys_linear = pd.read_csv('/home/shilpa/glory_mat_analysis/noaaoisst_lastday_included/combined_linear_glorys.csv')
sorted_gl = glorys_linear.sort_values('time_index')
sorted_gl['time'] = pd.to_datetime(sorted_gl['time'])

# noaa full
noaa_1981_2023 = pd.read_csv('/home/shilpa/glory_mat_analysis/means_std_median_polyareas_noaaoisst/noaa_%eez_spatial_extent.csv')
sorted_n_1981_2023 = noaa_1981_2023.sort_values('time_index')
sorted_n_1981_2023['time'] = pd.to_datetime(sorted_n_1981_2023['time'])

# noaa full macroscale
noaa_macro_scale = pd.read_csv('/home/shilpa/glory_mat_analysis/means_std_median_polyareas_noaaoisst/noaaoisst_mean_std_mhw_prop/%eez_spatialextent_25_more.csv')
noaa_macro_scale_1981_2023 = noaa_macro_scale.sort_values('time_index')
noaa_macro_scale_1981_2023['time'] = pd.to_datetime(noaa_macro_scale_1981_2023['time'])


# glorys full macroscale
glorys_macro_scale = pd.read_csv('/media/shilpa/Expansion/glorys_spatial_macroscale.csv')
glorys_macro_scale_1993_2019 = glorys_macro_scale.sort_values('time_index')
glorys_macro_scale_1993_2019['time'] = pd.to_datetime(glorys_macro_scale_1993_2019['time'])


# replace zeros with nans
df_replacedgl = sorted_gl.replace(0, np.nan)
df_replacedn_1981_2023 = sorted_n_1981_2023.replace(0, np.nan)
noaa_macro_scale_1981_2023_replaced = noaa_macro_scale_1981_2023.replace(0,np.nan)
glorys_macro_scale_1993_2019_replaced = glorys_macro_scale_1993_2019.replace(0,np.nan)

# lets get the coastal information
directory = '/media/shilpa/Expansion/coastal_mhw'
coastal_files = [file for file in os.listdir(directory) if file.startswith("Wallis and Futuna_spatial")]

# Read all CSV files and concatenate into a single DataFrame
dfs = []
for file in coastal_files:
    # Handle file names with spaces by adding double quotes around the file path
    file_path = os.path.join(directory, file)
    df = pd.read_csv(file_path)
    dfs.append(df)

# Concatenate all DataFrames
concatenated_df = pd.concat(dfs, ignore_index=True)

# Assuming the column containing presence and absence data is named 'presence'
# Group by 'time' (assuming it's the column name) and calculate the percentage of presence for each day
coastal_percentage = concatenated_df.groupby('time')['spatial_extent'].mean() * 100


presence_df = pd.DataFrame({'time': coastal_percentage.index, 'coastal_percentage': coastal_percentage.values})


start_date = '1993-01-01'
end_date = '2023-10-24'

# Create datetime range
date_range = pd.date_range(start=start_date, end=end_date)
print(len(date_range))


presence_df = pd.DataFrame({'time': coastal_percentage.index,'time_dt':date_range.values, 'coastal_percentage': coastal_percentage.values})


coastal = presence_df.replace(0,np.nan)

fig, axs = plt.subplots(6, 1,sharex=True, figsize=(19, 12))# constrained_layout=True)

ax = axs[0]
ax.plot(df_replacedgl.time,df_replacedgl['Wallis_and_Futuna_%spatial_extent'], color='cyan',alpha = 0.7,linestyle='solid',label='% Wallis and Futuna GLORYS [1993-2019] ')
ax.plot(df_replacedn_1981_2023.time,df_replacedn_1981_2023['Wallis_and_Futuna_%spatial_extent'], color='black',alpha = 0.5,linestyle='solid',label='% Wallis and Futuna NOAAOISST [1981-2023]')
ax.plot(noaa_macro_scale_1981_2023_replaced.time,noaa_macro_scale_1981_2023_replaced['Wallis_and_Futuna_%spatial_extent'], color='magenta',alpha = 0.5,linestyle='solid',label='% Wallis and Futuna NOAAOISST Macroscale [1981-2023]')
ax.plot(glorys_macro_scale_1993_2019_replaced.time,glorys_macro_scale_1993_2019_replaced['Wallis_and_Futuna_%spatial_extent'], color='green',alpha = 0.6,linestyle='solid',label='% Wallis and Futuna GLORYS Macroscale [1993-2019]')
ax.plot(coastal.time_dt,coastal['coastal_percentage'], color='blue',alpha = 0.4,linestyle='dashed',label='% Wallis and Futuna coastal GLORYS [1993-2023]')
ax.legend()
ax.set_ylim(0, 100)

ax = axs[1]
ax.plot(df_replacedgl.time,df_replacedgl['Wallis_and_Futuna_%spatial_extent'], color='cyan',alpha = 0.7,linestyle='solid',label='% Wallis and Futuna GLORYS [1993-2019] ')
ax.plot(df_replacedn_1981_2023.time,df_replacedn_1981_2023['Wallis_and_Futuna_%spatial_extent'], color='black',alpha = 0.5,linestyle='solid',label='% Wallis and Futuna NOAAOISST [1981-2023]')
#ax.plot(noaa_macro_scale_1981_2023_replaced.time,noaa_macro_scale_1981_2023_replaced['New_Caledonia_%spatial_extent'], color='magenta',alpha = 0.5,linestyle='solid',label='% New Caledonia NOAAOISST Macroscale [1981-2023]')
#ax.plot(glorys_macro_scale_1993_2019_replaced.time,glorys_macro_scale_1993_2019_replaced['New_Caledonia_%spatial_extent'], color='green',alpha = 0.5,linestyle='solid',label='% New Caledonia GLORYS Macroscale [1993-2019]')
#ax.plot(coastal.time_dt,coastal['coastal_percentage'], color='blue',alpha = 0.5,linestyle='dashed',label='% New Caledonia coastal GLORYS [1993-2023]')
ax.legend()
ax.set_ylim(0, 100)

ax = axs[2]
#ax.plot(df_replacedgl.time,df_replacedgl['New_Caledonia_%spatial_extent'], color='cyan',alpha = 0.7,linestyle='solid',label='% New Caledonia GLORYS [1993-2019] ')
#ax.plot(df_replacedn_1981_2023.time,df_replacedn_1981_2023['New_Caledonia_%spatial_extent'], color='black',alpha = 0.5,linestyle='solid',label='% New Caledonia NOAAOISST [1981-2023]')
ax.plot(noaa_macro_scale_1981_2023_replaced.time,noaa_macro_scale_1981_2023_replaced['Wallis_and_Futuna_%spatial_extent'], color='magenta',alpha = 0.5,linestyle='solid',label='% Wallis and Futuna NOAAOISST Macroscale [1981-2023]')
ax.plot(glorys_macro_scale_1993_2019_replaced.time,glorys_macro_scale_1993_2019_replaced['Wallis_and_Futuna_%spatial_extent'], color='green',alpha = 0.6,linestyle='solid',label='% Wallis and Futuna GLORYS Macroscale [1993-2019]')
#ax.plot(coastal.time_dt,coastal['coastal_percentage'], color='blue',alpha = 0.5,linestyle='dashed',label='% New Caledonia coastal GLORYS [1993-2023]')
ax.legend()
ax.set_ylim(0, 100)

ax = axs[3]
#ax.plot(df_replacedgl.time,df_replacedgl['New_Caledonia_%spatial_extent'], color='cyan',alpha = 0.7,linestyle='solid',label='% New Caledonia GLORYS [1993-2019] ')
#ax.plot(df_replacedn_1981_2023.time,df_replacedn_1981_2023['New_Caledonia_%spatial_extent'], color='black',alpha = 0.5,linestyle='solid',label='% New Caledonia NOAAOISST [1981-2023]')
ax.plot(noaa_macro_scale_1981_2023_replaced.time,noaa_macro_scale_1981_2023_replaced['Wallis_and_Futuna_%spatial_extent'], color='magenta',alpha = 0.5,linestyle='solid',label='% Wallis and Futuna NOAAOISST Macroscale [1981-2023]')
ax.plot(glorys_macro_scale_1993_2019_replaced.time,glorys_macro_scale_1993_2019_replaced['Wallis_and_Futuna_%spatial_extent'], color='green',alpha = 0.6,linestyle='solid',label='% Wallis and Futuna GLORYS Macroscale [1993-2019]')
ax.plot(coastal.time_dt,coastal['coastal_percentage'], color='blue',alpha = 0.4,linestyle='dashed',label='% Wallis and Futuna coastal GLORYS [1993-2023]')
ax.legend()
ax.set_ylim(0, 100)


ax = axs[4]
#ax.plot(df_replacedgl.time,df_replacedgl['New_Caledonia_%spatial_extent'], color='cyan',alpha = 0.7,linestyle='solid',label='% New Caledonia GLORYS [1993-2019] ')
#ax.plot(df_replacedn_1981_2023.time,df_replacedn_1981_2023['New_Caledonia_%spatial_extent'], color='black',alpha = 0.5,linestyle='solid',label='% New Caledonia NOAAOISST [1981-2023]')
#ax.plot(noaa_macro_scale_1981_2023_replaced.time,noaa_macro_scale_1981_2023_replaced['New_Caledonia_%spatial_extent'], color='magenta',alpha = 0.5,linestyle='solid',label='% New Caledonia NOAAOISST Macroscale [1981-2023]')
ax.plot(glorys_macro_scale_1993_2019_replaced.time,glorys_macro_scale_1993_2019_replaced['Wallis_and_Futuna_%spatial_extent'], color='green',alpha = 0.6,linestyle='solid',label='% Wallis and Futuna GLORYS Macroscale [1993-2019]')
ax.plot(coastal.time_dt,coastal['coastal_percentage'], color='blue',alpha = 0.4,linestyle='dashed',label='% Wallis and Futuna coastal GLORYS [1993-2023]')
ax.legend()
ax.set_ylim(0, 100)

ax = axs[5]
#ax.plot(df_replacedgl.time,df_replacedgl['New_Caledonia_%spatial_extent'], color='cyan',alpha = 0.7,linestyle='solid',label='% New Caledonia GLORYS [1993-2019] ')
#ax.plot(df_replacedn_1981_2023.time,df_replacedn_1981_2023['New_Caledonia_%spatial_extent'], color='black',alpha = 0.5,linestyle='solid',label='% New Caledonia NOAAOISST [1981-2023]')
ax.plot(noaa_macro_scale_1981_2023_replaced.time,noaa_macro_scale_1981_2023_replaced['Wallis_and_Futuna_%spatial_extent'], color='magenta',alpha = 0.5,linestyle='solid',label='% Wallis and Futuna NOAAOISST Macroscale [1981-2023]')
#ax.plot(glorys_macro_scale_1993_2019_replaced.time,glorys_macro_scale_1993_2019_replaced['New_Caledonia_%spatial_extent'], color='green',alpha = 0.5,linestyle='solid',label='% New Caledonia GLORYS Macroscale [1993-2019]')
ax.plot(coastal.time_dt,coastal['coastal_percentage'], color='blue',alpha = 0.4,linestyle='dashed',label='% Wallis and Futuna coastal GLORYS [1993-2023]')
ax.legend()
ax.set_ylim(0, 100)

fig.supylabel('% EEZ MHW spatial extent')
fig.savefig('/home/shilpa/glory_mat_analysis/coastal_points/timeseries_%eez_spatial_extent_wf.png')
