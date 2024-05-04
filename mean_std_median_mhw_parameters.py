import xarray as xr
import numpy as np

#ratedecline

file_path_data = '/home/datawork-lead/datarmor-only/shilpa/ratedecline_noaaoisst_for_polyareas.nc'

ds_data = xr.open_dataset(file_path_data)
ds2 = ds_data.where(ds_data != 0, np.nan) #replace 0s with nans

file_path_mask = xr.open_dataset('/home/datawork-lead/datarmor-only/shilpa/noaaoisst_full_mhw_areas_attempt2.nc')
file_path_mask_slice  = file_path_mask.sel(lat=slice(-34.875,-2.625),lon=slice(145.125,209.875))

#polygon_size_mask = file_path_mask_slice['mhw_spatial_extent_area'] > 50

#ds2['duration_pacific'] = ds2['duration_pacific'].where(polygon_size_mask)

time_slice_data = ds2.sel(time=slice('1993-01-01','2022-12-31'))

var_arr = np.array(ds2.ratedecline_arr[:,:,:])

#mean_var_arr = np.zeros((130,260))
#std_var_arr = np.zeros((130,260))

processed_data = []

for j in range(130):
    for k in range(260):
        time_series = var_arr[:,j,k]

        # Remove NaNs and keep one value for each group of same values
        cleaned_series = []
        last_val = None
        for val in time_series:
            if np.isnan(val):
                continue
            if val != last_val:
                cleaned_series.append(val)
                last_val = val

        processed_data.append(cleaned_series)

# Convert to a NumPy array for easier mean and std deviation calculation
processed_data_np = np.array([np.nanmean(group) if group else np.nan for group in processed_data])
processed_data_std = np.array([np.nanstd(group) if group else np.nan for group in processed_data])

# Reshape back to lat-lon dimensions for mean and std deviation
mean_values = processed_data_np.reshape(130, 260)
std_values = processed_data_std.reshape(130, 260)

def size_selector(value):
	file_path_data = '/home/datawork-lead/datarmor-only/shilpa/ratedecline_noaaoisst_for_polyareas.nc'

	ds_data = xr.open_dataset(file_path_data)
	ds2 = ds_data.where(ds_data != 0, np.nan) #replace 0s with nans

	polygon_size_mask = file_path_mask_slice['mhw_spatial_extent_area'] > value

	ds2['ratedecline_arr'] = ds2['ratedecline_arr'].where(polygon_size_mask)

	time_slice_data = ds2.sel(time=slice('1993-01-01','2022-12-31'))

	var_arr = np.array(ds2.ratedecline_arr[:,:,:])

	processed_data = []

	for j in range(130):
    		for k in range(260):
        		time_series = var_arr[:,j,k]

        # Remove NaNs and keep one value for each group of same values
        		cleaned_series = []
        		last_val = None
        		for val in time_series:
            			if np.isnan(val):
                			continue
            			if val != last_val:
                			cleaned_series.append(val)
                			last_val = val

        		processed_data.append(cleaned_series)

			
# Convert to a NumPy array for easier mean and std deviation calculation
	processed_data_mean = np.array([np.nanmean(group) if group else np.nan for group in processed_data])
	processed_data_std = np.array([np.nanstd(group) if group else np.nan for group in processed_data])
	processed_data_median = np.array([np.nanmedian(group) if group else np.nan for group in processed_data])

# Reshape back to lat-lon dimensions for mean and std deviation
	mean_values = processed_data_mean.reshape(130, 260)
	std_values = processed_data_std.reshape(130, 260)
	median_values = processed_data_median.reshape(130, 260)

    return mean_values, std_values, median_values

mean_values_0_1_rem,std_values_0_1_rem,median_values_0_1_rem = size_selector(1)
mean_values_0_10_rem,std_values_0_10_rem,median_values_0_10_rem = size_selector(10)
mean_values_0_25_rem,std_values_0_25_rem,median_values_0_25_rem = size_selector(25)
mean_values_0_50_rem,std_values_0_50_rem,median_values_0_50_rem = size_selector(50)
mean_values_0_100_rem,std_values_0_100_rem,median_values_0_100_rem = size_selector(100)




longi = np.array(ds2.lon.values)
lati = np.array(ds2.lat.values)


xrds = xr.Dataset(
       coords = dict(lon = longi,lat = lati),
       data_vars = dict(mean_ratedecline = (['lat','lon'], mean_values),std_ratedecline = (['lat','lon'], std_values),median_ratedecline = (['lat','lon'], median_values),
                   mean_ratedecline_0_1_rem = (['lat','lon'], mean_values_0_1_rem),std_ratedecline_0_1_rem = (['lat','lon'], std_values_0_1_rem),median_ratedecline_0_1_rem = (['lat','lon'], median_values_0_1_rem),
                   mean_ratedecline_0_10_rem = (['lat','lon'], mean_values_0_10_rem),std_ratedecline_0_10_rem = (['lat','lon'], std_values_0_10_rem),median_ratedecline_0_10_rem = (['lat','lon'], median_values_0_10_rem),
                   mean_ratedecline_0_25_rem = (['lat','lon'], mean_values_0_25_rem),std_ratedecline_0_25_rem = (['lat','lon'], std_values_0_25_rem),median_ratedecline_0_25_rem = (['lat','lon'], median_values_0_25_rem),
                   mean_ratedecline_0_50_rem = (['lat','lon'], mean_values_0_50_rem),std_ratedecline_0_50_rem = (['lat','lon'], std_values_0_50_rem),median_ratedecline_0_50_rem = (['lat','lon'], median_values_0_50_rem),
                   mean_ratedecline_0_100_rem = (['lat','lon'], mean_values_0_100_rem),std_ratedecline_0_100_rem = (['lat','lon'], std_values_0_100_rem),median_ratedecline_0_100_rem = (['lat','lon'], median_values_0_100_rem)
                    ))
xrds.to_netcdf('/home/datawork-lead/datarmor-only/shilpa/polyareas_removed_noaa_93_22_mean_std_median_ratedecline.nc')


