import xarray as xr
import numpy as np

#get duration

file_path_mask = xr.open_dataset('/home/datawork-lead/datarmor-only/shilpa/noaaoisst_full_mhw_areas_attempt2.nc')
file_path_mask_slice  = file_path_mask.sel(lat=slice(-34.875,-2.625),lon=slice(145.125,209.875))
spatial = file_path_mask_slice['mhw_spatial_extent_area']

def size_selector(upper_value, lower_value):
	file_path_data = '/home/datawork-lead/datarmor-only/shilpa/noaaoisst_mhw_output_1981_2023_aoi_clim_1993_2019.nc'

	ds_data = xr.open_dataset(file_path_data)
	ds2 = ds_data['duration_pacific'].where(ds_data['duration_pacific'] != 0, np.nan) #replace 0s with nans

	lower_value = lower_value  # Example lower bound
	upper_value = upper_value  # Example upper bound

	polygon_size_mask = (file_path_mask_slice['mhw_spatial_extent_area'] > lower_value) & \
	(file_path_mask_slice['mhw_spatial_extent_area'] < upper_value)
	masked_spatial = file_path_mask_slice['mhw_spatial_extent_area'].where(polygon_size_mask)
	ds2 = ds2.where(polygon_size_mask)
  
	#time_slice_data = ds2.sel(time=slice('1993-01-01','2022-12-31'))
	var_arr = np.array(ds2[:,:,:])

	processed_data = []

	for j in range(130):
		for k in range(260):
			time_series = var_arr[:,j,k]

        #  keep one value for each group of same values
			cleaned_series = []
			last_val = np.nan
			for val in time_series:
				if np.isnan(val):
					if not np.isnan(last_val):
						last_val = val
				else:
					if  val != 0 and val != last_val:
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

def size_selector_hotseason(upper_value, lower_value):
        file_path_data = '/home/datawork-lead/datarmor-only/shilpa/noaaoisst_mhw_output_1981_2023_aoi_clim_1993_2019.nc'

        ds_data = xr.open_dataset(file_path_data)
        ds2 = ds_data['duration_pacific'].where(ds_data['duration_pacific'] != 0, np.nan) #replace 0s with nans

        lower_value = lower_value  # Example lower bound
        upper_value = upper_value  # Example upper bound

        polygon_size_mask = (file_path_mask_slice['mhw_spatial_extent_area'] > lower_value) & \
        (file_path_mask_slice['mhw_spatial_extent_area'] < upper_value)
        masked_spatial = file_path_mask_slice['mhw_spatial_extent_area'].where(polygon_size_mask)
        ds2 = ds2.where(polygon_size_mask)

        #time_slice_data = ds2.sel(time=slice('1993-01-01','2022-12-31'))
        # Condition 1: November to December
        #condition1 = time_slice_data.time.dt.month >= 11
        condition1 = ds2.time.dt.month >= 11

        # Condition 2: January to April
        #condition2 = time_slice_data.time.dt.month <= 4
        condition2 = ds2.time.dt.month <= 4

        # Combine the two conditions
        combined_condition = condition1 | condition2

        # Filter the original time-sliced data using the combined condition
        filtered_data = ds2.sel(time=combined_condition)

        var_arr = np.array(filtered_data[:,:,:])

        #var_arr = np.array(ds2[:,:,:])

        processed_data = []

        for j in range(130):
                for k in range(260):
                        time_series = var_arr[:,j,k]

        #  keep one value for each group of same values
                        cleaned_series = []
                        last_val = np.nan
                        for val in time_series:
                                if np.isnan(val):
                                        if not np.isnan(last_val):
                                                last_val = val
                                else:
                                        if  val != 0 and val != last_val:
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

def size_selector_coldseason(upper_value, lower_value):
        file_path_data = '/home/datawork-lead/datarmor-only/shilpa/noaaoisst_mhw_output_1981_2023_aoi_clim_1993_2019.nc'

        ds_data = xr.open_dataset(file_path_data)
        ds2 = ds_data['duration_pacific'].where(ds_data['duration_pacific'] != 0, np.nan) #replace 0s with nans

        lower_value = lower_value  # Example lower bound
        upper_value = upper_value  # Example upper bound

        polygon_size_mask = (file_path_mask_slice['mhw_spatial_extent_area'] > lower_value) & \
        (file_path_mask_slice['mhw_spatial_extent_area'] < upper_value)
        masked_spatial = file_path_mask_slice['mhw_spatial_extent_area'].where(polygon_size_mask)
        ds2 = ds2.where(polygon_size_mask)

        #time_slice_data = ds2.sel(time=slice('1993-01-01','2022-12-31'))
        filtered_data = ds2.sel(time=ds2.time.dt.month.isin(range(5, 10))) # cold season is frm may to october

        var_arr = np.array(filtered_data[:,:,:])

        processed_data = []

	for j in range(130):
		for k in range(260):
			time_series = var_arr[:,j,k]

        #  keep one value for each group of same values
			cleaned_series = []
			last_val = np.nan
			for val in time_series:
				if np.isnan(val):
					if not np.isnan(last_val):
						last_val = val
				else:
					if  val != 0 and val != last_val:
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


values1 = [0,25,700]
values2 = [0,700]
data_vars = {}
variable = 'duration_pacific'
for i in range(len(values)-1):
	upper_value = values1[i+1]
	lower_value = values1[i]
	mean_values,std_values,median_values = size_selector(upper_value=upper_value,lower_value=lower_value)
	var_name = f"mean_{variable}_{lower_value}_{upper_value}"
	data_vars[var_name] = (['lat', 'lon'], mean_values)
	var_name = f"std_{variable}_{lower_value}_{upper_value}"
	data_vars[var_name] = (['lat', 'lon'], std_values)
	var_name = f"median_{variable}_{lower_value}_{upper_value}"
	data_vars[var_name] = (['lat', 'lon'], median_values)

	upper_value = values2[i+1]
	lower_value = values2[i]
	mean_values,std_values,median_values = size_selector(upper_value=upper_value,lower_value=lower_value)
	var_name = f"mean_{variable}_{lower_value}_{upper_value}"
	data_vars[var_name] = (['lat', 'lon'], mean_values)
	var_name = f"std_{variable}_{lower_value}_{upper_value}"
	data_vars[var_name] = (['lat', 'lon'], std_values)
	var_name = f"median_{variable}_{lower_value}_{upper_value}"
	data_vars[var_name] = (['lat', 'lon'], median_values)

	upper_value = values1[i+1]
	lower_value = values1[i]
	mean_values,std_values,median_values = size_selector_hotseason(upper_value=upper_value,lower_value=lower_value)
	var_name = f"mean_{variable}_{lower_value}_{upper_value}_hot"
	data_vars[var_name] = (['lat', 'lon'], mean_values)
	var_name = f"std_{variable}_{lower_value}_{upper_value}_hot"
	data_vars[var_name] = (['lat', 'lon'], std_values)
	var_name = f"median_{variable}_{lower_value}_{upper_value}_hot"
	data_vars[var_name] = (['lat', 'lon'], median_values)

	upper_value = values2[i+1]
	lower_value = values2[i]
	mean_values,std_values,median_values = size_selector_hotseason(upper_value=upper_value,lower_value=lower_value)
	var_name = f"mean_{variable}_{lower_value}_{upper_value}_hot"
	data_vars[var_name] = (['lat', 'lon'], mean_values)
	var_name = f"std_{variable}_{lower_value}_{upper_value}_hot"
	data_vars[var_name] = (['lat', 'lon'], std_values)
	var_name = f"median_{variable}_{lower_value}_{upper_value}_hot"
        data_vars[var_name] = (['lat', 'lon'], median_values)

	upper_value = values1[i+1]
	lower_value = values1[i]
	mean_values,std_values,median_values = size_selector_coldseason(upper_value=upper_value,lower_value=lower_value)
	var_name = f"mean_{variable}_{lower_value}_{upper_value}_cold"
	data_vars[var_name] = (['lat', 'lon'], mean_values)
	var_name = f"std_{variable}_{lower_value}_{upper_value}_cold"
	data_vars[var_name] = (['lat', 'lon'], std_values)
	var_name = f"median_{variable}_{lower_value}_{upper_value}_cold"
	data_vars[var_name] = (['lat', 'lon'], median_values)

	upper_value = values2[i+1]
	lower_value = values2[i]
	mean_values,std_values,median_values = size_selector_coldseason(upper_value=upper_value,lower_value=lower_value)
	var_name = f"mean_{variable}_{lower_value}_{upper_value}_cold"
	data_vars[var_name] = (['lat', 'lon'], mean_values)
	var_name = f"std_{variable}_{lower_value}_{upper_value}_cold"
	data_vars[var_name] = (['lat', 'lon'], std_values)
	var_name = f"median_{variable}_{lower_value}_{upper_value}_cold"
	data_vars[var_name] = (['lat', 'lon'], median_values)


longi = np.array(file_path_mask_slice.lon.values)
lati = np.array(file_path_mask_slice.lat.values)


xrds = xr.Dataset(
       coords = dict(lon = longi,lat = lati),
       data_vars = data_vars)

xrds.to_netcdf('/home/datawork-lead/datarmor-only/shilpa/noaaoisst_1981_2023_mean_std_median_duration.nc')


                                                           
