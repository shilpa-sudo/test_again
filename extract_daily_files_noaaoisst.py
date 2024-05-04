#!/usr/bin/env python
# coding: utf-8

# In[1]:


from datetime import datetime, timedelta
import xarray as xr
import os


# In[2]:


#from datetime import datetime, timedelta

# Define the start and end dates
start_date_str = "1981-09-01"
end_date_str = "2023-06-26"

# Convert the start and end date strings to datetime objects
start_date = datetime.strptime(start_date_str, "%Y-%m-%d")
end_date = datetime.strptime(end_date_str, "%Y-%m-%d")

# Create an empty list to store the dates
date_list = []

# Generate dates within the specified range and add them to the list
current_date = start_date
while current_date <= end_date:
    date_list.append(current_date.strftime("%Y-%m-%d"))
    current_date += timedelta(days=1)

# Print the list of dates
print(date_list)


# In[ ]:


#import xarray as xr
#import os

# List of dates you want to extract (replace with your dates)
dates = date_list

# Input NetCDF file
input_file = "/home/datawork-lead/datarmor-only/shilpa/noaaoisst_mhw_output_1981_2023_aoi_clim_1993_2019.nc"

# Output directory (create if it doesn't exist)
output_dir = "/home/datawork-lead/datarmor-only/shilpa/dailyfiles_noaaoisst_1981_2023_aoi_clim_1993_2019"
os.makedirs(output_dir, exist_ok=True)

# Open the NetCDF file using xarray
ds = xr.open_dataset(input_file)

# Convert dates to a list of datetime objects
date_objs = [xr.coding.times.decode_cf_datetime(date) for date in dates]

# Iterate through the list of dates
for date_obj in date_objs:
    # Extract data for the current date
    date_str = date_obj.strftime("%Y-%m-%d")
    subset = ds.sel(time=date_str)

    if not subset.empty:
        # Save the extracted data to a new NetCDF file
        output_file = os.path.join(output_dir, f"data_{date_str}.nc")
        subset.to_netcdf(output_file)
        print(f"Extracted data for {date_str} and saved as {output_file}")
    else:
        print(f"No data found for {date_str} in the input file.")

# Close the input NetCDF file
ds.close()

