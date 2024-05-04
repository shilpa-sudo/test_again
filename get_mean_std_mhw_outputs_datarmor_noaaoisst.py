#!/usr/bin/env python
# coding: utf-8

# In[2]:


import xarray as xr
import numpy as np


# In[ ]:





# In[ ]:


ds = xr.open_dataset('/home/datawork-lead/datarmor-only/shilpa/noaaoisst_mhw_output_test_ratedecline.nc')


# In[ ]:


duration = ds.duration_pacific
intensity_max = ds.intensity_max_pacific
intensity_mean = ds.intensity_mean_pacific
intensity_cumulative = ds.intensity_cumulative_pacific
rateonset = ds.rateonset_pacific
ratedecline = ds.ratedecline_arr
reaction_window = ds.reaction_window_days_pacific
coping_window = ds.coping_window_days_pacific
recovery_window = ds.recovery_window_days_pacific
severity = ds.severity_index_pacific
event_no = ds.mhw_event_no_pacificrest_method


# In[1]:


duration_array = duration.values
intensity_max_array = intensity_max.values
intensity_mean_array = intensity_mean.values
intensity_cumulative_array = intensity_cumulative.values
rateonset_array = rateonset.values
ratedecline_array = ratedecline.values
reaction_window_array = reaction_window.values
coping_window_array = coping_window.values
recovery_window_array = recovery_window.values
severity_array = severity.values
event_no_array = event_no.values


# In[ ]:


mean_duration_arr = np.zeros((144,588))
std_duration_arr = np.zeros((144,588))

mean_intensity_max_arr = np.zeros((144,588))
std_intensity_max_arr = np.zeros((144,588))

mean_intensity_mean_arr = np.zeros((144,588))
std_intensity_mean_arr = np.zeros((144,588))

mean_intensity_cumulative_arr = np.zeros((144,588))
std_intensity_cumulative_arr = np.zeros((144,588))

mean_rateonset_arr= np.zeros((144,588)) 
std_rateonset_arr = np.zeros((144,588))

mean_ratedecline_arr = np.zeros((144,588))
std_ratedecline_arr= np.zeros((144,588)) 

mean_reaction_window_arr = np.zeros((144,588))
std_reaction_window_arr = np.zeros((144,588))

mean_coping_window_arr = np.zeros((144,588))
std_coping_window_arr = np.zeros((144,588))  

mean_recovery_window_arr = np.zeros((144,588))
std_recovery_window_arr = np.zeros((144,588))

mean_severity_arr = np.zeros((144,588))
std_severity_arr = np.zeros((144,588))

#mean_event_no_arr = np.zeros((141,604))  
#std_event_no_arr = np.zeros((141,604))  

max_event_no_arr = np.zeros((144,588)) 


# In[ ]:


for i in range(144):
    for j in range(588):
        
        nonzero_duration_indices = np.nonzero(duration_array[:,i,j])
        nonzero_duration = duration_array[nonzero_duration_indices]
        mean_duration = np.mean(nonzero_duration)
        mean_duration_arr[i,j] = mean_duration
        std_duration = np.std(nonzero_duration)
        std_duration_arr[i,j] = std_duration
        
        nonzero_intensity_max_indices = np.nonzero(intensity_max_array[:,i,j])
        nonzero_intensity_max = intensity_max_array[nonzero_intensity_max_indices]
        mean_intensity_max = np.mean(nonzero_intensity_max)
        mean_intensity_max_arr[i,j] = mean_intensity_max 
        std_intensity_max = np.std(nonzero_intensity_max)
        std_intensity_max_arr[i,j] = std_intensity_max
        
        nonzero_intensity_mean_indices = np.nonzero(intensity_mean_array[:,i,j]) 
        nonzero_intensity_mean = intensity_mean_array[nonzero_intensity_mean_indices]
        mean_intensity_mean = np.mean(nonzero_intensity_mean)
        mean_intensity_mean_arr[i,j] = mean_intensity_mean
        std_intensity_mean = np.std(nonzero_intensity_mean)
        std_intensity_mean_arr[i,j] =  std_intensity_mean 
        
        nonzero_intensity_cumulative_indices = np.nonzero(intensity_cumulative_array[:,i,j])
        nonzero_intensity_cumulative = intensity_cumulative_array[nonzero_intensity_cumulative_indices]
        mean_intensity_cumulative = np.mean(nonzero_intensity_cumulative)
        mean_intensity_cumulative_arr[i,j] = mean_intensity_cumulative
        std_intensity_cumulative = np.std(nonzero_intensity_cumulative)
        std_intensity_cumulative_arr[i,j] = std_intensity_cumulative
        
        nonzero_rateonset_indices = np.nonzero(rateonset_array[:,i,j])
        nonzero_rateonset = rateonset_array[nonzero_rateonset_indices]
        mean_rateonset = np.mean(nonzero_rateonset)
        mean_rateonset_arr[i,j] = mean_rateonset
        std_rateonset = np.std(nonzero_rateonset)
        std_rateonset_arr[i,j] = std_rateonset
        
        nonzero_ratedecline_indices = np.nonzero(ratedecline_array[:,i,j])
        nonzero_ratedecline = ratedecline_array[nonzero_ratedecline_indices]
        mean_ratedecline = np.mean(nonzero_ratedecline)
        mean_ratedecline_arr[i,j] = mean_ratedecline
        std_ratedecline = np.std(nonzero_ratedecline)
        std_ratedecline_arr[i,j] = std_ratedecline
        
        nonzero_reaction_window_indices = np.nonzero(reaction_window_array[:,i,j])
        nonzero_reaction_window = reaction_window_array[nonzero_reaction_window_indices]
        mean_reaction_window = np.mean(nonzero_reaction_window)
        mean_reaction_window_arr[i,j] = mean_reaction_window
        std_reaction_window = np.std(nonzero_reaction_window)
        std_reaction_window_arr[i,j] = std_reaction_window
        
        nonzero_coping_window_indices = np.nonzero(coping_window_array[:,i,j])
        nonzero_coping_window = coping_window_array[nonzero_coping_window_indices]
        mean_coping_window = np.mean(nonzero_coping_window)
        mean_coping_window_arr[i,j] = mean_coping_window
        std_coping_window = np.std(nonzero_coping_window)
        std_coping_window_arr[i,j] = std_coping_window
        
                                                   
        nonzero_recovery_window_indices = np.nonzero(recovery_window_array[:,i,j])
        nonzero_recovery_window = recovery_window_array[nonzero_recovery_window_indices]
        mean_recovery_window = np.mean(nonzero_recovery_window)
        mean_recovery_window_arr[i,j] = mean_recovery_window
        std_recovery_window = np.std(nonzero_recovery_window)
        std_recovery_window_arr[i,j] = std_recovery_window 
        
        nonzero_severity_indices = np.nonzero(severity_array[:,i,j])
        nonzero_severity = severity_array[nonzero_severity_indices]
        mean_severity = np.mean(nonzero_severity)
        mean_severity_arr[i,j] = mean_severity
        std_severity = np.std(nonzero_severity)
        std_severity_arr[i,j] = std_severity
        
        nonzero_event_no_indices = np.nonzero(event_no_array[:,i,j])
        nonzero_event_no = event_no_array[nonzero_event_no_indices]  
        max_event_no = np.max(nonzero_event_no)
        max_event_no_arr[i,j] = max_event_no
        #std_event_no = np.std(nonzero_event_no)
        #std_event_no_arr[i,j] = std_event_no


# In[ ]:


longi = np.array(ds.lon.values)
lati = np.array(ds.lat.values)


# In[ ]:


xrds = xr.Dataset(
       coords = dict(lon = longi,lat = lati),
       data_vars = dict(
                   mean_duration_pacifc = (['lat','lon'], mean_duration_arr),
                   std_duration_pacific = (['lat','lon'], std_duration_arr),  
           
                   mean_intensity_max_pacific = (['lat','lon'], mean_intensity_max_arr),  
                   std_intensity_max_pacific = (['lat','lon'], std_intensity_max_arr),  
                   
                   mean_intensity_mean_pacific = (['lat','lon'],mean_intensity_mean_arr), 
                   std_intensity_mean_pacific = (['lat','lon'],std_intensity_mean_arr),  
                  
                   mean_intensity_cumulative_pacifc = (['lat','lon'], mean_intensity_cumulative_arr),
                   std_intensity_cumulative_pacific = (['lat','lon'], std_intensity_cumulative_arr),  
           
                   mean_rateonset_pacific = (['lat','lon'], mean_rateonset_arr),  
                   std_rateonset_pacific = (['lat','lon'], std_rateonset_arr),  
                   
                   mean_ratedecline_pacific = (['lat','lon'],mean_ratedecline_arr), 
                   std_ratedecline_pacific = (['lat','lon'],std_ratedecline_arr),  
                  
                   mean_reaction_window_pacifc = (['lat','lon'], mean_reaction_window_arr),
                   std_reaction_window_pacific = (['lat','lon'], std_reaction_window_arr),  
           
                   mean_coping_window_pacific = (['lat','lon'], mean_coping_window_arr),  
                   std_coping_window_pacific = (['lat','lon'], std_coping_window_arr),  
                   
                   mean_recovery_window_pacific = (['lat','lon'],mean_recovery_window_arr), 
                   std_recovery_window_pacific = (['lat','lon'],std_recovery_window_arr),  
                  
                   mean_severity_pacific = (['lat','lon'], mean_severity_arr),
                   std_severity_pacific = (['lat','lon'], std_severity_arr),  
           
                   max_event_no_pacific = (['lat','lon'], max_event_no_arr)
                  
           
                       ))
             


# In[ ]:


xrds.to_netcdf('/home2/datawork/slal/noaaoisst_output_mean_std.nc')


# In[ ]:





# In[ ]:




