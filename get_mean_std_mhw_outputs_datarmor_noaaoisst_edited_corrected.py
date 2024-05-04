#!/usr/bin/env python
# coding: utf-8

# In[2]:


import xarray as xr
import numpy as np


# In[ ]:





# In[ ]:


ds = xr.open_dataset('/home/datawork-lead/datarmor-only/shilpa/noaaoisst_mhw_output_test_ratedecline.nc')


# In[ ]:


#duration = ds.duration_pacific
#intensity_max = ds.intensity_max_pacific
#intensity_mean = ds.intensity_mean_pacific
#intensity_cumulative = ds.intensity_cumulative_pacific
#rateonset = ds.rateonset_pacific
#ratedecline = ds.ratedecline_arr
#reaction_window = ds.reaction_window_days_pacific
#coping_window = ds.coping_window_days_pacific
#recovery_window = ds.recovery_window_days_pacific
#severity = ds.severity_index_pacific
#event_no = ds.mhw_event_no_pacific





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


max_event_no_arr = np.zeros((144,588)) 


# In[ ]:


for i in range(144):
    for j in range(588):
        
        
        duration = ds.duration_pacific[:,i,j][ds.duration_pacific[:,i,j] > 0]
        #if duration != []:
	#print (len(duration))
        #if (len(duration)) != 0:

        mean_duration = np.mean(duration)
        mean_duration_arr[i,j] = mean_duration
        std_duration = np.std(duration)
        std_duration_arr[i,j] = std_duration

        intensity_max = ds.intensity_max_pacific[:,i,j][ds.intensity_max_pacific[:,i,j] > 0]
        mean_intensity_max = np.mean(intensity_max)
        mean_intensity_max_arr[i,j] = mean_intensity_max 
        std_intensity_max = np.std(intensity_max)
        std_intensity_max_arr[i,j] = std_intensity_max
        
        
        intensity_mean = ds.intensity_mean_pacific[:,i,j][ds.intensity_mean_pacific[:,i,j] > 0]
        mean_intensity_mean = np.mean(intensity_mean)
        mean_intensity_mean_arr[i,j] = mean_intensity_mean
        std_intensity_mean = np.std(intensity_mean)
        std_intensity_mean_arr[i,j] =  std_intensity_mean 
        
        intensity_cumulative = ds.intensity_cumulative_pacific[:,i,j][ds.intensity_cumulative_pacific[:,i,j] > 0]
        mean_intensity_cumulative = np.mean(intensity_cumulative)
        mean_intensity_cumulative_arr[i,j] = mean_intensity_cumulative
        std_intensity_mean = np.std(intensity_cumulative)
        std_intensity_mean_arr[i,j] =  std_intensity_cumulative
       
        rateonset = ds.rateonset_pacific[:,i,j][ds.rateonset_pacific[:,i,j] > 0]
        mean_rateonset = np.mean(rateonset)
        mean_rateonset_arr[i,j] = mean_rateonset
        std_rateonset = np.std(rateonset)
        std_rateonset_arr[i,j] =  std_rateonset

        ratedecline = ds.ratedecline_arr[:,i,j][ds.ratedecline_arr[:,i,j] > 0]
        mean_ratedecline = np.mean(ratedecline)
        mean_ratedecline_arr[i,j] = mean_ratedecline
        std_ratedecline = np.std(ratedecline)
        std_ratedecline_arr[i,j] =  std_ratedecline
        
        reaction_window = ds.reaction_window_days_pacific[:,i,j][ds.reaction_window_days_pacific[:,i,j] > 0]
        mean_reaction_window = np.mean(reaction_window)
        mean_reaction_window_arr[i,j] = mean_reaction_window
        std_reaction_window = np.std(reaction_window)
        std_reaction_window_arr[i,j] = std_reaction_window
        
        coping_window = ds.coping_window_days_pacific[:,i,j][ds.coping_window_days_pacific[:,i,j] > 0]
        mean_coping_window = np.mean(coping_window)
        mean_coping_window_arr[i,j] = mean_coping_window
        std_coping_window = np.std(coping_window)
        std_coping_window_arr[i,j] = std_coping_window
        
        recovery_window = ds.recovery_window_days_pacific[:,i,j][ds.recovery_window_days_pacific[:,i,j] > 0]
        mean_recovery_window = np.mean(recovery_window)
        mean_recovery_window_arr[i,j] = mean_recovery_window
        std_recovery_window = np.std(recovery_window)
        std_recovery_window_arr[i,j] = std_recovery_window 
        
        severity = ds.severity_index_pacific[:,i,j][ds.severity_index_pacific[:,i,j] > 0]
        mean_severity = np.mean(severity)
        mean_severity_arr[i,j] = mean_severity
        std_severity = np.std(severity)
        std_severity_arr[i,j] = std_severity
        
        mhw_event_no = ds.mhw_event_no_pacific[:,i,j]
        max_event_no = np.max(mhw_event_no)
        max_event_no_arr[i,j] = max_event_no
       


# In[ ]:


longi = np.array(ds.lon.values)
lati = np.array(ds.lat.values)


# In[ ]:


xrds = xr.Dataset(
       coords = dict(lon = longi,lat = lati),
       data_vars = dict(
                   mean_duration_pacifc = (['lat','lon'], mean_duration_arr),
                   std_duration_pacific = (['lat','lon'], std_duration_arr) ,
           
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


xrds.to_netcdf('/home/datawork-lead/datarmor-only/shilpa/noaaoisst_output_mean_std_all_parameters.nc')


# In[ ]:





# In[ ]:




#!/usr/bin/env python
# coding: utf-8

# In[2]:


#import xarray as xr
#import numpy as np


# In[ ]:





# In[ ]:

