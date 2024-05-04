#!/usr/bin/env python
# coding: utf-8

# In[1]:


# country daily mean intensity plus minus 2sd, percent cover, max intensity on datarmor


# In[2]:


# standard imports

import numpy as np
from cartopy import config
import cartopy.crs as ccrs
import xarray as xr
import pandas as pd
import geopandas as gpd
import glob as glob
from shapely.geometry import Polygon

# In[ ]:


eez = gpd.read_file("/home2/datahome/slal/World_EEZ_v11_20191118/eez_v11.shp") # get the EEZs


# In[ ]:


# get eezs
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


# In[ ]:


df = pd.read_csv('/home2/datahome/slal/lonlist_clip_eez.csv') 
dfarrno = df.to_numpy()
#dfarrno.shape
nlon = dfarrno[:,1]


# In[ ]:


spatial = "/home/datawork-lead/datarmor-only/shilpa/spatial_extent_noaa_aoi_masking.nc"
#intensity = "/home2/datawork/slal/full_intensity_spatial_noaaoisst_aoi_lastdayincluded.nc"
severity_index = "/home/datawork-lead/datarmor-only/shilpa/aoi_mhw_severity_noaaoisst.nc"

ds_spatial = xr.open_dataset(spatial)
#ds_intensity = xr.open_dataset(intensity)
ds_severity_index = xr.open_dataset(severity_index)



# In[ ]:


def country_percent_area(country):    
    edgar_gdf_spatial.crs = country.crs
    country_data = gpd.clip(edgar_gdf_spatial, country)
    percent_area_cov_country = (np.nanmean(country_data.spatial_extent))*100
    return percent_area_cov_country





def country_percentage_severity_index_0_1(country): 
    lower_threshold = 0
    upper_threshold = 1
    edgar_gdf_severity.crs = country.crs
    country_data = gpd.clip(edgar_gdf_severity, country)
    country_data_filtered = country_data[(country_data['mhw_severity'] != 0)]
    filtered_data = country_data_filtered[(country_data_filtered['mhw_severity'] > lower_threshold)&(country_data_filtered['mhw_severity'] < upper_threshold)]
    percentage_severity_index_0_1 = ((len(filtered_data))/(len(country_data.mhw_severity)))*100
    return percentage_severity_index_0_1

def country_percentage_severity_index_1_2(country): 
    lower_threshold = 1
    upper_threshold = 2
    edgar_gdf_severity.crs = country.crs
    country_data = gpd.clip(edgar_gdf_severity, country)
    country_data_filtered = country_data[(country_data['mhw_severity'] != 0)]
    filtered_data = country_data_filtered[(country_data_filtered['mhw_severity'] >= lower_threshold)&(country_data_filtered['mhw_severity'] < upper_threshold)]
    percentage_severity_index_1_2 = ((len(filtered_data))/(len(country_data.mhw_severity)))*100
    return percentage_severity_index_1_2

def country_percentage_severity_index_2_3(country): 
    lower_threshold = 2
    upper_threshold = 3
    edgar_gdf_severity.crs = country.crs
    country_data = gpd.clip(edgar_gdf_severity, country)
    country_data_filtered = country_data[(country_data['mhw_severity'] != 0)]
    filtered_data = country_data_filtered[(country_data_filtered['mhw_severity'] >= lower_threshold)&(country_data_filtered['mhw_severity'] < upper_threshold)]
    percentage_severity_index_2_3 = ((len(filtered_data))/(len(country_data.mhw_severity)))*100
    return percentage_severity_index_2_3
    

def country_percentage_severity_index_3_4(country): 
    lower_threshold = 3
    upper_threshold = 4
    edgar_gdf_severity.crs = country.crs
    country_data = gpd.clip(edgar_gdf_severity, country)
    country_data_filtered = country_data[(country_data['mhw_severity'] != 0)]
    filtered_data = country_data_filtered[(country_data_filtered['mhw_severity'] >= lower_threshold)&(country_data_filtered['mhw_severity'] < upper_threshold)]
    percentage_severity_index_3_4 = ((len(filtered_data))/(len(country_data.mhw_severity)))*100
    return percentage_severity_index_3_4


def country_percentage_severity_index_4_more(country): 
    lower_threshold = 4
    #upper_threshold = 4
    edgar_gdf_severity.crs = country.crs
    country_data = gpd.clip(edgar_gdf_severity, country)
    country_data_filtered = country_data[(country_data['mhw_severity'] != 0)]
    filtered_data = country_data_filtered[(country_data_filtered['mhw_severity'] >= lower_threshold)]
    percentage_severity_index_4_more = ((len(filtered_data))/(len(country_data.mhw_severity)))*100
    return percentage_severity_index_4_more
    
    

def country_mean_severity(country):
    edgar_gdf_severity.crs = country.crs
    country_data = gpd.clip(edgar_gdf_severity, country)
    filtered_data = country_data[(country_data['mhw_severity'] > 0)]
    mean_filtered_data = (np.nanmean(filtered_data.mhw_severity)) 
    return mean_filtered_data

def country_std_severity(country):
    edgar_gdf_severity.crs = country.crs
    country_data = gpd.clip(edgar_gdf_severity, country)
    filtered_data = country_data[(country_data['mhw_severity'] > 0)]
    std_filtered_data = (np.nanstd(filtered_data.spatial_intensity))
    return std_filtered_data



# In[ ]:


datessfull = pd.date_range(start='9/1/1981', end='6/26/2023')
print(len(datessfull))
#princ(len(datessfull))


# In[ ]:




for i in range (9375,15274):
    print('i =',i)
    edgar_spatial = ds_spatial.spatial_extent_pacific[i,:,:].to_dataframe()
    edgar_spatial['newlons'] = nlon
    edgar_spatial = edgar_spatial.reset_index()
    edgar_gdf_spatial = gpd.GeoDataFrame(edgar_spatial, geometry=gpd.points_from_xy(edgar_spatial.newlons, edgar_spatial.lats))
    
    
    edgar_severity = ds_severity_index.mhw_severity[i,:,:].to_dataframe()
    edgar_severity['newlons'] = nlon
    edgar_severity = edgar_severity.reset_index()
    edgar_gdf_severity = gpd.GeoDataFrame(edgar_severity, geometry=gpd.points_from_xy(edgar_severity.newlons, edgar_severity.lats))
    
    fiji_area = country_percent_area(fiji)
    fiji_sevrity_0_1_percent = country_percentage_severity_index_0_1(fiji)
    fiji_sevrity_1_2_percent = country_percentage_severity_index_1_2(fiji)
    fiji_sevrity_2_3_percent = country_percentage_severity_index_2_3(fiji)
    fiji_sevrity_3_4_percent = country_percentage_severity_index_3_4(fiji)
    fiji_sevrity_4_more_percent = country_percentage_severity_index_4_more(fiji)
    fiji_sevrity_mean = country_mean_severity(fiji)
    fiji_sevrity_std = country_std_severity(fiji)



    ncd_area = country_percent_area(ncd)
    ncd_sevrity_0_1_percent = country_percentage_severity_index_0_1(ncd)
    ncd_sevrity_1_2_percent = country_percentage_severity_index_1_2(ncd)
    ncd_sevrity_2_3_percent = country_percentage_severity_index_2_3(ncd)
    ncd_sevrity_3_4_percent = country_percentage_severity_index_3_4(ncd)
    ncd_sevrity_4_more_percent = country_percentage_severity_index_4_more(ncd)
    ncd_sevrity_mean = country_mean_severity(ncd)
    ncd_sevrity_std = country_std_severity(ncd)

   
    vanuatu_area = country_percent_area(vanuatu)
    vanuatu_sevrity_0_1_percent = country_percentage_severity_index_0_1(vanuatu)
    vanuatu_sevrity_1_2_percent = country_percentage_severity_index_1_2(vanuatu)
    vanuatu_sevrity_2_3_percent = country_percentage_severity_index_2_3(vanuatu)
    vanuatu_sevrity_3_4_percent = country_percentage_severity_index_3_4(vanuatu)
    vanuatu_sevrity_4_more_percent = country_percentage_severity_index_4_more(vanuatu)
    vanuatu_sevrity_mean = country_mean_severity(vanuatu)
    vanuatu_sevrity_std = country_std_severity(vanuatu)


    wf_area = country_percent_area(wf)
    wf_sevrity_0_1_percent = country_percentage_severity_index_0_1(wf)
    wf_sevrity_1_2_percent = country_percentage_severity_index_1_2(wf)
    wf_sevrity_2_3_percent = country_percentage_severity_index_2_3(wf)
    wf_sevrity_3_4_percent = country_percentage_severity_index_3_4(wf)
    wf_sevrity_4_more_percent = country_percentage_severity_index_4_more(wf)
    wf_sevrity_mean = country_mean_severity(wf)
    wf_sevrity_std = country_std_severity(wf)

    
    solo_area = country_percent_area(solo)
    solo_sevrity_0_1_percent = country_percentage_severity_index_0_1(solo)
    solo_sevrity_1_2_percent = country_percentage_severity_index_1_2(solo)
    solo_sevrity_2_3_percent = country_percentage_severity_index_2_3(solo)
    solo_sevrity_3_4_percent = country_percentage_severity_index_3_4(solo)
    solo_sevrity_4_more_percent = country_percentage_severity_index_4_more(solo)
    solo_sevrity_mean = country_mean_severity(solo)
    solo_sevrity_std = country_std_severity(solo)

    
    tonga_area = country_percent_area(tonga)
    tonga_sevrity_0_1_percent = country_percentage_severity_index_0_1(tonga)
    tonga_sevrity_1_2_percent = country_percentage_severity_index_1_2(tonga)
    tonga_sevrity_2_3_percent = country_percentage_severity_index_2_3(tonga)
    tonga_sevrity_3_4_percent = country_percentage_severity_index_3_4(tonga)
    tonga_sevrity_4_more_percent = country_percentage_severity_index_4_more(tonga)
    tonga_sevrity_mean = country_mean_severity(tonga)
    tonga_sevrity_std = country_std_severity(tonga)

    
    samoa_area = country_percent_area(samoa)
    samoa_sevrity_0_1_percent = country_percentage_severity_index_0_1(samoa)
    samoa_sevrity_1_2_percent = country_percentage_severity_index_1_2(samoa)
    samoa_sevrity_2_3_percent = country_percentage_severity_index_2_3(samoa)
    samoa_sevrity_3_4_percent = country_percentage_severity_index_3_4(samoa)
    samoa_sevrity_4_more_percent = country_percentage_severity_index_4_more(samoa)
    samoa_sevrity_mean = country_mean_severity(samoa)
    samoa_sevrity_std = country_std_severity(samoa)

    
    tuvalu_area = country_percent_area(tuvalu)
    tuvalu_sevrity_0_1_percent = country_percentage_severity_index_0_1(tuvalu)
    tuvalu_sevrity_1_2_percent = country_percentage_severity_index_1_2(tuvalu)
    tuvalu_sevrity_2_3_percent = country_percentage_severity_index_2_3(tuvalu)
    tuvalu_sevrity_3_4_percent = country_percentage_severity_index_3_4(tuvalu)
    tuvalu_sevrity_4_more_percent = country_percentage_severity_index_4_more(tuvalu)
    tuvalu_sevrity_mean = country_mean_severity(tuvalu)
    tuvalu_sevrity_std = country_std_severity(tuvalu)


    niue_area = country_percent_area(niue)
    niue_sevrity_0_1_percent = country_percentage_severity_index_0_1(niue)
    niue_sevrity_1_2_percent = country_percentage_severity_index_1_2(niue)
    niue_sevrity_2_3_percent = country_percentage_severity_index_2_3(niue)
    niue_sevrity_3_4_percent = country_percentage_severity_index_3_4(niue)
    niue_sevrity_4_more_percent = country_percentage_severity_index_4_more(niue)
    niue_sevrity_mean = country_mean_severity(niue)
    niue_sevrity_std = country_std_severity(niue)

    cooks_area = country_percent_area(cooks)
    cooks_sevrity_0_1_percent = country_percentage_severity_index_0_1(cooks)
    cooks_sevrity_1_2_percent = country_percentage_severity_index_1_2(cooks)
    cooks_sevrity_2_3_percent = country_percentage_severity_index_2_3(cooks)
    cooks_sevrity_3_4_percent = country_percentage_severity_index_3_4(cooks)
    cooks_sevrity_4_more_percent = country_percentage_severity_index_4_more(cooks)
    cooks_sevrity_mean = country_mean_severity(cooks)
    cooks_sevrity_std = country_std_severity(cooks)

    tokelau_area = country_percent_area(tokelau)
    tokelau_sevrity_0_1_percent = country_percentage_severity_index_0_1(tokelau)
    tokelau_sevrity_1_2_percent = country_percentage_severity_index_1_2(tokelau)
    tokelau_sevrity_2_3_percent = country_percentage_severity_index_2_3(tokelau)
    tokelau_sevrity_3_4_percent = country_percentage_severity_index_3_4(tokelau)
    tokelau_sevrity_4_more_percent = country_percentage_severity_index_4_more(tokelau)
    tokelau_sevrity_mean = country_mean_severity(tokelau)
    tokelau_sevrity_std = country_std_severity(tokelau)

    
    amsam_area = country_percent_area(amsam)
    amsam_sevrity_0_1_percent = country_percentage_severity_index_0_1(amsam)
    amsam_sevrity_1_2_percent = country_percentage_severity_index_1_2(amsam)
    amsam_sevrity_2_3_percent = country_percentage_severity_index_2_3(amsam)
    amsam_sevrity_3_4_percent = country_percentage_severity_index_3_4(amsam)
    amsam_sevrity_4_more_percent = country_percentage_severity_index_4_more(amsam)
    amsam_sevrity_mean = country_mean_severity(amsam)
    amsam_sevrity_std = country_std_severity(amsam)
 


    d = {'time_index': i, 'time' : datessfull[i],
         
    'Fiji_%spatial_extent':fiji_area,
    'Fiji_sevrity_0_1_percent' : fiji_sevrity_0_1_percent,
    'Fiji_sevrity_1_2_percent' : fiji_sevrity_1_2_percent,
    'Fiji_sevrity_2_3_percent' : fiji_sevrity_2_3_percent,
    'Fiji_sevrity_3_4_percent' : fiji_sevrity_3_4_percent,
    'Fiji_sevrity_4_more_percent' : fiji_sevrity_4_more_percent,
    'Fiji_sevrity_mean' : fiji_sevrity_mean,
    'Fiji_sevrity_std'  : fiji_sevrity_std,                   
         
    'New_Caledonia_%spatial_extent':ncd_area,
    'New_Caledonia_sevrity_0_1_percent' : ncd_sevrity_0_1_percent,
    'New_Caledonia_sevrity_1_2_percent' : ncd_sevrity_1_2_percent,
    'New_Caledonia_sevrity_2_3_percent' : ncd_sevrity_2_3_percent,
    'New_Caledonia_sevrity_3_4_percent' : ncd_sevrity_3_4_percent,
    'New_Caledonia_sevrity_4_more_percent' : ncd_sevrity_4_more_percent,
    'New_Caledonia_sevrity_mean' : ncd_sevrity_mean,
    'New_Caledonia_sevrity_std'  : ncd_sevrity_std,      
     
    'Vanuatu_spatial_%extent':vanuatu_area,
    'Vanuatu_sevrity_0_1_percent' : vanuatu_sevrity_0_1_percent,
    'Vanuatu_sevrity_1_2_percent' : vanuatu_sevrity_1_2_percent,
    'Vanuatu_sevrity_2_3_percent' : vanuatu_sevrity_2_3_percent,
    'Vanuatu_sevrity_3_4_percent' : vanuatu_sevrity_3_4_percent,
    'Vanuatu_sevrity_4_more_percent' : vanuatu_sevrity_4_more_percent,
    'Vanuatu_sevrity_mean' : vanuatu_sevrity_mean,
    'Vanuatu_sevrity_std'  : vanuatu_sevrity_std, 
         
         
    'Wallis_and_Futuna_%spatial_extent':wf_area,
    'Wallis_and_Futuna_sevrity_0_1_percent' : wf_sevrity_0_1_percent,
    'Wallis_and_Futuna_sevrity_1_2_percent' : wf_sevrity_1_2_percent,
    'Wallis_and_Futuna_sevrity_2_3_percent' : wf_sevrity_2_3_percent,
    'Wallis_and_Futuna_sevrity_3_4_percent' : wf_sevrity_3_4_percent,
    'Wallis_and_Futuna_sevrity_4_more_percent' : wf_sevrity_4_more_percent,
    'Wallis_and_Futuna_sevrity_mean' : wf_sevrity_mean,
    'Wallis_and_Futuna_sevrity_std'  : wf_sevrity_std,  

         
         
    'Solomon_Islands_%spatial_extent':solo_area,
    'Solomon_Islands_sevrity_0_1_percent' : solo_sevrity_0_1_percent,
    'Solomon_Islands_sevrity_1_2_percent' : solo_sevrity_1_2_percent,
    'Solomon_Islands_sevrity_2_3_percent' : solo_sevrity_2_3_percent,
    'Solomon_Islands_sevrity_3_4_percent' : solo_sevrity_3_4_percent,
    'Solomon_Islands_sevrity_4_more_percent' : solo_sevrity_4_more_percent,
    'Solomon_Islands_sevrity_mean' : solo_sevrity_mean,
    'Solomon_Islands_sevrity_std'  : solo_sevrity_std,

         
    'Tonga_%spatial_extent':tonga_area,
    'Tonga_sevrity_0_1_percent' : tonga_sevrity_0_1_percent,
    'Tonga_sevrity_1_2_percent' : tonga_sevrity_1_2_percent,
    'Tonga_sevrity_2_3_percent' : tonga_sevrity_2_3_percent,
    'Tonga_sevrity_3_4_percent' : tonga_sevrity_3_4_percent,
    'Tonga_sevrity_4_more_percent' : tonga_sevrity_4_more_percent,
    'Tonga_sevrity_mean' : tonga_sevrity_mean,
    'Tonga_sevrity_std'  : tonga_sevrity_std,

    
    'Samoa_%spatial_extent':samoa_area,
    'Samoa_sevrity_0_1_percent' : samoa_sevrity_0_1_percent,
    'Samoa_sevrity_1_2_percent' : samoa_sevrity_1_2_percent,
    'Samoa_sevrity_2_3_percent' : samoa_sevrity_2_3_percent,
    'Samoa_sevrity_3_4_percent' : samoa_sevrity_3_4_percent,
    'Samoa_sevrity_4_more_percent' : samoa_sevrity_4_more_percent,
    'Samoa_sevrity_mean' : samoa_sevrity_mean,
    'Samoa_sevrity_std'  : samoa_sevrity_std,

         
    'Tuvalu_%spatial_extent':tuvalu_area,
    'Tuvalu_sevrity_0_1_percent' : tuvalu_sevrity_0_1_percent,
    'Tuvalu_sevrity_1_2_percent' : tuvalu_sevrity_1_2_percent,
    'Tuvalu_sevrity_2_3_percent' : tuvalu_sevrity_2_3_percent,
    'Tuvalu_sevrity_3_4_percent' : tuvalu_sevrity_3_4_percent,
    'Tuvalu_sevrity_4_more_percent' : tuvalu_sevrity_4_more_percent,
    'Tuvalu_sevrity_mean' : tuvalu_sevrity_mean,
    'Tuvalu_sevrity_std'  : tuvalu_sevrity_std,


         
    'Niue_%spatial_extent':niue_area,
    'Niue_sevrity_0_1_percent' : niue_sevrity_0_1_percent,
    'Niue_sevrity_1_2_percent' : niue_sevrity_1_2_percent,
    'Niue_sevrity_2_3_percent' : niue_sevrity_2_3_percent,
    'Niue_sevrity_3_4_percent' : niue_sevrity_3_4_percent,
    'Niue_sevrity_4_more_percent' : niue_sevrity_4_more_percent,
    'Niue_sevrity_mean' : niue_sevrity_mean,
    'Niue_sevrity_std'  : niue_sevrity_std,

         
    'Cook_Islands_%spatial_extent':cooks_area,
    'Cook_Islands_sevrity_0_1_percent' : cooks_sevrity_0_1_percent,
    'Cook_Islands_sevrity_1_2_percent' : cooks_sevrity_1_2_percent,
    'Cook_Islands_sevrity_2_3_percent' : cooks_sevrity_2_3_percent,
    'Cook_Islands_sevrity_3_4_percent' : cooks_sevrity_3_4_percent,
    'Cook_Islands_sevrity_4_more_percent' : cooks_sevrity_4_more_percent,
    'Cook_Islands_sevrity_mean' : cooks_sevrity_mean,
    'Cook_Islands_sevrity_std'  : cooks_sevrity_std,


  

    'Tokelau_%spatial_extent':tokelau_area,
    'Tokelau_sevrity_0_1_percent' : tokelau_sevrity_0_1_percent,
    'Tokelau_sevrity_1_2_percent' : tokelau_sevrity_1_2_percent,
    'Tokelau_sevrity_2_3_percent' : tokelau_sevrity_2_3_percent,
    'Tokelau_sevrity_3_4_percent' : tokelau_sevrity_3_4_percent,
    'Tokelau_sevrity_4_more_percent' : tokelau_sevrity_4_more_percent,
    'Tokelau_sevrity_mean' : tokelau_sevrity_mean,
    'Tokelau_sevrity_std'  : tokelau_sevrity_std,


   
  
         
    'American_Samoa_%spatial_extent':amsam_area,
    'American_Samoa_sevrity_0_1_percent' : amsam_sevrity_0_1_percent,
    'American_Samoa_sevrity_1_2_percent' : amsam_sevrity_1_2_percent,
    'American_Samoa_sevrity_2_3_percent' : amsam_sevrity_2_3_percent,
    'American_Samoa_sevrity_3_4_percent' : amsam_sevrity_3_4_percent,
    'American_Samoa_sevrity_4_more_percent' : amsam_sevrity_4_more_percent,
    'American_Samoa_sevrity_mean' : amsam_sevrity_mean,
    'American_Samoa_sevrity_std'  : amsam_sevrity_std


    
      }

    dfnew = pd.DataFrame(data=d,index=[0])
    dfnew.to_csv(f"/home2/datawork/slal/country_daily_noaaoisst_1981_2023/percent_mhwstate_severity_countries_noaaoisst_{i}.csv")



# In[ ]:


#fnames = glob.glob('/home2/datawork/slal/country_daily_noaaoisst/percent_mhwstate_max_mean_*')


# In[ ]:


#df_concat2 = pd.concat([pd.read_csv(f) for f in fnames], ignore_index = True)
#df_concat2


# In[ ]:


#df_concat2.to_csv('/home2/datawork/slal/country_daily_noaaoisst/fullset_percent_mhwstate_max_mean.csv')


# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:




