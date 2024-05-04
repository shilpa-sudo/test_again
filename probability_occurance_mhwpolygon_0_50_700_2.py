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

var = xr.open_dataset('/home/shilpa/glory_mat_analysis/means_std_median_polyareas_noaaoisst/noaaoisst_mean_std_mhw_prop/probability_of_occurance_polyareas_specific_sizes_range_noaa_93_22_new_0_25_700.nc')

#prob of occurance only_mhw_days

var_0_25 = np.array(var.prob_occu_of_sp_size_poly_only_mhw_days_0_25)
var_25_700 = np.array(var.prob_occu_of_sp_size_poly_only_mhw_days_25_700)
#var_10_25 = np.array(var.prob_occu_of_sp_size_poly_only_mhw_days_10_25)
#var_25_50 = np.array(var.prob_occu_of_sp_size_poly_only_mhw_days_25_50)
#var_50_100 = np.array(var.prob_occu_of_sp_size_poly_only_mhw_days_50_100)
#var_100_200 = np.array(var.prob_occu_of_sp_size_poly_only_mhw_days_100_200)
#var_200_500 = np.array(var.prob_occu_of_sp_size_poly_only_mhw_days_200_500)
#var_500_700 = np.array(var.prob_occu_of_sp_size_poly_only_mhw_days_500_700)



print (np.nanmean(var_0_25),np.nanmin(var_0_25),np.nanmax(var_0_25))
print (np.nanmean(var_25_700),np.nanmin(var_25_700),np.nanmax(var_25_700))
#print (np.nanmean(var_10_25),np.nanmin(var_10_25),np.nanmax(var_10_25))
#print (np.nanmean(var_25_50),np.nanmin(var_25_50),np.nanmax(var_25_50))
#print (np.nanmean(var_50_100),np.nanmin(var_50_100),np.nanmax(var_50_100))
#print (np.nanmean(var_100_200),np.nanmin(var_100_200),np.nanmax(var_100_200))
#print (np.nanmean(var_200_500),np.nanmin(var_200_500),np.nanmax(var_200_500))
#print (np.nanmean(var_500_700),np.nanmin(var_500_700),np.nanmax(var_500_700))


#min_mean=min([np.nanmin(var_0_1),np.nanmin(var_1_10),np.nanmin(var_10_25),np.nanmin(var_25_50),np.nanmin(var_50_100),np.nanmin(var_100_200),np.nanmin(var_200_500),np.nanmin(var_500_700)])
#max_mean=max([np.nanmax(var_0_1),np.nanmax(var_1_10),np.nanmax(var_10_25),np.nanmax(var_25_50),np.nanmax(var_50_100),np.nanmax(var_100_200),np.nanmax(var_200_500),np.nanmax(var_500_700)])




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
lat = var.lat
lon = var.lon

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


variablename = 'Probability of occurance (Any MHW day)'
units_ = 'sqr degs'
units = 'Probability'

nrows = 2
ncols = 2

fig, axs = plt.subplots(nrows=nrows,ncols=ncols,subplot_kw={'projection': proj_crs},figsize=(20,20))

uneven_levels = [0,0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1,0.12,0.14,0.16,0.18,0.2,0.22,0.24,0.26,0.28,0.3,0.32,0.34,0.36,0.38,0.4,0.42,0.44,0.46,0.48,0.5,0.55,0.6,0.65,0.7,0.75,1]
easy_plot_unevencolorbars(lon=xx,lat=yy, data=var_0_25,
                  uneven_levels=uneven_levels,name=f'{variablename} 0-25 {units_}',units=f'{units}',i=0,j=0)

easy_plot_unevencolorbars(lon=xx,lat=yy, data=var_25_700,
                  uneven_levels=uneven_levels,name=f'{variablename} 25-700 {units_}',units=f'{units}',i=0,j=1)

easy_plot_unevencolorbars(lon=xx,lat=yy, data=var_0_25,
                  uneven_levels=uneven_levels,name=f'{variablename} 0-25 {units_}',units=f'{units}',i=1,j=0)

easy_plot_unevencolorbars(lon=xx,lat=yy, data=var_25_700,
                  uneven_levels=uneven_levels,name=f'{variablename} 25-700 {units_}',units=f'{units}',i=1,j=1)




axs[0, 0].text(-0.2, 1.05, 'a', transform=axs[0, 0].transAxes, fontsize=20, fontweight='bold')
axs[0, 1].text(-0.2, 1.05, 'b', transform=axs[0, 1].transAxes, fontsize=20, fontweight='bold')
axs[1, 0].text(-0.2, 1.05, 'C', transform=axs[1, 0].transAxes, fontsize=20, fontweight='bold')
axs[1, 1].text(-0.2, 1.05, 'd', transform=axs[1, 1].transAxes, fontsize=20, fontweight='bold')
#axs[2, 0].text(-0.2, 1.05, 'e', transform=axs[2, 0].transAxes, fontsize=20, fontweight='bold')
#axs[2, 1].text(-0.2, 1.05, 'f', transform=axs[2, 1].transAxes, fontsize=20, fontweight='bold')
#axs[3, 0].text(-0.2, 1.05, 'g', transform=axs[3, 0].transAxes, fontsize=20, fontweight='bold')
#axs[3, 1].text(-0.2, 1.05, 'h', transform=axs[3, 1].transAxes, fontsize=20, fontweight='bold')


plt.tight_layout()
plt.savefig('prob_occurance_mhw_days_polygon_size_range_0_25_700.png',bbox_inches='tight')
