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

years_list = list(range(1993, 2020))
print(years_list)

filelist = []
for year in years_list:
	file =  f'/media/shilpa/Expansion/regridded_glorys_files_for_coastal/regridded_glorys_{year}_0m.nc'
	filelist.append(file)

filelist2 = ['/media/shilpa/Expansion/glorys_20200101_20210630_0.49m.nc', '/media/shilpa/Expansion/glorys_20210701_20231024_0.49m.nc']

df_lon = pd.read_csv('/media/shilpa/Expansion/regridded_glorys_files_for_coastal/glorys_lon.csv')
df_lat = pd.read_csv('/media/shilpa/Expansion/regridded_glorys_files_for_coastal/glorys_lat.csv')

df_lon['adjusted_longitude'] = np.where(df_lon['longitude'] < 0, df_lon['longitude'] + 360, df_lon['longitude'])

longi = df_lon.adjusted_longitude
lati = df_lat.latitude

latarr = (lati).tolist()
print(len(latarr))
lonarr = (longi).tolist()
print(len(lonarr))

get_coords = pd.read_csv('/media/shilpa/Expansion/regridded_glorys_files_for_coastal/coastal_points/tokelau_coastalpoints.csv')

for index, row in get_coords.iterrows():  # for each lat and lon value in newcal coast
    temparray  = []
    timearray = []
    
    target_latitude = row['latitude']
    target_longitude = row['longitude']
    
    # Adjust lon if it's less than 0
    if target_longitude < 0:
        target_longitude += 360

    print(target_latitude,target_longitude)

    country_name = row['Country']

    target_latitude_ind = min(range(len(latarr)), key=lambda i:abs((latarr[i]) - (target_latitude)))
    target_longitude_ind = min(range(len(lonarr)), key=lambda i:abs((lonarr[i]) - (target_longitude)))
	
    def _preprocess(xx, target_latitude_ind, target_longitude_ind):
        return xx.isel(deptht=[0], y=[target_latitude_ind], x=[target_longitude_ind])
    
    ds = xr.open_mfdataset(filelist, combine='nested', concat_dim="time_counter", preprocess=lambda ds: _preprocess(ds, target_latitude_ind, target_longitude_ind))
    da = ds['votemper']
    numpyarr = da.values  # Alternative to da.as_numpy()
    temp = numpyarr[:, 0, 0, 0]
    temparray.append(temp)
	#print(numpyarr.shape)

    time = ds.time_counter
    strtime = time.data[:].astype(str)
    pdtime = pd.to_datetime(strtime)
    t = np.array([date.toordinal(x) for x in pdtime])
    timearray.append(t)

    def _preprocess2(xx, target_latitude_ind, target_longitude_ind):
        return xx.isel(depth=[0], latitude =[target_latitude_ind], longitude =[target_longitude_ind])
    ds2 = xr.open_mfdataset(filelist2, combine='nested', concat_dim="time", preprocess=lambda ds2: _preprocess2(ds2, target_latitude_ind, target_longitude_ind))
    da2 = ds2['thetao']
    da2_ = da2.ffill(dim='time')  #forward fill to give same numerical value to nans
    numpyarr2 = da2_.values  # Alternative to da.as_numpy()
    temp2 = numpyarr2[:, 0, 0, 0]
    temparray.append(temp2)
    time2 = ds2.time
    strtime2 = time2.data[:].astype(str)
    pdtime2 = pd.to_datetime(strtime2)
    t2 = np.array([date.toordinal(x) for x in pdtime2])
    timearray.append(t2)

    concatenated_array_temp = np.concatenate((temparray[0] , temparray[1]), axis=0)
    concatenated_array_time = np.concatenate((timearray[0] , timearray[1]), axis=0)

    t_ = concatenated_array_time
    temp_ = concatenated_array_temp
    print(t_.shape, temp_.shape)

    mhws, clim = mhw.detect(t_, temp_)

    events = mhws['n_events']
	#print(events)

    if events == 0:
        print('no events at this location')

    else:
        date_start = mhws['date_start']
        date_end = mhws['date_end']
        index_start = mhws['index_start']
        index_end = mhws['index_end']
        index_peak = mhws['index_peak']
        date_peak = mhws['date_peak']
        duration = mhws['duration']
        intensity_max = mhws['intensity_max']
        #intensity_max_max = mhws['intensity_max_max']
        intensity_mean = mhws['intensity_mean']
        intensity_var = mhws['intensity_var']
        intensity_cumulative = mhws['intensity_cumulative']
        category = mhws['category']
        rateonset = mhws['rate_onset']
        ratedecline = mhws['rate_decline']
        lat = [target_latitude] * events
        lon = [target_longitude] * events
        country = [country_name] * events
        reaction_window_days =  [x - y for x, y in zip((mhws['index_peak']),(mhws['index_start']))]
        coping_window_days =  [x - y for x, y in zip((mhws['index_end']),(mhws['index_peak']))]
        recovery_window_days =   [(x - y for x, y in zip((mhws['index_start'][1:]),(mhws['index_end'][:-1])))]

        d = {'country': country,'latitude':lat, 'longitude':lon,'index_start': index_start, 'index_end': index_end,'index_peak': index_peak,
                 'date_start': date_start, 'date_end': date_end, 'date_peak': date_peak,
                 'duration' : duration, 'intensity_max': intensity_max, 'intensity_mean' : intensity_mean,
                 'intensity_cumulative':intensity_cumulative, 'category':category, 'rateonset': rateonset,
                 'ratedecline':ratedecline, 'reaction_window_days':reaction_window_days,'coping_window_days':coping_window_days}
                 #'recovery_window_days':recovery_window_days}

        df = pd.DataFrame(data=d)
        df.to_csv('/media/shilpa/Expansion/coastal_tokelau/%s_all_mhws_at_%s_%s.csv'%(country_name,target_latitude,target_longitude))

    		# spatial extent script comes here
        index_start = mhws['index_start']
        index_end = mhws['index_end']
        z = np.zeros(len(t_))
        for i in range(events):
            z[index_start[i]:index_end[i]]=1

        spatial_extent = z

        lat = [target_latitude] * (len(t_))
        lon = [target_longitude] * (len(t_))
        country = [country_name] * (len(t_))

        d3 = {'country':country,'latitude':lat, 'longitude':lon,
            'time' : t_, 'spatial_extent' : spatial_extent }
        df3 = pd.DataFrame(data=d3)
        df3.to_csv('/media/shilpa/Expansion/coastal_tokelau/%s_spatial_extent_mhws_at_%s_%s.csv'%(country_name,target_latitude,target_longitude))



        mhwBlock = mhw.blockAverage(t_, mhws)

        total_days = mhwBlock['total_days']
        total_icum = mhwBlock['total_icum']

        mean_per_year,trend, dtrend = mhw.meanTrend(mhwBlock,alpha=0.05)
    		#mean_per_year
        count_mean_per_year = mean_per_year['count']
        duration_mean_per_year = mean_per_year['duration']
        intensity_mean_mean_per_year = mean_per_year['intensity_max']
        intensity_max_max_mean_per_year = mean_per_year['intensity_max_max']
        intensity_mean_mean_per_year = mean_per_year['intensity_mean']
        intensity_var_mean_per_year = mean_per_year['intensity_var']
        intensity_cumulative_mean_per_year = mean_per_year['intensity_cumulative']
        rate_onset_mean_per_year = mean_per_year['rate_onset']
        rate_decline_mean_per_year = mean_per_year['rate_decline']
        total_days_mean_per_year = mean_per_year['total_days']
        total_icum_mean_per_year = mean_per_year['total_icum']

    		#trend_per_decade
        count_trend_per_decade = (10*trend['count'])
        duration_trend_per_decade = (10*trend['duration'])
        intensity_max_trend_per_decade = (10*trend['intensity_max'])
        intensity_max_max_trend_per_decade = (10*trend['intensity_max_max'])
        intensity_mean_trend_per_decade = (10*trend['intensity_mean'])
        intensity_var_trend_per_decade = (10*trend['intensity_var'])
        intensity_cumulative_trend_per_decade = (10*trend['intensity_cumulative'])
        rate_onset_trend_per_decade = (10*trend['rate_onset'])
        rate_decline_trend_per_decade = (10*trend['rate_decline'])
        total_days_trend_per_decade = (10*trend['total_days'])
        total_icum_trend_per_decade = (10*trend['total_icum'])



		#significance of trend
        if (np.abs(trend['count']) > dtrend['count']) == False:
            sig_trend_per_decade_count =  0
        elif (np.abs(trend['count']) > dtrend['count'])== True:
            sig_trend_per_decade_count =  1
    			#print(sig_trend_per_decade_count)

        if (np.abs(trend['duration']) > dtrend['duration']) == False:
            sig_trend_per_decade_duration =  0
        elif (np.abs(trend['duration']) > dtrend['duration'])== True:
            sig_trend_per_decade_duration =  1
    			#print(sig_trend_per_decade_duration)
        if (np.abs(trend['intensity_max']) > dtrend['intensity_max']) == False:
            sig_trend_per_decade_intensity_max =  0
        elif (np.abs(trend['intensity_max']) > dtrend['intensity_max'])== True:
            sig_trend_per_decade_intensity_max =  1

        if (np.abs(trend['intensity_max_max']) > dtrend['intensity_max_max']) == False:
            sig_trend_per_decade_intensity_max_max =  0
        elif (np.abs(trend['intensity_max_max']) > dtrend['intensity_max_max'])== True:
            sig_trend_per_decade_intensity_max_max =  1

        if (np.abs(trend['intensity_mean']) > dtrend['intensity_mean']) == False:
            sig_trend_per_decade_intensity_mean =  0
        elif (np.abs(trend['intensity_mean']) > dtrend['intensity_mean'])== True:
            sig_trend_per_decade_intensity_mean =  1

        if (np.abs(trend['intensity_var']) > dtrend['intensity_var']) == False:
            sig_trend_per_decade_intensity_var =  0
        elif (np.abs(trend['intensity_var']) > dtrend['intensity_var'])== True:
            sig_trend_per_decade_intensity_var =  1
        if (np.abs(trend['intensity_cumulative']) > dtrend['intensity_cumulative']) == False:
            sig_trend_per_decade_intensity_cumulative =  0
        elif (np.abs(trend['intensity_cumulative']) > dtrend['intensity_cumulative'])== True:
            sig_trend_per_decade_intensity_cumulative =  1

        if (np.abs(trend['rate_onset']) > dtrend['rate_onset']) == False:
            sig_trend_per_decade_rate_onset =  0
        elif (np.abs(trend['rate_onset']) > dtrend['rate_onset'])== True:
            sig_trend_per_decade_rate_onset =  1

        if (np.abs(trend['rate_decline']) > dtrend['rate_decline']) == False:
            sig_trend_per_decade_rate_decline =  0
        elif (np.abs(trend['rate_decline']) > dtrend['rate_decline'])== True:
            sig_trend_per_decade_rate_decline =  1


        if (np.abs(trend['total_days']) > dtrend['total_days']) == False:
            sig_trend_per_decade_total_days =  0
        elif (np.abs(trend['total_days']) > dtrend['total_days'])== True:
            sig_trend_per_decade_total_days =  1

        if (np.abs(trend['total_icum']) > dtrend['total_icum']) == False:
            sig_trend_per_decade_total_icum =  0
        elif (np.abs(trend['total_icum']) > dtrend['total_icum'])== True:
            sig_trend_per_decade_total_icum =  1

        mean_duration = statistics.mean(duration)
        std_duration = statistics.stdev(duration)

        mean_intensity_mean = statistics.mean(intensity_mean)
        std_intensity_mean = statistics.stdev(intensity_mean)

        mean_intensity_var = statistics.mean(intensity_var)
        std_intensity_var = statistics.stdev(intensity_var)

        mean_intensity_max = statistics.mean(intensity_max)
        std_intensity_max = statistics.stdev(intensity_max)

    		#mean_intensity_max_max = statistics.mean(intensity_max_max)
    		#std_intensity_max_max = statistics.stdev(intensity_max_max)

        mean_intensity_cumulative = statistics.mean(intensity_cumulative)
        std_intensity_cumulative = statistics.stdev(intensity_cumulative)

        mean_rate_onset = statistics.mean(rateonset)
        std_rate_onset = statistics.stdev(rateonset)

        mean_rate_decline = statistics.mean(ratedecline)
        std_rate_decline = statistics.stdev(ratedecline)

        mean_total_days = statistics.mean(total_days)
        std_total_days = statistics.stdev(total_days)

        mean_total_icum = statistics.mean(total_icum)
        std_total_icum= statistics.stdev(total_icum)

        d2 = {'country':country_name,'latitude':[target_latitude], 'longitude':[target_longitude],
        'mean_duration' : mean_duration, 'std_duration' : std_duration,
        'mean_intensity_mean' : mean_intensity_mean, 'std_intensity_mean' : std_intensity_mean,
        'mean_intensity_var' : mean_intensity_var, 'std_intensity_var' : std_intensity_var,
        'mean_intensity_max' : mean_intensity_max, 'std_intensity_max' : std_intensity_max,
        'mean_intensity_cumulative' : mean_intensity_cumulative, 'std_intensity_cumulative' : std_intensity_cumulative,
        'mean_rate_onset' : mean_rate_onset, 'std_rate_onset' : std_rate_onset,
        'mean_rate_decline' : mean_rate_decline, 'std_rate_decline' : std_rate_decline,
        'mean_total_days' : mean_total_days, 'std_total_days' : std_total_days,
        'mean_total_icum' : mean_total_icum, 'std_total_icum' : std_total_icum,

        'sig_trend_per_decade_count' : sig_trend_per_decade_count,
        'sig_trend_per_decade_duration' :sig_trend_per_decade_duration,
        'sig_trend_per_decade_intensity_max' :sig_trend_per_decade_intensity_max,
        'sig_trend_per_decade_intensity_max_max' : sig_trend_per_decade_intensity_max_max,
        'sig_trend_per_decade_intensity_mean' : sig_trend_per_decade_intensity_mean,
        'sig_trend_per_decade_intensity_var' :sig_trend_per_decade_intensity_var,
        'sig_trend_per_decade_intensity_cumulative' : sig_trend_per_decade_intensity_cumulative,
        'sig_trend_per_decade_rate_onset' :sig_trend_per_decade_rate_onset,
        'sig_trend_per_decade_rate_decline' : sig_trend_per_decade_rate_decline,
        'sig_trend_per_decade_total_days' : sig_trend_per_decade_total_days,
        'sig_trend_per_decade_total_icum' :sig_trend_per_decade_total_icum,

        'sum_total_days' :np.sum(total_days),
        'sum_total_icum' : np.sum(total_icum),
        'number_of_events' :events,

        'count_trend_per_decade' : count_trend_per_decade,
        'duration_trend_per_decade' : duration_trend_per_decade,
        'intensity_max_trend_per_decade' :  intensity_max_trend_per_decade,
        'intensity_max_max_trend_per_decade' : intensity_max_max_trend_per_decade,
        'intensity_mean_trend_per_decade' :   intensity_mean_trend_per_decade,
        'intensity_var_trend_per_decade' : intensity_var_trend_per_decade,
        'intensity_cumulative_trend_per_decade' :  intensity_cumulative_trend_per_decade,
        'rate_onset_trend_per_decade' : rate_onset_trend_per_decade,
        'rate_decline_trend_per_decade' :  rate_decline_trend_per_decade,
        'total_days_trend_per_decade' : total_days_trend_per_decade,
        'total_icum_trend_per_decade' : total_icum_trend_per_decade,


        'count_mean_per_year' : count_mean_per_year,
        'duration_mean_per_year' : duration_mean_per_year,
        'intensity_mean_mean_per_year' : intensity_mean_mean_per_year,
        'intensity_max_max_mean_per_year' : intensity_max_max_mean_per_year,
        'intensity_mean_mean_per_year' : intensity_mean_mean_per_year,
        'intensity_var_mean_per_year' : intensity_var_mean_per_year,
        'intensity_cumulative_mean_per_year' : intensity_cumulative_mean_per_year,
        'rate_onset_mean_per_year' : rate_onset_mean_per_year,
        'rate_decline_mean_per_year' : rate_decline_mean_per_year,
        'total_days_mean_per_year' : total_days_mean_per_year,
        'total_icum_mean_per_year' : total_icum_mean_per_year

        }

        df2 = pd.DataFrame(data=d2)
        df2.to_csv('/media/shilpa/Expansion/coastal_tokelau/%s_mean_trends_mhws_at_%s_%s.csv'%(country_name,target_latitude,target_longitude))

