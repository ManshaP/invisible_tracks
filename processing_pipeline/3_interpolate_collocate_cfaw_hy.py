# Name: Collocate the advected trajectories with MODIS satellite data
#
# Description: Python code using the cis library to construct the trajectory at the time of the satellite overpass and collocate to the satellite imagery
#
# Inputs: advected emissions of the last 24h, MODIS data
#
# Output: collocated trajectories and cloud properties
#
# Modes: normal/null experiment, Aqua/Terra
# 
# Libraries: numpy, cis, xarray, pandas, pickle, calendar, datetime
#-------------------------------------------------------------------------------------------


from __future__ import division, unicode_literals, print_function  # for compatibility with Python 2 and 3
import numba
import numpy as np
import pandas as pd
from pandas import DataFrame, Series  # for convenience
import time
import sys
import tables
import pims
import pickle
import trackpy as tp
import xarray as xr
from cis.data_io.ungridded_data import UngriddedDataList
import os
os.environ['CIS_PLUGIN_HOME'] = '/home/users/pete_nut/plugins/'
from cis import read_data, read_data_list, get_variables
from modis_tools import delete_no_overlap, UngriddedData_from_data_frame
import glob
import calendar
from multiprocessing import Pool
import datetime
import matplotlib.pyplot as plt
import time

#np.unique(t1.particle.values)
#interp.loc[0].iloc[-1].x
def get_offset(dataframe, particle_ids, dist):
    offsets=np.zeros((int(particle_ids.max())+1,2))
    for i in particle_ids:
        i=int(i)
        x_start=dataframe.loc[i].iloc[0].longitude
        #print(dataframe.loc[i])
        y_start=dataframe.loc[i].iloc[0].latitude
        x_end=dataframe.loc[i].iloc[-1].longitude
        y_end=dataframe.loc[i].iloc[-1].latitude
        x_dir=x_end-x_start
        y_dir=y_end-y_start
        veclen=np.sqrt(x_dir**2+y_dir**2)
        normal_vec_x=-y_dir*dist/veclen
        normal_vec_y=x_dir*dist/veclen
        if normal_vec_x<0:
            offsets[i]=[-normal_vec_x,-normal_vec_y]
        else: offsets[i]=[normal_vec_x,normal_vec_y]
    return offsets


pd.options.mode.chained_assignment = None 

def collocate_day(i):

    if i+1<10: day='0{}'.format(i+1)
    else: day=str(i+1)

    # get the emissions before wind advection
    print('working on day:', day, month, year)
    if year=='2019': t1=pd.read_csv('/gws/nopw/j04/eo_shared_data_vol2/scratch/pete_nut/emissions_tracked/{}/2019_{}_{}'.format(year, month, day))
    else: t1=pd.read_csv('/gws/nopw/j04/eo_shared_data_vol2/scratch/pete_nut/emissions_tracked/{}/{}_{}'.format(year, month, day))

    day_before = pd.to_datetime(year+ month+ day, format='%Y%m%d') - pd.Timedelta('1d')
    try:
        if year=='2019': 
            t1_b=pd.read_csv('/gws/nopw/j04/eo_shared_data_vol2/scratch/pete_nut/emissions_tracked/{}/2019_{}_{}'.format(str(day_before.year), str(day_before.month).zfill(2), str(day_before.day).zfill(2)))   
        else: 
            t1_b=pd.read_csv('/gws/nopw/j04/eo_shared_data_vol2/scratch/pete_nut/emissions_tracked/{}/{}_{}'.format(str(day_before.year), str(day_before.month).zfill(2), str(day_before.day).zfill(2)))    
 
        t1_b.particle=t1_b.particle+10000
        t1=pd.concat( [t1_b, t1])
    except FileNotFoundError: None

    t1=t1.set_index(pd.to_datetime(t1['hour']))
    t1=t1.drop(columns=['hour']).sort_index()
    
    if year=='2019':
        infile = open('/home/users/pete_nut/IV_shiptracks/modis_tiles/atlantic_modis_swaths_intersect_{}_{}.pkl'.format(year, month),'rb')
    else: 
        infile = open('/home/users/pete_nut/IV_shiptracks/modis_tiles/modis_swaths_intersect_{}_{}.pkl'.format(year, month),'rb')
    swath_dict = pickle.load(infile)
    infile.close()
    today_modis=[]
    if terra: swaths_flat=[item for subl in swath_dict['terra'] for item in subl]
    else: swaths_flat=[item for subl in swath_dict['aqua'] for item in subl]
    for filename in swaths_flat:
        split_f=filename.split('/')
        if null: 
            if (split_f[6]==str(day_before.year).zfill(2)) & (split_f[7]==str(day_before.month).zfill(2)) & (split_f[8]==str(day_before.day).zfill(2)):
                today_modis.append(filename)
        else:
            if (split_f[6]==year) & (split_f[7]==month) & (split_f[8]==day):
                today_modis.append(filename)
    print(len(today_modis),'modis files to collocate on day', day, '/',month,'/',year)
    modis_daily=today_modis
    collect_overpasses=pd.DataFrame()

    for filename in modis_daily:
        print(filename)
        t_0=time.time()
        overpass_hour_min=filename.split('.')[2]
        print(overpass_hour_min)


        # advect polluted trajectories with the wind

        time_of_overpass=pd.to_datetime(year+month+day+overpass_hour_min,  format='%Y%m%d%H%M')
        t_first_counted_emis=pd.to_datetime(year+month+day+'0000',  format='%Y%m%d%H%M')

        t1_adv=t1.copy().loc[slice(time_of_overpass-pd.Timedelta('24h'), time_of_overpass),:]

        t1_adv=t1_adv.sort_values('particle')
        for i, t in enumerate(np.unique(t1_adv.index)):

            print(i+1, '/',len(np.unique(t1_adv.index)))

            start_time=t1_adv.index.where(t1_adv.index==t).dropna()[0]
            if (time_of_overpass - start_time).total_seconds()<0 : continue
            else: None

            try:
                if ((year=='2019') or (year=='2018') or (year=='2017')):
                    filloc = '/gws/nopw/j04/eo_shared_data_vol2/scratch/pete_nut/hysplit_adv_20m_non_isob/{}_{}_{}_{}:{}'
                else: 
                    filloc = '/gws/nopw/j04/eo_shared_data_vol1/satellite/.modistracks/hysplit_adv_20m_non_isob/{}_{}_{}_{}:{}'
                #this csv file has the locations where emissions from year, month, day, start_time end up at any given moment within 24h in 5min steps

                traj_read=pd.read_csv(filloc.format(                
                year, str(start_time.month).zfill(2), str(start_time.day).zfill(2), str(start_time.hour).zfill(2), str(start_time.minute).zfill(2)))
                traj_read.jday=pd.to_datetime(traj_read.jday)

                traj_read=traj_read[traj_read.jday==time_of_overpass]

                if len(traj_read) != len(t1_adv[t1_adv.index==t]):
                    for i in range(1,len(t1_adv[t1_adv.index==t])+1):
                        if not np.isin(i, np.unique(traj_read.particle)) : 
                            print('missing particle', i)
                            new_row = {'particle':i, 'jday':time_of_overpass, 'lat':np.nan, 'lon':np.nan ,'alt':np.nan ,'timestep':np.nan }
                            #append row to the dataframe
                            traj_read=traj_read.append(new_row, ignore_index=True)
                            traj_read=traj_read.sort_values('particle')
                        else: None
                else: None
                t1_adv.latitude[t1_adv.index==t]=traj_read.lat.values
                t1_adv.longitude[t1_adv.index==t]=traj_read.lon.values

            except (FileNotFoundError,): #ValueError
                t1_adv.latitude[t1_adv.index==t]=np.nan
                t1_adv.longitude[t1_adv.index==t]=np.nan
                print('no hysplit data for {}'.format(start_time))


        interp_adv=t1_adv.groupby('particle').resample('5T').interpolate()
        interp_adv['overpass']=time_of_overpass


        # clean up
        interp_adv=interp_adv.drop(columns=['mass', 'y', 'x'])

        # subset the emission traj to the observation window in both time and space 
        mask=[(interp_adv['latitude']<observation_window[0][1])&
          (interp_adv['longitude']<observation_window[1][1])&
          (interp_adv['latitude']>observation_window[0][0])&
          (interp_adv['longitude']>observation_window[1][0])]

        interp_adv=interp_adv[mask[0]]#.loc[(slice(None),slice(t_first_counted_emis, time_of_overpass)),:]
        # sometimes, this leaves no tracks to look at
        if (len(interp_adv.dropna()) == 0) : 
            print('no overlapping data')
            continue


    #get modis data

        print('loading modis file') 
        mod_lr=read_data_list(filename,
                            ['Cloud_Effective_Radius',
                            'Cloud_Effective_Radius_Uncertainty', 
                            'Cloud_Phase_Optical_Properties',
                            'Cloud_Water_Path',
                            'Cloud_Optical_Thickness',
                            # 'cloud_top_pressure_1km'
                            'cloud_top_temperature_1km'
                            ], product='MOD06_HACK')

        print('done loading modis file',filename)
        # want liquid clouds with low uncertainty in r_eff estimate
        sat_mask=(mod_lr[2].data==2) & (mod_lr[1].data<=10)
        xy=UngriddedDataList([mod_lr[0][sat_mask], mod_lr[3][sat_mask], mod_lr[4][sat_mask], mod_lr[5][sat_mask]])

        modis_coord_box=[[mod_lr.coords('latitude')[0].points.min(), 
           mod_lr.coords('latitude')[0].points.max()],
           [mod_lr.coords('longitude')[0].points.min(),           
           mod_lr.coords('longitude')[0].points.max()]]

        mask_hysplit_from_modis=[(interp_adv['latitude']<modis_coord_box[0][1])&
          (interp_adv['longitude']<modis_coord_box[1][1])&
          (interp_adv['latitude']>modis_coord_box[0][0])&
          (interp_adv['longitude']>modis_coord_box[1][0])]

        interp_adv=interp_adv[mask_hysplit_from_modis[0]]
        # sometimes, this leaves no tracks to look at
        if (len(interp_adv.dropna()) == 0) : 
            print('no overlapping data')
            continue

        # get shifted trajectories for unpolluted counterfactual
        offsets=get_offset(interp_adv, np.unique(interp_adv.particle.values), 0.3)

        shifted_adv=interp_adv.copy()
        shifted_adv2=interp_adv.copy()

        shifted=[shifted_adv,interp_adv,shifted_adv2]
        n_traj = len(shifted)
        k = n_traj//2

        for i in np.unique(interp_adv.particle.values):
            i=int(i)
            for j, current_traj in enumerate(shifted):
                current_traj.loc[i].longitude = interp_adv.loc[i].longitude.values + (j-k) * offsets[i][0]
                current_traj.loc[i].latitude = interp_adv.loc[i].latitude.values + (j-k) * offsets[i][1]

        for j, current_traj in enumerate(shifted):
            current_traj.index = current_traj.index.droplevel(0)


        # wrangle datasets (need to be merged before collocating to save time)
        concatenated=pd.concat(shifted)

        # create cis UngriddedData object from DataFrame
        tracks_ungridded=UngriddedData_from_data_frame(
            concatenated,
            ['signal'])

        print('collocate emissions and modis')
        # collocate the modis and emission data and append to collection dataframe
        try:
            coll = xy.collocated_onto(tracks_ungridded[0], 'box', fill_value=-999,h_sep='10km')
            # 0=reff, 2=npoints, 3=LWP, 6=COT, 9=CTT
            coldat=UngriddedDataList([coll[0], coll[2], coll[3], coll[6], coll[9]]).as_data_frame()
            print(coldat)
            coldat=coldat.rename(
                columns={'Cloud Particle Effective Radius two-channel retrieval using band 7(2.1um) and either band 1(0.65um), 2(0.86um), or 5(1.2um)  (specified in Quality_Assurance_1km)from best points: not failed in any way, not marked for clear sky restoral':'r_eff',
                        'Number of points used to calculate the mean of Cloud Particle Effective Radius two-channel retrieval using band 7(2.1um) and either band 1(0.65um), 2(0.86um), or 5(1.2um)  (specified in Quality_Assurance_1km)from best points: not failed in any way, not marked for clear sky restoral':'npoints',
                        'Column Water Path two-channel retrieval using band 7(2.1um) and either band 1(0.65um), 2(0.86um), or 5(1.2um) (specified in Quality_Assurance_1km)from best points: not failed in any way, not marked for clear sky restoral':'LWP',
                        'Cloud Optical Thickness two-channel retrieval using band 7(2.1um) and either band 1(0.65um), 2(0.86um), or 5(1.2um)  (specified in Quality_Assurance_1km)from best points: not failed in any way, not marked for clear sky restoral':'COT',
                        'Cloud Top Temperature at 1-km resolution from LEOCAT, Temperature from Ancillary Data at Retrieved Cloud Top Pressure Level':'CTT'
                        }
            )

            middle = np.shape(coldat)[0]//n_traj
            saveframe=coldat[k*middle :(k+1) * middle].copy()

            saveframe['signal']=concatenated.signal[:middle].values
            saveframe['overpass']=concatenated.overpass[:middle].values
            saveframe['particle']=concatenated.particle[:middle].values
            print('number of paths:',n_traj)
            for n in range(n_traj):
                if n==n_traj//2: continue
                traj_name=str(n+1)
                saveframe['latitude_'+traj_name]=coldat.latitude[n*middle:(n+1)*middle]
                saveframe['longitude_'+traj_name]=coldat.longitude[n*middle:(n+1)*middle]
                saveframe['r_eff_'+traj_name]=coldat.r_eff[n*middle:(n+1)*middle]
                saveframe['LWP_'+traj_name]=coldat.LWP[n*middle:(n+1)*middle]
                saveframe['COT_'+traj_name]=coldat.COT[n*middle:(n+1)*middle]
                saveframe['CTT_'+traj_name]=coldat.CTT[n*middle:(n+1)*middle]
                saveframe['npoints_'+traj_name]=coldat.npoints[n*middle:(n+1)*middle]







            collect_overpasses=collect_overpasses.append(saveframe)
            print(len(saveframe.dropna()), 'data points collocated')
        except ValueError: 
            print('no data for'+filename)
            continue
        t_1=time.time()
        print('time elapsed for this modis file:', t_1-t_0)
    if terra:
        try:collect_overpasses.to_hdf(savepath_t.format(year, month,day), key='df')

        except tables.exceptions.HDF5ExtError: 
            time.sleep(60)
            collect_overpasses.to_hdf(savepath_t.format(year, month,day), key='df')
    else:
        try: collect_overpasses.to_hdf(savepath_a.format(year, month,day), key='df')
        except tables.exceptions.HDF5ExtError:
            time.sleep(60)
            collect_overpasses.to_hdf(savepath_t.format(year, month,day), key='df')



nullexp = sys.argv[4]
null = nullexp == 'null'
if null: 
    savepath_t="/gws/nopw/j04/eo_shared_data_vol1/satellite/.modistracks/null/{}/terra_{}{}.h5"
    savepath_a="/gws/nopw/j04/eo_shared_data_vol1/satellite/.modistracks/null/{}/aqua_{}{}.h5"
else: 
    savepath_t="/gws/nopw/j04/eo_shared_data_vol1/satellite/.modistracks/{}/ctt_terra_{}{}.h5"
    savepath_a="/gws/nopw/j04/eo_shared_data_vol1/satellite/.modistracks/{}/ctt_aqua_{}{}.h5"
print('null =',null)
year=sys.argv[3]

observation_window=[[-50., 50.],[ -90.,20.]]

index=int(sys.argv[1])
months=['00','01','02','03','04','05','06','07','08','09','10','11','12']
month=months[index]
terra = (sys.argv[2]=='Terra')
print('Terra=',terra)


with Pool(16) as p:
    p.map(collocate_day, range(calendar.monthrange(int(year), int(month))[1]) )
