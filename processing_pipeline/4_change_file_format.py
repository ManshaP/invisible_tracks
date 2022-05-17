# Name: Post-processing the data
#
# Description: Python code adding derived data fields and summarising the outputs
#
# Inputs: collocated trajectories and cloud properties
#
# Output: trajectories and cloud properties monthly, aqua and terra, with added droplet number concentration, EIS, boolean masks for the regions of interest 
# 
# Modes: normal, null experiment
#
# Libraries: numpy, xarray, pandas, calendar, datetime
#-------------------------------------------------------------------------------------------


from __future__ import division, unicode_literals, print_function  # for compatibility with Python 2 and 3
# import numba
import glob
import numpy as np
import sys
import pandas as pd
from pandas import DataFrame, Series  # for convenience
import datetime
import calendar
from global_land_mask import globe
from multiprocessing import Pool
import xarray as xr
import tables

year=2019

def condi(h5,observation_window): 
        cond = ((h5['latitude_1']<observation_window[0][1])&
        (h5['longitude_1']<observation_window[1][1])&
        (h5['latitude_1']>observation_window[0][0])&
        (h5['longitude_1']>observation_window[1][0])&
        (h5['latitude_3']<observation_window[0][1])&
        (h5['longitude_3']<observation_window[1][1])&
        (h5['latitude_3']>observation_window[0][0])&
        (h5['longitude_3']>observation_window[1][0])&
        (h5['latitude']<observation_window[0][1])&
        (h5['longitude']<observation_window[1][1])&
        (h5['latitude']>observation_window[0][0])&
        (h5['longitude']>observation_window[1][0]))
        return cond

def subsume_month(month):

    h5=pd.DataFrame()
    print(month)
    # filename = '/gws/nopw/j04/eo_shared_data_vol2/scratch/pete_nut/colloc_lines/atlantic/aqua_cth_{}{}{}.h5'
    # filename_terra = '/gws/nopw/j04/eo_shared_data_vol2/scratch/pete_nut/colloc_lines/atlantic/terra_cth_{}{}{}.h5'
    filename = '/gws/nopw/j04/eo_shared_data_vol1/satellite/.modistracks/{}/ctt_aqua_{}{}.h5'
    filename_terra = '/gws/nopw/j04/eo_shared_data_vol1/satellite/.modistracks/{}/ctt_terra_{}{}.h5'
    # filename = '/gws/nopw/j04/eo_shared_data_vol1/satellite/.modistracks/null/{}/aqua_{}{}.h5'
    # filename_terra='/gws/nopw/j04/eo_shared_data_vol1/satellite/.modistracks/null/{}/terra_{}{}.h5'

    # filename = '/gws/nopw/j04/eo_shared_data_vol2/scratch/pete_nut/colloc_lines/null/{}/aqua_{}{}.h5'
    # filename_terra = '/gws/nopw/j04/eo_shared_data_vol2/scratch/pete_nut/colloc_lines/null/{}/terra_{}{}.h5'
    EIS = xr.open_dataset('/gws/nopw/j04/eo_shared_data_vol1/satellite/.modistracks/reanalysisx/globEIS_{}'.format(year))
    # EIS = xr.open_dataset('/gws/nopw/j04/eo_shared_data_vol2/scratch/pete_nut/reanalysisx/globEIS_{}'.format(year))
    print(EIS.sel(time=['0h']))
    EIS = EIS.assign_coords(longitude=(((EIS.longitude + 180) % 360) - 180))
    for i in range(calendar.monthrange(int(year), int(month))[1]):
        
        if i+1<10: day='0{}'.format(i+1)
        else: day=str(i+1)
        print(month, day)
        try:
            new = pd.read_hdf(filename.format(year,month, day), key='df')
            # drop some pure nan lines from the files
            new = new.dropna(how='all', subset=['LWP','LWP_1', 'LWP_3'])
            new['terra']=0
            new.overpass=pd.to_datetime(new.overpass)
            interp_EIS = EIS.interp(time=xr.DataArray(new.overpass, dims='obs'), longitude=xr.DataArray(new.longitude, dims='obs'), latitude=xr.DataArray(new.latitude, dims='obs'))

            new['EIS'] = interp_EIS.t.values
            h5=pd.concat([h5, new])

        except (KeyError,FileNotFoundError, tables.exceptions.HDF5ExtError, AttributeError):
            print('no aqua data for:',day,'/',month,year)

        try:
            new = pd.read_hdf(filename_terra.format(year,month, day), key='df')
            new = new.dropna(how='all', subset=['LWP','LWP_1', 'LWP_3'])
            new['terra']=1
            new.overpass=pd.to_datetime(new.overpass)
            interp_EIS = EIS.interp(time=xr.DataArray(new.overpass, dims='obs'), longitude=xr.DataArray(new.longitude, dims='obs'), latitude=xr.DataArray(new.latitude, dims='obs'))
            new['EIS'] = interp_EIS.t.values
            h5=pd.concat([h5, new])

        except (KeyError, FileNotFoundError, tables.exceptions.HDF5ExtError,AttributeError):
            print('no terra data for:',day,'/',month, year)
    h5['ocean'] = np.logical_not(globe.is_land(h5.latitude.values, h5.longitude.values)) & np.logical_not(globe.is_land(h5.latitude.values+1, h5.longitude.values)) & np.logical_not(globe.is_land(h5.latitude.values, h5.longitude.values+1)) & np.logical_not(globe.is_land(h5.latitude.values-1, h5.longitude.values)) & np.logical_not(globe.is_land(h5.latitude.values, h5.longitude.values-1)) 
    h5=h5.set_index(pd.to_datetime(h5.index.values.astype("datetime64[ns]")))
    h5.overpass=pd.to_datetime(h5.overpass)
    # add Nd calculated from COT and effective radius
    h5['Nd'] = (1.37e-5 * (h5.r_eff * 1e-6)**(-5/2)*h5.COT**(1/2) / 1e6)
    h5['Nd_1'] = (1.37e-5 * (h5.r_eff_1 * 1e-6)**(-5/2)* h5.COT_1 **(1/2) / 1e6)
    h5['Nd_3'] = (1.37e-5 * (h5.r_eff_3 * 1e-6)**(-5/2)* h5.COT_3 **(1/2) / 1e6)
    #  add markers for the regions of interest
    h5['chil'] = condi(h5, [[-30., -17.],[ -82.5,-72.5]])
    h5['azor'] = condi(h5,[[30, 50],[ -40,-10]])
    h5['cver'] = condi(h5,[[10., 30.],[ -50,-20]])
    h5['ango'] = condi(h5,[[-30., -10.],[ -10.,15.]])
    # add the difference between emission time and satellite overpass (how long the emissions have been advected for)
    h5['hours_diff']=pd.to_timedelta((h5.overpass-h5.index).values).seconds//3600
    
    print('writing')
    h5.to_csv('/gws/nopw/j04/eo_shared_data_vol1/satellite/.modistracks/{}/ctp_n_{}'.format(year,month))
    # h5.to_csv('/gws/nopw/j04/eo_shared_data_vol1/satellite/.modistracks/null/{}/null_{}'.format(year,month))

    print(month, 'done')


with Pool(12) as p:
    p.map(subsume_month, ['01','02','03','04','05','06','07','08','09','10','11','12'])
