from __future__ import division, unicode_literals, print_function  # for compatibility with Python 2 and 3
import numba
import sys
import numpy as np
import pandas as pd
import xarray as xr
import trackpy as tp
import calendar
import trackpy.diag
from multiprocessing import Pool


trackpy.diag.performance_report()
def x_to_lon(x):
    lon=x/100-90
    return lon
def y_to_lat(y):
    lat=-y/100
    return lat


def x_coarse_to_lon(x):
    lon=x/20-90
    return lon
def y_coarse_to_lat(y):
    lat=50-y/20
    return lat


def track_emissions_coarse(i):
    if i+1<10: day='0{}'.format(i+1)
    else: day=str(i+1)
    
    # find features and trajectories
    print('working on day:', day, month, year)
    monthly_emis=xr.open_mfdataset('/gws/nopw/j04/eo_shared_data_vol2/scratch/pete_nut/Emissions/'+year+'/SOx_*.nc', chunks={'latitude': 10, 'longitude': 10})
    monthly_emis=monthly_emis.SOx.loc[year+'-'+month]
    frames=monthly_emis
    f=tp.batch(frames[(int(day)-1)*24:int(day)*24].values, 7)

    #search_range argument as usual, along with an adaptive_stop argument, and optionally an adaptive_step.
    pred = tp.predict.NearestVelocityPredict()
    t = pred.link_df(f, search_range=20, adaptive_stop=4,adaptive_step=0.98,  memory=3)
    t1 = tp.filter_stubs(t,3)
    
    # Compare the number of particles in the unfiltered and filtered data.
    print('Before:', t['particle'].nunique())
    print('After:', t1['particle'].nunique())

    #clean up the dataset
    t1['hour']=pd.to_datetime(year+month+day)+pd.to_timedelta(t1.frame.values, unit='hours')
    t1['latitude']=y_coarse_to_lat(t1.y)
    t1['longitude']=x_coarse_to_lon(t1.x)
    t1=t1.set_index(t1['hour'])
    t1=t1.drop(columns=['hour','frame','size','ecc','raw_mass','ep']).sort_index()
    t1.to_csv('/gws/nopw/j04/eo_shared_data_vol2/scratch/pete_nut/emissions_tracked/{}/{}_{}'.format(year, month, day))
    
index=int(sys.argv[1])
year=sys.argv[2]
months=['00','01','02','03','04','05','06','07','08','09','10','11','12']
month=months[index]

with Pool(12) as p:
        p.map(track_emissions_coarse, range(calendar.monthrange(int(year), int(month))[1]) )
