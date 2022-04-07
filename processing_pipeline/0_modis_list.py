import pickle

from modis_tools import delete_no_overlap, UngriddedData_from_data_frame
import glob
import calendar
from multiprocessing import Pool
import sys

observation_window=[[-50., 50.],[ -90.,20.]]

months=['01','02','03','04','05','06','07','08','09','10','11','12']

def find_modis_granules_aqua(year, month, i):
    day=str(i+1).zfill(2)
    modis_dir = '/neodc/modis/data/MYD06_L2/collection61/{}/{}/{}/'
    modis_daily = sorted(glob.glob(modis_dir.format(year, month, day) + '*.hdf'))
    
    try: 
        end_day = glob.glob(modis_dir.format(year, month, day) + '*.2100.*.hdf')[0]
        end_day = modis_daily.index(end_day)
    except IndexError: end_day = -1 #glob.glob(modis_dir.format(year, month, day) + '*.2300.*.hdf')[0]
    # end_day = modis_daily.index(end_day)
    try: 
        start_time=glob.glob(modis_dir.format(year, month, day) + '*.0850.*.hdf')[0]
        start_time=modis_daily.index(start_time)
    except IndexError: 
        start_time=0 #glob.glob(modis_dir.format(year, month, day) + '*.0800.*.hdf')[0]    
    
    modis_daily = modis_daily[start_time:end_day]
    modis_daily=delete_no_overlap(modis_daily, observation_window)
    return modis_daily

def find_modis_granules_terra(year, month, i):
    day=str(i+1).zfill(2)
    modis_dir = '/neodc/modis/data/MOD06_L2/collection61/{}/{}/{}/'
    
    modis_daily = sorted(glob.glob(modis_dir.format(year, month, day) + '*.hdf'))
    try: 
        end_day = glob.glob(modis_dir.format(year, month, day) + '*.1800.*.hdf')[0]
        end_day = modis_daily.index(end_day)
    except IndexError: end_day = -1 #glob.glob(modis_dir.format(year, month, day) + '*.2000.*.hdf')[0]
    
    try: 
        start_time=glob.glob(modis_dir.format(year, month, day) + '*.0550.*.hdf')[0]
        start_time=modis_daily.index(start_time)
    except IndexError: start_time=0 #glob.glob(modis_dir.format(year, month, day) + '*.0500.*.hdf')[0]

    modis_daily = modis_daily[start_time:end_day]
    modis_daily=delete_no_overlap(modis_daily, observation_window)
    return modis_daily


year=sys.argv[2]
index=int(sys.argv[1])
months=['00','01','02','03','04','05','06','07','08','09','10','11','12']
month=months[index]

modis_list={'aqua':[],'terra':[]}

for i in range(calendar.monthrange(int(year), int(month))[1]):
    print(i)
    modis_list['aqua'].append(find_modis_granules_aqua(year, month, i))
    modis_list['terra'].append(find_modis_granules_terra(year, month, i))

f = open("../modis_tiles/modis_swaths_intersect_{}_{}.pkl".format(year, month),"wb")
pickle.dump(modis_list,f)
f.close()