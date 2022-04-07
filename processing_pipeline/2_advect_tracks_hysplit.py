# Name: HYSPLIT Lagrangian Trajectory Calculation
#
# Description: Simple python code to run HYSPLIT and output trajectory file
#
# Modes:
# FORWARD trajectory, BACKWARD trajectory, and ALLWARDS trajectory (combines them)
#
# Inputs:
# latitude, longitude, height, trajectory duration, outputpath, HYSPLIT executable path, meteorology path
#
# Output: tdump file
#
# Libraries: numpy, os and datetime
#
# Notes: This version outputs a trajectory location every 15 minutes and uses
# NOAA ARL generated meteorological files for RP, GDAS, ERA5 or MERRA2.
# For example of running code see RUN HYSPLIT trajectory and make modifications as necessary
#
# Current example produces a 12 hour backward and forward trajectory from 45 N 120 W on June 12th 2016 at 15:30 UTC
# The two trajectories are "stitched" together. It is straightforward to run just a forward or back trajectory by negating the 'A' in trajHours variable.
#-------------------------------------------------------------------------------------------
import numpy as np
import os
import datetime
import time
import pandas as pd
import glob
import sys
import calendar
from multiprocessing import Pool


def create_hysplit_control_file_multi(year,month,day,hour,minute,lats,lons,MLHeight,trajHours,path,workingpath,metpath):

    nTraj = len(lats)
    print('YEAR: ',year)
    print('MONTH: ',month)
    print('DAY: ',day)
    print('HOUR: ',hour)
    print('MINUTE: ',minute)
    print('HYSPLIT Path: ',path)
    print('Number of trajectories: ',nTraj)
    #print('Height: ',MLHeight)
    #print('Longitude: ',lons)
    #print('Latitude: ',lats)
    
    #Fetch reanalysis files
    if metpath.find("RP") > 0:  #RP reanalysis path used
        metfiles = fetch_rp_reanalysis(month,year,metpath)

    if metpath.find("GDAS") > 0:  #GDAS reanalysis path used
        metfiles = fetch_gdas_reanalysis(day,month,year,metpath)

    if metpath.find("ERA5") > 0:  #ERA5 reanalysis path
        metfiles = fetch_era5_reanalysis(day,month,year,metpath,trajHours)
    
    if metpath.find("MERRA2") > 0:  #ERA5 reanalysis path
        metfiles = fetch_merra2_reanalysis(day,month,year,metpath,trajHours)
    
    nmetfiles = len(metfiles)
    
    #---------------------------------------------------------------
    # Create HYSPLIT CONTROL FILE
    #---------------------------------------------------------------
    fName = workingpath+'CONTROL'
    print('CONTROL FILE: ',fName)

    f = open(fName, 'w')
    f.writelines( str(year-2000).zfill(2)+' '+str(month).zfill(2)+' '+str(day).zfill(2)+' '+str(hour).zfill(2) +' '+str(minute).zfill(2)+"\n" ) #time
    f.writelines( str(nTraj) +"\n") #number of trajectories
    #starting location of each trajectory
    for i in range(nTraj):
        f.writelines("{:6.2f}".format((float(lats[i])))+'  '+"{:6.2f}".format((float(lons[i])))+'  '+"{:6.2f}".format(float(MLHeight[i])) +"\n")
    f.writelines( str(int(trajHours)) +"\n")
    f.writelines( str(0) + "\n") #non-isobaric trajectory
    f.writelines( str(10000.0) + "\n") #top of model in meters
    
    #Meteorological files
    f.writelines( str(nmetfiles) + "\n") #number of meteorologial files
    for i in range(nmetfiles):
        if len(os.path.dirname(metfiles[0])) == 0:
            f.writelines(metpath + "\n")
            f.writelines(metfiles[i] + "\n")
        if len(os.path.dirname(metfiles[0])) > 0:
            f.writelines(os.path.dirname(metfiles[i])+'/' + "\n")
            f.writelines(os.path.basename(metfiles[i]) + "\n")
    f.writelines(workingpath + "\n")
    f.writelines('tdump' + "\n")
    f.close()

    errval=0
    return errval


#---------------------------------------------------------------
#Function to fetch ERA5 reanalysis files
#---------------------------------------------------------------
def fetch_era5_reanalysis(day,month,year,metpath,trajHours):
    if trajHours > 0: t0 = datetime.datetime(year,month,day)
    if trajHours > 0: t1 = t0+datetime.timedelta(hours=trajHours)
    if trajHours < 0: t1 = datetime.datetime(year,month,day)
    if trajHours < 0: t0 = t1-datetime.timedelta(hours=abs(trajHours))
    nDays = int(np.ceil(((t1-t0).total_seconds())/86400.))+1
    metfiles = []
    metjday  = []
    for file in os.listdir(metpath):
        if file.endswith(".ARL"):
            metfiles.append(file)
            metjday.append(datetime.datetime(int(file[5:9]),int(file[9:11]),int(file[11:13]),0,0))
    #Sort by time
    sortID = np.asarray(np.argsort(metjday))
    metfiles = np.asarray(metfiles)
    metfiles = metfiles[sortID]
    fct = sortID.shape[0]
    metyears = [metjday[sortID[i]].year for i in range(len(metjday))]
    metmonths= [metjday[sortID[i]].month for i in range(len(metjday))]
    metdays= [metjday[sortID[i]].day for i in range(len(metjday))]
    metjday2 = [metjday[sortID[i]] for i in range(len(metjday))]
    index, = np.where( (np.asarray(metjday2) >= t0) & (np.asarray(metjday2) <= t1) )
    indices = [index[0]-1,index[len(index)-1]+1]
    return metfiles[index[0]-1:index[len(index)-1]+2]




#Function to run HYSPLIT particle trajectory code
def traj(latPt,lonPt,hPt,year,month,day,hour,minute,trajHours,path_hysplit_code,metpath,outpath,setupFile,ASCFile,tName):
    cwd = os.getcwd()
    path_working = outpath+'working/'
    os.system('cp '+setupFile+' '+path_working+'/SETUP.CFG')
    os.system('cp '+ASCFile+' '+path_working+'/ASCDATA.CFG')
    os.chdir(path_working)

    lat  = np.asarray(latPt)
    lon  = np.asarray(lonPt)
    MLHeight = np.asarray(hPt)
    create_hysplit_control_file_multi(int(year),int(month),int(day),int(hour),int(minute),lat,lon,MLHeight,trajHours,path_hysplit_code,path_working,metpath)
    
    #Run HYSPLIT code
    os.system(path_hysplit_code+'/exec/hyts_std')
    os.system('chmod 777 '+path_working+'tdump')
    print('mv '+path_working+'tdump '+outpath+tName)
    os.system('mv '+path_working+'tdump '+outpath+tName)
    os.chdir(cwd)

def run_trajectory(latPt,lonPt,yyyy,mm,dd,hh,mn,hPt,trajHours,outpath,path_hysplit_code,metpath,setupFile,ASCFile,tName):
    os.system('mkdir -p '+outpath)
    path_working = outpath+'working/'
    os.system('mkdir -p '+path_working)
    
    #Run both forward and backward trajectories
    if trajHours.find('A') == 0:
        #Backward
        BACK_trajHours = int( (-1.)*abs( int( trajHours[1:] ) ) )
        BACK_outpath = outpath + '/bwd/'
        os.system('mkdir -p '+BACK_outpath)
        path_working = BACK_outpath+'working/'
        os.system('mkdir -p '+path_working)
        BACK_tName = 'tdump_traj_BACKWARD_test.txt'
        traj(latPt,lonPt,hPt,yyyy,mm,dd,hh,mn,BACK_trajHours,path_hysplit_code,metpath,BACK_outpath,setupFile,ASCFile,BACK_tName)
        
        #Forward
        FOR_trajHours = int( abs( int( trajHours[1:] ) ) )
        FOR_outpath = outpath + '/fwd/'
        os.system('mkdir -p '+FOR_outpath)
        path_working = FOR_outpath+'working/'
        os.system('mkdir -p '+path_working)
        FOR_tName = 'tdump_traj_FORWARD_test.txt'
        traj(latPt,lonPt,hPt,yyyy,mm,dd,hh,mn,FOR_trajHours,path_hysplit_code,metpath,FOR_outpath,setupFile,ASCFile,FOR_tName)
        
        #Combine files
        cFile = outpath+tName
        fw = open(cFile,'w')
        
        file1 = open(BACK_outpath+BACK_tName,'r')
        BWD_Lines = file1.readlines()
        file1.close()
        file1 = open(FOR_outpath+FOR_tName,'r')
        FWD_Lines = file1.readlines()
        file1.close()
        
        fwdI = 0
        for i in range(len(FWD_Lines)):
            fw.write(FWD_Lines[i])
            if FWD_Lines[i].find('FORWARD') >= 0:
                fw.write(FWD_Lines[i+1])
                fw.write(FWD_Lines[i+2])
                fwdI = i+2
                break
        bwdI = 0
        for i in range(len(BWD_Lines)):
            if BWD_Lines[i].find('BACKWARD') >= 0:
                bwdI= i + 2
                break
        #Fetch lines
        for i in reversed(range(bwdI+1,len(BWD_Lines))):
            fw.write(BWD_Lines[i])
        for i in range(fwdI+2,len(FWD_Lines)):
            fw.write(FWD_Lines[i])
        fw.close()
        print('combined_file: '+cFile)
    else:
        traj(latPt,lonPt,hPt,yyyy,mm,dd,hh,mn,int(trajHours),path_hysplit_code,metpath,outpath,setupFile,ASCFile,tName)



def read_trajectory_file(TRAJECTORY_FILE,fileFlag=0, rmfile=False):
    file1 = open(TRAJECTORY_FILE,'r')
    L = file1.readlines()
    file1.close()
    traj = {'particle':[],
#             'year':[],'month':[],'day':[],'hour':[],'minute':[],
            'jday':[],'lat':[],'lon':[],'alt':[],'timestep':[]}
    for i in range(len(L)):
        #if (L[i].find('BACKWARD') > 0) | (L[i].find('FORWARD') > 0):
        if (L[i].find('PRESSURE') > 0):
            break
    for i in range(i+1,len(L)):
        traj['particle'].append    ( int( (L[i].split())[0] ) )
#         traj['year'].append    ( int( (L[i].split())[2] )+2000 )
#         traj['month'].append   ( int( (L[i].split())[3] ) )
#         traj['day'].append     ( int( (L[i].split())[4] ) )
#         traj['hour'].append    ( int( (L[i].split())[5] ) )
#         traj['minute'].append  ( int( (L[i].split())[6] ) )
        traj['timestep'].append( float( L[i].split()[8] ) )
        traj['lat'].append( float( (L[i].split())[9] ) )
        traj['lon'].append( float( (L[i].split())[10] ) )
        traj['alt'].append( float( (L[i].split())[11] ) )
        traj['jday'].append( datetime.datetime( int( (L[i].split())[2] )+2000, int( (L[i].split())[3] ), int( (L[i].split())[4] ), int( (L[i].split())[5] ), int( (L[i].split())[6] )))
        
    return traj

def hysplit_traj(latInit, lonInit, jday):#, overpass_time):
    # lat and lon arrays, time and overpass_time datetime
    t0=time.time()

    outpath='/home/users/pete_nut/test/'
    os.system('mkdir -p '+outpath)
    path_hysplit_code = '/home/users/pete_nut/hysplit.v5.0.0_CentOS/'
    metpath = '/gws/nopw/j04/eo_shared_data_vol1/reanalysis/ARL_noaa/reanalysis_data/ERA5/'
    setupFile = '/home/users/pete_nut/IV_shiptracks/matt_traj_code/SETUP_traj_lowcloud_15min.CFG'
    ASCFile  = '/home/users/pete_nut/IV_shiptracks/matt_traj_code/ASCDATA_JASMIN.CFG'

    year   = str( jday.year  ).zfill(4)
    month  = str( jday.month ).zfill(2)
    day    = str( jday.day ).zfill(2)
    hour   = str( jday.hour ).zfill(2)
    minute = str( jday.minute ).zfill(2)
    height = np.full(len(latInit),'20').flatten() # not midway BL, rather 20m height |midway through the PBL (for now)-

    trajHours=24
    trajHours=str(trajHours)
    
    tName = 'tdump_traj_FORWARD_v0.txt'
    EXP = year+month+day+hour+minute
    rootPath = outpath+EXP+'/'
    run_trajectory(latInit,lonInit,year,month,day,hour,minute,height,
                   trajHours,rootPath,path_hysplit_code,metpath,setupFile,ASCFile,tName)
    
    trajectory = read_trajectory_file( rootPath + tName )
    os.system('rm -r '+rootPath)

    t1=time.time()
    print('time elapsed = ', t1-t0)
    return trajectory


def advect_emissions_day(i, month, year):
    if i+1<10: day='0{}'.format(i+1)
    else: day=str(i+1)

    # find features and trajectories
    print('working on day:', day, month, year)
    t1=pd.read_csv('/gws/nopw/j04/eo_shared_data_vol2/scratch/pete_nut/emissions_tracked/{}/{}_{}'.format(year, month, day))
    t1=t1.set_index(pd.to_datetime(t1['hour']))
    t1=t1.drop(columns=['hour']).sort_index()
    for i, t in enumerate(np.unique(t1.index)):
        print(i+1, '/',len(np.unique(t1.index)))
        lat_init=t1.sort_values('particle').latitude.where(t1.sort_values('particle').index==t).dropna()
        lon_init=t1.sort_values('particle').longitude.where(t1.sort_values('particle').index==t).dropna()
        start_time=t1.sort_values('particle').index.where(t1.sort_values('particle').index==t).dropna()[0]
        try: 
            adv_pos = hysplit_traj(lat_init, lon_init, start_time)
        except FileNotFoundError: continue
        adv_pos = pd.DataFrame.from_dict(adv_pos)
        adv_pos.to_csv('/gws/nopw/j04/eo_shared_data_vol1/satellite/.modistracks/hysplit_adv_20m_non_isob/{}_{}_{}_{}:{}'.format(
        year, month, day, str(start_time.hour).zfill(2), str(start_time.minute).zfill(2)))

def parallelise_advection(i):
    index=int(sys.argv[1])
    month=months[index]
    advect_emissions_day(i, month, year)

months=['00','01','02','03','04','05','06','07','08','09','10','11','12']
year = sys.argv[2]
index = int(sys.argv[1])
month = months[index]

with Pool(12) as p:
    p.map(parallelise_advection, range(calendar.monthrange(int(year), int(month))[1]) )

