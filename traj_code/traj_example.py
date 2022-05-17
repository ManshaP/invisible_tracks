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
    f.writelines( str(year-2000).zfill(2)+' '+str(month).zfill(2)+' '+str(day).zfill(2)+' '+str(hour).zfill(2) + "\n" ) #time
    f.writelines( str(nTraj) +"\n") #number of trajectories
    #starting location of each trajectory
    for i in range(nTraj):
        f.writelines("{:6.2f}".format((float(lats[i])))+'  '+"{:6.2f}".format((float(lons[i])))+'  '+"{:6.2f}".format(float(MLHeight[i])) +"\n")
    f.writelines( str(int(trajHours)) +"\n")
    f.writelines( str(1) + "\n") #isobaric trajectory
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
#Function to fetch RP reanalysis files
#---------------------------------------------------------------
def fetch_rp_reanalysis(month,year,metpath):
    metfiles = []
    metjday  = []
    for file in os.listdir(metpath):
        if file.endswith(".gbl"):
            metfiles.append(file)
            metjday.append(datetime.datetime(int(file[2:6]),int(file[6:8]),15,0,0))

    #Sort by time
    sortID = np.asarray(np.argsort(metjday))
    metfiles = np.asarray(metfiles)
    metfiles = metfiles[sortID]
    fct = sortID.shape[0]
    metyears = [metjday[sortID[i]].year for i in range(len(metjday))]
    metmonths= [metjday[sortID[i]].month for i in range(len(metjday))]

    #Select 3 met files closest in time to the input time
    nmetfiles = 3
    index = (np.where((np.asarray(metyears) == year) & (np.asarray(metmonths) == month))[0])[0]
    file0 = metfiles[ index-1] #prior month
    file1 = metfiles[ index  ]
    file2 = metfiles[ index+1] #next month
    metfiles = [file0,file1,file2]
    return metfiles


#---------------------------------------------------------------
#Function to fetch GDAS reanalysis files
#---------------------------------------------------------------
def fetch_gdas_reanalysis(day,month,year,metpath):
    monSTR = ['jan','feb','mar','apr','may','jun','jul','aug','sep','oct','nov','dec']
    monVal = [  1  ,  2  ,  3  ,  4  ,  5  ,  6  , 7   , 8   , 9   , 10  , 11  , 12  ]
    weekSTR = ['w1' ,'w2' ,'w3' ,'w4','w5']
    weekVal = [ 1   , 7   , 14   , 21  , 28]
    metfiles = []
    metjday  = []
    for file in os.listdir(metpath):
        if file.startswith("gdas1."):
            #print(file)
            monthName = file[6:9]
            TMPyear = int(file[9:11])+2000
            week = file[12:14]
            mValID=monSTR.index(monthName)
            TMPmonth=monVal[mValID]
            wValID=weekSTR.index(week)
            TMPday=weekVal[wValID]
            metfiles.append(file)
            metjday.append(datetime.datetime(TMPyear,TMPmonth,TMPday,0,0))
            #print(monthName,' ',week,' ',year,' ',month,' ',day)
            #print(datetime.datetime(year,month,day,0,0))


    #Sort by time
    sortID = np.asarray(np.argsort(metjday))
    metfiles = np.asarray(metfiles)
    metfiles = metfiles[sortID]
    fct = sortID.shape[0]
    metyears = [metjday[sortID[i]].year for i in range(len(metjday))]
    metmonths= [metjday[sortID[i]].month for i in range(len(metjday))]
    metdays  = [metjday[sortID[i]].day for i in range(len(metjday))]
    metjdays = [metjday[sortID[i]] for i in range(len(metjday))]

    #Select 5 met files closest in time to the input time
    nmetfiles = 5
    jdayRef= datetime.datetime(year,month,day)
    diffTime = np.asarray([(jdayRef-metjdays[i]).total_seconds()  for i in range(len(metjday))])
    index = (np.where( abs(diffTime) == min(abs(diffTime)) ))[0][0]
    file0 = metfiles[ index-2] #prior weeks
    file1 = metfiles[ index-1] #prior week
    file2 = metfiles[ index]   
    file3 = metfiles[ index+1] #next week
    file4 = metfiles[ index+2] #next weeks
    metfiles = [file0,file1,file2,file3,file4]
    return metfiles


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


#---------------------------------------------------------------
#Function to fetch MERRA2 reanalysis files
#---------------------------------------------------------------
def fetch_merra2_reanalysis(day,month,year,metpath,trajHours):
    if trajHours > 0: t0 = datetime.datetime(year,month,day)
    if trajHours > 0: t1 = t0+datetime.timedelta(hours=trajHours)
    if trajHours < 0: t1 = datetime.datetime(year,month,day)
    if trajHours < 0: t0 = t1-datetime.timedelta(hours=abs(trajHours))
    nDays = int(np.ceil(((t1-t0).total_seconds())/86400.))+1
    metfiles = []
    metjday  = []
    for iD in range(nDays):
        t = t0 + datetime.timedelta(days=iD)
        f = metpath + str(t.year).zfill(4) +'/'+ str(t.month).zfill(2) +'/'+ str(t.day).zfill(2) + '/' + 'MERRA2_'+str(t.year).zfill(4)+str(t.month).zfill(2)+str(t.day).zfill(2)+'.ARL'
        if os.path.exists(f): 
            metfiles.append(f)
            metjday.append( t )
        print(t.year,t.month,t.day)
    return metfiles



#Function to run HYSPLIT particle trajectory code
def traj(latPt,lonPt,hPt,year,month,day,hour,minute,trajHours,path_hysplit_code,metpath,outpath,setupFile,ASCFile,tName):
    cwd = os.getcwd()
    path_working = outpath+'working/'
    os.system('cp '+setupFile+' '+path_working+'/SETUP.CFG')
    os.system('cp '+ASCFile+' '+path_working+'/ASCDATA.CFG')
    os.chdir(path_working)
    
    lat  = np.asarray([latPt])
    lon  = np.asarray([lonPt])
    MLHeight = np.asarray([hPt])
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
    traj = {'n':[],'year':[],'month':[],'day':[],'hour':[],'minute':[],'jday':[],'lat':[],'lon':[],'alt':[],'timestep':[]}
    for i in range(len(L)):
        if (L[i].find('BACKWARD') > 0) | (L[i].find('FORWARD') > 0):
            break
    for i in range(i+3,len(L)):
        traj['n'].append    ( int( (L[i].split())[0] ) )
        traj['year'].append    ( int( (L[i].split())[2] )+2000 )
        traj['month'].append   ( int( (L[i].split())[3] ) )
        traj['day'].append     ( int( (L[i].split())[4] ) )
        traj['hour'].append    ( int( (L[i].split())[5] ) )
        traj['minute'].append  ( int( (L[i].split())[6] ) )
        traj['timestep'].append( float( L[i].split()[8] ) )
        traj['lat'].append( float( (L[i].split())[9] ) )
        traj['lon'].append( float( (L[i].split())[10] ) )
        traj['alt'].append( float( (L[i].split())[11] ) )
        traj['jday'].append( datetime.datetime( int( (L[i].split())[2] )+2000, int( (L[i].split())[3] ), int( (L[i].split())[4] ), int( (L[i].split())[5] ), int( (L[i].split())[6] )))
        
    return traj






#-----------------------------------------------------------------------------
# RUN HYSPLIT trajectory
#-----------------------------------------------------------------------------
outpath='/home/users/pete_nut/test/'
os.system('mkdir -p '+outpath)
path_hysplit_code = '/home/users/pete_nut/hysplit.v5.0.0_CentOS/'
metpath = '/gws/nopw/j04/eo_shared_data_vol1/reanalysis/ARL_noaa/reanalysis_data/ERA5/'
setupFile = '/home/users/pete_nut/IV_shiptracks/matt_traj_code/SETUP_traj_lowcloud_15min.CFG'
ASCFile  = '/home/users/pete_nut/IV_shiptracks/matt_traj_code/ASCDATA_JASMIN.CFG'

lonInit = '-120.'
latInit = '45.'
year=2016
month=6
day=12
hour=15
minute=30
jday = datetime.datetime(year,month,day,hour,minute)

year   = str( jday.year  ).zfill(4)
month  = str( jday.month ).zfill(2)
day    = str( jday.day ).zfill(2)
hour   = str( jday.hour ).zfill(2)
minute = str( jday.minute ).zfill(2)
height = '0.5' #midway through the PBL (for now)
trajHours = 'A12'
if trajHours.find('A') >= 0: tName = 'tdump_traj_ALLWARD_v0.txt'
if trajHours.find('A') < 0:
    if int(trajHours) > 0: tName = 'tdump_traj_FORWARD_v0.txt'
if trajHours.find('A') < 0:
    if int(trajHours) < 0: tName = 'tdump_traj_BACKWARD_v0.txt'

EXP = year+month+day
rootPath = outpath+EXP+'/'
run_trajectory(latInit,lonInit,year,month,day,hour,minute,height,trajHours,rootPath,path_hysplit_code,metpath,setupFile,ASCFile,tName)
trajectory = read_trajectory_file( rootPath + tName )
print('timestep, latitude, longitude, time')
for i in range(len(trajectory['lat'])): print(i,trajectory['timestep'][i],trajectory['lat'][i],trajectory['lon'][i],trajectory['jday'][i])
