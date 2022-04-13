# invisible_tracks

Code to track ship trajectories from AIS data, advect them with ERA5 winds using HYSPLIT, and collocate them to MODIS satellite data.

## Data set production
Below are examples of which keywords the python scripts expect to run the analysis for the month of 01/2016. 
In order to run it, maps of emission locations are needed as well as a HYSPLIT installation, meteorological (e.g. ERA5) data, as well as MODIS cloud product data. The filepaths that need replacing on a different system are listed for each python script.  

`python 0_modis_list.py 1 2016`
files in: MODIS files: `/neodc/modis/data/MOD06_L2/collection61/`\
files out: MODIS that overlaps with Region of Interest\
filepath0: `/home/users/pete_nut/IV_shiptracks/modis_tiles/atlantic_modis_swaths_intersect_{}_{}.pkl`

`python 1_track_emissions.py 1 2016`

files in: Emissions
filepath1: `/gws/nopw/j04/eo_shared_data_vol2/scratch/pete_nut/Emissions/'+year+'/SOx_*.nc`\
files out: tracked Emissions\
filepath2: `/gws/nopw/j04/eo_shared_data_vol2/scratch/pete_nut/emissions_tracked/{}/{}_{}`

`python 2_advect_tracks_hysplit.py 1 2016`

files in: tracked Emissions
filepath2: `/gws/nopw/j04/eo_shared_data_vol2/scratch/pete_nut/emissions_tracked/{}/{}_{}`\
files out: hysplit-advected Emissions\
filepath3: `/gws/nopw/j04/eo_shared_data_vol1/satellite/.modistracks/hysplit_adv_20m_non_isob/{}_{}_{}_{}:{}`\
this script also requires a hysplit installation and meteorology files, as well as a setup file and an ASCFile:

    outpath='/home/users/pete_nut/test/'
    os.system('mkdir -p '+outpath)
    path_hysplit_code = '/home/users/pete_nut/hysplit.v5.0.0_CentOS/'
    metpath = '/gws/nopw/j04/eo_shared_data_vol1/reanalysis/ARL_noaa/reanalysis_data/ERA5/'
    setupFile = '/home/users/pete_nut/IV_shiptracks/matt_traj_code/SETUP_traj_lowcloud_15min.CFG'
    ASCFile  = '/home/users/pete_nut/IV_shiptracks/matt_traj_code/ASCDATA_JASMIN.CFG'

`python 3_interpolate_collocate_cfaw_hy.py 1 Terra 2016 nonull`

files in:\
tracked Emissions\
hysplit-advected Emissions\
MODIS tiles that overlap with Region of Interest

files out: monthly daily collocated MODIS and emissions\
filepath4: `/gws/nopw/j04/eo_shared_data_vol1/satellite/.modistracks/{}/ctt_aqua_{}{}.h5`





