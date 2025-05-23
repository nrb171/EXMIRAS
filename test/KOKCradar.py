# KOKC case study
# this script reads in the RADAR data from KTLX and extracts the reflectivity and ZDR profiles
# and writes them to a netcdf file
import cartopy.crs as ccrs
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import numpy as np
np.round_ = np.round

from metpy.calc import azimuth_range_to_lat_lon
from metpy.cbook import get_test_data
from metpy.io import Level2File
from metpy.plots import add_metpy_logo, add_timestamp, USCOUNTIES
from metpy.units import units
from metpy.plots import ctables
import pyart

import os
import sys
import datetime
import netCDF4 as nc

#### USER INPUTS ###########################################################

query_lat = 35.3892  #  latitude of target
query_lon = -97.6005  #  longitude of target
dx, dy = -2.9275e+04, 6.2106e+03  # dx/y in meters from radar to target

############################################################################
# list all files in directory
directory = '/h/eol/nbarron/export/radardata/'
# os.chdir(directory)
files = os.listdir(directory)
files = [file for file in files if file.endswith('.gz')]
files.sort()

dts = np.zeros(len(files),dtype=datetime.datetime)
dts2 = np.zeros(len(files),dtype=datetime.datetime)
refProfile = np.zeros([len(files),18])
zdrProfile = np.zeros([len(files),18])
heightProfile = np.zeros([len(files),18])
for file,ifile in zip(files,range(0,len(files))):
        

    ###########################################
    # Open the file

    radar = pyart.io.read(directory+files[ifile])

    dts[ifile] = (datetime.datetime.strptime(files[ifile], 'KTLX%Y%m%d_%H%M%S_V06.gz') - datetime.datetime(1970,1,1)).total_seconds()
    dts2[ifile] = datetime.datetime.strptime(files[ifile], 'KTLX%Y%m%d_%H%M%S_V06.gz')

    

    for ii in range(0,18):
        try:
            sweep = ii
            x,y,z=radar.get_gate_x_y_z(sweep)
            lat,lon,a = radar.get_gate_lat_lon_alt(sweep)
            radarAtSweep = radar.extract_sweeps([ii])

            diff = np.sqrt((x - (dx))**2 + (y - (6.2106e+03))**2)
            idx = np.where(diff == np.min(diff))
            idx

            refProfile[ifile,ii] = 10*np.log10(10**(radarAtSweep.fields['reflectivity']['data'][
                slice(int(idx[0][0]-1), int(idx[0][0]+1)), 
                slice(int(idx[1][0]-1), int(idx[1][0]+1))]/10).mean())
            zdrProfile[ifile,ii] = 10*np.log10(10**(radarAtSweep.fields['differential_reflectivity']['data'][
                slice(int(idx[0][0]-1), int(idx[0][0]+1)), 
                slice(int(idx[1][0]-1), int(idx[1][0]+1))]/10).mean())
            heightProfile[ifile,ii] = z[idx]
        except IndexError:
            refProfile[ifile,ii] = np.nan
            zdrProfile[ifile,ii] = np.nan
            heightProfile[ifile,ii] = np.nan

    print('Finished with file '+files[ifile])

## some sweeps are repeated, average those
uHeights = np.unique(heightProfile[0,:])
uHeights.sort()

refProfile2 = np.zeros([len(files),len(uHeights)])
zdrProfile2 = np.zeros([len(files),len(uHeights)])
for ii in range(0,len(uHeights)):
    mask = heightProfile[0,:] == uHeights[ii]
    refProfile2[:,ii] = 10*np.log10(np.nanmean(10**(refProfile[:,mask]/10),axis=1))
    zdrProfile2[:,ii] = 10*np.log10(np.nanmean(10**(zdrProfile[:,mask]/10),axis=1))

# write data to netcdf file
ncfile = nc.Dataset('KTLXProfiles.nc','w')
ncfile.createDimension('time',len(files))
ncfile.createDimension('height',len(uHeights))

ncfile.createVariable('time','f8',('time',))
ncfile.createVariable('height','f8',('height',))
ncfile.createVariable('refProfile','f8',('time','height',))
ncfile.createVariable('zdrProfile','f8',('time','height',))

ncfile.variables['time'][:] = dts
ncfile.variables['height'][:] = uHeights
ncfile.variables['refProfile'][:] = refProfile2
ncfile.variables['zdrProfile'][:] = zdrProfile2

ncfile.close()
