#--------------------------------------------------------------------
# Example script to demonstrate how to create snow mass budget ASCII files
# from standard monthly-mean SIMIP diaganostics
# including applying the mask, area weigting, and summation
#
# Ed Blockley, Feb 2019
#--------------------------------------------------------------------

import numpy as np
import glob
from netCDF4 import Dataset
from netCDF4 import MFDataset
import datetime as dt
import xarray as xr
import cftime
import os

#--------------------------------------------------------------------
# define and count input files
# NB specific to our file naming
datafile='NCAR_CESM2f09g17_MUSHY_001_snow.txt'
case = 'b.e21.B1850.f09_g17.CMIP6-piControl.001_snow1'
case2 = 'b.e21.B1850.f09_g17.CMIP6-piControl.001_bl99snow1'
path = '/gpfs/fs1/p/cesm/pcwg/dbailey/archive/'
fileGlob = path+case+'/ice/proc/tseries/month_1/*sidmasslat*.nc'
files = sorted(glob.glob(fileGlob))

# derive dates from time_bounds array

ds = xr.open_mfdataset(files)

times = ds.time_bounds.values

leftbounds_yr = [x[0].timetuple()[0] for x in times]
leftbounds_mo = [x[0].timetuple()[1] for x in times]

print(leftbounds_yr)

#--------------------------------------------------------------------
# define sea ice budget variables
budgetVars =   [
                  'sndmasssnf',    
                  'sndmassmelt', 
                  'sndmasssi',    
                  'sndmasswindrif',
                  'sndmassubl',
                  'sndmassdyn',    
                  'total'
               ]

#--------------------------------------------------------------------
# read in static fields: mask and grid-cell-areas
#
# mask
fh = Dataset('arctic_region_mask_gx1v7.nc', mode='r')
mask = fh.variables['mask'][:,:]  
fh.close()
# areas
maskFile = '/glade/p/cesm/omwg/grids/gx1v7_grid.nc'
fh = Dataset(maskFile, mode='r')
tarea = fh.variables['TAREA'][:,:]  
tlat = fh.variables['TLAT'][:,:]  
fh.close()

tarea = tarea*1.0e-4

#--------------------------------------------------------------------
# define output file and populate header
#
# define formats
title_format = "%1s %4s %6s %14s %12s %14s %14s %14s %16s %14s %14s %12s"
data_format = "%6i %6i %14.5e %12.5e %14.5e %14.5e %14.5e %16.5e %14.5e %14.5e %12.5e"
headers = ('#','Year','Month','Area (Km**2)','Mass (Kg)',
                  'sndmasssnf',    
                  'sndmassmelt', 
                  'sndmasssi',    
                  'sndmasswindrif',
                  'sndmasssubl',
                  'sndmassdyn',   
                  'total')

# create header
data_fileh = open(datafile,'w')
data_fileh.write('# Contact: David Bailey dbailey@ucar.edu')
data_fileh.write("\n")
data_fileh.write('# Corresponding HIST file: n/a')
data_fileh.write("\n")
data_fileh.write('# Components of the Arctic snow on sea ice mass budget (Kg s-1):')
data_fileh.write("\n")
data_fileh.write(title_format % headers)
data_fileh.write("\n")

#--------------------------------------------------------------------
# loop over monthly files, calculate budgets, mass and area
# then output to ascii file

files = sorted(glob.glob(path+case+'/ice/proc/tseries/month_1/*sndmasssnf*.nc'))
fh1 = MFDataset(files)
files = sorted(glob.glob(path+case+'/ice/proc/tseries/month_1/*sndmassmelt*.nc'))
fh2 = MFDataset(files)
files = sorted(glob.glob(path+case+'/ice/proc/tseries/month_1/*sidmasssi*.nc'))
fh3 = MFDataset(files)
files = sorted(glob.glob(path+case+'/ice/proc/tseries/month_1/*sndmassubl*.nc'))
fh4 = MFDataset(files)
files = sorted(glob.glob(path+case+'/ice/proc/tseries/month_1/*aice.*.nc'))
fh5 = MFDataset(files)
files = sorted(glob.glob(path+case2+'/ice/proc/tseries/month_1/*aice.*.nc'))
fh8 = MFDataset(files)
files = sorted(glob.glob(path+case+'/ice/proc/tseries/month_1/*hs.*.nc'))
fh6 = MFDataset(files)

files = sorted(glob.glob(path+case+'/ice/proc/tseries/day_1/*hs_d.*.nc'))
fh7 = MFDataset(files)

days = [0.,31.,28.,31.,30.,31.,30.,31.,31.,30.,31.,30.,31.]
juldays = [-1,30,58,89,119,150,180,211,242,272,303,333,364]

time = fh8.variables['time']
ntimes = len(time)

# reset/zero the budget for this month
thisBudget = np.ma.masked_all([len(budgetVars)],dtype=float)

#
# calculate budget components

rhoi = 917.
rhos = 330.
dt = 1800.

for n in range(0,ntimes):


   aice1 = fh5.variables['aice'][n,:,:]
   aice2 = fh8.variables['aice'][n,:,:]

   mask = np.where((aice1 > 0.15) & (aice2 > 0.15) & (mask < 1.0e10),mask,0.0)

   thisVar1 = fh1.variables[budgetVars[0]][n,:,:]
   thisBudget[0] = np.sum(thisVar1*tarea*mask/rhos,dtype=float)
   thisVar2 = fh2.variables[budgetVars[1]][n,:,:]
   thisBudget[1] = np.sum(-thisVar2*tarea*mask,dtype=float)
   # convert snow-ice formation
   thisVar3 = fh3.variables['sidmasssi'][n,:,:]
   thisBudget[2] = np.sum(-thisVar3*tarea*mask*rhos/rhoi,dtype=float)
   # we don't have wind drifting
   thisBudget[3] = 0.
   # Fix sublimation
   thisVar5 = fh4.variables[budgetVars[4]][n,:,:]
   thisBudget[4] = np.sum(thisVar5*tarea*mask/rhos/dt,dtype=float)

   #
   #
   # calculate mass and area
   thisVar = fh6.variables['hs'][n,:,:] 
   total_mass = np.sum(thisVar*tarea*mask*rhos,dtype=float)

   #Estimate dsnowmassdt

   first_day = juldays[leftbounds_mo[n]-1] + 1 + (leftbounds_yr[n]-1)*365
   last_day = juldays[leftbounds_mo[n]] + (leftbounds_yr[n]-1)*365
#  first_day = juldays[leftbounds_mo[n]-1] + 1 + (leftbounds_yr[n]-1850)*365
#  last_day = juldays[leftbounds_mo[n]] + (leftbounds_yr[n]-1850)*365
#  first_day = juldays[leftbounds_mo[n]-1] + 1 + (leftbounds_yr[n]-2015)*365
#  last_day = juldays[leftbounds_mo[n]] + (leftbounds_yr[n]-2015)*365

   thisVar_first = fh7.variables['hs_d'][first_day,:,:] 
   mass_d_first = np.sum(thisVar_first*tarea*mask*rhos,dtype=float)

   thisVar_last = fh7.variables['hs_d'][last_day,:,:] 
   mass_d_last = np.sum(thisVar_last*tarea*mask*rhos,dtype=float)

   thisBudget[5] = (mass_d_last-mass_d_first)/((days[leftbounds_mo[n]]-1)*86400.)-np.sum(thisBudget[0:5])

   # sum all for total
   thisBudget[-1] = np.sum(thisBudget[0:-1])

   siconc = fh5.variables['aice'][n,:,:]
   total_area = np.sum(siconc*tarea*mask,dtype=float)
   total_area = total_area / 1e6    # convert from m^2 to km^2 
   #
   #
   # add row to ascii file
   data_fileh.write(data_format % (leftbounds_yr[n],
                    leftbounds_mo[n],
                    total_area,
                    total_mass,
                    thisBudget[0],
                    thisBudget[1],
                    thisBudget[2],
                    thisBudget[3],
                    thisBudget[4],
                    thisBudget[5],
                    thisBudget[6]))
   data_fileh.write("\n")

# close output file
data_fileh.close()

