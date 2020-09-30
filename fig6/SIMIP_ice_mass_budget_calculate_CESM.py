#--------------------------------------------------------------------
# Example script to demonstrate how to create ice mass budget ASCII files
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
datafile='NCAR_CESM2f09g17_BL99_001_ice.txt'
case = 'b.e21.B1850.f09_g17.CMIP6-piControl.001_bl99snow1'
case2 = 'b.e21.B1850.f09_g17.CMIP6-piControl.001_snow1'
path = '/glade/p/cesm/pcwg/dbailey/archive/'
fileGlob = path+case+'/ice/proc/tseries/month_1/*sidmasslat*.nc'
files = sorted(glob.glob(fileGlob))

# derive dates from time_bounds array

ds = xr.open_mfdataset(files)

times = ds.time_bounds.values

print(times)

leftbounds_yr = [x[0].timetuple()[0] for x in times]
leftbounds_mo = [x[0].timetuple()[1] for x in times]

#--------------------------------------------------------------------
# define sea ice budget variables
budgetVars =   [
                  'sidmassgrowthbot',
                  'sidmassgrowthwat',
                  'sidmassmelttop',
                  'sidmassmeltbot',
                  'sidmasslat',
                  'sidmasssi',
                  'sidmassevapsubl',
                  'sidmassdyn',
                  'total'   
               ]

#--------------------------------------------------------------------
# read in static fields: mask and grid-cell-areas
#
# mask
fh = Dataset('arctic_region_mask_gx1v7.nc')
mask = fh.variables['mask'][:,:]  
fh.close()
# areas
maskFile = '/glade/p/cesm/omwg/grids/gx1v7_grid.nc'
fh = Dataset(maskFile)
tarea = fh.variables['TAREA'][:,:]  
tlat = fh.variables['TLAT'][:,:]  
fh.close()

tarea = tarea*1.0e-4

#--------------------------------------------------------------------
# define output file and populate header
#
# define formats
title_format = "%1s %4s %6s %14s %12s %17s %17s %15s %15s %12s %12s %15s %12s %12s"
data_format = "%6i %6i %14.5e %12.5e %17.5e %17.5e %15.5e %15.5e %12.5e %12.5e %15.5e %12.5e %12.5e"
headers = ('#','Year','Month','Area (Km**2)', 'Mass (Kg)', 
                  'sidmassgrowthbot',    
                  'sidmassgrowthwat',
                  'sidmassmelttop',    
                  'sidmassmeltbot',    
                  'sidmasslat',    
                  'sidmasssi',    
                  'sidmassevapsubl',   
                  'sidmassdyn',   
                  'total')

# create header
data_fileh = open(datafile,'w')
data_fileh.write('# Contact: David Bailey dbailey@ucar.edu')
data_fileh.write("\n")
data_fileh.write('# Corresponding HIST file: n/a')
data_fileh.write("\n")
data_fileh.write('# Components of the Arctic sea ice mass budget (Kg s-1):')
data_fileh.write("\n")
data_fileh.write(title_format % headers)
data_fileh.write("\n")

#--------------------------------------------------------------------
# loop over monthly files, calculate budgets, mass and area
# then output to ascii file
# open netcdf files

files = sorted(glob.glob(path+case+'/ice/proc/tseries/month_1/*sidmassgrowthbot*.nc'))
fh1 = MFDataset(files)
files = sorted(glob.glob(path+case+'/ice/proc/tseries/month_1/*sidmassgrowthwat*.nc'))
fh2 = MFDataset(files)
files = sorted(glob.glob(path+case+'/ice/proc/tseries/month_1/*sidmassmelttop*.nc'))
fh3 = MFDataset(files)
files = sorted(glob.glob(path+case+'/ice/proc/tseries/month_1/*sidmassmeltbot*.nc'))
fh4 = MFDataset(files)
files = sorted(glob.glob(path+case+'/ice/proc/tseries/month_1/*sidmasslat*.nc'))
fh5 = MFDataset(files)
files = sorted(glob.glob(path+case+'/ice/proc/tseries/month_1/*sidmasssi*.nc'))
fh6 = MFDataset(files)
files = sorted(glob.glob(path+case+'/ice/proc/tseries/month_1/*sidmassevapsubl*.nc'))
fh7 = MFDataset(files)
files = sorted(glob.glob(path+case+'/ice/proc/tseries/month_1/*sidmassdyn*.nc'))
fh8 = MFDataset(files)
files = sorted(glob.glob(path+case+'/ice/proc/tseries/month_1/*aice.*.nc'))
fh9 = MFDataset(files)
files = sorted(glob.glob(path+case2+'/ice/proc/tseries/month_1/*aice.*.nc'))
fh11 = MFDataset(files)
files = sorted(glob.glob(path+case+'/ice/proc/tseries/month_1/*hi.*.nc'))
fh10 = MFDataset(files)

time = fh9.variables['time']
ntimes = len(time)

#
# reset/zero the budget for this month
thisBudget = np.ma.masked_all([len(budgetVars)],dtype=float)
#
#
# calculate budget components

rhoi = 917.
dt = 1800.

for n in range(0,ntimes):
   aice1 = fh9.variables['aice'][n,:,:]
   aice2 = fh11.variables['aice'][n,:,:]

#  mask = np.where((aice1 > 0.15) & (aice2 > 0.15) & (mask < 1.0e10),mask,0.0)

   thisVar1 = fh1.variables[budgetVars[0]][n,:,:]
   thisBudget[0] = np.sum(thisVar1*tarea*mask,dtype=float)
   thisVar2 = fh2.variables[budgetVars[1]][n,:,:]
   thisBudget[1] = np.sum(thisVar2*tarea*mask,dtype=float)
   thisVar3 = fh3.variables[budgetVars[2]][n,:,:]
   thisBudget[2] = np.sum(-thisVar3*tarea*mask,dtype=float)
   thisVar4 = fh4.variables[budgetVars[3]][n,:,:]
   thisBudget[3] = np.sum(-thisVar4*tarea*mask,dtype=float)
   thisVar5 = fh5.variables[budgetVars[4]][n,:,:]
   thisBudget[4] = np.sum(-thisVar5*tarea*mask,dtype=float)
   thisVar6 = fh6.variables[budgetVars[5]][n,:,:]
   thisBudget[5] = np.sum(thisVar6*tarea*mask,dtype=float)
   thisVar7 = fh7.variables[budgetVars[6]][n,:,:]
   thisBudget[6] = np.sum(thisVar7*tarea*mask,dtype=float)/rhoi/dt
   thisVar8 = fh8.variables[budgetVars[7]][n,:,:]
   thisBudget[7] = np.sum(thisVar8*tarea*mask,dtype=float)
   #
   # sum all for total
   thisBudget[-1] = np.sum(thisBudget[0:-1])
   #
   # calculate mass and area
   thisVar = fh10.variables['hi'][n,:,:] 
   total_mass = np.sum(thisVar*tarea*mask,dtype=float)*rhoi
   thisVar = fh9.variables['aice'][n,:,:]
   total_area = np.sum(thisVar*tarea*mask,dtype=float)
   total_area = total_area / 1e6    # convert from m^2 to km^2 
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
                    thisBudget[6],
                    thisBudget[7],
                    thisBudget[8]))
   data_fileh.write("\n")

# close output file
data_fileh.close()

