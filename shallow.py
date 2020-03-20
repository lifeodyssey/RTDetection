import netCDF4 as nc4
import numpy as np
from collections  import OrderedDict
from geo_Collection import geo_web as gs
from QAAV6 import  QAAv6
import os
import glob
import warnings
warnings.filterwarnings("ignore")#用来去掉Warning,
import pandas as pd
import matplotlib.pyplot as plt
os.chdir('I:/Nagoya University/Project/Seto/MODIS')
path='I:/Nagoya University/Project/Seto/MODIS/'
datalist=glob.glob('20180723*.nc')
f="F:/GEBCO_2019_10_Mar_2020_543a6365a94f/GEBCO_2019_10_Mar_2020_543a6365a94f/gebco.nc"
# print(filename
ele = nc4.Dataset(f, 'r')
minlat = 32.5
minlon = 130.5
maxlat = 35
maxlon = 136

#Full Seto-Inland Sea
L=len(datalist)
for a1 in range(L):
    file=datalist[a1]
    filename = path + file
    #print(filename)
    nc_file = nc4.Dataset(filename, 'r')
    f= nc4.Dataset(f, 'r')
    lone = np.asarray(f.variables['lon'])
    late = np.asarray(f.variables['lat'])
    print(nc_file.time_coverage_end)
    lon=nc_file.groups['navigation_data'].variables['longitude'][:]
    lat=nc_file.groups['navigation_data'].variables['latitude'][:]
    variables=nc_file.groups['geophysical_data'].variables

    # x = np.arange(minlon, maxlon, 0.01) # 1 km grid,
    # y = np.arange(maxlat, minlat, -0.01)
    for i in variables:
        var=variables[i][:]
        np.where(var<=0,var,np.nan)
        if i!='l2_flags':
            var_re, grid = gs.swath_resampling(var, lon, lat, lone, late, 3000)  # 1 km grid
            #var_re=var_re.filled()
            variables[i] = var_re
        else:
            #var_re = var_re.filled()
            variables[i]=var

    lons=grid.lons
    lats=grid.lats
    kd=variables['Kd_490']

    ele = np.abs(np.asarray(f.variables['elevation']))
    # lonr=np.tile(lone,(562,1))
    # latr=np.tile(late,(1229,1)).T
    #
    # ere=gs.swath_resampling(ele, lonr, latr, x, y, 1000)
    r=ele*np.sqrt(kd)
    rr=r
    rr[np.where(r<=1)]=1
    rr[np.where(r > 1)] = 0
    gs.plot_geo_image(rr, lone, late,)
