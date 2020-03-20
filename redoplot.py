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
os.chdir('I:/Nagoya University/Project/Seto/GOCI')
path='I:/Nagoya University/Project/Seto/GOCI/'
datalist=glob.glob('*.nc')
# minlat = 32.5
# minlon = 130.5
# maxlat = 35
# maxlon = 136
#Full Seto-Inland Sea
# minlat = 32.65
# minlon = 131.39
# maxlat = 33.48
# maxlon = 132.77
#UWAJIMA
# minlat = 34.18
# minlon = 134.00
# maxlat = 34.88
# maxlon = 135.59
#OSAKA
#Specific for 20150715
minlat = 33.07428
minlon = 131.927472
maxlat = 33.473592
maxlon = 132.625198
L=len(datalist)
#area of full seto-inland sea
#for j in range(len(datalist)):
for a1 in range(L):
    file=datalist[a1]
    filename = path + file
    #print(filename)
    nc_file = nc4.Dataset(filename, 'r')
    print(nc_file.time_coverage_end)
    lon=nc_file.groups['navigation_data'].variables['longitude'][:]
    lat=nc_file.groups['navigation_data'].variables['latitude'][:]
    variables=nc_file.groups['geophysical_data'].variables

    x = np.arange(minlon, maxlon, 0.01) # 1 km grid,
    y = np.arange(maxlat, minlat, -0.01)
    for i in variables:
        var=variables[i][:]
        np.where(var<=0,var,np.nan)
        if i!='l2_flags':
            var_re, grid = gs.swath_resampling(var, lon, lat, x, y, 1500)  # 1 km grid
            #var_re=var_re.filled()
            variables[i] = var_re
        else:
            #var_re = var_re.filled()
            variables[i]=var

    lons=grid.lons
    lats=grid.lats
    Rrs_412=variables['Rrs_412']
    Rrs_443 = variables['Rrs_443']
    #Rrs_469 = variables['Rrs_469']
    Rrs_490 = variables['Rrs_490']
    #Rrs_531 = variables['Rrs_531']
    #Rrs_547 = variables['Rrs_547']
    Rrs_555 = variables['Rrs_555']
    #Rrs_645 = variables['Rrs_645']
    Rrs_660 = variables['Rrs_660']
    Rrs_680= variables['Rrs_680']
    #nflh=variables['nflh']
    chl=variables['chlor_a']
    gs.plot_geo_image(Rrs_412, lons, lats, log10=False, title='Rrs412' + nc_file.time_coverage_end[0:13],
                      save_image='Rrs412' + nc_file.time_coverage_end[0:13])
    gs.plot_geo_image(Rrs_443, lons, lats, log10=False, title='Rrs443' + nc_file.time_coverage_end[0:13],
                      save_image='Rrs443' + nc_file.time_coverage_end[0:13])
    gs.plot_geo_image(Rrs_490, lons, lats, log10=False, title='Rrs490' + nc_file.time_coverage_end[0:13],
                      save_image='Rrs490' + nc_file.time_coverage_end[0:13])
    gs.plot_geo_image(Rrs_555, lons, lats, log10=False, title='Rrs555' + nc_file.time_coverage_end[0:13],
                      save_image='Rrs555' + nc_file.time_coverage_end[0:13])
    gs.plot_geo_image(Rrs_660, lons, lats, log10=False, title='Rrs660' + nc_file.time_coverage_end[0:13],
                      save_image='Rrs660' + nc_file.time_coverage_end[0:13])
    gs.plot_geo_image(Rrs_680, lons, lats, log10=False, title='Rrs680' + nc_file.time_coverage_end[0:13],
                      save_image='Rrs680' + nc_file.time_coverage_end[0:13])
    # gs.plot_geo_image(Rrs_555, lons, lats, log10=False, title='Rrs555' + nc_file.time_coverage_end[0:13],
    #                   save_image='Rrs555' + nc_file.time_coverage_end[0:13])
    # gs.plot_geo_image(Rrs_645, lons, lats, log10=False, title='Rrs645' + nc_file.time_coverage_end[0:13],
    #                   save_image='Rrs645' + nc_file.time_coverage_end[0:13])
    # gs.plot_geo_image(Rrs_667, lons, lats, log10=False, title='Rrs667' + nc_file.time_coverage_end[0:13],
    #                   save_image='Rrs667' + nc_file.time_coverage_end[0:13])
    # gs.plot_geo_image(Rrs_678, lons, lats, log10=False, title='Rrs678' + nc_file.time_coverage_end[0:13],
    #                   save_image='Rrs678' + nc_file.time_coverage_end[0:13])
    print(print('percent: {:.2%}'.format((a1 + 1) / L)))