
import netCDF4 as nc4
import numpy as np

from geo_Collection import geo_web as gs

import os
import glob
os.chdir('I:/Nagoya University/Project/Seto/requested_files')
path='I:/Nagoya University/Project/Seto/requested_files/'
datalist=glob.glob('*.nc')
minlat = 32.5
minlon = 130.5
maxlat = 35
maxlon = 136
#area of full seto-inland sea
for j in range(len(datalist)):
    file=datalist[j]
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
            var_re, grid = gs.swath_resampling(var, lon, lat, x, y, 1000)  # 1 km grid
            variables[i] = var_re
        else:
            variables[i]=var

    lons=grid.lons
    lats=grid.lats
    time=nc_file.time_coverage_end[0:13]
    #TODO can't batch? currently save as npz
    np.savez(path+time+'.npz',lons,lats,variables,time)
    p=j/len(datalist)
    nc_file.close()
    print('percent: {:.2%}'.format(p))
# open a netCDF file to write
# time=nc_file.date_created[0:10]
# nfname="I:/Nagoya University/Project/Seto/"+time+'.nc'
# nfname=nfname.replace('-','')
# ncout = nc4.Dataset(nfname, 'w', format='NETCDF4')
#
# # define axis size
# ncout.createDimension('time', None)  # unlimited
# ncout.createDimension('lat', len(y))
# ncout.createDimension('lon', len(x))
#
# # create time axis
# time = ncout.createVariable('time', np.dtype('double').char, ('time',))
# time.long_name = 'time'
# time.units = nc_file.date_created
# time.calendar = 'standard'
# time.axis = 'T'
#
# # create latitude axis
# lat = ncout.createVariable('lat', np.dtype('double').char, ('lat'))
# lat.standard_name = 'latitude'
# lat.long_name = 'latitude'
# lat.units = 'degrees_north'
# lat.axis = 'Y'
#
# # create longitude axis
# lon = ncout.createVariable('lon', np.dtype('double').char, ('lon'))
# lon.standard_name = 'longitude'
# lon.long_name = 'longitude'
# lon.units = 'degrees_east'
# lon.axis = 'X'
#
# # create variable array
# vout = ncout.createVariable('t2m', np.dtype('OrderDict').char, ('time', 'lat', 'lon'))
# vout.long_name = 'variables'
#
# # copy axis from original dataset
# time[:] = time
# lon[:] = x[:]
# lat[:] = y[:]
# vout= variables
#
# # close files
# nc_file.close()
# ncout.close()

#gs.plot_geo_image(chla,x,y,caxis=[0.1,70])
