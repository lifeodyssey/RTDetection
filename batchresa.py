
import netCDF4 as nc4
import numpy as np

from geo_Collection import geo_web as gs
import pandas as pd
from collections import OrderedDict
import json
import geopandas as gpd
import h5py
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
#for j in range(len(datalist)):
for j in range(1):
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
            #var_re=var_re.filled()
            variables[i] = var_re
        else:
            #var_re = var_re.filled()
            variables[i]=var

    lons=grid.lons
    lats=grid.lats
    #time=nc_file.time_coverage_end
    variables['lons']=lons
    variables['lats']=lats
    #variables['time']=time
    #I haven't solve that how to write a single string to netcdf

    # #JSON not support ndarray
    # #gdf = gpd.GeoDataFrame(variables)
    #df = pd.DataFrame([k.values() for k in variables], columns=variables.keys())

    #TODO can't batch? currently save as json


    #fname=path+time
    # with open(fname + '.json', 'a') as f:
    #     json.dump(variables, f, ensure_ascii=False)
    #     f.write('\n')
    #np.savez(path+time+'.npz',lons,lats,variables,time)
    #NPZ have problem of reading
    #hf=h5py.File(fname+'.h5','w')
    #hf.create_dataset('modis',data=variables)
    #hf.close()
    #p=j/len(datalist)
    #nc_file.close()
    #print('percent: {:.2%}'.format(p))
    #open a netCDF file to write
    #time=nc_file.date_created[0:10]
    T=nc_file.time_coverage_end[0:13]
    nfname="I:/Nagoya University/Project/Seto/"+str(T)+'.nc'
    nfname=nfname.replace('-','')
    ncout = nc4.Dataset(nfname, 'w', format='NETCDF4')

    # define axis size
    #ncout.createDimension('time', None)  # unlimited
    ncout.createDimension('lat', len(y))
    ncout.createDimension('lon', len(x))

    # # create time axis
    # time = ncout.createVariable('time', np.dtype('double').char, ('time',))
    # time.long_name = 'time'
    # time.units = nc_file.date_created
    # time.calendar = 'standard'
    # time.axis = 'T'

    # create latitude axis
    lat = ncout.createVariable('lat', np.dtype('float64').char, ('lat'))
    lat.standard_name = 'latitude'
    lat.long_name = 'latitude'
    lat.units = 'degrees_north'
    lat.axis = 'Y'

    # create longitude axis
    lon = ncout.createVariable('lon', np.dtype('float64').char, ('lon'))
    lon.standard_name = 'longitude'
    lon.long_name = 'longitude'
    lon.units = 'degrees_east'
    lon.axis = 'X'

    # create variable array
    group=ncout.createGroup('geovar')

    for key, value in variables.items():
        var=group.createVariable(key+'_re', np.dtype('float32').char, ('lon', 'lat'))
        var.long_name =key

        var=value


    group.close
    #存储还是有问题
    #先不管了

    # copy axis from original dataset
    #time= nc_file.time_coverage_end
    lon = lons
    lat= lats
    #vout= variables

    # close files
    nc_file.close()
    ncout.close()
    p=j/len(datalist)
    #nc_file.close()
    print('percent: {:.2%}'.format(p))

#gs.plot_geo_image(chla,x,y,caxis=[0.1,70])
