from mpl_toolkits.basemap import Basemap
import netCDF4
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import pyresample
from geo_Collection import geo_web as gs
from shapely.geometry import Polygon, mapping
import rasterio as rio
from rasterio.mask import mask
from rasterio.plot import show
import earthpy as et
minlat = 32.5
minlon = 130.5
maxlat = 35
maxlon = 136
#area of full seto-inland sea
filename="I:/Nagoya University/Project/Seto/A2019152042500.L2_LAC_OC.nc"
nc_file = netCDF4.Dataset(filename, 'r')
lon=nc_file.groups['navigation_data'].variables['longitude'][:]
lat=nc_file.groups['navigation_data'].variables['latitude'][:]
lonindex = np.where((lon >= minlon) & (lon <= maxlon))
latindex = np.where((lat >= minlat) & (lat <= maxlat))
lon = np.array(lon[lonindex])
lat = np.array(lat[latindex])
minlati = np.min(latindex)
maxlati = np.max(latindex)
minloni = np.min(lonindex)
maxloni = np.max(lonindex)
chlor_a=nc_file.groups['geophysical_data'].variables['chlor_a'][minlati:maxlati, minloni:maxloni]
#crop out
#lons,lats=gs.creategrid(minlon,maxlon,minlat,maxlat,0.01)#1km resolution
#nlat,nlon,nchla=gs.geointerp(lats,lons,chlor_a,0.01)
x = np.arange(minlon, maxlon, 0.01) # 1 km grid, full resolution did not work but 500 m was OK!
y = np.arange(maxlat, minlat, -0.01)
#shape=np.zeros((lon.shape[0],lat.shape[0]),dtype='uint8')
#lon=lon.reshape(shape)
#lat=lat.reshape(shape)
result,grid=gs.swath_resampling(chlor_a,lon,lat,x,y,1000)
