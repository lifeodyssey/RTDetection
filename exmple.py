from mpl_toolkits.basemap import Basemap
import netCDF4
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.colors import LogNorm

# filename="I:/Nagoya University/Project/HAB/Data/A20191019Ariake.hdf"
# print(filename)
# nc_file = netCDF4.Dataset(filename, 'r')
# print(nc_file)
#nc_var = nc_file.variables.keys()
#print(nc_var)

# set latitude, Longitude of Image
# nc_lat  = nc_file.variables['Latitude']
# nc_lon  = nc_file.variables['Longitude']
# lat = np.array(nc_lat)
# lon = np.array(nc_lon)
# print(lat,lon)
#
# south = lat.min()
# north = lat.max()
# west = lon.min()
# east = lon.max()
# llat=lat.shape[0]
# llon=lat.shape[1]
#print (south, north, west, east,llat,llon)

# DNchl = nc_file.variables['chlor_a-Conv inv-remapped to \'MODIS\''][:]
# slope = nc_file.variables['chlor_a-Conv inv-remapped to \'MODIS\''].Slope
# intercept = nc_file.variables['chlor_a-Conv inv-remapped to \'MODIS\''].Intercept
#print(f'Slope: {slope}\nIntercept: {intercept}\nDNChl: {DNchl}')
#print(DNchl.max(), DNchl.min(), DNchl.mean())
#print('min=, max=, mean=', chl.min(),chl.max(),chl.mean())
#print('data shape=',chl.shape)
#print('ndim=',chl.ndim)
#print(DNchl.dtype.name)
#print(chl.size)
#print(chl.shape)

#nlw565 = nc_file.variables\
#['Image_data/NWLR_565 inv-remapped to \'Grid\' inv-remapped to \'Linear\'']
#nlw670 = nc_file.variables\
#['Image_data/NWLR_670 inv-remapped to \'Grid\' inv-remapped to \'Linear\'']

# DNchl=np.where((DNchl < 0), 255 + DNchl, DNchl)
# print('DNchl=', DNchl.max(), DNchl.min(), DNchl.mean())
# chl = 10**(DNchl*slope+intercept)
# print('chl=', chl)
# print(chl.max(),chl.min(), chl.mean())

# Draw chl image
# num_x=1
# num_y=1
#fig = plt.figure(figsize=[llon/100, llat/100], dpi=300)
fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)


# setup Lambert Conformal basemap.
ax = Basemap(projection='mill',
            llcrnrlat = 32.5,
            llcrnrlon = 130.8,
            urcrnrlat = 35,
            urcrnrlon = 136,
            resolution='f')
# draw coastlines.
ax.drawcoastlines()
# draw a boundary around the map, fill the background.
# this background will end up being the ocean color, since
# the continents will be drawn on top.
ax.drawmapboundary()
# fill continents, set lake color same as ocean color.
ax.fillcontinents()
ax.drawcountries()
ax.drawrivers(color='#0000ff')
parallels=np.arange(32.5,35,0.5)
ax.drawparallels(parallels,fontsize=10,linewidth=0.25,dashes=[7,15],
                color='k',labels=[1,0,1,1])
##dash：http://bbs.06climate.com/home.php?mod=space&uid=8934&do=blog&id=3284 参考
##linewidth：线的宽度
meridians = np.arange(131,136,0.5)
ax.drawmeridians(meridians,fontsize=10,dashes=[7,15],
                linewidth=0.3, color='k',labels=[1,1,0,1])
# x,y=ax(lon,lat)
#
# im = ax.pcolor(x,y,np.squeeze(chl), cmap='jet',norm=LogNorm(vmin=0.1, vmax=70))
#ax.tick_params(labelsize=10)
#Draw colorbar
#
# cbar =ax.colorbar(im, location='right',pad="10%")
# # Set colorbar label
# unit = 'Chl-a (mg m-3)'
# cbar.set_label(unit, rotation=270, labelpad=10.0, fontsize=10)
# cbar.ax.tick_params(labelsize=10)
#plt.show()