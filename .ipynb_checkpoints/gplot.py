import netCDF4 as nc4
from mpl_toolkits.basemap import Basemap
import numpy as np

import matplotlib.pyplot as plt

from matplotlib.colors import LogNorm

#data need to be reprojected by seadas using nearst meathod to WGS84
#先放着，去写QAA
nc1=nc4.Dataset('I:/Nagoya University/Project/HAB/satellite/GOCI/G2019246051643.L2_COMS_OC_reprojected.nc', 'r')

lon=np.array(nc1.variables['lon'])
lat=np.array(nc1.variables['lat'])

latmax=32.7
latmin=32
lonmax=130.7
lonmin=130

lonindex=np.where((lon>=lonmin)&(lon<=lonmax))
latindex=np.where((lat>=latmin)&(lat<=latmax))
lon=np.array(lon[lonindex])
lat=np.array(lat[latindex])

minlati=np.min(latindex)
maxlati=np.max(latindex)
minloni=np.min(lonindex)
maxloni=np.max(lonindex)
# Rrs412=np.array(nc1.variables['Rrs_412'])[minlati:maxlati,minloni:maxloni]
# Rrs443=np.array(nc1.variables['Rrs_443'])[minlati:maxlati,minloni:maxloni]
# Rrs469=np.array(nc1.variables['Rrs_469'])[minlati:maxlati,minloni:maxloni]
# Rrs488=np.array(nc1.variables['Rrs_488'])[minlati:maxlati,minloni:maxloni]
# Rrs531=np.array(nc1.variables['Rrs_531'])[minlati:maxlati,minloni:maxloni]
# Rrs547=np.array(nc1.variables['Rrs_547'])[minlati:maxlati,minloni:maxloni]
# Rrs555=np.array(nc1.variables['Rrs_555'])[minlati:maxlati,minloni:maxloni]
# Rrs645=np.array(nc1.variables['Rrs_645'])[minlati:maxlati,minloni:maxloni]
# Rrs667=np.array(nc1.variables['Rrs_667'])[minlati:maxlati,minloni:maxloni]
# Rrs678=np.array(nc1.variables['Rrs_678'])[minlati:maxlati,minloni:maxloni]
chl_a=np.array(nc1.variables['chlor_a'])[minlati:maxlati,minloni:maxloni]
# #


fig = plt.figure()

ax = fig.add_subplot(1, 1, 1)

#setup Lambert Conformal basemap.
ax = Basemap(projection='merc',
llcrnrlat = latmin,
llcrnrlon = lonmin,
urcrnrlat = latmax,
urcrnrlon = lonmax,
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
# ax.drawrivers(color='#0000ff')
parallels=np.arange(latmin,latmax+0.01,0.2)
ax.drawparallels(parallels,fontsize=10,linewidth=0.25,dashes=[7,15],
color='k',labels=[1,0,1,1])
##dash：http://bbs.06climate.com/home.php?mod=space&uid=8934&do=blog&id=3284 参考
##linewidth：线的宽度
meridians = np.arange(lonmin,lonmax+0.01,0.2)
ax.drawmeridians(meridians,fontsize=10,dashes=[7,15],
linewidth=0.3, color='k',labels=[1,1,0,1])




#
lons,lats=np.meshgrid(lon,lat)#把网格画出来QAQ，一直忘了这个
#
x,y=ax(lons,lats)
im = ax.pcolor(x,y,chl_a, cmap='jet',norm=LogNorm(vmin=0.1, vmax=70))
#ax.tick_params(labelsize=10)
#Draw colorbar
#WIM里面用的standard inverse究竟是啥
cbar =ax.colorbar(im, location='right',pad="10%")
# Set colorbar label
unit = 'Chlorophyll-a (mg m-3)'
cbar.set_label(unit, rotation=270, labelpad=10.0, fontsize=10)
cbar.ax.tick_params(labelsize=10)
plt.show()
plt.savefig('G2019246051643.L2_COMS_OC_reprojectedchla.png')