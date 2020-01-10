import netCDF4 as nc4
from mpl_toolkits.basemap import Basemap
import numpy as np
import pandas as pd

from matplotlib.colors import LogNorm
import os
import glob
import matplotlib.pyplot as plt
from QAAV6 import QAAv6
nc1=nc4.Dataset('I:/Nagoya University/Project/HAB/satellite/modis/A2019248042500.L2_LAC_OC.nc_subset_reprojected.nc','r')
lon = np.array(nc1.variables['lon'])
lat = np.array(nc1.variables['lat'])

latmax = 32.7
latmin = 32
lonmax = 130.7
lonmin = 130

lonindex = np.where((lon >= lonmin) & (lon <= lonmax))
latindex = np.where((lat >= latmin) & (lat <= latmax))
lon = np.array(lon[lonindex])
lat = np.array(lat[latindex])
minlati = np.min(latindex)
maxlati = np.max(latindex)
minloni = np.min(lonindex)
maxloni = np.max(lonindex)
Rrs412=np.array(nc1.variables['Rrs_412'])[minlati:maxlati,minloni:maxloni]
Rrs443=np.array(nc1.variables['Rrs_443'])[minlati:maxlati,minloni:maxloni]
Rrs469=np.array(nc1.variables['Rrs_469'])[minlati:maxlati,minloni:maxloni]
Rrs488=np.array(nc1.variables['Rrs_488'])[minlati:maxlati,minloni:maxloni]
Rrs531=np.array(nc1.variables['Rrs_531'])[minlati:maxlati,minloni:maxloni]
Rrs547=np.array(nc1.variables['Rrs_547'])[minlati:maxlati,minloni:maxloni]
Rrs555=np.array(nc1.variables['Rrs_555'])[minlati:maxlati,minloni:maxloni]
Rrs645=np.array(nc1.variables['Rrs_645'])[minlati:maxlati,minloni:maxloni]
Rrs667=np.array(nc1.variables['Rrs_667'])[minlati:maxlati,minloni:maxloni]
Rrs678=np.array(nc1.variables['Rrs_678'])[minlati:maxlati,minloni:maxloni]

chl_a = np.array(nc1.variables['chlor_a'])[minlati:maxlati, minloni:maxloni]
# chla10index=np.where((chl_a>=10))
# chla10=chl_a[chla10index[0],chla10index[1]]
# lat10=lat[chla10index[0]]
# lon10=lon[chla10index[1]]
# Rrs412_10=Rrs412[chla10index[0],chla10index[1]]
# Rrs443_10=Rrs443[chla10index[0],chla10index[1]]
# Rrs469_10=Rrs443[chla10index[0],chla10index[1]]
# Rrs488_10=Rrs488[chla10index[0],chla10index[1]]
# Rrs531_10=Rrs531[chla10index[0],chla10index[1]]
# Rrs547_10=Rrs547[chla10index[0],chla10index[1]]
# Rrs555_10=Rrs555[chla10index[0],chla10index[1]]
# Rrs645_10=Rrs645[chla10index[0],chla10index[1]]
# Rrs667_10=Rrs667[chla10index[0],chla10index[1]]
# Rrs678_10=Rrs678[chla10index[0],chla10index[1]]
#还要再去掉Rrs是负值的地方
Rrsindex=np.where((Rrs412>0)&(Rrs443>0)&(Rrs469>0)&(Rrs488>0)&(Rrs531>0)&(Rrs547>0)&(Rrs555>0)&(Rrs645>0)&(Rrs667>0)&(Rrs678>0))
chla10=chl_a[Rrsindex[0],Rrsindex[1]]
lat10=lat[Rrsindex[0]]
lon10=lon[Rrsindex[1]]
Rrs412_10=Rrs412[Rrsindex[0],Rrsindex[1]]
Rrs443_10=Rrs443[Rrsindex[0],Rrsindex[1]]
Rrs469_10=Rrs469[Rrsindex[0],Rrsindex[1]]
Rrs488_10=Rrs488[Rrsindex[0],Rrsindex[1]]
Rrs531_10=Rrs531[Rrsindex[0],Rrsindex[1]]
Rrs547_10=Rrs547[Rrsindex[0],Rrsindex[1]]
Rrs555_10=Rrs555[Rrsindex[0],Rrsindex[1]]
Rrs645_10=Rrs645[Rrsindex[0],Rrsindex[1]]
Rrs667_10=Rrs667[Rrsindex[0],Rrsindex[1]]
Rrs678_10=Rrs678[Rrsindex[0],Rrsindex[1]]
aph=np.ones([len(chla10),10])
bbp=np.ones([len(chla10),10])
a=np.ones([len(chla10),10])
adg=np.ones([len(chla10),10])
rrs=np.ones([len(chla10),10])
u=np.ones([len(chla10),10])
for j in range(len(chla10)):
    [rrs[j,:],u[j,:],a[j,:],adg[j,:],aph[j,:],bbp[j,:]]=QAAv6(np.array([Rrs412_10[j],
                                                                        Rrs443_10[j],
                                                                        Rrs469_10[j],
                                                                        Rrs488_10[j],
                                                                        Rrs531_10[j],
                                                                        Rrs547_10[j],
                                                                        Rrs555_10[j],
                                                                        Rrs645_10[j],
                                                                        Rrs667_10[j],
                                                                        Rrs678_10[j]]))


Rrs= np.vstack([np.array([Rrs412_10,
                        Rrs443_10,
                        Rrs469_10,
                        Rrs488_10,
                        Rrs531_10,
                        Rrs547_10,
                        Rrs555_10,
                        Rrs645_10,
                        Rrs667_10,
                        Rrs678_10])])

wave=[412, 443,469, 488,531,547, 555, 645,667, 678]
from errobarplot import errorplot
plt.figure(1)
plt.subplot(4,2,1)

errorplot(wave,a,'a')

plt.subplot(4,2,2)

errorplot(wave,aph,'aph')

plt.subplot(4,2,3)

errorplot(wave,adg,'adg')

plt.subplot(4,2,4)

errorplot(wave,bbp,'bbp')

plt.subplot(4,2,5)

errorplot(wave,Rrs,'Rrs')

plt.subplot(4,2,6)

errorplot(wave,rrs,'rrs')

plt.subplot(4,2,7)

errorplot(wave,u,'u')

plt.show()
#
# fig = plt.figure(2)
#
#
#
# # setup Lambert Conformal basemap.
# ax = Basemap(projection='merc',
#              llcrnrlat=latmin,
#              llcrnrlon=lonmin,
#              urcrnrlat=latmax,
#              urcrnrlon=lonmax,
#              resolution='f')
# # draw coastlines.
# ax.drawcoastlines()
# # draw a boundary around the map, fill the background.
# # this background will end up being the ocean color, since
# # the continents will be drawn on top.
# ax.drawmapboundary()
# # fill continents, set lake color same as ocean color.
# ax.fillcontinents()
# ax.drawcountries()
# # ax.drawrivers(color='#0000ff')
# parallels = np.arange(latmin, latmax + 0.01, 0.2)
# ax.drawparallels(parallels, fontsize=10, linewidth=0.25, dashes=[7, 15],
#                  color='k', labels=[1, 0, 1, 1])
# ##dash：http://bbs.06climate.com/home.php?mod=space&uid=8934&do=blog&id=3284 参考
# ##linewidth：线的宽度
# meridians = np.arange(lonmin, lonmax + 0.01, 0.2)
# ax.drawmeridians(meridians, fontsize=10, dashes=[7, 15],
#                  linewidth=0.3, color='k', labels=[1, 1, 0, 1])
#
# #
# lons, lats = np.meshgrid(lon, lat)  # 把网格画出来QAQ，一直忘了这个
# #
# x, y = ax(lons, lats)
# im = ax.pcolor(x, y, chl_a, cmap='jet', norm=LogNorm(vmin=0.1, vmax=70))
# # ax.tick_params(labelsize=10)
# # Draw colorbar
# # WIM里面用的standard inverse究竟是啥
# cbar = ax.colorbar(im, location='right', pad="10%")
# # Set colorbar label
# unit = 'Chlorophyll-a (mg m-3)'
# cbar.set_label(unit, rotation=270, labelpad=10.0, fontsize=10)
# cbar.ax.tick_params(labelsize=10)
# plt.show()