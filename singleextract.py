import netCDF4 as nc4
from mpl_toolkits.basemap import Basemap
import numpy as np
import pandas as pd

from matplotlib.colors import LogNorm
import os
import glob
import matplotlib.pyplot as plt
from QAAV6GOCI import QAAv6
nc1=nc4.Dataset('I:/Nagoya University/Project/HAB/satellite/GOCI/wanyine484.nc','r')
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
Rrs490=np.array(nc1.variables['Rrs_490'])[minlati:maxlati,minloni:maxloni]
Rrs555=np.array(nc1.variables['Rrs_555'])[minlati:maxlati,minloni:maxloni]
Rrs660=np.array(nc1.variables['Rrs_660'])[minlati:maxlati,minloni:maxloni]
Rrs680=np.array(nc1.variables['Rrs_680'])[minlati:maxlati,minloni:maxloni]

chl_a = np.array(nc1.variables['chlor_a'])[minlati:maxlati, minloni:maxloni]
chla10index=np.where((chl_a>=10))
chla10=chl_a[chla10index[0],chla10index[1]]
lat10=lat[chla10index[0]]
lon10=lon[chla10index[1]]
Rrs412_10=Rrs412[chla10index[0],chla10index[1]]
Rrs443_10=Rrs443[chla10index[0],chla10index[1]]
Rrs490_10=Rrs490[chla10index[0],chla10index[1]]
Rrs555_10=Rrs555[chla10index[0],chla10index[1]]
Rrs660_10=Rrs660[chla10index[0],chla10index[1]]
Rrs680_10=Rrs680[chla10index[0],chla10index[1]]
#还要再去掉Rrs是负值的地方
Rrsindex=np.where((Rrs412_10>0)&(Rrs443_10>0)&(Rrs490_10>0)&(Rrs555_10>0)&(Rrs660_10>0)&(Rrs680_10>0))
chla10=chla10[Rrsindex]
lat10=lat10[Rrsindex]
lon10=lon10[Rrsindex]
Rrs412_10=Rrs412_10[Rrsindex]
Rrs443_10=Rrs443_10[Rrsindex]
Rrs490_10=Rrs490_10[Rrsindex]
Rrs555_10=Rrs555_10[Rrsindex]
Rrs660_10=Rrs660_10[Rrsindex]
Rrs680_10=Rrs680_10[Rrsindex]
aph=np.ones([len(chla10),6])
bbp=np.ones([len(chla10),6])
a=np.ones([len(chla10),6])
adg=np.ones([len(chla10),6])
rrs=np.ones([len(chla10),6])
u=np.ones([len(chla10),6])
for j in range(len(chla10)):
    [rrs[j,:],u[j,:],a[j,:],adg[j,:],aph[j,:],bbp[j,:]]=np.abs(QAAv6(np.array([Rrs412_10[j],Rrs443_10[j],Rrs490_10[j],Rrs555_10[j],Rrs660_10[j],Rrs680_10[j]])))


result = np.vstack([np.array(lon10),
                        np.array(lat10),
                        np.array(chla10),
                        np.array(Rrs412_10),
                        np.array(Rrs443_10),
                        np.array(Rrs490_10),
                        np.array(Rrs555_10),
                        np.array(Rrs660_10),
                        np.array(Rrs680_10)])
a_mean=[np.mean(a[:,0]),np.mean(a[:,1]),np.mean(a[:,2]),np.mean(a[:,3]),np.mean(a[:,4]),np.mean(a[:,5])]
chla_mean=np.mean(chla10)
aph_mean=[np.mean(aph[:,0]),np.mean(aph[:,1]),np.mean(aph[:,2]),np.mean(aph[:,3]),np.mean(aph[:,4]),np.mean(aph[:,5])]
adg_mean=[np.mean(adg[:,0]),np.mean(adg[:,1]),np.mean(adg[:,2]),np.mean(adg[:,3]),np.mean(adg[:,4]),np.mean(adg[:,5])]
bbp_mean=[np.mean(bbp[:,0]),np.mean(bbp[:,1]),np.mean(bbp[:,2]),np.mean(bbp[:,3]),np.mean(bbp[:,4]),np.mean(bbp[:,5])]
Rrs_mean=[np.mean(Rrs412_10),np.mean(Rrs443_10),np.mean(Rrs490_10),np.mean(Rrs555_10),np.mean(Rrs660_10),np.mean(Rrs680_10)]
rrs_mean=[np.mean(rrs[:,0]),np.mean(rrs[:,1]),np.mean(rrs[:,2]),np.mean(rrs[:,3]),np.mean(rrs[:,4]),np.mean(rrs[:,5])]
u_mean=[np.mean(u[:,0]),np.mean(u[:,1]),np.mean(u[:,2]),np.mean(u[:,3]),np.mean(u[:,4]),np.mean(u[:,5])]
result=np.concatenate((result.T,aph,bbp,adg,a,rrs,u),axis=1)
#np.savetxt('G2019248011640.csv', result, delimiter=',')
wave=[412, 443, 490, 555, 660, 680]
plt.figure(1)
plt.subplot(4,2,1)
plt.plot(wave,a_mean/chla_mean,color="r", linestyle="-", marker="^", linewidth=1)
plt.xlabel("wavelength")
plt.ylabel("a")
plt.title("a")

plt.subplot(4,2,2)
plt.plot(wave,aph_mean/chla_mean,color="r", linestyle="-", marker="^", linewidth=1)
plt.xlabel("wavelength")
plt.ylabel("aph")
plt.title("aph")

plt.subplot(4,2,3)
plt.plot(wave,adg_mean/chla_mean,color="r", linestyle="-", marker="^", linewidth=1)
plt.xlabel("wavelength")
plt.ylabel("adg")
plt.title("adg")

plt.subplot(4,2,4)
plt.plot(wave,bbp_mean/chla_mean,color="r", linestyle="-", marker="^", linewidth=1)
plt.xlabel("wavelength")
plt.ylabel("bbp")
plt.title("bbp")

plt.subplot(4,2,5)
plt.plot(wave,Rrs_mean,color="r", linestyle="-", marker="^", linewidth=1)
plt.xlabel("wavelength")
plt.ylabel("Rrs")
plt.title("Rrs")

plt.subplot(4,2,6)
plt.plot(wave,rrs_mean,color="r", linestyle="-", marker="^", linewidth=1)
plt.xlabel("wavelength")
plt.ylabel("rrs")
plt.title("rrs")

plt.subplot(4,2,7)
plt.plot(wave,u_mean,color="r", linestyle="-", marker="^", linewidth=1)
plt.xlabel("wavelength")
plt.ylabel("u")
plt.title("u")

plt.show()

fig = plt.figure(2)



# setup Lambert Conformal basemap.
ax = Basemap(projection='merc',
             llcrnrlat=latmin,
             llcrnrlon=lonmin,
             urcrnrlat=latmax,
             urcrnrlon=lonmax,
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
parallels = np.arange(latmin, latmax + 0.01, 0.2)
ax.drawparallels(parallels, fontsize=10, linewidth=0.25, dashes=[7, 15],
                 color='k', labels=[1, 0, 1, 1])
##dash：http://bbs.06climate.com/home.php?mod=space&uid=8934&do=blog&id=3284 参考
##linewidth：线的宽度
meridians = np.arange(lonmin, lonmax + 0.01, 0.2)
ax.drawmeridians(meridians, fontsize=10, dashes=[7, 15],
                 linewidth=0.3, color='k', labels=[1, 1, 0, 1])

#
lons, lats = np.meshgrid(lon, lat)  # 把网格画出来QAQ，一直忘了这个
#
x, y = ax(lons, lats)
im = ax.pcolor(x, y, chl_a, cmap='jet', norm=LogNorm(vmin=0.1, vmax=70))
# ax.tick_params(labelsize=10)
# Draw colorbar
# WIM里面用的standard inverse究竟是啥
cbar = ax.colorbar(im, location='right', pad="10%")
# Set colorbar label
unit = 'Chlorophyll-a (mg m-3)'
cbar.set_label(unit, rotation=270, labelpad=10.0, fontsize=10)
cbar.ax.tick_params(labelsize=10)
plt.show()