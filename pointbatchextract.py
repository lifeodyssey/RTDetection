import netCDF4 as nc4
from mpl_toolkits.basemap import Basemap
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import os
import glob
from QAAV6GOCI import QAAv6
os.chdir('I:/Nagoya University/Project/HAB/satellite/GOCI/')
path='I:/Nagoya University/Project/HAB/satellite/GOCI/'
datalist=glob.glob('*248*_reprojected.nc')
for i in range(len(datalist)):
    file=datalist[i]


    filename=path+file
    print(filename)
    nc1=nc4.Dataset(filename,'r')
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
    Rrs412 = np.array(nc1.variables['Rrs_412'])[minlati:maxlati, minloni:maxloni]
    Rrs443 = np.array(nc1.variables['Rrs_443'])[minlati:maxlati, minloni:maxloni]
    Rrs490 = np.array(nc1.variables['Rrs_490'])[minlati:maxlati, minloni:maxloni]
    Rrs555 = np.array(nc1.variables['Rrs_555'])[minlati:maxlati, minloni:maxloni]
    Rrs660 = np.array(nc1.variables['Rrs_660'])[minlati:maxlati, minloni:maxloni]
    Rrs680 = np.array(nc1.variables['Rrs_680'])[minlati:maxlati, minloni:maxloni]

    chl_a = np.array(nc1.variables['chlor_a'])[minlati:maxlati, minloni:maxloni]
    chla10index = np.where((chl_a >= 10))
    chla10 = chl_a[chla10index[0], chla10index[1]]
    lat10 = lat[chla10index[0]]
    lon10 = lon[chla10index[1]]
    Rrs412_10 = Rrs412[chla10index[0], chla10index[1]]
    Rrs443_10 = Rrs443[chla10index[0], chla10index[1]]
    Rrs490_10 = Rrs490[chla10index[0], chla10index[1]]
    Rrs555_10 = Rrs555[chla10index[0], chla10index[1]]
    Rrs660_10 = Rrs660[chla10index[0], chla10index[1]]
    Rrs680_10 = Rrs680[chla10index[0], chla10index[1]]

    Rrsindex = np.where(
        (Rrs412_10 > 0) & (Rrs443_10 > 0) & (Rrs490_10 > 0) & (Rrs555_10 > 0) & (Rrs660_10 > 0) & (Rrs680_10 > 0))
    chla10 = chla10[Rrsindex]
    lat10 = lat10[Rrsindex]
    lon10 = lon10[Rrsindex]
    Rrs412_10 = Rrs412_10[Rrsindex]
    Rrs443_10 = Rrs443_10[Rrsindex]
    Rrs490_10 = Rrs490_10[Rrsindex]
    Rrs555_10 = Rrs555_10[Rrsindex]
    Rrs660_10 = Rrs660_10[Rrsindex]
    Rrs680_10 = Rrs680_10[Rrsindex]

    aph=np.ones([len(chla10),6])
    bbp=np.ones([len(chla10),6])



    for j in range(len(chla10)):
        [aph[j,:],bbp[j,:]]=QAAv6(np.array([Rrs412_10[j],Rrs443_10[j],Rrs490_10[j],Rrs555_10[j],Rrs660_10[j],Rrs680_10[j]]))


    result = np.vstack([np.array(lon10),
                        np.array(lat10),
                        np.array(chla10),
                        np.array(Rrs412_10),
                        np.array(Rrs443_10),
                        np.array(Rrs490_10),
                        np.array(Rrs555_10),
                        np.array(Rrs660_10),
                        np.array(Rrs680_10)])
    result=np.concatenate((result.T,aph,bbp),axis=1)

    np.savetxt(file+'.csv', result, delimiter=',')


