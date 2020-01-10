import netCDF4 as nc4
from mpl_toolkits.basemap import Basemap
import numpy as np
import pandas as pd

from matplotlib.colors import LogNorm
import os
import glob
import matplotlib.pyplot as plt
from QAAV6GOCI import QAAv6
##Non-blooming(NB)
ncnb=nc4.Dataset('I:/Nagoya University/Project/HAB/satellite/GOCI/non bloom/G2019193061642.L2_COMS_OCnc_reprojected.nc','r')
lonnb = np.array(ncnb.variables['lon'])
latnb = np.array(ncnb.variables['lat'])

latmax = 32.7
latmin = 32
lonmax = 130.7
lonmin = 130

lonindexnb = np.where((lonnb >= lonmin) & (lonnb <= lonmax))
latindexnb = np.where((latnb >= latmin) & (latnb <= latmax))
lon = np.array(lonnb[lonindexnb])
lat = np.array(latnb[latindexnb])
minlatinb = np.min(latindexnb)
maxlatinb = np.max(latindexnb)
minloninb = np.min(lonindexnb)
maxloninb = np.max(lonindexnb)
Rrs412nb=np.array(ncnb.variables['Rrs_412'])[minlatinb:maxlatinb,minloninb:maxloninb]
Rrs443nb=np.array(ncnb.variables['Rrs_443'])[minlatinb:maxlatinb,minloninb:maxloninb]
Rrs490nb=np.array(ncnb.variables['Rrs_490'])[minlatinb:maxlatinb,minloninb:maxloninb]
Rrs555nb=np.array(ncnb.variables['Rrs_555'])[minlatinb:maxlatinb,minloninb:maxloninb]
Rrs660nb=np.array(ncnb.variables['Rrs_660'])[minlatinb:maxlatinb,minloninb:maxloninb]
Rrs680nb=np.array(ncnb.variables['Rrs_680'])[minlatinb:maxlatinb,minloninb:maxloninb]

chl_anb = np.array(ncnb.variables['chlor_a'])[minlatinb:maxlatinb,minloninb:maxloninb]
# chla10indexnb=np.where((chl_anb>=10))
# chla10nb=chl_anb[chla10indexnb[0],chla10indexnb[1]]
# lat10nb=latnb[chla10indexnb[0]]
# lon10nb=lonnb[chla10indexnb[1]]
# Rrs412_10nb=Rrs412nb[chla10indexnb[0],chla10indexnb[1]]
# Rrs443_10nb=Rrs443nb[chla10indexnb[0],chla10indexnb[1]]
# Rrs490_10nb=Rrs490nb[chla10indexnb[0],chla10indexnb[1]]
# Rrs555_10nb=Rrs555nb[chla10indexnb[0],chla10indexnb[1]]
# Rrs660_10nb=Rrs660nb[chla10indexnb[0],chla10indexnb[1]]
# Rrs680_10nb=Rrs680nb[chla10indexnb[0],chla10indexnb[1]]
#还要再去掉Rrs是负值的地方
Rrsindexnb=np.where((Rrs412nb>0)&(Rrs443nb>0)&(Rrs490nb>0)&(Rrs555nb>0)&(Rrs660nb>0)&(Rrs680nb>0))
chla10nb=chl_anb[Rrsindexnb]
lat10nb=latnb[Rrsindexnb[0]]
lon10nb=lonnb[Rrsindexnb[1]]
Rrs412_10nb=Rrs412nb[Rrsindexnb[0],Rrsindexnb[1]]
Rrs443_10nb=Rrs443nb[Rrsindexnb[0],Rrsindexnb[1]]
Rrs490_10nb=Rrs490nb[Rrsindexnb[0],Rrsindexnb[1]]
Rrs555_10nb=Rrs555nb[Rrsindexnb[0],Rrsindexnb[1]]
Rrs660_10nb=Rrs660nb[Rrsindexnb[0],Rrsindexnb[1]]
Rrs680_10nb=Rrs680nb[Rrsindexnb[0],Rrsindexnb[1]]
aphnb=np.ones([len(chla10nb),6])
bbpnb=np.ones([len(chla10nb),6])
anb=np.ones([len(chla10nb),6])
adgnb=np.ones([len(chla10nb),6])
rrsnb=np.ones([len(chla10nb),6])
unb=np.ones([len(chla10nb),6])
for j in range(len(chla10nb)):
    [rrsnb[j,:],unb[j,:],anb[j,:],adgnb[j,:],aphnb[j,:],bbpnb[j,:]]=QAAv6(
        np.array(
            [Rrs412_10nb.T[j],
             Rrs443_10nb.T[j],
             Rrs490_10nb.T[j],
             Rrs555_10nb.T[j],
             Rrs660_10nb.T[j],
             Rrs680_10nb.T[j]]))


# result = np.vstack([np.array(lon10),
#                         np.array(lat10),
#                         np.array(chla10),
#                         np.array(Rrs412_10),
#                         np.array(Rrs443_10),
#                         np.array(Rrs490_10),
#                         np.array(Rrs555_10),
#                         np.array(Rrs660_10),
#                         np.array(Rrs680_10)])

a_meannb=[np.mean(anb[:,0]),np.mean(anb[:,1]),np.mean(anb[:,2]),np.mean(anb[:,3]),np.mean(anb[:,4]),np.mean(anb[:,5])]
a_stdnb=[np.std(anb[:,0]),np.std(anb[:,1]),np.std(anb[:,2]),np.std(anb[:,3]),np.std(anb[:,4]),np.std(anb[:,5])]
#chla_meannb=np.mean(chla10nb)
aph_meannb=[np.mean(aphnb[:,0]),np.mean(aphnb[:,1]),np.mean(aphnb[:,2]),np.mean(aphnb[:,3]),np.mean(aphnb[:,4]),np.mean(aphnb[:,5])]
aph_stdnb=[np.std(aphnb[:,0]),np.std(aphnb[:,1]),np.std(aphnb[:,2]),np.std(aphnb[:,3]),np.std(aphnb[:,4]),np.std(aphnb[:,5])]

adg_meannb=[np.mean(adgnb[:,0]),np.mean(adgnb[:,1]),np.mean(adgnb[:,2]),np.mean(adgnb[:,3]),np.mean(adgnb[:,4]),np.mean(adgnb[:,5])]
adg_stdnb=[np.std(adgnb[:,0]),np.std(adgnb[:,1]),np.std(adgnb[:,2]),np.std(adgnb[:,3]),np.std(adgnb[:,4]),np.std(adgnb[:,5])]

bbp_meannb=[np.mean(bbpnb[:,0]),np.mean(bbpnb[:,1]),np.mean(bbpnb[:,2]),np.mean(bbpnb[:,3]),np.mean(bbpnb[:,4]),np.mean(bbpnb[:,5])]
bbp_stdnb=[np.std(bbpnb[:,0]),np.std(bbpnb[:,1]),np.std(bbpnb[:,2]),np.std(bbpnb[:,3]),np.std(bbpnb[:,4]),np.std(bbpnb[:,5])]

Rrs_meannb=[np.mean(Rrs412_10nb),np.mean(Rrs443_10nb),np.mean(Rrs490_10nb),np.mean(Rrs555_10nb),np.mean(Rrs660_10nb),np.mean(Rrs680_10nb)]
Rrs_stdnb=[np.std(Rrs412_10nb),np.std(Rrs443_10nb),np.std(Rrs490_10nb),np.std(Rrs555_10nb),np.std(Rrs660_10nb),np.std(Rrs680_10nb)]

rrs_meannb=[np.mean(rrsnb[:,0]),np.mean(rrsnb[:,1]),np.mean(rrsnb[:,2]),np.mean(rrsnb[:,3]),np.mean(rrsnb[:,4]),np.mean(rrsnb[:,5])]
rrs_stdnb=[np.std(rrsnb[:,0]),np.std(rrsnb[:,1]),np.std(rrsnb[:,2]),np.std(rrsnb[:,3]),np.std(rrsnb[:,4]),np.std(rrsnb[:,5])]

u_meannb=[np.mean(unb[:,0]),np.mean(unb[:,1]),np.mean(unb[:,2]),np.mean(unb[:,3]),np.mean(unb[:,4]),np.mean(unb[:,5])]
u_stdnb=[np.std(unb[:,0]),np.std(unb[:,1]),np.std(unb[:,2]),np.std(unb[:,3]),np.std(unb[:,4]),np.std(unb[:,5])]
#result=np.concatenate((result.T,aph,bbp,adg,a,rrs,u),axis=1)
#np.savetxt('G2019248011640.csv', result, delimiter=',')

##Blooming

nc1=nc4.Dataset('I:/Nagoya University/Project/HAB/satellite/GOCI/G2019248011640.L2_COMS_OC_reprojected.nc','r')
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
# chla10index=np.where((chl_a>=10))
# chla10=chl_a[chla10index[0],chla10index[1]]
# lat10=lat[chla10index[0]]
# lon10=lon[chla10index[1]]
# Rrs412_10=Rrs412[chla10index[0],chla10index[1]]
# Rrs443_10=Rrs443[chla10index[0],chla10index[1]]
# Rrs490_10=Rrs490[chla10index[0],chla10index[1]]
# Rrs555_10=Rrs555[chla10index[0],chla10index[1]]
# Rrs660_10=Rrs660[chla10index[0],chla10index[1]]
# Rrs680_10=Rrs680[chla10index[0],chla10index[1]]
#还要再去掉Rrs是负值的地方
Rrsindex=np.where((Rrs412>0)&(Rrs443>0)&(Rrs490>0)&(Rrs555>0)&(Rrs660>0)&(Rrs680>0))
chla10=chl_a[Rrsindex]
lat10=lat[Rrsindex[0]]
lon10=lon[Rrsindex[1]]
Rrs412_10=Rrs412[Rrsindex]
Rrs443_10=Rrs443[Rrsindex]
Rrs490_10=Rrs490[Rrsindex]
Rrs555_10=Rrs555[Rrsindex]
Rrs660_10=Rrs660[Rrsindex]
Rrs680_10=Rrs680[Rrsindex]
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
a_std=[np.std(a[:,0]),np.std(a[:,1]),np.std(a[:,2]),np.std(a[:,3]),np.std(a[:,4]),np.std(a[:,5])]
#chla_meannb=np.mean(chla10nb)
aph_mean=[np.mean(aph[:,0]),np.mean(aph[:,1]),np.mean(aph[:,2]),np.mean(aph[:,3]),np.mean(aph[:,4]),np.mean(aph[:,5])]
aph_std=[np.std(aph[:,0]),np.std(aph[:,1]),np.std(aph[:,2]),np.std(aph[:,3]),np.std(aph[:,4]),np.std(aph[:,5])]

adg_mean=[np.mean(adg[:,0]),np.mean(adg[:,1]),np.mean(adg[:,2]),np.mean(adg[:,3]),np.mean(adg[:,4]),np.mean(adg[:,5])]
adg_std=[np.std(adg[:,0]),np.std(adg[:,1]),np.std(adg[:,2]),np.std(adg[:,3]),np.std(adg[:,4]),np.std(adg[:,5])]

bbp_mean=[np.mean(bbp[:,0]),np.mean(bbp[:,1]),np.mean(bbp[:,2]),np.mean(bbp[:,3]),np.mean(bbp[:,4]),np.mean(bbp[:,5])]
bbp_std=[np.std(bbp[:,0]),np.std(bbp[:,1]),np.std(bbp[:,2]),np.std(bbp[:,3]),np.std(bbp[:,4]),np.std(bbp[:,5])]

Rrs_mean=[np.mean(Rrs412_10),np.mean(Rrs443_10),np.mean(Rrs490_10),np.mean(Rrs555_10),np.mean(Rrs660_10),np.mean(Rrs680_10)]
Rrs_std=[np.std(Rrs412_10),np.std(Rrs443_10),np.std(Rrs490_10),np.std(Rrs555_10),np.std(Rrs660_10),np.std(Rrs680_10)]

rrs_mean=[np.mean(rrs[:,0]),np.mean(rrs[:,1]),np.mean(rrs[:,2]),np.mean(rrs[:,3]),np.mean(rrs[:,4]),np.mean(rrs[:,5])]
rrs_std=[np.std(rrs[:,0]),np.std(rrs[:,1]),np.std(rrs[:,2]),np.std(rrs[:,3]),np.std(rrs[:,4]),np.std(rrs[:,5])]

u_mean=[np.mean(u[:,0]),np.mean(u[:,1]),np.mean(u[:,2]),np.mean(u[:,3]),np.mean(u[:,4]),np.mean(u[:,5])]
u_std=[np.std(u[:,0]),np.std(u[:,1]),np.std(u[:,2]),np.std(u[:,3]),np.std(u[:,4]),np.std(u[:,5])]
#result=np.concatenate((result.T,aph,bbp,adg,a,rrs,u),axis=1)
#np.savetxt('G2019248011640.csv', result, delimiter=',')
wave=[412, 443, 490, 555, 660, 680]
plt.figure(1)

p1,=plt.errorbar(wave,a_mean,yerr=a_std,color="r", linestyle="-", marker="^", linewidth=1)
p2,=plt.errorbar(wave,a_meannb,yerr=a_stdnb,color="b", linestyle="-", marker="o", linewidth=1)
plt.xlabel("wavelength")
plt.ylabel("a")
plt.title("a")
plt.legend([p1,p2],["blooming water","non-blooming water"],loc='upper left')
plt.show()

plt.figure(2)
p3,=plt.errorbar(wave,aph_mean,yerr=aph_std,color="r", linestyle="-", marker="^", linewidth=1)
p4,=plt.errorbar(wave,aph_meannb,yerr=aph_stdnb,color="b", linestyle="-", marker="^", linewidth=1)
plt.xlabel("wavelength")
plt.ylabel("aph")
plt.title("aph")
plt.legend([p3,p4],["blooming water","non-blooming water"],loc='upper left')
plt.show()
#TODO 有一点问题
plt.figure(3)
p5,=plt.errorbar(wave,adg_mean,yerr=adg_std,color="r", linestyle="-", marker="^", linewidth=1)
p6,=plt.errorbar(wave,adg_meannb,yerr=adg_stdnb,color="b", linestyle="-", marker="^", linewidth=1)
plt.xlabel("wavelength")
plt.ylabel("adg")
plt.title("adg")
plt.legend([p5,p6],["blooming water","non-blooming water"],loc='upper left')
plt.show()

plt.figure(4)
p7,=plt.errorbar(wave,bbp_mean,yerr=bbp_std,color="r", linestyle="-", marker="^", linewidth=1)
p8,=plt.errorbar(wave,bbp_meannb,yerr=bbp_stdnb,color="b", linestyle="-", marker="^", linewidth=1)
plt.xlabel("wavelength")
plt.ylabel("bbp")
plt.title("bbp")
plt.legend([p7,p8],["blooming water","non-blooming water"],loc='upper left')
plt.show()

plt.figure(5)
p9,=plt.errorbar(wave,Rrs_mean,yerr=Rrs_std,color="r", linestyle="-", marker="^", linewidth=1)
p10,=plt.errorbar(wave,Rrs_meannb,yerr=Rrs_stdnb,color="b", linestyle="-", marker="^", linewidth=1)
plt.xlabel("wavelength")
plt.ylabel("Rrs")
plt.title("Rrs")
plt.legend([p9,p10],["blooming water","non-blooming water"],loc='upper left')
plt.show()

plt.figure(6)
p11,=plt.errorbar(wave,rrs_mean,yerr=rrs_std,color="r", linestyle="-", marker="^", linewidth=1)
p12,=plt.errorbar(wave,rrs_meannb,yerr=rrs_std,color="b", linestyle="-", marker="^", linewidth=1)
plt.xlabel("wavelength")
plt.ylabel("rrs")
plt.title("rrs")
plt.legend([p11,p12],["blooming water","non-blooming water"],loc='upper left')
plt.show()

plt.figure(7)
p13,=plt.errorbar(wave,u_mean,yerr=u_std,color="r", linestyle="-", marker="^", linewidth=1)
p14,=plt.errorbar(wave,u_meannb,yerr=u_stdnb,color="b", linestyle="-", marker="^", linewidth=1)
plt.xlabel("wavelength")
plt.ylabel("u")
plt.title("u")
plt.legend([p13,p14],["blooming water","non-blooming water"],loc='upper left')
plt.show()


