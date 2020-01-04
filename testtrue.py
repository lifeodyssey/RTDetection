import netCDF4 as nc4
from mpl_toolkits.basemap import Basemap
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import os
import glob
import matplotlib.pyplot as plt
from QAAV6GOCI import QAAv6
import xlrd
checkfile=xlrd.open_workbook('I:/Nagoya University/Project/QAAv6/【Data】Ise-Bay_Dataset_MODIS_Rrs_matchup.xlsx')
booksheet = checkfile.sheet_by_name('aph&Rrs')
j=booksheet.nrows
RRS412=np.array(booksheet.col_values(20), dtype=np.float32)[1:j]
RRS443=np.array(booksheet.col_values(21), dtype=np.float32)[1:j]
RRS490=np.array(booksheet.col_values(23), dtype=np.float32)[1:j]
RRS555=np.array(booksheet.col_values(26), dtype=np.float32)[1:j]
RRS660=np.array(booksheet.col_values(28), dtype=np.float32)[1:j]
RRS680=np.array(booksheet.col_values(29), dtype=np.float32)[1:j]
n=j-1
rrs=np.ones([n,6])
u=np.ones([n,6])
a=np.ones([n,6])
adg=np.ones([n,6])
aph=np.ones([n,6])
bbp=np.ones([n,6])
for i in range(n):
    [rrs[i,:],u[i,:],a[i,:],adg[i,:],aph[i,:],bbp[i,:]]=QAAv6(np.array([RRS412[i],RRS443[i],RRS490[i],RRS555[i],RRS660[i],RRS680[i]]))

a_mean=[np.mean(a[:,0]),np.mean(a[:,1]),np.mean(a[:,2]),np.mean(a[:,3]),np.mean(a[:,4]),np.mean(a[:,5])]

aph_mean=[np.mean(aph[:,0]),np.mean(aph[:,1]),np.mean(aph[:,2]),np.mean(aph[:,3]),np.mean(aph[:,4]),np.mean(aph[:,5])]
adg_mean=[np.mean(adg[:,0]),np.mean(adg[:,1]),np.mean(adg[:,2]),np.mean(adg[:,3]),np.mean(adg[:,4]),np.mean(adg[:,5])]
bbp_mean=[np.mean(bbp[:,0]),np.mean(bbp[:,1]),np.mean(bbp[:,2]),np.mean(bbp[:,3]),np.mean(bbp[:,4]),np.mean(bbp[:,5])]

rrs_mean=[np.mean(rrs[:,0]),np.mean(rrs[:,1]),np.mean(rrs[:,2]),np.mean(rrs[:,3]),np.mean(rrs[:,4]),np.mean(rrs[:,5])]
u_mean=[np.mean(u[:,0]),np.mean(u[:,1]),np.mean(u[:,2]),np.mean(u[:,3]),np.mean(u[:,4]),np.mean(u[:,5])]
wave=[412, 443, 490, 555, 660, 680]
plt.figure(1)
plt.subplot(3,2,1)
plt.plot(wave,a_mean,color="r", linestyle="-", marker="^", linewidth=1)
plt.xlabel("wavelength")
plt.ylabel("a")
plt.title("a")

plt.subplot(3,2,2)
plt.plot(wave,aph_mean,color="r", linestyle="-", marker="^", linewidth=1)
plt.xlabel("wavelength")
plt.ylabel("aph")
plt.title("aph")

plt.subplot(3,2,3)
plt.plot(wave,adg_mean,color="r", linestyle="-", marker="^", linewidth=1)
plt.xlabel("wavelength")
plt.ylabel("adg")
plt.title("adg")

plt.subplot(3,2,4)
plt.plot(wave,bbp_mean,color="r", linestyle="-", marker="^", linewidth=1)
plt.xlabel("wavelength")
plt.ylabel("bbp")
plt.title("bbp")



plt.subplot(3,2,5)
plt.plot(wave,rrs_mean,color="r", linestyle="-", marker="^", linewidth=1)
plt.xlabel("wavelength")
plt.ylabel("rrs")
plt.title("rrs")

plt.subplot(3,2,6)
plt.plot(wave,u_mean,color="r", linestyle="-", marker="^", linewidth=1)
plt.xlabel("wavelength")
plt.ylabel("u")
plt.title("u")

plt.show()
from sklearn import linear_model
aph412=np.array(booksheet.col_values(9), dtype=np.float32)[1:j]
aph443=np.array(booksheet.col_values(10), dtype=np.float32)[1:j]
aph490=np.array(booksheet.col_values(12), dtype=np.float32)[1:j]
aph555=np.array(booksheet.col_values(15), dtype=np.float32)[1:j]
aph660=np.array(booksheet.col_values(17), dtype=np.float32)[1:j]
aph680=np.array(booksheet.col_values(18), dtype=np.float32)[1:j]

f2=plt.figure(2)
regr = linear_model.LinearRegression()
regr.fit(aph[:,0].reshape(-1,1),aph412.reshape(-1,1))

plt.title('aph412')
plt.xlabel('QAAv6 derived')
plt.ylabel('in-situ')
plt.scatter(aph[:,0],aph412)
p412=regr.predict(aph[:,0].reshape(-1,1))
r2_412=regr.score(aph[:,0].reshape(-1,1),aph412.reshape(-1,1))
plt.plot(aph[:,0],p412,'r','-')
plt.text(min(aph[:,0]),max(p412),'r2='+str(r2_412))
plt.plot()
plt.show()

f3=plt.figure(3)
regr.fit(aph[:,1].reshape(-1,1),aph443.reshape(-1,1))

plt.title('aph443')
plt.xlabel('QAAv6 derived')
plt.ylabel('in-situ')
plt.scatter(aph[:,1],aph443)
p443=regr.predict(aph[:,1].reshape(-1,1))
r2_443=regr.score(aph[:,1].reshape(-1,1),aph443.reshape(-1,1))
plt.plot(aph[:,1],p443,'r','-')
plt.text(min(aph[:,1]),max(p443),'r2='+str(r2_443))
plt.plot()
plt.show()

f4=plt.figure(4)
regr.fit(aph[:,2].reshape(-1,1),aph490.reshape(-1,1))

plt.title('aph490')
plt.xlabel('QAAv6 derived')
plt.ylabel('in-situ')
plt.scatter(aph[:,2],aph490)
p490=regr.predict(aph[:,2].reshape(-1,1))
r2_490=regr.score(aph[:,2].reshape(-1,1),aph490.reshape(-1,1))
plt.plot(aph[:,2],p490,'r','-')
plt.text(min(aph[:,2]),max(p490),'r2='+str(r2_490))
plt.plot()
plt.show()

f5=plt.figure(5)
regr.fit(aph[:,3].reshape(-1,1),aph555.reshape(-1,1))

plt.title('aph555')
plt.xlabel('QAAv6 derived')
plt.ylabel('in-situ')
plt.scatter(aph[:,3],aph555)
p555=regr.predict(aph[:,3].reshape(-1,1))
r2_555=regr.score(aph[:,3].reshape(-1,1),aph555.reshape(-1,1))
plt.plot(aph[:,3],p555,'r','-')
plt.text(min(aph[:,3]),max(p555),'r2='+str(r2_555))
plt.plot()
plt.show()

f6=plt.figure(6)
regr.fit(aph[:,4].reshape(-1,1),aph660.reshape(-1,1))

plt.title('aph660')
plt.xlabel('QAAv6 derived')
plt.ylabel('in-situ')
plt.scatter(aph[:,4],aph660)
p660=regr.predict(aph[:,4].reshape(-1,1))
r2_660=regr.score(aph[:,4].reshape(-1,1),aph660.reshape(-1,1))
plt.plot(aph[:,4],p660,'r','-')
plt.text(min(aph[:,4]),max(p660),'r2='+str(r2_660))
plt.plot()
plt.show()

f7=plt.figure(7)

regr.fit(aph[:,5].reshape(-1,1),aph680.reshape(-1,1))

plt.title('aph680')
plt.xlabel('QAAv6 derived')
plt.ylabel('in-situ')
plt.scatter(aph[:,5],aph680)
p680=regr.predict(aph[:,5].reshape(-1,1))
r2_680=regr.score(aph[:,5].reshape(-1,1),aph680.reshape(-1,1))
plt.plot(aph[:,5],p680,'r','-')
plt.text(min(aph[:,5]),max(p680),'r2='+str(r2_680))
plt.plot()
plt.show()
#
# result=np.concatenate((aph,bbp,adg,a,rrs,u),axis=1)
# np.savetxt('validation.csv', result, delimiter=',')
