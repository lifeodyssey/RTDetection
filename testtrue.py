import netCDF4 as nc4
from mpl_toolkits.basemap import Basemap
import numpy as np
import pandas as pd

from matplotlib.colors import LogNorm
import os
import glob
import matplotlib.pyplot as plt
from QAAV6 import QAAv6
import xlrd
checkfile=xlrd.open_workbook('I:/Nagoya University/Project/Dataset/Ariake.xlsx')
booksheet = checkfile.sheet_by_name('Test')
j=booksheet.nrows
RRS412=np.array(booksheet.col_values(12)[1:j], dtype=np.float32)
RRS443=np.array(booksheet.col_values(13)[1:j], dtype=np.float32)
RRS469=np.array(booksheet.col_values(14)[1:j], dtype=np.float32)
RRS490=np.array(booksheet.col_values(15)[1:j], dtype=np.float32)
RRS531=np.array(booksheet.col_values(17)[1:j], dtype=np.float32)
RRS547=np.array(booksheet.col_values(18)[1:j], dtype=np.float32)
RRS555=np.array(booksheet.col_values(18)[1:j], dtype=np.float32)
RRS645=np.array(booksheet.col_values(21)[1:j], dtype=np.float32)
RRS667=np.array(booksheet.col_values(22)[1:j], dtype=np.float32)
RRS678=np.array(booksheet.col_values(23)[1:j], dtype=np.float32)
n=j-1
rrs=np.ones([n,10])
u=np.ones([n,10])
a=np.ones([n,10])
adg=np.ones([n,10])
aph=np.ones([n,10])
bbp=np.ones([n,10])
for i in range(n):
    [rrs[i,:],u[i,:],a[i,:],adg[i,:],aph[i,:],bbp[i,:]]=QAAv6(np.array([RRS412[i],RRS443[i],RRS469[i],RRS490[i],RRS531[i],RRS547[i],RRS555[i],RRS645[i],RRS667[i],RRS678[i]]))


wave=[412, 443, 469,488,531,547,555,645,667,678]
from errobarplot import errorplot
plt.figure(1)
plt.subplot(3,2,1)
errorplot(wave,a,'a')

plt.subplot(3,2,2)
errorplot(wave,aph,'aph')

plt.subplot(3,2,3)
errorplot(wave,adg,'adg')

plt.subplot(3,2,4)
errorplot(wave,bbp,'bbp')



plt.subplot(3,2,5)
errorplot(wave,rrs,'rrs')

plt.subplot(3,2,6)
errorplot(wave,u,'u')

plt.show()
from sklearn import linear_model
aph412=np.array(booksheet.col_values(25)[1:j], dtype=np.float32)
aph443=np.array(booksheet.col_values(26)[1:j], dtype=np.float32)
aph469=np.array(booksheet.col_values(27)[1:j], dtype=np.float32)
aph488=np.array(booksheet.col_values(28)[1:j], dtype=np.float32)
aph531=np.array(booksheet.col_values(30)[1:j], dtype=np.float32)
aph547=np.array(booksheet.col_values(31)[1:j], dtype=np.float32)
aph555=np.array(booksheet.col_values(31)[1:j], dtype=np.float32)
aph645=np.array(booksheet.col_values(34)[1:j], dtype=np.float32)
aph667=np.array(booksheet.col_values(35)[1:j], dtype=np.float32)
aph678=np.array(booksheet.col_values(36)[1:j], dtype=np.float32)

f2=plt.figure(2)

from Linearplot import linearplot
linearplot(aph[:,0],aph412,'aph412')

f3=plt.figure(3)
linearplot(aph[:,1],aph443,'aph443')

f4=plt.figure(4)
linearplot(aph[:,2],aph469,'aph469')

f5=plt.figure(5)
linearplot(aph[:,3],aph488,'aph488')

f6=plt.figure(6)
linearplot(aph[:,4],aph531,'aph531')

f7=plt.figure(7)
linearplot(aph[:,5],aph547,'aph547')

f8=plt.figure(8)
linearplot(aph[:,6],aph555,'aph555')

f9=plt.figure(9)
linearplot(aph[:,7],aph645,'aph645')

f10=plt.figure(10)
linearplot(aph[:,8],aph667,'aph667')

f11=plt.figure(11)
linearplot(aph[:,9],aph678,'aph678')

# result=np.concatenate((aph,bbp,adg,a,rrs,u),axis=1)
# np.savetxt('validation.csv', result, delimiter=',')
