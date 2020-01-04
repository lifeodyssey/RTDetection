# read WIM output of MODIS for Ariake to estimate HAB by Fengs algorithm

import netCDF4
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.colors import LogNorm

filename=\
"I:/Nagoya University/Project/HAB/Data/A20190918Ariake.hdf"
print(filename)
nc_file = netCDF4.Dataset(filename, 'r')
print(nc_file)
#nc_var = nc_file.variables.keys()
#print(nc_var)

# set latitude, Longitude of Image
nc_lat  = nc_file.variables['Latitude']
nc_lon  = nc_file.variables['Longitude']
lat = np.array(nc_lat)
lon = np.array(nc_lon)        
print(lat,lon)

south = lat.min()
north = lat.max()
west = lon.min()
east = lon.max()
llat=lat.shape[0]
llon=lat.shape[1]
#print (south, north, west, east,llat,llon)

DNchl = nc_file.variables['chlor_a-Conv inv-remapped to \'MODIS\''][:]
slope = nc_file.variables['chlor_a-Conv inv-remapped to \'MODIS\''].Slope
intercept = nc_file.variables['chlor_a-Conv inv-remapped to \'MODIS\''].Intercept
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

DNchl=np.where((DNchl < 0), 255 + DNchl, DNchl)
print('DNchl=', DNchl.max(), DNchl.min(), DNchl.mean())
chl = 10**(DNchl*slope+intercept)
print('chl=', chl)
print(chl.max(),chl.min(), chl.mean())

# Draw chl image
num_x=1
num_y=1
#fig = plt.figure(figsize=[llon/100, llat/100], dpi=300)
fig = plt.figure()
ax = fig.add_subplot(num_y, num_x, 1)
im = ax.imshow(chl,extent=(lon.min(), lon.max(), lat.min(), lat.max()),
         interpolation='nearest', cmap='jet',norm=LogNorm(vmin=0.1, vmax=70))
ax.tick_params(labelsize=10)
# Draw colorbar
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.1)
cbar = plt.colorbar(im, cax=cax)
# Set colorbar label
unit = 'Chl-a (mg m-3)'
cbar.set_label(unit, rotation=270, labelpad=10.0, fontsize=10)
cbar.ax.tick_params(labelsize=10)
plt.show()
plt.savefig('A20190905chl.png')

# Calculate Rrs555
DN555 = nc_file.variables['Rrs_555 inv-remapped to \'MODIS\''][:]
slope555 = nc_file.variables['Rrs_555 inv-remapped to \'MODIS\''].Slope
intercept555 = nc_file.variables['Rrs_555 inv-remapped to \'MODIS\''].Intercept
Rrs555 = DN555*slope555+intercept555
Rrs555[Rrs555<=0]=np.nan

#print('DN565=', DN565[0:1])
#print('slope565=', slope565)
#print('intercept565=', intercept565)
#print('nlw565=', nlw565)
#print('nlw565=', nlw565.max(), nlw565.min(), nlw565.mean())

# Draw Rrs555 image
num_x=1
num_y=1
fig = plt.figure()
ax = fig.add_subplot(num_y, num_x, 1)
im = ax.imshow(Rrs555,extent=(lon.min(), lon.max(), lat.min(), lat.max()),
         interpolation='nearest', cmap='jet',vmin=0, vmax=0.01)
ax.tick_params(labelsize=10)
# Draw colorbar
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.1)
cbar = plt.colorbar(im, cax=cax)
unit = 'Rrs555'
cbar.set_label(unit, rotation=270, labelpad=10.0, fontsize=10)
cbar.ax.tick_params(labelsize=10)
#plt.show()
plt.savefig('A20190905Rrs555.png')

# Calculate Rrs488, Rrs645, Rrs667, Rrs678
DN488 = nc_file.variables['Rrs_488 inv-remapped to \'MODIS\''][:]
slope488 = nc_file.variables['Rrs_488 inv-remapped to \'MODIS\''].Slope
intercept488 = nc_file.variables['Rrs_488 inv-remapped to \'MODIS\''].Intercept
Rrs488 = DN488*slope488+intercept488

DN645 = nc_file.variables['Rrs_645 inv-remapped to \'MODIS\''][:]
slope645 = nc_file.variables['Rrs_645 inv-remapped to \'MODIS\''].Slope
intercept645 = nc_file.variables['Rrs_645 inv-remapped to \'MODIS\''].Intercept
Rrs645 = DN645*slope645+intercept645

DN667 = nc_file.variables['Rrs_667 inv-remapped to \'MODIS\''][:]
slope667 = nc_file.variables['Rrs_667 inv-remapped to \'MODIS\''].Slope
intercept667 = nc_file.variables['Rrs_667 inv-remapped to \'MODIS\''].Intercept
Rrs667 = DN667*slope667+intercept667

DN678 = nc_file.variables['Rrs_678 inv-remapped to \'MODIS\''][:]
slope678 = nc_file.variables['Rrs_678 inv-remapped to \'MODIS\''].Slope
intercept678 = nc_file.variables['Rrs_678 inv-remapped to \'MODIS\''].Intercept
Rrs678 = DN678*slope678+intercept678

# Calculate Water Type Parameters
diff488_555 = Rrs488-Rrs555
SS645=(Rrs645-Rrs555)-(Rrs667-Rrs555)*(645-555)/(667-555)
#SS645=np.ones((llat,llon))  #mixed銈掋仾銇忋仚
Bbpi555=0.26*Rrs555*Rrs667/(Rrs555-Rrs667)
RBR=(Rrs678)/(Rrs667)
#DHindex=Bbpi555+0.0033*RBR-0.0053 old
DHindex=Bbpi555-0.0019*RBR**(-2.261)

# Drwo Image of diff488_555
num_x=1
num_y=1
fig = plt.figure()
ax = fig.add_subplot(num_y, num_x, 1)
im = ax.imshow(diff488_555,extent=(lon.min(), lon.max(), lat.min(), lat.max()),
         interpolation='nearest', cmap='seismic',vmin=-0.01, vmax=0.01)
ax.tick_params(labelsize=10)
# Draw colorbar
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.1)
cbar = plt.colorbar(im, cax=cax)
unit = 'difference Rrs488-Rrs555'
cbar.set_label(unit, rotation=270, labelpad=10.0, fontsize=10)
cbar.ax.tick_params(labelsize=10)
#plt.show()

# Draw Image of Rrs555-0.007
num_x=1
num_y=1
fig = plt.figure()
ax = fig.add_subplot(num_y, num_x, 1)
im = ax.imshow(Rrs555-0.007,extent=(lon.min(), lon.max(), lat.min(), lat.max()),
         interpolation='nearest', cmap='seismic',vmin=-0.01, vmax=0.01)
ax.tick_params(labelsize=10)
# Draw colorbar
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.1)
cbar = plt.colorbar(im, cax=cax)
unit = 'Rrs555-0.007'
cbar.set_label(unit, rotation=270, labelpad=10.0, fontsize=10)
cbar.ax.tick_params(labelsize=10)
#plt.show()

# define image
Image=np.zeros((llat,llon))
Image[(diff488_555>=0)]=1 #clear
Image[(diff488_555<0) & (Rrs555>=0.007)]=4 #turbid
Image[(diff488_555<0) & (Rrs555<0.007) & (SS645<0)]=3 #mixed
Image[(diff488_555<0) & (Rrs555<0.007) & (SS645>=0) & (DHindex<0)]=6  #Diatom
Image[(diff488_555<0) & (Rrs555<0.007) & (SS645>=0) & (DHindex>=0)]=5 #Raphid
#Image[(diff488_555<0) & (Rrs555<0.007) & (chl<10)]=3 #mixed
#Image[(diff488_555<0) & (Rrs555<0.007) & (chl>=0) & (DHindex<0)]=6  #Diatom
#Image[(diff488_555<0) & (Rrs555<0.007) & (chl>=0) & (DHindex>=0)]=5 #Raphid

# Draw Water Type Image
num_x=1
num_y=1
fig = plt.figure()
ax = fig.add_subplot(num_y, num_x, 1)
im = ax.imshow(Image,extent=(lon.min(), lon.max(), lat.min(), lat.max()),
         interpolation='nearest', cmap='jet', vmin=0, vmax=6)
ax.tick_params(labelsize=10)
# Draw colorbar
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.1)
cbar = plt.colorbar(im, cax=cax)
unit = 'Water type'
cbar.set_label(unit, rotation=270, labelpad=10.0, fontsize=10)
cbar.ax.tick_params(labelsize=10)
#plt.show()
plt.savefig('A20190905wt.png')

# Scattering Diagram of Rrs555 and SS645
fig=plt.figure(figsize=[8,6])
ax=fig.add_axes([0.2,0.2,0.75,0.75])
ax.scatter(Rrs555[Image==1],SS645[Image==1], marker=".", c='blue',label='clear')
ax.scatter(Rrs555[Image==4],SS645[Image==4], marker=".", c='orange',label='turbid')
ax.scatter(Rrs555[Image==3],SS645[Image==3], marker=".", c='green',label='mixed')
ax.scatter(Rrs555[Image==6],SS645[Image==6], marker="o", c='brown',label='diatom')
ax.scatter(Rrs555[Image==5],SS645[Image==5], marker="o", c='red',label='raphid')
x=np.array([0,0.014]) #draw line
y=np.array([0,0])
ax.plot(x,y)
ax.legend(loc=(0.7, 0.05), fontsize=15)
ax.set_xlim([0,0.014])
ax.set_ylim([-0.002,0.002])
ax.tick_params(labelsize=20)
ax.set_xlabel('Rrs555',fontsize=20)
ax.set_ylabel('SS645',fontsize=20)
#plt.show
plt.savefig('A20190905Rrs555_SS.png')

# Scattering Diagram of RBR and Bbpi555
fig=plt.figure(figsize=[8,6])
ax=fig.add_axes([0.2,0.2,0.75,0.75])
ax.scatter(RBR[Image==1],Bbpi555[Image==1], marker=".", c='blue',label='clear')
ax.scatter(RBR[Image==4],Bbpi555[Image==4], marker=".", c='orange',label='turbid')
ax.scatter(RBR[Image==3],Bbpi555[Image==3], marker=".", c='green',label='mixed')
ax.scatter(RBR[Image==6],Bbpi555[Image==6], marker="o", c='brown',label='diatom')
ax.scatter(RBR[Image==5],Bbpi555[Image==5], marker="o", c='red',label='raphid')
#x=np.array([0,0.0053/0.0033]) #draw line
#y=np.array([0.0053,0])
x = np.linspace(0,3,1000)
y=0.0019*x**(-2.261)
ax.plot(x,y)
ax.legend(loc=(0.7, 0.5), fontsize=15)
ax.set_xlim([0,3])
ax.set_ylim([0,0.004])
ax.tick_params(labelsize=20)
ax.set_yticks(np.arange(0, 0.005, 0.001))
ax.set_xlabel('RBR',fontsize=20)
ax.set_ylabel('Bbpi555',fontsize=20)
#plt.show
plt.savefig('A20190905RBR_Bbpi555.png')

# Scattering Diagram of Chl and Rrs555
fig=plt.figure(figsize=[8,6])
ax=fig.add_axes([0.2,0.2,0.75,0.75])
ax.scatter(chl[Image==1],Rrs555[Image==1], marker=".", c='blue',label='clear')
ax.scatter(chl[Image==4],Rrs555[Image==4], marker=".", c='orange',label='turbid')
ax.scatter(chl[Image==3],Rrs555[Image==3], marker=".", c='green',label='mixed')
ax.scatter(chl[Image==6],Rrs555[Image==6], marker="o", c='brown',label='diatom')
ax.scatter(chl[Image==5],Rrs555[Image==5], marker="o", c='red',label='raphid')
ax.set_xscale('log')
#x=np.array([0,0.0053/0.0033]) #draw line
#y=np.array([0.0053,0])
#ax.plot(x,y)
ax.legend(loc=(0.1, 0.5), fontsize=15)
ax.set_xlim([0.1,100])
ax.set_ylim([0,0.014])
ax.tick_params(labelsize=20)
#ax.set_yticks(np.arange(0, 0.005, 0.001))
ax.set_xlabel('Chl',fontsize=20)
ax.set_ylabel('Rrs555',fontsize=20)
#plt.show
plt.savefig('A20190905chl_Rrs555.png')
