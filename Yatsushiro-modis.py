# read WIM output of MODIS for Ariake to estimate HAB by Slope difference algorithm
from mpl_toolkits.basemap import Basemap
import netCDF4
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.colors import LogNorm
import os
import glob
import warnings
warnings.filterwarnings("ignore")#用来去掉Warning,
os.chdir('I:/Nagoya University/Project/HAB/Data')
path='I:/Nagoya University/Project/HAB/Data/'
datalist=glob.glob('A*'+'hdf')
for i in range(len(datalist)):
    file=datalist[i]


    filename=path+file
    print(filename)

    nc_file = netCDF4.Dataset(filename, 'r')
    #print(nc_file)
    #nc_var = nc_file.variables.keys()
    #print(nc_var)

    # set latitude, Longitude of Image
    nc_lat  = nc_file.variables['Latitude']
    nc_lon  = nc_file.variables['Longitude']
    lat = np.array(nc_lat)
    lon = np.array(nc_lon)
    #print(lat,lon)

    south = lat.min()
    north = lat.max()
    west = lon.min()
    east = lon.max()
    llat=lat.shape[0]
    llon=lat.shape[1]

    # Calculate Rrs 412 443 488 555 645 667 678

    # DN412 = nc_file.variables['Rrs_412 inv-remapped to \'MODIS\''][:]
    # slope412 = nc_file.variables['Rrs_412 inv-remapped to \'MODIS\''].Slope
    # intercept412 = nc_file.variables['Rrs_412 inv-remapped to \'MODIS\''].Intercept
    # Rrs412 = DN412*slope412+intercept412
    # Rrs412[Rrs412<=0]=np.nan
    #
    # DN443 = nc_file.variables['Rrs_443 inv-remapped to \'MODIS\''][:]
    # slope443 = nc_file.variables['Rrs_443 inv-remapped to \'MODIS\''].Slope
    # intercept443 = nc_file.variables['Rrs_443 inv-remapped to \'MODIS\''].Intercept
    # Rrs443 = DN443*slope443+intercept443
    # Rrs443[Rrs443<=0]=np.nan

    DN488 = nc_file.variables['Rrs_488 inv-remapped to \'MODIS\''][:]
    slope488 = nc_file.variables['Rrs_488 inv-remapped to \'MODIS\''].Slope
    intercept488 = nc_file.variables['Rrs_488 inv-remapped to \'MODIS\''].Intercept
    Rrs488 = DN488*slope488+intercept488
    Rrs488[Rrs488<=0]=np.nan

    DN555 = nc_file.variables['Rrs_555 inv-remapped to \'MODIS\''][:]
    slope555 = nc_file.variables['Rrs_555 inv-remapped to \'MODIS\''].Slope
    intercept555 = nc_file.variables['Rrs_555 inv-remapped to \'MODIS\''].Intercept
    Rrs555 = DN555*slope555+intercept555
    Rrs555[Rrs555<=0]=np.nan

    DN645 = nc_file.variables['Rrs_645 inv-remapped to \'MODIS\''][:]
    slope645 = nc_file.variables['Rrs_645 inv-remapped to \'MODIS\''].Slope
    intercept645 = nc_file.variables['Rrs_645 inv-remapped to \'MODIS\''].Intercept
    Rrs645 = DN645*slope645+intercept645
    Rrs645[Rrs645<=0]=np.nan

    DN667 = nc_file.variables['Rrs_667 inv-remapped to \'MODIS\''][:]
    slope667 = nc_file.variables['Rrs_667 inv-remapped to \'MODIS\''].Slope
    intercept667 = nc_file.variables['Rrs_667 inv-remapped to \'MODIS\''].Intercept
    Rrs667 = DN667*slope667+intercept667
    Rrs667[Rrs667<=0]=np.nan

    DN678 = nc_file.variables['Rrs_678 inv-remapped to \'MODIS\''][:]
    slope678 = nc_file.variables['Rrs_678 inv-remapped to \'MODIS\''].Slope
    intercept678 = nc_file.variables['Rrs_678 inv-remapped to \'MODIS\''].Intercept
    Rrs678 = DN678*slope678+intercept678
    Rrs555[Rrs678<=0]=np.nan

    #calculate slope
    # slope412_443=(Rrs412-Rrs443)/(412-443)
    # slope443_488=(Rrs443-Rrs488)/(443-488)
    slope488_555=(Rrs488-Rrs555)/(488-555)
    slope555_645=(Rrs555-Rrs645)/(555-645)
    slope645_667=(Rrs645-Rrs667)/(645-667)
    slope667_678=(Rrs667-Rrs678)/(667-678)


    # Calculate Water Type Parameters
    diff488_555 = Rrs488-Rrs555



    # define image
    Image=np.zeros((llat,llon))
    Image[(diff488_555>=0)]=2 #clear
    Image[(diff488_555<0) & (Rrs555>=0.007)]=6 #turbid

    Image[(diff488_555<0) & (Rrs555<0.007) & (slope667_678<0.00002) & (slope555_645<-4.331111111111111e-05)]=3  #Diatom
    Image[(diff488_555<0) & (Rrs555<0.007) & (slope667_678>0.00002) & (slope555_645>=-4.331111111111111e-05)]=5 #Raphid


    # Draw Water Type Image
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)

    ax = Basemap(projection='mill',
                 llcrnrlat=32,
                 llcrnrlon=130,
                 urcrnrlat=33.2,
                 urcrnrlon=131,
                 resolution='f')
    ax.drawcoastlines()
    ax.fillcontinents()
    ax.drawcountries()
    ax.drawrivers(color='#0000ff')
    parallels = np.arange(32, 33.2, 0.4)
    ax.drawparallels(parallels, fontsize=10, linewidth=0.25, dashes=[7, 15],
                     color='k', labels=[1, 0, 1, 1])
    ##dash：http://bbs.06climate.com/home.php?mod=space&uid=8934&do=blog&id=3284 参考
    ##linewidth：线的宽度
    meridians = np.arange(130, 131, 0.2)
    ax.drawmeridians(meridians, fontsize=10, dashes=[7, 15],
                     linewidth=0.3, color='k', labels=[1, 1, 0, 1])
    x, y = ax(lon, lat)

    im = ax.pcolor(x,y,np.squeeze(Image),cmap='jet', vmin=0, vmax=6)

    cbar = ax.colorbar(im,location="right",size="5%",pad=0.1)
    unit = 'Water type'
    cbar.set_label(unit, rotation=270, labelpad=10.0, fontsize=10)
    cbar.ax.tick_params(labelsize=10)
    #plt.show()
    plt.savefig(file+'555_645.png')
