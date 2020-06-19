import netCDF4 as nc4
import numpy as np
from collections import OrderedDict
from geo_Collection import geo_web as gs
from QAAV6 import QAAv6
import os
import glob
import datetime
import warnings
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
from matplotlib import colors
from matplotlib.colors import ListedColormap
import earthpy.plot as ep
import math
from scipy.interpolate import RectSphereBivariateSpline
from matplotlib.colors import LogNorm

import concurrent.futures
import math
from deco import *
import multiprocessing as mp

# from multiprocess import Pool
warnings.filterwarnings("ignore")  # 用来去掉Warning,
os.environ['PROJ_LIB'] = '/Users/zhenjia/opt/anaconda3/share/proj'
import pandas as pd
import matplotlib.pyplot as plt

# os.chdir('I:/Nagoya University/Project/Seto/MODIS')
# path='I:/Nagoya University/Project/Seto/MODIS/'
os.chdir('/Users/zhenjia/Desktop/Project/Seto/MODIS/l1/l2_lac')
path = '/Users/zhenjia/Desktop/Project/Seto/MODIS/l1/l2_lac/'
datalist = glob.glob('2008*.nc')
# 20110715 Uwajima
# 20180723 Osaka
minlat = 32.5
minlon = 130
maxlat = 35
maxlon = 136
# area of full seto-inland sea
# # # for a1 in range(1):
# minlat = 32.32
# minlon = 131.39
# maxlat = 33.48
# maxlon = 132.77

# # bun go si i to
# minlat = 33.88
# minlon = 132.85
# maxlat = 34.73
# maxlon = 133.837
# middle
# minlat = 34.18
# minlon = 134.00
# maxlat = 34.88
# maxlon = 135.59
# # Osaka
L = len(datalist)
initime = datetime.datetime.now()
#
# # plt.subplot(2,3,3)
# # with concurrent.futures.ProcessPoolExecutor() as executor:
# a1 = np.arange(len(datalist))
#
#
# #
#
#
# def main(a1):
file = datalist[0]
filename = path + file
nc_file = nc4.Dataset(filename, 'r')
print(nc_file.time_coverage_end)
lon = nc_file.groups['navigation_data'].variables['longitude'][:]
lat = nc_file.groups['navigation_data'].variables['latitude'][:]
variables = nc_file.groups['geophysical_data'].variables

x = np.arange(minlon, maxlon, 0.01)  # 1 km grid,
y = np.arange(maxlat, minlat, -0.01)

for i in variables:
    var = variables[i][:]
    np.where(var <= 0, var, np.nan)
    if i != 'l2_flags':
        var_re, grid = gs.swath_resampling(var, lon, lat, x, y, 1000)  # 1 km grid
        # var_re=var_re.filled()
        if np.ma.is_mask(var_re):
            var_re = var_re.filled(np.nan)
            var_re[var_re == -32767.0] = np.nan
        variables[i] = var_re
    else:
        # var_re = var_re.filled()
        variables[i] = var

lons = grid.lons
lats = grid.lats
Rrs_412 = variables['Rrs_412']
Rrs_443 = variables['Rrs_443']
Rrs_469 = variables['Rrs_469']
Rrs_488 = variables['Rrs_488']
Rrs_531 = variables['Rrs_531']
Rrs_547 = variables['Rrs_547']
Rrs_555 = variables['Rrs_555']
Rrs_645 = variables['Rrs_645']
Rrs_667 = variables['Rrs_667']
Rrs_678 = variables['Rrs_678']
Rrs_748 = variables['Rrs_748']
# F0 = [172.912, 187.622, 205.878, 194.933, 185.747, 186.539, 183.869, 157.811, 152.255, 148.052, 128.065]
# nlw412 = Rrs_412 * F0[0]
# nlw443 = Rrs_443 * F0[1]
# nlw469 = Rrs_469 * F0[2]
# nlw488 = Rrs_488 * F0[3]
# nlw531 = Rrs_531 * F0[4]
# nlw547 = Rrs_547 * F0[5]
# nlw555 = Rrs_555 * F0[6]
# nlw645 = Rrs_645 * F0[7]
# nlw667 = Rrs_667 * F0[8]
# nlw678 = Rrs_678 * F0[9]
# nlw748 = Rrs_748 * F0[10]
# # nlw412 = variables['nLw_412']
# # nlw443 = variables['nLw_443']
# # nlw469 = variables['nLw_469']
# # nlw488 = variables['nLw_488']
# # nlw531 = variables['nLw_531']
# # nlw547 = variables['nLw_547']
# # nlw555 = variables['nLw_555']
# # nlw645 = variables['nLw_645']
# # nlw667 = variables['nLw_667']
# # nlw678 = variables['nLw_678']
# # nlw748 = variables['nLw_748']
#
# # Rrs_412 = variables['Rrs_412'].filled(np.nan)
# # Rrs_443 = variables['Rrs_443'].filled(np.nan)
# # Rrs_469 = variables['Rrs_469'].filled(np.nan)
# # Rrs_488 = variables['Rrs_488'].filled(np.nan)
# # Rrs_531 = variables['Rrs_531'].filled(np.nan)
# # Rrs_547 = variables['Rrs_547'].filled(np.nan)
# # Rrs_555 = variables['Rrs_555'].filled(np.nan)
# # Rrs_645 = variables['Rrs_645'].filled(np.nan)
# # Rrs_667 = variables['Rrs_667'].filled(np.nan)
# # Rrs_678 = variables['Rrs_678'].filled(np.nan)
# # Rrs_748 = variables['Rrs_748'].filled(np.nan)
# # nlw412 = variables['nLw_412'].filled(np.nan)
# # nlw443 = variables['nLw_443'].filled(np.nan)
# # nlw469 = variables['nLw_469'].filled(np.nan)
# # nlw488 = variables['nLw_488'].filled(np.nan)
# # nlw531 = variables['nLw_531'].filled(np.nan)
# # nlw547 = variables['nLw_547'].filled(np.nan)
# # nlw555 = variables['nLw_555'].filled(np.nan)
# # nlw645 = variables['nLw_645'].filled(np.nan)
# # nlw667 = variables['nLw_667'].filled(np.nan)
# # nlw678 = variables['nLw_678'].filled(np.nan)
# # nlw748 = variables['nLw_748'].filled(np.nan)
Rrs = np.array([
    Rrs_412,
    Rrs_443,
    Rrs_469,
    Rrs_488,
    Rrs_531,
    Rrs_547,
    Rrs_555,
    Rrs_645,
    Rrs_667,
    Rrs_678])
# nflh = nlw678 - (70 / 81) * nlw667 - (11 / 81) * nlw748
chl = variables['chlor_a']

wt = np.ma.array(np.zeros(np.shape(chl)))
wt[chl > 10] = 3
wt[(chl < 10) & (chl > 5)] = 2
wt[(chl < 5)] = 1
wt[chl.mask] = 0


def plot_WT_image(sds: np.ma.array, lon: np.ndarray, lat: np.ndarray, title: str = None,
                  lon_range: list = None, lat_range: list = None, save_image: str = None,
                  dpi: int = 400):
    if len(lon.shape) == 1:
        print('MeshGridding...')
        lon, lat = np.meshgrid(lon, lat)

    lon_0 = (lon.min() + lon.max()) / 2
    lat_0 = (lat.min() + lat.max()) / 2

    print(f'Lat: [{lat.min():.3f}, {lat.max():.3f}] | '
          f'Lon: [{lon.min():.3f}, {lon.max():.3f}] | '
          f'SDS: [{sds.min():.3f}, {sds.max():.3f}]')

    if (lon_range is not None) and (lat_range is not None):
        m = Basemap(llcrnrlon=min(lon_range), llcrnrlat=min(lat_range),
                    urcrnrlon=max(lon_range), urcrnrlat=max(lat_range),
                    resolution='f', lon_0=lon_0, lat_0=lat_0, projection='tmerc')
    else:
        m = Basemap(llcrnrlon=lon.min(), llcrnrlat=lat.min(),
                    urcrnrlon=lon.max(), urcrnrlat=lat.max(),
                    resolution='f', lon_0=lon_0, lat_0=lat_0, projection='tmerc')
    x2d, y2d = m(lon, lat)

    fig = plt.figure(figsize=(12, 12 * m.aspect))
    ax = fig.add_axes([0.08, 0.1, 0.7, 0.7], facecolor='w')

    if (lon_range is not None) and (lat_range is not None):
        parallels = np.arange(min(lat_range), max(lat_range), 0.5)
        meridians = np.arange(min(lon_range), max(lon_range), 0.5)
    else:
        parallels = np.arange(lat.min(), lat.max(), 0.5)
        meridians = np.arange(lon.min(), lon.max(), 0.5)
    labels = ['No Data',
              'Other Water',
              'Other Phytoplankton',
              'K.mikimotoi Bloom']
    colors = ['w',
              'deepskyblue',
              'g',
              'r'
              ]
    cmap = ListedColormap(colors)
    p = m.pcolor(x2d, y2d, sds, cmap=cmap)

    if title is not None:
        plt.title(title, fontsize=24)

    plt.plot([], [], label=labels[0], color=colors[0])
    plt.plot([], [], label=labels[1], color=colors[1])
    plt.plot([], [], label=labels[2], color=colors[2])
    plt.plot([], [], label=labels[3], color=colors[3])
    plt.legend()

    m.drawcoastlines()
    m.drawcountries()
    m.fillcontinents(color='lightgray')
    m.drawrivers()

    m.drawmeridians(meridians, fontsize=10, linewidth=0.25, dashes=[7, 15],
                    color='k', labels=[1, 0, 1, 1])
    m.drawparallels(parallels, fontsize=10, dashes=[7, 15],
                    linewidth=0.3, color='k', labels=[1, 1, 0, 1])
    plt.show()

    if save_image is not None:
        plt.savefig(save_image, dpi=dpi, facecolor='w', edgecolor='w', orientation='portrait')
        plt.show()
        plt.close()


plot_WT_image(wt, lons, lats)
