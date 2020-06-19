import netCDF4 as nc4
import numpy as np
from mpl_toolkits.basemap import Basemap
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from geo_Collection import geo_self as gs
from matplotlib.collections import PatchCollection
import os
import glob
import warnings
import cmocean
from matplotlib.patches import Polygon

warnings.filterwarnings("ignore")  # 用来去掉Warning,
import pandas as pd
import matplotlib.pyplot as plt

file = "/Users/zhenjia/Downloads/GEBCO_2020_16_Jun_2020_c364c2ffbd27/gebco.nc"
# print(filename
ele = nc4.Dataset(file, 'r')
minlat = 32.5
minlon = 130.5
maxlat = 35
maxlon = 136

# North of Oita
# minlon=130.9
# maxlon=131.75
# minlat=33.58
# maxlat=34.08

NO = {'minlon': 130.9,
      'maxlon': 131.75,
      'minlat': 33.54,
      'maxlat': 34.08}

# East of oita
# minlon=131.48
# maxlon=132.15
# maxlat=33.454
# minlat=32.74

EO = {'minlon': 131.48,
      'maxlon': 132.15,
      'minlat': 32.74,
      'maxlat': 33.454}

# Uwajima coast
# minlon=131.80
# maxlon=132.60
# minlat=32.91
# maxlat=33.36

UC = {'minlon': 132.20,
      'maxlon': 132.60,
      'minlat': 32.80,
      'maxlat': 33.36}

# osaka region
# minlon=134.02
# maxlon=135.446
# maxlat=34.775
# minlat=34.20
OR = {'minlon': 134.02,
      'maxlon': 135.45,
      'minlat': 34.20,
      'maxlat': 34.78}

lat = np.asarray(ele.variables['lat'])
lon = np.asarray(ele.variables['lon'])
bathymetry = -np.asarray(ele.variables['elevation'])
lon, lat = np.meshgrid(lon, lat)

lon_0 = (lon.min() + lon.max()) / 2
lat_0 = (lat.min() + lat.max()) / 2

m = Basemap(llcrnrlon=minlon, llcrnrlat=minlat,
            urcrnrlon=maxlon, urcrnrlat=maxlat,
            resolution='f', lon_0=lon_0, lat_0=lat_0, projection='merc')

x2d, y2d = m(lon, lat)

fig = plt.figure(figsize=(12, 12 * m.aspect))
ax = fig.add_axes([0.08, 0.1, 0.7, 0.7], facecolor='white')


cmn, cmx = 0, 200

cmap = cmocean.cm.deep
p = m.pcolor(x2d, y2d, bathymetry, vmin=cmn, vmax=cmx, cmap=cmap)
# plt.title(title, fontsize=24)

cb = m.colorbar(p, location="right", size="5%", pad=0.1)  # draw colorbar

cb.ax.tick_params(labelsize=10)
label = 'Bathymetry(m)'
cb.set_label(label, labelpad=10.0, fontsize=10)


m.drawcoastlines()
m.drawcountries()
m.fillcontinents(color='lightgray')
m.drawrivers()

parallels = np.arange(minlat, maxlat, 0.5)
meridians = np.arange(minlon, maxlon, 1)
###The four positions are [left, right, top, bottom]


m.drawmeridians(meridians, fontsize=10, labels=[0, 0, 0, 1])
m.drawparallels(parallels, fontsize=10, labels=[1, 0, 0, 0])

dict = {
    'EO': EO,
    'NO': NO,
    'UC': UC,
    'OR': OR
}


def plot_rectangle(bmap, lonmin, lonmax, latmin, latmax):
    xs = [lonmin, lonmax, lonmax, lonmin, lonmin]
    ys = [latmin, latmin, latmax, latmax, latmin]
    bmap.plot(xs, ys, latlon=True,color='r')


for key in dict:
    area = dict[key]
    plot_rectangle(m, area['minlon'], area['maxlon'], area['minlat'], area['maxlat'])

plt.show()
