import netCDF4 as nc4
import numpy as np
from collections import OrderedDict
from geo_Collection import geo_web as gs
from geo_Collection import geo_self as gss
from QAAV6GOCI import QAAv6
import os
import glob
import datetime
import warnings
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
os.chdir('/Users/zhenjia/Downloads/')
path = '/Users/zhenjia/Downloads/'
datalist = glob.glob('G2018*.nc')
# 20150715 Uwajima
# 20180723 Osaka
# minlat = 32
# minlon = 130
# maxlat = 35
# maxlon = 136
# area of full seto-inland sea
# # # for a1 in range(1):
# minlat = 32.32
# minlon = 131.39
# maxlat = 33.48
# maxlon = 132.77

# # bun go si i to
minlat = 33.88
minlon = 132.85
maxlat = 34.8
maxlon = 135
# middle
# minlat = 34.18
# minlon = 134.00
# maxlat = 34.88
# maxlon = 135.59
# # Osaka
l = np.arange(len(datalist))


def bbpp(i):
    file = datalist[i]
    filename = path + file
    nc_file = nc4.Dataset(filename, 'r')
    print(nc_file.time_coverage_end)
    lon = nc_file.groups['navigation_data'].variables['longitude'][:]
    lat = nc_file.groups['navigation_data'].variables['latitude'][:]
    variables = nc_file.groups['geophysical_data'].variables

    x = np.arange(minlon, maxlon, 0.005)  # 1 km grid,
    y = np.arange(maxlat, minlat, -0.005)

    for i in variables:
        var = variables[i][:]
        np.where(var <= 0, var, np.nan)
        if i != 'l2_flags':
            var_re, grid = gs.swath_resampling(var, lon, lat, x, y, 2500)  # 1 km grid
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
    Rrs_490 = variables['Rrs_490']

    Rrs_555 = variables['Rrs_555']

    Rrs_660 = variables['Rrs_660']
    Rrs_680 = variables['Rrs_680']

    Rrs = np.array([
        Rrs_412,
        Rrs_443,

        Rrs_490,

        Rrs_555,
        Rrs_660,
        Rrs_680])
    # nflh = nlw678 - (70 / 81) * nlw667 - (11 / 81) * nlw748
    chl = variables['chlor_a']

    def rrs_to_qaa(Rrs):
        bbp555QAA = np.zeros(np.shape(chl))

        for m in range(len(y)):
            for n in range(len(x)):
                # for k in range(10):
                # if (Rrs[k, i, j] < 0) or (Rrs[k, i, j] == np.nan):
                # bbp[k] = np.nan
                # else:
                bbp = QAAv6(Rrs[:, m, n])
                bbp555QAA[m, n] = bbp[3]

        return bbp555QAA

    bbp555qaa = rrs_to_qaa(Rrs)
    bbpindex = 0.37 * (Rrs_555 * Rrs_660) / (Rrs_555 - Rrs_660)
    chlanorm_bbp555qaa = bbp555qaa / chl
    chlanorm_bbpindex = bbpindex / chl
    gs.plot_geo_image(chlanorm_bbp555qaa, lons, lats, log10=False,
                      title='chlanorm_bbp555qaa' + nc_file.time_coverage_end[0:13],
                      save_image='chlanorm_bbp555qaa' + nc_file.time_coverage_end[0:13])
    gs.plot_geo_image(chlanorm_bbpindex, lons, lats, log10=False,
                      title='chlanorm_bbpindex' + nc_file.time_coverage_end[0:13],
                      save_image='chlanorm_bbpindex' + nc_file.time_coverage_end[0:13])
    gs.plot_geo_image(bbp555qaa, lons, lats, log10=False, title='bbp555qaa' + nc_file.time_coverage_end[0:13],
                      caxis=[0, 0.03], save_image='bbp555qaa' + nc_file.time_coverage_end[0:13])
    gs.plot_geo_image(bbpindex, lons, lats, log10=False, title='bbp555index' + nc_file.time_coverage_end[0:13],
                      save_image='bbp555index' + nc_file.time_coverage_end[0:13])


if __name__ == '__main__':
    mp.set_start_method("spawn")
    pool = mp.Pool(processes=8)
    pool.map(bbpp, l)
