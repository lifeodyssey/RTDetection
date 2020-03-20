import netCDF4 as nc4
import numpy as np

from geo_Collection import geo_web as gs
import pandas as pd
from collections import OrderedDict
import json
import geopandas as gpd
import h5py
import os
import glob
os.chdir('I:/Nagoya University/Project/Seto/GOCI')
path='I:/Nagoya University/Project/Seto/GOCI/'
datalist=glob.glob('*.nc')
L=len(datalist)
for j in range(L):
    file=datalist[j]
    filename = path + file
    #print(filename)
    nc_file = nc4.Dataset(filename, 'r')
    T=nc_file.time_coverage_end[0:13]
    T = T.replace('-', '')
    nc_file.close()
    os.rename(filename,path+str(T)+'.nc')
    print('percent: {:.2%}'.format(j+1/L))
