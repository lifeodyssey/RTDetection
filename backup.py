import pandas as pd
import netCDF4 as nc4
import numpy as np

from geo_Collection import geo_self as gs

import os
import glob
import datetime
import warnings
import multiprocessing as mp
from geo_Collection import K_malgorithm as km
import pandas as pd
from datetime import datetime

os.environ['PROJ_LIB'] = '/Users/zhenjia/opt/anaconda3/share/proj'

obs500 = pd.read_excel('/Users/zhenjia/Desktop/Field_zheng2_OR_seperate.xlsx', sheet_name='500')
warnings.filterwarnings("ignore")  # 用来去掉Warning,

import matplotlib.pyplot as plt

# os.chdir('I:/Nagoya University/Project/Seto/MODIS')
# path='I:/Nagoya University/Project/Seto/MODIS/'
os.chdir('/Users/zhenjia/Desktop/Project/Seto/MODIS/l2')
path = '/Users/zhenjia/Desktop/Project/Seto/MODIS/l2/'
datalist = glob.glob('*.nc')

#
# def datasorter(datalist, obs):
#     list=datalist
#     l = len(datalist)
#     for i in range(0, l, ):
#         T = datalist[i]
#         ddd = T[1:8]
#         parsed = datetime.strptime(ddd, '%Y%j')
#         day = datetime.strftime(parsed, '%Y%m%d')
#         if day not in obs['DATE']:
#             list.pop(i)
#     return list
#
#
# datalist = datasorter(datalist, obs500)
obs500['DATE'] = obs500['DATE'].astype('str')
l = len(datalist)
for i in range(l - 1, -1, -1):
    T = datalist[i]
    ddd = T[1:8]
    parsed = datetime.strptime(ddd, '%Y%j')
    day = datetime.strftime(parsed, '%Y%m%d')
    if day in obs500['DATE'].values:
        print(day)
    else:
        datalist.pop(i)

minlat = 32.5
minlon = 131
maxlat = 35
maxlon = 136

# area = [
#     'EO',
#     'NO',
#     'UC',
#     'OR']
#
area = [
    'OB',
    'HN']

methods = ['bbp',
           'Eko',
           'Shang',
           'RBDKBBI'
           ]

for a in area:
    for b in methods:
        obs500[a + b] = ''
df = obs500.set_index('DATE')

#
# mgr = mp.Manager()
# ns = mgr.Namespace()
# ns.df = obs500


# def main(q):
# def apply_parallel(df, func, args={}):
initime = datetime.now()
l2 = len(datalist)
for q in range(l2):
    file = datalist[q]
    filename = path + file
    nc_file = nc4.Dataset(filename, 'r')
    ddd = file[1:8]
    parsed = datetime.strptime(ddd, '%Y%j')
    day = datetime.strftime(parsed, '%Y%m%d')
    print('processing')
    print(day)

    lon = nc_file.groups['navigation_data'].variables['longitude'][:]
    lat = nc_file.groups['navigation_data'].variables['latitude'][:]
    variables = nc_file.groups['geophysical_data'].variables

    x = np.arange(minlon, maxlon, 0.01)  # 1 km grid,
    y = np.arange(maxlat, minlat, -0.01)

    for i in variables:
        var = variables[i][:]
        np.where(var <= 0, var, np.nan)
        if i != 'l2_flags':
            var_re, grid = gs.swath_resampling(var, lon, lat, x, y, 5000)  # 1 km grid
            # var_re=var_re.filled()
            if np.ma.is_mask(var_re):
                var_re = var_re.filled(np.nan)
                var_re[var_re == -32767.0] = np.nan
            variables[i] = var_re
        else:
            # var_re = var_re.filled()
            variables[i] = var
    nc_file.close()
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

    F0 = [172.912, 187.622, 205.878, 194.933, 185.747, 186.539, 183.869, 157.811, 152.255, 148.052, 128.065]
    nlw412 = Rrs_412 * F0[0]
    nlw443 = Rrs_443 * F0[1]
    nlw469 = Rrs_469 * F0[2]
    nlw488 = Rrs_488 * F0[3]
    nlw531 = Rrs_531 * F0[4]
    nlw547 = Rrs_547 * F0[5]
    nlw555 = Rrs_555 * F0[6]
    nlw645 = Rrs_645 * F0[7]
    nlw667 = Rrs_667 * F0[8]
    nlw678 = Rrs_678 * F0[9]
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
    nflh = variables['nflh']
    chl = variables['chlor_a']


    def rrs_to_qaa(Rrs):
        bbp555QAA = np.zeros(np.shape(chl))
        a443 = np.zeros(np.shape(chl))
        for m in range(len(y)):
            for n in range(len(x)):
                # for k in range(10):
                # if (Rrs[k, i, j] < 0) or (Rrs[k, i, j] == np.nan):
                # bbp[k] = np.nan
                # else:
                bbp, a = gs.QAAv6MODIS(Rrs[:, m, n])
                bbp555QAA[m, n] = bbp[6]
                a443[m, n] = a[1]
        return bbp555QAA, a443


    bbp555QAA, a443 = rrs_to_qaa(Rrs)

    Imagebbp = km.bbpratio(chl, bbp555QAA, nflh, Rrs_555)

    ImageRBDKBBI = km.RBDKBBI(nlw667, nlw678, chl)

    Imageshang = km.shang(Rrs_443, Rrs_488, Rrs_531, Rrs_555, nflh, a443)

    ImageEKO = km.siswanto(x, y, chl, nlw412, nlw443, nlw488, nlw547, nlw555, nlw645, nlw667)

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
    # OR = {'minlon': 134.02,
    #       'maxlon': 135.45,
    #       'minlat': 34.20,
    #       'maxlat': 34.78}
    OB = {'minlon': 134.89,
          'maxlon': 135.48,
          'minlat': 34.27,
          'maxlat': 34.73}  # Osaka bay

    HN = {'minlon': 134.0,
          'maxlon': 134.88,
          'minlat': 34.8,
          'maxlat': 34.20}  # Harima nada
    # areadict = {
    #     'EO': EO,
    #     'NO': NO,
    #     'UC': UC,
    #     'OR': OR
    # }
    # areadict = {
    #     'EO': EO,
    #     'NO': NO,
    #     'UC': UC,
    #     'OB': OB,
    #     'HN':HN
    # }
    areadict = {
        'OB': OB,
        'HN':HN
    }



    methoddict = {'bbp': Imagebbp,
                  'Eko': ImageEKO,
                  'Shang': Imageshang,
                  'RBDKBBI': ImageRBDKBBI
                  }

    for key in areadict:
        area = areadict[key]
        latmin = area['minlat']
        latmax = area['maxlat']
        lonmin = area['minlon']
        lonmax = area['maxlon']
        for key2 in methoddict:
            result = methoddict[key2]
            latindex = np.arange(len(y))
            lonindex = np.arange(len(x))
            latindexmin = np.max(latindex[y >= latmin])
            latindexmax = np.min(latindex[y >= latmax])
            lonindexmin = np.min(lonindex[x >= lonmin])
            lonindexmax = np.min(lonindex[x >= lonmax])

            sub_result = result[latindexmax:latindexmin, lonindexmin:lonindexmax]
            # lat2 = lat[latindexmax:latindexmin]
            # lon2 = lon[lonindexmin:lonindexmax]
            # gs.plot_WT_image(sub_result, lon2, lat2)
            if 3 in sub_result:
                df.loc[day, key + key2] = 1

            if (2 in sub_result) & (3 not in sub_result):
                df.loc[day, key + key2] = 0

            if (1 in sub_result) & (3 not in sub_result):
                df.loc[day, key + key2] = 0

    print('File processed...')
    print(q / l2)
    time2 = datetime.now()
    print('Time consuming...')
    print(time2-initime)

df.to_excel('result_t4_OR_seprated_back.xlsx')

## todo there are some problem in the output
## todo possible multiprocessing code
# import multiprocessing
#
# def apply_parallel(df, func, args={}):
#     """ Multiprocessing apply for Dataframe """
#     cores = multiprocessing.cpu_count()
#     if args: func = partial(func, **args)
#
#     df_split = numpy.array_split(df, cores)
#
#     with multiprocessing.Pool(cores) as pool:
#         results = pool.map(func, df_split)
#         try:
#             df = pandas.concat(results)
#         except ValueError:
#             # result could be a list of Nones
#             pass
#     return df
