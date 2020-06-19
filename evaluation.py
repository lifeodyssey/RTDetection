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

import ctypes

warnings.filterwarnings("ignore")  # 用来去掉Warning,
os.environ['PROJ_LIB'] = '/Users/zhenjia/opt/anaconda3/share/proj'

obs500 = pd.read_excel('/Users/zhenjia/Desktop/Field_zheng.xlsx', sheet_name='500')
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

area = [
    'EO',
    'NO',
    'UC',
    'OR']

methods = ['bbp',
           'Eko',
           'Shang',
           'RBDKBBI'
           ]

for a in area:
    for b in methods:
        obs500[a + b] = ''
obs500.set_index('DATE')

#
# mgr = mp.Manager()
# ns = mgr.Namespace()
# ns.df = obs500
df = obs500

# def main(q):
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
    OR = {'minlon': 134.02,
          'maxlon': 135.45,
          'minlat': 34.20,
          'maxlat': 34.78}
    areadict = {
        'EO': EO,
        'NO': NO,
        'UC': UC,
        'OR': OR
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
            latindexmin = np.max(latindex[y <= latmin])
            latindexmax = np.min(latindex[y >= latmax])
            lonindexmin = np.min(lonindex[x >= lonmin])
            lonindexmax = np.min(lonindex[x >= lonmax])

            sub_result = result[latindexmin:latindexmax, lonindexmin:lonindexmax]
            # lat2 = lat[latindexmax:latindexmin]
            # lon2 = lon[lonindexmin:lonindexmax]
            # gs.plot_WT_image(sub_result, lon2, lat2)
            if 3 in sub_result:
                df.loc[day, key + key2] = 1
                print(day)
                print(key+key2)
                print(1)
            if 2 in sub_result:
                df.loc[day, key + key2] = 0
                print(day)
                print(key+key2)
                print(0)
            if 1 in sub_result:
                df.loc[day, key + key2] = 0
                print(day)
                print(key+key2)
                print(0)

    print('File processed...')
    print(q / l2)

df.to_excel('result_t.xlsx')

# if __name__ == '__main__':
#
#     l2 = np.arange(len(datalist))
#     pool = mp.Pool(processes=10)
#     pool.map(main,l2)
#     pool.close()
#     pool.join()
#     ns.df.to_excel('result_t.xlsx')

#
# def perf_measure(y_actual, y_hat):
#     TP = 0
#     FP = 0
#     TN = 0
#     FN = 0
#
#     for i in range(len(y_hat)):
#         if y_actual[i] == y_hat[i] == 1:
#             TP += 1
#         if y_hat[i] == 1 and y_actual[i] != y_hat[i]:
#             FP += 1
#         if y_actual[i] == y_hat[i] == 0:
#             TN += 1
#         if y_hat[i] == 0 and y_actual[i] != y_hat[i]:
#             FN += 1
#
#     return TP, FP, TN, FN
#
#
# def perf(TP, FP, TN, FN):
#     # Sensitivity, hit rate, recall, or true positive rate
#     try:
#         TPR = TP / (TP + FN)
#     except:
#         TPR = 0
#     # # Specificity or true negative rate
#     # TNR = TN / (TN + FP)
#     # # Precision or positive predictive value
#     PPV = TP / (TP + FP)
#     # # Negative predictive value
#     # NPV = TN / (TN + FN)
#     # # Fall out or false positive rate
#     # FPR = FP / (FP + TN)
#     # # False negative rate
#     # FNR = FN / (TP + FN)
#     # # False discovery rate
#     # FDR = FP / (TP + FP)
#
#     # Overall accuracy
#     ACC = (TP + TN) / (TP + FP + FN + TN)
#     try:
#         f1 = (2 * PPV * TPR) / (PPV + TPR)
#     except:
#         f1 = 0
#     return ACC
#
#
# def perf_by_region(sheetname):
#     sheet = pd.read_excel('/Users/zhenjia/Desktop/Project/Seto/result 6/Field.xlsx', sheet_name=sheetname)
#     Oitaest = perf_measure(sheet['East of Oita'], sheet['East of Oita.1'])
#     Oitanorth = perf_measure(sheet['North of Oita'], sheet['North of Oita.1'])
#     Uwajimaco = perf_measure(sheet['Uwajima coast'], sheet['Uwajima coast.1'])
#     Osaka = perf_measure(sheet['Osaka'], sheet['Osaka.1'])
#
#     P_EO = perf(Oitaest[0], Oitaest[1], Oitaest[2], Oitaest[3])
#     P_NO = perf(Oitanorth[0], Oitanorth[1], Oitanorth[2], Oitanorth[3])
#     P_UC = perf(Uwajimaco[0], Uwajimaco[1], Uwajimaco[2], Uwajimaco[3])
#     P_OS = perf(Osaka[0], Osaka[1], Osaka[2], Osaka[3])
#     return P_EO, P_NO, P_UC, P_OS
#
#
# list = ['sisiwanto', 'rbd-kbbi', 'shang_need', 'bbp ratio']
# l = len(list)
# for i in range(l):
#     print(list[i])
#     print('P_EO, P_NO, P_UC,P_OS')
#     P_EO, P_NO, P_UC, P_OS = perf_by_region(list[i])
#     print(P_EO, P_NO, P_UC, P_OS)
