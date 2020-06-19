import netCDF4 as nc4
import numpy as np


from geo_Collection import geo_self as gs

import os
import glob
import datetime
import warnings

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
datalist = glob.glob('*.nc')
# 20150715 Uwajima
# 20180723 Osaka


minlat = 32.5
minlon = 131
maxlat = 35
maxlon = 136
# area of full seto-inland sea

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

L = len(datalist)
initime = datetime.datetime.now()

# plt.subplot(2,3,3)
# with concurrent.futures.ProcessPoolExecutor() as executor:
a1 = np.arange(len(datalist))

# 想一下我这个代码要怎么跑


def main(a1):
    file = datalist[a1]
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
            var_re, grid = gs.swath_resampling(var, lon, lat, x, y, 5000)  # 1 km grid
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
    nlw748 = Rrs_748 * F0[10]
    # nlw412 = variables['nLw_412']
    # nlw443 = variables['nLw_443']
    # nlw469 = variables['nLw_469']
    # nlw488 = variables['nLw_488']
    # nlw531 = variables['nLw_531']
    # nlw547 = variables['nLw_547']
    # nlw555 = variables['nLw_555']
    # nlw645 = variables['nLw_645']
    # nlw667 = variables['nLw_667']
    # nlw678 = variables['nLw_678']
    # nlw748 = variables['nLw_748']

    # Rrs_412 = variables['Rrs_412'].filled(np.nan)
    # Rrs_443 = variables['Rrs_443'].filled(np.nan)
    # Rrs_469 = variables['Rrs_469'].filled(np.nan)
    # Rrs_488 = variables['Rrs_488'].filled(np.nan)
    # Rrs_531 = variables['Rrs_531'].filled(np.nan)
    # Rrs_547 = variables['Rrs_547'].filled(np.nan)
    # Rrs_555 = variables['Rrs_555'].filled(np.nan)
    # Rrs_645 = variables['Rrs_645'].filled(np.nan)
    # Rrs_667 = variables['Rrs_667'].filled(np.nan)
    # Rrs_678 = variables['Rrs_678'].filled(np.nan)
    # Rrs_748 = variables['Rrs_748'].filled(np.nan)
    # nlw412 = variables['nLw_412'].filled(np.nan)
    # nlw443 = variables['nLw_443'].filled(np.nan)
    # nlw469 = variables['nLw_469'].filled(np.nan)
    # nlw488 = variables['nLw_488'].filled(np.nan)
    # nlw531 = variables['nLw_531'].filled(np.nan)
    # nlw547 = variables['nLw_547'].filled(np.nan)
    # nlw555 = variables['nLw_555'].filled(np.nan)
    # nlw645 = variables['nLw_645'].filled(np.nan)
    # nlw667 = variables['nLw_667'].filled(np.nan)
    # nlw678 = variables['nLw_678'].filled(np.nan)
    # nlw748 = variables['nLw_748'].filled(np.nan)
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
    nflh = nlw678 - (70 / 81) * nlw667 - (11 / 81) * nlw748
    chl = variables['chlor_a']
    # gs.plot_geo_image(nflh, lons, lats, log10=False,
    #                   title='nflh' + '\n' + nc_file.time_coverage_end[0:13],
    #                   save_image='nflh' + nc_file.time_coverage_end[0:13])
    # gs.plot_geo_image(nlw678, lons, lats, log10=False,
    #                   title='nlw678' + '\n' + nc_file.time_coverage_end[0:13],
    #                   save_image='nlw678' + nc_file.time_coverage_end[0:13])
    # gs.plot_geo_image(nlw667, lons, lats, log10=False,
    #                   title='nlw667' + '\n' + nc_file.time_coverage_end[0:13],
    #                   save_image='nlw667' + nc_file.time_coverage_end[0:13])
    # gs.plot_geo_image(nlw748, lons, lats, log10=False,
    #                   title='nlw748' + '\n' + nc_file.time_coverage_end[0:13],
    #                   save_image='nlw748' + nc_file.time_coverage_end[0:13])

    nlw443_412slope = (nlw443 - nlw412) / (443 - 412)
    nlw488_443slope = (nlw488 - nlw443) / (488 - 443)
    nlw547_488slope = (nlw547 - nlw488) / (547 - 488)
    nlw547_412slope = (nlw547 - nlw412) / (547 - 412)
    nlw488_412slope = (nlw488 - nlw412) / (488 - 412)
    ImageEKO = np.ones(np.shape(chl))
    # TODO There must be some problem with EKO's method!
    ImageEKO[chl.mask] = 0
    for i in range(len(y)):
        for j in range(len(x)):
            if nlw547[i, j] > max(nlw412[i, j], nlw443[i, j], nlw488[i, j], nlw555[i, j], nlw645[i, j], nlw667[i, j]):
                if (not ((nlw443_412slope[i, j] > nlw488_443slope[i, j]) & (
                        nlw488_443slope[i, j] > nlw547_488slope[i, j]))):
                    if (not (nlw547_412slope[i, j]) > -0.0003 * (
                            (np.log(chl[i, j])) ** 2) + 0.0024 * np.log(chl[i, j]) - 0.00005):
                        if not (nlw547[i, j] > 0.8):
                            if not (abs(nlw547_488slope[i, j] - nlw488_412slope[i, j]) > 0.006):
                                ImageEKO[i, j] = 3
                            else:
                                ImageEKO[i, j] = 2

    gs.plot_WT_image(ImageEKO, lons, lats,
                      title='Result by SSD Siswanto et al 2013' + nc_file.time_coverage_end[0:13],
                      save_image='Result by SSD Siswanto et al 2013 ' + nc_file.time_coverage_end[0:13])

    def bbpratio(chl, bbp555QAA, bbpmorel):
        # 0 nodata
        # 1 other water
        # 2 other phytoplankton
        # 3k,bloom

        Imagebbp = np.ones(np.shape(chl))

        Imagebbp[chl.mask] = 0
        Imagebbp[(chl > 1.5) & (bbp555QAA < bbpmorel)] = 3  # Karenia
        return Imagebbp

    def nlwratio(chl, bbpmorel, nlw555):
        Imagenlw = np.ones(np.shape(chl))
        # #TODO 有待商议
        Imagenlw[chl.mask] = 0
        Imagenlw[(chl > 1.5) & (nlw555 < bbpmorel)] = 3
        return Imagenlw

    def SS(nlw488, nlw443, nlw531):
        SS_90 = nlw488 - nlw443 - (nlw531 - nlw443) * ((488 - 443) / (531 - 443))
        ImageSS = np.ones(np.shape(chl))
        ImageSS[chl.mask] = 0
        ImageSS[SS_90 < 0] = 3
        return ImageSS

    def RBDKBBI(RBD, KBBI):
        ImageRBD = np.ones(np.shape(chl))
        ImageKBBI = np.ones(np.shape(chl))
        ImageRBD[chl.mask] = 0
        ImageKBBI[chl.mask] = 0
        ImageRBD[RBD > 0.015] = 3
        ImageKBBI[RBD > 0.015] = 2
        ImageKBBI[(RBD > 0.015) & (KBBI > 0.3 * RBD)] = 3
        return ImageRBD, ImageKBBI

    # plt.savefig('subplot' + nc_file.time_coverage_end[0:13])
    def shang(Rrs_443, Rrs_488, Rrs_531, Rrs_555, nflh, a443):
        BI = ((Rrs_488 - Rrs_443) / (488 - 443)) / ((Rrs_555 - Rrs_531) / (555 - 531))
        ImageBI = np.ones(np.shape(Rrs_443))

        # Finally, a criterion based on the combination of FLH, a(443), and BI is proposed below:
        # When FLH is doubled over the background level and a(443)》>=0.5 m21,
        # if 0.0 < BI <= 0.3, it suggests a dinoflagellate bloom;
        # if 0.3 < BI <= 1.0, it suggests a diatom bloom.
        # don`t have standard nflh product in the very coastal area
        ImageBI[nflh.mask] = 0
        ImageBI[(nflh > 0.02) & (a443 >= 0.5) & ((BI > 0.3) & (BI <= 1))] = 2
        ImageBI[(nflh > 0.02) & (a443 >= 0.5) & ((BI > 0.0) & (BI <= 0.3))] = 3
        return ImageBI, BI

    def shen(Rrs_555, Rrs_667, Rrs_748):
        RDI = (1 / Rrs_667 - 1 / Rrs_555) * Rrs_748
        ImageShen = np.ones(np.shape(Rrs_555))
        ImageShen[RDI > 0.16] = 3
        # R_slope = np.arctan(100 * (1 - ((Rrs_555 / Rrs_667) / (555 - 667))))
        # ImageShen[(RDI > 0.14) & (R_slope < 0.5)] = 5
        return ImageShen

    def rrs_to_qaa(Rrs):
        bbp555QAA = np.zeros(np.shape(chl))
        a443 = np.zeros(np.shape(chl))
        for m in range(len(y)):
            for n in range(len(x)):
                # for k in range(10):
                # if (Rrs[k, i, j] < 0) or (Rrs[k, i, j] == np.nan):
                # bbp[k] = np.nan
                # else:
                bbp, a = QAAv6(Rrs[:, m, n])
                bbp555QAA[m, n] = bbp[6]
                a443[m, n] = a[1]
        return bbp555QAA, a443


    bbp555QAA, a443 = rrs_to_qaa(Rrs)
    # ImageShen = shen(Rrs_555, Rrs_667, Rrs_748, )
    # gs.plot_geo_image(ImageShen, lons, lats, log10=False,
    #                   title='Result by Shen et al(2019) near Osaka' + '\n' + nc_file.time_coverage_end[0:13], caxis=[0, 6],
    #                   save_image='Result by shen et al(2019) near Osaka' + nc_file.time_coverage_end[0:13])
    ImageBI, BI = shang(Rrs_443, Rrs_488, Rrs_531, Rrs_555, nflh, a443)
    gs.plot_WT_image(ImageBI, lons, lats,
                      title='Result by BI Shang et al 2014 ' + '\n' + nc_file.time_coverage_end[0:13],
                      save_image='Result by BI shang et al 2014 ' + nc_file.time_coverage_end[0:13])
    # gs.plot_geo_image(ImageBI, lons, lats, log10=False,
    #                   title='Result by BI Shang et al 2014 near Osaka' + '\n' + nc_file.time_coverage_end[0:13],
    #                   caxis=[0, 6],
    #                   save_image='Result by BI shang et al 2014 near Osaka' + nc_file.time_coverage_end[0:13])
    # bbpmorel = 0.3 * (chl ** 0.62) * (0.002 + 0.02 * (0.5 - 0.25 * np.log10(chl)))

    #
    bbpmorel = 0.3 * (chl ** 0.62) * (0.002 + 0.02 * (0.5 - 0.25 * np.log10(chl)))
    Imagebbp = bbpratio(chl, bbp555QAA, bbpmorel)
    gs.plot_WT_image(Imagebbp, lons, lats,
                      title='Result by bbp ratio Cannizzaro et al 2009 ' + '\n' + nc_file.time_coverage_end[0:13],
                      save_image='Result by bbp ratio Cannizzaro et al 2009 ' + nc_file.time_coverage_end[0:13])
    # Imagenlw = nlwratio(chl, bbpmorel, nlw555)
    # gs.plot_geo_image(Imagenlw, lons, lats, log10=False,
    #                   title='Result by nlw ratio near Osaka' + nc_file.time_coverage_end[0:13],
    #                   caxis=[0, 6], save_image='Result by nlw ratio near Osaka' + nc_file.time_coverage_end[0:13])
    # ImageSS = SS(nlw488, nlw443, nlw531)
    # gs.plot_geo_image(ImageSS, lons, lats, log10=False,
    #                   title='Result by ss490 near Osaka' + '\n' + nc_file.time_coverage_end[0:13], caxis=[0, 6],
    #                   save_image='Result by ss490 near Osaka' + nc_file.time_coverage_end[0:13])
    RBD = nlw678 - nlw667
    KBBI = (nlw678 - nlw667) / (nlw678 + nlw667)
    ImageRBD, ImageKBBI = RBDKBBI(RBD, KBBI, )
    # # plt.subplot(2,3,4)
    gs.plot_WT_image(ImageKBBI, lons, lats,
                      title='Result by RBD-KBBI Amin et al 2009 ' + '\n' + nc_file.time_coverage_end[0:13],
                      save_image='Result by RBD-KBBI Amin et al 2009 ' + nc_file.time_coverage_end[0:13])


    nlw443_412slope = (nlw443 - nlw412) / (443 - 412)
    nlw488_443slope = (nlw488 - nlw443) / (488 - 443)
    nlw547_488slope = (nlw547 - nlw488) / (547 - 488)
    nlw547_412slope = (nlw547 - nlw412) / (547 - 412)
    nlw488_412slope = (nlw488 - nlw412) / (488 - 412)


if __name__ == '__main__':
    mp.set_start_method("spawn")
    pool = mp.Pool(processes=8)
    pool.map(main, a1)

end7 = datetime.datetime.now()
print(end7 - initime)
# todo add nfl and rrs to bbp ratio
# todo specify the region
# todo evaluation script
# todo step by step performance
# todo enhancement?
# todo atomospher correction influence?