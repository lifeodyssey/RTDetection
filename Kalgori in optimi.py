import netCDF4 as nc4
import numpy as np
from collections import OrderedDict
from geo_Collection import geo_web as gs
from QAAV6 import QAAv6
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
os.chdir('/Users/zhenjia/Desktop/Project/Seto/MODIS/l1/l2_lac')
path = '/Users/zhenjia/Desktop/Project/Seto/MODIS/l1/l2_lac/'
datalist = glob.glob('20080727*.nc')
minlat = 32.5
minlon = 130.5
maxlat = 35
maxlon = 136
# area of full seto-inland sea
# for a1 in range(1):
L = len(datalist)
initime = datetime.datetime.now()

# plt.subplot(2,3,3)
# with concurrent.futures.ProcessPoolExecutor() as executor:
a1 = np.arange(len(datalist))


def total(a1):
    start = datetime.datetime.now()

    def rrs_to_qaa(Rrs):
        bbp555QAA = np.zeros(np.shape(lats))
        a443 = np.zeros(np.shape(lats))
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

    def bbpratio(lats, chl, bbp555QAA, bbpmorel):

        Imagebbp = np.zeros(np.shape(lats))

        Imagebbp[(chl > 1.5) & (bbp555QAA < bbpmorel)] = 5  # Karenia
        return Imagebbp

    def nlwratio(lats, chl, bbpmorel, nlw555):
        Imagenlw = np.zeros(np.shape(lats))
        # #TODO 有待商议
        Imagenlw[(chl > 1.5) & (nlw555 < bbpmorel)] = 5
        return Imagenlw

    def SS(nlw488, nlw443, nlw531, lats):

        SS_90 = nlw488 - nlw443 - (nlw531 - nlw443) * ((488 - 443) / (531 - 443))
        ImageSS = np.zeros(np.shape(lats))
        ImageSS[SS_90 < 0] = 5
        return ImageSS

    def RBDKBBI(RBD, KBBI, lats):
        ImageRBD = np.zeros(np.shape(lats))
        ImageKBBI = np.zeros(np.shape(lats))
        ImageRBD[RBD > 0.015] = 5
        ImageKBBI[(RBD > 0.015) & (KBBI > 0.3 * RBD)] = 5
        return ImageRBD, ImageKBBI

    def EKO(nlw412, nlw443, nlw488, nlw555, nlw645, nlw667, nlw547, chl, lats, x, y):
        nlw443_412slope = (nlw443 - nlw412) / (443 - 412)
        nlw488_443slope = (nlw488 - nlw443) / (488 - 443)
        nlw547_488slope = (nlw547 - nlw488) / (547 - 488)
        nlw547_412slope = (nlw547 - nlw412) / (547 - 412)
        nlw488_412slope = (nlw488 - nlw412) / (488 - 412)
        ImageEKO0 = np.zeros(np.shape(lats))
        # # TODO There must be some problem with EKO's method!
        # for i in range(len(y)):
        #     for j in range(len(x)):
        #
        #         if nlw547[i, j] > max(nlw412[i, j], nlw443[i, j], nlw488[i, j], nlw555[i, j], nlw645[i, j],
        #                               nlw667[i, j]):
        #             if not ((nlw443_412slope[i, j] > nlw488_443slope[i, j]) &
        #                     (nlw488_443slope[i, j] > nlw547_488slope[i, j])):
        #                 if not nlw547_412slope[i, j] > -0.0003 * (np.square(np.log(chl[i, j]))) + 0.0024 * np.log(
        #                         chl[i, j]) - 0.00005:
        #                     if nlw547[i, j] > 0.8:
        #                         ImageEKO[i, j] = 3
        #                     else:
        #                         if abs(nlw547_488slope[i, j] - nlw488_412slope[i, j]) > 0.006:
        #                             ImageEKO[i, j] = 3
        #                         else:
        #                             ImageEKO[i, j] = 5
        # for i in range(len(y)):
        #     for j in range(len(x)):
        #         if nlw547[i, j] > max(nlw412[i, j], nlw443[i, j], nlw488[i, j], nlw555[i, j], nlw645[i, j],
        #                               nlw667[i, j]):
        #             if (not ((nlw443_412slope[i, j] > nlw488_443slope[i, j]) & (
        #                     nlw488_443slope[i, j] > nlw547_488slope[i, j]))):
        #                 if (not (nlw547_412slope[i, j]) > -0.0003 * (
        #                         (np.log(chl[i, j])) ** 2) + 0.0024 * np.log(chl[i, j]) - 0.00005):
        #                     if nlw547[i, j] > 0.8:
        #                         ImageEKO[i, j] = 3
        #                     else:
        #                         if abs(nlw547_488slope[i, j] - nlw488_412slope[i, j]) > 0.006:
        #                             ImageEKO[i, j] = 3
        #                         else:
        #                             ImageEKO[i, j] = 5
        #
        # return ImageEKO
        for bi in range(len(y)):
            for j in range(len(x)):
                if nlw547[bi, j] > max(nlw412[bi, j], nlw443[bi, j], nlw488[bi, j], nlw555[bi, j], nlw645[bi, j],
                                       nlw667[bi, j]):
                    if (not ((nlw443_412slope[bi, j] > nlw488_443slope[bi, j]) & (
                            nlw488_443slope[bi, j] > nlw547_488slope[bi, j]))):
                        if (not (nlw547_412slope[bi, j]) > -0.0003 * (
                                (np.log(chl[bi, j])) ** 2) + 0.0024 * np.log(chl[bi, j]) - 0.00005):
                            if not (nlw547[bi, j] > 0.8):
                                if not (abs(nlw547_488slope[bi, j] - nlw488_412slope[bi, j]) > 0.006):
                                    ImageEKO0[bi, j] = 5

        return ImageEKO0

    def shang(Rrs_443, Rrs_488, Rrs_531, Rrs_555, nflh, a443):
        BI = ((Rrs_488 - Rrs_443) / (488 - 443)) / ((Rrs_555 - Rrs_531) / (555 - 531))
        ImageBI = np.zeros(np.shape(lats))

        # Finally, a criterion based on the combination of FLH, a(443), and BI is proposed below:
        # When FLH is doubled over the background level and a(443)》>=0.5 m21,
        # if 0.0 < BI <= 0.3, it suggests a dinoflagellate bloom;
        # if 0.3 < BI <= 1.0, it suggests a diatom bloom.
        # don`t have standard nflh product in the very coastal area
        ImageBI[(nflh > 0.02) & (a443 >= 0.5) & ((BI > 0.3) & (BI <= 1))] = 3
        ImageBI[(nflh > 0.02) & (a443 >= 0.5) & ((BI > 0.0) & (BI <= 0.3))] = 5
        return ImageBI

    def shen(Rrs_555, Rrs_667, Rrs_748, lats):
        RDI = (1 / Rrs_667 - 1 / Rrs_555) * Rrs_748
        ImageShen = np.zeros(np.shape(lats))
        ImageShen[RDI > 0.16] = 5
        # R_slope = np.arctan(100 * (1 - ((Rrs_555 / Rrs_667) / (555 - 667))))
        # ImageShen[(RDI > 0.14) & (R_slope < 0.5)] = 5
        return ImageShen

    file = datalist[a1]
    filename = path + file
    # print(filename)
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
                 var_re.mask = np.ma.nomask
                 var_re[var_re == -32767.0] = np.nan
            variables[i] = var_re
        else:
            # var_re = var_re.filled()
            variables[i] = var

    lons = grid.lons
    lats = grid.lats
    end1 = datetime.datetime.now()
    # print('Resample')
    # print(end1 - start)

    # Lambda=np.array([412, 443, 469,488,531,547,555,645,667,678])
    # Rrs_412=variables['Rrs_412'].filled()
    # Rrs_443 = variables['Rrs_443'].filled()
    # Rrs_469 = variables['Rrs_469'].filled()
    # Rrs_488 = variables['Rrs_488'].filled()
    # Rrs_531 = variables['Rrs_531'].filled()
    # Rrs_547 = variables['Rrs_547'].filled()
    # Rrs_555 = variables['Rrs_555'].filled()
    # Rrs_645 = variables['Rrs_645'].filled()
    # Rrs_667 = variables['Rrs_667'].filled()
    # Rrs_678= variables['Rrs_678'].filled()
    # nflh=variables['nflh'].filled()
    # chl=variables['chlor_a'].filled()
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
    nflh = nlw678 - (70 / 81) * nlw667 - (11 / 81) * nlw748
    chl = variables['chlor_a']
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

    # bbpmorel = 0.3 * (chl ** 0.62) * (0.002 + 0.02 * (0.5 - 0.25 * np.log10(chl)))
    bbp555QAA, a443 = rrs_to_qaa(Rrs)
    #
    # Imagebbp=bbpratio(lats,chl,bbp555QAA,bbpmorel)
    # gs.plot_geo_image(Imagebbp, lons, lats, log10=False, title='Result by bbpratio' + nc_file.time_coverage_end[0:13],
    #                   caxis=[0, 6],
    #                   save_image='Result by bbp ratio' + nc_file.time_coverage_end[0:13])
    # end2 = datetime.datetime.now()
    # print('bbp ratio')
    # print(end2 - end1)
    #
    # #nLw ratio
    # #Estimates the relationship between water leaving radiance (Lw(551)) and CHL.
    # # If Lw(551) is lower than the backscattering calculated at λ = 550 nm (bbp MOREL(550); Morel, 1988)as function of CHL then the pixel is classified as a K. brevis bloom.
    # F0 = [172.912, 187.622, 205.878, 194.933, 185.747, 186.539, 183.869, 157.811, 152.255, 148.052]
    # nlw555 = Rrs_555 * F0[6]
    # Imagenlw=nlwratio(lats,chl,bbpmorel,nlw555)
    # gs.plot_geo_image(Imagenlw, lons, lats, log10=False,
    #                   title='Result by nlw ratio' + nc_file.time_coverage_end[0:13],
    #                   caxis=[0, 6],save_image='Result by nlw ratio' + nc_file.time_coverage_end[0:13])
    # end3 = datetime.datetime.now()
    # print('nlw ratio')
    # print(end3 - end2)
    # todo test eko indivisually
    # todo find chla thresholds
    # todo using seadas derive all of the band as well as rayleigh corrected reflectance
    # todo write shen,shang
    # todo optimize code
    # todo find other paper related to this

    #
    #
    # #Spectral shape at 490 nm (SS_490)
    # #Calculates the spectral shape at 490 nm (SS_490)
    # # using nLw: SS_90 = nLw(488) − nLw(443) − (nLw(531) − nLw(443)) × ((488 − 443) ÷ (531 − 443))
    # # A negative SS_490 is indicative of a K. brevis bloom.
    # nlw488 = Rrs_488 * F0[3]
    # nlw443 = Rrs_443 * F0[1]
    # nlw531 = Rrs_531 * F0[4]
    # ImageSS=SS(nlw488,nlw443,nlw531,lats)
    # gs.plot_geo_image(ImageSS, lons, lats, log10=False,
    #                   title='Result by ss490' + '\n' + nc_file.time_coverage_end[0:13],caxis=[0, 6],
    #                   save_image='Result by ss490' + nc_file.time_coverage_end[0:13])

    # end4 = datetime.datetime.now()
    # print('SS')
    # print(end4 - end3)
    # gs.plot_geo_image(ImageEKO,lonmlat)
    # #RBD
    # #Takes advantage of the high fluorescence properties of dinoflagellate blooms.
    # # RBD = nLw(678) − nLw(667) RBD N 0.15 Wm−2 μm−1 sr−1
    # # indicates a dinoflagellate bloom but not specifically a K. brevis bloom.
    # nlw667 = Rrs_667 * F0[8]
    # nlw678 = Rrs_678 * F0[9]
    # RBD = nlw678 - nlw667
    # KBBI = (nlw678 - nlw667) / (nlw678 + nlw667)
    # ImageRBD,ImageKBBI=RBDKBBI(RBD,KBBI,lats)
    # # plt.subplot(2,3,4)
    # gs.plot_geo_image(ImageRBD, lons, lats, log10=False,
    #                   title='Result by RBD' + '\n' + nc_file.time_coverage_end[0:13],caxis=[0, 6],
    #                   save_image='Result by RBD' + nc_file.time_coverage_end[0:13])
    # end5 = datetime.datetime.now()

    #
    # #RBD–KBBI
    # #Combines the high fluorescence and low backscattering properties of K. brevis blooms in 2 main equations:
    # # RBD = nLw(678) − nLw(667) and KBBI = (nLw(678) − nLw(667)) / (nLw(678) + nLw(667))
    # # A K. brevis bloom should meet the following criteria:
    # # RBD N 0.15 Wm−2 μm−1 sr−1,and KBBI N 0.3*RBD
    # KBBI = (nlw678 - nlw667) / (nlw678 + nlw667)
    #
    # gs.plot_geo_image(ImageKBBI, lons, lats, log10=False,
    #                   title='Result by RBD-KBBI' + '\n' + nc_file.time_coverage_end[0:13],caxis=[0, 6],
    #                   save_image='Result by RBD-KBBI' + nc_file.time_coverage_end[0:13])
    # end6 = datetime.datetime.now()
    # print('RBD-KBBI')
    # print(end6 - end5)
    # EKO
    # nlw Peak at 547
    # NOT！nlw443-412 slope>nlw488-443slope & nlw488-443>nlw547-488slope
    # NOT! nlw547-412 slope>-0.0003[ln(chla)}^2+0.0024[ln(chla)]-0.00005
    # NOT! nlw547>0.8
    # NOT! nlw488-412slope-nlw547-488slope>0.06
    # Lambda=np.array([412, 443, 469,488,531,547,555,645,667,678])
    # The most complex one
    # How much is chla affected?
    # nlw412 = Rrs_412 * F0[0]
    # nlw469 = Rrs_469 * F0[2]
    # nlw547 = Rrs_547 * F0[5]
    # nlw645 = Rrs_645 * F0[7]

    # gs.plot_geo_image(ImageEKO, lons, lats, log10=False,
    # title='EKO1',
    # caxis=[0.1, 6])

    #     if(not(nlw547_412slope[i,j])>-0.0003*((np.log(chl[i,j]))**2+0.0024*np.log(chl[i,j]))-0.00005):
    #         if(not(nlw547[i,j]>0.8)):
    #             if(not(abs(nlw547_488slope[i,j]-nlw488_412slope[i,j])>0.006)):
    #                 ImageEKO[i,j]=5
    # plt.subplot(2,3,6)
    ImageEKO = EKO(nlw412, nlw443, nlw488, nlw547, nlw555, nlw645, nlw667, chl, lats, x, y)
    gs.plot_geo_image(ImageEKO, lons, lats, log10=False,
                      title='Result by Eko Siswanto et al(2013)' + '\n' + nc_file.time_coverage_end[0:13], caxis=[0, 6],
                      save_image='Result by Eko Siswanto et al(2013)' + nc_file.time_coverage_end[0:13])
    # plt.savefig('subplot' + nc_file.time_coverage_end[0:13])
    end7 = datetime.datetime.now()
    print('EKO')
    print(end7 - start)

    # nFLH

    # shang 2014,diatom and dino
    # BI =R488-R443/488-443
    # /r555-r531/555-531
    # ImageBI=shang(Rrs_443,Rrs_488,Rrs_531,Rrs_555,nflh,a443)
    # gs.plot_geo_image(ImageBI, lons, lats, log10=False,
    #                   title='Result by Shang et al(2014)' + '\n' + nc_file.time_coverage_end[0:13],caxis=[0, 6],
    #                   save_image='Result by shang et al(2013)' + nc_file.time_coverage_end[0:13])
    # end8 = datetime.datetime.now()
    # print('SHANG')
    # print(end8 - end7)
    # 　shen 2019
    # 　red tide detection index:rdi=(1/Rrs1-1/Rrs2)*Rrs3
    # (a) the λ1, λ2, and λ3 of Eq. (1) were set using the MERIS bands of 665, 560, and 753 nm,
    # (b) the λ1, λ2, and λ3 were set using the MODIS bands of 667, 555, and 748 nm, and
    # (c) the λ1, λ2, and λ3 were set using the GOCI bands of 660, 555, and 745 nm.
    # 　0.14~0.16
    # R_slope =tan−1 (100×(1 – (Rrs(λ2)/ Rrs(λ1))) / (λ2–λ1))
    # λ1 represents the baseline wavelength in the range of 550–570 nm, and λ2 is proposed to be a wavelength in the range of 570–670 nm.
    # rrs slope 0.5 original paper do not included
    # ImageEKO=np.zeros(np.shape(lats))
    # RDI=((1/Rrs_667)-(1/Rrs_555))*RRS_7
    # ImageShen=shen(Rrs_555,Rrs_667,Rrs_748,lats)
    # gs.plot_geo_image(ImageShen, lons, lats, log10=False,
    #                   title='Result by Shen et al(2019)' + '\n' + nc_file.time_coverage_end[0:13],caxis=[0, 6],
    #                   save_image='Result by shen et al(2019)' + nc_file.time_coverage_end[0:13])
    # end9 = datetime.datetime.now()
    # print('shen')
    # print(end9-end8)
    #
    # print('one loop')
    # print(end9 - start)

    print(print('percent: {:.2%}'.format((a1 + 1) / L)))


pool = mp.Pool(processes=8)
pool.map(total, a1)
endt = datetime.datetime.now()
print('all time')
print(endt - initime)
# 2008-07-26T03:54:59.771Z
# Resample
# 0:00:02.266570
# Lat: [32.510, 35.000] | Lon: [130.500, 135.990] | SDS: [0.000, 5.000]
# bbp ratio
# 0:03:18.745898
# Lat: [32.510, 35.000] | Lon: [130.500, 135.990] | SDS: [0.000, 5.000]
# nlw ration
# 0:00:38.067859
# Lat: [32.510, 35.000] | Lon: [130.500, 135.990] | SDS: [0.000, 5.000]
# SS
# 0:00:37.835756
# Lat: [32.510, 35.000] | Lon: [130.500, 135.990] | SDS: [0.000, 5.000]
# RBD
# 0:00:38.416747
# Lat: [32.510, 35.000] | Lon: [130.500, 135.990] | SDS: [0.000, 5.000]
# RBD-KBBI
# 0:00:39.070265
# Lat: [32.510, 35.000] | Lon: [130.500, 135.990] | SDS: [0.000, 0.000]
# EKO
# 0:00:43.163418
