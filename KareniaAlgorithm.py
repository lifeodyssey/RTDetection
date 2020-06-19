import netCDF4 as nc4
import numpy as np
from collections  import OrderedDict
from geo_Collection import geo_web as gs
from QAAV6 import  QAAv6
import os

os.environ['PROJ_LIB'] = '/Users/zhenjia/opt/anaconda3/share/proj'

import glob
import warnings
warnings.filterwarnings("ignore")#用来去arning,
import pandas as pd
import matplotlib.pyplot as plt
os.chdir('/Users/zhenjia/Desktop/Project/Seto/MODIS')
path='/Users/zhenjia/Desktop/Project/Seto/MODIS/'
datalist=glob.glob('20180725*.nc')
# minlat = 32.5
# minlon = 130.5
# maxlat = 35
# maxlon = 136
#Full Seto-Inland Sea
# minlat = 32.65
# minlon = 131.39
# maxlat = 33.48
# maxlon = 132.77
#UWAJIMA
# minlat = 34.18
# minlon = 134.00
# maxlat = 34.88
# maxlon = 135.59
#OSAKA
#Specific for 20150715
minlat = 33.07428
minlon = 131.927472
maxlat = 33.473592
maxlon = 132.625198
L=len(datalist)
#area of full seto-inland sea
#for j in range(len(datalist)):
for a1 in range(L):
    file=datalist[a1]
    filename = path + file
    #print(filename)
    nc_file = nc4.Dataset(filename, 'r')
    print(nc_file.time_coverage_end)
    lon=nc_file.groups['navigation_data'].variables['longitude'][:]
    lat=nc_file.groups['navigation_data'].variables['latitude'][:]
    variables=nc_file.groups['geophysical_data'].variables

    x = np.arange(minlon, maxlon, 0.01) # 1 km grid,
    y = np.arange(maxlat, minlat, -0.01)
    for i in variables:
        var=variables[i][:]
        np.where(var<=0,var,np.nan)
        if i!='l2_flags':
            var_re, grid = gs.swath_resampling(var, lon, lat, x, y, 3000)  # 1 km grid
            #var_re=var_re.filled()
            variables[i] = var_re
        else:
            #var_re = var_re.filled()
            variables[i]=var

    lons=grid.lons
    lats=grid.lats
    #lons=x
    #lats=y
    #Lambda=np.array([412, 443, 469,488,531,547,555,645,667,678])
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
    Rrs_412=variables['Rrs_412']
    Rrs_443 = variables['Rrs_443']
    Rrs_469 = variables['Rrs_469']
    Rrs_488 = variables['Rrs_488']
    Rrs_531 = variables['Rrs_531']
    Rrs_547 = variables['Rrs_547']
    Rrs_555 = variables['Rrs_555']
    Rrs_645 = variables['Rrs_645']
    Rrs_667 = variables['Rrs_667']
    Rrs_678= variables['Rrs_678']
    #nflh=variables['nflh']
    chl=variables['chlor_a']
    # chl[chl<0]=np.nan
    # Rrs_412[Rrs_412<0]=np.nan
    # Rrs_443[Rrs_443 < 0] = np.nan
    # Rrs_469[Rrs_469 < 0] = np.nan
    # Rrs_488[Rrs_488 < 0] = np.nan
    # Rrs_531[Rrs_531 < 0] = np.nan
    # Rrs_547[Rrs_547 < 0] = np.nan
    # Rrs_555[Rrs_555 < 0] = np.nan
    # Rrs_645[Rrs_645 < 0] = np.nan
    # Rrs_667[Rrs_678 < 0] = np.nan
    # Rrs_678[Rrs_678 < 0] = np.nan
    gs.plot_geo_image(Rrs_412, lons, lats, log10=False, title='Rrs412' + nc_file.time_coverage_end[0:13],
                      save_image='Rrs412' + nc_file.time_coverage_end[0:13])
    gs.plot_geo_image(Rrs_443, lons, lats, log10=False, title='Rrs443' + nc_file.time_coverage_end[0:13],
                      save_image='Rrs443' + nc_file.time_coverage_end[0:13])
    gs.plot_geo_image(Rrs_469, lons, lats, log10=False, title='Rrs469' + nc_file.time_coverage_end[0:13],
                      save_image='Rrs469' + nc_file.time_coverage_end[0:13])
    gs.plot_geo_image(Rrs_488, lons, lats, log10=False, title='Rrs488' + nc_file.time_coverage_end[0:13],
                      save_image='Rrs488' + nc_file.time_coverage_end[0:13])
    gs.plot_geo_image(Rrs_531, lons, lats, log10=False, title='Rrs531' + nc_file.time_coverage_end[0:13],
                      save_image='Rrs531' + nc_file.time_coverage_end[0:13])
    gs.plot_geo_image(Rrs_547, lons, lats, log10=False, title='Rrs547' + nc_file.time_coverage_end[0:13],
                      save_image='Rrs547' + nc_file.time_coverage_end[0:13])
    gs.plot_geo_image(Rrs_555, lons, lats, log10=False, title='Rrs555' + nc_file.time_coverage_end[0:13],
                      save_image='Rrs555' + nc_file.time_coverage_end[0:13])
    gs.plot_geo_image(Rrs_645, lons, lats, log10=False, title='Rrs645' + nc_file.time_coverage_end[0:13],
                      save_image='Rrs645' + nc_file.time_coverage_end[0:13])
    gs.plot_geo_image(Rrs_667, lons, lats, log10=False, title='Rrs667' + nc_file.time_coverage_end[0:13],
                      save_image='Rrs667' + nc_file.time_coverage_end[0:13])
    gs.plot_geo_image(Rrs_678, lons, lats, log10=False, title='Rrs678' + nc_file.time_coverage_end[0:13],
                      save_image='Rrs678' + nc_file.time_coverage_end[0:13])



#
# #CHL anamaly
# #Description and implementation
# #Difference in CHL between a single image and the mean over 2 months ending 2 weeks prior to the image.
# # A CHL anomaly of N1 mg m−3 is flagged as a bloom. Not specificto K. brevis blooms.
#
# #bbp ration
# # bbp morel=0.3 * CHL^0.62 *( 0.002 + 0.02 *( 0.5−0.25 * log10(CHL))
# # If CHL is greater than 1.5 mg m−3 and bbp (551) is less than bbp calculated with the Morel (1988) model ,then the pixel is flagged as a K. brevis bloom.
# # Bbp can be calculated using either Carder et al. (1999) or the QAA algorithm (Lee et al., 2002).
    bbp555QAA=np.zeros(np.shape(lats))
    bbp=np.ones(10)
    bbpmorel=0.3*(chl**0.62)*(0.002+0.02*(0.5-0.25*np.log10(chl)))
    for i in range(len(y)):
        for j in range(len(x)):
            Rrs=np.array([
                Rrs_412[i,j],
                Rrs_443[i, j],
                Rrs_469[i, j],
                Rrs_488[i, j],
                Rrs_531[i, j],
                Rrs_547[i, j],
                Rrs_555[i, j],
                Rrs_645[i, j],
                Rrs_667[i, j],
                Rrs_678[i, j]])
            for k in range(10):
                if((Rrs[k]<0)or(Rrs[k]==np.nan) ):
                    bbp[k]=np.nan
                else:
                    bbp = QAAv6(Rrs)
            bbp555QAA[i,j]=bbp[6]

    Imagebbp=np.zeros(np.shape(lats))

    Imagebbp[(chl>1.5) & (bbp555QAA<bbpmorel)]=5 #Karenia
    gs.plot_geo_image(Imagebbp,lons,lats,log10=False,title='Result by Cannizzaro et al(2004,2008,2009), bbpratio '+nc_file.time_coverage_end[0:13],caxis=[0.1,6],save_image='Result by Cannizzaro et al(2004,2008,2009) in Uwajima'+nc_file.time_coverage_end[0:13])
# TODO UwaJIMA再跑一次
#
# #nLw ratio
# #Estimates the relationship between water leaving radiance (Lw(551)) and CHL.
# # If Lw(551) is lower than the backscattering calculated at λ = 550 nm (bbp MOREL(550); Morel, 1988)as function of CHL then the pixel is classified as a K. brevis bloom.
    F0=[172.912,187.622,205.878,194.933,185.747,186.539,183.869,157.811,152.255,148.052]
    nlw555=Rrs_555*F0[6]
# #     Imagenlw=np.zeros(np.shape(lats))
# #TODO 有待商议
# #     Imagenlw[(chl>1.5)&(nlw555<bbpmorel)]=5
# #     gs.plot_geo_image(Imagenlw, lons, lats, log10=False,
# #                       title='Result by Carvalho et al(2008,2011)' + nc_file.time_coverage_end[0:13],
# #                       caxis=[0.1, 6], save_image='Result by Carvalho et al(2008,2011)' + nc_file.time_coverage_end[0:13])
#
#
# #Spectral shape at 490 nm (SS_490)
# #Calculates the spectral shape at 490 nm (SS_490)
# # using nLw: SS_90 = nLw(488) − nLw(443) − (nLw(531) − nLw(443)) × ((488 − 443) ÷ (531 − 443))
# # A negative SS_490 is indicative of a K. brevis bloom.
    nlw488=Rrs_488*F0[3]
    nlw443=Rrs_443*F0[1]
    nlw531=Rrs_531*F0[4]
    SS_90 = nlw488-nlw443-(nlw531-nlw443)*((488-443)/(531-443))
    ImageSS=np.zeros(np.shape(lats))
    ImageSS[SS_90<0]=5
    gs.plot_geo_image(ImageSS, lons, lats, log10=False,
                      title='Result by Tominson et al(2009)' +'\n'+ nc_file.time_coverage_end[0:13],
                      caxis=[0, 6], save_image='Result by Tominson et al(2009)in ,SS' + nc_file.time_coverage_end[0:13])
    #gs.plot_geo_image(ImageEKO,lonmlat)
# #RBD
# #Takes advantage of the high fluorescence properties of dinoflagellate blooms.
# # RBD = nLw(678) − nLw(667) RBD N 0.15 Wm−2 μm−1 sr−1
# # indicates a dinoflagellate bloom but not specifically a K. brevis bloom.
    nlw667=Rrs_667*F0[8]
    nlw678=Rrs_678*F0[9]
    RBD=nlw678-nlw667
    ImageRBD=np.zeros(np.shape(lats))
    ImageRBD[RBD>0.015]=5
    gs.plot_geo_image(ImageRBD, lons, lats, log10=False,
                      title='Result by Amin, Gilerson et al(2009,RBD)' + '\n' + nc_file.time_coverage_end[0:13],
                      caxis=[0.1, 6], save_image='Result by Amin, Gilerson et al(2009,RBD)in Uwajima' + nc_file.time_coverage_end[0:13])
#
#
# #RBD–KBBI
# #Combines the high fluorescence and low backscattering properties of K. brevis blooms in 2 main equations:
# # RBD = nLw(678) − nLw(667) and KBBI = (nLw(678) − nLw(667)) / (nLw(678) + nLw(667))
# # A K. brevis bloom should meet the following criteria:
# # RBD N 0.15 Wm−2 μm−1 sr−1,and KBBI N 0.3*RBD
    KBBI=(nlw678-nlw667)/(nlw678+nlw667)
    ImageKBBI = np.zeros(np.shape(lats))
    ImageKBBI[(RBD > 0.015)&(KBBI>0.3*RBD)] = 5
    gs.plot_geo_image(ImageKBBI,lons, lats, log10=False,
                      title='Result by Amin, Gilerson et al(2009,RBD-KBBI)' + '\n' + nc_file.time_coverage_end[0:13],
                      caxis=[0.1, 6], save_image='Result by Amin, Gilerson et al(2009,RBD-KBBI)in Uwajima' + nc_file.time_coverage_end[0:13])
#EKO
#nlw Peak at 547
#NOT！nlw443-412 slope>nlw488-443slope & nlw488-443>nlw547-488slope
#NOT! nlw547-412 slope>-0.0003[ln(chla)}^2+0.0024[ln(chla)]-0.00005
#NOT! nlw547>0.8
#NOT! nlw488-412slope-nlw547-488slope>0.06
# Lambda=np.array([412, 443, 469,488,531,547,555,645,667,678])
# The most complex one
# How much is chla affected?
    nlw412=Rrs_412*F0[0]
    nlw469=Rrs_469*F0[2]
    nlw547=Rrs_547*F0[5]
    nlw645=Rrs_645*F0[7]
    nlw443_412slope=(nlw443-nlw412)/(443-412)
    nlw488_443slope=(nlw488-nlw443)/(488-443)
    nlw547_488slope=(nlw547-nlw488)/(547-488)
    nlw547_412slope=(nlw547-nlw412)/(547-412)
    nlw488_412slope=(nlw488-nlw412)/(488-412)
    ImageEKO=np.zeros(np.shape(lats))
    # TODO There must be some problem with EKO's method!
    for i in range(len(y)):
        for j in range(len(x)):
            if nlw547[i, j]>max(nlw412[i, j], nlw443[i, j], nlw488[i, j], nlw555[i, j], nlw645[i, j], nlw667[i, j]):
                if (not ((nlw443_412slope[i, j] > nlw488_443slope[i, j]) & (
                        nlw488_443slope[i, j] > nlw547_488slope[i, j]))):
                    if (not (nlw547_412slope[i, j]) > -0.0003 * (
                            (np.log(chl[i, j])) ** 2) + 0.0024 * np.log(chl[i, j]) - 0.00005):
                        if (not (nlw547[i, j] > 0.8)):
                            if (not (abs(nlw547_488slope[i, j] - nlw488_412slope[i, j]) > 0.006)):
                                ImageEKO[i, j] = 5
    #gs.plot_geo_image(ImageEKO, lons, lats, log10=False,
                                   #title='EKO1',
                                   #caxis=[0.1, 6])


                #     if(not(nlw547_412slope[i,j])>-0.0003*((np.log(chl[i,j]))**2+0.0024*np.log(chl[i,j]))-0.00005):
                #         if(not(nlw547[i,j]>0.8)):
                #             if(not(abs(nlw547_488slope[i,j]-nlw488_412slope[i,j])>0.006)):
                #                 ImageEKO[i,j]=5
    gs.plot_geo_image(ImageEKO,lons, lats, log10=False,
                       title='Result by Eko Siswanto et al(2013)' + '\n' + nc_file.time_coverage_end[0:13],
                       caxis=[0, 6], save_image='Result by Eko Siswanto et al(2013)in Uwajima' + nc_file.time_coverage_end[0:13])

    print(print('percent: {:.2%}'.format((a1+1)/L)))

    #nFLH
    #　shen 2019
    #　red tide detection index:rdi=(1/Rrs1-1/Rrs2)*Rrs3
    # (a) the λ1, λ2, and λ3 of Eq. (1) were set using the MERIS bands of 665, 560, and 753 nm,
    # (b) the λ1, λ2, and λ3 were set using the MODIS bands of 667, 555, and 748 nm, and
    # (c) the λ1, λ2, and λ3 were set using the GOCI bands of 660, 555, and 745 nm.
    #　0.14~0.16
    # R_slope =tan−1 (100×(1 – (Rrs(λ2)/ Rrs(λ1))) / (λ2–λ1))
    # λ1 represents the baseline wavelength in the range of 550–570 nm, and λ2 is proposed to be a wavelength in the range of 570–670 nm.
    #rrs slope 0.5 original paper do not included
