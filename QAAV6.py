from scipy.interpolate import Akima1DInterpolator as Akima
import numpy as np
from deco import *
import math as m


# wavelengths = {
#     'OLI'   : [442.98, 482.49, 561.33, 654.61],
#     'MSI'   : [443.93, 496.54, 560.01, 664.45],
#     'OLCI'  : [411.3999939, 442.63000488, 490.07998657, 510.07000732, 560.05999756, 619.97998047, 664.85998535, 673.61999512, 681.15002441], # 9 band insitu
#     'OLCI2' : [400, 412.5, 442.5, 490, 510, 560, 620, 665, 673.75, 681.25, 708.75, 753.75, 761.25, 764.375, 767.5, 778.75], # 16 band LUT
#     'VI'    : [412.49, 444.17, 486.81, 549.99, 670.01],
#     'AER'   : [412, 442, 490, 530, 551, 668],
#     'MOSIA' : [412, 443, 469,488,531,547,555,645,667,678],
#     'GOCI'  : [412, 443, 490, 555, 660, 680],
#     'SGLI'  : [380,412,443,490,530,565,673.5]
#}
#@concurrent
def QAAv6MODIS(Rrs):
    #acording to Lee,QAAv6
    #write by Zhou,20191130
    #Input data need to be an arrary that contain 10 bands Rrs of MODIS,from short wavelength to long wavelength in a certain station
    #Output is a tuple, first array is aph,second array is bbp
    #use as import QAAV6
# B1: 412 nm 0
# B2: 443 nm 1
# B3: 469 nm 2
# B4: 488 nm 3
# B5: 531 nm 4
# B6: 547 nm 5
# B7: 555 nm 6
# B8: 645 nm 7
# B9: 667 nm 8
# B10: 678 nm 9

    Lambda=np.array([412, 443, 469,488,531,547,555,645,667,678])
    nbands=np.shape(Lambda)[0]

    IOPw=np.array([[0.003344468,0.004572564],
    [0.00244466,0.00635],
    [0.001910803,0.010483637],
    [0.001609567,0.014361745],
    [0.00111757,0.043747657],
    [0.000983055,0.053262848],
    [0.000923288,0.0595],
    [0.000482375,0.325],
    [0.00041731,0.433497423],
    [0.00038884,0.457440162]])

    #     if(Rrs=np.nan):
    #     return np.nan
    # else:
    # bbw from Morel (1974).aw  from Pope and Fry (1997)
    bbp = np.ones(10)
    adg = np.ones(10)
    if(np.nan in Rrs):
        bbp[:]=np.nan
    else:
        bw=IOPw[:,0]#backscaterring of pure water
        aw=IOPw[:,1]#absorption of pure water
        rrs = Rrs / (0.52 + 1.7 * Rrs)
        g0 = 0.089
        g1 = 0.1245
        u = (-g0 + ((g0 ** 2) + 4 * g1 * rrs) ** 0.5) / (2 * g1)

        aph = np.ones(10)#adg is the absorption of CDOM and NAP
        if Rrs[6]<0.0015:#select 555 as reference
            r=550
            p1=(rrs[1] + rrs[3])
            p2 = rrs[6] + 5 * (((rrs[8]) ** 2)) / (rrs[3])
            x = np.log10(p1 / p2)
            ar = aw[6] + np.power(10, (-1.146 - 1.366 * x - 0.469 * (x ** 2)))# step 2
            bbpr=((u[6]*ar)/(1-u[6]))-bw[6]#step3
        else:
            r=670
            p1 = Rrs[8] / (Rrs[1] + Rrs[3])
            p2 = 0.39 * (p1 ** 1.14)
            ar = (aw[8]) + p2  # step2
            bbpr = (u[8] * ar / (1 - (u[8])) - (bw[8]))  # step3
        eta=2*(1-1.2*np.exp(-0.9*(rrs[1]/rrs[6]))) #step4

        zeta = 0.74 + 0.2 / (0.8 + rrs[1] / rrs[6])#step 7&8
        S = 0.015 + 0.002 / (0.6 + rrs[1] / rrs[6])
        xi = np.exp(S * (442.5 - 415.5))
        for i in range(nbands):
            bbp[i]= bbpr * np.power(r/Lambda[i], eta)#step5
        a = ((1 - u) * (bw + bbp)) / u#step6
        for i in range(nbands):

            ag443=((a[0]-zeta*a[1])/(xi-zeta))-((aw[0]-zeta*aw[1])/(xi-zeta))
            adg[i]=ag443*np.exp(-S*(Lambda[i]-443))
            aph[i]=a[i]-adg[i]-aw[i]
        return bbp,a


def QAAv6GOCI(Rrs):
    #acording to Lee,QAAv6
    #write by Zhou,20191130
    #Input data need to be an arrary that contain 10 bands Rrs of GOCI,from short wavelength to long wavelength in a certain station
    #Output is a list, first array is aph,second array is bbp
    #use as import QAAV6
# B1: 412 nm 0
# B2: 443 nm 1
# B3: 490 nm 2
# B4: 555 nm 3
# B5: 660 nm 4
# B6: 680 nm 5
    Lambda=np.array([412, 443, 490, 555, 660, 680])
    nbands=np.shape(Lambda)[0]

    IOPw=np.array([[0.003344468,0.004572564412],
    [0.00244466,0.00635443],
    [0.001609567,0.014361745488],
    [0.000923288,0.0595555],
    [0.00041731,0.433497423667],
    [0.00038884,0.457440162678]])
    # bbw from Morel (1974).aw  from Pope and Fry (1997)
    bw=IOPw[:,0]#backscaterring of pure water
    aw=IOPw[:,1]#absorption of pure water
    rrs=Rrs/(0.52+1.7*Rrs)
    g0=0.089
    g1=0.1245
    u=(-g0+((g0**2)+4*g1*rrs)**0.5)/(2*g1)

    bbp=np.ones(6)
    adg=np.ones(6)
    aph=np.ones(6)#adg is the absorption of CDOM and NAP
    if ((Rrs[5])/2)<0.0015:#select 555 as reference#用680代替
        r=550
        p1=(rrs[1]+rrs[2])
        p2=rrs[3]+5*(((rrs[5])**2))/(rrs[2])
        x=np.log10(p1/p2)
        ar=aw[3]+np.power(10,(-1.146-1.366*x-0.469*(x**2)))# step 2
        bbpr=((u[3]*ar)/(1-u[3]))-bw[3]#step3
        #TODO bbp每个都多了0.05左右
    else:
        r=670
        p1=Rrs[5]/(Rrs[1]+Rrs[2])
        p2=0.39*(p1**1.14)
        ar=(aw[5])+p2 #step2
        bbpr=(u[5]*ar/(1-(u[5]))-(bw[5]))# step3
    eta=2*(1-1.2*np.exp(-0.9*(rrs[1]/rrs[3]))) #step4

    zeta = 0.74 + (0.2 / (0.8 + (rrs[1] / rrs[3])))#step 7&8
    S = 0.015 + 0.002 / (0.6 + (rrs[1] / rrs[3]))
    xi = np.exp(S * ((442.5 - 415.5)))

    for i in range(nbands):
        bbp[i]= bbpr * np.power(r/Lambda[i], eta)#step5

    a = ((1 - u) * (bw + bbp)) / u
#TODO 现在确定了 就是a算错了 别的都没有问题
    for i in range(nbands):
        ag443=((a[0]-zeta*a[1])/(xi-zeta))-((aw[0]-zeta*aw[1])/(xi-zeta))
        adg[i]=ag443*np.exp(-S*(Lambda[i]-443))
        aph[i]=a[i]-adg[i]-aw[i]

    return bbp





















