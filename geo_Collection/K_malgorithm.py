import numpy as np


def bbpratio(chl, bbp555QAA, nflh, Rrs555):
    # 0 nodata
    # 1 other water
    # 2 other phytoplankton
    # 3k,bloom
    bbpmorel = 0.3 * (chl ** 0.62) * (0.002 + 0.02 * (0.5 - 0.25 * np.log10(chl)))
    Imagebbp = np.ones(np.shape(chl))

    Imagebbp[chl.mask] = 0
    Imagebbp[(nflh > 0.02) & (Rrs555 < 0.007)] = 2
    Imagebbp[(nflh > 0.02) & (Rrs555 < 0.007) & (bbp555QAA < bbpmorel)] = 3  # Karenia
    return Imagebbp


def nlwratio(chl, bbpmorel, nlw555):
    Imagenlw = np.ones(np.shape(chl))
    # #TODO 有待商议
    Imagenlw[chl.mask] = 0
    Imagenlw[(chl > 1.5) & (nlw555 < bbpmorel)] = 3
    return Imagenlw


def SS(nlw488, nlw443, nlw531, chl):
    SS_90 = nlw488 - nlw443 - (nlw531 - nlw443) * ((488 - 443) / (531 - 443))
    ImageSS = np.ones(np.shape(chl))
    ImageSS[chl.mask] = 0
    ImageSS[SS_90 < 0] = 3
    return ImageSS


def RBDKBBI(nlw667, nlw678, chl):
    RBD = nlw678 - nlw667
    KBBI = (nlw678 - nlw667) / (nlw678 + nlw667)

    ImageKBBI = np.ones(np.shape(chl))

    ImageKBBI[chl.mask] = 0

    ImageKBBI[RBD > 0.015] = 2
    ImageKBBI[(RBD > 0.015) & (KBBI > 0.3 * RBD)] = 3
    return ImageKBBI


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
    return ImageBI


def shen(Rrs_555, Rrs_667, Rrs_748):
    RDI = (1 / Rrs_667 - 1 / Rrs_555) * Rrs_748
    ImageShen = np.ones(np.shape(Rrs_555))
    ImageShen[RDI > 0.16] = 3
    # R_slope = np.arctan(100 * (1 - ((Rrs_555 / Rrs_667) / (555 - 667))))
    # ImageShen[(RDI > 0.14) & (R_slope < 0.5)] = 5
    return ImageShen


def siswanto(x, y, chl, nlw412, nlw443, nlw488, nlw547, nlw555, nlw645, nlw667):
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
    return ImageEKO
