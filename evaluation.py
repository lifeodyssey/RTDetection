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
#todo 读取文件，计算出各个区域match的case，算出月份
result500=pd.read_excel('/Users/zhenjia/Desktop/resul_500.xlsx',sheet_name='Sheet3')


def perf_measure(y_actual, y_hat):
    TP = 0
    FP = 0
    TN = 0
    FN = 0

    for i in range(len(y_hat)):
        if y_actual[i] == y_hat[i] == 1:
            TP += 1
        if y_hat[i] == 1 and y_actual[i] != y_hat[i]:
            FP += 1
        if y_actual[i] == y_hat[i] == 0:
            TN += 1
        if y_hat[i] == 0 and y_actual[i] != y_hat[i]:
            FN += 1

    return TP, FP, TN, FN


def perf(TP, FP, TN, FN):
    # Sensitivity, hit rate, recall, or true positive rate
    try:
        TPR = TP / (TP + FN)
    except:
        TPR = 0
    # # Specificity or true negative rate
    # TNR = TN / (TN + FP)
    # # Precision or positive predictive value
    PPV = TP / (TP + FP)
    # # Negative predictive value
    # NPV = TN / (TN + FN)
    # # Fall out or false positive rate
    # FPR = FP / (FP + TN)
    # # False negative rate
    # FNR = FN / (TP + FN)
    # # False discovery rate
    # FDR = FP / (TP + FP)

    # Overall accuracy
    ACC = (TP + TN) / (TP + FP + FN + TN)
    try:
        f1 = (2 * PPV * TPR) / (PPV + TPR)
    except:
        f1 = 0
    return ACC


def perf_by_region(sheetname):
    sheet = pd.read_excel('/Users/zhenjia/Desktop/Project/Seto/result 6/Field.xlsx', sheet_name=sheetname)
    Oitaest = perf_measure(sheet['East of Oita'], sheet['East of Oita.1'])
    Oitanorth = perf_measure(sheet['North of Oita'], sheet['North of Oita.1'])
    Uwajimaco = perf_measure(sheet['Uwajima coast'], sheet['Uwajima coast.1'])
    Osaka = perf_measure(sheet['Osaka'], sheet['Osaka.1'])

    P_EO = perf(Oitaest[0], Oitaest[1], Oitaest[2], Oitaest[3])
    P_NO = perf(Oitanorth[0], Oitanorth[1], Oitanorth[2], Oitanorth[3])
    P_UC = perf(Uwajimaco[0], Uwajimaco[1], Uwajimaco[2], Uwajimaco[3])
    P_OS = perf(Osaka[0], Osaka[1], Osaka[2], Osaka[3])
    return P_EO, P_NO, P_UC, P_OS


list = ['sisiwanto', 'rbd-kbbi', 'shang_need', 'bbp ratio']
l = len(list)
for i in range(l):
    print(list[i])
    print('P_EO, P_NO, P_UC,P_OS')
    P_EO, P_NO, P_UC, P_OS = perf_by_region(list[i])
    print(P_EO, P_NO, P_UC, P_OS)
