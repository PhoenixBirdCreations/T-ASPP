#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 16 11:38:53 2022

@author: mberbel
"""

import numpy as np
import matplotlib.pyplot as plt

#%%

def getData(table):
    data = np.loadtxt(table)
    N = len(data)
    rho = np.zeros(N); P = np.zeros(N); epsAd = np.zeros(N);
    for i in range (0,N):
        rho[i] = data[i][0]
        P[i] = data[i][1]
        epsAd[i] = (data[i][2]/data[i][0]-1)
    fend = {}
    fend['rho'] = rho
    fend['P'] = P
    fend['epsAd'] = epsAd
    return fend

#Derivatives of Lagrangian interpolator
def der3_middle(x0, x1, x2, y0, y1, y2):
    x12 = x1 - x2
    x01 = x0 - x1
    x02 = x0 - x2
    return y0*x12/(x01*x02)+y1*(1.0/x12-1.0/x01)-y2*x01/(x02*x12)

def der3_initial(x0, x1, x2, y0, y1, y2):
    x12 = x1 - x2
    x01 = x0 - x1
    x02 = x0 - x2
    return y0*(x01+x02)/(x01*x02)-y1*x02/(x01*x12)+y2*x01/(x02*x12)

def der3_last(x0, x1, x2, y0, y1, y2):
    x12 = x1 - x2
    x01 = x0 - x1
    x02 = x0 - x2
    return -y0*x12/(x01*x02)+y1*x02/(x12*x01)-y2*(x02+x12)/(x02*x12)

def calculateTermoQuantities(data):
    N = len(data['P'])
    P = data['P']
    rho = data['rho']
    
    #Classic sound speed
    csc = np.zeros(N)
    csc[0] = der3_initial(rho[0], rho[1], rho[2], P[0], P[1], P[2])
    for i in range(1,N-1):
        csc[i] = der3_middle(rho[i-1], rho[i], rho[i+1], P[i-1], P[i], P[i+1])
    csc[N-1] = der3_last(rho[N-3], rho[N-2], rho[N-1], P[N-3], P[N-2], P[N-1])
    
    data['csc'] = csc
    
    
    #Adiabatic exponent
    LnP = np.log(data['P']); Lnrho = np.log(data['rho'])
    cGd = np.zeros(N)
    cGd[0] = der3_initial(Lnrho[0], Lnrho[1], Lnrho[2], LnP[0], LnP[1], LnP[2])
    for i in range(1,N-1):
        cGd[i] = der3_middle(Lnrho[i-1], Lnrho[i], Lnrho[i+1], LnP[i-1], LnP[i], LnP[i+1])
    cGd[N-1] = der3_last(Lnrho[N-3], Lnrho[N-2], Lnrho[N-1], LnP[N-3], LnP[N-2], LnP[N-1])

    data['adexp'] = cGd
    
    
    #Fundamental derivative
    csc = np.where(csc<0, 1e-13, csc) #fix floating point 0 values for log
    Lncsc = 0.5*np.log(csc)
    G = np.zeros(N)
    G[0] = 1 + der3_initial(Lnrho[0], Lnrho[1], Lnrho[2], Lncsc[0], Lncsc[1], Lncsc[2])
    for i in range(1,N-1):
        G[i] = 1 + der3_middle(Lnrho[i-1], Lnrho[i], Lnrho[i+1], Lncsc[i-1], Lncsc[i], Lncsc[i+1])
    G[N-1] = 1 + der3_last(Lnrho[N-3], Lnrho[N-2], Lnrho[N-1], Lncsc[N-3], Lncsc[N-2], Lncsc[N-1])
    
    data['G'] = G
    
    
    #Relativistic sound speed
    data['cs'] = csc/(1 + data['epsAd'] + data['P']/data['rho'])
    
    
    #Grüneisen coefficient
    eps = data['epsAd']
    Gru = np.zeros(N)
    Gru[0] = der3_initial(eps[0], eps[1], eps[2], P[0], P[1], P[2])/rho[0]
    for i in range(1,N-1):
        Gru[i] = der3_middle(eps[i-1], eps[i], eps[i+1], P[i-1], P[i], P[i+1])/rho[i]
    Gru[N-1] = der3_last(eps[N-3], eps[N-2], eps[N-1], P[N-3], P[N-2], P[N-1])/rho[N-1]
    
    cs = data['cs']
    for i in range (0,N):
        Gru[i] = Gru[i]/(Gru[i] - cs[i])
        
    data['Gru'] = Gru
    
    
    return data    

def locateBadPressure(data):
    P = data['P']
    indexes = []
    N = len(P)
    for i in range(0,N-1):
        if (P[i]>P[i+1]):
            print("Pressure decreasing at line ", i + 1)
            indexes.append(i + 1)
    return indexes

def fixBadPressure(data, indexes): #linear interpolation. Suggest changing to log
    P = data['P']
    rho = data['rho']
    for j in indexes:
        data['P'][j] = P[j-1] + (P[j+1] - P[j-1])/(rho[j+1] - rho[j-1])*(rho[j] - rho[j-1])
    return data

def locateBadEpsilon(data):
    eps = data['epsAd']
    indexes = []
    indexes = np.where(eps<0)[0]
    j=-1
    if (len(indexes)>0):
        print("Negative epsilon found until line", indexes[-1])
        j = indexes[-1]
    return j

def fixBadEpsilon(data, j):
    P = data['P']
    rho = data['rho']
    if j==0:
        data['epsAd'][0] = min(P[0]/rho[0], data['epsAd'][1]/2)
    data['epsAd'][0] = P[0]/rho[0]
    for i in range(1,j+1):
        data['epsAd'][i] = data['epsAd'][i-1]+P[i]/rho[i]**2*(rho[i]-rho[i-1])
    return data

def exportData(imin, data, name):
    file=open(name,"w")
    print("#1-rho(g/cm3) 2-G 3-adExp 4-P(g/cm3) 5-eps(ad) 6-cs^2(r) 7-cs^2(c) 8-GruCoeff ", file=file)
    for i in range (imin, len(data['G'])):
        print('{0:1.9e}'.format(data['rho'][i]),"\t", '{0:1.9e}'.format(data['G'][i]),"\t",'{0:1.9e}'.format(data['adexp'][i]),"\t", '{0:1.9e}'.format(data['P'][i]),"\t",'{0:1.9e}'.format(data['epsAd'][i]),"\t", '{0:1.9e}'.format(data['cs'][i]),"\t",'{0:1.9e}'.format(data['csc'][i]),"\t", '{0:1.9e}'.format(data['Gru'][i]),file=file)
    file.close()
   

def doEverything(path,eos):
    dic=getData(path+eos+'\\'+eos+'.table')
    
    index_bad_P = locateBadPressure(dic)
    if (len(index_bad_P)>0):
        dic = fixBadPressure(dic, index_bad_P)
    max_index_bad_eps = locateBadEpsilon(dic)
    if (max_index_bad_eps>=0):
        dic = fixBadEpsilon(dic, max_index_bad_eps)

    dic = calculateTermoQuantities(dic)
    
    exportData(0, dic, eos+'.txt')

#%%

path="/home/mberbel/Escritorio/PhDCopy/95_coldEoSTablesArizona/"
EOS=['BGN1H1', 'ALF1', 'ALF2']

#From Arizona tables with header
#
#   rho [g/cm^3]  P [g/cm^3]   eps [g/cm^3]
#

for eos in EOS:
    print(eos)
    doEverything(path, eos)
