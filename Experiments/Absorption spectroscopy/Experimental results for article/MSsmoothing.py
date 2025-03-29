# -*- coding: utf-8 -*-
"""
MSsmoothing.py

Package for smoothing according to Schmid et al. 2022 "Why and How Savitzky−Golay Filters Should Be Replaced"
https://doi.org/10.1021/acsmeasuresciau.1c00054

@author: Aurelien Legoupil
"""

import numpy as np

CorrCoeffsMS_dic={
    2:[],
    4:[],
    6:[[0.001717576, 0.02437382, 1.64375]],
    8:[[0.0043993373, 0.088211164, 2.359375], #j=0
       [0.006146815, 0.024715371, 3.6359375]],#j=1
    10:[[0.0011840032, 0.04219344, 2.746875], #j=0
        [0.0036718843, 0.12780383, 2.7703125]]#j=1    
}

def windowsMS(x,alpha):
    '''
    w_a(x) from (4)
    '''
    return np.exp(-alpha*x**2) + np.exp(-alpha*(x+2)**2) + np.exp(-alpha*(x-2)**2) - 2*np.exp(-alpha) - np.exp(-9*alpha)

def kernelMS(n,m,alpha=4):
    '''
    Function that returns the convolution kernel array a (eq (3) schmid et al)
    Inputs:
        - m : half size of kernel (kernel size: 2m+1
        - n : degree
        - alpha : steepness parameter. value 4 from 2.2 in schmid et al
    '''
    KAPPA = [ a + b/(c-m)**3 for [a,b,c] in CorrCoeffsMS_dic[n] ]
    X = [i/(m+1) for i in range(-m,m+1)]
    nu = 2 - ((n/2)%2)
    A = np.array([windowsMS(x,alpha) * (
                    np.sinc((n+4)/2*np.pi*x)
                 + np.sum([ kappa*x*np.sin((2*j+nu)*np.pi*x) for j,kappa in enumerate(KAPPA) ])
    ) for x in X])
    A /= np.sum(A)
    return A

def edgeWeights(n,m):
    '''
    Hann-square weights for linear fit at the edges from (17) and (18)
    Inputs:
        - m : half size of kernel (kernel size: 2m+1
        - n : degree
    '''
    beta = 0.7 + 0.14*np.exp(-0.6*(n-4))
    fitlengthD = (m+1)*beta/(1.5+0.5*n)
    fitlength = int(np.floor(fitlengthD))
    return np.array([ (np.cos(np.pi/2*i/fitlengthD))**2  for i in range(fitlength+1) ])

def fitWeighted(xData,yData,weights):
    sumWeights = np.sum(weights)
    sumX = np.sum(xData*weights)
    sumY = np.sum(yData*weights)
    sumX2 = np.sum(xData**2*weights)
    sumXY = np.sum(xData*yData*weights)
    varX2 = sumX2*sumWeights - sumX**2
    if varX2==0:
        slope=0
    else:
        slope = (sumXY*sumWeights - sumX*sumY)/varX2
    offset = (sumY-slope*sumX)/sumWeights
    return offset,slope

def extendData(data,m,fitWeights):
    '''
    Extends the data by weighted linear extrapolation, for smoothing to the ends
    '''
    extData=np.zeros(len(data)+2*m)
    extData[m:len(data)+m]=data
    fitX = np.arange(len(fitWeights))
    fitY = data[:len(fitWeights)]
    offset,slope = fitWeighted(fitX,fitY,fitWeights)
    extData[:m] = [offset + slope*x for x in range(-m+1,1)]
    fitY = data[-len(fitWeights):][::-1]
    offset,slope = fitWeighted(fitX,fitY,fitWeights)
    extData[-m:]= [offset + slope*x for x in range(-m+1,1)][::-1]
    return extData

def smoothMS(data,n,m):
    '''
    MS data smoothing (see Schmid et al) with correction and linear extrapolation for edge smoothing
    Inputs:
        - data : row vector on which to apply smoothing
        - m : half size of kernel (kernel size: 2m+1
        - n : degree
    '''
    kernel = kernelMS(n,m)
    fitWeights = edgeWeights(n,m)
    extData = extendData(data,m,fitWeights)
    smoothedExtData = np.convolve(extData, kernel, mode="same")
    smoothedData = smoothedExtData[m:len(data)+m]
    return smoothedData