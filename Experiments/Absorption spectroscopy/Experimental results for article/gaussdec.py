# -*- coding: utf-8 -*-
"""
gaussdec.py

@author: Aurelien Legoupil
"""

#%% Imports
import numpy as np
from scipy.optimize import curve_fit

#%% Catalog of defect peaks
# Nameofdefect, peak number, OA peak eV , FWHM eV, sources (list)

#GR91 = Griscom, 1991
#GI19 = Girard, 2019
#KA22 = Kashaykin, 2022
#MO20 = Morana, 2020

#{"peak":,"FWHM":,"sources":["GI19"]}
#FWHM=1e-6 when unknown 

DICOPEAKS = {
    "ODC-II":[
        {"peak":5.05,"FWHM":0.32,"sources":["GI19"]},
        {"peak":3.15,"FWHM":0.30,"sources":["GI19"]},
        {"peak":6.9,"FWHM":0.4,"sources":["GI19"]}
        ],
    "ODC-I":[
        {"peak":7.6,"FWHM":{"min":0.5,"max":0.6},"sources":["GI19"]}
        ],
    "NBOHC":[
        {"peak":1.97,"FWHM":0.17,"sources":["GI19"]},
        {"peak":4.8,"FWHM":1.0,"sources":["GI19"]},
        {"peak":6.4,"FWHM":1.7,"sources":["GI19"]}
        ],
    "E-prime":[
        {"peak":5.8,"FWHM":0.7,"sources":["GI19"]}
        ],
    "STH1":[
        {"peak":2.61,"FWHM":1.2,"sources":["GI19"]},
        {"peak":1.88,"FWHM":{"min":0.2,"max":0.5},"sources":["GI19"]}
        ],
    "STH2":[
        {"peak":2.16,"FWHM":{"min":0.3,"max":0.6},"sources":["GI19"]},
        {"peak":1.63,"FWHM":{"min":0.3,"max":0.7},"sources":["GI19"]}
        ],
    "POL":[
        {"peak":3.8,"FWHM":0.2,"sources":["GI19"]},
        {"peak":4.2,"FWHM":0.6,"sources":["GI19"]},
        {"peak":7.3,"FWHM":0.2,"sources":["GI19"]},
        {"peak":7.5,"FWHM":0.1,"sources":["GI19"]}
        ],
    "POR":[
        {"peak":2.02,"FWHM":None,"sources":["GI19"]},
        {"peak":4.08,"FWHM":None,"sources":["GI19"]},
        {"peak":5.02,"FWHM":None,"sources":["GI19"]},
        {"peak":2.0,"FWHM":None,"sources":["GI19"]},
        {"peak":4.8,"FWHM":None,"sources":["GI19"]}
        ],
    "STE":[
        {"peak":3.7,"FWHM":None,"sources":["GI19"]},
        {"peak":4.6,"FWHM":None,"sources":["GI19"]},
        {"peak":6.4,"FWHM":None,"sources":["GI19"]}
        ],
    "STEX":[
        {"peak":4.2,"FWHM":1.16,"sources":["GI19"]},
        {"peak":5.3,"FWHM":0.78,"sources":["GI19"]}
        ],
    "Composite 1-eV band":[
        {"peak":1.2,"FWHM":0.56,"sources":["KA22","MO20"]},
        {"peak":0.93,"FWHM":0.42,"sources":["KA22","MO20"]}
        ]
    }

DICOIMPUR = {
    "O3":[
        {"peak":4.8,"FWHM":{"min":0.8,"max":0.86},"sources":["GI19"]}
        ],
    "O2":[
        {"peak":0.97,"FWHM":0.013,"sources":["GI19"]},
        {"peak":1.62,"FWHM":0.013,"sources":["GI19"]}
        ],
    "Cl0":[
        {"peak":3.26,"FWHM":None,"sources":["GI19"]},
        {"peak":3.65,"FWHM":None,"sources":["GI19"]}
        ],
    "Cl2":[
        {"peak":3.78,"FWHM":0.6,"sources":["GI19"]},
        {"peak":2.3,"FWHM":None,"sources":["GI19"]}
        ],
    "H(I)":[], #not observed
    "LTIRA":[] #non gaussian, see Girard 2019, Dianov Chernov 1989  
    }

DICOLTIRA = {"D":0.39, #eV (even if they say eV^-1?) #From Ka22 and Dianov-Chernov 89
             "S":0.11, #eV (even if they say eV^-1?)
             "R":0.43, #eV (even if they say eV^-1?)
             "Ep":0.65 #eV
             }

#%% Helper functions

def Sigma(FWHM):
    return FWHM/(2*np.sqrt(2*np.log(2)))

eV = 1.60218e-19 #J
h_Planck = 6.62607015e-34 #m^2.kg/s
c_light = 299792458 #m/s
def eV2WLnm(E):
    return h_Planck*c_light/(E*eV*1e-9) #in nm
def WLnm2eV(lmbda):
    return h_Planck*c_light/(lmbda*1e-9*eV)


#%% Functions for bands

def LTIRAfunc(lmbda,magnitude,D,S,R,Ep,model="Kubo-Greenwood"):
    E = WLnm2eV(lmbda)
    if magnitude==0:
        return 0
    if model=="Kubo-Greenwood":
        Q=D+S/(D-S)
        B0=magnitude
        return B0*E*np.exp(-E/S)*(np.exp(E/Q)-np.exp(Ep/Q))
    elif model=="exp":
        C=magnitude
        return (C*np.exp(-E/R))
    else:
        raise NameError("You fucked up")
    
def gaussband(lmbda,magnitude,peak_eV,FWHM_eV):
    E = WLnm2eV(lmbda)
    sigma = Sigma(FWHM_eV)
    return magnitude/(sigma*np.sqrt(2*np.pi)) * np.exp(-(E-peak_eV)**2/(2*sigma**2))

#%% Function for total absorption fitting

def totabs(lmbda,
           ltira_magnitude,ltira_D,ltira_S,ltira_R,ltira_Ep,
           *paramsgauss):
    
    assert len(paramsgauss)%3==0
    n_gaussians=int(len(paramsgauss)/3)
    
    total = LTIRAfunc(lmbda, ltira_magnitude, ltira_D, ltira_S, ltira_R, ltira_Ep)
    #Add gaussian bands:
    for i in range(n_gaussians):
        magnitude = paramsgauss[i*3]
        peak_eV = paramsgauss[i*3+1]
        FWHM_eV = paramsgauss[i*3+2]
        total += gaussband(lmbda, magnitude, peak_eV, FWHM_eV)
    
    return total

#%%

def decompose(lmbdadata,ydata,LTIRA=False,GaussList=[],fixedreldelta=1e-6,defautfwhm_minmaxinit=(1e-6,2,0.1),relmaxfev=100,defaultmag=1e2):
    #GaussList must be a list of keys in DICOPEAKS or DICOIMPUR, setting these gaussians as available for decomposition
    #Possible improvements if necessary after decomposition attempt:
        #give margin for ltira parameters and other fixed parameters (in percent?)
        #modify the optimization ranges for peaks with unknown FWHM
    
    #Construct the intial guess and parameter ranges for curve fitting
    
    d_minus = 1-fixedreldelta/2
    d_plus  = 1+fixedreldelta/2
    
    low_bounds=[]
    high_bounds=[]
    if LTIRA:
        initial_guesses =[defaultmag,DICOLTIRA["D"],DICOLTIRA["S"],DICOLTIRA["R"],DICOLTIRA["Ep"]]
        low_bounds=[0,DICOLTIRA["D"]*d_minus,DICOLTIRA["S"]*d_minus,DICOLTIRA["R"]*d_minus,DICOLTIRA["Ep"]*d_minus]
        high_bounds=[np.inf,DICOLTIRA["D"]*d_plus,DICOLTIRA["S"]*d_plus,DICOLTIRA["R"]*d_plus,DICOLTIRA["Ep"]*d_plus]
    else:
        initial_guesses =[0,DICOLTIRA["D"],DICOLTIRA["S"],DICOLTIRA["R"],DICOLTIRA["Ep"]]
        low_bounds=[0-fixedreldelta/2,DICOLTIRA["D"]*d_minus,DICOLTIRA["S"]*d_minus,DICOLTIRA["R"]*d_minus,DICOLTIRA["Ep"]*d_minus]
        high_bounds=[0+fixedreldelta/2,DICOLTIRA["D"]*d_plus,DICOLTIRA["S"]*d_plus,DICOLTIRA["R"]*d_plus,DICOLTIRA["Ep"]*d_plus]
    
    for KeyGauss in GaussList:
        if KeyGauss in DICOPEAKS.keys():
            PEAKS=DICOPEAKS[KeyGauss]
        elif KeyGauss in DICOIMPUR.keys():
            PEAKS=DICOIMPUR[KeyGauss]
        else:
            PEAKS=[]
        for dicpeak in PEAKS:
            if dicpeak["FWHM"]==None:
                initial_guesses+=[defaultmag,dicpeak["peak"],defautfwhm_minmaxinit[2]]
                low_bounds+=[0,dicpeak["peak"]*d_minus,defautfwhm_minmaxinit[0]]
                high_bounds+=[np.inf,dicpeak["peak"]*d_plus,defautfwhm_minmaxinit[1]]
            elif type(dicpeak["FWHM"])==float:
                initial_guesses+=[defaultmag,dicpeak["peak"],dicpeak["FWHM"]]
                low_bounds+=[0,dicpeak["peak"]*d_minus,dicpeak["FWHM"]*d_minus]
                high_bounds+=[np.inf,dicpeak["peak"]*d_plus,dicpeak["FWHM"]*d_plus]
            elif type(dicpeak["FWHM"])==dict:
                initial_guesses+=[defaultmag,dicpeak["peak"],(dicpeak["FWHM"]["min"]+dicpeak["FWHM"]["max"])/2]
                low_bounds+=[0,dicpeak["peak"]*d_minus,dicpeak["FWHM"]["min"]]
                high_bounds+=[np.inf,dicpeak["peak"]*d_plus,dicpeak["FWHM"]["max"]]
  
    #Apply curve fitting
    popt, pcov = curve_fit(totabs, lmbdadata, ydata, p0=initial_guesses, bounds=(low_bounds,high_bounds),max_nfev=relmaxfev*len(initial_guesses))
    
    return popt,pcov







